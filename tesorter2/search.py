"""
HMM search engine with IQR-balanced parallelism.

Core search function: legacy_search() runs nobias hmmsearch against all
models with IQR-based parallel strategy selection (normal models use
parallel='queries', M^2 outliers use parallel='targets').

Also contains pass1_screen and pass2_search for the deprecated two-pass
mode, retained for potential future use.

Results are captured via TopHits.write() to BytesIO and parsed from
the standard HMMER domain table text format.
"""

import io
import logging
import multiprocessing
import statistics
from collections import defaultdict

import pyfastx
import pyhmmer
import pyhmmer.easel as easel

from .hmm import AMINO_ALPHABET

log = logging.getLogger(__name__)


def tophits_to_domtbl(top_hits, header=False):
    """
    Write TopHits to HMMER domain table format string via BytesIO.
    """
    buf = io.BytesIO()
    top_hits.write(buf, format="domains", header=header)
    return buf.getvalue().decode()


def parse_domtbl_line(line):
    """
    Parse a single HMMER domain table line into a dict.
    Returns None for comment/empty lines.
    """
    if not line or line.startswith("#"):
        return None

    cols = line.split()
    if len(cols) < 22:
        return None

    return {
        "target_name":  cols[0],
        "target_acc":   cols[1],
        "target_len":   int(cols[2]),
        "query_name":   cols[3],
        "query_acc":    cols[4],
        "query_len":    int(cols[5]),
        "evalue":       float(cols[6]),
        "score":        float(cols[7]),
        "bias":         float(cols[8]),
        "dom_num":      int(cols[9]),
        "dom_of":       int(cols[10]),
        "c_evalue":     float(cols[11]),
        "i_evalue":     float(cols[12]),
        "dom_score":    float(cols[13]),
        "dom_bias":     float(cols[14]),
        "hmm_from":     int(cols[15]),
        "hmm_to":       int(cols[16]),
        "ali_from":     int(cols[17]),
        "ali_to":       int(cols[18]),
        "env_from":     int(cols[19]),
        "env_to":       int(cols[20]),
        "acc":          float(cols[21]),
        "description":  " ".join(cols[22:]),
    }


def parse_domtbl_text(text):
    """Parse multi-line domain table text into a list of hit dicts."""
    hits = []
    for line in text.strip().split("\n"):
        hit = parse_domtbl_line(line)
        if hit is not None:
            hits.append(hit)
    return hits


def build_sequence_block(aa_fasta, alphabet=None):
    """
    Build a DigitalSequenceBlock from a FASTA file.
    """
    if alphabet is None:
        alphabet = AMINO_ALPHABET

    seqs = []
    for rec in pyfastx.Fasta(aa_fasta, build_index=True):
        ts = easel.TextSequence(
            name=rec.name.encode("ascii"),
            sequence=str(rec.seq),
        )
        seqs.append(ts.digitize(alphabet))

    return easel.DigitalSequenceBlock(alphabet, seqs)


def _collect_hits(top_hits_iter):
    """Collect all domain hits from an hmmsearch iterator."""
    all_hits = []
    for top_hits in top_hits_iter:
        if len(top_hits) == 0:
            continue
        text = tophits_to_domtbl(top_hits, header=False)
        all_hits.extend(parse_domtbl_text(text))
    return all_hits


def _normalize_nucl_hit(hit):
    """
    Put an nhmmer hit into the convention the rest of the pipeline expects.

    nhmmer reports descending coordinates for minus-strand hits. Sort them
    ascending and encode the strand in the target suffix, the same way
    bath_search does, so _parse_frame_info recovers the base sequence and
    strand. Frame is always 1 here: nucleotide models are not translated.
    """
    af, at = hit["ali_from"], hit["ali_to"]
    ef, et = hit["env_from"], hit["env_to"]
    hit["target_name"] = f"{hit['target_name']}|{'fwd1' if af <= at else 'rev1'}"
    hit["ali_from"], hit["ali_to"] = sorted((af, at))
    hit["env_from"], hit["env_to"] = sorted((ef, et))
    return hit


# E-value fields reported per hit; all share the same search-space scaling.
_EVALUE_FIELDS = ("evalue", "c_evalue", "i_evalue")


def legacy_search_nucl(hmms, seq_block, megabases=None, n_workers=0):
    """
    Single-pass nobias search of DNA models against nucleotide sequence.

    Uses nhmmer, the DNA-DNA tool, rather than hmmsearch. hmmsearch only
    scores the strand it is handed, so it misses every minus-strand copy, and
    it rejects sequences over pyhmmer's 100k-residue limit. nhmmer scans both
    strands in one pass and handles long targets.

    bias_filter is off to match legacy_search (TEsorter's --nobias).

    No model-size partitioning here, unlike legacy_search. That split exists to
    choose between pyhmmer's parallel='queries' and parallel='targets', a switch
    pyhmmer.nhmmer does not expose; nhmmer spreads queries across its threads on
    its own. Hand-rolled target chunking on top of that competes with nhmmer's
    own parallelism rather than adding to it.

    n_workers is passed through as pyhmmer's ``cpus``. It is stated explicitly
    rather than left at pyhmmer's default of 0 (auto-select via
    psutil.cpu_count), so the thread count is the one the caller asked for with
    -p instead of whatever the machine reports -- which on a shared or
    cgroup-limited node is not the same number. Pass 0 to restore auto-select.

    E-values need care. Z does not mean the same thing to the two engines:
    hmmsearch's -Z is a number of comparisons, while nhmmer's -Z is a database
    size in megabases. On top of that, pyhmmer's nhmmer applies Z twice, giving
    E = P * Z^2 where the nhmmer binary gives E = P * megabases.

    Both formulas agree at Z=1, where each returns P, the E-value against a
    1 Mb database. So ask for Z=1 and scale by the real database size here.
    That reproduces the nhmmer binary exactly and keeps working if pyhmmer
    fixes the double-scaling.
    """
    # Caller may pin the search space. A staged cascade must, so that E-values
    # do not improve simply because an upstream stage removed sequences.
    if megabases is None:
        megabases = sum(len(s) for s in seq_block) / 1e6

    # Longest-first, for the same reason legacy_search's partitions are sorted:
    # nhmmer dispatches queries to threads in order, so the expensive models
    # should start first (see _partition_hmms_by_size).
    # domZ is pinned alongside Z. Without it the conditional E-value is
    # normalised against however many domains this particular target set
    # produced, so c_evalue moves when the input is split even though every
    # other E-value holds -- the protein path already pins both (Z=Z, domZ=Z).
    hits = _collect_hits(pyhmmer.nhmmer(
        sorted(hmms, key=lambda h: -h.M), seq_block,
        bias_filter=False,
        Z=1, domZ=1, E=1e10,
        cpus=n_workers,
    ))
    for h in hits:
        for field in _EVALUE_FIELDS:
            h[field] *= megabases
    return [_normalize_nucl_hit(h) for h in hits]


def _partition_hmms_by_size(hmms):
    """
    Partition HMMs into normal and outlier groups by M^2 cost.

    Uses a 2x IQR outlier definition: any model with
    M^2 > Q3 + 2.0 * IQR is an outlier. Outliers benefit from
    parallel='targets' (parallelize over sequences), while normal
    models are faster with parallel='queries' (parallelize over models).

    Returns:
        (normal_hmms, outlier_hmms) -- two lists
    """
    if len(hmms) < 4:
        return hmms, []

    m2_values = sorted(h.M ** 2 for h in hmms)
    q1 = statistics.median(m2_values[:len(m2_values) // 2])
    q3 = statistics.median(m2_values[(len(m2_values) + 1) // 2:])
    iqr = q3 - q1
    threshold = q3 + 2.0 * iqr

    # Longest-first within each group. pyhmmer hands queries to its worker
    # threads in the order given, so starting with the expensive models keeps
    # the short ones as tail filler instead of leaving one long model running
    # alone after everything else has finished (LPT scheduling). Measured on
    # AnnoSINE (87 models, M 95-703) against 3000 targets: 16% faster at
    # cpus=16, 12% at cpus=4, with an identical hit set.
    normal = sorted((h for h in hmms if h.M ** 2 <= threshold),
                    key=lambda h: -h.M)
    outliers = sorted((h for h in hmms if h.M ** 2 > threshold),
                      key=lambda h: -h.M)

    return normal, outliers


def legacy_search(hmms, seq_block, optimized=None):
    """
    Legacy mode: single-pass search with bias filter OFF on all sequences
    against all models. Equivalent to TEsorter's --nobias behavior, just
    faster (hmmsearch instead of hmmscan, pyhmmer instead of subprocess).

    Models are partitioned by M^2 cost: normal models use parallel='queries'
    (faster when models are similarly sized), outliers use parallel='targets'
    (avoids one huge model monopolizing a thread).

    If optimized dict is provided, uses pre-built OptimizedProfiles for
    faster searching (skips per-call profile optimization).

    Z is set to len(hmms) so E-values match hmmscan convention.
    """
    Z = len(hmms)

    # Use optimized profiles if available, otherwise raw HMMs
    def _get_search_objs(hmm_list):
        if optimized:
            return [optimized[h.name] for h in hmm_list if h.name in optimized]
        return hmm_list

    # Single model: always use parallel=targets
    if len(hmms) == 1:
        log.info(f"    1 model (parallel=targets)")
        return _collect_hits(pyhmmer.hmmsearch(
            _get_search_objs(hmms), seq_block,
            bias_filter=False,
            Z=Z, domZ=Z, E=1e10, domE=1e10,
            parallel="targets",
        ))

    normal, outliers = _partition_hmms_by_size(hmms)

    all_hits = []

    if normal:
        log.info(f"    {len(normal)} normal models (parallel=queries)")
        results_iter = pyhmmer.hmmsearch(
            _get_search_objs(normal), seq_block,
            bias_filter=False,
            Z=Z, domZ=Z, E=1e10, domE=1e10,
            parallel="queries",
        )
        all_hits.extend(_collect_hits(results_iter))

    if outliers:
        # pyhmmer returns HMM names as bytes, so this log line has to decode
        # before joining -- it only fires when a database has outlier-sized
        # models, which is why it went unnoticed.
        names = [h.name.decode() if isinstance(h.name, (bytes, bytearray))
                 else str(h.name) for h in outliers]
        log.info(f"    {len(outliers)} outlier models (parallel=targets): "
                 f"{', '.join(n[:40] for n in names[:3])}"
                 f"{'...' if len(names) > 3 else ''}")
        results_iter = pyhmmer.hmmsearch(
            _get_search_objs(outliers), seq_block,
            bias_filter=False,
            Z=Z, domZ=Z, E=1e10, domE=1e10,
            parallel="targets",
        )
        all_hits.extend(_collect_hits(results_iter))

    return all_hits


def pass1_screen(hmms, seq_block, F1=0.02, F2=1e-3, F3=1e-5):
    """
    Pass 1: coarse screen to identify plausible sequence-model pairs.

    Runs hmmsearch with bias filter ON. Filter thresholds F1/F2/F3
    control the MSV/Viterbi/Forward stages respectively. HMMER defaults
    are F1=0.02, F2=1e-3, F3=1e-5. Relaxing F1 (e.g. to 0.1) captures
    more hits in compositionally biased sequences at modest runtime cost.
    Reporting thresholds are maximally permissive -- E-value filtering
    is a downstream concern.

    Z is set to len(hmms) so E-values match hmmscan convention.
    """
    Z = len(hmms)

    normal, outliers = _partition_hmms_by_size(hmms)

    all_hits = []

    if normal:
        all_hits.extend(_collect_hits(pyhmmer.hmmsearch(
            normal, seq_block,
            bias_filter=True, F1=F1, F2=F2, F3=F3,
            Z=Z, domZ=Z, E=1e10, domE=1e10,
            parallel="queries",
        )))

    if outliers:
        all_hits.extend(_collect_hits(pyhmmer.hmmsearch(
            outliers, seq_block,
            bias_filter=True, F1=F1, F2=F2, F3=F3,
            Z=Z, domZ=Z, E=1e10, domE=1e10,
            parallel="targets",
        )))

    seq_models = defaultdict(set)
    for hit in all_hits:
        seq_models[hit["target_name"]].add(hit["query_name"])

    return all_hits, dict(seq_models)


# ---------------------------------------------------------------------------
# Pass 2: M^2-aware cost-balanced work scheduling
# ---------------------------------------------------------------------------

def _invert_seq_models(seq_models):
    """Invert {seq_name: set(model_names)} to {model_name: set(seq_names)}."""
    model_seqs = defaultdict(set)
    for seq_name, models in seq_models.items():
        for model in models:
            model_seqs[model].add(seq_name)
    return dict(model_seqs)


def _plan_work_units(model_seqs, hmms_dict, seq_lens, n_workers):
    """
    Plan cost-balanced work units for pass 2.

    Cost of searching a model against a set of sequences:
        M^2 * total_sequence_length

    Target number of work units is n_workers * 8: enough granularity for
    good load balancing (expensive units finish, workers grab cheap ones),
    without excessive per-unit dispatch overhead.

    Expensive models whose cost exceeds the target unit budget get their
    sequences split into chunks. Cheap models are packed together.
    Units are returned sorted descending by cost so expensive jobs start
    first.

    Returns:
        list of (model_names, seq_names, estimated_cost) tuples,
        sorted descending by cost
    """
    model_costs = {}
    for model_name, seqs in model_seqs.items():
        if model_name not in hmms_dict:
            continue
        m = hmms_dict[model_name].M
        total_seqlen = sum(seq_lens.get(s, 0) for s in seqs)
        model_costs[model_name] = m * m * total_seqlen

    if not model_costs:
        return []

    total_cost = sum(model_costs.values())
    target_n_units = max(n_workers * 8, 1)
    target_cost_per_unit = total_cost / target_n_units

    sorted_models = sorted(model_costs.keys(),
                           key=lambda x: model_costs[x], reverse=True)

    work_units = []  # (model_names, seq_names, cost)

    # Phase 1: expensive models -- split sequences into chunks
    remaining_models = []
    for model_name in sorted_models:
        cost = model_costs[model_name]
        m_sq = hmms_dict[model_name].M ** 2

        if cost > target_cost_per_unit and len(model_seqs[model_name]) > 1:
            # How many chunks based on cost ratio
            n_chunks = max(1, round(cost / target_cost_per_unit))
            n_chunks = min(n_chunks, len(model_seqs[model_name]))

            seqs = list(model_seqs[model_name])
            seq_with_len = sorted(
                [(s, seq_lens.get(s, 0)) for s in seqs],
                key=lambda x: x[1], reverse=True,
            )

            # Greedy bin-packing by sequence length
            chunks = [[] for _ in range(n_chunks)]
            chunk_costs = [0.0] * n_chunks

            for seq_name, slen in seq_with_len:
                min_idx = chunk_costs.index(min(chunk_costs))
                chunks[min_idx].append(seq_name)
                chunk_costs[min_idx] += m_sq * slen

            for chunk, ccost in zip(chunks, chunk_costs):
                if chunk:
                    work_units.append(([model_name], chunk, ccost))
        else:
            remaining_models.append(model_name)

    # Phase 2: pack cheap models together into batch units
    if remaining_models:
        current_models = []
        current_seqs = set()
        current_cost = 0.0

        for model_name in remaining_models:
            m_sq = hmms_dict[model_name].M ** 2
            new_seqs = model_seqs[model_name]
            add_cost = m_sq * sum(seq_lens.get(s, 0) for s in new_seqs)

            current_models.append(model_name)
            current_seqs |= new_seqs
            current_cost += add_cost

            if current_cost >= target_cost_per_unit:
                work_units.append((list(current_models),
                                   list(current_seqs),
                                   current_cost))
                current_models = []
                current_seqs = set()
                current_cost = 0.0

        if current_models:
            work_units.append((list(current_models),
                               list(current_seqs),
                               current_cost))

    # Sort descending by cost -- expensive jobs first
    work_units.sort(key=lambda x: x[2], reverse=True)

    return work_units


# Module-level worker state, initialized per-process
_worker_full_block = None
_worker_name_to_idx = None
_worker_hmms_dict = None
_worker_optimized = None
_worker_alphabet = None


def _init_worker(hmm_path, alphabet_str, seq_fasta):
    """
    Initialize worker: load HMMs, build OptimizedProfiles, and build
    full sequence block + index. Profile optimization happens once here
    instead of on every hmmsearch call.
    """
    global _worker_full_block, _worker_name_to_idx
    global _worker_hmms_dict, _worker_optimized, _worker_alphabet

    if alphabet_str == "amino":
        _worker_alphabet = easel.Alphabet.amino()
    elif alphabet_str == "DNA":
        _worker_alphabet = easel.Alphabet.dna()
    else:
        _worker_alphabet = easel.Alphabet.rna()

    from .hmm import load_hmms, build_optimized_profiles
    all_hmms = load_hmms(hmm_path)
    _worker_hmms_dict = {h.name: h for h in all_hmms}
    _worker_optimized = build_optimized_profiles(all_hmms, _worker_alphabet)

    _worker_full_block = build_sequence_block(seq_fasta, _worker_alphabet)
    _worker_name_to_idx = {
        _worker_full_block[i].name: i
        for i in range(len(_worker_full_block))
    }


def _run_work_unit(args):
    """
    Execute a work unit using pre-initialized worker state.

    Builds a per-unit sequence subset block (sub-millisecond) so that
    large models only scan their assigned chunk of sequences.

    Args:
        tuple of (model_names, seq_names, Z, valid_pairs_list)

    Returns:
        list of hit dicts
    """
    model_names, seq_names, Z, valid_pairs_list = args

    needed_hmms = [_worker_optimized[n] for n in model_names
                   if n in _worker_optimized]

    if not needed_hmms:
        return []

    # Build subset sequence block for just this work unit's sequences
    seq_name_set = set(seq_names)
    subset_indices = [_worker_name_to_idx[n] for n in seq_name_set
                      if n in _worker_name_to_idx]
    subset_seqs = [_worker_full_block[i] for i in subset_indices]

    if not subset_seqs:
        return []

    sub_block = easel.DigitalSequenceBlock(_worker_alphabet, subset_seqs)

    valid_pairs = set(tuple(p) for p in valid_pairs_list)

    results_iter = pyhmmer.hmmsearch(
        needed_hmms, sub_block,
        bias_filter=False,
        Z=Z,
        domZ=Z,
        E=1e10,
        domE=1e10,
    )

    all_hits = _collect_hits(results_iter)

    filtered = []
    for hit in all_hits:
        pair = (hit["target_name"], hit["query_name"])
        if pair in valid_pairs:
            filtered.append(hit)

    return filtered


def pass2_search(hmms_dict, seq_block, seq_models):
    """
    Pass 2 (single-process): sensitive search with bias filter disabled.

    Simple fallback for small workloads or single-processor mode.
    """
    needed_models = set()
    for models in seq_models.values():
        needed_models |= models

    if not needed_models:
        return []

    needed_hmms = [hmms_dict[name] for name in needed_models
                   if name in hmms_dict]

    Z = len(hmms_dict)

    results_iter = pyhmmer.hmmsearch(
        needed_hmms, seq_block,
        bias_filter=False,
        Z=Z,
        domZ=Z,
        E=1e10,
        domE=1e10,
    )

    all_hits = _collect_hits(results_iter)

    filtered = []
    for hit in all_hits:
        seq_name = hit["target_name"]
        model_name = hit["query_name"]
        if seq_name in seq_models and model_name in seq_models[seq_name]:
            filtered.append(hit)

    return filtered


def pass2_search_parallel(hmm_path, seq_fasta, seq_models, hmms_dict,
                          alphabet, n_workers=4):
    """
    Pass 2 (parallel): M^2-aware cost-balanced sensitive search.

    Each worker loads the full HMM database and sequence block at init.
    Per work unit, a subset DigitalSequenceBlock is built (sub-ms cost)
    so large models only scan their assigned sequence chunk.

    Work is planned so that expensive models (large M^2 * total_seqlen)
    get their sequences split across multiple workers, while cheap models
    are packed together into shared units.
    """
    needed_models = set()
    for models in seq_models.values():
        needed_models |= models

    if not needed_models:
        return []

    # Build seq_lens from the FASTA index
    fa = pyfastx.Fasta(seq_fasta, build_index=True)
    seq_lens = {rec.name: len(rec) for rec in fa}

    # Invert: model -> sequences that need it
    model_seqs = _invert_seq_models(seq_models)
    model_seqs = {m: s for m, s in model_seqs.items() if m in hmms_dict}

    # Plan work units
    work_units = _plan_work_units(model_seqs, hmms_dict, seq_lens, n_workers)

    if not work_units:
        return []

    Z = len(hmms_dict)

    # Log the plan
    total_cost = sum(cost for _, _, cost in work_units)
    log.info(f"    Work plan: {len(work_units)} units across {n_workers} workers")
    for i, (models, seqs, cost) in enumerate(work_units[:5]):
        pct = 100 * cost / total_cost if total_cost else 0
        log.info(f"      unit {i}: {len(models)} model(s), "
                 f"{len(seqs)} seqs, {pct:.1f}% of cost")
    if len(work_units) > 5:
        log.info(f"      ... and {len(work_units) - 5} more units")

    # Build valid pairs for each work unit
    work_args = []
    for model_names, seq_names, _ in work_units:
        valid_pairs = []
        model_set = set(model_names)
        for seq_name, models in seq_models.items():
            for m in models:
                if m in model_set:
                    valid_pairs.append((seq_name, m))
        work_args.append((model_names, seq_names, Z, valid_pairs))

    # Alphabet string for serialization to workers
    if alphabet == easel.Alphabet.amino():
        alphabet_str = "amino"
    elif alphabet == easel.Alphabet.dna():
        alphabet_str = "DNA"
    else:
        alphabet_str = "RNA"

    all_hits = []
    with multiprocessing.Pool(
        processes=n_workers,
        initializer=_init_worker,
        initargs=(hmm_path, alphabet_str, seq_fasta),
    ) as pool:
        for result in pool.imap_unordered(_run_work_unit, work_args):
            all_hits.extend(result)

    return all_hits
