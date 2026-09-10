"""
hierarchical_search.py — staged multi-engine search cascade.

A cascade runs several search engines in sequence over one database. Each stage
searches only the sequences no earlier stage resolved, so a fast engine strips
the easy majority and the slower, more sensitive engines see a shrinking
remainder:

    stage 0  engine A  -> classified (done) | rejected (done) | remainder v
    stage 1  engine B  -> classified (done) | rejected (done) | remainder v
    stage 2  engine C  -> last recall pass over whatever is left

Two structural decisions:

**The cascade is per database.** Each database runs its own independent
cascade and reaches its own verdict; a sequence resolved by rexdb still runs
the full cascade under gydb. Cross-database deconfliction happens only once
every search is finished, which is where it already happens -- reconciliation
and BLAST pass-2 both treat "unclassified" as unclassified by *any* database
(blast_pass2 reads the per-database classifications table, not the reconciled
set). Running the cascade globally instead would let one database's confident
call suppress another's search and silently starve the vote. The cost is
redundant work on sequences some database already resolved, which is accepted.

**Engine order is data, not code.** Stages are a list; swapping the order or
inserting an engine is a configuration change (`--stages bath,hmmer`). Which
order is best is an empirical question about relative speed and recall, so the
code deliberately does not encode an answer.

Handoff between stages is a FASTA file, because that is the one interface every
engine shares -- BATH and nail are subprocesses that read a path, and pyhmmer
builds its sequence block from one. See subset.subset_fasta.

Pass/fail is decided per *sequence*, never per frame: a translated sequence has
six records but one verdict, matching where classification already draws the
line (results._parse_frame_info keys hits by base name).

Each stage classifies using only its own hits. A sequence reaching stage N was
not classifiable from stages 0..N-1, so there is nothing to carry forward, and
scores from different engines are not on a comparable scale -- pooling them
would produce a normalized score with no consistent meaning.

Known caveat -- BATH results are weakly input-grouping dependent
----------------------------------------------------------------
Two things make an engine's verdict depend on what else was in its input, and
only the first is fixed here:

  1. E-value scales with search-space size, so a subset yields better E-values
     for the same alignment. Fixed: search_space_mb is pinned to the full input
     and passed to every stage (bathsearch -Z, nhmmer's megabase scaling).
     Verified: with -Z pinned, common hits keep identical E-values across
     subsets; without it, none do.

  2. BATH truncates an alignment depending on how much sequence precedes the
     target in the input file. NOT fixed here -- it is a bathsearch bug.

     Minimal case, TIR database, model Ginger-TPase (M=153), one 4477 nt
     target. With <= 20,267 nt of other sequence ahead of it in the file BATH
     reports hmm 44-148, score 19.7. With 25,562 nt ahead it reports hmm
     123-148, score 13.6 -- the same alignment, truncated by 79 model
     positions. Normalized: 0.129 vs 0.089, straddling the classifier's
     min_norm_score of 0.1, so the sequence is classified in one case and not
     the other.

     It is the heuristic filter stage, not the DP: --max (all heuristics off)
     recovers the full alignment at every offset tested. --block_length, purely
     a threading/IO knob, flips it too (50000 recovers 19.7 on the same file
     the 262144 default truncates). BATH is otherwise deterministic -- byte
     identical output across repeat runs and across --cpu 1 vs 8 -- so this is
     reproducible, not a race. The alignment never crosses a sequence boundary;
     coordinates stay inside the target. What moves is how much of the model
     the filter lets through.

     Ruled out: a running RNG accumulator tiebreaking between two alignment
     paths. --seed changes nothing (13.6 across seeds 1, 2, 3, 7, 42, 99 and
     12345 on the same input), so the behaviour is seed-independent as well as
     thread-independent. The remaining candidate is long-target windowing --
     window placement derived from a running offset within the block, so a
     boundary landing inside the domain surfaces only part of it -- but that is
     unconfirmed and belongs upstream rather than here.

The cascade does not cause (2), but it does expose it: different stage orders
hand BATH inputs packed differently, so a cascade's output can differ by stage
order beyond the intended sensitivity difference. Any engine-comparison
analysis must therefore hold the BATH input file byte-identical -- a full run
and a subset run are not comparable for BATH, whatever Z is pinned to. If it
proves material in production the workaround is to run BATH once over the full
input and apply the cascade filter to its results afterwards, which is correct
but forfeits the work reduction that staging exists for.
"""

import logging
import os
import resource
import time

import numpy as np

from .hmm import (load_hmms, build_optimized_profiles, peek_alphabet,
                  AMINO_ALPHABET, DNA_ALPHABET)
from .search import build_sequence_block, legacy_search, legacy_search_nucl
from .sequence import translate_fasta
from .subset import subset_fasta, base_name
from .classifier import classify_sequences, store_classifications, DB_CONFIGS
from .results import store_legacy
from .search import _EVALUE_FIELDS
from . import bath_search, nail_search


log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Hit dicts -> the array layout the classifier consumes
# ---------------------------------------------------------------------------

def hits_to_arrays(hits):
    """Convert engine hit dicts into the structure classify_sequences expects.

    deconflict.load_hits builds this same layout by reading back from SQLite.
    The cascade classifies straight from memory instead: a stage's hits are
    already in hand, and a write-then-read round trip per stage per database
    would dominate the stages that resolve only a handful of sequences. Hits
    are still persisted -- that is provenance, not a data dependency.

    Field choices mirror load_hits exactly so both paths classify identically:
    score is the per-domain score, evalue the independent E-value, and
    model_len the query (model) length.
    """
    if not hits:
        return None

    target = np.array([h["target_name"] for h in hits])
    model = np.array([h["query_name"] for h in hits])
    score = np.array([h["dom_score"] for h in hits], dtype=np.float64)
    evalue = np.array([h["i_evalue"] for h in hits], dtype=np.float64)
    acc = np.array([h["acc"] for h in hits], dtype=np.float64)
    hmm_from = np.array([h["hmm_from"] for h in hits], dtype=np.int32)
    hmm_to = np.array([h["hmm_to"] for h in hits], dtype=np.int32)
    model_len = np.array([h["query_len"] for h in hits], dtype=np.int32)
    env_from = np.array([h["env_from"] for h in hits], dtype=np.int32)
    env_to = np.array([h["env_to"] for h in hits], dtype=np.int32)

    # Same derivations as deconflict._get_base_seq / _get_family.
    base = np.array([t.rsplit("|", 1)[0] if "|" in t else t for t in target])
    family = np.array([
        m.split(":")[1].split("-")[-1] if ":" in m else m.split("_")[0]
        for m in model])

    return {
        "target": target, "model": model, "score": score,
        "evalue": evalue, "acc": acc,
        "hmm_from": hmm_from, "hmm_to": hmm_to, "model_len": model_len,
        "env_from": env_from, "env_to": env_to,
        "base_seq": base, "family": family,
        "hmm_cov": 100.0 * (hmm_to - hmm_from + 1) / model_len,
        "norm_score": score / model_len,
    }


# ---------------------------------------------------------------------------
# Engines
# ---------------------------------------------------------------------------

class Engine:
    """A search tool, reduced to what a cascade needs to know about it.

    name          identifier used in --stages and stamped on every hit
    input_kind    which FASTA the driver materializes for this stage:
                  "aa" (six-frame translations), "aa_nostop" (the same with
                  stop codons rewritten as 'X', which is all nail can read),
                  or "nucl" (the nucleotide input, read directly).
    db_kinds      database alphabets it can search
    """

    name = None
    input_kind = None
    db_kinds = ()

    def search(self, db_path, fasta, db_name, workdir, n_workers,
               search_space_mb=None):
        """Return a list of hit dicts in the shared internal schema.

        search_space_mb pins the E-value search space to the FULL input rather
        than to this stage's subset. Engines whose E-values scale with target
        count (BATH, nhmmer) must honour it or their verdicts drift as upstream
        stages shrink the input; hmmsearch's Z is a model count and is already
        stage-invariant, so pyhmmer ignores it.
        """
        raise NotImplementedError


class _HmmCache:
    """Models are reloaded per (database, engine); parsing REXdb repeatedly is
    pure overhead when a cascade revisits it."""

    def __init__(self):
        self._hmms = {}
        self._optimized = {}

    def hmms(self, db_path):
        if db_path not in self._hmms:
            self._hmms[db_path] = load_hmms(db_path)
        return self._hmms[db_path]

    def optimized(self, db_path, alphabet):
        if db_path not in self._optimized:
            self._optimized[db_path] = build_optimized_profiles(
                self.hmms(db_path), alphabet=alphabet)
        return self._optimized[db_path]


_CACHE = _HmmCache()


class HmmerEngine(Engine):
    """pyHMMER hmmsearch over six-frame translations. The replication
    baseline: this is the path whose output matches TEsorter.

    Named "hmmer" on the command line -- it is HMMER's hmmsearch, run
    in-process through pyHMMER rather than as a subprocess, and the
    distinction is an implementation detail from a user's side.
    """

    name = "hmmer"
    input_kind = "aa"
    db_kinds = ("aa",)

    def search(self, db_path, fasta, db_name, workdir, n_workers,
               search_space_mb=None):
        # search_space_mb is unused: hmmsearch's Z is len(hmms), so E-values
        # do not move when the target set shrinks.
        hmms = _CACHE.hmms(db_path)
        optimized = _CACHE.optimized(db_path, AMINO_ALPHABET)
        block = build_sequence_block(fasta, AMINO_ALPHABET)
        log.info("    pyhmmer: %d models over %d frames", len(hmms), len(block))
        return legacy_search(hmms, block, optimized=optimized)


class NhmmerEngine(Engine):
    """pyHMMER nhmmer for DNA profile databases. Both strands in one pass;
    the only engine that can search a DNA database at all, since BATH and
    nail are protein-profile tools."""

    name = "nhmmer"
    input_kind = "nucl"
    db_kinds = ("dna",)

    def search(self, db_path, fasta, db_name, workdir, n_workers,
               search_space_mb=None):
        hmms = _CACHE.hmms(db_path)
        block = build_sequence_block(fasta, DNA_ALPHABET)
        log.info("    nhmmer: %d models over %d sequences",
                 len(hmms), len(block))
        return legacy_search_nucl(hmms, block, megabases=search_space_mb,
                                  n_workers=n_workers)


class BathEngine(Engine):
    """BATH bathsearch --fs: frameshift-aware translated search run directly
    on nucleotide sequence, so no translation and no coordinate back-mapping.

    report_evalue caps what reaches the tblout. The classifier's own cutoff is
    1e-3, so the default leaves one order of magnitude of sub-threshold
    evidence. A cascade stage that rejects on "no evidence at all" needs a far
    looser cap than that to have anything to judge -- hence the parameter.
    """

    name = "bath"
    input_kind = "nucl"
    db_kinds = ("aa",)

    def __init__(self, report_evalue=0.01, frameshift=True):
        self.report_evalue = report_evalue
        self.frameshift = frameshift

    def search(self, db_path, fasta, db_name, workdir, n_workers,
               search_space_mb=None):
        bath_hmm = bath_search.resolve_bath_db(db_path, db_name=db_name)
        tblout = os.path.join(workdir, "%s.bath.tblout" % db_name)
        bath_search.run_bathsearch(
            bath_hmm, fasta, tblout, n_workers=n_workers,
            frameshift=self.frameshift, report_evalue=self.report_evalue,
            z_megabases=search_space_mb)
        return bath_search.parse_bath_tblout(tblout)


class NailEngine(Engine):
    """nail: MMseqs2-seeded sparse approximation of HMMER's Forward/Backward.

    Reads the HMMER3 .hmm directly, so no conversion step. Amino-acid only.

    Reads the "aa_nostop" source: nail rejects '*' outright, so its targets
    are the shared translation with stops rewritten as 'X'. That rewrite is
    nail-local -- hmmer and bath keep the stop-bearing translation -- so a
    cascade mixing them is deliberately searching two variants of the same
    sequence. See build_stages for what that costs.

    Intended as a cheap first stage. Note that nail exposes no -Z equivalent,
    so search_space_mb cannot be honoured: its E-values scale with whatever
    target set it is handed. Harmless at stage 0, which sees the full input;
    placed later its verdicts would drift with how much upstream stages
    removed, and nothing in its CLI can correct for it. See nail_search.
    """

    name = "nail"
    input_kind = "aa_nostop"
    db_kinds = ("aa",)

    def __init__(self, report_evalue=10.0):
        self.report_evalue = report_evalue

    def search(self, db_path, fasta, db_name, workdir, n_workers,
               search_space_mb=None):
        if search_space_mb is not None:
            log.debug("    nail: search space cannot be pinned (no -Z)")
        model_lengths = {h.name.decode() if isinstance(h.name, bytes)
                         else h.name: h.M
                         for h in _CACHE.hmms(db_path)}
        return nail_search.run_and_parse(
            db_path, fasta, db_name, model_lengths, workdir,
            n_workers=n_workers, report_evalue=self.report_evalue)


ENGINES = {
    "hmmer": HmmerEngine,
    "nhmmer": NhmmerEngine,
    "bath": BathEngine,
    "nail": NailEngine,
}

# "pyhmmer" was the engine's name before the CLI settled on "hmmer"; accepted
# so existing --stages strings keep working.
_ENGINE_ALIASES = {"pyhmmer": "hmmer"}

# Applied to DNA databases regardless of the protein cascade: nhmmer is the
# only engine that can read a DNA profile, so there is nothing to stage.
DNA_STAGES = ("nhmmer",)

# The protein cascade every run gets unless --stages says otherwise: cheapest
# first, most sensitive last. nail strips the easy majority, hmmer is the
# replication baseline, bath is the frameshift-aware last recall pass over
# whatever the first two could not label.
DEFAULT_STAGES = ("nail", "hmmer", "bath")


def build_stages(names):
    """Instantiate engines from a list of names (e.g. ['bath', 'pyhmmer']).

    A cascade containing nail searches two variants of the same sequence, and
    that is accepted rather than refused. nail cannot read stop codons, so its
    targets are always '*' -> 'X' masked; hmmer and bath read the stop-bearing
    translation unless --mask-stops was asked for globally.

    The asymmetry is not cosmetic. 'X' is an unknown residue rather than a
    terminator, so a profile can align straight through a premature stop:
    masking took the TIR fixture from 15 classified to 17. Those two extra
    sequences exit at nail's stage on evidence hmmer would never have seen,
    which is exactly the read-through-degraded-copies behaviour nail is placed
    first to get cheaply -- but it also means a stage-to-stage comparison
    within one cascade is measuring the substitution as well as the engines.
    Anything comparing engine against engine must hold the input identical
    (--mask-stops on, or a stage list without nail).
    """
    stages = []
    for n in names:
        n = _ENGINE_ALIASES.get(n.strip(), n.strip())
        if n not in ENGINES:
            # A configuration error, not a bug: exit cleanly rather than
            # dumping a traceback at someone with a typo in --stages.
            hint = ""
            if n == "facet":
                hint = (" The facet engine was removed; use 'hmmer' for the "
                        "same databases.")
            raise SystemExit("Unknown engine %r; known engines: %s.%s"
                             % (n, ", ".join(sorted(ENGINES)), hint))
        stages.append(ENGINES[n]())

    return stages


def stages_for_input(names, seq_type):
    """Drop the stages this input cannot feed.

    Protein input skips six-frame translation, so there is no nucleotide
    source: BATH reads nucleotide sequence -- it translates internally, with
    frameshifts -- and cannot run at all. Dropping it from the default cascade
    is routine; being left with no stage is a configuration error.
    """
    names = [n.strip() for n in names]
    if seq_type != "prot":
        return names
    kept, dropped = [], []
    for n in names:
        engine = ENGINES.get(_ENGINE_ALIASES.get(n, n))
        target = dropped if (engine is not None
                             and engine.input_kind == "nucl") else kept
        target.append(n)
    if dropped:
        log.warning("Protein input: dropping %s -- it reads nucleotide "
                    "sequence, which -st prot does not provide",
                    ", ".join(dropped))
    if not kept:
        raise SystemExit(
            "No stage left to run: every engine in --stages (%s) reads "
            "nucleotide sequence, which -st prot does not provide. Use nail "
            "and/or hmmer." % ", ".join(dropped))
    return kept


# ---------------------------------------------------------------------------
# Cascade
# ---------------------------------------------------------------------------

CASCADE_EXITS_SCHEMA = """
CREATE TABLE IF NOT EXISTS cascade_exits (
    seq_id    TEXT NOT NULL,
    database  TEXT NOT NULL,
    stage     INTEGER NOT NULL,
    engine    TEXT NOT NULL,
    outcome   TEXT NOT NULL
)
"""


def _record_exits(conn, db_name, stage_idx, engine_name, seq_ids, outcome):
    """Log why each sequence left the cascade.

    This is the record the engine-comparison analyses read: which stage
    resolved a sequence, which rejected it, and which ran out of stages with it
    still unclassified. None of that is recoverable from the hits table alone,
    because "no hit" and "hit that failed a filter" look the same there.
    """
    if not seq_ids:
        return
    conn.execute(CASCADE_EXITS_SCHEMA)
    conn.executemany(
        "INSERT INTO cascade_exits (seq_id, database, stage, engine, outcome) "
        "VALUES (?, ?, ?, ?, ?)",
        [(s, db_name, stage_idx, engine_name, outcome) for s in seq_ids])
    conn.commit()


def _rejected(arrays, remaining, floor):
    """Sequences to drop as showing no usable evidence.

    A sequence is rejected when its best normalized domain score falls below
    `floor` -- including the case of no hits at all. Rejection is the only
    lossy step in the cascade: a rejected sequence never reaches the more
    sensitive engine behind it, so the floor is only safe once measurement
    shows it retains everything a later stage would have classified. It is
    therefore off (floor=None) unless explicitly configured.
    """
    if floor is None:
        return set()
    best = {}
    if arrays is not None:
        for seq, ns in zip(arrays["base_seq"], arrays["norm_score"]):
            if ns > best.get(seq, float("-inf")):
                best[seq] = ns
    return {s for s in remaining if best.get(s, float("-inf")) < floor}


def _cpu_times():
    """(user, sys) seconds for this process AND its reaped children.

    RUSAGE_CHILDREN is the load-bearing half: pyhmmer runs in-process threads,
    but BATH, nail and the MMseqs2 prefilter nail shells out to are all
    subprocesses, and RUSAGE_SELF sees none of their CPU. Reported per stage so
    an in-process engine and a subprocess engine can be compared at all.
    """
    me = resource.getrusage(resource.RUSAGE_SELF)
    kids = resource.getrusage(resource.RUSAGE_CHILDREN)
    return (me.ru_utime + kids.ru_utime, me.ru_stime + kids.ru_stime)


def _count_records(path):
    """Records in a FASTA. For an AA handoff this is frames, not sequences --
    six per sequence -- which is the number the engine actually searches."""
    n = 0
    try:
        with open(path, errors="replace") as fh:
            for line in fh:
                if line[0] == ">":
                    n += 1
    except OSError:
        return None
    return n


def _materialize(kind, remaining, all_names, sources, workdir, tag):
    """FASTA of `remaining` in the alphabet this engine reads.

    Returns the full source file unchanged when nothing has been resolved yet,
    which is the common case for stage 0 and avoids copying the whole input.
    """
    src = sources[kind]
    if len(remaining) == len(all_names):
        return src
    out = os.path.join(workdir, "%s.%s.fa" % (tag, kind))
    n_written, _ = subset_fasta(src, out, remaining, exclude=False)
    log.info("    handoff: %d sequences -> %s", len(remaining),
             os.path.basename(out))
    if n_written == 0:
        return None
    return out


def run_cascade_for_db(conn, db_name, db_path, db_kind, stages, sources,
                       all_names, outdir, n_workers=4, reject_floors=None,
                       compat_rounding=False, compat_voting=False,
                       search_space_mb=None, timing=None, min_clade_delta=0.0,
                       evalue_scale=1.0, seq_index=None):
    """Run one database's cascade. Returns its classification dicts.

    Every stage's hits land in legacy_hits stamped with engine and stage, and
    every sequence's exit is recorded in cascade_exits.
    """
    config = DB_CONFIGS.get(db_name)
    if config is None:
        log.warning("  No classifier config for %s, skipping", db_name)
        return []

    reject_floors = reject_floors or {}
    remaining = set(all_names)
    results = []

    for stage_idx, engine in enumerate(stages):
        if not remaining:
            log.info("  [%s] stage %d (%s): nothing left, stopping",
                     db_name, stage_idx, engine.name)
            break

        workdir = os.path.join(outdir, "cascade", db_name, "s%d" % stage_idx)
        os.makedirs(workdir, exist_ok=True)

        fasta = _materialize(engine.input_kind, remaining, all_names,
                             sources, workdir, db_name)
        if fasta is None:
            break

        n_seqs_in = len(remaining)
        n_recs_in = _count_records(fasta)
        log.info("  [%s] stage %d (%s): %d sequences in"
                 "%s",
                 db_name, stage_idx, engine.name, n_seqs_in,
                 ("" if n_recs_in in (None, n_seqs_in)
                  else " (%d %s records searched)"
                       % (n_recs_in,
                          "frame" if engine.input_kind.startswith("aa")
                          else "nucl")))
        t0 = time.time()
        u0, s0 = _cpu_times()
        hits = engine.search(db_path, fasta, db_name, workdir, n_workers,
                             search_space_mb=search_space_mb)

        # nail exposes no -Z, so its E-values are computed against whatever
        # target set it was handed. Every other engine honours the pinned
        # search space; nail is corrected here instead, which is exact because
        # an E-value scales linearly with search-space size. Under a
        # partitioned run evalue_scale is full_residues/partition_residues;
        # unpartitioned it is 1.0 and this is a no-op.
        if evalue_scale != 1.0 and engine.name == "nail" and hits:
            for h in hits:
                for field in _EVALUE_FIELDS:
                    if field in h:
                        h[field] *= evalue_scale
        wall = time.time() - t0
        u1, s1 = _cpu_times()
        user_s, sys_s = u1 - u0, s1 - s0
        log.info("    %d hits in %.1fs wall  (user %.1fs, sys %.1fs, "
                 "cpu %.1fs, %.0f%% of %d cores)",
                 len(hits), wall, user_s, sys_s, user_s + sys_s,
                 100.0 * (user_s + sys_s) / wall / max(n_workers, 1)
                 if wall > 0 else 0.0, n_workers)

        # A stage may see hits for sequences an earlier stage already resolved
        # when it re-reads a shared source file; keep only the ones still in
        # play so the stage's verdict cannot overwrite a settled one.
        hits = [h for h in hits if base_name(h["target_name"]) in remaining]

        if hits:
            store_legacy(conn, hits, db_name, engine=engine.name,
                         seq_index=seq_index,
                         stage=stage_idx)

        arrays = hits_to_arrays(hits)
        stage_results = []
        if arrays is not None:
            stage_results = classify_sequences(
                arrays, config,
                compat_rounding=compat_rounding,
                compat_voting=compat_voting,
                min_clade_delta=min_clade_delta)
            for r in stage_results:
                r["engine"] = engine.name
                r["stage"] = stage_idx

        labeled = {r["id"] for r in stage_results}
        remaining -= labeled
        _record_exits(conn, db_name, stage_idx, engine.name, labeled,
                      "classified")

        rejected = _rejected(arrays, remaining,
                             reject_floors.get(engine.name))
        if rejected:
            remaining -= rejected
            _record_exits(conn, db_name, stage_idx, engine.name, rejected,
                          "rejected")

        results.extend(stage_results)
        log.info("    %d classified, %d rejected, %d remaining",
                 len(labeled), len(rejected), len(remaining))

        if timing is not None:
            timing.append({
                "database": db_name, "stage": stage_idx,
                "engine": engine.name, "input_kind": engine.input_kind,
                "seqs_in": n_seqs_in,
                "records_searched": n_recs_in if n_recs_in is not None else "",
                "wall_s": round(wall, 2),
                "user_s": round(user_s, 2), "sys_s": round(sys_s, 2),
                "cpu_s": round(user_s + sys_s, 2),
                "hits": len(hits), "classified": len(labeled),
                "rejected": len(rejected), "remaining": len(remaining),
            })

    if remaining:
        _record_exits(conn, db_name, len(stages) - 1,
                      stages[-1].name if stages else "", remaining,
                      "unresolved")

    if results:
        store_classifications(conn, results, database=db_name,
                              mode="hierarchical")
    return results


def run_cascade(conn, input_fasta, db_paths, db_alphabets, outdir,
                protein_stages=DEFAULT_STAGES, n_workers=4,
                reject_floors=None, compat_rounding=False,
                compat_voting=False, aa_fasta=None, mask_stops=False,
                min_clade_delta=0.0, seq_type="nucl",
                search_space_mb=None, evalue_scale=1.0, seq_index=None):
    """Run the cascade for every database. Returns {db_name: [results]}.

    Databases are searched independently and reconciled by the caller -- this
    function is the whole search+classify half of a `sequences` run, and the
    cross-database logic after it is unchanged.

    seq_type="prot" means the input is already amino acid: no translation, no
    nucleotide-reading engine, and no DNA profile database.
    """
    # Pinned once, from the whole input, and reused by every stage of every
    # database. Without this a stage's E-values improve as upstream stages
    # remove sequences, so a cascade's verdicts would depend on how much came
    # before it rather than on the sequence itself.
    lengths = list(_iter_names(input_fasta))
    all_names = [name for name, _ in lengths]
    own_mb = sum(n for _, n in lengths) / 1e6
    if search_space_mb is None:
        search_space_mb = own_mb
        log.info("Search space pinned at %.3f Mb (%d sequences)",
                 search_space_mb, len(all_names))
    else:
        # Partitioned run: this worker holds a slice, but E-values must be
        # those of the whole input or a sequence's verdict would depend on
        # which partition it landed in.
        log.info("Search space pinned at %.3f Mb (whole input); this worker "
                 "holds %.3f Mb (%d sequences)",
                 search_space_mb, own_mb, len(all_names))

    protein_stages = stages_for_input(protein_stages, seq_type)

    # Every source is built once for the whole run and subset per stage, never
    # regenerated. Protein input already is what the AA engines read, so it is
    # its own "aa" source and there is no nucleotide source at all.
    sources = {"nucl": None if seq_type == "prot" else input_fasta,
               "aa": input_fasta if seq_type == "prot" else aa_fasta,
               "aa_nostop": None}
    kinds = {ENGINES[n].input_kind for n in protein_stages if n in ENGINES}
    needs_aa = bool(kinds & {"aa", "aa_nostop"}) and any(
        a == AMINO_ALPHABET for a in db_alphabets.values())
    if needs_aa and not sources["aa"]:
        sources["aa"] = os.path.join(outdir, "cascade_input.aa")
        log.info("Six-frame translating input for AA engines%s",
                 " (stops masked as X)" if mask_stops else "")
        translate_fasta(input_fasta, sources["aa"], mask_stops=mask_stops)

    # nail's stop-free targets, prepared once for the run rather than once per
    # database. The rewrite is a full copy of the translation -- minutes on a
    # large library -- and every nail stage of every database would otherwise
    # redo it. Under --mask-stops the shared translation already says X, so
    # this is the same file and nothing is written.
    if needs_aa and "aa_nostop" in kinds:
        if mask_stops:
            sources["aa_nostop"] = sources["aa"]
        else:
            sources["aa_nostop"], _ = nail_search.prepare_targets(
                sources["aa"],
                os.path.join(outdir, "cascade_input.nostop.aa"))

    per_db = {}
    timing = []
    for db_name, db_path in db_paths.items():
        is_dna = db_alphabets[db_name] == DNA_ALPHABET
        if is_dna and seq_type == "prot":
            log.warning("  Skipping %s: nhmmer is the only engine that reads "
                        "a DNA profile and it needs nucleotide input",
                        db_name)
            continue
        names = DNA_STAGES if is_dna else protein_stages
        stages = build_stages(names)

        # External tools are hard dependencies only if a stage uses them.
        if any(s.name == "bath" for s in stages):
            bath_search.require_binaries()
        if any(s.name == "nail" for s in stages):
            nail_search.require_binaries()

        log.info("--- Cascade: %s [%s] ---", db_name,
                 " -> ".join(s.name for s in stages))
        per_db[db_name] = run_cascade_for_db(
            conn, db_name, db_path, "dna" if is_dna else "aa", stages,
            sources, all_names, outdir, n_workers=n_workers,
            reject_floors=reject_floors, compat_rounding=compat_rounding,
            compat_voting=compat_voting, search_space_mb=search_space_mb,
            timing=timing, min_clade_delta=min_clade_delta,
            evalue_scale=evalue_scale, seq_index=seq_index)

    if timing:
        _write_timing(timing, os.path.join(outdir, "cascade_timing.tsv"))

    return per_db


def _write_timing(rows, path):
    """Per-stage cost table, and a summary on the log.

    A cascade's whole premise is that later stages see fewer sequences, so a
    single total for the run hides the thing being claimed. This records what
    each stage was actually handed and what it cost.
    """
    cols = ["database", "stage", "engine", "input_kind", "seqs_in",
            "records_searched", "wall_s", "user_s", "sys_s", "cpu_s",
            "hits", "classified", "rejected", "remaining"]
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    by_engine = {}
    for r in rows:
        e = by_engine.setdefault(r["engine"], {"wall": 0.0, "cpu": 0.0,
                                               "seqs": 0, "cls": 0})
        e["wall"] += r["wall_s"]; e["cpu"] += r["cpu_s"]
        e["seqs"] += r["seqs_in"]; e["cls"] += r["classified"]
    log.info("Cascade cost by engine (summed over databases):")
    for name, e in sorted(by_engine.items(), key=lambda kv: -kv[1]["wall"]):
        log.info("  %-7s %8.1fs wall  %9.1fs cpu  %8d seqs in  %7d classified",
                 name, e["wall"], e["cpu"], e["seqs"], e["cls"])
    log.info("Per-stage detail -> %s", path)


def _iter_names(fasta):
    """Stream (name, length) without building an index."""
    import pyfastx
    for name, seq in pyfastx.Fasta(fasta, build_index=False):
        yield name, len(seq)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    import argparse

    from .paths import get_db_dir
    from .pipeline import DB_ALIASES, resolve_db
    from .results import create_db, store_sequences, index_hits_tables, finalize_db
    from .classifier import (export_classification_tsv,
                             reconcile_classifications)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s")

    p = argparse.ArgumentParser(
        prog="tesorter2-hierarchical",
        description="Staged multi-engine search cascade: each stage searches "
                    "only what earlier stages left unresolved.")
    p.add_argument("sequence", help="Input FASTA")
    p.add_argument("-d", "--database", default=None,
                   help="Comma-separated database aliases or paths "
                        "[default: every bundled database except sine-so]")
    p.add_argument("-o", "--outdir", default=None,
                   help="Output directory [default: {input}.TEsorter2]")
    p.add_argument("--prefix", default=None, help="Output file prefix")
    p.add_argument("--db-dir", default=None, help="HMM database directory")
    p.add_argument("--stages", default=",".join(DEFAULT_STAGES),
                   help="Ordered protein-engine cascade, comma separated "
                        "(known: %s). DNA databases always use nhmmer, the "
                        "only engine that reads DNA profiles. [default: %s]"
                        % (", ".join(sorted(ENGINES)),
                           ",".join(DEFAULT_STAGES)))
    p.add_argument("-p", "--processors", type=int, default=4)
    p.add_argument(
        "--reject-floor", default=None,
        help="Per-engine normalized-score floor below which a sequence is "
             "dropped instead of passed on, as engine=value pairs "
             "(e.g. pyhmmer=0.05). Off by default: rejection is the only "
             "lossy step, so a floor should be set from measurement showing "
             "it retains what later stages would have classified.")
    p.add_argument(
        "--mask-stops", action="store_true", default=False,
        help="Translate stop codons as 'X' rather than '*', letting profiles "
             "align through premature stops. Off by default; see "
             "sequence.translate_fasta for what it changes.")
    p.add_argument("--compat-tesorter-voting", action="store_true",
                   default=False)
    p.add_argument("--compat-tesorter-rounding", action="store_true",
                   default=False)
    args = p.parse_args(argv)

    input_base = os.path.basename(args.sequence)
    outdir = args.outdir or "%s.TEsorter2" % input_base
    prefix = args.prefix or input_base
    os.makedirs(outdir, exist_ok=True)

    db_dir = get_db_dir(args.db_dir)
    if args.database:
        db_names = [s.strip() for s in args.database.split(",")]
    else:
        db_names = [k for k in DB_ALIASES if k != "sine-so"]
    db_paths = {n: resolve_db(n, db_dir=db_dir) for n in db_names}
    db_alphabets = {n: peek_alphabet(pth) for n, pth in db_paths.items()}

    reject_floors = {}
    if args.reject_floor:
        for item in args.reject_floor.split(","):
            engine, _, value = item.partition("=")
            reject_floors[engine.strip()] = float(value)

    db_out = os.path.join(outdir, "%s.db" % prefix)
    conn = create_db(db_out)
    store_sequences(conn, {n: ln for n, ln in _iter_names(args.sequence)})

    t0 = time.time()
    per_db = run_cascade(
        conn, args.sequence, db_paths, db_alphabets, outdir,
        protein_stages=[s.strip() for s in args.stages.split(",")],
        n_workers=args.processors, reject_floors=reject_floors,
        compat_rounding=args.compat_tesorter_rounding,
        compat_voting=args.compat_tesorter_voting,
        mask_stops=args.mask_stops)

    index_hits_tables(conn)

    for name, results in per_db.items():
        if results:
            tsv = os.path.join(outdir, "%s.%s.cls.tsv" % (prefix, name))
            export_classification_tsv(results, tsv)
            log.info("  %s: %d classified -> %s", name, len(results), tsv)

    reconciled = reconcile_classifications(per_db)
    log.info("Reconciled across %d databases: %d sequences",
             len(per_db), len(reconciled))
    if reconciled:
        combined = os.path.join(outdir, "%s.cls.tsv" % prefix)
        export_classification_tsv(reconciled, combined,
                                  include_secondary=True, include_so=True,
                                  include_lineage=True)
        log.info("  Combined -> %s", combined)

    finalize_db(conn)
    conn.close()
    log.info("Done in %.1fs; results database: %s", time.time() - t0, db_out)
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
