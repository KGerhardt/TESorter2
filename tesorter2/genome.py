"""
genome.py — Genome mode: scan whole-genome sequences for TE protein domains.

Mirrors TEsorter's `-genome` mode (TEsorter/app.py: genomeAnn /
_hmm2best(genome=True) / hmm2best(genome=True) / group_resolve_overlaps /
summary_genome), without its windowing. Each chromosome is searched whole for
TE protein-domain models; EACH domain occurrence is classified individually;
overlapping features are resolved; and a domain-level GFF3 + summary table are
emitted.

TEsorter cut the genome into overlapping windows because hmmscan could not
take a chromosome. Neither engine here needs that: BATH reuses nhmmer's
long-target code, splitting each target internally into 256 kb blocks padded
by the longest model's domain envelope, so a hit is never cut at a block edge
and none is reported twice. Both engines read the input as given and report
genomic coordinates directly.

Unlike element mode (classifier.classify_sequences), genome mode:
  * does NOT collapse domains into one classification per input sequence —
    every domain occurrence becomes its own GFF feature,
  * runs NO BLAST pass-2,
  * emits NO per-element .cls.tsv.

Protein-domain search is BATH-only. bathsearch --fs runs directly on the
nucleotide chromosomes and its tblout already yields nucleotide ali coords + strand
per domain occurrence (ascending-normalized, |fwd1/|rev1 suffix), so there is no
six-frame translation and no frame math. BATH's frameshift-aware alignment
recovers a domain broken by an indel as one hit with its true extent, which is
why genome mode uses it exclusively; --bath is implied here and has no effect.
BATH is consequently a hard requirement for genome mode -- see bath_search.
require_binaries, called before any search so a missing binary fails fast.

DNA-alphabet databases (e.g. AnnoSINE) are searched with a second engine:
  * nhmmer, shelled out to stock HMMER (nhmmer_search), not pyHMMER. SINEs are
    non-coding and carry no protein domain, so BATH cannot see them at all --
    bathconvert rejects a DNA profile outright ("Expect Amino Acid"). At genome
    scale the reference nhmmer is what the models were calibrated against and
    what published results compare to, and process startup is negligible against
    a whole assembly; short-sequence mode keeps the in-process pyHMMER path,
    where that tradeoff reverses.
    Like BATH, it reads the raw input and its hits carry genomic coordinates
    plus the |fwd1/|rev1 strand suffix the rest of the pipeline expects.
    A run may mix engines (rexdb via BATH + sine via nhmmer); both report
    nucleotide coordinates, so every feature's sub-sequence is nucleotide.
"""

import itertools
import logging
import os
import re
import sys
import time

from .hmm import AMINO_ALPHABET, DNA_ALPHABET
from .sequence import (open_input, clean_seq, revcomp,
                      parse_frame_suffix, load_sequences_dict)
from .so_map import so_gff_type
from .classifier import (parse_clade_rexdb, parse_clade_gydb, classify_element,
                        load_gydb_clade_map, DB_CONFIGS)
from . import bath_search, nhmmer_search


log = logging.getLogger(__name__)

# Order of TE orders for summary sorting (mirrors TEsorter app.py:56).
ORDERS = ['LTR', 'pararetrovirus', 'DIRS', 'Penelope', 'LINE', 'SINE',
          'TIR', 'Helitron', 'Maverick', 'mixture', 'Unknown', 'Total']

# Hit filters — match classifier.apply_filters defaults / TEsorter genome mode.
MIN_COV = 20.0
MAX_EVALUE = 1e-3
MIN_ACC = 0.5
MIN_NORM_SCORE = 0.1

# ---------------------------------------------------------------------------
# Overlap resolution (ports TEsorter app.py:843-885)
# ---------------------------------------------------------------------------

def _overlap_pct(a, b):
    """Percent overlap relative to the shorter feature. Features are tuples
    with start at index 3 and end at index 4."""
    ovl = max(0, min(a[4], b[4]) - max(a[3], b[3]))
    shorter = min((a[4] - a[3] + 1), (b[4] - b[3] + 1))
    return 100.0 * ovl / shorter if shorter > 0 else 0.0


def _resolve_sort_key(feat):
    """Canonical, run-order-independent ordering for the overlap sweep:
    start asc, end asc, score desc, then gid to break any remaining ties.

    The sweep below is greedy and order-sensitive, so without a total order its
    output depends on the order features were appended (e.g. which databases ran
    and in what sequence). Keying on gid (unique per feature) makes the result
    deterministic and identical regardless of co-searched databases -- so the
    SINE annotation is the same whether or not a protein database ran alongside."""
    return (feat[3], feat[4], -feat[5], feat[9])


def resolve_overlaps(features, max_ovl=20):
    """Within one chromosome (and one feature type), drop equal/overlapping
    features keeping the higher score (feature[5]). Faithful port of TEsorter
    resolve_overlaps, with a deterministic tie-broken sweep order."""
    last = None
    discards = []
    for feat in sorted(features, key=_resolve_sort_key):
        discard = None
        if last:
            if feat == last:
                pair = [last, feat]            # retain, discard (duplicate)
            elif _overlap_pct(feat, last) > max_ovl:
                pair = [feat, last] if feat[5] > last[5] else [last, feat]
            else:
                last = feat
                continue
            _, discard = pair
            discards.append(discard)
        if not last or discard != feat:
            last = feat
    return sorted(set(features) - set(discards), key=lambda x: x[3])


def group_resolve_overlaps(features):
    """Resolve overlaps per chromosome (feature[0]) AND per feature type
    (feature[2]).

    Overlaps are resolved within each feature type, not globally, so a
    protein-domain feature (CDS) and a non-coding-element feature (e.g.
    SINE_element from nhmmer) that share genomic coordinates coexist as
    independent annotation layers rather than evicting one another. This is
    deliberate: their scores live on different scales (nucleotide-HMM bits from
    nhmmer vs protein-HMM bits from HMMER/BATH), so a global 'keep the higher
    score' comparison across the boundary is not meaningful; and the DNA-element
    annotation must not depend on whether a protein database was co-searched.
    Within a single type, resolution is unchanged, so protein-only runs behave
    exactly as before (every feature is CDS -> one group per chromosome)."""
    resolved = []
    keyf = lambda x: (x[0], x[2])
    for (chrom, ftype), items in itertools.groupby(sorted(features, key=keyf),
                                                    key=keyf):
        items = list(items)
        log.info("  resolving overlaps in {} {} ({} features)".format(
            chrom, ftype, len(items)))
        resolved += resolve_overlaps(items)
    return resolved


# ---------------------------------------------------------------------------
# Classification helpers
# ---------------------------------------------------------------------------

def _fmt_cls(*args):
    """Join classification levels, dropping 'unknown' and duplicates
    (faithful port of TEsorter fmt_cls)."""
    values = []
    for arg in args:
        if arg == "unknown" or arg in values:
            continue
        values.append(arg)
    return "/".join(values)


def _parse_model(model, config):
    """Model name -> (gene, clade) for the configured parser."""
    parser = config["clade_parser"]
    if parser == "rexdb":
        _, _, clade, gene = parse_clade_rexdb(model)
        return gene, clade
    if parser == "gydb":
        return parse_clade_gydb(model)
    return "SINE", "SINE"


# ---------------------------------------------------------------------------
# Per-hit -> genomic feature
# ---------------------------------------------------------------------------

def _hit_to_feature(h, engine, config):
    """Turn one passing domain hit into a GFF feature tuple, classified on its
    own. Returns (feature_tuple, classification_string) or None.

    Feature tuple layout (indices used by resolve_overlaps): 0 chrom, 1 source,
    2 type, 3 start, 4 end, 5 score, 6 strand, 7 frame, 8 attributes,
    9 gid, 10 (chrom, start, end, strand) for sequence extraction.
    """
    model = h["query_name"]
    model_len = h["query_len"]
    if model_len <= 0:
        return None
    # Per-DOMAIN score (TEsorter genome uses domscore/tlen). A chromosome
    # holds many domains per model, so the sequence-level score (col 7) would
    # be inflated; dom_score is this domain alone. For BATH and nhmmer,
    # dom_score == score (one alignment per reported hit).
    dom_score = h.get("dom_score")
    if dom_score is None:
        dom_score = h["score"]
    norm_score = dom_score / model_len
    hmm_cov = 100.0 * (h["hmm_to"] - h["hmm_from"] + 1) / model_len
    if not (hmm_cov >= MIN_COV and h["evalue"] <= MAX_EVALUE
            and h["acc"] >= MIN_ACC and norm_score >= MIN_NORM_SCORE):
        return None

    gene, clade = _parse_model(model, config)

    # Both engines search the input chromosomes whole and yield ascending
    # nucleotide ali coords + strand via the |fwd1/|rev1 suffix (BATH from its
    # tblout, nhmmer from search._normalize_nucl_hit), so the coordinates are
    # already genomic.
    chrom, strand, _frame = parse_frame_suffix(h["target_name"])
    g_start, g_end = h["ali_from"], h["ali_to"]   # nucleotide, ascending
    if g_start > g_end:
        g_start, g_end = g_end, g_start
    frame_str = "."
    extract = (chrom, g_start, g_end, strand)

    order, superfamily, disp_clade, complete, _ = classify_element(
        [gene], [clade], [model], [norm_score], config)
    cls = _fmt_cls(order, superfamily, disp_clade)

    gid = "{}:{}-{}|{}".format(chrom, g_start, g_end, model)
    name = "{}-{}".format(clade, gene)
    so_name, so_id = so_gff_type(order, superfamily)
    # GFF type depends on what the hit represents. For BATH the feature is a
    # protein domain, not the element, so its type stays CDS:
    # typing it as the element's SO term would assert the domain *is* the
    # retrotransposon. The element's ontology term goes in Ontology_term instead.
    # An nhmmer hit, by contrast, IS the element (a SINE has no protein domain),
    # so it is typed directly with the element's SO term.
    gff_type = so_name if engine == "nhmmer" else "CDS"
    attr = ("ID={};Name={};Classification={};gene={};clade={};"
            "coverage={:.1f};evalue={:g};score={:.3f};"
            "Ontology_term={};so_name={}").format(
        gid, name, cls, gene, clade, hmm_cov, h["evalue"], norm_score,
        so_id, so_name)

    feature = (chrom, "TEsorter2", gff_type, g_start, g_end,
               round(norm_score, 3), strand, frame_str, attr,
               gid, extract)
    return feature, cls


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def _extract_domain_seq(feature, seq_cache):
    """Nucleotide sub-sequence for one feature (revcomp on minus). The
    extraction tuple is (chrom, start, end, strand), 1-based genomic."""
    key, start, end, strand = feature[10]
    seq = seq_cache.get(key, "")
    sub = seq[start - 1:end]
    if strand == "-":
        sub = revcomp(sub)
    return sub


def write_outputs(features, outdir, prefix, input_fasta):
    """Write {prefix}.dom.gff3, {prefix}.dom.fna, and the summary table.

    Every feature is nucleotide and in genomic coordinates, whichever engine
    found it, so each is read straight out of the input genome.

    Returns (gff_path, [seq_path], summary_path)."""
    gff_path = os.path.join(outdir, "{}.dom.gff3".format(prefix))
    summary_path = os.path.join(outdir, "{}.genome.summary.tsv".format(prefix))
    seq_path = os.path.join(outdir, "{}.dom.fna".format(prefix))

    genome = load_sequences_dict(input_fasta) if features else {}

    with open(gff_path, "w") as fg, open(seq_path, "w") as fs:
        fg.write("##gff-version 3\n")
        for feat in features:
            fg.write("\t".join(map(str, feat[:9])) + "\n")
            sub = _extract_domain_seq(feat, genome)
            fs.write(">{} {}\n{}\n".format(feat[9], feat[8], sub))

    summarize(features, summary_path)
    return gff_path, [seq_path], summary_path


def summarize(features, summary_path):
    """Collapse classified features into a per-classification stats table
    (Order/Superfamily/Clade/Number/Total_length/Mean_length). Ports the
    spirit of TEsorter summary_genome, robust to orders outside ORDERS."""
    d_stats = {}
    for feat in sorted(features, key=lambda x: (x[0], x[3])):
        # Classification is stored in the attribute string.
        m = re.search(r'Classification=([^;]+)', feat[8])
        if not m:
            continue
        cls = tuple(m.group(1).split("/"))
        length = feat[4] - feat[3] + 1
        # expand to Order, Order/SF, Order/SF/Clade, and Total
        xcls = []
        if len(cls) == 3:
            xcls += [cls[:1], cls[:2]]
        elif len(cls) == 2:
            xcls += [cls[:1]]
        xcls += [cls, ("Total",)]
        for c in xcls:
            d_stats.setdefault(c, []).append(length)

    def order_key(item):
        c = item[0]
        try:
            oi = ORDERS.index(c[0])
        except ValueError:
            oi = len(ORDERS)
        return (oi, c)

    lines = [["Order", "Superfamily", "Clade", "Number",
              "Total_length", "Mean_length"]]
    for cls, lens in sorted(d_stats.items(), key=order_key):
        padded = [""] * (len(cls) - 1) + [cls[-1]] + [""] * (3 - len(cls))
        n, total = len(lens), sum(lens)
        lines.append(padded + [str(n), str(total), str(round(total / n, 1))])

    with open(summary_path, "w") as f:
        for line in lines:
            f.write("\t".join(line) + "\n")

    # Echo to stdout (TEsorter prints the summary).
    log.info("Summary of classifications:")
    for line in lines:
        sys.stdout.write("\t".join(line) + "\n")
    sys.stdout.flush()


# ---------------------------------------------------------------------------
# Engine dispatch
# ---------------------------------------------------------------------------

def _search_bath(input_fasta, db_path, db_name, n_workers, outdir):
    """bathsearch --fs directly on the input chromosomes, whole."""
    return bath_search.run_and_parse(
        db_path, input_fasta, db_name, n_workers=n_workers, outdir=outdir)


def _search_nhmmer(db_path, db_name, input_fasta, outdir, n_workers):
    """DNA-DNA search of a nucleotide model set against the RAW input genome,
    both strands in one pass, via stock HMMER's nhmmer binary.

    Genome mode shells out rather than using pyHMMER here. At genome scale the
    reference implementation is what the models were calibrated against and
    what published results are comparable to; pyHMMER's LongTargetsPipeline is
    a reimplementation whose E-value scaling already needed a workaround
    (search.legacy_search_nucl pins Z=1 and rescales by hand because pyHMMER
    applies Z twice). Process startup is irrelevant against a genome, so the
    tradeoff that favours in-process search for short sequences reverses here.

    nhmmer handles long targets, so the target is the whole chromosome and
    hits carry genomic coordinates plus a |fwd1/|rev1 strand suffix."""
    log.info("  Searching DNA models (stock nhmmer, both strands, raw input)")
    hits = nhmmer_search.run_and_parse(
        db_path, input_fasta, db_name,
        os.path.join(outdir, "nhmmer"), n_workers=n_workers)
    log.info("  {} raw nucleotide hits".format(len(hits)))
    return hits


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def run_genome(input_fasta, db_paths, db_alphabets, outdir, prefix,
               n_workers=4):
    """Run genome mode end-to-end. db_paths/db_alphabets keyed by db name.

    Amino-acid databases find TE protein domains with BATH; DNA-alphabet
    databases (e.g. sine) find non-coding elements directly with nhmmer. A
    single run may mix both; both engines map back to the same genomic feature
    space and the outputs are merged."""
    t_start = time.time()

    aa_dbs = [n for n, a in db_alphabets.items() if a == AMINO_ALPHABET]
    nucl_dbs = [n for n, a in db_alphabets.items() if a == DNA_ALPHABET]
    if not aa_dbs and not nucl_dbs:
        raise SystemExit(
            "Genome mode needs at least one searchable database (amino-acid for "
            "TE protein domains, or DNA for non-coding elements like sine).")

    # BATH is the only protein engine, so check for it before searching. A
    # DNA-only run never calls BATH and must not require it.
    if aa_dbs:
        bath_search.require_binaries()
    if nucl_dbs:
        nhmmer_search.require_binaries()

    all_features = []
    for name in aa_dbs + nucl_dbs:
        db_path = db_paths[name]
        config = DB_CONFIGS.get(name)
        if config is None:
            log.warning("  No classifier config for {}, skipping".format(name))
            continue
        if config["clade_parser"] == "gydb":
            config = dict(config)
            config["_clade_map"] = load_gydb_clade_map()

        engine = "nhmmer" if db_alphabets[name] == DNA_ALPHABET else "bath"
        log.info("--- Genome search: {} ({}) [{}] ---".format(
            name, os.path.basename(db_path), engine))
        t0 = time.time()
        if engine == "nhmmer":
            hits = _search_nhmmer(db_path, name, input_fasta, outdir, n_workers)
        else:
            hits = _search_bath(input_fasta, db_path, name, n_workers, outdir)
        log.info("  Search done in {:.1f}s".format(time.time() - t0))

        n_pass = 0
        for h in hits:
            result = _hit_to_feature(h, engine, config)
            if result is not None:
                all_features.append(result[0])
                n_pass += 1
        log.info("  {} features pass filters".format(n_pass))

    log.info("Resolving overlapping features across {} raw features".format(
        len(all_features)))
    resolved = group_resolve_overlaps(all_features)
    log.info("  {} features after overlap resolution".format(len(resolved)))

    # Sort for stable output: chromosome, then start.
    resolved.sort(key=lambda x: (x[0], x[3]))

    gff, seqfs, summ = write_outputs(
        resolved, outdir, prefix, input_fasta)
    log.info("  Feature GFF3: {}".format(gff))
    for p in seqfs:
        log.info("  Feature seqs: {}".format(p))
    log.info("  Summary:      {}".format(summ))
    log.info("Genome mode done in {:.1f}s".format(time.time() - t_start))
    return gff, seqfs, summ
