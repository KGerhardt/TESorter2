"""
Main pipeline for TE classification.

Orchestrates: FASTA ingestion -> alphabet detection -> optional translation
-> HMM search -> classification -> pass-2 similarity search -> SQLite + TSV
output.
"""

import argparse
import logging
import os
import time

from .paths import get_db_dir
from .hmm import peek_alphabet, needs_translation, load_hmms, AMINO_ALPHABET, DNA_ALPHABET
from .search import build_sequence_block, legacy_search, legacy_search_nucl
from .sequence import translate_fasta, open_input
from .results import (create_db, store_sequences, store_legacy, store_facet,
                     index_hits_tables, finalize_db,
                     FACET_STAGE_VERIFIED, FACET_STAGE_CROSS_FAMILY,
                     FACET_STAGE_LEGACY_FALLBACK)
from .facet_classify import (facet_classify, facet_classify_v2,
                            export_classifications_tsv)
from .cross_family import find_missing_families, search_missing, search_missing_v2
from .classifier import (classify_sequences, export_classification_tsv,
                       store_classifications, reconcile_classifications,
                       DB_CONFIGS)
from .blast_pass2 import blast_pass2
from .minimap import minimap2_version
from . import bath_search


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)

# Database directory: the databases ship inside the package. Override with
# --db-dir or TESORTER2_DB to point at a custom collection.
DB_DIR = get_db_dir()

DB_ALIASES = {
    "rexdb":    "REXdb_protein_database_viridiplantae_v4.0_plus_metazoa_v3.1.hmm",
    "gydb":     "GyDB2.hmm",
    "line":     "Kapitonov_et_al.GENE.LINE.hmm",
    "tir":      "Yuan_and_Wessler.PNAS.TIR.hmm",
    "sine":     "AnnoSINE_core.hmm",
    "sine-so":  "SINE_SO.hmm",
}


def resolve_db(name, db_dir=None):
    """Resolve a database name or alias to an absolute path."""
    if os.path.isfile(name):
        return os.path.abspath(name)
    if name in DB_ALIASES:
        base = db_dir or DB_DIR
        path = os.path.join(base, DB_ALIASES[name])
        if os.path.isfile(path):
            return os.path.abspath(path)
        raise FileNotFoundError(
            f"Database alias '{name}' -> {path} not found under {base}")
    raise FileNotFoundError(f"Database '{name}' not found (not a file or known alias)")


def parse_args():
    parser = argparse.ArgumentParser(
        prog="tesorter2",
        description="Fast TE classification using HMM profile databases. "
                    "Default mode produces results identical to TEsorter.",
    )
    parser.add_argument(
        "sequence",
        help="Input TE/LTR sequences in FASTA format",
    )
    parser.add_argument(
        "-d", "--database",
        default="rexdb",
        help="Comma-separated list of database names or paths "
             f"(aliases: {', '.join(DB_ALIASES.keys())}) [default: rexdb]",
    )
    parser.add_argument(
        "--max-search",
        action="store_true",
        default=False,
        help="Search against all known databases",
    )
    parser.add_argument(
        "--genome",
        action="store_true",
        default=False,
        help="Genome mode: input is genome sequence(s). Windows the genome, "
             "detects TE protein domains throughout, classifies each domain "
             "individually, and emits a feature-level GFF3 + summary table "
             "(no per-element .cls.tsv, no BLAST pass-2). Works with the "
             "default HMMER engine and with --bath. DNA databases (e.g. sine) "
             "are searched with nhmmer to find non-coding elements too.",
    )
    parser.add_argument(
        "--win-size",
        type=float,
        default=1e6,
        help="Genome mode: window size for chunking [default: 1e6]",
    )
    parser.add_argument(
        "--win-ovl",
        type=float,
        default=1e5,
        help="Genome mode: window overlap [default: 1e5]",
    )
    parser.add_argument(
        "--facet",
        action="store_true",
        default=False,
        help="Facet mode: sub-HMM pre-screen -> verify top hit per family "
             "-> cross-family completion -> legacy fallback. Faster on AA "
             "databases with 99.8%% post-filter recall. DNA databases "
             "automatically use default mode.",
    )
    # Deprecated modes -- retained for backward compatibility, hidden from help
    parser.add_argument("--quick", action="store_true", default=False,
                        help=argparse.SUPPRESS)
    parser.add_argument("--iterative", action="store_true", default=False,
                        help=argparse.SUPPRESS)
    parser.add_argument("--two-pass", action="store_true", default=False,
                        help=argparse.SUPPRESS)
    parser.add_argument("--pass-1-only", action="store_true", default=False,
                        help=argparse.SUPPRESS)
    parser.add_argument(
        "-o", "--outdir",
        default=None,
        help="Output directory [default: {input}.TEsorter2]",
    )
    parser.add_argument(
        "--db-dir",
        default=None,
        help="Directory containing the HMM databases "
             "[default: the databases bundled with the package]",
    )
    parser.add_argument(
        "--dna-engine",
        choices=("nhmmer", "hmmsearch"),
        default="nhmmer",
        help="Engine for DNA databases. nhmmer scans both strands and handles "
             "long targets; hmmsearch only scores the strand it is given and "
             "is limited to 100k-residue sequences. [default: nhmmer]",
    )
    parser.add_argument(
        "--prefix",
        default=None,
        help="Output file prefix [default: basename of input]",
    )
    parser.add_argument(
        "-p", "--processors",
        type=int,
        default=4,
        help="Processors to use [default: 4]",
    )
    parser.add_argument("--F1", type=float, default=0.02,
                        help=argparse.SUPPRESS)  # deprecated, two-pass only
    parser.add_argument(
        "--emit-bath",
        action="store_true",
        default=False,
        help="Emit routed FASTA partitions for the BATH aligner. "
             "Output goes to {outdir}/BATHwater/ directory.",
    )
    parser.add_argument(
        "--bath",
        action="store_true",
        default=False,
        help="Use the BATH aligner (frameshift-aware translated nucleotide "
             "search) instead of pyhmmer/HMMER for amino-acid databases. "
             "Runs bathsearch --fs on the raw nucleotide input (no six-frame "
             "translation). BATH covers protein profiles only, so DNA "
             "databases (AnnoSINE) still go through nhmmer. "
             "Set BATH_BIN_DIR if bathsearch/bathconvert are not on the "
             "default path.",
    )
    parser.add_argument(
        "--compat-tesorter-voting",
        action="store_true",
        default=False,
        help="Decide the clade by raw domain-count plurality, replicating "
             "TEsorter exactly. Default is a score-weighted clade vote (summed "
             "normalized domain score), which resolves sibling-clade swaps that "
             "TEsorter breaks arbitrarily on tied vote counts.",
    )
    parser.add_argument(
        "--include-sine-so",
        action="store_true",
        default=False,
        help="Include the SINE_SO model (M=4176) in AnnoSINE searches. "
             "Excluded by default: SINE_SO costs 71%% of AnnoSINE's compute "
             "but produces <1 filterable hit per 100k sequences.",
    )
    parser.add_argument(
        "--compat-tesorter-rounding",
        action="store_true",
        default=False,
        help="Round normalized scores to 2 decimal places before threshold "
             "comparison, replicating a TEsorter rounding bug. Use only for "
             "exact result reproduction against old TEsorter output.",
    )
    parser.add_argument(
        "--compat-tesorter-output",
        action="store_true",
        default=False,
        help="Emit the combined .cls.tsv in TEsorter's 7-column format. "
             "Default emits an 8th SecondaryHits column listing all "
             "per-database classifications and their summed normalized "
             "scores in descending order of evidence strength.",
    )

    # pass-2 options
    parser.add_argument(
        "-dp2", "--disable-pass2",
        action="store_true", default=False,
        help="Skip pass-2 similarity search (HMM-only classification)",
    )
    parser.add_argument(
        "-rule", "--pass2-rule",
        default="80-80-80", type=str, metavar="I-C-L",
        help="Pass-2 threshold as identity-coverage-length. For the blast "
             "aligner: pident, qcovs, and alignment-length filters (80-80-80 "
             "matches TEsorter2 master). For minimap2: I drives "
             "classify_ltr_paf_fast --min-pid; C drives BOTH --min-qcov and "
             "--min-tcov; L is parsed for grammar compatibility but is not "
             "consumed by classify_ltr_paf_fast [default: %(default)s]",
    )
    parser.add_argument(
        "--pass2-classified-fasta",
        default=None, type=str, metavar="FASTA",
        help="Optional FASTA of previously-classified elements to augment "
             "the pass-2 target database. Headers must be like "
             ">id#Order/Superfamily/Clade",
    )
    parser.add_argument(
        "--minimap2-extra",
        default="", type=str, metavar="STR",
        help="Extra flags passed through to minimap2 (advanced) "
             "[default: empty]",
    )
    parser.add_argument(
        "--pass2-aligner",
        choices=["blast", "minimap2"], default="blast",
        help="Aligner for the pass-2 similarity search. 'blast' (default) "
             "reproduces TEsorter2 master's blastn pass-2 (qcovs + "
             "alignment-length filter, clade=unknown); 'minimap2' uses the "
             "PAF qcov+tcov path and inherits the best target's full "
             "classification. Both share the same -rule and the "
             "--pass2-classified-fasta external-pool merge.",
    )
    parser.add_argument(
        "--blast-task",
        choices=["megablast", "dc-megablast"], default="megablast",
        help="blastn -task for the 'blast' pass-2 aligner. megablast "
             "(default) is fastest and tuned for near-identical matches; "
             "dc-megablast uses discontiguous seeds — slower but more "
             "sensitive to diverged/cross-species matches. Ignored when "
             "--pass2-aligner=minimap2 [default: %(default)s]",
    )

    parser.add_argument(
        "--no-tesorter-outputs",
        action="store_true",
        default=False,
        help="Skip the TEsorter-compatible companion files (.cls.lib, "
             ".cls.pep, .dom.gff3, .dom.tsv, .dom.faa), which are written "
             "alongside each per-database .cls.tsv by default. Writing them "
             "costs extra I/O and holds the input FASTA in memory.",
    )
    parser.add_argument(
        "-nolib", "--no-library",
        action="store_true",
        default=False,
        help="Do not write the RepeatMasker library (.cls.lib).",
    )
    parser.add_argument(
        "-norc", "--no-reverse",
        action="store_true",
        default=False,
        help="Do not reverse-complement minus-strand sequences when writing "
             "the .cls.lib library.",
    )
    return parser.parse_args()


def run_database_legacy(db_path, seq_block, db_name, conn, alphabet=None,
                        facet_fallback=False, dna_engine="nhmmer"):
    """
    Exhaustive single-pass nobias search against all models.

    DNA-alphabet databases go through nhmmer, which scans both strands;
    amino-acid databases go through hmmsearch. Pass dna_engine="hmmsearch" to
    force the old single-strand behaviour on DNA databases.

    When facet_fallback=False (default), this is a true default-mode run and
    hits go to legacy_hits. When True, this call is the DNA-alphabet branch
    of a facet-mode run (facet mode is AA-only); those hits go to facet_hits
    as a legacy-fallback-stage write so the facet and legacy tables stay
    fully disjoint.
    """
    log.info(f"Loading HMMs from {db_name}")
    t0 = time.time()
    hmms = load_hmms(db_path)
    use_nhmmer = alphabet == DNA_ALPHABET and dna_engine == "nhmmer"
    optimized = None
    if not use_nhmmer:
        # nhmmer's pipeline builds its own profiles; optimizing here would be
        # wasted work.
        from .hmm import build_optimized_profiles
        optimized = build_optimized_profiles(hmms, alphabet=alphabet)
    t1 = time.time()
    log.info(f"  Loaded {len(hmms)} models in {t1 - t0:.1f}s")

    t2 = time.time()
    if use_nhmmer:
        log.info(f"  nhmmer search: bias filter OFF, both strands, all models")
        hits = legacy_search_nucl(hmms, seq_block)
    else:
        log.info(f"  Legacy search: bias filter OFF, all models, all sequences")
        hits = legacy_search(hmms, seq_block, optimized=optimized)
    t3 = time.time()
    log.info(f"  {len(hits)} hits in {t3 - t2:.1f}s")

    if facet_fallback:
        store_facet(conn, hits, db_name, stage=FACET_STAGE_LEGACY_FALLBACK)
    else:
        store_legacy(conn, hits, db_name)
    return len(hits)


def run_database(db_path, seq_block, seq_fasta, db_name, alphabet, conn,
                 pass1_only=False, n_workers=4, F1=0.02):
    """
    Run the two-pass search for a single database.

    Args:
        db_path: path to HMM database file
        seq_block: DigitalSequenceBlock (amino or nucl as appropriate)
        seq_fasta: path to the FASTA file (for parallel worker init)
        db_name: short name for this database (for tagging results)
        alphabet: easel.Alphabet for this database
        conn: sqlite3 connection for storing results
        pass1_only: if True, skip pass 2
        n_workers: number of worker processes for pass 2

    Returns:
        tuple of (pass1_hit_count, pass2_hit_count)
    """
    log.info(f"Loading HMMs from {db_name}")
    t0 = time.time()
    hmms = load_hmms(db_path)
    hmms_dict = {hmm.name: hmm for hmm in hmms}
    t1 = time.time()
    log.info(f"  Loaded {len(hmms)} models in {t1 - t0:.1f}s")

    # Pass 1
    log.info(f"  Pass 1: coarse screen (bias filter ON)")
    t2 = time.time()
    p1_hits, seq_models = pass1_screen(hmms, seq_block, F1=F1)
    t3 = time.time()
    n_seqs = len(seq_models)
    n_pairs = sum(len(v) for v in seq_models.values())
    log.info(f"  Pass 1: {len(p1_hits)} hits, {n_seqs} seqs with signal, "
             f"{n_pairs} seq-model pairs in {t3 - t2:.1f}s")

    store_pass1(conn, p1_hits, db_name)

    if pass1_only:
        log.info(f"  --pass-1-only: skipping pass 2 for {db_name}")
        return len(p1_hits), 0

    # Pass 2
    needed_models = set()
    for models in seq_models.values():
        needed_models |= models
    log.info(f"  Pass 2: sensitive search (bias filter OFF) on "
             f"{len(needed_models)} models")

    t4 = time.time()
    if n_workers > 1:
        p2_hits = pass2_search_parallel(
            db_path, seq_fasta, seq_models, hmms_dict,
            alphabet, n_workers=n_workers,
        )
    else:
        p2_hits = pass2_search(hmms_dict, seq_block, seq_models)
    t5 = time.time()
    log.info(f"  Pass 2: {len(p2_hits)} hits in {t5 - t4:.1f}s")

    store_pass2(conn, p2_hits, db_name)

    return len(p1_hits), len(p2_hits)


def main():
    args = parse_args()

    if args.bath and args.facet:
        raise SystemExit(
            "--bath and --facet are incompatible: facet mode is pyhmmer-only. "
            "Run BATH without --facet.")

    if args.genome and args.facet:
        raise SystemExit(
            "--genome and --facet are incompatible: facet mode is element-level. "
            "Run genome mode without --facet.")

    # Resolve output directory and prefix
    input_base = os.path.basename(args.sequence)
    outdir = args.outdir or f"{input_base}.TEsorter2"
    prefix = args.prefix or input_base
    os.makedirs(outdir, exist_ok=True)

    db_path_out = os.path.join(outdir, f"{prefix}.db")
    aa_fasta = os.path.join(outdir, f"{prefix}.aa")

    # Resolve databases
    if args.max_search:
        db_names = [k for k in DB_ALIASES.keys() if k != "sine-so"]
        if args.include_sine_so:
            db_names.append("sine-so")
    else:
        db_names = [s.strip() for s in args.database.split(",")]

    db_dir = get_db_dir(args.db_dir)
    db_paths = {}
    for name in db_names:
        db_paths[name] = resolve_db(name, db_dir=db_dir)

    log.info(f"Input: {args.sequence}")
    log.info(f"Databases: {', '.join(db_names)}")
    log.info(f"Output directory: {outdir}")
    log.info(f"File prefix: {prefix}")

    # Check which alphabets we need
    any_amino = False
    any_nucl = False
    db_alphabets = {}
    for name, path in db_paths.items():
        alphabet = peek_alphabet(path)
        db_alphabets[name] = alphabet
        if alphabet == AMINO_ALPHABET:
            any_amino = True
        else:
            any_nucl = True
        log.info(f"  {name}: {alphabet}, translate={'yes' if alphabet == AMINO_ALPHABET else 'no'}")

    # Genome mode: dispatch to the domain-level genome scanner and stop here.
    # It windows the genome, finds/classifies each TE protein domain, and emits
    # a GFF3 + summary -- no element classification, reconcile, or BLAST pass-2.
    if args.genome:
        from . import genome

        log.info("Genome mode")
        genome.run_genome(
            args.sequence, db_paths, db_alphabets, outdir, prefix,
            use_bath=args.bath, n_workers=args.processors,
            win_size=args.win_size, win_ovl=args.win_ovl)
        return

    # Create results database
    conn = create_db(db_path_out)

    # Read input and store sequence metadata
    t_start = time.time()
    fa = open_input(args.sequence)
    nucl_lengths = {rec.name: len(rec) for rec in fa}
    store_sequences(conn, nucl_lengths)
    log.info(f"Input: {len(nucl_lengths)} sequences")

    # Translate if any database needs amino acid sequences. Under --bath the
    # amino-acid databases are searched by bathsearch directly against the
    # nucleotide input (frameshift-aware translated search), so no six-frame
    # translation is needed for them.
    aa_block = None
    if any_amino and not args.bath:
        # Remove stale index if present
        for f in [aa_fasta + ".fxi"]:
            if os.path.exists(f):
                os.remove(f)
        log.info("Six-frame translating input sequences")
        t0 = time.time()
        translate_fasta(args.sequence, aa_fasta)
        t1 = time.time()
        log.info(f"  Translation done in {t1 - t0:.1f}s -> {aa_fasta}")
        aa_block = build_sequence_block(aa_fasta, AMINO_ALPHABET)
        log.info(f"  Built amino acid sequence block: {len(aa_block)} frames")

    # Build nucleotide block if needed
    nucl_block = None
    if any_nucl:
        log.info("Building nucleotide sequence block")
        nucl_block = build_sequence_block(args.sequence, DNA_ALPHABET)
        log.info(f"  Built nucleotide sequence block: {len(nucl_block)} sequences")

    # Per-DB mode tracks which search path was taken; determines whether the
    # classification step loads from legacy_hits or facet_hits.
    db_modes = {}

    # Run each database
    for name in db_names:
        path = db_paths[name]
        alphabet = db_alphabets[name]

        if alphabet == AMINO_ALPHABET:
            seq_block = aa_block
            seq_fasta = aa_fasta
        else:
            seq_block = nucl_block
            seq_fasta = args.sequence

        log.info(f"--- Searching {name} ({os.path.basename(path)}) ---")

        two_pass = args.two_pass or args.pass_1_only or args.emit_bath

        if args.bath and alphabet == AMINO_ALPHABET:
            t_b0 = time.time()
            hits = bath_search.run_and_parse(
                path, args.sequence, name,
                n_workers=args.processors, outdir=outdir)
            store_legacy(conn, hits, name)
            db_modes[name] = "default"
            log.info(f"  BATH search: {len(hits)} hits in {time.time() - t_b0:.1f}s")
        elif args.facet and alphabet != DNA_ALPHABET:
            t_f0 = time.time()
            classifications, f_verified, f_cross, f_legacy = facet_classify_v2(
                path, seq_block, seq_fasta, alphabet,
                n_workers=args.processors,
                checkpoint_dir=outdir)
            t_f1 = time.time()
            n_primary = sum(1 for c in classifications if not c.get("is_secondary"))
            log.info(f"  Facet mode: {n_primary} assignments, "
                     f"{len(f_verified)} verified, "
                     f"{len(f_cross)} cross-family, "
                     f"{len(f_legacy)} legacy hits in {t_f1 - t_f0:.1f}s")
            if f_verified:
                store_facet(conn, f_verified, name, stage=FACET_STAGE_VERIFIED)
            if f_cross:
                store_facet(conn, f_cross, name, stage=FACET_STAGE_CROSS_FAMILY)
            if f_legacy:
                store_facet(conn, f_legacy, name, stage=FACET_STAGE_LEGACY_FALLBACK)
            db_modes[name] = "facet"
            # Export facet-tier classifications (facet-specific TSV)
            cls_tsv = os.path.join(outdir, f"{prefix}.{name}.classifications.tsv")
            export_classifications_tsv(classifications, cls_tsv)
            log.info(f"  Classifications: {cls_tsv}")
        elif args.facet and alphabet == DNA_ALPHABET:
            log.info(f"  DNA database: using legacy search (facets AA-only)")
            run_database_legacy(path, seq_block, name, conn, alphabet=alphabet,
                                facet_fallback=True, dna_engine=args.dna_engine)
            db_modes[name] = "facet"
        else:
            run_database_legacy(path, seq_block, name, conn, alphabet=alphabet,
                                dna_engine=args.dna_engine)
            db_modes[name] = "default"

    # Build hits-table indexes now that all HMM hits are written and
    # before the classification phase starts reading from them.
    log.info("Indexing hits tables")
    t_idx0 = time.time()
    index_hits_tables(conn)
    log.info(f"  Indexed in {time.time() - t_idx0:.1f}s")

    # --- Classification ---
    log.info("--- Classification ---")
    from .deconflict import load_hits
    per_db_results = {}

    for name in db_names:
        config = DB_CONFIGS.get(name)
        if config is None:
            log.warning(f"  No classifier config for {name}, skipping")
            continue

        mode = db_modes.get(name, "default")
        hits_table = "facet_hits" if mode == "facet" else "legacy_hits"
        hits = load_hits(db_path_out, table=hits_table, database=name)
        if hits is None:
            continue

        log.info(f"  Classifying {name} ({mode} mode, {hits_table})")
        results = classify_sequences(hits, config,
                                     compat_rounding=args.compat_tesorter_rounding,
                                     compat_voting=args.compat_tesorter_voting)

        # Store and export per-database classification (TEsorter format)
        store_classifications(conn, results, database=name, mode=mode)
        cls_tsv = os.path.join(outdir, f"{prefix}.{name}.cls.tsv")
        export_classification_tsv(results, cls_tsv)
        log.info(f"    {len(results)} classified -> {cls_tsv}")

        # TEsorter-compatible companion files, named {prefix}.{db}.* to sit
        # beside the per-database .cls.tsv. The domain-level files encode
        # six-frame translated coordinates, so they are only written when this
        # database was searched on the translated block; --bath and
        # DNA-alphabet databases get the library alone.
        if not args.no_tesorter_outputs and results:
            from .tesorter_output import generate_all_outputs, domain_keys
            from .classifier import select_domain_indices
            six_frame = (aa_block is not None
                         and db_alphabets[name] == AMINO_ALPHABET)
            keep = domain_keys(hits, select_domain_indices(
                hits, config, compat_rounding=args.compat_tesorter_rounding))
            generate_all_outputs(
                conn, os.path.join(outdir, f"{prefix}.{name}"), name,
                args.sequence,
                aa_fasta if six_frame else None,
                nucl_lengths, results,
                seq_type="nucl",
                no_reverse=args.no_reverse,
                no_library=args.no_library,
                hits_table=hits_table,
                domain_files=six_frame,
                keep=keep,
            )

        per_db_results[name] = results

    # Reconcile across databases via hierarchical weighted vote
    reconciled = reconcile_classifications(per_db_results)
    all_classifications = {r["id"]: r for r in reconciled}
    log.info(f"  Reconciled across {len(per_db_results)} databases: "
             f"{len(reconciled)} sequences")

    # --- pass-2 similarity search ---
    all_results = list(reconciled)
    if (not args.pass_1_only and not args.disable_pass2
            and all_classifications):
        try:
            p2_id, p2_cov, p2_len = args.pass2_rule.split("-")
            p2_id = float(p2_id)
            p2_cov = float(p2_cov)
            p2_len = float(p2_len)
        except ValueError:
            raise SystemExit(
                f"--pass2-rule must be I-C-L (three numbers separated by '-'), "
                f"got {args.pass2_rule!r}"
            )

        log.info(f"--- pass-2 ({args.pass2_aligner}) ---")
        if args.pass2_aligner == "minimap2":
            minimap2_version()
        blast_cls = blast_pass2(
            args.sequence, conn,
            hmm_classifications=all_classifications,
            seq_type="nucl",
            n_processors=args.processors,
            min_identity=p2_id,
            min_coverage=p2_cov,
            min_length=p2_len,
            outdir=outdir,
            pass2_classified_fasta=args.pass2_classified_fasta,
            minimap2_extra=args.minimap2_extra,
            aligner=args.pass2_aligner,
            blast_task=args.blast_task,
        )

        if blast_cls:
            run_mode = "facet" if args.facet else "default"
            store_classifications(conn, blast_cls, database="blast_pass2",
                                  mode=run_mode)
            all_results = list(reconciled) + blast_cls

    # Export combined classification
    if all_results:
        combined_tsv = os.path.join(outdir, f"{prefix}.cls.tsv")
        # Lineage rides on the combined file only: the per-database .cls.tsv
        # stays at TEsorter's exact 7 columns, which is what EDTA consumes.
        export_classification_tsv(
            all_results, combined_tsv,
            include_secondary=not args.compat_tesorter_output,
            include_so=not args.compat_tesorter_output,
            include_lineage=not args.compat_tesorter_output,
        )
        log.info(f"  Combined: {len(all_results)} classified -> {combined_tsv}")

    # TODO: rewrite exports with numpy for large datasets
    # Flat file exports temporarily disabled -- results are in the SQLite db
    # # Export flat files
    # log.info("Exporting results")
    #
    # # Determine which table to export from
    # if not two_pass:
    #     hit_table = "legacy_hits"
    # else:
    #     hit_table = "pass2_hits"
    #
    # def outpath(filename):
    #     return os.path.join(outdir, filename)
    #
    # if two_pass:
    #     p1_tsv = outpath(f"{prefix}.pass1.tsv")
    #     export_tsv(conn, p1_tsv, table="pass1_hits")
    #     log.info(f"  Raw pass-1 hits: {p1_tsv}")
    #
    # if not args.pass_1_only:
    #     raw_tsv = outpath(f"{prefix}.{'pass2' if two_pass else 'legacy'}.tsv")
    #     export_tsv(conn, raw_tsv, table=hit_table)
    #     log.info(f"  Raw hits: {raw_tsv}")
    #
    #     best_tsv = outpath(f"{prefix}.best.tsv")
    #     export_best_hits_tsv(conn, best_tsv, nucl_lengths=nucl_lengths,
    #                          table=hit_table)
    #     log.info(f"  Best hits: {best_tsv}")
    #
    #     domains_tsv = outpath(f"{prefix}.domains.tsv")
    #     export_all_domains_tsv(conn, domains_tsv, nucl_lengths=nucl_lengths,
    #                            table=hit_table)
    #     log.info(f"  All domains: {domains_tsv}")
    #
    #     if any_amino:
    #         dom_faa = outpath(f"{prefix}.domains.faa")
    #         export_domain_sequences(conn, dom_faa, aa_fasta,
    #                                 nucl_lengths=nucl_lengths,
    #                                 table=hit_table)
    #         log.info(f"  Domain sequences: {dom_faa}")

    # Build remaining indexes (classifications, blast_hits) and truncate WAL
    log.info("Finalizing database")
    t_fin0 = time.time()
    finalize_db(conn)
    log.info(f"  Finalized in {time.time() - t_fin0:.1f}s")

    t_end = time.time()
    log.info(f"Done in {t_end - t_start:.1f}s")
    log.info(f"Results database: {db_path_out}")

    conn.close()


if __name__ == "__main__":
    main()
