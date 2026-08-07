"""
tesorter_compat.py — TEsorter-compatible CLI entry point.

Accepts the same command-line arguments as TEsorter and produces
output files in the same format and naming convention. Downstream
tools that depend on TEsorter's output structure can use TEsorter2
as a drop-in replacement.

Usage:
    TEsorter input.fasta -db rexdb -p 4
    TEsorter input.fasta -db rexdb -p 4 --facet
"""

import argparse
import logging
import os
import sys

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)

# Map TEsorter -db aliases to our aliases
DB_MAP = {
    "rexdb":         "rexdb",
    "rexdb-plant":   "rexdb",      # v4.0 plant-only (we use v4+metazoa)
    "rexdb-metazoa": "rexdb",      # metazoa-only (subset, use full)
    "gydb":          "gydb",
    "rexdb-pnas":    "tir",
    "rexdb-line":    "line",
    "sine":          "sine",
    # Deprecated v3 aliases -> map to current equivalents
    "rexdb-v3":       "rexdb",
    "rexdb-plantv3":  "rexdb",
    "rexdb-metazoav3":"rexdb",
}


def parse_args():
    parser = argparse.ArgumentParser(
        prog="TEsorter",
        description="TEsorter2 running in TEsorter-compatible mode. "
                    "Lineage-level classification of transposable elements "
                    "using conserved protein domains.",
    )
    parser.add_argument("sequence", type=str,
                        help="Input TE/LTR sequences in FASTA format [required]")
    parser.add_argument("-db", "--hmm-database", type=str, default="rexdb",
                        choices=list(DB_MAP.keys()),
                        help="The database name used [default: rexdb]")
    parser.add_argument("--db-hmm", type=str, default=None,
                        help="Custom HMM database file (overrides -db)")
    parser.add_argument("-st", "--seq-type", type=str, default="nucl",
                        choices=["nucl", "prot"],
                        help="'nucl' for DNA or 'prot' for protein [default: nucl]")
    parser.add_argument("-pre", "--prefix", type=str, default=None,
                        help="Output prefix [default: '{input}.{db}']")
    parser.add_argument("-p", "--processors", type=int, default=4,
                        help="Processors to use [default: 4]")
    parser.add_argument("-tmp", "--tmp-dir", type=str, default=None,
                        help="Directory for temporary files")
    parser.add_argument("-cov", "--min-coverage", type=float, default=20,
                        help="Minimum HMM coverage (0-100) [default: 20]")
    parser.add_argument("-eval", "--max-evalue", type=float, default=1e-3,
                        help="Maximum E-value [default: 1e-3]")
    parser.add_argument("-prob", "--min-probability", type=float, default=0.5,
                        help="Minimum posterior probability [default: 0.5]")
    parser.add_argument("-score", "--min-score", type=float, default=0.1,
                        help="Minimum normalized score [default: 0.1]")
    parser.add_argument("-dp2", "--disable-pass2", action="store_true",
                        default=False,
                        help="Do not run pass-2 classification")
    parser.add_argument("-rule", "--pass2-rule", type=str, default="80-80-80",
                        metavar="I-C-L",
                        help="Pass-2 threshold identity-coverage-length. "
                             "blast: pident, qcovs, and alignment-length "
                             "filters. minimap2: I drives "
                             "classify_ltr_paf_fast --min-pid; C drives both "
                             "--min-qcov and --min-tcov; L unused "
                             "[default: 80-80-80]")
    parser.add_argument("--pass2-aligner", choices=["blast", "minimap2"],
                        default="blast",
                        help="Aligner for the pass-2 similarity search. "
                             "blast (default) reproduces TEsorter2 master's "
                             "blastn pass-2; minimap2 uses the PAF qcov+tcov "
                             "path [default: blast]")
    parser.add_argument("--pass2-classified-fasta", type=str, default=None,
                        metavar="FASTA",
                        help="Optional FASTA of previously-classified elements "
                             "to augment pass-2 target DB. Headers must be "
                             "like >id#Order/Superfamily/Clade")
    parser.add_argument("--minimap2-extra", type=str, default="",
                        metavar="STR",
                        help="Extra flags passed through to minimap2")
    parser.add_argument("-nolib", "--no-library", action="store_true",
                        default=False,
                        help="Do not generate RepeatMasker library file")
    parser.add_argument("-norc", "--no-reverse", action="store_true",
                        default=False,
                        help="Do not reverse complement minus-strand sequences")
    parser.add_argument("-nocln", "--no-cleanup", action="store_true",
                        default=False,
                        help="Do not clean up temporary directory")

    # TEsorter2 extension
    parser.add_argument("--facet", action="store_true", default=False,
                        help="Use facet pre-screening for faster search "
                             "(AA databases only)")

    return parser.parse_args()


def main():
    args = parse_args()

    # Resolve database
    if args.db_hmm:
        db_arg = args.db_hmm
        db_name = os.path.splitext(os.path.basename(args.db_hmm))[0]
    else:
        db_arg = DB_MAP.get(args.hmm_database, args.hmm_database)
        db_name = args.hmm_database

    # Resolve prefix
    if args.prefix is None:
        prefix = f"{os.path.basename(args.sequence)}.{db_name}"
    else:
        prefix = args.prefix

    # Build output directory from prefix
    outdir = os.path.dirname(prefix) or "."
    file_prefix = os.path.basename(prefix)

    log.info(f"TEsorter2 (TEsorter-compatible mode)")
    log.info(f"Input: {args.sequence}")
    log.info(f"Database: {db_name} -> {db_arg}")
    log.info(f"Prefix: {prefix}")

    # Import and run the pipeline
    # We need to set up sys.path if running from the src directory
    src_dir = os.path.dirname(os.path.realpath(__file__))
    if src_dir not in sys.path:
        sys.path.insert(0, src_dir)

    from .hmm import peek_alphabet, load_hmms, AMINO_ALPHABET, DNA_ALPHABET
    from .search import build_sequence_block, legacy_search
    from .sequence import translate_fasta, open_input, clean_seq
    from .results import (create_db, store_sequences, store_legacy, store_facet,
                         index_hits_tables, finalize_db,
                         FACET_STAGE_VERIFIED, FACET_STAGE_CROSS_FAMILY,
                         FACET_STAGE_LEGACY_FALLBACK)
    from .classifier import (classify_sequences, export_classification_tsv,
                            store_classifications, DB_CONFIGS)
    from .blast_pass2 import blast_pass2
    from .deconflict import load_hits

    import time
    t_start = time.time()

    # Resolve database path
    from .pipeline import resolve_db, DB_ALIASES
    try:
        db_path = resolve_db(db_arg)
    except FileNotFoundError:
        if args.db_hmm:
            db_path = args.db_hmm
        else:
            log.error(f"Database not found: {db_arg}")
            sys.exit(1)

    alphabet = peek_alphabet(db_path)

    # Create output database
    db_out = f"{prefix}.db"
    conn = create_db(db_out)

    # Read input
    fa = open_input(args.sequence)
    nucl_lengths = {rec.name: len(rec) for rec in fa}
    store_sequences(conn, nucl_lengths)
    log.info(f"Input: {len(nucl_lengths)} sequences")

    # Translate if needed
    if alphabet == AMINO_ALPHABET and args.seq_type == "nucl":
        aa_fasta = f"{prefix}.aa"
        if os.path.exists(aa_fasta + ".fxi"):
            os.remove(aa_fasta + ".fxi")
        translate_fasta(args.sequence, aa_fasta)
        seq_block = build_sequence_block(aa_fasta, AMINO_ALPHABET)
        seq_fasta = aa_fasta
    elif alphabet == DNA_ALPHABET:
        seq_block = build_sequence_block(args.sequence, DNA_ALPHABET)
        seq_fasta = args.sequence
    else:
        seq_block = build_sequence_block(args.sequence, AMINO_ALPHABET)
        seq_fasta = args.sequence

    # Search
    if args.facet and alphabet == AMINO_ALPHABET:
        from .facet_classify import facet_classify, export_classifications_tsv
        from .cross_family import find_missing_families, search_missing

        classifications, verified_hits, legacy_hits = facet_classify(
            db_path, seq_block, seq_fasta, alphabet,
            n_workers=args.processors)

        if verified_hits:
            store_facet(conn, verified_hits, db_arg, stage=FACET_STAGE_VERIFIED)

        # Cross-family
        hmms = load_hmms(db_path)
        hmms_dict = {h.name: h for h in hmms}
        missing, _ = find_missing_families(classifications, hmms_dict)
        if missing:
            cf_hits = search_missing(missing, hmms_dict, seq_block, alphabet)
            if cf_hits:
                store_facet(conn, cf_hits, db_arg,
                            stage=FACET_STAGE_CROSS_FAMILY)

        if legacy_hits:
            store_facet(conn, legacy_hits, db_arg,
                        stage=FACET_STAGE_LEGACY_FALLBACK)
        run_mode = "facet"
    else:
        hmms = load_hmms(db_path)
        hits = legacy_search(hmms, seq_block)
        store_legacy(conn, hits, db_arg)
        run_mode = "default"

    # Classify
    config = DB_CONFIGS.get(db_arg)
    if config is None:
        # Try to match by name
        for k, v in DB_CONFIGS.items():
            if k in db_arg or db_arg in k:
                config = v
                break

    # Build hits-table indexes before the classification reads
    index_hits_tables(conn)

    if config:
        hits_table = "facet_hits" if run_mode == "facet" else "legacy_hits"
        db_hits = load_hits(db_out, table=hits_table, database=db_arg)
        if db_hits:
            # Drop-in TEsorter mode: keep raw domain-count clade voting so
            # output matches TEsorter exactly. The score-weighted vote is the
            # default only in the main pipeline.py CLI.
            results = classify_sequences(db_hits, config, compat_voting=True)
            store_classifications(conn, results, database=db_arg, mode=run_mode)

            # Export TEsorter-format cls.tsv
            cls_out = f"{prefix}.cls.tsv"
            export_classification_tsv(results, cls_out)
            log.info(f"Classification: {len(results)} sequences -> {cls_out}")

            # pass-2 similarity search
            if not args.disable_pass2 and args.seq_type == "nucl":
                try:
                    p2_id, p2_cov, p2_len = args.pass2_rule.split("-")
                    p2_id = float(p2_id)
                    p2_cov = float(p2_cov)
                    p2_len = float(p2_len)
                except ValueError:
                    log.error(f"--pass2-rule must be I-C-L, got {args.pass2_rule!r}")
                    sys.exit(1)

                hmm_cls = {r["id"]: r for r in results}
                blast_cls = blast_pass2(
                    args.sequence, conn,
                    hmm_classifications=hmm_cls,
                    seq_type="nucl",
                    n_processors=args.processors,
                    min_identity=p2_id,
                    min_coverage=p2_cov,
                    min_length=p2_len,
                    outdir=args.tmp_dir or os.path.dirname(prefix) or ".",
                    pass2_classified_fasta=args.pass2_classified_fasta,
                    minimap2_extra=args.minimap2_extra,
                    aligner=args.pass2_aligner,
                )

                if blast_cls:
                    store_classifications(conn, blast_cls, database="blast_pass2",
                                          mode=run_mode)
                    all_results = results + blast_cls
                    export_classification_tsv(all_results, cls_out)
                    log.info(f"minimap2 pass-2: {len(blast_cls)} additional -> {cls_out}")

    # Generate TEsorter-format output files
    if config and results:
        from .tesorter_output import generate_all_outputs
        all_cls = results + (blast_cls if not args.disable_pass2 and 'blast_cls' in dir() else [])
        generate_all_outputs(
            conn, prefix, db_arg, args.sequence,
            aa_fasta if alphabet == AMINO_ALPHABET and args.seq_type == "nucl" else None,
            nucl_lengths, all_cls,
            seq_type=args.seq_type,
            no_reverse=args.no_reverse,
            no_library=args.no_library,
        )

    finalize_db(conn)

    t_end = time.time()
    log.info(f"Done in {t_end - t_start:.1f}s")
    log.info(f"Results: {db_out}")

    conn.close()


if __name__ == "__main__":
    main()
