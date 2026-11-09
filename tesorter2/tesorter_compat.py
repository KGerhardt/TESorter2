"""
tesorter_compat.py -- TEsorter's own pipeline: the `TEsorter` command.

Same command line as TEsorter (cli.py) and the same output files under the
same names: one database (rexdb by default), HMMER on the six-frame
translation, TEsorter's filters and domain-count clade vote, BLAST pass-2,
and .cls.tsv/.cls.lib/.cls.pep/.dom.gff3/.dom.tsv/.dom.faa. Pipelines that
call TEsorter (EDTA, for one) can use it unchanged. -genome hands over to
TEsorter2's genome mode, which writes TEsorter's genome-mode file names.

Usage:
    TEsorter input.fasta -db rexdb-plant -p 4
"""

import logging
import os
import shutil
import sys

from . import cli

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)


def main(argv=None):
    from .pipeline import DB_ALIASES
    args = cli.parse_args(argv, prog="TEsorter", classic=True,
                          aliases=DB_ALIASES)
    return run_classic(args)


def run_classic(args):
    names, label = cli.database_names(args)
    if cli.ALL_DBS in names or len(names) > 1:
        raise SystemExit("TEsorter's pipeline searches one database; run "
                         "TEsorter2 to search several")
    if not names:
        names, label = ["rexdb"], "rexdb"
    db_arg = names[0]
    db_name = label
    if args.force_write_hmmscan:
        log.info("-fw: TEsorter2 keeps no hmmscan output to reuse; searching")

    prefix = os.path.join(*cli.output_paths(args, label, args.sequence))
    outdir = os.path.dirname(prefix) or "."
    os.makedirs(outdir, exist_ok=True)
    tmp_existed = os.path.isdir(args.tmp_dir)
    os.makedirs(args.tmp_dir, exist_ok=True)
    filters = dict(min_cov=args.min_coverage, max_evalue=args.max_evalue,
                   min_acc=args.min_probability, min_norm_score=args.min_score)

    log.info("TEsorter2 (TEsorter pipeline)")
    log.info(f"Input: {args.sequence}")
    log.info(f"Database: {db_name} -> {db_arg}")
    log.info(f"Prefix: {prefix}")

    from .hmm import peek_alphabet, load_hmms, AMINO_ALPHABET, DNA_ALPHABET
    from .search import build_sequence_block, legacy_search
    from .sequence import translate_fasta, open_input, clean_seq
    from .results import (create_db, store_sequences, store_legacy,
                         finalize_db)
    from .classifier import (classify_sequences, export_classification_tsv,
                            store_classifications, DB_CONFIGS)
    from .blast_pass2 import blast_pass2
    from .deconflict import load_hits

    import time
    t_start = time.time()

    # Resolve database path
    from .pipeline import resolve_db
    from .paths import get_db_dir
    try:
        db_path = resolve_db(db_arg, db_dir=get_db_dir(args.db_dir))
    except FileNotFoundError as e:
        raise SystemExit(str(e))
    if os.path.isfile(db_arg):
        db_arg = db_name               # a custom file goes by its short name

    alphabet = peek_alphabet(db_path)

    if args.genome:
        from . import genome
        if args.win_size is not None or args.win_ovl is not None:
            log.info("-win_size/-win_ovl ignored: TEsorter2 searches each "
                     "chromosome whole")
        if DB_CONFIGS.get(db_arg) is None:
            from .pipeline import _register_custom_db
            _register_custom_db(db_arg)
        gff, _, _ = genome.run_genome(
            args.sequence, {db_arg: db_path}, {db_arg: alphabet}, outdir,
            os.path.basename(prefix), n_workers=args.processors,
            filters=filters, workdir=args.tmp_dir)
        if args.mask:
            cli.mask_input(args.sequence, gff, prefix, args.mask)
        if args.cleanup and not tmp_existed:
            shutil.rmtree(args.tmp_dir, ignore_errors=True)
        return

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

    # No index build before the read below. This mode searches one database,
    # so `WHERE database = ?` matches essentially every row and an index on
    # that column only adds a rowid lookup per row. finalize_db builds them.
    results = []
    blast_cls = []
    if config is None:
        log.warning(f"No classifier config for {db_arg}; nothing classified")
    if config:
        db_hits = load_hits(db_out, table="legacy_hits", database=db_arg)
        if db_hits:
            # Drop-in TEsorter mode: keep raw domain-count clade voting so
            # output matches TEsorter exactly. The score-weighted vote is the
            # default only in the main pipeline.py CLI.
            results = classify_sequences(
                db_hits, config, compat_voting=True,
                filters=filters)
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
                    outdir=args.tmp_dir,
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
        all_cls = results + blast_cls
        generate_all_outputs(
            conn, prefix, db_arg, args.sequence,
            aa_fasta if alphabet == AMINO_ALPHABET and args.seq_type == "nucl" else None,
            nucl_lengths, all_cls,
            seq_type=args.seq_type,
            no_reverse=args.no_reverse,
            no_library=args.no_library,
        )

    if args.mask:
        cli.mask_input(args.sequence, f"{prefix}.dom.gff3", prefix, args.mask,
                       seq_type=args.seq_type)

    finalize_db(conn)

    # -nocln keeps pass-2 scratch; by default remove it, as TEsorter does.
    if args.cleanup:
        shutil.rmtree(os.path.join(args.tmp_dir, "blast_pass2"),
                      ignore_errors=True)
        if not tmp_existed:
            shutil.rmtree(args.tmp_dir, ignore_errors=True)
        for f in (f"{prefix}.aa", f"{prefix}.aa.fxi"):
            if os.path.exists(f):
                os.remove(f)

    t_end = time.time()
    log.info(f"Done in {t_end - t_start:.1f}s")
    log.info(f"Results: {db_out}")

    conn.close()


if __name__ == "__main__":
    main()
