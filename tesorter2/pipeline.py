"""
Main pipeline for TE classification.

Orchestrates: FASTA ingestion -> alphabet detection -> the staged search
cascade (hierarchical_search) -> cross-database reconciliation -> BLAST pass-2
-> SQLite + TSV output.

Engine selection lives entirely in --stages. There is no second way in: the
cascade runs for every `sequences` run, and a single-engine run is the stage
list with one name in it.
"""

import argparse
import logging
import os
import sys
import time

from .paths import get_db_dir
from .hmm import peek_alphabet, AMINO_ALPHABET
from .sequence import open_input, looks_nucleotide
from .results import (create_db, store_sequences,
                     index_hits_tables, finalize_db)
from .classifier import (export_classification_tsv, store_classifications,
                       reconcile_classifications)
from .blast_pass2 import blast_pass2
from .subset import subset_fasta
from .hierarchical_search import ENGINES, DEFAULT_STAGES


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
    # AnnoSINE ships plant models; AnnoSINE v2 adds a separate animal set, and
    # RepBase's SINEs are ~68% animal. Kept as its own database rather than
    # merged into AnnoSINE_core.hmm so the plant and animal model sets stay
    # separable, and because a user searching plant genomes should be able to
    # drop it. On by default, like every alias except sine-so.
    "sine-animals": "AnnoSINE_animals.hmm",
    "sine-so":  "SINE_SO.hmm",
    # Pfam-derived TE domains, for families the curated databases do not model:
    # L1 ORF1p/ORF2p flanking domains, DNA-transposon DNA-binding and
    # dimerisation domains, Zisupton/Academ/Dada transposases, the
    # Helitron/Academ helicases. Model names are rewritten into the REXdb
    # convention (see database/kenjis_pfam_TE.manifest.tsv for the Pfam mapping), so no
    # new parser is needed. NOT on by default: Pfam families are defined across
    # a whole fold rather than a TE superfamily, so breadth has to be priced
    # before it can be trusted in a default run.
    "pfam-te":  "kenjis_pfam_TE.hmm",
    # Dfam curated, cut into four collections by measured cost per correct
    # label against RepBase. Every model clears a 95% order-accuracy gate
    # (one-sided Wilson bound), and the collections are disjoint: the aliases
    # below stack, so a model stored in more than one file would be searched
    # more than once.
    "dfam-core":     "Dfam_curated_core.hmm.gz",
    "dfam-extended": "Dfam_curated_extended.hmm.gz",
    "dfam-deep":     "Dfam_curated_deep.hmm.gz",
    "dfam-complete": "Dfam_curated_complete.hmm.gz",
}

# Hierarchical aliases: naming a level implies the cheaper levels beneath it.
# `dfam-deep` is 106 models on its own and means nothing without the 346 in
# front of it, so asking for depth asks for everything up to that depth. Kept
# as an expansion over disjoint files rather than as four cumulative files,
# which would store the core models four times.
DB_EXPANSIONS = {
    "dfam-extended": ("dfam-core", "dfam-extended"),
    "dfam-deep":     ("dfam-core", "dfam-extended", "dfam-deep"),
    "dfam-complete": ("dfam-core", "dfam-extended", "dfam-deep",
                      "dfam-complete"),
}

# Databases that run AFTER the primary arm has classified what it can, against
# only what it left. They are a fallback layer, not a competitor: the protein
# cascade is 97.8% accurate at order level on what it classifies, so letting a
# DNA model outvote it would trade a good call for a worse one. Running them
# second also keeps them out of the cross-database vote entirely.
SUBORDINATE_DBS = ("dfam-core", "dfam-extended", "dfam-deep", "dfam-complete")

# Aliases that exist but are not searched unless named with -d.
#
# pfam-te and dfam-core are now default-ON. dfam-core costs about a third of
# the AnnoSINE search that already ships and recovers ~5,000 RepBase sequences
# the protein cascade misses, at 26 wrong labels; the deeper collections cost
# 3x, 6x and 12x that and are opt-in.
NON_DEFAULT_DBS = ("sine-so", "dfam-extended", "dfam-deep", "dfam-complete")


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


# Subcommands. The mode is the verb: what the input *is* decides almost
# everything downstream -- which engines can run, whether the input is
# windowed, whether there is a per-element classification at all -- so it reads
# better as a mode than as a flag hung off a single command.
MODES = ("sequences", "genome")


def _add_shared(p):
    """Options meaningful in both modes."""
    p.add_argument("sequence", help="Input FASTA")
    p.add_argument(
        "-d", "--database", default=None,
        help="Comma-separated list of database names or paths "
             f"(aliases: {', '.join(DB_ALIASES.keys())}). "
             "[default: every bundled database except sine-so]")
    p.add_argument("--max-search", action="store_true", default=False,
                   help="Search against all known databases. This is the "
                        "default when -d is omitted; the flag still forces the "
                        "full set alongside an explicit -d.")
    p.add_argument("-o", "--outdir", default=None,
                   help="Output directory [default: {input}.TEsorter2]")
    p.add_argument("--prefix", default=None,
                   help="Output file prefix [default: basename of input]")
    p.add_argument("--db-dir", default=None,
                   help="Directory containing the HMM databases "
                        "[default: the databases bundled with the package]")
    p.add_argument("-p", "--processors", type=int, default=4,
                   help="Processors to use [default: 4]")
    p.add_argument("--include-sine-so", action="store_true", default=False,
                   help="Include the SINE_SO model (M=4176) in AnnoSINE "
                        "searches. Excluded by default: it costs 71%% of "
                        "AnnoSINE's compute for <1 filterable hit per 100k "
                        "sequences.")
    p.add_argument("--mask-stops", action="store_true", default=False,
                   help="Six-frame translate stop codons as 'X' (unknown "
                        "residue) rather than '*' (terminator), letting a "
                        "profile align through a premature stop instead of "
                        "being cut short by it. Recovers domains in degraded "
                        "copies, but covers only the in-frame part of what "
                        "BATH does and its false-positive cost is not yet "
                        "characterized.")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="tesorter2",
        description="Fast TE classification using HMM profile databases.")
    sub = parser.add_subparsers(dest="mode", metavar="{sequences,genome}")

    # ---- tesorter2 sequences ----
    seq = sub.add_parser(
        "sequences", help="Classify pre-extracted TE sequences",
        description="Classify pre-extracted TE/LTR sequences. Emits a "
                    "per-element classification per database, a reconciled "
                    "combined call, and a BLAST pass-2 over whatever no "
                    "database classified.")
    _add_shared(seq)
    seq.add_argument(
        "-st", "--seq-type", choices=("nucl", "prot"), default="nucl",
        help="Alphabet of the input sequences. 'prot' takes amino acids "
             "directly: no six-frame translation, no nucleotide-reading engine "
             "(BATH), and DNA profile databases are skipped since nothing can "
             "search them. [default: nucl]")
    seq.add_argument(
        "--stages", default=",".join(DEFAULT_STAGES),
        help="Ordered, comma-separated engine cascade for amino-acid "
             f"databases ({', '.join(sorted(k for k in ENGINES if k != 'nhmmer'))}). "
             "Each stage searches only what earlier stages left unresolved, "
             "per database, so a fast engine strips the easy majority and the "
             "sensitive ones see a shrinking remainder. One name runs that "
             "engine alone -- this is the only way to choose an engine. DNA "
             "databases always use nhmmer, the only engine that reads DNA "
             f"profiles. [default: {','.join(DEFAULT_STAGES)}]")
    seq.add_argument("--compat-tesorter-voting", action="store_true",
                     default=False,
                     help="Decide the clade by raw domain-count plurality, "
                          "replicating TEsorter exactly. Default is a "
                          "score-weighted clade vote.")
    seq.add_argument("--compat-tesorter-rounding", action="store_true",
                     default=False,
                     help="Round normalized scores to 2 dp before threshold "
                          "comparison, replicating a TEsorter rounding bug.")
    seq.add_argument("--compat-tesorter-output", action="store_true",
                     default=False,
                     help="Emit the combined .cls.tsv in TEsorter's 7-column "
                          "format, without SecondaryHits/SO/lineage.")
    seq.add_argument("--no-tesorter-outputs", action="store_true",
                     default=False,
                     help="Skip the TEsorter-compatible companion files "
                          "(.cls.lib, .cls.pep, .dom.gff3, .dom.tsv, .dom.faa).")
    seq.add_argument("-nolib", "--no-library", action="store_true",
                     default=False,
                     help="Do not write the RepeatMasker library (.cls.lib).")
    seq.add_argument("-norc", "--no-reverse", action="store_true",
                     default=False,
                     help="Do not reverse-complement minus-strand sequences "
                          "when writing the .cls.lib library.")
    seq.add_argument(
        "--min-clade-delta", type=float, default=0.0, metavar="BITS",
        help="Demote a clade call to 'mixture' when the vote could be "
             "overturned by a per-hit score offset smaller than BITS. The "
             "clade vote sums score/model_len, so a uniform offset favours "
             "clades carried by short models; this is the margin the call can "
             "absorb. Only the Clade field is affected -- Order and "
             "Superfamily are never touched. 0 (the default) never demotes, "
             "leaving classification "
             "bit-for-bit unchanged. ~5.6 is the measured pyhmmer-to-BATH "
             "offset, so that is the scale at which a call is engine-dependent.")
    seq.add_argument(
        "--emit-clade-delta", action="store_true", default=False,
        help="Append a CladeDelta column to the per-database .cls.tsv with "
             "that margin in bits ('inf' = only one clade in play). "
             "Diagnostic only; does not change any call.")
    # ---- tesorter2 genome ----
    gen = sub.add_parser(
        "genome", help="Annotate TE domains across a whole assembly",
        description="Scan whole-genome sequence for TE domains. Windows the "
                    "genome, classifies each domain occurrence individually, "
                    "and emits a feature-level GFF3 + summary. No per-element "
                    "classification and no BLAST pass-2. Protein profiles are "
                    "searched with BATH (required); DNA profiles with stock "
                    "HMMER nhmmer.")
    _add_shared(gen)
    gen.add_argument("--win-size", type=float, default=1e6,
                     help="Window size for chunking [default: 1e6]")
    gen.add_argument("--win-ovl", type=float, default=1e5,
                     help="Window overlap [default: 1e5]")

    args = parser.parse_args(_normalize_argv(argv))
    if args.mode is None:
        parser.error("a mode is required: %s" % " | ".join(MODES))
    return args


def _normalize_argv(argv):
    """Accept the pre-subcommand form, so existing scripts keep working.

    `tesorter2 in.fa --genome ...` becomes `tesorter2 genome in.fa ...`, and a
    bare `tesorter2 in.fa ...` becomes `tesorter2 sequences in.fa ...`.
    """
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in MODES or argv[0] in ("-h", "--help"):
        return argv
    if argv[0].startswith("-"):
        return argv                      # let argparse report it
    if "--genome" in argv:
        log.warning("`--genome` is deprecated; use `tesorter2 genome ...`")
        return ["genome"] + [a for a in argv if a != "--genome"]
    log.warning("Assuming `tesorter2 sequences ...`; name the mode explicitly.")
    return ["sequences"] + argv


def main():
    args = parse_args()

    genome_mode = args.mode == "genome"

    # Resolve output directory and prefix
    input_base = os.path.basename(args.sequence)
    outdir = args.outdir or f"{input_base}.TEsorter2"
    prefix = args.prefix or input_base
    os.makedirs(outdir, exist_ok=True)

    db_path_out = os.path.join(outdir, f"{prefix}.db")

    # Resolve databases. Omitting -d searches everything: each database
    # classifies independently and reconciliation resolves them afterwards, so
    # more databases is more evidence rather than more ambiguity.
    if args.max_search:
        # Literally everything, which is what the flag says. The default set
        # below is deliberately smaller.
        db_names = list(DB_ALIASES.keys())
    elif not args.database:
        db_names = [k for k in DB_ALIASES.keys() if k not in NON_DEFAULT_DBS]
        if args.include_sine_so:
            db_names.append("sine-so")
    else:
        db_names = [s.strip() for s in args.database.split(",")]

    # Expand hierarchical aliases, preserving order and dropping duplicates so
    # naming both dfam-core and dfam-deep searches core exactly once.
    expanded = []
    for name in db_names:
        for part in DB_EXPANSIONS.get(name, (name,)):
            if part not in expanded:
                expanded.append(part)
    if expanded != db_names:
        added = [n for n in expanded if n not in db_names]
        if added:
            log.info("Hierarchical aliases pulled in: %s", ", ".join(added))
    db_names = expanded

    # A database the user NAMED must exist -- a typo or a missing file should
    # stop the run rather than silently narrow it. One that is merely on by
    # default must not: the Dfam collections are a large optional download, and
    # an installation without them should still run, just without that layer.
    explicit = bool(args.database) and not args.max_search
    db_dir = get_db_dir(args.db_dir)
    db_paths = {}
    absent = []
    for name in db_names:
        try:
            db_paths[name] = resolve_db(name, db_dir=db_dir)
        except FileNotFoundError:
            if explicit:
                raise
            absent.append(name)
    if absent:
        log.warning("Not installed, skipping: %s", ", ".join(absent))
        db_names = [n for n in db_names if n not in absent]
    if not db_names:
        raise SystemExit(
            f"No databases available under {db_dir}. Point --db-dir or "
            f"TESORTER2_DB at a database directory.")

    log.info(f"Input: {args.sequence}")
    log.info(f"Databases: {', '.join(db_names)}")
    log.info(f"Output directory: {outdir}")
    log.info(f"File prefix: {prefix}")

    # Database alphabet decides the engines: AA databases run the --stages
    # cascade, DNA databases nhmmer.
    db_alphabets = {}
    for name, path in db_paths.items():
        alphabet = peek_alphabet(path)
        db_alphabets[name] = alphabet
        log.info(f"  {name}: {alphabet}, translate="
                 f"{'yes' if alphabet == AMINO_ALPHABET else 'no'}")

    # Genome mode: dispatch to the domain-level genome scanner and stop here.
    # It windows the genome, finds/classifies each TE protein domain, and emits
    # a GFF3 + summary -- no element classification, reconcile, or BLAST pass-2.
    if genome_mode:
        from . import genome

        log.info("Genome mode")
        genome.run_genome(
            args.sequence, db_paths, db_alphabets, outdir, prefix,
            n_workers=args.processors,
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

    # A mismatched -st is otherwise silent: protein profiles against
    # nucleotide letters simply find nothing, which reads as a clean run.
    if args.seq_type == "prot" and looks_nucleotide(args.sequence):
        log.warning("-st prot, but the input reads as nucleotide sequence. "
                    "Drop -st prot to six-frame translate it.")
    elif args.seq_type == "nucl" and not looks_nucleotide(args.sequence):
        log.warning("The input does not read as nucleotide sequence. Pass "
                    "-st prot if it is amino acid.")

    # Protein input feeds only the amino-acid databases: nhmmer is the one
    # engine that reads a DNA profile and it needs nucleotide sequence.
    if args.seq_type == "prot":
        dna_dbs = [n for n in db_names if db_alphabets[n] != AMINO_ALPHABET]
        if dna_dbs:
            log.warning("-st prot: skipping DNA profile database(s) %s -- "
                        "nothing can search them with protein input",
                        ", ".join(dna_dbs))
            db_names = [n for n in db_names if n not in dna_dbs]
            db_paths = {n: db_paths[n] for n in db_names}
            db_alphabets = {n: db_alphabets[n] for n in db_names}
        if not db_names:
            raise SystemExit(
                "Nothing to search: every requested database holds DNA "
                "profiles, which protein input cannot feed.")

    # --- Search and classify: the staged cascade ---
    # One cascade per database (--stages picks the engines, DNA databases
    # always nhmmer). Everything after it -- cross-database reconciliation,
    # BLAST pass-2, the combined export -- is shared and engine-agnostic.
    from . import hierarchical_search

    # Two rounds. The primary databases compete with each other and are
    # reconciled by weighted vote; the subordinate ones (Dfam) then run against
    # only what is still unclassified. Splitting the rounds is what makes them
    # subordinate: they never enter the vote, so they cannot outweigh a protein
    # call, and they never see a sequence the primary arm already resolved.
    sub_names = [n for n in db_names if n in SUBORDINATE_DBS]
    primary_names = [n for n in db_names if n not in SUBORDINATE_DBS]
    if sub_names and not primary_names:
        # Nothing to fall through from; run them as the primary arm instead of
        # searching an empty remainder.
        log.info("Only subordinate databases requested; running them directly")
        primary_names, sub_names = sub_names, []
    if sub_names:
        log.info("Primary: %s", ", ".join(primary_names))
        log.info("Subordinate (searched only on what the primary arm leaves): "
                 "%s", ", ".join(sub_names))

    per_db_results = hierarchical_search.run_cascade(
        conn, args.sequence,
        {n: db_paths[n] for n in primary_names},
        {n: db_alphabets[n] for n in primary_names}, outdir,
        protein_stages=[x.strip() for x in args.stages.split(",")],
        n_workers=args.processors,
        compat_rounding=args.compat_tesorter_rounding,
        compat_voting=args.compat_tesorter_voting,
        mask_stops=args.mask_stops,
        min_clade_delta=args.min_clade_delta,
        seq_type=args.seq_type)

    log.info("Indexing hits tables")
    index_hits_tables(conn)

    # Reconcile across databases via hierarchical weighted vote
    reconciled = reconcile_classifications(per_db_results)
    all_classifications = {r["id"]: r for r in reconciled}
    log.info(f"  Reconciled across {len(per_db_results)} databases: "
             f"{len(reconciled)} sequences")

    # --- Subordinate round ---
    # Ahead of BLAST pass-2 on purpose: this is profile evidence against curated
    # models, pass-2 is nucleotide similarity to sequences the run itself just
    # labelled. The weaker evidence should see only what the stronger could not
    # reach.
    if sub_names:
        remaining = [n for n in nucl_lengths if n not in all_classifications]
        log.info("--- Subordinate databases: %s ---", ", ".join(sub_names))
        if not remaining:
            log.info("  Nothing unclassified; skipping")
        else:
            sub_fa = os.path.join(outdir, "subordinate_input.fa")
            n_written, _ = subset_fasta(args.sequence, sub_fa,
                                        set(remaining), exclude=False)
            log.info("  %d sequences unclassified by the primary arm",
                     n_written)
            sub_per_db = hierarchical_search.run_cascade(
                conn, sub_fa,
                {n: db_paths[n] for n in sub_names},
                {n: db_alphabets[n] for n in sub_names}, outdir,
                protein_stages=[x.strip() for x in args.stages.split(",")],
                n_workers=args.processors,
                compat_rounding=args.compat_tesorter_rounding,
                compat_voting=args.compat_tesorter_voting,
                mask_stops=args.mask_stops,
                min_clade_delta=args.min_clade_delta,
                seq_type=args.seq_type)
            index_hits_tables(conn)
            # Reconciled among themselves only -- the collections are disjoint
            # slices of one library, so two of them hitting the same sequence
            # is ordinary agreement, not a cross-database conflict.
            sub_reconciled = reconcile_classifications(sub_per_db)
            gained = 0
            for r in sub_reconciled:
                if r["id"] not in all_classifications:
                    all_classifications[r["id"]] = r
                    reconciled.append(r)
                    gained += 1
            per_db_results.update(sub_per_db)
            log.info("  Recovered %d sequences the primary arm left "
                     "unclassified", gained)

    for name, results in per_db_results.items():
        if not results:
            continue
        cls_tsv = os.path.join(outdir, f"{prefix}.{name}.cls.tsv")
        export_classification_tsv(
            results, cls_tsv, include_clade_delta=args.emit_clade_delta)
        log.info(f"  {name}: {len(results)} classified -> {cls_tsv}")

        # Companion files: the RepeatMasker library only. The domain-level
        # files (.dom.gff3/.dom.tsv/.dom.faa, .cls.pep) encode coordinates on
        # the six-frame translation, and a cascade's hits for one database can
        # come from several engines -- BATH reports nucleotide coordinates,
        # nail its own -- so there is no single frame to write them against.
        if not args.no_tesorter_outputs:
            from .tesorter_output import generate_all_outputs
            generate_all_outputs(
                conn, os.path.join(outdir, f"{prefix}.{name}"), name,
                args.sequence, None, nucl_lengths, results,
                seq_type=args.seq_type,
                no_reverse=args.no_reverse,
                no_library=args.no_library,
                hits_table="legacy_hits",
                domain_files=False,
                keep=None,
            )

    # --- BLAST pass-2 ---
    all_results = list(reconciled)
    if all_classifications:
        log.info("--- BLAST pass-2 ---")
        blast_cls = blast_pass2(
            args.sequence, conn,
            hmm_classifications=all_classifications,
            seq_type=args.seq_type,
            n_processors=args.processors,
            outdir=outdir,
        )

        if blast_cls:
            store_classifications(conn, blast_cls, database="blast_pass2")
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
