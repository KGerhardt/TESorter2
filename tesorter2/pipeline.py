"""
Main pipeline for TE classification.

Orchestrates: FASTA ingestion -> alphabet detection -> the staged search
cascade (hierarchical_search) -> cross-database reconciliation -> BLAST pass-2
-> SQLite + TSV output.

Engine selection lives entirely in --stages. There is no second way in: the
cascade runs for every element-mode run, and a single-engine run is the stage
list with one name in it.
"""

import logging
import os
import time

from . import cli
from .paths import get_db_dir
from .hmm import peek_alphabet, AMINO_ALPHABET, DNA_ALPHABET
from .sequence import open_input, looks_nucleotide
from .results import create_db, store_sequences, finalize_db
from .classifier import (export_classification_tsv, store_classifications,
                       reconcile_classifications, vote_tier, DB_CONFIGS)
from .blast_pass2 import blast_pass2
from .minimap import minimap2_version
from .subset import subset_fasta
from .hierarchical_search import DEFAULT_STAGES


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
    "sine":     "AnnoSINE_core.hmm.gz",
    # AnnoSINE ships plant models; AnnoSINE v2 adds a separate animal set, and
    # RepBase's SINEs are ~68% animal. Kept as its own database rather than
    # merged into AnnoSINE_core.hmm so the plant and animal model sets stay
    # separable, and because a user searching plant genomes should be able to
    # drop it. On by default.
    "sine-animals": "AnnoSINE_animals.hmm.gz",
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
NON_DEFAULT_DBS = ("dfam-extended", "dfam-deep", "dfam-complete")


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


def parse_args(argv=None):
    """TEsorter's command line plus TEsorter2's options (see cli.py)."""
    return cli.parse_args(argv, prog="TEsorter2", aliases=DB_ALIASES)


def _short_name(path):
    """The name a custom HMM file goes by: its basename, minus .gz/.hmm."""
    name = os.path.basename(path)
    for ext in (".gz", ".hmm"):
        name = name.removesuffix(ext)
    return name


def _register_custom_db(name):
    """Give a custom HMM file a classifier config.

    Model names decide how a hit is read, so the config of the bundled
    database whose name the file's name contains is used (TEsorter's --db-hmm
    behaves the same way), falling back to REXdb's naming convention.
    """
    if name in DB_CONFIGS:
        return
    for key, config in DB_CONFIGS.items():
        if key in name.lower():
            DB_CONFIGS[name] = config
            log.info("Custom database %s: read with the %s config", name, key)
            return
    DB_CONFIGS[name] = DB_CONFIGS["rexdb"]
    log.warning("Custom database %s: no bundled database name in it; reading "
                "its model names with the REXdb convention", name)


def _resolve_stages(stages_arg):
    """The protein-stage list for this run.

    An explicit --stages is taken as given: a missing engine there is a
    configuration error and fails in the cascade's own binary check. The
    default cascade instead drops any stage whose external engine is not
    installed, with a warning, so a plain conda install (which cannot ship
    BATH) still runs. hmmer is pyhmmer and always available.
    """
    if stages_arg is not None:
        return [x.strip() for x in stages_arg.split(",") if x.strip()]

    from . import bath_search, nail_search
    checks = {"bath": bath_search.require_binaries,
              "nail": nail_search.require_binaries}
    kept = []
    for name in DEFAULT_STAGES:
        check = checks.get(name)
        if check is not None:
            try:
                check()
            except SystemExit:
                log.warning("Default cascade: skipping the %s stage, %s is "
                            "not installed (see the README). Pass --stages "
                            "to choose engines explicitly.", name,
                            {"bath": "BATH", "nail": "nail/MMseqs2"}[name])
                continue
        kept.append(name)
    log.info("Protein cascade stages: %s", ",".join(kept))
    return kept


def main(argv=None):
    args = parse_args(argv)

    genome_mode = args.genome
    if args.force_write_hmmscan:
        log.info("-fw: TEsorter2 keeps no hmmscan output to reuse; searching")
    if args.win_size is not None or args.win_ovl is not None:
        log.info("-win_size/-win_ovl ignored: TEsorter2 searches each "
                 "chromosome whole")
    filters = dict(min_cov=args.min_coverage, max_evalue=args.max_evalue,
                   min_acc=args.min_probability,
                   min_norm_score=args.min_score)

    # Parse the pass-2 rule up front so a typo fails before the search runs.
    p2_id = p2_cov = p2_len = None
    if not genome_mode:
        try:
            p2_id, p2_cov, p2_len = (float(x) for x in
                                     args.pass2_rule.split("-"))
        except ValueError:
            raise SystemExit(
                f"--pass2-rule must be I-C-L (three numbers separated by '-'), "
                f"got {args.pass2_rule!r}")

    # Databases. -db takes TEsorter's names as well as TEsorter2's.
    db_names, db_label = cli.database_names(args)
    explicit = bool(db_names) and cli.ALL_DBS not in db_names

    # Output names follow TEsorter: files go beside -pre, which defaults to
    # {input}.{db} (or {input}.TEsorter2 for several databases), in the
    # working directory. An extended run writes into the database it extends,
    # and by default beside it under its name.
    extending = args.extend
    if extending and not db_names:
        raise SystemExit("--extend needs -db: the databases to add")
    if extending and args.prefix is None:
        outdir = args.outdir or os.path.dirname(os.path.abspath(extending))
        prefix = os.path.basename(extending).removesuffix(".db")
    else:
        outdir, prefix = cli.output_paths(args, db_label, args.sequence)
    os.makedirs(outdir, exist_ok=True)
    # Intermediates live under -tmp, one directory per round.
    tmp_dir = args.tmp_dir
    tmp_existed = os.path.isdir(tmp_dir)
    os.makedirs(tmp_dir, exist_ok=True)

    db_path_out = extending or os.path.join(outdir, f"{prefix}.db")

    # Omitting -db searches everything: each database classifies
    # independently and reconciliation resolves them afterwards, so more
    # databases is more evidence rather than more ambiguity.
    if cli.ALL_DBS in db_names:
        # Literally everything, which is what it says. The default set below
        # is deliberately smaller.
        db_names = list(DB_ALIASES.keys()) + [
            n for n in db_names if n != cli.ALL_DBS and n not in DB_ALIASES]
    elif not db_names:
        db_names = [k for k in DB_ALIASES.keys() if k not in NON_DEFAULT_DBS]

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
    # A custom HMM file goes by its short name, so its outputs and its row in
    # the results database read like a bundled one's.
    db_dir = get_db_dir(args.db_dir)
    db_paths = {}
    absent = []
    named = []
    for name in db_names:
        key = name
        if name not in DB_ALIASES and os.path.isfile(name):
            key = _short_name(name)
            _register_custom_db(key)
        try:
            db_paths[key] = resolve_db(name, db_dir=db_dir)
            named.append(key)
        except FileNotFoundError:
            if explicit:
                raise SystemExit(
                    f"Database '{name}' not found: not a file, and not one of "
                    f"{', '.join(dict.fromkeys([*cli.TESORTER_DB_NAMES, *DB_ALIASES]))}"
                    f" under {db_dir}")
            absent.append(name)
    if absent:
        log.warning("Not installed, skipping: %s", ", ".join(absent))
    db_names = named
    if not db_names:
        raise SystemExit(
            f"No databases available under {db_dir}. Point --db-dir or "
            f"TESORTER2_DB at a database directory.")

    log.info(f"Input: {args.sequence}")
    log.info(f"Databases: {', '.join(db_names)}")
    log.info(f"Output prefix: {os.path.join(outdir, prefix)}")
    log.info(f"Temporary directory: {tmp_dir}")

    # Database alphabet decides the engines: AA databases run the --stages
    # cascade, DNA databases nhmmer.
    db_alphabets = {}
    for name, path in db_paths.items():
        alphabet = peek_alphabet(path)
        db_alphabets[name] = alphabet
        log.info(f"  {name}: {alphabet}, translate="
                 f"{'yes' if alphabet == AMINO_ALPHABET else 'no'}")

    # Genome mode: dispatch to the domain-level genome scanner and stop here.
    # It searches each chromosome whole, classifies each TE domain, and emits
    # a GFF3 + summary -- no element classification, reconcile, or BLAST pass-2.
    if genome_mode:
        from . import genome

        log.info("Genome mode")
        gff, _, _ = genome.run_genome(
            args.sequence, db_paths, db_alphabets, outdir, prefix,
            n_workers=args.processors, filters=filters, workdir=tmp_dir)
        if args.mask:
            cli.mask_input(args.sequence, gff, os.path.join(outdir, prefix),
                           args.mask)
        if args.cleanup:
            _remove_tmp(tmp_dir, tmp_existed)
        return

    if args.mask:
        log.warning("-mask needs domain coordinates, which the engine cascade "
                    "does not produce for pre-extracted sequences (its hits "
                    "come from several engines with different coordinate "
                    "systems). Use -genome, or the TEsorter command, for masked "
                    "output.")

    # Create results database. --memory-db keeps it in RAM until the end; see
    # backup_to for why that is not just a speed trick.
    if extending:
        from . import extend
        conn = extend.open_existing(db_path_out, memory=args.memory_db)
        log.info("Extending %s", db_path_out)
    else:
        conn = create_db(":memory:" if args.memory_db else db_path_out)
    if args.memory_db:
        log.info("Results database in memory; written to %s at the end",
                 db_path_out)

    # Read input and store sequence metadata
    t_start = time.time()
    # Wall and CPU per pipeline stage, for stage_timing.tsv (--emit-timing).
    # CPU is this process plus its reaped children: partition workers are
    # reaped when their pool shuts down at the end of a round, and nail, BATH
    # and blastn are subprocesses, so all of it is inside the stage that ran it.
    from .hierarchical_search import _cpu_times
    stage_rows = []

    def _cpu_now():
        return sum(_cpu_times())

    c_start = _cpu_now()

    def _timed(name, cores, fn, *fa, **fk):
        t0, c0 = time.time(), _cpu_now()
        out = fn(*fa, **fk)
        wall, cpu = time.time() - t0, _cpu_now() - c0
        stage_rows.append((name, wall, cpu, cores))
        log.info("  %s: %.1fs wall, %.1fs cpu over %d cores (%.0f%% utilised)",
                 name, wall, cpu, cores,
                 100 * cpu / (wall * cores) if wall else 0)
        return out
    fa = open_input(args.sequence)
    nucl_lengths = {rec.name: len(rec) for rec in fa}
    if extending:
        extend.check_compatible(conn, nucl_lengths)
    else:
        store_sequences(conn, nucl_lengths)
    log.info(f"Input: {len(nucl_lengths)} sequences")

    # An extended run searches only what the earlier run did not. Naming a
    # database it already has is not an error: dfam-extended pulls in
    # dfam-core, which is exactly the one to skip.
    old_per_db, old_is_dna = {}, {}
    if extending:
        done = extend.searched_databases(conn)
        skipped = [n for n in db_names if n in done]
        if skipped:
            log.info("Already in %s, not searched again: %s",
                     os.path.basename(db_path_out), ", ".join(skipped))
        db_names = [n for n in db_names if n not in done]
        db_paths = {n: db_paths[n] for n in db_names}
        db_alphabets = {n: db_alphabets[n] for n in db_names}
        if not db_names:
            raise SystemExit("--extend: every requested database is already "
                             "in %s; nothing to add" % db_path_out)
        log.info("Adding: %s", ", ".join(db_names))
        # Only now is the database touched. The sequences are already stored
        # (check_compatible matched them), so nothing is written for them.
        extend.begin(conn)
        old_per_db, old_is_dna = _timed(
            "reload", 1, extend.reload_results, conn, len(nucl_lengths),
            compat_rounding=args.compat_tesorter_rounding,
            compat_voting=args.compat_tesorter_voting,
            min_clade_delta=args.min_clade_delta, filters=filters)
        log.info("  Rebuilt the calls of %d stored database(s): %s",
                 len(old_per_db), ", ".join(
                     "%s %d" % (n, len(r)) for n, r in old_per_db.items()))

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
    # Stored calls split the same way. A run that had only subordinate
    # databases ran them as its primary arm, so they stay primary here.
    old_sub = {n: r for n, r in old_per_db.items() if n in SUBORDINATE_DBS}
    old_primary = {n: r for n, r in old_per_db.items() if n not in old_sub}
    if old_sub and not old_primary and not primary_names:
        old_primary, old_sub = old_sub, {}
    if sub_names and not primary_names and not old_primary:
        # Nothing to fall through from; run them as the primary arm instead of
        # searching an empty remainder. (Extending a run that had a primary
        # arm, that arm is what they fall through from.)
        log.info("Only subordinate databases requested; running them directly")
        primary_names, sub_names = sub_names, []
    if sub_names:
        log.info("Primary: %s", ", ".join(primary_names))
        log.info("Subordinate (searched only on what the primary arm leaves): "
                 "%s", ", ".join(sub_names))

    # Partitioning is on by default, one group per processor. There is no
    # group-count knob: the groups are balanced on residues to 0.00% and
    # over-decomposing measured worse, not better -- 256 groups over 64 cores
    # made the engines do 20% more total work for the same hits, because each
    # invocation reloads its models. One group per worker is the whole design.
    n_groups = args.processors if args.partition else 0
    if n_groups and n_groups > len(nucl_lengths):
        log.warning("%d processors exceeds the %d input sequences; "
                    "partitioning into %d groups",
                    n_groups, len(nucl_lengths), len(nucl_lengths))
        n_groups = len(nucl_lengths)

    protein_stages = _resolve_stages(args.stages)
    # Convert each protein database to BATH format once, here, before the
    # partition workers start: otherwise every worker converts the same
    # database at the same time (16 workers made 42 conversions of 5).
    if "bath" in protein_stages and args.seq_type != "prot":
        from . import bath_search
        for name, path in db_paths.items():
            if db_alphabets[name] == AMINO_ALPHABET:
                bath_search.resolve_bath_db(path, db_name=name)

    cascade_kwargs = dict(
        protein_stages=protein_stages,
        compat_rounding=args.compat_tesorter_rounding,
        compat_voting=args.compat_tesorter_voting,
        mask_stops=args.mask_stops,
        min_clade_delta=args.min_clade_delta,
        seq_type=args.seq_type,
        filters=filters)

    # A hit's `id` is its sequence's position in the input, so the map is the
    # same for every worker and every round.
    seq_index = {name: i for i, name in enumerate(nucl_lengths)}

    # Every intermediate a round writes lives under tmp_dir/<round>, whichever
    # search path runs it, so cleanup is one directory per round and nothing
    # intermediate is ever written loose beside the results.
    def _workdir(tag):
        path = os.path.join(tmp_dir, tag)
        os.makedirs(path, exist_ok=True)
        return path

    def _search(fasta, names, tag, subset=None):
        """One search round, partitioned or not, over the named databases."""
        paths = {n: db_paths[n] for n in names}
        alphas = {n: db_alphabets[n] for n in names}
        # queue_search.py (the old --partition-queue) is a dead end and is not
        # wired in; see the note at the top of that module.
        if n_groups and n_groups > 1:
            from . import partition
            # Sequence lengths are already in hand from reading the input, so
            # planning needs neither an index query nor an index build.
            lengths = (nucl_lengths if subset is None
                       else {k: nucl_lengths[k] for k in subset})
            return partition.run_partitioned(
                conn, fasta, paths, alphas, _workdir(tag),
                n_groups, args.processors,
                dict(cascade_kwargs, n_workers=1),
                seq_index=seq_index, lengths=lengths, cleanup=args.cleanup,
                timing=args.emit_timing)
        return hierarchical_search.run_cascade(
            conn, fasta, paths, alphas, _workdir(tag),
            n_workers=args.processors, seq_index=seq_index, **cascade_kwargs)

    per_db_results, rounds_run = {}, []
    if primary_names:
        per_db_results = _timed("primary_search", args.processors, _search,
                                args.sequence, primary_names, "primary")
        rounds_run.append("primary")
    # An extended run's stored primary calls vote alongside the new ones, as
    # they would have had every database been searched together.
    per_db_results.update(old_primary)
    for n in old_per_db:
        db_alphabets.setdefault(
            n, DNA_ALPHABET if old_is_dna.get(n) else AMINO_ALPHABET)

    # No index build here. reconcile_classifications works off per_db_results,
    # the objects the cascade just returned, and never reads the database; the
    # legacy_hits indexes are built once at the end by finalize_db.
    # Reconcile across databases via hierarchical weighted vote
    tiers = {n: vote_tier(n, db_alphabets[n]) for n in per_db_results}
    reconciled = reconcile_classifications(per_db_results, tiers=tiers)
    all_classifications = {r["id"]: r for r in reconciled}
    log.info(f"  Reconciled across {len(per_db_results)} databases: "
             f"{len(reconciled)} sequences")

    # --- Subordinate round ---
    # Ahead of BLAST pass-2 on purpose: this is profile evidence against curated
    # models, pass-2 is nucleotide similarity to sequences the run itself just
    # labelled. The weaker evidence should see only what the stronger could not
    # reach.
    # Stored Dfam calls first: they were made on the earlier remainder, and a
    # new collection searches only what they and the primary arm leave.
    if old_sub:
        gained = 0
        for r in reconcile_classifications(old_sub):
            if r["id"] not in all_classifications:
                all_classifications[r["id"]] = r
                reconciled.append(r)
                gained += 1
        per_db_results.update(old_sub)
        log.info("  Stored subordinate calls (%s): %d sequences",
                 ", ".join(old_sub), gained)

    if sub_names:
        remaining = [n for n in nucl_lengths if n not in all_classifications]
        log.info("--- Subordinate databases: %s ---", ", ".join(sub_names))
        if not remaining:
            log.info("  Nothing unclassified; skipping")
        else:
            sub_fa = os.path.join(_workdir("subordinate"), "input.fa")
            n_written, _ = subset_fasta(args.sequence, sub_fa,
                                        set(remaining), exclude=False)
            log.info("  %d sequences unclassified by the primary arm%s",
                     n_written, " and the stored Dfam calls" if old_sub else "")
            sub_per_db = _timed("subordinate_search", args.processors,
                                _search, sub_fa, sub_names, "subordinate",
                                subset=remaining)
            rounds_run.append("subordinate")
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

    # One file per database, TEsorter's format, when there is more than one:
    # with a single database the combined file below is that database's
    # call (plus pass-2), exactly what TEsorter wrote under the same name.
    single_db = len(per_db_results) <= 1
    for name, results in per_db_results.items():
        if not results or single_db:
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
        if not args.no_tesorter_outputs and not args.no_library:
            from .tesorter_output import generate_cls_lib
            generate_cls_lib(args.sequence,
                             os.path.join(outdir, f"{prefix}.{name}"),
                             results, args.no_reverse)

    # --- BLAST pass-2 ---
    all_results = list(reconciled)
    if extending:
        # The earlier pass-2 searched a different split of hit and no-hit
        # sequences; redo it over everything.
        extend.clear_pass2(conn)
    if all_classifications and not args.disable_pass2:
        log.info(f"--- pass-2 ({args.pass2_aligner}) ---")
        if args.pass2_aligner == "minimap2":
            minimap2_version()
        blast_cls = _timed(
            "blast_pass2", args.processors, blast_pass2,
            args.sequence, conn,
            hmm_classifications=all_classifications,
            seq_type=args.seq_type,
            n_processors=args.processors,
            min_identity=p2_id,
            min_coverage=p2_cov,
            min_length=p2_len,
            outdir=tmp_dir,
            pass2_classified_fasta=args.pass2_classified_fasta,
            minimap2_extra=args.minimap2_extra,
            aligner=args.pass2_aligner,
            blast_task=args.blast_task,
            timing=stage_rows,
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
    # TEsorter's RepeatMasker library, over the combined calls and pass-2.
    if not args.no_library:
        from .tesorter_output import generate_cls_lib
        generate_cls_lib(args.sequence, os.path.join(outdir, prefix),
                         all_results, args.no_reverse)

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
    if args.memory_db:
        t_db = time.time()
        from .results import backup_to
        backup_to(conn, db_path_out)
        log.info("  Results database written to disk in %.1fs",
                 time.time() - t_db)
    log.info(f"  Finalized in {time.time() - t_fin0:.1f}s")

    conn.close()

    # Only once the database is closed and on disk: a run that dies before this
    # point keeps every intermediate for diagnosis.
    from .partition import remove_path
    if args.emit_timing:
        _collect_timing(outdir, rounds_run, tmp_dir)
        stage_rows.append(("total", time.time() - t_start,
                           _cpu_now() - c_start, args.processors))
        _write_stage_timing(outdir, stage_rows, n_groups or 0)
    else:
        # A timed run's table, left in a reused outdir, would read as this one's.
        remove_path(os.path.join(outdir, "cascade_timing.tsv"))
        remove_path(os.path.join(outdir, "stage_timing.tsv"))
    if args.cleanup:
        for d in ("primary", "subordinate", "blast_pass2"):
            remove_path(os.path.join(tmp_dir, d))
        _remove_tmp(tmp_dir, tmp_existed)
        log.info("Removed intermediate files (-nocln keeps them)")

    t_end = time.time()
    log.info(f"Done in {t_end - t_start:.1f}s")
    log.info(f"Results database: {db_path_out}")


def _write_stage_timing(outdir, rows, n_groups):
    """stage_timing.tsv: one row per pipeline stage.

    `cores` is what the stage was allotted (-p), so cpu_s / (wall_s * cores)
    is the share of the allotment it used. `partition_groups` is 0 for an
    unpartitioned run. `total` runs from argument parsing to the end, after
    finalize; the BLAST pass-2 sub-steps (blast_makeblastdb, blast_search)
    sit inside blast_pass2, so the rows do not sum.
    """
    with open(os.path.join(outdir, "stage_timing.tsv"), "w") as fh:
        fh.write("stage\twall_s\tcpu_s\tcores\tutilisation\tpartition_groups\n")
        for name, wall, cpu, cores in rows:
            util = cpu / (wall * cores) if wall and cores else 0.0
            fh.write(f"{name}\t{wall:.2f}\t{cpu:.2f}\t{cores}\t{util:.4f}"
                     f"\t{n_groups}\n")


def _remove_tmp(tmp_dir, existed):
    """Remove the temporary directory if this run created it, as TEsorter does.

    A -tmp the user pointed at an existing directory is theirs: only what the
    run wrote in it is removed (by the caller), never the directory.
    """
    if not existed:
        import shutil
        shutil.rmtree(tmp_dir, ignore_errors=True)


def _collect_timing(outdir, rounds, workroot=None):
    """One cascade_timing.tsv at the top of outdir, every round and group.

    Each round writes its own table under outdir/<round>; a partitioned round
    already carries a `group` column, an unpartitioned one gets an empty one so
    every run has the same layout. `round` tells the primary arm from the
    subordinate one, which previously wrote to the same path and overwrote it.
    """
    top = os.path.join(outdir, "cascade_timing.tsv")
    if os.path.exists(top):
        os.remove(top)                 # a previous run's, in a reused outdir
    out_rows, header = [], None
    for rnd in rounds:
        path = os.path.join(workroot or outdir, rnd, "cascade_timing.tsv")
        if not os.path.exists(path):
            continue
        with open(path) as fh:
            cols = fh.readline().rstrip("\n").split("\t")
            grouped = cols[0] == "group"
            if header is None:
                header = ["round"] + (cols if grouped else ["group"] + cols)
            for line in fh:
                if not line.strip():
                    continue
                vals = line.rstrip("\n").split("\t")
                out_rows.append([rnd] + (vals if grouped else [""] + vals))
    if header is None:
        return
    with open(top, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in out_rows:
            fh.write("\t".join(r) + "\n")


if __name__ == "__main__":
    main()
