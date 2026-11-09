"""
cli.py -- the one command line shared by `TEsorter2`, `tesorter2` and
`TEsorter`.

It is TEsorter's command line, flag for flag, plus TEsorter2's own options.
A TEsorter command keeps working when only the executable name changes:

    TEsorter  in.fa -db rexdb-plant -p 8 -pre out
    TEsorter2 in.fa -db rexdb-plant -p 8 -pre out

There is no mode to name. As in TEsorter, element classification is the
default and `-genome` switches to whole-genome domain annotation.

`TEsorter2` with no -db searches every default database with the TEsorter2
engine cascade. `TEsorter` runs TEsorter's own single-database pipeline
(rexdb by default) and writes TEsorter's output files.

Accepted but left out of --help, so existing scripts keep working: the
2.0.0 subcommand words (`tesorter2 sequences ...`, `tesorter2 genome ...`),
`-d`/`--database` for -db, `--max-search` for `-db all`, and `--genome`.
"""

import argparse
import os
import sys
import uuid

from . import __version__


# TEsorter's -db names -> TEsorter2 database aliases. TEsorter2 bundles one
# REXdb file (viridiplantae v4.0 + metazoa v3.1), so the plant-only,
# metazoa-only and v3 names all search it; the name the user gave still
# names the output files, as it did in TEsorter.
TESORTER_DB_NAMES = {
    "rexdb":           "rexdb",
    "rexdb-plant":     "rexdb",
    "rexdb-metazoa":   "rexdb",
    "rexdb-v3":        "rexdb",
    "rexdb-plantv3":   "rexdb",
    "rexdb-metazoav3": "rexdb",
    "gydb":            "gydb",
    "rexdb-pnas":      "tir",
    "rexdb-line":      "line",
    "sine":            "sine",
}

# 2.0.0 subcommands, accepted as an optional first word.
MODES = ("sequences", "genome")

ALL_DBS = "all"


def _db_help(aliases):
    return ("Database(s) to search: a name, an HMM file, or a comma-separated "
            "list of them, or 'all' for every bundled database. Accepts "
            "TEsorter's names (%s) and TEsorter2's (%s). "
            "[default: TEsorter2 searches every bundled database except "
            "dfam-extended, dfam-deep and dfam-complete; TEsorter searches "
            "rexdb]"
            % (", ".join(TESORTER_DB_NAMES), ", ".join(aliases)))


def build_parser(prog="TEsorter2", classic=False, aliases=()):
    """The full parser, shared by both executables. `classic` (the
    `TEsorter` command) hides the options its pipeline does not use."""
    parser = argparse.ArgumentParser(
        prog=prog,
        description="Lineage-level classification of transposable elements "
                    "using conserved protein domains (TEsorter2 %s). Accepts "
                    "every TEsorter argument." % __version__)
    parser.add_argument("-v", "--version", action="version",
                        version="%s %s" % (prog, __version__))
    parser.add_argument("sequence",
                        help="Input TE/LTR or genome sequences in FASTA format "
                             "[required]")

    # ---- TEsorter's options, same names and defaults ----
    parser.add_argument("-db", "--hmm-database", dest="database",
                        default=None, metavar="DB", help=_db_help(aliases))
    # 2.0.0 spelling, kept working but not advertised.
    parser.add_argument("-d", "--database", dest="database",
                        help=argparse.SUPPRESS)
    parser.add_argument("--db-hmm", default=None, metavar="HMM",
                        help="Custom HMM database file, used instead of -db "
                             "[default: None]")
    parser.add_argument("-st", "--seq-type", choices=("nucl", "prot"),
                        default="nucl",
                        help="'nucl' for DNA or 'prot' for protein. 'prot' "
                             "takes amino acids directly: no six-frame "
                             "translation, no BATH, DNA profile databases "
                             "skipped [default: nucl]")
    parser.add_argument("-pre", "--prefix", default=None,
                        help="Output prefix, may include a directory "
                             "[default: '{input}.{db}' for one database, "
                             "'{input}.TEsorter2' for several]")
    parser.add_argument("-fw", "--force-write-hmmscan", action="store_true",
                        default=False,
                        help="Accepted for TEsorter compatibility. TEsorter2 "
                             "keeps no hmmscan output to reuse, so it always "
                             "searches.")
    parser.add_argument("-p", "--processors", type=int, default=4,
                        help="Processors to use [default: 4]")
    parser.add_argument("-tmp", "--tmp-dir", default=None,
                        help="Directory for temporary files "
                             "[default: ./tmp-{random id}]")
    parser.add_argument("-cov", "--min-coverage", type=float, default=20,
                        help="Minimum coverage of a protein domain's HMM "
                             "(0-100) [default: 20]")
    parser.add_argument("-eval", "--max-evalue", type=float, default=1e-3,
                        help="Maximum E-value for protein domains "
                             "[default: 1e-3]")
    parser.add_argument("-prob", "--min-probability", type=float, default=0.5,
                        help="Minimum posterior probability for protein "
                             "domains (0-1) [default: 0.5]")
    parser.add_argument("-score", "--min-score", type=float, default=0.1,
                        help="Minimum score per HMM position for protein "
                             "domains (0-2) [default: 0.1]")
    parser.add_argument("-mask", nargs="+", default=None,
                        choices=("soft", "hard"),
                        help="Write the input with its TE domains masked "
                             "(soft: lowercase; hard: N, or X for protein) "
                             "as {prefix}.softmasked / {prefix}.hardmasked. "
                             "Needs domain coordinates: -genome, or the TEsorter "
                             "command "
                             "[default: None]")
    parser.add_argument("-nocln", "--no-cleanup", dest="cleanup",
                        action="store_false", default=True,
                        help="Keep the temporary directory and intermediate "
                             "files. A run that fails keeps them anyway.")
    parser.add_argument("-cite", "--citation", action=_Citation, nargs=0,
                        help="Print the citation and exit")

    elem = parser.add_argument_group(
        "ELEMENT mode (default)",
        "Input TE/LTR sequences to classify them to clade level.")
    elem.add_argument("-dp2", "--disable-pass2", action="store_true",
                      default=False,
                      help="Do not further classify the unclassified sequences "
                           "by similarity. The TEsorter command always skips it "
                           "for -st prot [default: off]")
    elem.add_argument("-rule", "--pass2-rule", default="80-80-80",
                      metavar="I-C-L",
                      help="Pass-2 classifying rule, identity-coverage-length. "
                           "For the blast aligner: pident, qcovs and "
                           "alignment length. For minimap2: I is --min-pid, C "
                           "is both --min-qcov and --min-tcov, L is unused "
                           "[default: 80-80-80]")
    elem.add_argument("-nolib", "--no-library", action="store_true",
                      default=False,
                      help="Do not write the RepeatMasker library (.cls.lib)")
    elem.add_argument("-norc", "--no-reverse", action="store_true",
                      default=False,
                      help="Do not reverse-complement minus-strand sequences "
                           "in the library")

    gen = parser.add_argument_group(
        "GENOME mode",
        "Input genome sequences to detect TE protein domains throughout the "
        "whole genome. Protein databases are searched with BATH (required), "
        "DNA databases with nhmmer.")
    gen.add_argument("-genome", action="store_true", default=False,
                     help="Input is genome sequences [default: off]")
    gen.add_argument("--genome", dest="genome", action="store_true",
                     help=argparse.SUPPRESS)
    gen.add_argument("-win_size", type=int, default=None,
                     help="Accepted for TEsorter compatibility and ignored: "
                          "TEsorter2 searches each chromosome whole, so there "
                          "are no windows")
    gen.add_argument("-win_ovl", type=int, default=None,
                     help="Accepted for TEsorter compatibility and ignored "
                          "(see -win_size)")

    # ---- TEsorter2 options ----
    new = parser.add_argument_group("TEsorter2 options")
    if classic:
        # The TEsorter command runs TEsorter's pipeline, which reads only
        # these; the rest are accepted (one parser) but not advertised.
        used = {"-o", "--db-dir", "--pass2-aligner",
                "--pass2-classified-fasta", "--minimap2-extra"}
        add = new.add_argument

        def _add(*flags, **kw):
            if not used.intersection(flags):
                kw["help"] = argparse.SUPPRESS
            return add(*flags, **kw)
        new.add_argument = _add
    new.add_argument("--max-search", action="store_true", default=False,
                     help=argparse.SUPPRESS)            # = -db all
    new.add_argument("-o", "--outdir", default=None,
                     help="Directory for the output files; a relative -pre is "
                          "placed inside it [default: the directory of -pre, "
                          "else the working directory]")
    new.add_argument("--db-dir", default=None,
                     help="Directory containing the HMM databases "
                          "[default: the databases bundled with the package]")
    new.add_argument("--stages", default=None,
                     help="Ordered, comma-separated engine cascade for "
                          "amino-acid databases (hmmer, nail, bath). Each "
                          "stage searches only what earlier stages left "
                          "unresolved, so a fast engine strips the easy "
                          "majority. One name runs that engine alone. DNA "
                          "databases always use nhmmer. [default: "
                          "nail,hmmer,bath; a stage whose engine is not "
                          "installed is skipped with a warning, while an "
                          "explicit --stages list fails instead]")
    new.add_argument("--mask-stops", action="store_true", default=False,
                     help="Translate stop codons as 'X' rather than '*', so a "
                          "profile can align through a premature stop. "
                          "Recovers domains in degraded copies; its "
                          "false-positive cost is not yet characterized.")
    new.add_argument("--no-partition", dest="partition", action="store_false",
                     default=True,
                     help="Search the input as one unit instead of one "
                          "balanced group per processor (partitioning is "
                          "4.6x faster on 135k sequences over 64 cores)")
    new.add_argument("--memory-db", action="store_true", default=False,
                     help="Build the results database in RAM and write it to "
                          "disk at the end. Faster on parallel filesystems; "
                          "costs memory and loses the results if the run "
                          "dies.")
    new.add_argument("--extend", default=None, metavar="DB",
                     help="Add the -db databases to a finished run's results "
                          "database (DB) instead of starting a new one. "
                          "Databases it already searched are skipped, and the "
                          "vote, pass-2 and every output are redone over old "
                          "and new together. Give the same input and "
                          "classification options as the original run. "
                          "Outputs are named after DB unless -pre says "
                          "otherwise.")
    new.add_argument("--emit-timing", action="store_true", default=False,
                     help="Write cascade_timing.tsv and stage_timing.tsv "
                          "(wall and CPU time per stage). Diagnostic only.")
    new.add_argument("--compat-tesorter-voting", action="store_true",
                     default=False,
                     help="Decide the clade by raw domain-count plurality, as "
                          "TEsorter does. Default is a score-weighted vote.")
    new.add_argument("--compat-tesorter-rounding", action="store_true",
                     default=False,
                     help="Round normalized scores to 2 dp before the "
                          "threshold comparison, as TEsorter does.")
    new.add_argument("--compat-tesorter-output", action="store_true",
                     default=False,
                     help="Write the combined .cls.tsv in TEsorter's 7 "
                          "columns, without SecondaryHits/SO/Lineage.")
    new.add_argument("--no-tesorter-outputs", action="store_true",
                     default=False,
                     help="Skip the per-database RepeatMasker libraries")
    new.add_argument("--min-clade-delta", type=float, default=0.0,
                     metavar="BITS",
                     help="Demote a clade call to 'mixture' when a per-hit "
                          "score offset smaller than BITS could overturn the "
                          "vote. Only Clade is affected. 0 never demotes "
                          "[default: 0]")
    new.add_argument("--emit-clade-delta", action="store_true", default=False,
                     help="Append a CladeDelta column (that margin, in bits) "
                          "to the per-database .cls.tsv. Diagnostic only.")
    new.add_argument("--pass2-aligner", choices=("blast", "minimap2"),
                     default="blast",
                     help="Aligner for pass-2. 'blast' uses blastn (blastp "
                          "for -st prot) with clade=unknown; 'minimap2' "
                          "inherits the best target's full classification "
                          "[default: blast]")
    new.add_argument("--blast-task", choices=("megablast", "dc-megablast"),
                     default="megablast",
                     help="blastn -task for the blast pass-2 aligner "
                          "[default: megablast]")
    new.add_argument("--pass2-classified-fasta", default=None,
                     metavar="FASTA",
                     help="FASTA of previously classified elements added to "
                          "the pass-2 targets; headers like "
                          ">id#Order/Superfamily/Clade")
    new.add_argument("--minimap2-extra", default="", metavar="STR",
                     help="Extra flags passed through to minimap2")
    return parser


class _Citation(argparse.Action):
    """-cite prints and exits, before the positional is required."""

    def __call__(self, parser, namespace, values, option_string=None):
        print(CITATION)
        parser.exit()


def normalize_argv(argv):
    """Strip a 2.0.0 subcommand word; `genome` becomes -genome."""
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] in MODES:
        word = argv.pop(0)
        if word == "genome" and "-genome" not in argv and "--genome" not in argv:
            argv.append("-genome")
    return argv


def parse_args(argv=None, prog="TEsorter2", classic=False, aliases=()):
    parser = build_parser(prog=prog, classic=classic, aliases=aliases)
    args = parser.parse_args(normalize_argv(argv))
    args.classic = classic

    # TEsorter: protein input has no pass-2 and nothing to reverse-complement.
    # TEsorter2 itself runs a blastp pass-2 on protein input.
    if args.seq_type == "prot" and args.classic:
        args.disable_pass2 = True
        args.no_reverse = True
    if args.genome and args.seq_type == "prot":
        parser.error("-genome needs nucleotide input (-st nucl)")
    if args.tmp_dir is None:
        args.tmp_dir = "./tmp-%s" % uuid.uuid1()
        args.tmp_given = False
    else:
        args.tmp_given = True
    return args


def database_names(args):
    """(searched names, label) from -db/--db-hmm/--max-search.

    `searched` holds TEsorter2 aliases or file paths in order, without
    duplicates; `label` is the name for the default prefix: the name the user
    gave when there is one database, else None.
    """
    if args.db_hmm:
        given = [args.db_hmm]
    elif args.max_search:
        given = [ALL_DBS]
    elif args.database:
        given = [s.strip() for s in args.database.split(",") if s.strip()]
    else:
        given = []
    label = None
    if len(given) == 1 and given[0] != ALL_DBS:
        label = (os.path.splitext(os.path.basename(given[0]))[0]
                 if os.path.isfile(given[0]) else given[0])
        if label.endswith(".hmm"):            # name.hmm.gz
            label = label[:-4]
    names = []
    for g in given:
        name = TESORTER_DB_NAMES.get(g, g)
        if name not in names:
            names.append(name)
    return names, label


def output_paths(args, label, input_path):
    """(outdir, file prefix) following TEsorter's -pre convention."""
    prefix = args.prefix
    if prefix is None:
        prefix = "%s.%s" % (os.path.basename(input_path), label or "TEsorter2")
    if args.outdir and not os.path.isabs(prefix):
        prefix = os.path.join(args.outdir, prefix)
    outdir = os.path.dirname(prefix) or "."
    return outdir, os.path.basename(prefix)


def mask_input(input_fasta, gff3, prefix, types, seq_type="nucl"):
    """TEsorter's -mask: write {prefix}.{soft,hard}masked from a domain GFF3."""
    import logging
    import pyfastx
    log = logging.getLogger(__name__)

    regions = {}
    if os.path.exists(gff3):
        with open(gff3) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 5:
                    continue
                regions.setdefault(f[0], []).append((int(f[3]), int(f[4])))
    gap = "X" if seq_type == "prot" else "N"
    for kind in types:
        out = "%s.%smasked" % (prefix, kind)
        masked = total = 0
        with open(out, "w") as fo:
            for rec in pyfastx.Fasta(input_fasta, build_index=False,
                                     full_name=True):
                name, seq = rec[0], list(rec[1])
                sid = name.split()[0]
                ivs = regions.get(sid) or regions.get(sid.split("#")[0]) or []
                hit = bytearray(len(seq))
                for s, e in ivs:
                    lo, hi = max(min(s, e) - 1, 0), min(max(s, e), len(seq))
                    hit[lo:hi] = b"\x01" * (hi - lo)
                for i in range(len(seq)):
                    if hit[i]:
                        seq[i] = seq[i].lower() if kind == "soft" else gap
                masked += sum(hit)
                total += len(seq)
                fo.write(">%s\n" % name)
                s = "".join(seq)
                for i in range(0, len(s), 60):
                    fo.write(s[i:i + 60] + "\n")
        log.info("%s-masked %s -> %s: %d / %d (%.2f%%)", kind, input_fasta,
                 out, masked, total, 100.0 * masked / total if total else 0)


CITATION = """\
If you use TEsorter2, please cite it (paper in preparation):
Mahecha D, Gerhardt K, et al. TEsorter2. https://github.com/KGerhardt/TEsorter2

TEsorter2 builds on TEsorter:
Zhang RG, Li GL, Wang XL et al. TEsorter: an accurate and fast method to
classify LTR retrotransposons in plant genomes. Horticulture Research, 2022,
9: uhac017. https://doi.org/10.1093/hr/uhac017

If you use the REXdb database (rexdb, rexdb-plant, rexdb-metazoa), cite:
Neumann P, Novak P, Hostakova N et al. Systematic survey of plant
LTR-retrotransposons elucidates phylogenetic relationships of their polyprotein
domains and provides a reference for element classification. Mobile DNA, 2019,
10: 1. https://doi.org/10.1186/s13100-018-0144-1

If you use GyDB (gydb), cite:
Llorens C, Futami R, Covelli L et al. The Gypsy Database (GyDB) of mobile
genetic elements: release 2.0. Nucleic Acids Research, 2011, 39: D70-D74.
https://doi.org/10.1093/nar/gkq1061

If you use AnnoSINE (sine, sine-animals), cite:
Li Y, Jiang N, Sun Y. AnnoSINE: a short interspersed nuclear elements
annotation tool for plant genomes. Plant Physiology, 2022, 188: 955-970.
https://doi.org/10.1093/plphys/kiab524

If you use the LINE/RT database (line, rexdb-line), cite:
Kapitonov VV, Tempel S, Jurka J. Simple and fast classification of non-LTR
retrotransposons based on phylogeny of their RT domain protein sequences.
Gene, 2009, 448: 207-213. https://doi.org/10.1016/j.gene.2009.07.019

If you use the DNA/TIR database (tir, rexdb-pnas), cite:
Yuan YW, Wessler SR. The catalytic domain of all eukaryotic cut-and-paste
transposase superfamilies. PNAS, 2011, 108: 7884-7889.
https://doi.org/10.1073/pnas.1104208108

If you use the Dfam collections (dfam-*), cite:
Storer J, Hubley R, Rosen J, Wheeler TJ, Smit AF. The Dfam community resource
of transposable element families, sequence models, and genome annotations.
Mobile DNA, 2021, 12: 2. https://doi.org/10.1186/s13100-020-00230-y

If BATH ran (the default cascade and -genome), cite:
Krause GR, Shands W, Wheeler TJ. Sensitive and error-tolerant annotation of
protein-coding DNA with BATH. Bioinformatics Advances, 2024, 4: vbae088.
https://doi.org/10.1093/bioadv/vbae088
"""
