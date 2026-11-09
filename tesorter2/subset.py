"""
subset.py — write the subset of a FASTA that a cascade stage hands on.

Each stage of the search cascade labels some sequences and passes the rest to
the next engine. That handoff is a file: pyhmmer, BATH and nail all take a
FASTA path, so a FASTA on disk is the one interchange format they share, and
subsetting is how a stage narrows the work for the stage behind it.

Two properties matter:

  * Streaming. pyfastx.Fasta(build_index=False) iterates without building a
    .fxi sidecar, so a pass over a 387 MB library costs one read and no index.
    A cascade subsets once per stage per database; paying to index (or to hold
    the library in RAM) each time would cost more than the filtering saves.

  * Frame awareness. A translated FASTA carries six records per input sequence
    ({id}|fwd1 .. {id}|rev3), but pass/fail is decided per *sequence*, not per
    frame — classification already collapses frames, keying hits by base name
    (results._parse_frame_info). Membership is therefore tested on the base
    name so all six frames of a sequence move together. Names are normalized on
    both sides, so a caller may pass frame-suffixed names (straight out of a
    hits table) or bare sequence IDs interchangeably.

Nucleotide FASTAs are unaffected: a name with no |fwd/|rev suffix is its own
base name, so the same call works for the BATH/nhmmer inputs.
"""

import argparse
import sys

import pyfastx

from .sequence import parse_frame_suffix


def base_name(name):
    """Sequence ID with any six-frame suffix stripped.

    Safe on IDs that themselves contain '|' (e.g.
    'GCA_008658365|TE_00000554_INT#LTR/unknown'): parse_frame_suffix only
    strips the final field when it looks like fwd<n>/rev<n>.
    """
    return parse_frame_suffix(name)[0]


def subset_fasta(input_fasta, output_fasta, names, exclude=True, wrap=0):
    """Copy input_fasta to output_fasta, filtering records by base name.

    Args:
        input_fasta: FASTA to read (nucleotide or six-frame translated)
        output_fasta: FASTA to write
        names: iterable of sequence names; frame suffixes are stripped, so
            either bare IDs or frame-suffixed names work
        exclude: True (default) writes the records NOT in `names` — the
            cascade's usual direction, where `names` is what this stage
            resolved and the remainder goes to the next engine. False writes
            only the records that ARE in `names`.
        wrap: line width for sequence output; 0 (default) writes each sequence
            on one line, matching translate_fasta

    Returns:
        (n_written, n_skipped)
    """
    keys = {base_name(n) for n in names}
    n_written = n_skipped = 0

    with open(output_fasta, "w") as out:
        # build_index=False streams (name, seq) tuples and writes no .fxi.
        for name, seq in pyfastx.Fasta(input_fasta, build_index=False):
            if (base_name(name) in keys) == exclude:
                n_skipped += 1
                continue
            if wrap > 0:
                body = "\n".join(seq[i:i + wrap]
                                 for i in range(0, len(seq), wrap))
            else:
                body = seq
            out.write(">{}\n{}\n".format(name, body))
            n_written += 1

    return n_written, n_skipped


def read_names(path):
    """One name per line; blank lines and '#' comments ignored."""
    names = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                names.append(line)
    return names


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="tesorter2-subset",
        description="Subset a FASTA by sequence name, collapsing six-frame "
                    "suffixes so all frames of a sequence move together.")
    p.add_argument("input_fasta", help="FASTA to read")
    p.add_argument("names", help="File of sequence names, one per line")
    p.add_argument("output_fasta", help="FASTA to write")
    p.add_argument(
        "--keep", action="store_true", default=False,
        help="Write the records named in NAMES. Default writes the records "
             "NOT named — the cascade handoff, where NAMES is what the "
             "current stage resolved and the rest goes to the next engine.")
    p.add_argument("--wrap", type=int, default=0,
                   help="Wrap sequence lines at N characters [default: 0, "
                        "one line per sequence]")
    args = p.parse_args(argv)

    names = read_names(args.names)
    n_written, n_skipped = subset_fasta(
        args.input_fasta, args.output_fasta, names,
        exclude=not args.keep, wrap=args.wrap)
    sys.stderr.write(
        "{}: {} records written, {} filtered out ({} names)\n".format(
            args.output_fasta, n_written, n_skipped, len(names)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
