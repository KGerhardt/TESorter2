"""
nail_search.py — Run nail as a search engine over six-frame translations.

nail (https://github.com/TravisWheelerLab/nail) seeds with MMseqs2 and then
computes a sparse approximation of HMMER's Forward/Backward, recovering most of
HMMER's recall at a fraction of the runtime. It reads a HMMER3 .hmm directly as
its query, so unlike BATH there is no conversion step.

It is amino-acid only, which makes it a candidate for the *first* stage of a
short-sequence cascade: strip the sequences it can label confidently, hand the
rest to a slower and more sensitive engine. It cannot search DNA profile
databases at all, and it has no role in genome mode.

Two integration details that are not in nail's documentation:

**Stop codons must be removed.** nail rejects '*' outright ("byte 42 is not in
the compressed amino acid alphabet"), and a six-frame translation is full of
them -- 24,072 in a 60-sequence test library. Targets are therefore rewritten
with '*' -> 'X' before the search. The rewrite is nail-local: the shared
translated FASTA keeps its stop codons for the pyhmmer path, which handles them
natively.

That substitution is NOT neutral, and anything comparing nail against another
engine has to control for it. X is an unknown residue rather than a stop, so a
profile can align straight through a premature stop instead of being cut short
by it. Measured on the TIR fixture: pyhmmer classifies 15 sequences on
stop-bearing targets and 17 on the same targets with '*' -> 'X', and one
sequence goes from no hit at all to score 20.5 over hmm 59-148. Handing nail
X-substituted targets while pyhmmer reads stop-bearing ones therefore flatters
nail by exactly that margin: on identical X-substituted input nail produces no
call pyhmmer does not, which is what nail's own benchmarks predict (it recovers
most of the MMseqs2-to-HMMER recall gap, so it sits at or below HMMER).

The gain from X-substitution is a separate finding worth its own experiment: it
is a cheap partial substitute for BATH's read-through-stops behaviour, though
only the in-frame part -- it cannot change frame at an indel, which is the rest
of what BATH does.

**Model length is not in the output.** nail reports query start/end but never
the profile's length, which the classifier needs for coverage and normalized
score. It is taken from the HMM file instead, keyed by model name.

Caveat -- there is no way to pin nail's search space. nail exposes no -Z
equivalent, so its E-values scale with the size of the target set it is handed,
exactly as BATH's do without -Z. As a *first* cascade stage that is harmless,
since stage 0 sees the whole input. Placed later, its verdicts would shift with
however much upstream stages removed, and nothing in its CLI can correct for
that.
"""

import logging
import os
import shutil
import subprocess

log = logging.getLogger(__name__)

_INSTALL_HELP = (
    "install it with 'cargo install nail' (and ensure MMseqs2 is on PATH, "
    "which nail shells out to for seeding), then put it on PATH or set "
    "NAIL_BIN")

# The leaf row of nail's wrapped three-line header. Parsing is positional --
# the column names span two lines, so there is nothing clean to match on -- but
# this is checked first so a format change fails loudly instead of silently
# mapping the wrong fields, which is how the BATH parser went wrong.
_EXPECTED_COLUMNS = ["target", "query", "start", "end", "start", "end",
                     "score", "bias", "evalue", "frac"]
_N_COLUMNS = 10

# Positions within a data row, matching _EXPECTED_COLUMNS.
_TARGET, _QUERY = 0, 1
_TARGET_START, _TARGET_END = 2, 3
_QUERY_START, _QUERY_END = 4, 5
_SCORE, _BIAS, _EVALUE = 6, 7, 8


def _bin():
    """Absolute path to the nail binary.

    NAIL_BIN wins if set (and fails loudly rather than falling through), then
    PATH, then cargo's default install location -- 'cargo install nail' puts it
    in ~/.cargo/bin, which is not always on PATH in a non-login shell.
    """
    explicit = os.environ.get("NAIL_BIN")
    if explicit:
        if os.path.isfile(explicit):
            return explicit
        raise FileNotFoundError("NAIL_BIN is set to %r, which is not a file."
                                % explicit)

    found = shutil.which("nail")
    if found:
        return found

    cargo = os.path.join(os.path.expanduser("~"), ".cargo", "bin", "nail")
    if os.path.isfile(cargo):
        return cargo

    raise FileNotFoundError("nail binary not found on PATH -- %s."
                            % _INSTALL_HELP)


def require_binaries():
    """Fail early, and as a configuration error rather than a traceback.

    Checks MMseqs2 too: nail invokes it for seeding and fails mid-run without
    it, which surfaces as an opaque seeding error rather than a missing
    dependency.
    """
    try:
        _bin()
    except FileNotFoundError as exc:
        raise SystemExit(str(exc))
    if not (os.environ.get("MMSEQS_BIN") or shutil.which("mmseqs")):
        raise SystemExit(
            "MMseqs2 not found on PATH. nail uses 'mmseqs search' for seeding "
            "and cannot run without it; install it (it is on bioconda) or set "
            "MMSEQS_BIN.")


_STOP_TO_X = bytes.maketrans(b"*", b"X")


def _scan(fasta):
    """One binary read pass: {name: residue count} and the total stop count.

    Lengths are needed because nail's output does not report target length and
    the hit schema does. Binary mode avoids decoding several hundred MB of
    sequence that is only ever compared byte-wise.
    """
    lengths, stops, name = {}, 0, None
    with open(fasta, "rb") as f:
        for line in f:
            if line[:1] == b">":
                name = line[1:].split()[0].decode("utf-8", "replace")
                lengths[name] = 0
            elif name is not None:
                n = len(line) - (1 if line.endswith(b"\n") else 0)
                lengths[name] += n
                stops += line.count(b"*", 0, n)
    return lengths, stops


def prepare_targets(in_fasta, out_fasta):
    """Return (path for nail to search, {name: length}).

    nail rejects '*', so a translated FASTA carrying stop codons has to be
    rewritten with 'X'. Two things keep that cheap:

    A file with no stop codons is handed to nail untouched and no copy is
    written. That is the normal cascade case -- --mask-stops is mandatory
    whenever nail shares a cascade with another amino-acid engine, so the
    shared translated FASTA already says X. Copying a 300 MB file to change
    nothing took over three minutes on the RepBase LTR set.

    When a rewrite is needed it runs in binary with bytes.translate, which is a
    C-level table substitution, rather than decoding and re-encoding every line.
    Header lines are passed through untouched so an id containing '*' survives.
    """
    lengths, stops = _scan(in_fasta)
    if stops == 0:
        log.info("    nail: no stop codons, searching the input directly")
        return in_fasta, lengths
    with open(in_fasta, "rb") as fin, open(out_fasta, "wb") as fout:
        for line in fin:
            fout.write(line if line[:1] == b">" else line.translate(_STOP_TO_X))
    log.info("    nail: rewrote %d stop codons as X", stops)
    return out_fasta, lengths


def run_nail(hmm_path, target_fasta, tbl_out, workdir, n_workers=4,
             report_evalue=10.0):
    """Run `nail search`, writing tabular output. Returns tbl_out."""
    cmd = [_bin(), "search",
           "-t", str(max(1, n_workers)),
           "-X",                      # allow overwriting a previous run
           "-E", str(report_evalue),
           "--tbl-out", tbl_out,
           "--tmp-dir", os.path.join(workdir, "nail_tmp")]
    mmseqs = os.environ.get("MMSEQS_BIN")
    if mmseqs:
        cmd += ["--mmseqs-path", mmseqs]
    cmd += [hmm_path, target_fasta]
    log.info("    nail: %s", " ".join(cmd))
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    return tbl_out


def _check_header(rule_line, name_line, path):
    """Validate the tabular layout before trusting positional parsing."""
    groups = [g for g in rule_line.lstrip("#").split() if set(g) == {"-"}]
    names = name_line.lstrip("#").split()
    if len(groups) != _N_COLUMNS or names != _EXPECTED_COLUMNS:
        raise ValueError(
            "%s: unexpected nail tabular layout (%d columns %r); this parser "
            "reads columns positionally and expects %d %r. nail's format has "
            "changed -- update nail_search before trusting these results."
            % (path, len(groups), names, _N_COLUMNS, _EXPECTED_COLUMNS))


def parse_nail_tbl(tbl_path, model_lengths, target_lengths):
    """Parse nail tabular output into hit dicts in the shared internal schema.

    model_lengths maps model name -> profile length (from the HMM file; nail
    does not report it). target_lengths maps frame name -> residue count.
    """
    hits = []
    checked = False
    prev = ""
    with open(tbl_path) as f:
        for line in f:
            if line.startswith("#"):
                if not checked and set(line.strip()) <= set("#- ") and "-" in line:
                    _check_header(line, prev, tbl_path)
                    checked = True
                prev = line
                continue
            if not checked:
                raise ValueError(
                    "%s: no nail column header found; refusing to parse "
                    "positionally without confirming the layout." % tbl_path)

            p = line.split()
            if len(p) != _N_COLUMNS:
                continue

            target = p[_TARGET]
            model = p[_QUERY]
            model_len = model_lengths.get(model)
            if model_len is None:
                # A hit against a model that is not in the HMM we searched
                # cannot be scored for coverage; skipping is safer than
                # inventing a length.
                log.warning("    nail: no model length for %r, skipping", model)
                continue

            t_start, t_end = int(p[_TARGET_START]), int(p[_TARGET_END])
            q_start, q_end = int(p[_QUERY_START]), int(p[_QUERY_END])
            evalue = float(p[_EVALUE])
            score = float(p[_SCORE])
            bias = float(p[_BIAS])

            hits.append({
                "target_name": target,
                "target_len": target_lengths.get(target, t_end),
                "query_name": model,
                "query_len": model_len,
                "evalue": evalue,
                "score": score,
                "bias": bias,
                "dom_num": 1,
                "dom_of": 1,
                "c_evalue": evalue,
                "i_evalue": evalue,
                "dom_score": score,
                "dom_bias": bias,
                "hmm_from": q_start,
                "hmm_to": q_end,
                "ali_from": t_start,
                "ali_to": t_end,
                # nail reports one alignment per hit and no separate envelope.
                "env_from": t_start,
                "env_to": t_end,
                # No posterior-probability column, as with BATH; 1.0 makes the
                # classifier's min_acc filter a no-op for nail hits.
                "acc": 1.0,
            })
    return hits


def run_and_parse(hmm_path, aa_fasta, db_name, model_lengths, workdir,
                  n_workers=4, report_evalue=10.0):
    """Sanitize, search, and parse. Returns hit dicts."""
    os.makedirs(workdir, exist_ok=True)
    target = os.path.join(workdir, "%s.nail_targets.faa" % db_name)
    target, target_lengths = prepare_targets(aa_fasta, target)
    tbl = os.path.join(workdir, "%s.nail.tbl" % db_name)
    run_nail(hmm_path, target, tbl, workdir, n_workers=n_workers,
             report_evalue=report_evalue)
    hits = parse_nail_tbl(tbl, model_lengths, target_lengths)
    log.info("    nail: %d hits for %s", len(hits), db_name)
    return hits
