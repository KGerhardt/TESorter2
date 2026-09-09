"""
nhmmer_search.py — Run stock HMMER's nhmmer binary for DNA profile databases.

Genome mode searches DNA profiles (AnnoSINE) by shelling out to HMMER's own
nhmmer rather than going through pyHMMER's LongTargetsPipeline. On genome-scale
targets the reference implementation is what the models were calibrated against
and what the field's results are comparable to; pyHMMER's binding is a
reimplementation of the same pipeline and needed a workaround for its E-value
scaling (see search.legacy_search_nucl, which pins Z=1 and rescales by hand
because pyHMMER applies Z twice). Short-sequence mode keeps the in-process
pyHMMER path, where per-call process startup would dominate.

nhmmer's --tblout is a flat 16-column table, one row per hit, with an explicit
strand column -- unlike pyHMMER's objects, which report descending coordinates
for minus-strand hits and leave the strand to be inferred.

Two fields the pipeline needs are absent from the table and are supplied here:

  query_len   the profile's length, read from the HMM file and keyed by model
              name. Coverage and normalized score both need it.
  acc         nhmmer's table carries no per-residue posterior column, so acc is
              set to 1.0 and the classifier's min_acc filter becomes a no-op for
              these hits -- the same convention bath_search and nail_search use.
              This costs nothing in practice: on the rice/sine fixture the
              accuracy filter rejects no hits at any threshold.
"""

import logging
import os
import shutil
import subprocess

from .hmm import open_hmm_text

log = logging.getLogger(__name__)

# Column positions in `nhmmer --tblout`. Validated against the header before
# use, so a format change fails loudly rather than mapping the wrong fields.
_EXPECTED_HEADER = ["target", "name", "accession", "query", "name",
                    "accession", "hmmfrom", "hmm", "to", "alifrom", "ali",
                    "to", "envfrom", "env", "to", "sq", "len", "strand",
                    "E-value", "score", "bias", "description", "of", "target"]
_N_COLUMNS = 16
_TARGET, _QUERY = 0, 2
_HMM_FROM, _HMM_TO = 4, 5
_ALI_FROM, _ALI_TO = 6, 7
_ENV_FROM, _ENV_TO = 8, 9
_SQ_LEN, _STRAND, _EVALUE, _SCORE, _BIAS = 10, 11, 12, 13, 14


def _bin():
    """Path to the nhmmer binary. NHMMER_BIN wins, then PATH."""
    explicit = os.environ.get("NHMMER_BIN")
    if explicit:
        if os.path.isfile(explicit):
            return explicit
        raise FileNotFoundError(
            "NHMMER_BIN is set to %r, which is not a file." % explicit)
    found = shutil.which("nhmmer")
    if found:
        return found
    raise FileNotFoundError(
        "nhmmer not found on PATH. It ships with HMMER, which is on bioconda "
        "('conda install -c bioconda hmmer'); or set NHMMER_BIN.")


def require_binaries():
    """Fail early and as a configuration error, not a traceback."""
    try:
        _bin()
    except FileNotFoundError as exc:
        raise SystemExit(str(exc))


def model_lengths(hmm_path):
    """{model name: profile length} straight from the HMM file."""
    lengths, name = {}, None
    with open_hmm_text(hmm_path) as f:
        for line in f:
            if line.startswith("NAME"):
                name = line.split(None, 1)[1].strip()
            elif line.startswith("LENG") and name is not None:
                lengths[name] = int(line.split()[1])
                name = None
    return lengths


def run_nhmmer(hmm_path, fasta, tblout, n_workers=4, dfam_e=None):
    """Run nhmmer, writing --tblout. Returns tblout.

    --nobias matches the rest of the pipeline (TEsorter's --nobias); -E is left
    permissive so the classifier's own 1e-3 cutoff, not the reporting
    threshold, decides what survives.
    """
    cmd = [_bin(), "--cpu", str(max(1, n_workers)), "--nobias",
           "-E", str(dfam_e if dfam_e is not None else 10.0),
           "--tblout", tblout, "-o", os.devnull, hmm_path, fasta]
    log.info("    nhmmer: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)
    return tblout


def _check_header(line, path):
    got = line.lstrip("#").split()
    if got[:len(_EXPECTED_HEADER)] != _EXPECTED_HEADER:
        raise ValueError(
            "%s: unexpected nhmmer --tblout layout %r; this parser reads "
            "columns positionally and expects %r. Update nhmmer_search before "
            "trusting these results." % (path, got, _EXPECTED_HEADER))


def parse_tblout(tblout, lengths):
    """Parse nhmmer --tblout into hit dicts in the shared internal schema.

    Coordinates are normalized the same way search._normalize_nucl_hit does:
    ascending, with the strand encoded as a |fwd1 / |rev1 target suffix, so
    results._parse_frame_info recovers base sequence and strand unchanged.
    """
    hits, checked = [], False
    with open(tblout, errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                if not checked and "target name" in line:
                    _check_header(line, tblout)
                    checked = True
                continue
            if not checked:
                raise ValueError(
                    "%s: no column header found; refusing to parse nhmmer "
                    "output by guessed offsets." % tblout)
            p = line.split()
            if len(p) < _N_COLUMNS:
                continue

            model = p[_QUERY]
            mlen = lengths.get(model)
            if mlen is None:
                log.warning("    nhmmer: no model length for %r, skipping", model)
                continue

            strand = p[_STRAND]
            ali = sorted((int(p[_ALI_FROM]), int(p[_ALI_TO])))
            env = sorted((int(p[_ENV_FROM]), int(p[_ENV_TO])))
            evalue, score, bias = float(p[_EVALUE]), float(p[_SCORE]), float(p[_BIAS])

            hits.append({
                "target_name": "%s|%s" % (p[_TARGET],
                                          "fwd1" if strand == "+" else "rev1"),
                "target_len": int(p[_SQ_LEN]),
                "query_name": model,
                "query_len": mlen,
                "evalue": evalue, "score": score, "bias": bias,
                "dom_num": 1, "dom_of": 1,
                "c_evalue": evalue, "i_evalue": evalue,
                "dom_score": score, "dom_bias": bias,
                "hmm_from": int(p[_HMM_FROM]), "hmm_to": int(p[_HMM_TO]),
                "ali_from": ali[0], "ali_to": ali[1],
                "env_from": env[0], "env_to": env[1],
                "acc": 1.0,
            })
    return hits


def run_and_parse(hmm_path, fasta, db_name, workdir, n_workers=4):
    """Run stock nhmmer over `fasta` and return parsed hit dicts."""
    os.makedirs(workdir, exist_ok=True)
    tblout = os.path.join(workdir, "%s.nhmmer.tblout" % db_name)
    run_nhmmer(hmm_path, fasta, tblout, n_workers=n_workers)
    hits = parse_tblout(tblout, model_lengths(hmm_path))
    log.info("    nhmmer: %d hits for %s", len(hits), db_name)
    return hits
