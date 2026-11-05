"""
bath_search.py — Run the BATH aligner as a drop-in replacement for the
pyhmmer/HMMER search on amino-acid databases.

BATH (bathsearch) performs a frameshift-aware translated search of protein
HMMs directly against *nucleotide* sequences, so unlike the HMMER path it does
NOT need six-frame translation of the input. This module:

  1. Resolves a BATH-format HMM for a given repo HMMER3 database (reusing a
     cached / pre-converted .bath.hmm when available, else running bathconvert).
  2. Runs bathsearch (--fs) against the raw nucleotide FASTA.
  3. Parses the tblout into hit dicts matching the schema that
     results.store_legacy expects, so the existing classifier runs unchanged.

BATH's tblout has no posterior-probability ("acc") column, so acc is set to
1.0 — the min_acc>=0.5 filter becomes a no-op for BATH hits, and effective
filtering happens on E-value / coverage / normalized score (the same fields
TEsorter's plant benchmark filtered on).
"""

import logging
import os
import shutil
import subprocess
import tempfile

log = logging.getLogger(__name__)

# BATH publishes no release artifacts (the repo has tags but zero releases), so
# a source build is the only install route.
_BUILD_HELP = (
    "build it with:\n"
    "  git clone https://github.com/TravisWheelerLab/BATH.git\n"
    "  cd BATH && git clone -b BATH https://github.com/TravisWheelerLab/easel\n"
    "  autoconf && ./configure && make\n"
    "then add BATH/src to PATH or set BATH_BIN_DIR to it")

# Fallback for the original development environment. Only consulted when
# BATH_BIN_DIR is unset and the binary is not on PATH.
_LEGACY_BATH_BIN_DIR = "/anvil/projects/x-bio250374/daniel/software/BATH_2/BATH/src"


def _bin(name):
    """Absolute path to a BATH binary.

    Resolution order: BATH_BIN_DIR (an explicit setting always wins, and fails
    loudly rather than silently falling through), then PATH, then the legacy
    development directory.
    """
    bin_dir = os.environ.get("BATH_BIN_DIR")
    if bin_dir:
        path = os.path.join(bin_dir, name)
        if os.path.isfile(path):
            return path
        raise FileNotFoundError(
            f"BATH binary '{name}' not found in BATH_BIN_DIR ({bin_dir}).")

    found = shutil.which(name)
    if found:
        return found

    path = os.path.join(_LEGACY_BATH_BIN_DIR, name)
    if os.path.isfile(path):
        return path

    raise FileNotFoundError(
        f"BATH binary '{name}' not found on PATH. BATH is not packaged on "
        f"conda and ships no prebuilt binaries -- {_BUILD_HELP}.")


def require_binaries():
    """Fail early if BATH is missing.

    Genome mode searches protein profiles exclusively with BATH, so this runs
    before any expensive setup rather than letting the run window a genome and
    then die at the first search. Re-raised as SystemExit: a missing external
    dependency is a user-facing configuration problem, not a crash, so it
    should print the instructions rather than a traceback.
    """
    for name in ("bathsearch", "bathconvert"):
        try:
            _bin(name)
        except FileNotFoundError as exc:
            raise SystemExit(str(exc))


# There used to be a KNOWN_CONVERTED shortcut here mapping db aliases to
# pre-converted BATH HMMs under a specific user's Anvil project directory. It
# was removed: those paths exist on that machine, so the shortcut silently won
# over the bundled databases and searched a third party's copy instead of the
# one shipped with this package -- including after GyDB2 was corrected in
# place. A database substituted without the caller knowing is worse than a
# conversion that takes a minute. The cache / convert path below is
# self-contained and always derives from the database actually being used.


def _cache_path(hmm_path):
    """Where the converted BATH HMM for hmm_path should live.

    Beside the source HMM when that directory is writable, which keeps the
    conversion next to what it was derived from. The bundled databases ship
    inside the package, though, and an installed package is frequently
    read-only (site-packages, a conda env, a shared module tree) -- so fall
    back to a user cache directory rather than failing to convert.

    The fallback key includes the source's size and mtime so a database that
    is updated in place does not silently keep serving a stale conversion.
    """
    directory = os.path.dirname(hmm_path) or "."
    if os.access(directory, os.W_OK):
        return hmm_path + ".bath.hmm"

    root = os.environ.get("TESORTER2_CACHE") or os.path.join(
        os.path.expanduser("~"), ".cache", "tesorter2", "bath")
    os.makedirs(root, exist_ok=True)
    st = os.stat(hmm_path)
    return os.path.join(root, "%s.%d.%d.bath.hmm" % (
        os.path.basename(hmm_path), st.st_size, int(st.st_mtime)))


def resolve_bath_db(hmm_path, db_name=None):
    """Return a path to a BATH-format HMM for the given repo HMMER3 database.

    Resolution order:
      1. Cached conversion (see _cache_path: beside the source HMM when that
         directory is writable, else under a user cache dir) — reuse if present.
      2. Run bathconvert on hmm_path, writing the cache, and use that.
    """
    cache = _cache_path(hmm_path)
    if os.path.isfile(cache):
        log.info(f"  BATH HMM (cached): {cache}")
        return cache

    log.info(f"  Converting {os.path.basename(hmm_path)} to BATH format -> {cache}")
    # Write to a temp file first so an interrupted convert never leaves a
    # half-written cache that a later run would happily reuse.
    fd, tmp = tempfile.mkstemp(suffix=".bath.hmm", dir=os.path.dirname(cache) or ".")
    os.close(fd)
    try:
        subprocess.run([_bin("bathconvert"), tmp, hmm_path], check=True)
        shutil.move(tmp, cache)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    return cache


def run_bathsearch(bath_hmm, nucl_fasta, tblout_path, n_workers=4,
                   frameshift=True, report_evalue=0.01, z_megabases=None):
    """Run bathsearch, writing a parseable tblout. Returns tblout_path.

    frameshift=True enables --fs (BATH's frameshift-aware mode), which is the
    whole reason to prefer BATH over HMMER for degraded elements.
    report_evalue caps the reported hits to keep the tblout small; the real
    TEsorter cutoff (1e-3) is applied later by the classifier's filters.

    z_megabases pins the search space (bathsearch -Z, in megabases) instead of
    letting BATH derive it from whatever was passed in. This matters whenever
    the input is a subset: E-value scales with search-space size, so searching
    30 of 60 sequences yields better E-values for the same alignment, and hits
    sitting near the 1e-3 cutoff flip purely on how much an upstream stage
    removed. A staged cascade must therefore pin Z to the full input so every
    stage's E-values mean the same thing and match a single full-input run.
    (search.legacy_search_nucl pins nhmmer's Z for the same reason; hmmsearch
    needs no equivalent, as its Z is a model count, not a target count.)
    """
    cmd = [_bin("bathsearch")]
    if frameshift:
        cmd.append("--fs")
    cmd += ["--cpu", str(max(1, n_workers)),
            "-E", str(report_evalue)]
    if z_megabases is not None:
        cmd += ["-Z", repr(float(z_megabases))]
    cmd += ["--tblout", tblout_path,
            "-o", os.devnull,
            bath_hmm, nucl_fasta]
    log.info(f"  bathsearch: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    return tblout_path


# Multi-token column names, longest first so greedy matching prefers the
# longest phrase. Everything else in the header is a single token.
_HEADER_PHRASES = (
    ("description", "of", "target"),
    ("hit", "id"), ("target", "name"), ("query", "name"),
    ("hmm", "len"), ("hmm", "from"), ("hmm", "to"),
    ("seq", "len"), ("ali", "from"), ("ali", "to"),
    ("env", "from"), ("env", "to"),
)

# Columns the hit dict cannot be built without.
_REQUIRED_COLUMNS = ("target name", "query name", "hmm len", "hmm from",
                     "hmm to", "seq len", "ali from", "ali to", "e-value",
                     "score")


def _tblout_columns(name_line):
    """Map column name -> token index from the tblout's column-name line.

    BATH's column set is not stable across versions. 1.x reports env
    coordinates; 2.0 drops them, prepends a numeric "hit ID", and appends
    PID / shifts / stops. Reading by name keeps one parser correct for both,
    and for whatever the next release moves.

    Matching is done on the token stream rather than by slicing at the rule
    line's dash positions: column values are right-aligned and headers are not
    reliably within their own dash run, so position-slicing shears multi-word
    names apart ("env fro", "m    env t") and then fails silently. Greedy
    longest-phrase matching is whitespace-independent instead.

    Every column before the free-text description is a single whitespace-free
    token in the data rows, so token index == column index.
    """
    tokens = [t.lower() for t in name_line.lstrip("#").split()]
    cols, i, idx = {}, 0, 0
    while i < len(tokens):
        for phrase in _HEADER_PHRASES:
            if tuple(tokens[i:i + len(phrase)]) == phrase:
                cols.setdefault(" ".join(phrase), idx)
                i += len(phrase)
                break
        else:
            # Two columns are both "accession"; neither is read, and
            # first-wins keeps the map unambiguous.
            cols.setdefault(tokens[i], idx)
            i += 1
        idx += 1
    return cols


def parse_bath_tblout(tblout_path):
    """Parse a BATH tblout into hit dicts for results.store_legacy.

    Columns are located by name (see _tblout_columns) rather than by fixed
    offset. When the engine reports no env coordinates -- BATH 2.0 does not --
    the ali coordinates stand in, which is what the downstream geometry wants
    anyway: for BATH the reported alignment IS the domain envelope.
    """
    hits = []
    cols = None
    with open(tblout_path) as f:
        for line in f:
            if line.startswith("#"):
                # The column-name line is the one naming the target column.
                if cols is None:
                    found = _tblout_columns(line)
                    if "target name" in found:
                        missing = [c for c in _REQUIRED_COLUMNS
                                   if c not in found]
                        if missing:
                            raise ValueError(
                                f"{tblout_path}: BATH tblout header is missing "
                                f"required column(s) {missing}. Parsed header: "
                                f"{sorted(found)}")
                        cols = found
                continue
            if cols is None:
                # Never fall back to guessing offsets: a mis-read header would
                # silently yield zero hits, which looks like "BATH found
                # nothing" rather than a parse failure.
                raise ValueError(
                    f"{tblout_path}: no column header found; refusing to parse "
                    f"BATH output by guessed offsets.")

            p = line.split()

            def col(name, cast=str, default=None):
                idx = cols.get(name)
                if idx is None or idx >= len(p):
                    return default
                return cast(p[idx])

            seq_name = col("target name")
            hmm_name = col("query name")
            hmm_len = col("hmm len", int)
            if seq_name is None or hmm_name is None or hmm_len is None:
                continue
            hmm_from = col("hmm from", int)
            hmm_to = col("hmm to", int)
            seq_len = col("seq len", int)
            ali_from = col("ali from", int)
            ali_to = col("ali to", int)
            # BATH 2.0 reports no envelope; the alignment is the envelope.
            env_from = col("env from", int, ali_from)
            env_to = col("env to", int, ali_to)
            evalue = col("e-value", float)
            score = col("score", float)
            bias = col("bias", float, 0.0)
            if None in (hmm_from, hmm_to, seq_len, ali_from, ali_to,
                        evalue, score):
                continue

            # BATH reports descending coords on the minus strand. Encode strand
            # in the target suffix (matching results._parse_frame_info's
            # |fwd1/|rev1 convention) and normalize coords to ascending so the
            # downstream overlap / ordering geometry matches the HMMER path.
            strand_suffix = "fwd1" if ali_from <= ali_to else "rev1"
            ali_lo, ali_hi = sorted((ali_from, ali_to))
            env_lo, env_hi = sorted((env_from, env_to))

            hits.append({
                "target_name": f"{seq_name}|{strand_suffix}",
                "target_len": seq_len,
                "query_name": hmm_name,
                "query_len": hmm_len,
                "evalue": evalue,
                "score": score,
                "bias": bias,
                "dom_num": 1,
                "dom_of": 1,
                "c_evalue": evalue,
                "i_evalue": evalue,
                "dom_score": score,
                "dom_bias": bias,
                "hmm_from": hmm_from,
                "hmm_to": hmm_to,
                "ali_from": ali_lo,
                "ali_to": ali_hi,
                "env_from": env_lo,
                "env_to": env_hi,
                # BATH tblout has no posterior column; 1.0 makes the acc filter
                # a no-op (see module docstring).
                "acc": 1.0,
            })
    return hits


def run_and_parse(hmm_path, nucl_fasta, db_name, n_workers=4, outdir=None,
                  frameshift=True):
    """Resolve the BATH HMM, run bathsearch, and return parsed hit dicts."""
    bath_hmm = resolve_bath_db(hmm_path, db_name=db_name)

    tmp_tblout = None
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        tblout = os.path.join(outdir, f"{db_name}.bath.tblout")
    else:
        fd, tblout = tempfile.mkstemp(suffix=".bath.tblout")
        os.close(fd)
        tmp_tblout = tblout

    try:
        run_bathsearch(bath_hmm, nucl_fasta, tblout, n_workers=n_workers,
                       frameshift=frameshift)
        hits = parse_bath_tblout(tblout)
    finally:
        if tmp_tblout and os.path.exists(tmp_tblout):
            os.remove(tmp_tblout)
    log.info(f"  BATH: {len(hits)} hits for {db_name}")
    return hits
