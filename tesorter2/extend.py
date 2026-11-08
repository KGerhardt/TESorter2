"""Add databases to a finished run's results database (`--extend`).

A run that searched the core set can take dfam-extended, or any other
database, without searching the core set again. The new databases' hits go
into the existing database through its own schema, the earlier databases'
calls are rebuilt from the hits already stored there, and everything that
depends on the whole set of calls -- the cross-database vote, the Dfam round's
own vote, BLAST pass-2 and the combined export -- is redone over old and new
together.

Rebuilding the calls from the stored hits, rather than reading back the
classifications table, is deliberate: the vote needs the per-call scores and
domain evidence, and only the hits carry them. Classification is per sequence
and deterministic, so classifying a database's stored hits stage by stage
gives back exactly the calls the run made.

What an extended run does NOT reproduce, by design: in a fresh run every Dfam
collection searches the same remainder and they vote among themselves, so a
deeper collection can join a call a shallower one made. Extended later, the
new collection searches only what is still unclassified, and the earlier
calls stand. That costs a few calls and saves searching what is already
called.

SQLite does only what it does well here: whole-table reads, appends, and a
numeric range on `id`. Everything else is Python.
"""

import logging
import math
import os
import sqlite3
from collections import defaultdict

import numpy as np

from .classifier import classify_sequences, DB_CONFIGS
from .hierarchical_search import DNA_STAGES
from .results import _INSERT_COLS, prepare_db

log = logging.getLogger(__name__)

# Rows of legacy_hits held in memory at once while the calls are rebuilt.
# The reload walks the table in ranges of sequence position (`id`), one
# sequential scan per range, so peak memory is one range's hits rather than
# the whole table's (8.4 M rows on Repbase).
RELOAD_ROWS = 2_000_000

_RELOAD_COLS = ("database, stage, engine, target_name, query_name, dom_score, "
                "i_evalue, acc, hmm_from, hmm_to, query_len, env_from, env_to, "
                "base_seq, domain_type")


def open_existing(db_path, memory=False):
    """Open a finished run's results database for extending, untouched.

    Nothing is changed until `begin` is called, so a run refused by the
    checks leaves the database exactly as it was, indexes and all. Under
    --memory-db the whole file is copied into RAM and written back at the
    end, the same as a fresh --memory-db run, and a run that dies leaves the
    file on disk as it was.
    """
    if not os.path.exists(db_path):
        raise SystemExit(f"--extend: {db_path} does not exist")
    if not memory:
        return sqlite3.connect(db_path)
    conn = sqlite3.connect(":memory:")
    src = sqlite3.connect(db_path)
    try:
        src.backup(conn)
    finally:
        src.close()
    return conn


def begin(conn):
    """Once the checks pass: bulk-insert tuning, and the finished run's
    indexes dropped so the new hits append to bare tables; finalize_db
    rebuilds them, as for a fresh run."""
    return prepare_db(conn)


def _tables(conn):
    return {r[0] for r in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table'")}


def check_compatible(conn, names):
    """Refuse a database this version cannot extend, or a different input.

    The new hits are inserted through the existing table, so its columns must
    be the ones this version writes. And every hit carries its sequence's
    position in the input, so the input must be the one the run was made on.
    """
    if "legacy_hits" not in _tables(conn):
        raise SystemExit("--extend: no legacy_hits table; not a TEsorter2 "
                         "sequences-mode results database")
    have = {r[1] for r in conn.execute("PRAGMA table_info(legacy_hits)")}
    want = {c.strip() for c in _INSERT_COLS.split(",")} | {"id"}
    missing = sorted(want - have)
    if missing:
        raise SystemExit("--extend: the results database predates this "
                         "version (legacy_hits lacks %s); rerun it from "
                         "scratch" % ", ".join(missing))
    stored = {r[0] for r in conn.execute("SELECT name FROM sequences")}
    if stored != set(names):
        raise SystemExit(
            "--extend: the input does not match the run being extended "
            "(%d sequences stored, %d in the input, %d in common)"
            % (len(stored), len(names), len(stored & set(names))))


def searched_databases(conn):
    """Every database the earlier run searched.

    hit_seqs names the ones that hit anything; cascade_exits also names the
    ones that hit nothing, since every searched database records its
    unresolved sequences there.
    """
    found = set()
    tables = _tables(conn)
    for table in ("hit_seqs", "cascade_exits"):
        if table in tables:
            found.update(r[0] for r in conn.execute(
                f"SELECT database FROM {table}"))
    return found


def reload_results(conn, n_sequences, compat_rounding=False,
                   compat_voting=False, min_clade_delta=0.0, filters=None):
    """Rebuild each searched database's calls from its stored hits.

    Returns ({database: [call, ...]}, {database: is_dna}). Calls are
    classified per (database, cascade stage), exactly as the cascade made
    them, and stamped with the engine and stage that produced them.

    The table is read in ranges of sequence position. A sequence's hits all
    carry its position, so a range holds whole sequences and classifying range
    by range is the same as classifying all at once; it only bounds memory.
    """
    n_rows = conn.execute("SELECT MAX(rowid) FROM legacy_hits").fetchone()[0]
    if not n_rows:
        return {}, {}
    n_ranges = max(1, math.ceil(n_rows / RELOAD_ROWS))
    step = math.ceil(max(n_sequences, 1) / n_ranges)
    log.info("  Reloading stored hits in %d range(s)", n_ranges)

    per_db = defaultdict(list)
    is_dna = {}
    for k in range(n_ranges):
        lo, hi = k * step, (k + 1) * step
        # The first range also takes id -1, "position unknown"; the last takes
        # everything above, so no row falls outside every range.
        if n_ranges == 1:
            rows = conn.execute(f"SELECT {_RELOAD_COLS} FROM legacy_hits")
        elif k == 0:
            rows = conn.execute(f"SELECT {_RELOAD_COLS} FROM legacy_hits "
                                "WHERE id < ?", (hi,))
        elif k == n_ranges - 1:
            rows = conn.execute(f"SELECT {_RELOAD_COLS} FROM legacy_hits "
                                "WHERE id >= ?", (lo,))
        else:
            rows = conn.execute(f"SELECT {_RELOAD_COLS} FROM legacy_hits "
                                "WHERE id >= ? AND id < ?", (lo, hi))
        groups = defaultdict(list)
        engines = {}
        for r in rows:
            key = (r[0], r[1])
            groups[key].append(r[3:])
            engines[key] = r[2]
        for (db, stage), g in groups.items():
            engine = engines[(db, stage)]
            is_dna[db] = engine in DNA_STAGES
            if DB_CONFIGS.get(db) is None:
                log.warning("  No classifier config for stored database %s; "
                            "its hits are kept but it cannot vote", db)
                continue
            results = classify_sequences(
                _arrays(g), DB_CONFIGS[db], compat_rounding=compat_rounding,
                compat_voting=compat_voting, min_clade_delta=min_clade_delta,
                filters=filters)
            for r in results:
                r["engine"] = engine
                r["stage"] = stage
            per_db[db].extend(results)
        del groups
    return dict(per_db), is_dna


def _arrays(rows):
    """The layout classify_sequences reads; mirrors hits_to_arrays."""
    (target, model, score, evalue, acc, hmm_from, hmm_to, model_len,
     env_from, env_to, base, family) = zip(*rows)
    hmm_from = np.array(hmm_from, dtype=np.int32)
    hmm_to = np.array(hmm_to, dtype=np.int32)
    model_len = np.array(model_len, dtype=np.int32)
    score = np.array(score, dtype=np.float64)
    return {
        "target": np.array(target), "model": np.array(model),
        "score": score, "evalue": np.array(evalue, dtype=np.float64),
        "acc": np.array(acc, dtype=np.float64),
        "hmm_from": hmm_from, "hmm_to": hmm_to, "model_len": model_len,
        "env_from": np.array(env_from, dtype=np.int32),
        "env_to": np.array(env_to, dtype=np.int32),
        "base_seq": np.array(base), "family": np.array(family),
        "hmm_cov": 100.0 * (hmm_to - hmm_from + 1) / model_len,
        "norm_score": score / model_len,
    }


def clear_pass2(conn):
    """Drop the earlier run's BLAST pass-2 so it can be redone.

    Pass-2 searches whatever has no HMM hit against everything that has one,
    so new hits change both sides and the old result cannot be kept. blast_hits
    is pass-2's alone and is dropped whole. Its calls share the
    classifications table with the per-database calls, so that table is read
    back, filtered in Python and rewritten under its own DDL.
    """
    conn.execute("DROP TABLE IF EXISTS blast_hits")
    if "classifications" in _tables(conn):
        ddl = conn.execute("SELECT sql FROM sqlite_master WHERE type='table' "
                           "AND name='classifications'").fetchone()[0]
        cols = [r[1] for r in conn.execute(
            "PRAGMA table_info(classifications)")]
        db_pos = cols.index("database")
        rows = conn.execute("SELECT * FROM classifications").fetchall()
        keep = [r for r in rows if r[db_pos] != "blast_pass2"]
        conn.execute("DROP TABLE classifications")
        conn.execute(ddl)
        conn.executemany("INSERT INTO classifications VALUES (%s)"
                         % ", ".join("?" * len(cols)), keep)
        log.info("  Cleared %d BLAST pass-2 calls from the earlier run",
                 len(rows) - len(keep))
    conn.commit()
