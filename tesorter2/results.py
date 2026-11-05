"""
SQLite persistence for HMM search results.

Stores domain hits from pass 1 (coarse screen) and pass 2 (sensitive)
in separate tables. Exports clean tab-separated flat files.
"""

import sqlite3

from .sequence import parse_frame_suffix, aa_to_nucl_coords, load_sequences_dict


_HIT_COLUMNS = """
    id          INTEGER NOT NULL,
    database    TEXT NOT NULL,
    target_name TEXT NOT NULL,
    base_seq    TEXT NOT NULL,
    strand      TEXT NOT NULL,
    frame       INTEGER NOT NULL,
    target_len  INTEGER NOT NULL,
    query_name  TEXT NOT NULL,
    domain_type TEXT NOT NULL,
    query_len   INTEGER NOT NULL,
    evalue      REAL NOT NULL,
    score       REAL NOT NULL,
    bias        REAL NOT NULL,
    dom_num     INTEGER NOT NULL,
    dom_of      INTEGER NOT NULL,
    c_evalue    REAL NOT NULL,
    i_evalue    REAL NOT NULL,
    dom_score   REAL NOT NULL,
    dom_bias    REAL NOT NULL,
    hmm_from    INTEGER NOT NULL,
    hmm_to      INTEGER NOT NULL,
    ali_from    INTEGER NOT NULL,
    ali_to      INTEGER NOT NULL,
    env_from    INTEGER NOT NULL,
    env_to      INTEGER NOT NULL,
    acc         REAL NOT NULL,
    search_mode INTEGER NOT NULL DEFAULT 0,
    engine      TEXT NOT NULL DEFAULT '',
    stage       INTEGER NOT NULL DEFAULT 0
"""
# `id` is the hit's sequence position in the original input, not a row
# identifier and deliberately not the rowid. Making it INTEGER PRIMARY KEY
# aliased it onto the rowid, which meant inserting a row dictated where it
# landed in the B-tree: a partitioned run merging groups in completion order
# wrote non-monotonic rowids and paid page splits instead of appending. The
# rowid is now arbitrary and never touched, so every insert appends and merge
# order stops mattering. `id` carries what is actually useful -- where the
# sequence sat in the input -- and is stable no matter how the input was split.
#
# engine/stage are the cascade's provenance: which tool produced a hit and at
# which position in the per-database cascade. A single-engine run leaves engine
# at the tool's own name and stage at 0. They are on the hit rows rather than
# inferred later because a cascade run mixes engines within one database, and
# scores from different engines are not on a comparable scale -- any analysis
# that compares them has to be able to tell them apart.

SCHEMA = f"""
CREATE TABLE IF NOT EXISTS sequences (
    name        TEXT PRIMARY KEY,
    length      INTEGER NOT NULL
);

CREATE TABLE IF NOT EXISTS legacy_hits (
    {_HIT_COLUMNS}
);

-- Which (sequence, database) pairs produced any hit at all: the only thing
-- BLAST pass-2 needs in order to pick its targets, and the only question
-- anything asks of the hit rows during a run.
--
-- It exists so that nobody has to ask sqlite for it. The old route was
-- `SELECT DISTINCT base_seq, database FROM legacy_hits`, which is a scan of
-- every hit in the run -- 8.38 M rows to learn about 350 k facts, 26 s of a
-- 1,066 s run. store_legacy fills this table instead, deduping in Python
-- over rows it already holds, and pass-2 reads it back with a bare SELECT
-- and dedupes in Python again. No DISTINCT, no WHERE, no index: sqlite only
-- ever appends to this table and hands the whole of it back.
--
-- Deliberately not keyed. A primary key would dedupe across store_legacy
-- calls for us, but only by making every insert an index probe, and it would
-- make partition.merge fall back to INSERT OR IGNORE. Repeats across calls
-- are cheap to carry and free to drop on read.
CREATE TABLE IF NOT EXISTS hit_seqs (
    base_seq    TEXT NOT NULL,
    database    TEXT NOT NULL
);
"""

# Every index here is built at the end of a run, by finalize_db. None of them
# exists while the pipeline is working. They are for querying the finished
# database afterwards.
#
# _HITS_TABLE_INDEXES used to be built mid-run, after the search phase and
# before classification, on the belief that classification read legacy_hits
# through them. It does not: reconcile_classifications takes the Python
# objects the cascade already returned and never touches the connection.
# On a 135k-sequence run that build cost 44 s and the reconciliation it was
# meant to serve cost 0.03 s.
#
# Worse, the two things that DO read legacy_hits during a run were both made
# slower by having them. `EXPLAIN QUERY PLAN` on this schema:
#
#   SELECT * FROM legacy_hits WHERE database = ?
#       SEARCH legacy_hits USING INDEX idx_leg_db (database=?)
#   SELECT DISTINCT base_seq, database FROM legacy_hits
#       SCAN legacy_hits USING INDEX idx_leg_baseseq
#       USE TEMP B-TREE FOR DISTINCT
#
# `database` holds eight distinct values, so the first plan pulls a third of
# the table through random rowid lookups where a sequential scan would do.
# The second scans an index and then does a rowid lookup per row to fetch
# `database`, turning one sequential pass into a table's worth of random
# reads. Random access is exactly the pattern the parallel filesystem a
# cluster run writes to handles worst.
#
# The trade is not free in one direction: without an index, tesorter_output
# scans the whole table once per database instead of seeking, which costs the
# small databases (tir is 0.6% of rows) what it saves the large ones (rexdb is
# a third). That read totalled 13.6 s across eight databases with the indexes
# present, against 44 s to build them, so the net is clearly positive -- and
# under --memory-db the scans are at RAM speed and the question disappears.
#
# The split between the two lists is kept because they cover different tables:
# classifications and blast_hits are created on demand and may not exist.
_HITS_TABLE_INDEXES = [
    "CREATE INDEX IF NOT EXISTS idx_leg_target  ON legacy_hits(target_name)",
    "CREATE INDEX IF NOT EXISTS idx_leg_baseseq ON legacy_hits(base_seq)",
    "CREATE INDEX IF NOT EXISTS idx_leg_query   ON legacy_hits(query_name)",
    "CREATE INDEX IF NOT EXISTS idx_leg_domain  ON legacy_hits(domain_type)",
    "CREATE INDEX IF NOT EXISTS idx_leg_evalue  ON legacy_hits(i_evalue)",
    "CREATE INDEX IF NOT EXISTS idx_leg_db      ON legacy_hits(database)",
    "CREATE INDEX IF NOT EXISTS idx_leg_engine  ON legacy_hits(engine)",
]

_FINAL_INDEXES = [
    "CREATE INDEX IF NOT EXISTS idx_cls_seq     ON classifications(seq_id)",
    "CREATE INDEX IF NOT EXISTS idx_cls_db      ON classifications(database)",
    "CREATE INDEX IF NOT EXISTS idx_cls_mode    ON classifications(mode)",
    "CREATE INDEX IF NOT EXISTS idx_blast_q     ON blast_hits(qseqid)",
    "CREATE INDEX IF NOT EXISTS idx_blast_s     ON blast_hits(sseqid)",
]


def _index_name(stmt):
    # "CREATE INDEX IF NOT EXISTS <name>  ON <table>(<col>)"
    after_exists = stmt.split("EXISTS", 1)[1].strip()
    return after_exists.split()[0]


_ALL_DEFERRED_INDEX_NAMES = [
    _index_name(s) for s in (_HITS_TABLE_INDEXES + _FINAL_INDEXES)
]

_INSERT_COLS = (
    "database, target_name, base_seq, strand, frame, target_len, "
    "query_name, domain_type, query_len, "
    "evalue, score, bias, dom_num, dom_of, "
    "c_evalue, i_evalue, dom_score, dom_bias, "
    "hmm_from, hmm_to, ali_from, ali_to, "
    "env_from, env_to, acc, search_mode, engine, stage"
)

# Derived from _INSERT_COLS rather than hardcoded: the two must always agree,
# and a stale count fails at insert time with an unhelpful message.
_INSERT_PLACEHOLDERS = ", ".join(["?"] * len(_INSERT_COLS.split(",")))

# Where base_seq sits in the tuples store_legacy inserts. The id column is
# prepended to _INSERT_COLS, so everything shifts by one. Derived for the same
# reason _INSERT_PLACEHOLDERS is.
_BASE_SEQ_POS = 1 + [c.strip() for c in _INSERT_COLS.split(",")].index("base_seq")


def create_db(db_path):
    """Open or create the results database, apply bulk-insert tuning,
    and drop any deferred indexes so subsequent inserts are not slowed
    by per-row B-tree maintenance.

    A reused database (a second pipeline invocation against the same
    .db, e.g. writing a second mode into a companion file) still gets
    the full benefit of deferred indexing: any indexes built by a
    prior finalize_db are dropped here and rebuilt at the end of this
    run by finalize_db.
    """
    return prepare_db(sqlite3.connect(db_path))


def prepare_db(conn):
    """create_db's work on a connection already open: bulk-insert tuning, the
    schema, and the deferred indexes dropped. Separate so a database copied
    into memory (`--extend --memory-db`) gets the same treatment."""
    # WAL trades a tiny durability window (losing only the last transaction
    # on a power cut) for much faster concurrent writes and amortized fsync.
    # synchronous=NORMAL skips the inner fsync between WAL writes and
    # leaves the one at checkpoint time; a rebuildable pipeline output
    # does not need FULL.
    conn.execute("PRAGMA journal_mode = WAL")
    conn.execute("PRAGMA synchronous = NORMAL")
    # 256 MB page cache — keeps the working set in memory during multi-
    # million-row bulk inserts, avoiding constant fault-back from disk.
    conn.execute("PRAGMA cache_size = -262144")
    conn.execute("PRAGMA temp_store = MEMORY")
    conn.executescript(SCHEMA)
    for name in _ALL_DEFERRED_INDEX_NAMES:
        conn.execute(f"DROP INDEX IF EXISTS {name}")
    conn.commit()
    return conn


def backup_to(conn, db_path):
    """Copy an in-memory results database to disk.

    sqlite's Online Backup API, not `VACUUM INTO`: backup copies pages as they
    are, carrying the indexes across already built, and writes them
    sequentially. VACUUM INTO repacks and rebuilds, which is the work being
    avoided.

    Building the database in memory matters because index creation is
    random-write, which is the worst case for a parallel filesystem like the
    scratch a cluster run writes to. Under `--memory-db` the whole results
    database lives in RAM and reaches disk in one sequential pass at the end.
    """
    dest = sqlite3.connect(db_path)
    try:
        conn.commit()
        conn.backup(dest)
        dest.commit()
    finally:
        dest.close()


def index_hits_tables(conn):
    """Build the legacy_hits indexes.

    Called by finalize_db, at the end of a run. Do not call it during one:
    nothing in a pipeline run reads legacy_hits in a way these indexes help,
    and two things read it in ways they actively hurt. The comment on
    _HITS_TABLE_INDEXES has the query plans.
    """
    for stmt in _HITS_TABLE_INDEXES:
        conn.execute(stmt)
    conn.commit()


def finalize_db(conn):
    """Build every index and checkpoint the WAL.

    Call once at the end of the pipeline, after classifications and
    blast_hits are fully populated. Truncating the WAL into the main
    database file leaves the output as a single self-contained .db with
    no sidecar -wal file.

    This is the only place indexes are built. A run works against an
    unindexed legacy_hits on purpose; see the comment on _HITS_TABLE_INDEXES.

    Tables classifications and blast_hits are created on-demand by
    their writer functions; if a pipeline run skipped those phases
    the table may not exist and its index is silently skipped.
    """
    index_hits_tables(conn)
    tables = {r[0] for r in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table'")}
    for stmt in _FINAL_INDEXES:
        target_table = stmt.split(" ON ", 1)[1].split("(", 1)[0].strip()
        if target_table in tables:
            conn.execute(stmt)
    conn.commit()
    conn.execute("PRAGMA wal_checkpoint(TRUNCATE)")


def store_sequences(conn, nucl_lengths):
    """
    Store sequence metadata.

    Args:
        conn: sqlite3 connection
        nucl_lengths: dict of {seq_name: length}
    """
    conn.executemany(
        "INSERT OR REPLACE INTO sequences (name, length) VALUES (?, ?)",
        nucl_lengths.items(),
    )
    conn.commit()


def _parse_frame_info(target_name):
    """Extract base_seq, strand, frame from target name."""
    parts = target_name.rsplit("|", 1)
    if len(parts) != 2:
        return target_name, ".", 0

    base_seq, suffix = parts
    if suffix.startswith("fwd"):
        return base_seq, "+", int(suffix[-1]) - 1
    elif suffix.startswith("rev"):
        return base_seq, "-", int(suffix[-1]) - 1
    else:
        return target_name, ".", 0


def _parse_domain_type(query_name):
    """Extract domain type from model name."""
    if ":" in query_name:
        return query_name.split(":")[1].split("-")[-1]
    return query_name.split("_")[0]


def _hits_to_rows(hits, db_name, search_mode=0, engine="", stage=0,
                  seq_index=None):
    """Convert hit dicts to insert-ready tuples.

    engine/stage stamp cascade provenance onto every row; they default to the
    single-engine case so existing callers are unchanged.
    """
    rows = []
    for h in hits:
        base_seq, strand, frame = _parse_frame_info(h["target_name"])
        domain_type = _parse_domain_type(h["query_name"])

        rows.append((
            -1 if seq_index is None else seq_index.get(base_seq, -1),
            db_name,
            h["target_name"],
            base_seq,
            strand,
            frame,
            h["target_len"],
            h["query_name"],
            domain_type,
            h["query_len"],
            h["evalue"],
            h["score"],
            h["bias"],
            h["dom_num"],
            h["dom_of"],
            h["c_evalue"],
            h["i_evalue"],
            h["dom_score"],
            h["dom_bias"],
            h["hmm_from"],
            h["hmm_to"],
            h["ali_from"],
            h["ali_to"],
            h["env_from"],
            h["env_to"],
            h["acc"],
            search_mode,
            engine,
            stage,
        ))
    return rows


def store_legacy(conn, hits, db_name, engine="", stage=0, seq_index=None):
    """Store search hits to legacy_hits, the only hits table.

    `seq_index` maps a sequence name to its position in the original input;
    that position is written to the `id` column. Without it the column is -1,
    which says "position unknown" rather than inventing an order.
    """
    rows = _hits_to_rows(hits, db_name, search_mode=0, engine=engine,
                         stage=stage, seq_index=seq_index)
    conn.executemany(
        f"INSERT INTO legacy_hits (id, {_INSERT_COLS}) "
        f"VALUES (?, {_INSERT_PLACEHOLDERS})", rows)
    # Summarise as the rows go in. One pass over a list already in memory,
    # and the set collapses a database's hits on a sequence to the single
    # fact pass-2 needs. In a partitioned run this is per worker, so it is
    # 64-way parallel and the parent never reads legacy_hits at all.
    conn.executemany(
        "INSERT INTO hit_seqs (base_seq, database) VALUES (?, ?)",
        {(r[_BASE_SEQ_POS], db_name) for r in rows})
    conn.commit()


def export_tsv(conn, tsv_path, table="pass2_hits", db_name=None):
    """
    Export raw domain hits to a tab-separated file.

    Args:
        conn: sqlite3 connection
        tsv_path: output file path
        table: "pass1_hits" or "pass2_hits"
        db_name: if set, filter to this database only
    """
    assert table in ("pass1_hits", "pass2_hits", "legacy_hits")

    columns = [
        "database", "target_name", "target_len", "query_name", "query_len",
        "evalue", "score", "bias", "dom_num", "dom_of",
        "c_evalue", "i_evalue", "dom_score", "dom_bias",
        "hmm_from", "hmm_to", "ali_from", "ali_to",
        "env_from", "env_to", "acc",
    ]

    query = f"SELECT {', '.join(columns)} FROM {table}"
    params = ()
    if db_name is not None:
        query += " WHERE database = ?"
        params = (db_name,)
    query += " ORDER BY database, target_name, query_name, dom_num"

    cursor = conn.execute(query, params)

    with open(tsv_path, "w") as f:
        f.write("\t".join(columns) + "\n")
        for row in cursor:
            f.write("\t".join(str(v) for v in row) + "\n")


def export_best_hits_tsv(conn, tsv_path, nucl_lengths=None,
                         table="pass2_hits", db_name=None):
    """
    Export best domain hits as a clean tab-separated file.

    One row per (sequence, model) pair -- the highest-scoring domain.
    For translated sequences, coordinates are converted to nucleotide
    positions. Columns are biologist-friendly.

    Args:
        conn: sqlite3 connection
        tsv_path: output file path
        nucl_lengths: dict of {seq_name: nucl_length} for coordinate
                      conversion. If None, AA coordinates are reported.
        table: "pass1_hits" or "pass2_hits"
        db_name: if set, filter to this database only
    """
    assert table in ("pass1_hits", "pass2_hits", "legacy_hits")

    db_filter = ""
    params = ()
    if db_name is not None:
        db_filter = "AND d1.database = ? AND d2.database = ?"
        params = (db_name, db_name)

    rows = conn.execute(f"""
        SELECT d1.database, d1.target_name, d1.target_len,
               d1.query_name, d1.query_len,
               d1.i_evalue, d1.dom_score, d1.dom_bias, d1.acc,
               d1.hmm_from, d1.hmm_to,
               d1.env_from, d1.env_to
        FROM {table} d1
        WHERE d1.dom_score = (
            SELECT MAX(d2.dom_score) FROM {table} d2
            WHERE d2.target_name = d1.target_name
            AND d2.query_name = d1.query_name
            {db_filter}
        )
        {f"AND d1.database = ?" if db_name else ""}
        ORDER BY d1.database, d1.target_name, d1.dom_score DESC
    """, params + ((db_name,) if db_name else ()))

    columns = [
        "database", "seq_id", "seq_len", "model", "model_len",
        "strand", "frame",
        "nuc_start", "nuc_end", "env_from_aa", "env_to_aa",
        "hmm_from", "hmm_to", "hmm_coverage",
        "evalue", "score", "bias", "accuracy",
    ]

    with open(tsv_path, "w") as f:
        f.write("\t".join(columns) + "\n")

        for row in rows:
            (db, target, target_len, model, model_len,
             evalue, score, bias, acc,
             hmm_from, hmm_to, env_from, env_to) = row

            seq_id, strand, frame = parse_frame_suffix(target)

            hmm_cov = round(100.0 * (hmm_to - hmm_from + 1) / model_len, 1)

            # Coordinate conversion
            if strand in ("+", "-") and nucl_lengths and seq_id in nucl_lengths:
                nuc_start, nuc_end = aa_to_nucl_coords(
                    env_from, env_to, strand, frame, nucl_lengths[seq_id])
                seq_len = nucl_lengths[seq_id]
            else:
                nuc_start, nuc_end = env_from, env_to
                seq_len = target_len
                if strand == ".":
                    frame = "."

            vals = [
                db, seq_id, seq_len, model, model_len,
                strand, frame,
                nuc_start, nuc_end, env_from, env_to,
                hmm_from, hmm_to, hmm_cov,
                f"{evalue:.2e}", f"{score:.1f}", f"{bias:.1f}", f"{acc:.3f}",
            ]
            f.write("\t".join(str(v) for v in vals) + "\n")


def export_all_domains_tsv(conn, tsv_path, nucl_lengths=None,
                           table="pass2_hits", db_name=None):
    """
    Export all domain hits as a clean tab-separated file.

    Every domain occurrence, not just the best per pair. Useful for
    seeing multi-domain architectures within a single sequence.

    Args:
        conn: sqlite3 connection
        tsv_path: output file path
        nucl_lengths: dict for coordinate conversion (optional)
        table: "pass1_hits" or "pass2_hits"
        db_name: if set, filter to this database only
    """
    assert table in ("pass1_hits", "pass2_hits", "legacy_hits")

    query = f"""
        SELECT database, target_name, target_len,
               query_name, query_len,
               dom_num, dom_of,
               i_evalue, dom_score, dom_bias, acc,
               hmm_from, hmm_to,
               env_from, env_to
        FROM {table}
    """
    params = ()
    if db_name is not None:
        query += " WHERE database = ?"
        params = (db_name,)
    query += " ORDER BY database, target_name, query_name, dom_num"

    rows = conn.execute(query, params)

    columns = [
        "database", "seq_id", "seq_len", "model", "model_len",
        "strand", "frame",
        "dom_num", "dom_of",
        "nuc_start", "nuc_end", "env_from_aa", "env_to_aa",
        "hmm_from", "hmm_to", "hmm_coverage",
        "evalue", "score", "bias", "accuracy",
    ]

    with open(tsv_path, "w") as f:
        f.write("\t".join(columns) + "\n")

        for row in rows:
            (db, target, target_len, model, model_len,
             dom_num, dom_of,
             evalue, score, bias, acc,
             hmm_from, hmm_to, env_from, env_to) = row

            seq_id, strand, frame = parse_frame_suffix(target)

            hmm_cov = round(100.0 * (hmm_to - hmm_from + 1) / model_len, 1)

            if strand in ("+", "-") and nucl_lengths and seq_id in nucl_lengths:
                nuc_start, nuc_end = aa_to_nucl_coords(
                    env_from, env_to, strand, frame, nucl_lengths[seq_id])
                seq_len = nucl_lengths[seq_id]
            else:
                nuc_start, nuc_end = env_from, env_to
                seq_len = target_len
                if strand == ".":
                    frame = "."

            vals = [
                db, seq_id, seq_len, model, model_len,
                strand, frame,
                dom_num, dom_of,
                nuc_start, nuc_end, env_from, env_to,
                hmm_from, hmm_to, hmm_cov,
                f"{evalue:.2e}", f"{score:.1f}", f"{bias:.1f}", f"{acc:.3f}",
            ]
            f.write("\t".join(str(v) for v in vals) + "\n")


def export_domain_sequences(conn, fasta_path, aa_fasta, nucl_lengths=None,
                            table="pass2_hits", db_name=None):
    """
    Export domain protein subsequences as FASTA.

    Extracts the amino acid subsequence for each best domain hit.
    Loads the translated FASTA into memory once for fast slicing.

    Args:
        conn: sqlite3 connection
        fasta_path: output FASTA path
        aa_fasta: path to the translated amino acid FASTA
        nucl_lengths: dict for nucleotide coordinate annotation (optional)
        table: "pass1_hits" or "pass2_hits"
        db_name: if set, filter to this database only
    """
    assert table in ("pass1_hits", "pass2_hits", "legacy_hits")

    seqs = load_sequences_dict(aa_fasta)

    db_filter = ""
    params = ()
    if db_name is not None:
        db_filter = "AND d1.database = ? AND d2.database = ?"
        params = (db_name, db_name)

    rows = conn.execute(f"""
        SELECT d1.target_name, d1.query_name, d1.query_len,
               d1.env_from, d1.env_to, d1.dom_score, d1.i_evalue, d1.acc,
               d1.hmm_from, d1.hmm_to
        FROM {table} d1
        WHERE d1.dom_score = (
            SELECT MAX(d2.dom_score) FROM {table} d2
            WHERE d2.target_name = d1.target_name
            AND d2.query_name = d1.query_name
            {db_filter}
        )
        {f"AND d1.database = ?" if db_name else ""}
        ORDER BY d1.target_name, d1.dom_score DESC
    """, params + ((db_name,) if db_name else ()))

    with open(fasta_path, "w") as f:
        for row in rows:
            (target, model, model_len,
             env_from, env_to, score, evalue, acc,
             hmm_from, hmm_to) = row

            if target not in seqs:
                continue

            subseq = seqs[target][env_from - 1:env_to]

            seq_id, strand, frame = parse_frame_suffix(target)
            hmm_cov = round(100.0 * (hmm_to - hmm_from + 1) / model_len, 1)

            # Build nucleotide coordinate annotation if available
            if strand in ("+", "-") and nucl_lengths and seq_id in nucl_lengths:
                nuc_start, nuc_end = aa_to_nucl_coords(
                    env_from, env_to, strand, frame, nucl_lengths[seq_id])
                coord_str = f"nuc={nuc_start}-{nuc_end};"
            else:
                coord_str = ""

            header = (f"{seq_id}|{model} "
                      f"{coord_str}"
                      f"env={env_from}-{env_to};"
                      f"score={score:.1f};"
                      f"evalue={evalue:.2e};"
                      f"hmm_cov={hmm_cov};"
                      f"acc={acc:.3f}")

            f.write(f">{header}\n{subseq}\n")


def query_best_hits(conn, table="pass2_hits", db_name=None):
    """
    Get the best domain hit per (target, query) pair by dom_score.

    Returns:
        list of sqlite3.Row objects
    """
    assert table in ("pass1_hits", "pass2_hits", "legacy_hits")
    conn.row_factory = sqlite3.Row

    where = ""
    params = ()
    if db_name is not None:
        where = "AND d1.database = ? AND d2.database = ?"
        params = (db_name, db_name)

    cursor = conn.execute(f"""
        SELECT * FROM {table} d1
        WHERE dom_score = (
            SELECT MAX(dom_score) FROM {table} d2
            WHERE d2.target_name = d1.target_name
            AND d2.query_name = d1.query_name
            {where}
        )
        {f"AND d1.database = ?" if db_name else ""}
        ORDER BY target_name, dom_score DESC
    """, params + ((db_name,) if db_name else ()))
    return cursor.fetchall()
