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
"""

# Index creation is split by pipeline phase.
#
# _HITS_TABLE_INDEXES: built once the search phase finishes writing all
#   HMM hits and before the deconfliction/classification phase begins.
#   Those phases read from legacy_hits via filters on
#   database / base_seq / domain_type, which benefit from an index.
#
# _FINAL_INDEXES: built at the very end of the pipeline, after
#   classifications and blast_hits are fully populated. Neither table
#   is read during pipeline execution; these indexes exist for post-run
#   interactive analysis.
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


def create_db(db_path):
    """Open or create the results database, apply bulk-insert tuning,
    and drop any deferred indexes so subsequent inserts are not slowed
    by per-row B-tree maintenance.

    A reused database (a second pipeline invocation against the same
    .db, e.g. writing a second mode into a companion file) still gets
    the full benefit of deferred indexing: any indexes built by a
    prior finalize_db are dropped here and rebuilt at the end of this
    run via index_hits_tables + finalize_db.
    """
    conn = sqlite3.connect(db_path)
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


def index_hits_tables(conn):
    """Build indexes on legacy_hits.

    Call after the HMM search phase is complete (all hits written) and
    before the classification/deconfliction phase begins. The subsequent
    reads filter by database / base_seq / domain_type and benefit from
    these indexes; inserts during the search phase do not.
    """
    for stmt in _HITS_TABLE_INDEXES:
        conn.execute(stmt)
    conn.commit()


def finalize_db(conn):
    """Build the remaining indexes and checkpoint the WAL.

    Call once at the end of the pipeline, after classifications and
    blast_hits are fully populated. Truncating the WAL into the main
    database file leaves the output as a single self-contained .db with
    no sidecar -wal file.

    Tables classifications and blast_hits are created on-demand by
    their writer functions; if a pipeline run skipped those phases
    the table may not exist and its index is silently skipped.
    """
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
