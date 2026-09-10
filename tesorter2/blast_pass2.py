"""
blast_pass2.py — BLAST-based pass-2 classification for HMM-unclassified sequences.

Sequences not classified by HMM search are BLASTed against classified sequences.
If a strong similarity match exists (80-80-80 rule by default), the unclassified
sequence inherits the classification of its best BLAST hit.

Key design:
  - Cross-database: one BLAST search against all classified sequences from all databases
  - Per-database reconstruction via filtering on classified_by
  - Parallel chunked BLAST with greedy bin-packing by sequence length
  - All results stored in SQLite for post-hoc threshold adjustment
"""

import logging
import os
import subprocess
import tempfile
import time
from collections import defaultdict

import pyfastx

log = logging.getLogger(__name__)


def _get_classified_ids(conn):
    """Get classified sequence IDs per database from classifier results.

    Returns:
        classified: {base_seq: set(databases)} — which databases classified each seq
    """
    classified = defaultdict(set)

    # Check which tables exist
    tables = {r[0] for r in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table'").fetchall()}

    if "legacy_hits" in tables:
        for row in conn.execute(
            "SELECT DISTINCT base_seq, database FROM legacy_hits"
        ):
            classified[row[0]].add(row[1])

    return dict(classified)


def split_classified_unclassified(input_fasta, classified_ids, outdir,
                                  n_chunks=4, seq_type="nucl"):
    """Split input into classified (BLAST db) and chunked unclassified (queries).

    One pass through the input. Unclassified sequences are bin-packed into
    n_chunks files by total sequence length for even BLAST parallelism.

    Args:
        input_fasta: path to input FASTA
        classified_ids: {seq_name: set(databases)}
        outdir: directory for output files
        n_chunks: number of query chunks
        seq_type: "nucl" or "prot" (determines which sequences to write)

    Returns:
        db_fasta: path to classified sequences FASTA (BLAST database)
        query_chunks: list of paths to unclassified sequence chunks
        db_seq_to_dbs: {seq_name: set(databases)} for DB sequences
    """
    os.makedirs(outdir, exist_ok=True)

    db_fasta = os.path.join(outdir, "blast_db.fa")
    chunk_paths = [os.path.join(outdir, f"blast_query_{i}.fa") for i in range(n_chunks)]

    # Open all handles
    db_handle = open(db_fasta, "w")
    chunk_handles = [open(p, "w") for p in chunk_paths]
    chunk_lengths = [0] * n_chunks

    # build_index=False: this walks the whole input once, in order, and never
    # looks a sequence up by name. An index here is a SQLite sidecar built over
    # the entire 388 Mb input for nothing -- the largest of the wasted index
    # builds, since every other one was over a stage subset.
    fa = pyfastx.Fasta(input_fasta, build_index=False)
    n_classified = 0
    n_unclassified = 0
    db_seq_to_dbs = {}

    for name, seq in fa:

        if name in classified_ids:
            db_handle.write(f">{name}\n{seq}\n")
            db_seq_to_dbs[name] = classified_ids[name]
            n_classified += 1
        else:
            # Bin-pack into lightest chunk
            min_idx = chunk_lengths.index(min(chunk_lengths))
            chunk_handles[min_idx].write(f">{name}\n{seq}\n")
            chunk_lengths[min_idx] += len(seq)
            n_unclassified += 1

    db_handle.close()
    for h in chunk_handles:
        h.close()

    # Remove empty chunks
    query_chunks = [p for p in chunk_paths if os.path.getsize(p) > 0]

    log.info(f"  Split: {n_classified} classified (DB), "
             f"{n_unclassified} unclassified -> {len(query_chunks)} chunks")

    return db_fasta, query_chunks, db_seq_to_dbs


def make_blast_db(db_fasta, seq_type="nucl"):
    """Run makeblastdb."""
    dbtype = seq_type
    cmd = f"makeblastdb -in {db_fasta} -dbtype {dbtype} -out {db_fasta}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        log.error(f"makeblastdb failed: {result.stderr}")
        raise RuntimeError(f"makeblastdb failed: {result.stderr}")
    log.info(f"  BLAST database built: {db_fasta}")


def run_blast_chunk(query_chunk, db_fasta, output, seq_type="nucl", ncpu=1):
    """Run BLAST on one query chunk."""
    app = "blastn" if seq_type == "nucl" else "blastp"
    outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs qcovhsp sstrand"
    cmd = (f"{app} -query {query_chunk} -db {db_fasta} -out {output} "
           f"-outfmt '{outfmt}' -num_threads {ncpu}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        log.warning(f"BLAST chunk failed: {result.stderr[:200]}")
    return output


def parse_blast_output(blast_out):
    """Parse BLAST outfmt 6 into hit dicts."""
    fields = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
              "qstart", "qend", "sstart", "send", "evalue", "bitscore",
              "qlen", "slen", "qcovs", "qcovhsp", "sstrand"]
    types = [str, str, float, int, int, int, int, int, int, int,
             float, float, int, int, float, float, str]

    hits = []
    if not os.path.exists(blast_out):
        return hits

    with open(blast_out) as f:
        for line in f:
            vals = line.strip().split("\t")
            if len(vals) < len(fields):
                continue
            hit = {}
            for field, typ, val in zip(fields, types, vals):
                hit[field] = typ(val)
            hits.append(hit)

    return hits


def store_blast_hits(conn, hits, db_seq_to_dbs):
    """Store BLAST hits in SQLite.

    Adds classified_by field indicating which databases classified
    each target sequence.
    """
    # Indexes on blast_hits are built by results.finalize_db at end of run.
    conn.execute("""
        CREATE TABLE IF NOT EXISTS blast_hits (
            qseqid      TEXT NOT NULL,
            sseqid      TEXT NOT NULL,
            pident      REAL NOT NULL,
            length      INTEGER NOT NULL,
            evalue      REAL NOT NULL,
            bitscore    REAL NOT NULL,
            qlen        INTEGER NOT NULL,
            slen        INTEGER NOT NULL,
            qcovs       REAL NOT NULL,
            classified_by TEXT NOT NULL
        )
    """)

    rows = []
    for h in hits:
        dbs = db_seq_to_dbs.get(h["sseqid"], set())
        classified_by = ",".join(sorted(dbs)) if dbs else "unknown"
        rows.append((
            h["qseqid"], h["sseqid"], h["pident"], h["length"],
            h["evalue"], h["bitscore"], h["qlen"], h["slen"],
            h["qcovs"], classified_by,
        ))

    conn.executemany(
        "INSERT INTO blast_hits VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        rows,
    )
    conn.commit()
    log.info(f"  Stored {len(rows)} BLAST hits")


def classify_from_blast(conn, classifications, database=None,
                        min_identity=80, min_coverage=80, min_length=80):
    """Classify unclassified sequences from BLAST hits.

    Args:
        conn: sqlite3 connection with blast_hits table
        classifications: dict of {seq_id: {order, superfamily, ...}} from
                         classifier.classify_sequences()
        database: if set, only accept BLAST targets classified by this database
        min_identity: minimum percent identity
        min_coverage: minimum query coverage
        min_length: minimum alignment length

    Returns:
        list of classification dicts for newly classified sequences
    """
    # Check table exists
    tables = {r[0] for r in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table'").fetchall()}
    if "blast_hits" not in tables:
        return []

    # Build filter
    where = "WHERE pident >= ? AND qcovs >= ? AND length >= ?"
    params = [min_identity, min_coverage, min_length]

    if database:
        where += " AND classified_by LIKE ?"
        params.append(f"%{database}%")

    # Best hit per query by bitscore
    rows = conn.execute(f"""
        SELECT qseqid, sseqid, pident, qcovs, length, bitscore
        FROM blast_hits
        {where}
        ORDER BY bitscore DESC
    """, params).fetchall()

    # Classified ID set for quick lookup
    classified_set = set(classifications.keys())

    best = {}
    for qid, sid, pident, qcovs, length, bitscore in rows:
        if qid in classified_set:
            continue  # already classified by HMM, skip
        if qid not in best:
            best[qid] = (sid, pident, qcovs, length, bitscore)

    # Inherit classification from best hit's target
    new_classifications = []
    no_source = 0
    for qid, (sid, pident, qcovs, length, bitscore) in best.items():
        if sid in classifications:
            source = classifications[sid]
            new_classifications.append({
                "id": qid,
                "order": source["order"],
                "superfamily": source["superfamily"],
                "clade": "unknown",
                "complete": "none",
                "strand": "?",
                "domains": "none",
                "blast_source": sid,
                "blast_pident": pident,
                "blast_qcovs": qcovs,
                "blast_bitscore": bitscore,
            })
        else:
            no_source += 1

    if no_source:
        log.info(f"    {no_source} BLAST hits to unclassified targets (skipped)")

    log.info(f"  BLAST pass-2: {len(new_classifications)} sequences classified "
             f"(from {len(best)} hits passing filters)")
    return new_classifications


def blast_pass2(input_fasta, conn, hmm_classifications=None,
                seq_type="nucl", n_processors=4,
                min_identity=80, min_coverage=80, min_length=80,
                outdir=None):
    """Full BLAST pass-2 pipeline.

    Args:
        input_fasta: path to input FASTA
        conn: sqlite3 connection with HMM results
        hmm_classifications: dict of {seq_id: {order, superfamily, clade, ...}}
                             from classifier.classify_sequences(). Targets inherit
                             classification from their best BLAST match.
        seq_type: "nucl" or "prot"
        n_processors: number of parallel BLAST processes
        min_identity: filter threshold
        min_coverage: filter threshold
        min_length: filter threshold
        outdir: output directory (default: tempdir)

    Returns:
        list of new classification dicts
    """
    t0 = time.time()

    # Get classified IDs from database
    classified_ids = _get_classified_ids(conn)
    if not classified_ids:
        log.info("  No classified sequences for BLAST pass-2")
        return []

    log.info(f"  BLAST pass-2: {len(classified_ids)} classified sequences as targets")

    # Split
    if outdir is None:
        outdir = tempfile.mkdtemp(prefix="tesorter2_blast_")

    blast_dir = os.path.join(outdir, "blast_pass2")
    db_fasta, query_chunks, db_seq_to_dbs = split_classified_unclassified(
        input_fasta, classified_ids, blast_dir, n_chunks=n_processors,
        seq_type=seq_type)

    if not query_chunks:
        log.info("  No unclassified sequences to search")
        return []

    # Build BLAST database
    make_blast_db(db_fasta, seq_type=seq_type)

    # Run parallel BLAST
    log.info(f"  Running {len(query_chunks)} BLAST processes")
    t1 = time.time()

    import multiprocessing
    blast_outputs = []
    args_list = []
    for i, chunk in enumerate(query_chunks):
        out = chunk + ".blastout"
        blast_outputs.append(out)
        args_list.append((chunk, db_fasta, out, seq_type, 1))

    with multiprocessing.Pool(n_processors) as pool:
        pool.starmap(run_blast_chunk, args_list)

    t2 = time.time()
    log.info(f"  BLAST search: {t2 - t1:.1f}s")

    # Parse and store
    all_hits = []
    for blast_out in blast_outputs:
        all_hits.extend(parse_blast_output(blast_out))

    log.info(f"  {len(all_hits)} total BLAST hits")

    if all_hits:
        store_blast_hits(conn, all_hits, db_seq_to_dbs)

    # Build HMM classifications lookup for inheritance
    hmm_cls = {}
    if hmm_classifications is not None:
        hmm_cls = hmm_classifications
    log.info(f"  {len(hmm_cls)} HMM classifications available for inheritance")

    # Classify
    new_cls = classify_from_blast(
        conn, hmm_cls,
        min_identity=min_identity,
        min_coverage=min_coverage,
        min_length=min_length,
    )

    t3 = time.time()
    log.info(f"  BLAST pass-2 total: {t3 - t0:.1f}s")

    return new_cls
