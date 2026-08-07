"""
blast_backend.py — TEsorter2 master's blastn pass-2 logic, lifted verbatim so
the `--pass2-aligner blast` path runs and post-processes identically to
TEsorter2 master.

The orchestrator `run_pass2_blast` is called by blast_pass2.blast_pass2() AFTER
the shared classified/unclassified split and the optional pass2_external merge,
so the external-pool augmentation (`--pass2-classified-fasta`) applies to both
aligner backends. Only the alignment + parse + classify is master-specific here.

master functions (make_blast_db, run_blast_chunk, parse_blast_output,
store_blast_hits, classify_from_blast) are copied from
origin/master:src/blast_pass2.py; run_blast_chunk additionally takes a
blast_task parameter (--blast-task). chunk_fasta + run_pass2_blast are new.
"""

import logging
import multiprocessing
import os
import subprocess
import time

import pyfastx

log = logging.getLogger(__name__)


# ---- lifted verbatim from origin/master:src/blast_pass2.py ----

def make_blast_db(db_fasta, seq_type="nucl"):
    """Run makeblastdb."""
    dbtype = seq_type
    cmd = f"makeblastdb -in {db_fasta} -dbtype {dbtype} -out {db_fasta}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        log.error(f"makeblastdb failed: {result.stderr}")
        raise RuntimeError(f"makeblastdb failed: {result.stderr}")
    log.info(f"  BLAST database built: {db_fasta}")


def run_blast_chunk(query_chunk, db_fasta, output, seq_type="nucl", ncpu=1,
                    blast_task="megablast"):
    """Run BLAST on one query chunk."""
    app = "blastn" if seq_type == "nucl" else "blastp"
    outfmt = ("6 qseqid sseqid pident length mismatch gapopen qstart qend "
              "sstart send evalue bitscore qlen slen qcovs qcovhsp sstrand")
    # The -task value is blastn-only; blastp would reject it, so gate it on
    # the blastn branch.
    task = f" -task {blast_task}" if app == "blastn" else ""
    cmd = (f"{app}{task} -query {query_chunk} -db {db_fasta} -out {output} "
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
    """Store BLAST hits in SQLite."""
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
    """Classify unclassified sequences from BLAST hits (master logic)."""
    tables = {r[0] for r in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table'").fetchall()}
    if "blast_hits" not in tables:
        return []

    where = "WHERE pident >= ? AND qcovs >= ? AND length >= ?"
    params = [min_identity, min_coverage, min_length]

    if database:
        where += " AND classified_by LIKE ?"
        params.append(f"%{database}%")

    rows = conn.execute(f"""
        SELECT qseqid, sseqid, pident, qcovs, length, bitscore
        FROM blast_hits
        {where}
        ORDER BY bitscore DESC
    """, params).fetchall()

    classified_set = set(classifications.keys())

    best = {}
    for qid, sid, pident, qcovs, length, bitscore in rows:
        if qid in classified_set:
            continue
        if qid not in best:
            best[qid] = (sid, pident, qcovs, length, bitscore)

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


# ---- new: chunking + orchestration ----

def chunk_fasta(qry_fasta, n_chunks, outdir):
    """Bin-pack sequences from qry_fasta into n_chunks files by total length.

    Mirrors master's split bin-packing, but operates on an already-written
    unclassified-query FASTA (the shared split in blast_pass2 produced it).
    Returns the list of non-empty chunk paths. Returns [] for a missing or
    empty input (pyfastx raises on empty files).
    """
    os.makedirs(outdir, exist_ok=True)
    if not os.path.exists(qry_fasta) or os.path.getsize(qry_fasta) == 0:
        return []
    chunk_paths = [os.path.join(outdir, f"blast_query_{i}.fa")
                   for i in range(max(1, n_chunks))]
    handles = []
    try:
        handles = [open(p, "w") for p in chunk_paths]
        lengths = [0] * len(chunk_paths)
        fa = pyfastx.Fasta(qry_fasta, build_index=True)
        for rec in fa:
            i = lengths.index(min(lengths))
            handles[i].write(f">{rec.name}\n{rec.seq}\n")
            lengths[i] += len(rec.seq)
    finally:
        for h in handles:
            h.close()
    return [p for p in chunk_paths if os.path.getsize(p) > 0]


def run_pass2_blast(qry_fasta, db_fasta, conn, classifications, db_seq_to_dbs,
                    n_processors, min_identity, min_coverage, min_length, work,
                    blast_task="megablast"):
    """blastn pass-2 over an already-prepared (db_fasta, qry_fasta) pair.

    Reproduces TEsorter2 master's: makeblastdb -> chunked parallel blastn ->
    outfmt6 parse -> SQLite -> classify_from_blast (qcovs+length filter, best by
    bitscore, clade=unknown). The I/C/L thresholds come from the run's -rule.
    """
    t0 = time.time()
    os.makedirs(work, exist_ok=True)

    if not os.path.exists(db_fasta) or os.path.getsize(db_fasta) == 0:
        log.info("  pass-2 target FASTA is empty or missing; skipping blastn")
        return []

    make_blast_db(db_fasta, seq_type="nucl")

    query_chunks = chunk_fasta(qry_fasta, n_processors, work)
    if not query_chunks:
        log.info("  No unclassified sequences to search")
        return []

    log.info(f"  Running {len(query_chunks)} BLAST processes")
    t1 = time.time()
    blast_outputs = []
    args_list = []
    for chunk in query_chunks:
        out = chunk + ".blastout"
        blast_outputs.append(out)
        args_list.append((chunk, db_fasta, out, "nucl", 1, blast_task))
    with multiprocessing.Pool(len(query_chunks)) as pool:
        pool.starmap(run_blast_chunk, args_list)
    t2 = time.time()
    log.info(f"  BLAST search: {t2 - t1:.1f}s")

    all_hits = []
    for blast_out in blast_outputs:
        all_hits.extend(parse_blast_output(blast_out))
    log.info(f"  {len(all_hits)} total BLAST hits")

    if all_hits:
        store_blast_hits(conn, all_hits, db_seq_to_dbs)

    new_cls = classify_from_blast(
        conn, classifications,
        min_identity=min_identity,
        min_coverage=min_coverage,
        min_length=min_length,
    )
    log.info(f"  BLAST pass-2 total: {time.time() - t0:.1f}s")
    return new_cls
