"""
blast_pass2.py — minimap2-based pass-2 classification for HMM-unclassified
sequences.

Pipeline:
  1. minimap2 (sensitivity-tuned flags from `minimap.run_minimap2`) writes a
     PAF for the unclassified queries against the classified-pool target.
  2. `classify_ltr_paf_fast.process_paf` reduces the PAF to one row per query:
        qname  pass/fail  pid  eff_qcov  eff_tcov  best_tname
     under the rule `--min-pid I  --min-qcov C  --min-tcov C` derived from
     `--pass2-rule I-C-L`. (Per benchmarking, 70-70-70 is recommended for
     this minimap2 path; the CLI default is 80-80-80 to mirror upstream's
     blastn pass-2. The L value is parsed for backwards-compat with the
     I-C-L grammar but is not consumed downstream.)
  3. Each `pass` row inherits the target's order/superfamily/clade.

The SQLite `blast_hits` table is preserved for post-run introspection but the
schema has been simplified to the columns classify_ltr_paf_fast emits.
"""

import logging
import os
import tempfile
import time
from collections import defaultdict

import pyfastx

from . import minimap
from . import pass2_external
from .classify_ltr_paf_fast import process_paf, format_row

log = logging.getLogger(__name__)


def _get_classified_ids(conn):
    """Get classified sequence IDs per database from classifier results."""
    classified = defaultdict(set)
    tables = {r[0] for r in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table'").fetchall()}
    if "legacy_hits" in tables:
        for row in conn.execute(
            "SELECT DISTINCT base_seq, database FROM legacy_hits"
        ):
            classified[row[0]].add(row[1])
    return dict(classified)


def split_classified_unclassified(input_fasta, classified_ids, outdir,
                                  seq_type="nucl"):
    """Split input into classified-pool FASTA (target) and unclassified-query
    FASTA. Nucleotide pools are uppercased and stripped of non-ATCG to keep
    inputs clean (minimap2 itself tolerates ambiguous bases).
    """
    os.makedirs(outdir, exist_ok=True)
    db_fasta = os.path.join(outdir, "blast_db.fa")
    qry_fasta = os.path.join(outdir, "blast_query.fa")

    nucl = (seq_type == "nucl")
    fa = pyfastx.Fasta(input_fasta, build_index=True)

    n_classified = 0
    n_unclassified = 0
    db_seq_to_dbs = {}

    with open(db_fasta, "w") as dbh, open(qry_fasta, "w") as qh:
        for rec in fa:
            name = rec.name
            seq = str(rec.seq)
            if nucl:
                seq = "".join(c for c in seq.upper() if c in "ATCG")
            if name in classified_ids:
                dbh.write(f">{name}\n{seq}\n")
                db_seq_to_dbs[name] = classified_ids[name]
                n_classified += 1
            else:
                qh.write(f">{name}\n{seq}\n")
                n_unclassified += 1

    log.info(f"  Split: {n_classified} classified (DB), "
             f"{n_unclassified} unclassified (query)")
    return db_fasta, qry_fasta, db_seq_to_dbs


def run_alignment(query_fa, target_fa, paf_out, ncpu=4,
                  preset="asm20", extra=""):
    """Run minimap2 with the sensitivity-tuned pass-2 flag set."""
    minimap.run_minimap2(
        query_fa=query_fa, target_fa=target_fa, paf_out=paf_out,
        ncpu=ncpu, preset=preset, extra=extra,
    )
    return paf_out


def classify_paf_to_tsv(paf_path, tsv_out, min_pid, min_qcov, min_tcov):
    """Run classify_ltr_paf_fast over the PAF, write the TSV alongside,
    and return the parsed rows: [(qname, pass_str, pid, qcov, tcov, best_tname), ...].
    """
    if not os.path.exists(paf_path) or os.path.getsize(paf_path) == 0:
        return []

    with open(paf_path) as fh:
        rows = process_paf(
            fh, min_pid=min_pid, min_qcov=min_qcov, min_tcov=min_tcov,
            fill_colinear=False,
            verbose=False,
        )

    with open(tsv_out, "w") as fh:
        for row in rows:
            fh.write(format_row(row) + "\n")

    n_pass = sum(1 for r in rows if r[1] == "pass")
    log.info(f"  classify_ltr_paf_fast: {n_pass}/{len(rows)} queries pass "
             f"(min_pid={min_pid:.3f} min_qcov={min_qcov:.3f} "
             f"min_tcov={min_tcov:.3f}) -> {tsv_out}")
    return rows


def store_blast_hits(conn, tsv_rows, db_seq_to_dbs):
    """Store classify_ltr_paf_fast rows in SQLite. One row per query."""
    conn.execute("""
        CREATE TABLE IF NOT EXISTS blast_hits (
            qseqid        TEXT NOT NULL,
            sseqid        TEXT NOT NULL,
            pident        REAL NOT NULL,
            qcovs         REAL NOT NULL,
            tcovs         REAL NOT NULL,
            passes_rule   INTEGER NOT NULL,
            classified_by TEXT NOT NULL
        )
    """)
    rows = []
    for qname, pass_str, pid, qcov, tcov, tname in tsv_rows:
        dbs = db_seq_to_dbs.get(tname, set())
        classified_by = ",".join(sorted(dbs)) if dbs else "unknown"
        rows.append((
            qname, tname,
            pid * 100.0, qcov * 100.0, tcov * 100.0,
            1 if pass_str == "pass" else 0,
            classified_by,
        ))
    conn.executemany(
        "INSERT INTO blast_hits VALUES (?, ?, ?, ?, ?, ?, ?)", rows,
    )
    conn.commit()
    log.info(f"  Stored {len(rows)} minimap2 best-hit rows")


def classify_from_blast(tsv_rows, classifications):
    """Inherit order/superfamily/clade for queries whose best target passed
    the rule. `classifications` maps target_id -> classification dict.
    """
    classified_set = set(classifications.keys())
    new_classifications = []
    no_source = 0
    n_pass = 0
    for qname, pass_str, pid, qcov, tcov, tname in tsv_rows:
        if pass_str != "pass":
            continue
        n_pass += 1
        if qname in classified_set:
            continue
        if tname not in classifications:
            no_source += 1
            continue
        source = classifications[tname]
        new_classifications.append({
            "id": qname,
            "order": source["order"],
            "superfamily": source["superfamily"],
            "clade": source.get("clade", "unknown"),
            "complete": "none",
            "strand": "?",
            "domains": "none",
            "blast_source": tname,
            "blast_pident": pid * 100.0,
            "blast_qcovs": qcov * 100.0,
            "blast_tcovs": tcov * 100.0,
            "blast_bitscore": 0.0,
        })

    if no_source:
        log.info(f"    {no_source} pass-2 hits to unclassified targets (skipped)")
    log.info(f"  pass-2: {len(new_classifications)} sequences classified "
             f"(from {n_pass} rule-passing queries)")
    return new_classifications


def blast_pass2(input_fasta, conn, hmm_classifications=None,
                seq_type="nucl", n_processors=4,
                min_identity=80, min_coverage=80, min_length=80,
                outdir=None,
                pass2_classified_fasta=None,
                preset="asm20", minimap2_extra="",
                aligner="blast", blast_task="megablast"):
    """Pass-2 similarity search (blastn by default, minimap2 opt-in).

    Args:
      min_identity: I from --pass2-rule I-C-L (percent, e.g. 80)
      min_coverage: C from --pass2-rule I-C-L (percent; blast applies it to
                    qcovs, minimap2 to qcov AND tcov)
      min_length:   L from --pass2-rule I-C-L (blast: minimum alignment
                    length; minimap2: parsed for backwards-compat with the
                    I-C-L grammar, not consumed by classify_ltr_paf_fast)
    """
    t0 = time.time()
    if aligner == "minimap2":
        minimap.check_minimap2()

    if seq_type != "nucl":
        log.warning("minimap2 pass-2 only supports nucleotide sequences; "
                    f"seq_type={seq_type!r} will be treated as nucl")
        seq_type = "nucl"

    classified_ids = _get_classified_ids(conn)
    if not classified_ids:
        log.info("  No classified sequences for pass-2")
        return []

    log.info(f"  pass-2 ({aligner}): {len(classified_ids)} classified sequences as targets")

    if outdir is None:
        outdir = tempfile.mkdtemp(prefix="tesorter2_minimap2_")
    work = os.path.join(outdir, "minimap2_pass2")

    db_fasta, qry_fasta, db_seq_to_dbs = split_classified_unclassified(
        input_fasta, classified_ids, work, seq_type=seq_type)

    if os.path.getsize(qry_fasta) == 0:
        log.info("  No unclassified sequences to search")
        return []

    if hmm_classifications is None:
        hmm_classifications = {}

    if pass2_classified_fasta:
        updated = pass2_external.update_classified_fasta_headers(
            pass2_classified_fasta, hmm_classifications, work
        )
        src = updated or pass2_classified_fasta
        pass2_external.extend_hmm_classifications_from_fasta(
            hmm_classifications, src, db_seq_to_dbs
        )
        merged_db = os.path.join(work, "pass2_db_merged.fa")
        pass2_external.merge_classified_fastas(
            merged_db, db_fasta, src, clean_nucl=True
        )
        db_fasta = merged_db

    if os.path.getsize(db_fasta) == 0:
        log.info("  pass-2 target FASTA is empty; skipping pass-2")
        return []

    if aligner == "blast":
        from . import blast_backend
        new_cls = blast_backend.run_pass2_blast(
            qry_fasta=qry_fasta, db_fasta=db_fasta, conn=conn,
            classifications=hmm_classifications,
            db_seq_to_dbs=db_seq_to_dbs,
            n_processors=n_processors,
            min_identity=min_identity,
            min_coverage=min_coverage,
            min_length=min_length,
            work=work,
            blast_task=blast_task,
        )
        log.info(f"  blast pass-2 total: {time.time() - t0:.1f}s")
        return new_cls

    paf_out = os.path.join(work, "pass2.paf")
    tsv_out = os.path.join(work, "pass2.tsv")

    log.info(f"  Running minimap2 -x {preset} with {n_processors} threads")
    t1 = time.time()
    run_alignment(
        query_fa=qry_fasta, target_fa=db_fasta, paf_out=paf_out,
        ncpu=n_processors, preset=preset, extra=minimap2_extra,
    )
    t2 = time.time()
    log.info(f"  minimap2 alignment: {t2 - t1:.1f}s")

    min_pid = min_identity / 100.0
    min_cov = min_coverage / 100.0
    tsv_rows = classify_paf_to_tsv(
        paf_out, tsv_out,
        min_pid=min_pid, min_qcov=min_cov, min_tcov=min_cov,
    )

    if tsv_rows:
        store_blast_hits(conn, tsv_rows, db_seq_to_dbs)

    log.info(f"  {len(hmm_classifications)} HMM classifications available for inheritance")

    new_cls = classify_from_blast(tsv_rows, hmm_classifications)

    t3 = time.time()
    log.info(f"  minimap2 pass-2 total: {t3 - t0:.1f}s")
    return new_cls
