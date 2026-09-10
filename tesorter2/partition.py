"""partition.py — split the input across independent single-threaded workers.

A `sequences` run spends most of its wall clock in two engines that parallelise
badly. Measured on 135,042 RepBase sequences over 64 cores: nhmmer 40.5% of the
run and nail 28.7%, against hmmer's 6.5% — and the whole run used about 26% of
its allocation. nhmmer is the clearest case, because `legacy_search_nucl`
deliberately leaves parallelism to pyhmmer, which spreads *queries* across
threads. With 87 models and 64 threads that is a ceiling of 87 with a long
tail, and it measured 13.4x.

This takes the other axis. The input is split into balanced groups and each
group is searched by its own worker, single-threaded, running the whole cascade
for every database. It composes with the existing parallelism rather than
competing with it: each worker asks for one thread, so nothing is
oversubscribed.

Why the split is clean
----------------------
The cascade's stage-exit decision is **per sequence** — a sequence is labelled
or handed to the next engine on its own evidence — so a group can run the whole
cascade independently and the results concatenate. Nothing crosses a group
boundary until cross-database reconciliation, which the caller does afterwards
on the merged set exactly as it did before.

Two things would otherwise make a verdict depend on which group a sequence
landed in, and both are corrected:

  * **E-values scale with search-space size.** Every engine but one honours a
    pinned search space (`bathsearch -Z`, nhmmer's megabase scaling; hmmsearch's
    Z is a model count and is unaffected), so workers are handed the *whole*
    input's size rather than their own.
  * **nail exposes no -Z.** Its E-values are computed against whatever target
    set it was given, so they are rescaled by `full_residues / group_residues`
    after the fact. That is exact rather than approximate, because the
    dependence is linear.

BATH is the known exception: it truncates an alignment depending on how much
sequence precedes the target, so a partitioned run is not byte-identical to an
unpartitioned one on the protein side. That is upstream behaviour, not
something this can correct.

Balance
-------
Greedy longest-first bin packing on sequence length. Measured on RepBase:
0.00% imbalance at every group count up to 256, because the longest sequence is
50,000 bp — 0.8% of a 64-way group. Balance is the binding constraint at small
inputs, not overhead: 200 sequences into 64 groups is 28% imbalanced and into
128 is 156%, while 1,000 into 128 is 9.9% and 5,000 into 128 is 0.1%.

Cost of a worker
----------------
Each worker loads every database itself. Measured with 128 workers loading at
once, that is 2.31 s for the default eight and at most 2.2x what a single
loader costs — there is no page-cache serialisation, so the models do not need
sharing through fork. Building the sequence block, by contrast, divides across
workers.
"""

import heapq
import logging
import multiprocessing
import os
import sqlite3

import pyfastx

log = logging.getLogger(__name__)

def sequence_lengths(fasta):
    """[(name, length)] straight out of pyfastx's index.

    The `.fxi` sidecar is a SQLite database -- one `seq` table with `chrom` and
    `slen` -- so every length comes back in a single query. Asking pyfastx for
    them one at a time instead (`len(fa[name])` over `fa.keys()`) is 135,042
    random-access lookups and measured **1,657 seconds** on RepBase against
    0.36 for the query, which was 51% of a partitioned run's wall clock and
    entirely serial.
    """
    pyfastx.Fasta(fasta, build_index=True)     # ensure the .fxi exists
    conn = sqlite3.connect(fasta + ".fxi")
    try:
        return list(conn.execute("SELECT chrom, slen FROM seq"))
    finally:
        conn.close()


def plan_groups(fasta, n_groups, lengths=None):
    """[[name, ...], ...] balanced on total residues, longest sequence first.

    Longest-first is what makes the packing tight: placing the big sequences
    while every bin is still nearly empty leaves the small ones to level the
    result, and no sequence can be split.

    `lengths` lets a caller hand over sizes it already holds -- the pipeline
    measures every sequence before searching -- so neither the index query nor
    an index build is needed for a temporary subset file.
    """
    if lengths is None:
        lengths = sequence_lengths(fasta)
    elif isinstance(lengths, dict):
        lengths = list(lengths.items())
    else:
        lengths = list(lengths)
    lengths.sort(key=lambda x: -x[1])

    bins = [(0, i) for i in range(n_groups)]
    heapq.heapify(bins)
    groups = [[] for _ in range(n_groups)]
    for name, size in lengths:
        total, idx = heapq.heappop(bins)
        groups[idx].append(name)
        heapq.heappush(bins, (total + size, idx))

    size_of = dict(lengths)
    sizes = [sum(size_of[n] for n in g) for g in groups]
    total = sum(sizes) or 1
    ideal = total / n_groups
    log.info("Partitioned %d sequences into %d groups: %.2f Mb each, "
             "imbalance %.2f%%", len(lengths), n_groups, ideal / 1e6,
             100 * (max(sizes) - ideal) / ideal if ideal else 0.0)
    return groups, sizes, total


def write_groups(fasta, groups, outdir):
    """One FASTA per group. Returns their paths, in group order."""
    os.makedirs(outdir, exist_ok=True)
    where = {}
    for i, g in enumerate(groups):
        for name in g:
            where[name] = i
    # Read once into memory, then write each group in a single pass. The
    # alternative -- 64 open handles taking interleaved small writes -- turns
    # one sequential read into scattered I/O across a parallel filesystem.
    buckets = [[] for _ in groups]
    for name, seq in pyfastx.Fasta(fasta, build_index=False):
        i = where.get(name)
        if i is not None:
            buckets[i].append(">%s\n%s\n" % (name, seq))
    paths = []
    for i, chunk in enumerate(buckets):
        path = os.path.join(outdir, "part%03d.fa" % i)
        with open(path, "w") as fh:
            fh.write("".join(chunk))
        buckets[i] = None
        paths.append(path)
    return paths


def _worker(job):
    """One group: its own results database, its own full cascade.

    Runs in a forked process, so it must not touch the parent's sqlite
    connection. It returns the per-database classification dicts and the path
    to its own database, which the parent merges.
    """
    (idx, fasta, db_paths, db_alphabets, outdir, kwargs,
     search_space_mb, evalue_scale, seq_index) = job
    from .results import create_db
    from . import hierarchical_search

    work = os.path.join(outdir, "part%03d" % idx)
    os.makedirs(work, exist_ok=True)
    db_path = os.path.join(work, "part.db")
    conn = create_db(db_path)
    try:
        from .results import store_sequences
        lengths = {name: len(seq) for name, seq in
                   pyfastx.Fasta(fasta, build_index=False)}
        store_sequences(conn, lengths)
        per_db = hierarchical_search.run_cascade(
            conn, fasta, db_paths, db_alphabets, work,
            search_space_mb=search_space_mb, evalue_scale=evalue_scale,
            seq_index=seq_index, **kwargs)
        conn.commit()
    finally:
        conn.close()
    return idx, per_db, db_path


def merge_one(conn, path):
    """Copy one partition's rows into the parent database.

    Generic over the schema: tables come from each partition's own
    sqlite_master, so one that only some partitions created (classifications is
    made on first write) still arrives.

    **An INTEGER PRIMARY KEY is dropped from the copy.** Each partition is its
    own database and numbers `legacy_hits.id` from 1, so copying the column
    verbatim makes every group after the first collide with the first -- and
    with `INSERT OR IGNORE` those rows vanish without an error. Inserting the
    other columns explicitly lets the parent assign fresh ids. `OR IGNORE` is
    kept only for tables keyed on something real, like `sequences.name`, which
    the parent has already written.

    The row counts are checked rather than trusted: a table that loses rows in
    the copy raises instead of leaving a quietly short database.
    """
    if not os.path.exists(path):
        return 0
    n_rows = 0
    conn.execute("ATTACH DATABASE ? AS part", (path,))
    try:
        tables = [r[0] for r in conn.execute(
            "SELECT name FROM part.sqlite_master WHERE type='table' "
            "AND name NOT LIKE 'sqlite_%'")]
        for t in tables:
            ddl = conn.execute(
                "SELECT sql FROM part.sqlite_master WHERE type='table' "
                "AND name=?", (t,)).fetchone()
            if ddl and ddl[0]:
                conn.execute(ddl[0].replace(
                    "CREATE TABLE", "CREATE TABLE IF NOT EXISTS", 1))
            info = list(conn.execute("PRAGMA part.table_info(%s)" % t))
            # (cid, name, type, notnull, default, pk)
            keyed = any(c[5] for c in info)
            want = conn.execute(
                "SELECT COUNT(*) FROM part.%s" % t).fetchone()[0]
            before = conn.execute(
                "SELECT COUNT(*) FROM main.%s" % t).fetchone()[0]
            # `SELECT *` returns the declared columns only, never the
            # rowid, so the parent assigns its own and every row appends --
            # whatever order the groups happen to finish in. OR IGNORE only
            # where the key is real data the parent may already hold, as
            # with sequences.name.
            verb = "INSERT OR IGNORE" if keyed else "INSERT"
            conn.execute("%s INTO main.%s SELECT * FROM part.%s"
                         % (verb, t, t))
            after = conn.execute(
                "SELECT COUNT(*) FROM main.%s" % t).fetchone()[0]
            if not keyed and after - before != want:
                raise RuntimeError("merging %s from %s copied %d of %d rows"
                                   % (t, path, after - before, want))
            n_rows += after - before
        conn.commit()
    except Exception:
        # Roll back first: a failed insert leaves the transaction open and
        # DETACH then fails with "database is locked", replacing the real
        # error with a misleading one.
        conn.rollback()
        raise
    finally:
        conn.execute("DETACH DATABASE part")
    return n_rows


def run_partitioned(conn, input_fasta, db_paths, db_alphabets, outdir,
                    n_groups, n_workers, cascade_kwargs, seq_index=None,
                    lengths=None):
    """Search the input in `n_groups` independent groups. Returns per_db.

    `seq_index` maps a sequence name to its position in the original input and
    is written to each hit's `id`. It is the same map for every worker: the
    position is a property of the input, not of the group a sequence landed in,
    so nothing has to be reserved, strided or renumbered.
    """
    groups, sizes, total_residues = plan_groups(input_fasta, n_groups,
                                                lengths=lengths)
    part_dir = os.path.join(outdir, "partitions")
    paths = write_groups(input_fasta, groups, part_dir)

    search_space_mb = total_residues / 1e6
    jobs = []
    for i, (p, size) in enumerate(zip(paths, sizes)):
        # nail's E-values improve in proportion to how much smaller its target
        # set is; this puts them back on the whole input's scale.
        scale = (total_residues / size) if size else 1.0
        jobs.append((i, p, db_paths, db_alphabets, part_dir, cascade_kwargs,
                     search_space_mb, scale, seq_index))

    log.info("Searching %d groups across %d workers, one thread each",
             n_groups, min(n_groups, n_workers))
    per_db = {}
    merged = 0
    # Each group is merged as it lands, not after the pool drains. The parent
    # is otherwise idle while the workers search, and doing all 64 copies at
    # the end concentrates every row of I/O into one burst against a filesystem
    # that already has 64 workers on it -- so all but the last merge happens in
    # the shadow of the groups still running.
    #
    # Yields stay unordered. Ordering them would let one slow group hold back
    # every merge behind it, and there is no reason to: `id` is the input
    # ordinal rather than the rowid, so a row's value no longer decides where
    # it lands and appends are appends whatever order groups finish in.
    with multiprocessing.Pool(min(n_groups, n_workers)) as pool:
        for idx, part, db_path in pool.imap_unordered(_worker, jobs):
            for name, results in part.items():
                per_db.setdefault(name, []).extend(results)
            merged += merge_one(conn, db_path)
            log.info("  group %d done (%s), merged",
                     idx, ", ".join("%s=%d" % (k, len(v))
                                    for k, v in part.items()) or "nothing")
    log.info("  Merged %d rows from %d partition databases", merged, n_groups)
    return per_db
