"""queue_search.py — (chunk x database) tasks drained from one queue per level.

`partition.py` gives each worker a chunk and has it run the whole cascade for
every database. That balances badly: groups are packed to 0.00% on residues but
cost tracks what survives the filters, and a 64-group run packed at only 70% --
327 s of workers idle behind slower groups. The unit of work is too big, and a
worker that draws an expensive chunk carries it alone.

Here the unit is one (chunk, database) search at one cascade level. Every such
task is independent within its level, so 256 chunks x 7 databases is ~1,800
tasks for 64 workers to drain, and a slow chunk is shared out across workers
rather than owned by one. The tail becomes one task instead of one chunk's
entire cascade.

The level barrier is real and is the whole point
------------------------------------------------
A level completes for a chunk only when every database has searched it. Then
the level's results are reconciled across databases, and whatever the vote
settles is finished -- no later engine of any database sees it again. Level N+1
is queued only after that. Retiring a sequence at level 0 is exactly what makes
nail worth running: it removes work from every database's later stages, not
just its own.

The manager owns the database
-----------------------------
Workers return hits and classifications rather than writing them. There is no
per-worker results database and no merge: the parent inserts as tasks land,
spread across the run instead of concentrated at the end. A task carries on the
order of a thousand hit rows, which is cheap to hand back.
"""

import collections
import logging
import multiprocessing
import os

from .hmm import AMINO_ALPHABET

log = logging.getLogger(__name__)

_CTX = {}


def _init(sources_by_chunk, db_paths, search_space_mb, evalue_scale,
          compat_rounding, compat_voting, min_clade_delta):
    """Per-worker constants. Inherited through fork, so this costs nothing."""
    _CTX.update(sources=sources_by_chunk, db_paths=db_paths,
                search_space_mb=search_space_mb, evalue_scale=evalue_scale,
                compat_rounding=compat_rounding, compat_voting=compat_voting,
                min_clade_delta=min_clade_delta)


def _translate_task(job):
    """Six-frame translate one chunk. Its own task so 64 chunks translate at
    once rather than the manager doing them in series."""
    from .sequence import translate_fasta
    from . import nail_search
    chunk, nucl, aa_path, nostop_path, mask_stops = job
    translate_fasta(nucl, aa_path, mask_stops=mask_stops)
    if nostop_path:
        if mask_stops:
            nostop_path = aa_path
        else:
            nostop_path, _ = nail_search.prepare_targets(aa_path, nostop_path)
    return chunk, aa_path, nostop_path


def _search_task(job):
    """One (chunk, database) search at one level."""
    from .hierarchical_search import (_materialize, hits_to_arrays,
                                      _EVALUE_FIELDS, ENGINES)
    from .classifier import classify_sequences, DB_CONFIGS
    from .subset import base_name
    import time

    chunk, db_name, level, engine_name, targets, workdir = job
    engine = ENGINES[engine_name]()
    os.makedirs(workdir, exist_ok=True)
    sources = _CTX["sources"][chunk]

    fasta = _materialize(engine.input_kind, set(targets), sources["all_names"],
                         sources, workdir, "%s_L%d" % (db_name, level))
    if fasta is None:
        return chunk, db_name, level, [], [], 0.0

    t0 = time.time()
    hits = engine.search(_CTX["db_paths"][db_name], fasta, db_name, workdir, 1,
                         search_space_mb=_CTX["search_space_mb"])
    wall = time.time() - t0

    scale = _CTX["evalue_scale"].get(chunk, 1.0)
    if scale != 1.0 and engine.name == "nail":
        for h in hits:
            for f in _EVALUE_FIELDS:
                if f in h:
                    h[f] *= scale

    hits = [h for h in hits if base_name(h["target_name"]) in targets]
    arrays = hits_to_arrays(hits)
    results = []
    if arrays is not None:
        results = classify_sequences(
            arrays, DB_CONFIGS[db_name],
            compat_rounding=_CTX["compat_rounding"],
            compat_voting=_CTX["compat_voting"],
            min_clade_delta=_CTX["min_clade_delta"])
        for r in results:
            r["engine"] = engine.name
            r["stage"] = level
    return chunk, db_name, level, hits, results, wall


def run_queued(conn, input_fasta, db_paths, db_alphabets, outdir, n_groups,
               n_workers, protein_stages, seq_type="nucl", mask_stops=False,
               compat_rounding=False, compat_voting=False, min_clade_delta=0.0,
               seq_index=None, lengths=None, dna_level=0):
    """Search in (chunk x database) tasks, one queue per cascade level."""
    from .hierarchical_search import (build_stages, stages_for_input,
                                      DNA_STAGES, DNA_ALPHABET, ENGINES,
                                      _record_exits)
    from .classifier import (reconcile_classifications, store_classifications,
                             DB_CONFIGS)
    from .results import store_legacy
    from . import partition, bath_search, nail_search

    groups, sizes, total_residues = partition.plan_groups(
        input_fasta, n_groups, lengths=lengths)
    part_dir = os.path.join(outdir, "partitions")
    paths = partition.write_groups(input_fasta, groups, part_dir)
    search_space_mb = total_residues / 1e6
    protein_stages = stages_for_input(protein_stages, seq_type)

    # Which databases run at which level.
    plans, all_stages = {}, {}
    for db_name in db_paths:
        is_dna = db_alphabets[db_name] == DNA_ALPHABET
        if is_dna and seq_type == "prot":
            log.warning("  Skipping %s: DNA profiles need nucleotide input",
                        db_name)
            continue
        if DB_CONFIGS.get(db_name) is None:
            log.warning("  No classifier config for %s, skipping", db_name)
            continue
        names = DNA_STAGES if is_dna else protein_stages
        stages = build_stages(names)
        if any(s.name == "bath" for s in stages):
            bath_search.require_binaries()
        if any(s.name == "nail" for s in stages):
            nail_search.require_binaries()
        all_stages[db_name] = stages
        plans[db_name] = ({dna_level: stages[0].name} if is_dna
                          else {i: s.name for i, s in enumerate(stages)})
    if not plans:
        return {}
    n_levels = max(max(p) for p in plans.values()) + 1

    kinds = {ENGINES[n].input_kind for n in protein_stages if n in ENGINES}
    needs_aa = bool(kinds & {"aa", "aa_nostop"}) and any(
        db_alphabets[d] == AMINO_ALPHABET for d in plans)

    # Per-chunk sources, and the nail correction each chunk needs.
    sources_by_chunk, evalue_scale = {}, {}
    for i, (path, size) in enumerate(zip(paths, sizes)):
        sources_by_chunk[i] = {
            "nucl": None if seq_type == "prot" else path,
            "aa": path if seq_type == "prot" else None,
            "aa_nostop": None,
            "all_names": list(groups[i]),
        }
        evalue_scale[i] = (total_residues / size) if size else 1.0

    pool = multiprocessing.Pool(
        min(n_groups, n_workers), initializer=_init,
        initargs=(sources_by_chunk, db_paths, search_space_mb, evalue_scale,
                  compat_rounding, compat_voting, min_clade_delta))
    try:
        if needs_aa:
            log.info("Translating %d chunks in parallel", n_groups)
            jobs = []
            for i, path in enumerate(paths):
                d = os.path.join(part_dir, "chunk%03d" % i)
                os.makedirs(d, exist_ok=True)
                jobs.append((i, path, os.path.join(d, "aa.fa"),
                             os.path.join(d, "aa.nostop.fa")
                             if "aa_nostop" in kinds else None, mask_stops))
            for chunk, aa, nostop in pool.imap_unordered(_translate_task, jobs):
                sources_by_chunk[chunk]["aa"] = aa
                sources_by_chunk[chunk]["aa_nostop"] = nostop
            # Workers forked before this, so re-seed them with the paths.
            pool.close(); pool.join()
            pool = multiprocessing.Pool(
                min(n_groups, n_workers), initializer=_init,
                initargs=(sources_by_chunk, db_paths, search_space_mb,
                          evalue_scale, compat_rounding, compat_voting,
                          min_clade_delta))

        remaining = {i: set(groups[i]) for i in range(n_groups)}
        per_db = collections.defaultdict(list)

        for level in range(n_levels):
            runners = [(db, p[level]) for db, p in plans.items() if level in p]
            live = [i for i in range(n_groups) if remaining[i]]
            if not runners or not live:
                continue
            # One workdir per TASK, not per (chunk, level). Databases at a
            # level run concurrently on different workers here, and nail takes
            # a --tmp-dir under the workdir: sharing it means seven nail
            # processes writing the same scratch directory at once. The chunked
            # scheduler runs a chunk's databases in series, so it never hit
            # this.
            jobs = [(i, db, level, eng, frozenset(remaining[i]),
                     os.path.join(part_dir, "chunk%03d" % i, "L%d" % level, db))
                    for i in live for db, eng in runners]
            log.info("--- Level %d: %d tasks (%d chunks x %d databases) ---",
                     level, len(jobs), len(live), len(runners))

            by_chunk = collections.defaultdict(dict)
            n_hits = 0
            for chunk, db_name, lv, hits, results, wall in \
                    pool.imap_unordered(_search_task, jobs):
                if hits:
                    store_legacy(conn, hits, db_name, engine=plans[db_name][lv],
                                 seq_index=seq_index, stage=lv)
                    n_hits += len(hits)
                if results:
                    by_chunk[chunk][db_name] = results
                    per_db[db_name].extend(results)

            settled_total = 0
            for chunk in live:
                if not by_chunk[chunk]:
                    continue
                settled = {r["id"] for r in
                           reconcile_classifications(by_chunk[chunk])}
                remaining[chunk] -= settled
                settled_total += len(settled)
            left = sum(len(remaining[i]) for i in range(n_groups))
            log.info("  Level %d: %d hits stored, %d sequences resolved, "
                     "%d remaining", level, n_hits, settled_total, left)
    finally:
        pool.close()
        pool.join()

    for db_name, results in per_db.items():
        if results:
            store_classifications(conn, results, database=db_name,
                                  mode="hierarchical")
    return dict(per_db)
