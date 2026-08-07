#!/usr/bin/env python3
"""Classify putative LTR-RTs from a minimap2 PAF -- CIGAR-free version.

This is a faster, memory-leaner sibling of classify_ltr_paf.py that uses ONLY
the standard PAF columns plus the dv:f tag (approximate per-base sequence
divergence). It does NOT walk the cg:Z CIGAR. The intent is to let you run
minimap2 WITHOUT -c, saving substantial runtime on large inputs:

    minimap2 -k15 -w5 -A2 -B3 -r1k,10k -s30 -m30 -N50 -p0.1 \\
        target.fa query.fa > out.paf       # NOTE: no -c, no cg:Z, dv:f tag

Approximation tradeoff: q- and t-intervals for overlap dedup are taken as
[qstart, qend] and [tstart, tend], i.e. the alignment SPAN, which lumps
gap-induced bases (insertions and deletions) in with the matches+mismatches.
For LTR-RTs at >=70% identity, gap content is typically 2-5%, so eff_qcov /
eff_tcov are inflated by that small amount. dv:f is used in place of de:f
(de:f is gap-compressed; dv:f is approximate per-base divergence -- close
enough at our 70% threshold).

If both dv:f and de:f are present (e.g. the user kept -c for some reason),
de:f is used. Otherwise dv:f. If neither is present, the line is skipped.

OUTPUT (TSV; --header to add a header):
    qname  pass  pid  eff_qcov  eff_tcov  best_tname

PARSIMONY RULE for pid (same as the CIGAR-aware sibling):
    1. Drop alignments fully encompassed by another on q (encompassing aln has
       more statistical power).
    2. Sort survivors by dv:f DESCENDING (densest first).
    3. Greedy allocation, densest-first: each alignment's mutations are first
       absorbed into already-claimed (denser) overlap; leftovers spill into
       its unique region.
    4. pid = 1 - sum(unique_mutations) / sum(unique_q_spans).

Per-alignment "mutations" weight: dv * qspan (where qspan = qend - qstart).

Best-target per query: passes-rule first, then joint_score = pid*qcov*tcov,
then tname alphabetical.
"""

import argparse
import sys
from collections import defaultdict


def merge_intervals(ivs):
    """Merge overlapping/adjacent intervals; return sorted, disjoint list."""
    if not ivs:
        return []
    ivs = sorted(ivs)
    out = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= out[-1][1]:
            if e > out[-1][1]:
                out[-1][1] = e
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def interval_total(ivs):
    return sum(e - s for s, e in ivs)


def interval_difference(A, B):
    """Return A \\ B: parts of A not covered by B."""
    A = merge_intervals(A)
    B = merge_intervals(B)
    out = []
    bi = 0
    for s, e in A:
        cur = s
        while bi < len(B) and B[bi][1] <= cur:
            bi += 1
        j = bi
        while cur < e and j < len(B):
            bs, be = B[j]
            if bs >= e:
                break
            if bs > cur:
                out.append((cur, min(bs, e)))
            cur = max(cur, be)
            j += 1
        if cur < e:
            out.append((cur, e))
    return out


def interval_intersection_length(A, B):
    A = merge_intervals(A)
    B = merge_intervals(B)
    i = j = 0
    total = 0
    while i < len(A) and j < len(B):
        s = max(A[i][0], B[j][0])
        e = min(A[i][1], B[j][1])
        if s < e:
            total += e - s
        if A[i][1] < B[j][1]:
            i += 1
        else:
            j += 1
    return total


def is_fully_encompassed(small_iv, big_iv):
    """True iff every interval in small_iv is fully inside some interval in big_iv."""
    big_iv = merge_intervals(big_iv)
    for s, e in small_iv:
        contained = False
        for bs, be in big_iv:
            if bs <= s and e <= be:
                contained = True
                break
        if not contained:
            return False
    return True


def parse_paf_line(line, lineno, strict=False):
    """Return alignment record dict, or None if malformed and strict=False.

    Required: 12 standard columns + a dv:f or de:f tag.
    """
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 12:
        if strict:
            raise ValueError("line %d: only %d fields (need >=12)" % (lineno, len(fields)))
        return None
    dv = None
    de = None
    for tag in fields[12:]:
        if tag.startswith("dv:f:"):
            try:
                dv = float(tag[5:])
            except ValueError:
                pass
        elif tag.startswith("de:f:"):
            try:
                de = float(tag[5:])
            except ValueError:
                pass
    # Prefer the gap-compressed de:f if present (more accurate); fall back to dv:f.
    div = de if de is not None else dv
    if div is None:
        if strict:
            raise ValueError("line %d: neither dv:f nor de:f tag present" % lineno)
        return None
    try:
        return {
            "qname": fields[0],
            "qlen": int(fields[1]),
            "qstart": int(fields[2]),
            "qend": int(fields[3]),
            "strand": fields[4],
            "tname": fields[5],
            "tlen": int(fields[6]),
            "tstart": int(fields[7]),
            "tend": int(fields[8]),
            "div": div,
        }
    except (ValueError, IndexError):
        if strict:
            raise
        return None


def chained_intervals(alignments, max_gap=5000, gap_tol=0.20):
    """Group alignments into colinear chains; return per-chain (q,t) envelopes.

    Two adjacent alignments (sorted by qstart) join the same chain iff:
      1. Same strand.
      2. t-order consistent with strand (ascending t for +, descending t for -).
      3. Inner gap on each axis <= max_gap (bp).
      4. |q_gap - t_gap| / max(q_gap, t_gap) <= gap_tol  (synchronized indel).
    Otherwise a new chain begins.

    Returns ([(qmin, qmax), ...], [(tmin, tmax), ...]) -- one (q,t) envelope per
    chain. Multiple chains may overlap on q or t; the caller is responsible for
    merging across chains.
    """
    if not alignments:
        return [], []
    alns = sorted(alignments, key=lambda a: (a["qstart"], a["qend"]))
    chains = [[alns[0]]]
    for nxt in alns[1:]:
        cur = chains[-1][-1]
        join = cur["strand"] == nxt["strand"]
        if join:
            if cur["strand"] == "+":
                t_order_ok = nxt["tstart"] >= cur["tstart"]
                t_gap = max(0, nxt["tstart"] - cur["tend"])
            else:
                t_order_ok = nxt["tstart"] <= cur["tstart"]
                t_gap = max(0, cur["tstart"] - nxt["tend"])
            join = t_order_ok
        if join:
            q_gap = max(0, nxt["qstart"] - cur["qend"])
            big = max(q_gap, t_gap)
            if big > max_gap:
                join = False
            elif big > 0 and abs(q_gap - t_gap) / big > gap_tol:
                join = False
        if join:
            chains[-1].append(nxt)
        else:
            chains.append([nxt])
    q_ivs, t_ivs = [], []
    for chain in chains:
        q_ivs.append((min(a["qstart"] for a in chain),
                      max(a["qend"] for a in chain)))
        t_ivs.append((min(a["tstart"] for a in chain),
                      max(a["tend"] for a in chain)))
    return q_ivs, t_ivs


def compute_pair_metrics(alignments, qlen, tlen,
                         fill_colinear=False, max_gap=5000, gap_tol=0.20):
    """Compute (pid, eff_qcov, eff_tcov) for one (q,t) pair, span-based.

    `alignments` is a non-empty list of dicts: {qstart, qend, tstart, tend, div}.
    If fill_colinear is True, colinear-chain inner gaps on q and t are bridged
    in the qcov/tcov computation (pid logic is unchanged either way).
    """
    # Build per-alignment span intervals (no CIGAR -- one interval each).
    for a in alignments:
        a["q_iv"] = [(a["qstart"], a["qend"])]
        a["t_iv"] = [(a["tstart"], a["tend"])]
        a["qspan"] = a["qend"] - a["qstart"]

    # 1. Drop alignments whose q-span is fully encompassed by a STRICTLY longer one.
    survivors = []
    for i, A in enumerate(alignments):
        encompassed = False
        for j, B in enumerate(alignments):
            if i == j:
                continue
            if A["qspan"] < B["qspan"] and is_fully_encompassed(A["q_iv"], B["q_iv"]):
                encompassed = True
                break
        if not encompassed:
            survivors.append(A)
    if not survivors:
        survivors = [max(alignments, key=lambda a: a["div"])]

    # 2. Sort by div descending.
    survivors.sort(key=lambda a: -a["div"])

    # 3. Greedy allocation.
    claims = []
    claimed_q_all = []
    for X in survivors:
        X_q = X["q_iv"]
        X_total_mut = X["div"] * X["qspan"]
        mutations_in_overlap = 0.0
        for unique_iv, density in claims:
            mutations_in_overlap += density * interval_intersection_length(X_q, unique_iv)
        X_unique_iv = interval_difference(X_q, claimed_q_all)
        X_unique_len = interval_total(X_unique_iv)
        X_unique_mut = max(0.0, X_total_mut - mutations_in_overlap)
        if X_unique_len > 0:
            X_unique_mut = min(X_unique_mut, float(X_unique_len))
            density = X_unique_mut / X_unique_len
        else:
            X_unique_mut = 0.0
            density = 0.0
        claims.append((X_unique_iv, density))
        claimed_q_all = merge_intervals(claimed_q_all + X_unique_iv)

    total_unique_mut = sum(d * interval_total(iv) for iv, d in claims)
    total_unique_len = interval_total(claimed_q_all)
    pid = (1.0 - total_unique_mut / total_unique_len) if total_unique_len > 0 else 0.0

    # eff_qcov, eff_tcov: union of all alignment spans (use ALL alignments,
    # not just survivors -- encompassed ones add nothing new anyway). When
    # fill_colinear is on, inner gaps within colinear chains are bridged before
    # the union so a single element fragmented into HSPs counts as one span.
    if fill_colinear:
        q_chains, t_chains = chained_intervals(
            alignments, max_gap=max_gap, gap_tol=gap_tol)
        eff_qcov = interval_total(merge_intervals(q_chains)) / qlen if qlen > 0 else 0.0
        eff_tcov = interval_total(merge_intervals(t_chains)) / tlen if tlen > 0 else 0.0
    else:
        all_q_iv = [iv for a in alignments for iv in a["q_iv"]]
        all_t_iv = [iv for a in alignments for iv in a["t_iv"]]
        eff_qcov = interval_total(merge_intervals(all_q_iv)) / qlen if qlen > 0 else 0.0
        eff_tcov = interval_total(merge_intervals(all_t_iv)) / tlen if tlen > 0 else 0.0
    return pid, eff_qcov, eff_tcov


def process_paf(lines, min_pid=0.70, min_qcov=0.70, min_tcov=0.70,
                fill_colinear=False, max_gap=5000, gap_tol=0.20, verbose=False):
    """Stream PAF lines, group by (qname,tname), pick best target per query,
    return [(qname, pass_str, pid, qcov, tcov, best_tname), ...] sorted by qname.
    """
    pair_alns = defaultdict(list)
    pair_lengths = {}
    n_lines = 0
    n_skipped = 0
    for lineno, raw in enumerate(lines, start=1):
        if not raw.strip() or raw.startswith("#") or raw.startswith("["):
            continue
        rec = parse_paf_line(raw, lineno, strict=False)
        if rec is None:
            n_skipped += 1
            continue
        key = (rec["qname"], rec["tname"])
        pair_alns[key].append(rec)
        if key not in pair_lengths:
            pair_lengths[key] = (rec["qlen"], rec["tlen"])
        n_lines += 1
    if verbose:
        print("[classify_ltr_paf_fast] parsed %d alignments across %d (q,t) pairs"
              % (n_lines, len(pair_alns)), file=sys.stderr)
        if n_skipped:
            print("[classify_ltr_paf_fast] skipped %d malformed lines (missing dv:f/de:f or unparseable)"
                  % n_skipped, file=sys.stderr)

    per_query = defaultdict(list)
    for (qname, tname), alns in pair_alns.items():
        qlen, tlen = pair_lengths[(qname, tname)]
        pid, qcov, tcov = compute_pair_metrics(
            alns, qlen, tlen,
            fill_colinear=fill_colinear, max_gap=max_gap, gap_tol=gap_tol)
        passes = (pid >= min_pid) and (qcov >= min_qcov) and (tcov >= min_tcov)
        per_query[qname].append({
            "tname": tname,
            "pid": pid,
            "qcov": qcov,
            "tcov": tcov,
            "passes": passes,
            "joint": pid * qcov * tcov,
        })

    results = []
    n_pass = 0
    for qname in sorted(per_query):
        candidates = per_query[qname]
        candidates.sort(key=lambda c: (-int(c["passes"]), -c["joint"], c["tname"]))
        best = candidates[0]
        if best["passes"]:
            n_pass += 1
        results.append((
            qname,
            "pass" if best["passes"] else "fail",
            best["pid"],
            best["qcov"],
            best["tcov"],
            best["tname"],
        ))
    if verbose:
        print("[classify_ltr_paf_fast] %d/%d queries pass at pid>=%.3f qcov>=%.3f tcov>=%.3f"
              % (n_pass, len(results), min_pid, min_qcov, min_tcov), file=sys.stderr)
    return results


def format_row(row):
    qname, pass_str, pid, qcov, tcov, tname = row
    return "%s\t%s\t%.4f\t%.4f\t%.4f\t%s" % (qname, pass_str, pid, qcov, tcov, tname)


def main():
    ap = argparse.ArgumentParser(
        description="Classify putative LTR-RTs from a minimap2 PAF (CIGAR-free, "
                    "uses dv:f or de:f and standard PAF columns only).")
    ap.add_argument("paf", help="input PAF (use - for stdin)")
    ap.add_argument("-o", "--output", default="-",
                    help="output TSV path (default: stdout)")
    ap.add_argument("--min-pid", type=float, default=0.70)
    ap.add_argument("--min-qcov", type=float, default=0.70)
    ap.add_argument("--min-tcov", type=float, default=0.70)
    ap.add_argument("--fill-colinear-gaps", action="store_true",
                    help="bridge inner gaps within colinear HSP chains "
                         "(same strand, q-order matches t-order, q-gap ~ t-gap) "
                         "before computing eff_qcov/eff_tcov. Off by default.")
    ap.add_argument("--bridge-max-gap", type=int, default=5000,
                    help="max inner gap (bp) to bridge on either axis (default: 5000).")
    ap.add_argument("--bridge-gap-tol", type=float, default=0.20,
                    help="max relative mismatch |q_gap - t_gap| / max(q_gap, t_gap) "
                         "to treat as a synchronized indel (default: 0.20).")
    ap.add_argument("--header", action="store_true")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    in_fh = sys.stdin if args.paf == "-" else open(args.paf, "r")
    out_fh = sys.stdout if args.output == "-" else open(args.output, "w")
    if args.verbose:
        print("[classify_ltr_paf_fast] reading %s" % args.paf, file=sys.stderr)
    try:
        results = process_paf(
            in_fh,
            min_pid=args.min_pid, min_qcov=args.min_qcov, min_tcov=args.min_tcov,
            fill_colinear=args.fill_colinear_gaps,
            max_gap=args.bridge_max_gap, gap_tol=args.bridge_gap_tol,
            verbose=args.verbose,
        )
    finally:
        if args.paf != "-":
            in_fh.close()
    try:
        if args.header:
            out_fh.write("qname\tpass\tpid\teff_qcov\teff_tcov\tbest_tname\n")
        for row in results:
            out_fh.write(format_row(row) + "\n")
    finally:
        if args.output != "-":
            out_fh.close()
    if args.verbose:
        print("[classify_ltr_paf_fast] done", file=sys.stderr)


if __name__ == "__main__":
    main()
