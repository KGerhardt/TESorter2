"""
decompose_hmm.py — Decompose HMM models into high-scoring sub-model domains.

Standalone module with no project-specific dependencies. Takes pyhmmer HMM
objects in, returns pyhmmer HMM objects out. Works with DNA and amino acid
alphabets.

Core API:
    compute_log_odds(hmm)           → (M+1, K) position score matrix
    score_positions(hmm)            → length-M array of per-position scores
    find_domains(scores, window)    → [(start, end, score), ...]
    decompose_model(hmm, window)    → (name, M, domains, scores)
    splice_sub_hmm(hmm, start, end) → searchable sub-HMM pyhmmer object

Higher-level (takes file paths, loads HMMs internally):
    decompose_file(path, window)    → {name: [(start, end, score), ...]}
    build_sub_hmms_from_file(path)  → [(sub_hmm, parent_name, start, end)]

Usage:
    python decompose_hmm.py database.hmm --window 64 --overlap 0.33
"""

import argparse
import io
import sys

import numpy as np
import pyhmmer.easel as easel
import pyhmmer.plan7 as plan7

# ---------------------------------------------------------------------------
# Background frequencies for log-odds computation
# ---------------------------------------------------------------------------

DNA_BG = np.array([0.2869, 0.1899, 0.2161, 0.3072])  # A C G T

# Robinson-Robinson amino acid background (HMMER default)
AA_BG = np.array([
    0.0787945, 0.0151600, 0.0535222, 0.0668298,  # A C D E
    0.0397062, 0.0695071, 0.0229198, 0.0590092,  # F G H I
    0.0594422, 0.0963728, 0.0237718, 0.0414386,  # K L M N
    0.0482904, 0.0395639, 0.0540978, 0.0683364,  # P Q R S
    0.0540687, 0.0673417, 0.0114135, 0.0304133,  # T V W Y
])


# ---------------------------------------------------------------------------
# Core functions — operate on pyhmmer HMM objects, no file I/O
# ---------------------------------------------------------------------------

def compute_log_odds(hmm):
    """Compute log-odds emission scores for all match positions.

    Auto-detects DNA vs amino acid from the HMM's alphabet.

    Args:
        hmm: pyhmmer.plan7.HMM object

    Returns:
        (M+1, K) numpy array. Row 0 is unused (HMMER convention).
        K = 4 for DNA, 20 for amino acid.
    """
    M = hmm.M
    me = hmm.match_emissions

    if hmm.alphabet == easel.Alphabet.dna():
        bg = DNA_BG
        K = 4
    else:
        bg = AA_BG
        K = 20

    log_odds = np.zeros((M + 1, K))
    for i in range(1, M + 1):
        for j in range(K):
            p = me[i][j]
            log_odds[i, j] = np.log2(p / bg[j]) if p > 0 else -100
    return log_odds


def score_positions(hmm):
    """Compute per-position informativeness score.

    Returns array of length M with the max log-odds at each position.
    Higher = more informative (strongest base/residue is more distinct
    from background).
    """
    lo = compute_log_odds(hmm)
    return np.max(lo[1:], axis=1)  # skip row 0


def find_domains(position_scores, window_size, max_domains=None,
                 min_score_per_base=0.0, max_overlap_frac=0.0,
                 min_coverage=0.0, max_coverage=1.0):
    """Greedy window selection by total score with configurable overlap.

    Finds the highest-scoring non-overlapping (or partially overlapping)
    windows in a position score array. These are the most informative
    regions of the model.

    Args:
        position_scores: array of per-position scores, length M
        window_size: sub-model size in positions
        max_domains: optional cap on number of domains to return
        min_score_per_base: minimum average score per position to accept
        max_overlap_frac: maximum fraction of window that can overlap
                          with already-selected domains (0.0 = none,
                          0.33 = up to 33%)
        min_coverage: minimum fraction of parent model that must be
                      covered by selected domains (0.0 = no minimum).
                      Selection continues until this is met even if
                      score falls below min_score_per_base.
        max_coverage: maximum fraction of parent model to cover
                      (1.0 = no limit). Selection stops when reached.

    Returns:
        list of (start, end, total_score) sorted by model position.
        Positions are 0-based.
    """
    M = len(position_scores)
    if M < window_size:
        return [(0, M, float(np.sum(position_scores)))]

    max_overlap = int(window_size * max_overlap_frac)

    # Sliding window sums
    n_windows = M - window_size + 1
    window_scores = np.zeros(n_windows)
    window_scores[0] = np.sum(position_scores[:window_size])
    for i in range(1, n_windows):
        window_scores[i] = (window_scores[i - 1]
                            - position_scores[i - 1]
                            + position_scores[i + window_size - 1])

    # Greedy selection
    used = np.zeros(M, dtype=int)
    domains = []

    while True:
        if max_domains is not None and len(domains) >= max_domains:
            break

        # Check max coverage
        covered = int(np.sum(used > 0))
        if covered / M >= max_coverage:
            break

        valid_scores = window_scores.copy()
        for i in range(n_windows):
            overlap_count = int(np.sum(used[i:i + window_size] > 0))
            if overlap_count > max_overlap:
                valid_scores[i] = -np.inf

        best_idx = np.argmax(valid_scores)
        best_score = valid_scores[best_idx]

        if best_score == -np.inf:
            break

        # Allow low-scoring windows if we haven't met min_coverage
        coverage_so_far = covered / M
        if (best_score / window_size < min_score_per_base
                and coverage_so_far >= min_coverage):
            break

        start = best_idx
        end = best_idx + window_size
        domains.append((start, end, float(best_score)))
        used[start:end] += 1

    domains.sort(key=lambda d: d[0])
    return domains


def decompose_model(hmm, window_size, max_domains=None,
                    min_score_per_base=0.0, max_overlap_frac=0.0,
                    min_coverage=0.0, max_coverage=1.0):
    """Decompose one HMM into high-scoring sub-model domains.

    Args:
        hmm: pyhmmer.plan7.HMM object
        window_size: sub-model window size in positions
        max_domains: optional cap
        min_score_per_base: minimum average score to accept a window
        max_overlap_frac: overlap budget between windows
        min_coverage: minimum parent model coverage fraction
        max_coverage: maximum parent model coverage fraction

    Returns:
        (name, M, domains, position_scores) where domains is a list
        of (start, end, score) tuples with 0-based positions.
    """
    name = hmm.name if isinstance(hmm.name, str) else hmm.name.decode()
    M = hmm.M

    ps = score_positions(hmm)
    domains = find_domains(ps, window_size, max_domains,
                           min_score_per_base, max_overlap_frac,
                           min_coverage, max_coverage)

    return name, M, domains, ps


def splice_sub_hmm(src_hmm, start, end):
    """Splice a sub-HMM from a parent by copying profile probabilities.

    Copies match emissions, insert emissions, and transition probabilities
    for positions start:end from the parent model into a new HMM object.
    The result is a fully valid pyhmmer HMM that can be searched directly.

    Works with both DNA and amino acid models.

    Args:
        src_hmm: parent pyhmmer.plan7.HMM
        start: 0-based start position in parent model
        end: 0-based end position (exclusive)

    Returns:
        pyhmmer.plan7.HMM — a new HMM covering only the specified region
    """
    sub_M = end - start
    alphabet = src_hmm.alphabet
    K = alphabet.K  # 4 for DNA, 20 for amino

    name = src_hmm.name if isinstance(src_hmm.name, str) else src_hmm.name
    sub_name = f"{name}__sub_{start}-{end}"

    sub = plan7.HMM(alphabet, M=sub_M, name=sub_name)

    for i in range(1, sub_M + 1):
        src_i = int(start + i)
        for j in range(K):
            sub.match_emissions[i][j] = float(src_hmm.match_emissions[src_i][j])
            sub.insert_emissions[i][j] = float(src_hmm.insert_emissions[src_i][j])
        for j in range(7):
            sub.transition_probabilities[i][j] = float(
                src_hmm.transition_probabilities[src_i][j])

    for j in range(7):
        sub.transition_probabilities[0][j] = float(
            src_hmm.transition_probabilities[0][j])

    sub.set_composition()
    sub.consensus = src_hmm.consensus[start:end]

    # Write to text, inject required STATS fields, reload as valid HMM.
    # pyhmmer's ProfileConfig requires STATS LOCAL MSV/VITERBI/FORWARD
    # to be present or it returns eslEINVAL.
    buf = io.BytesIO()
    sub.write(buf)
    text = buf.getvalue().decode()

    if "STATS LOCAL MSV" not in text:
        lines = text.split("\n")
        new_lines = []
        for line in lines:
            new_lines.append(line)
            if line.startswith("MAP"):
                new_lines.append("STATS LOCAL MSV       -9.9014  0.70957")
                new_lines.append("STATS LOCAL VITERBI  -10.7224  0.70957")
                new_lines.append("STATS LOCAL FORWARD   -4.1637  0.70957")
        text = "\n".join(new_lines)

    buf2 = io.BytesIO(text.encode())
    with plan7.HMMFile(buf2) as hf:
        return next(hf)


# ---------------------------------------------------------------------------
# Convenience functions — load from files
# ---------------------------------------------------------------------------

def _load_hmms(hmm_path):
    """Load all HMMs from a file. Standalone version with no dependencies."""
    import gzip as _gzip
    opener = _gzip.open if str(hmm_path).endswith(".gz") else open
    with opener(hmm_path, "rb") as fh:
        return list(plan7.HMMFile(io.BytesIO(fh.read())))


def decompose_file(hmm_path, window_size=64, max_domains=None,
                   min_score_per_base=0.0, max_overlap_frac=0.0,
                   min_coverage=0.0, max_coverage=1.0):
    """Decompose all models in an HMM database file.

    Args:
        hmm_path: path to HMM file (single or multi-model)
        window_size: sub-model window size
        max_domains: optional cap per model
        min_score_per_base: minimum avg score to accept a window
        max_overlap_frac: overlap budget between windows
        min_coverage: minimum parent model coverage fraction
        max_coverage: maximum parent model coverage fraction

    Returns:
        dict of {model_name: [(start, end, score), ...]}
    """
    results = {}
    for hmm in _load_hmms(hmm_path):
        name, M, domains, _ = decompose_model(
            hmm, window_size, max_domains, min_score_per_base,
            max_overlap_frac, min_coverage, max_coverage)
        results[name] = domains
    return results


def build_sub_hmms_from_file(hmm_path, window_size=64, max_domains=None,
                             min_score_per_base=0.0, max_overlap_frac=0.0,
                             min_coverage=0.0, max_coverage=1.0):
    """Decompose all models in a file and build searchable sub-HMMs.

    Args:
        hmm_path: path to HMM file
        window_size: sub-model window size
        max_domains: optional cap per model
        min_score_per_base: minimum avg score to accept a window
        max_overlap_frac: overlap budget between windows
        min_coverage: minimum parent model coverage fraction
        max_coverage: maximum parent model coverage fraction

    Returns:
        list of (sub_hmm, parent_name, start, end) tuples.
        Models smaller than window_size are returned as-is.
    """
    sub_hmms = []

    for hmm in _load_hmms(hmm_path):
        name = hmm.name
        M = hmm.M

        if M <= window_size:
            sub_hmms.append((hmm, name, 0, M))
            continue

        _, _, domains, _ = decompose_model(
            hmm, window_size, max_domains, min_score_per_base,
            max_overlap_frac, min_coverage, max_coverage)

        for start, end, score in domains:
            sub = splice_sub_hmm(hmm, start, end)
            sub_hmms.append((sub, name, start, end))

    return sub_hmms


# Window tiers: try largest first, fall back to smaller
AA_WINDOW_TIERS = [96, 64, 48, 32]
DNA_WINDOW_TIERS = [192, 128, 64]


def build_sub_hmms_tiered(hmm_path, max_overlap_frac=0.5,
                          min_score_per_base=0.0,
                          min_coverage=0.0, max_coverage=1.0):
    """Decompose models using tiered window sizes.

    Auto-detects alphabet: uses AA tiers (96->32) for amino acid models,
    DNA tiers (192->64) for nucleotide models. For each model, tries
    the largest window that fits (window < model.M). Models shorter
    than the smallest tier are returned whole.

    Args:
        hmm_path: path to HMM file
        max_overlap_frac: overlap budget between windows
        min_score_per_base: minimum avg score per position
        min_coverage: minimum parent model coverage fraction
        max_coverage: maximum parent model coverage fraction

    Returns:
        list of (sub_hmm, parent_name, start, end) tuples
    """
    hmms = _load_hmms(hmm_path)

    # Detect alphabet from first model
    if hmms and hmms[0].alphabet == easel.Alphabet.dna():
        window_tiers = DNA_WINDOW_TIERS
        min_size = DNA_WINDOW_TIERS[-1]  # 64
    else:
        window_tiers = AA_WINDOW_TIERS
        min_size = AA_WINDOW_TIERS[-1]  # 32

    sub_hmms = []

    for hmm in hmms:
        name = hmm.name
        M = hmm.M

        # Too small to facet
        if M < min_size:
            sub_hmms.append((hmm, name, 0, M))
            continue

        # Try each tier: largest window smaller than the model
        best_domains = None

        for window in window_tiers:
            if window >= M:
                continue

            _, _, domains, _ = decompose_model(
                hmm, window, None, min_score_per_base,
                max_overlap_frac, min_coverage, max_coverage)

            if domains:
                best_domains = domains
                break

        if best_domains is None:
            sub_hmms.append((hmm, name, 0, M))
        else:
            for start, end, score in best_domains:
                sub = splice_sub_hmm(hmm, start, end)
                sub_hmms.append((sub, name, start, end))

    return sub_hmms


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Decompose HMM models into high-scoring sub-model domains")
    parser.add_argument("hmm_file",
                        help="HMM database file (single or multi-model)")
    parser.add_argument("--window", type=int, default=64,
                        help="Sub-model window size [default: 64]")
    parser.add_argument("--overlap", type=float, default=0.0,
                        help="Max overlap fraction between windows [default: 0.0]")
    parser.add_argument("--max-domains", type=int, default=None,
                        help="Max domains per model [default: no limit]")
    parser.add_argument("--min-score-per-base", type=float, default=0.0,
                        help="Minimum avg score/position to accept [default: 0.0]")
    parser.add_argument("--min-coverage", type=float, default=0.0,
                        help="Minimum parent model coverage fraction [default: 0.0]")
    parser.add_argument("--max-coverage", type=float, default=1.0,
                        help="Maximum parent model coverage fraction [default: 1.0]")
    parser.add_argument("--output", default=None,
                        help="Output TSV [default: stdout]")
    args = parser.parse_args()

    results = []
    for hmm in _load_hmms(args.hmm_file):
        name, M, domains, pos_scores = decompose_model(
            hmm, args.window, args.max_domains, args.min_score_per_base,
            args.overlap, args.min_coverage, args.max_coverage)
        results.append((name, M, domains, pos_scores))

    out = open(args.output, "w") if args.output else sys.stdout

    out.write("model\tM\tn_domains\tdomain_pos\tcoverage\t"
              "domains\tavg_score_per_pos\n")

    total_M = 0
    total_domain = 0

    for name, M, domains, _ in sorted(results, key=lambda r: -r[1]):
        n_dom = len(domains)
        domain_pos = sum(e - s for s, e, _ in domains)
        coverage = domain_pos / M if M > 0 else 0
        total_M += M
        total_domain += domain_pos

        dom_strs = [f"{s}-{e}({sc:.1f})" for s, e, sc in domains]
        avg_spb = (sum(sc for _, _, sc in domains) / domain_pos
                   if domain_pos else 0)

        out.write(f"{name}\t{M}\t{n_dom}\t{domain_pos}\t{coverage:.2f}\t"
                  f"{','.join(dom_strs)}\t{avg_spb:.2f}\n")

    if out is not sys.stdout:
        out.close()

    n_models = len(results)
    n_domains = sum(len(d) for _, _, d, _ in results)
    print(f"\n{n_models} models -> {n_domains} domains at window={args.window}"
          f" overlap={args.overlap:.0%}",
          file=sys.stderr)
    print(f"  Total model positions: {total_M:,}", file=sys.stderr)
    print(f"  Total domain positions: {total_domain:,} "
          f"({100 * total_domain / total_M:.1f}% coverage)",
          file=sys.stderr)
    print(f"  Avg domains/model: {n_domains / n_models:.1f}", file=sys.stderr)


if __name__ == "__main__":
    main()
