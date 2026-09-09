"""
classifier.py — Config-driven TE classification from HMM domain hits.

Replicates TEsorter's classification logic: domain deconfliction (hmm2best),
order/superfamily/clade assignment, completeness checking. Parameterized
by per-database configs rather than hardcoded per-database classes.

Two stages:
  1. hmm2best: select best domain hit per (sequence, domain_type),
     with database-specific remapping and overlap rules
  2. classify: assign order/superfamily/clade from domain architecture
"""

import logging
import os
import numpy as np
from collections import Counter, defaultdict

from .clade_harmonize import (load_harmonization, harmonize, harmonize_lineage,
                              load_scope, resolves)

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Database configs
# ---------------------------------------------------------------------------

REXDB_CONFIG = {
    "name": "rexdb",
    "domain_remap": {"aRH": "RH", "TPase": "INT"},
    "overlap_aware": True,
    "structures": {
        ("LTR", "Copia"):   ["GAG", "PROT", "INT", "RT", "RH"],
        ("LTR", "Gypsy"):   ["GAG", "PROT", "RT", "RH", "INT"],
        ("LTR", "Bel-Pao"): ["GAG", "PROT", "RT", "RH", "INT"],
    },
    "clade_parser": "rexdb",
    "clade_restrict": {"Copia", "Gypsy"},  # only these get clade names
}

GYDB_CONFIG = {
    "name": "gydb",
    "domain_remap": {},
    "overlap_aware": False,
    "structures": {
        ("LTR", "Copia"):          ["GAG", "AP", "INT", "RT", "RNaseH"],
        ("LTR", "Gypsy"):          ["GAG", "AP", "RT", "RNaseH", "INT"],
        ("LTR", "Pao"):            ["GAG", "AP", "RT", "RNaseH", "INT"],
        ("LTR", "Retroviridae"):   ["GAG", "AP", "RT", "RNaseH", "INT", "ENV"],
        ("LTR", "Caulimoviridae"): ["GAG", "AP", "RT", "RNaseH"],
    },
    "clade_parser": "gydb",
    "clade_restrict": None,
}

LINE_CONFIG = {
    "name": "line",
    "domain_remap": {},
    "overlap_aware": False,
    "structures": {},
    "clade_parser": "rexdb",
    "clade_restrict": None,
}

TIR_CONFIG = {
    "name": "tir",
    "domain_remap": {},
    "overlap_aware": False,
    "structures": {},
    "clade_parser": "rexdb",
    "clade_restrict": None,
}

SINE_CONFIG = {
    "name": "sine",
    "domain_remap": {},
    "overlap_aware": False,
    "structures": {},
    "clade_parser": "sine",
    "clade_restrict": None,
}

# Pfam-derived TE domains. The model names are written in the REXdb convention
# by build_pfam_te_db.py, so the REXdb parser reads them unchanged. No
# `structures`: Pfam families are single domains and the add-on makes no claim
# about domain architecture, so completeness is never asserted from it.
PFAMTE_CONFIG = {
    "name": "pfam-te",
    "domain_remap": {},
    "overlap_aware": False,
    "structures": {},
    "clade_parser": "rexdb",
    "clade_restrict": None,
}

# Dfam curated DNA consensus models, cut into cost tiers by build_dfam_dbs.py,
# which writes their names in the REXdb convention exactly as the Pfam TE
# add-on does -- so the REXdb parser reads them unchanged. No `structures`: a
# whole-element consensus is not a domain architecture, and these models make
# no claim about which domains an element carries. One config serves all four
# collections; they are disjoint slices of one library, not different
# databases.
DFAMDB_CONFIG = {
    "name": "dfam",
    "domain_remap": {},
    "overlap_aware": False,
    "structures": {},
    "clade_parser": "rexdb",
    "clade_restrict": None,
}

DB_CONFIGS = {
    "rexdb": REXDB_CONFIG,
    "gydb": GYDB_CONFIG,
    "line": LINE_CONFIG,
    "tir": TIR_CONFIG,
    "sine": SINE_CONFIG,
    "sine-animals": SINE_CONFIG,
    "sine-so": SINE_CONFIG,
    "pfam-te": PFAMTE_CONFIG,
    "dfam-core": DFAMDB_CONFIG,
    "dfam-extended": DFAMDB_CONFIG,
    "dfam-deep": DFAMDB_CONFIG,
    "dfam-complete": DFAMDB_CONFIG,
}

# GyDB clade map: loaded from GyDB2.hmm.info (ships alongside GyDB2.hmm)
from .paths import get_db_dir
from .so_map import lookup as so_lookup
_GYDB_INFO_PATH = os.path.join(get_db_dir(), "GyDB2.hmm.info")


def load_gydb_clade_map(info_path=None):
    """Load GyDB clade -> (order, superfamily) map from .info file.

    Replicates TEsorter's CladeInfo parser.
    """
    if info_path is None:
        info_path = _GYDB_INFO_PATH

    if not os.path.exists(info_path):
        log.warning(f"GyDB info file not found: {info_path}")
        return {}

    clade_aliases = {
        'Ty_(Pseudovirus)': 'pseudovirus',
        'Cer2-3': 'cer2-3',
        '412/Mdg1': '412_mdg1',
        'TF1-2': 'TF',
        'Micropia/Mdg3': 'micropia_mdg3',
        'CoDi-I': 'codi_I',
        'CoDi-II': 'codi_II',
        '17.6': '17_6',
    }

    clade_map = {}
    header = None

    for line in open(info_path):
        parts = line.strip().split('\t')
        if header is None:
            header = parts
            continue

        d = dict(zip(header, parts))
        clade = d.get('Clade', 'NA')
        if clade == 'NA':
            clade = d.get('Cluster_or_genus', '')

        superfamily = d.get('Family', '').split('/')[-1]
        if superfamily == 'Retroviridae':
            clade = d.get('Cluster_or_genus', '').replace('virus', 'viridae')
        if superfamily == 'Retrovirus':
            superfamily = 'Retroviridae'

        system = d.get('System', '')
        order = 'LTR' if system in {'LTR_retroelements', 'LTR_Retroelements',
                                      'LTR_retroid_elements'} else system

        # Register clade and all its aliases
        clade_map[clade] = (order, superfamily)
        if clade in clade_aliases:
            alias = clade_aliases[clade]
            clade_map[alias] = (order, superfamily)
            if alias == '412_mdg1':
                clade_map['412-mdg1'] = (order, superfamily)

        # Also register with underscores and lowercase
        clade_map[clade.replace('-', '_')] = (order, superfamily)
        clade_map[clade.lower()] = (order, superfamily)
        clade_map[clade.replace('-', '_').lower()] = (order, superfamily)

    # Special entries
    clade_map['ty1/copia'] = ('LTR', 'Copia')

    # Unknown clades
    for clade, order in {
        'retroelement': 'LTR', 'retroviridae': 'LTR',
        'B-type_betaretroviridae': 'LTR', 'D-type_betaretroviridae': 'LTR',
        'caulimoviruses': 'LTR', 'caulimoviridae_dom2': 'LTR',
        'errantiviridae': 'LTR', 'retropepsins': 'LTR',
        'VPX_retroviridae': 'LTR', 'cog5550': 'Unknown',
        'ddi': 'Unknown', 'dtg_ilg_template': 'Unknown',
        'saspase': 'Unknown', 'GIN1': 'Unknown',
        'shadow': 'Unknown', 'all': 'Unknown',
        'pepsins_A1a': 'Unknown', 'pepsins_A1b': 'Unknown',
    }.items():
        clade_map[clade] = (order, 'unknown')

    return clade_map


# ---------------------------------------------------------------------------
# Clade parsing
# ---------------------------------------------------------------------------

def parse_clade_rexdb(model_name):
    """Parse REXdb model name -> (order, superfamily, clade, gene).

    Format: Class_I/LTR/Ty1_copia/SIRE:Ty1-RT
    """
    if ":" not in model_name:
        return "Unknown", "unknown", model_name, model_name

    clade_path, domain = model_name.split(":", 1)
    gene = domain.split("-")[-1]

    if clade_path.startswith("Class_I/LTR/Ty1_copia"):
        order, superfamily = "LTR", "Copia"
    elif clade_path.startswith("Class_I/LTR/Ty3_gypsy"):
        order, superfamily = "LTR", "Gypsy"
    elif clade_path.startswith("Class_I/LTR/"):
        parts = clade_path.split("/")
        order, superfamily = parts[1], parts[2] if len(parts) > 2 else "unknown"
    elif clade_path.startswith("Class_I/"):
        parts = clade_path.split("/")
        order = parts[1]
        superfamily = parts[2] if len(parts) > 2 else "unknown"
    elif clade_path.startswith("Class_II/"):
        parts = clade_path.split("/")
        order = parts[2] if len(parts) > 2 else "unknown"
        superfamily = parts[3] if len(parts) > 3 else "unknown"
    elif clade_path.startswith("NA"):
        order, superfamily = "LTR", "Retrovirus"
    else:
        order, superfamily = "Unknown", "unknown"

    clade = clade_path.split("/")[-1] if "/" in clade_path else clade_path

    return order, superfamily, clade, gene


def parse_clade_gydb(model_name):
    """Parse GyDB model name -> (gene, clade).

    Format: AP_copia, RT_gypsy, GAG_lentiviridae
    """
    parts = model_name.split("_", 1)
    gene = parts[0]
    clade = parts[1] if len(parts) > 1 else model_name
    return gene, clade


# ---------------------------------------------------------------------------
# hmm2best: domain deconfliction
# ---------------------------------------------------------------------------

def hmm2best(hits, config, compat_rounding=False):
    """Select best domain hits per (base_seq, domain_type).

    Args:
        hits: dict from deconflict.load_hits()
        config: database config dict
        compat_rounding: if True, round norm_score to 2 decimal places
                         (replicates TEsorter rounding bug)

    Returns:
        numpy index array of selected hits
    """
    n = len(hits["score"])
    if n == 0:
        return np.array([], dtype=int)

    remap = config["domain_remap"]
    overlap_aware = config["overlap_aware"]

    # Apply domain remapping
    domain_types = hits["family"].copy()
    for old, new in remap.items():
        mask = domain_types == old
        domain_types[mask] = new

    # Normalized score
    norm_score = hits["score"] / hits["model_len"]
    if compat_rounding:
        norm_score = np.round(norm_score, 2)

    # Group by (base_seq, domain_type)
    # Process in score-descending order within each group
    best = {}  # (base_seq, domain_type) -> index

    sort_idx = np.argsort(-norm_score)

    for i in sort_idx:
        base = hits["base_seq"][i]
        dtype = domain_types[i]
        key = (base, dtype)

        if key not in best:
            best[key] = i
            continue

        if not overlap_aware:
            # Simple: already have a better score, skip
            continue

        # Overlap-aware logic (REXdb)
        curr = best[key]
        curr_score = norm_score[curr]
        new_score = norm_score[i]

        if new_score <= curr_score:
            continue

        # Only replace if same gene OR overlapping envelopes
        curr_model = hits["model"][curr]
        new_model = hits["model"][i]

        # Same gene = same domain prefix (e.g. both Ty1-RT)
        if ":" in curr_model and ":" in new_model:
            curr_gene = curr_model.split(":")[1]
            new_gene = new_model.split(":")[1]
            same_gene = curr_gene == new_gene
        else:
            same_gene = curr_model.split("_")[0] == new_model.split("_")[0]

        if same_gene:
            best[key] = i
            continue

        # Check envelope overlap
        curr_start, curr_end = hits["env_from"][curr], hits["env_to"][curr]
        new_start, new_end = hits["env_from"][i], hits["env_to"][i]
        if new_start <= curr_end and new_end >= curr_start:
            best[key] = i

    return np.array(list(best.values()), dtype=int)


def apply_filters(hits, indices, min_cov=20.0, max_evalue=1e-3,
                  min_acc=0.5, min_norm_score=0.1, compat_rounding=False):
    """Apply filters to selected hits.

    Uses full precision by default. --compat-tesorter-rounding rounds
    norm_score to 2 decimal places before threshold comparison.
    """
    idx = np.array(indices)
    cov = hits["hmm_cov"][idx]
    evalue = hits["evalue"][idx]
    acc = hits["acc"][idx]
    nscore = hits["score"][idx] / hits["model_len"][idx]
    if compat_rounding:
        nscore = np.round(nscore, 2)

    mask = (cov >= min_cov) & (evalue <= max_evalue) & (acc >= min_acc) & (nscore >= min_norm_score)
    return idx[mask]


def select_domain_indices(hits, config, compat_rounding=False):
    """Indices of the domain assignments a classification run actually uses.

    This is steps 1-2 of classify_sequences (hmm2best, then apply_filters).
    Exposed so the TEsorter-format domain exports can report the same one
    best model per region that the classification did, rather than every
    competing model in the raw hit table.
    """
    best_idx = hmm2best(hits, config, compat_rounding=compat_rounding)
    return apply_filters(hits, best_idx, compat_rounding=compat_rounding)


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

def classify_element(genes, clades, models, scores, config, compat_voting=False,
                     model_lens=None, min_clade_delta=0.0):
    """Classify a single TE element from its domain hits.

    Args:
        genes: list of gene names (domain types) in positional order
        clades: list of clade names corresponding to each gene
        models: list of full model names
        scores: list of per-domain normalized scores (dom_score / model_len),
                aligned with genes/clades/models. Used for score-weighted
                clade voting (the default). Ignored when compat_voting=True.
        config: database config dict
        compat_voting: if True, decide the clade by raw domain-count plurality
                       (exact TEsorter replication). Default False uses the
                       score-weighted vote.

    Returns:
        (order, superfamily, clade, complete, clade_delta)

        clade_delta is the offset margin from clade_delta_star -- always
        computed when model_lens is given, and reported whether or not
        min_clade_delta is in use.
    """
    parser = config["clade_parser"]

    if parser == "rexdb":
        return _classify_rexdb(genes, clades, models, scores, config,
                               compat_voting=compat_voting,
                               model_lens=model_lens,
                               min_clade_delta=min_clade_delta)
    elif parser == "gydb":
        return _classify_gydb(genes, clades, models, scores, config,
                              compat_voting=compat_voting,
                              model_lens=model_lens,
                              min_clade_delta=min_clade_delta)
    elif parser == "sine":
        return "SINE", "unknown", "unknown", "unknown", None
    else:
        return "Unknown", "unknown", "unknown", "unknown", None


def clade_delta_star(clades, scores, model_lens):
    """Smallest per-hit bit offset that would change the clade winner.

    The clade vote is an argmax over summed *normalized* score,
    W(c) = sum over that clade's hits of score/model_len. A per-hit offset of
    delta bits therefore does NOT shift every clade equally -- it is divided by
    model length, so

        W_delta(c) = W(c) + delta * K(c),   K(c) = sum of 1/model_len

    A clade carried by three short models gains far more from +delta than one
    carried by a single long model. That asymmetry is why a systematic
    inter-engine score offset moves clade calls at all: BATH scores ~5.6 bits
    above pyhmmer on shared pairs, which is enough to reorder GyDB clades (314
    models, clade encoded in the model name) while leaving superfamily intact.

    For the winner w and any rival c, the two are level when

        delta = (W(w) - W(c)) / (K(c) - K(w))

    and this returns the smallest such |delta| over all rivals: the largest
    offset the call can absorb in either direction and still stand. Big means
    the call is robust to rescoring, a different engine, or a database
    revision; small means it was decided by a margin no wider than the noise
    between engines.

    Returns inf when only one clade is in play (nothing to flip), and None if
    model lengths were not supplied.
    """
    if model_lens is None:
        return None
    weights = defaultdict(float)
    sens = defaultdict(float)
    for c, s, ln in zip(clades, scores, model_lens):
        weights[c] += s
        if ln:
            sens[c] += 1.0 / float(ln)
    if len(weights) < 2:
        return float("inf")
    win = max(weights, key=lambda x: weights[x])
    best = float("inf")
    for c in weights:
        if c == win:
            continue
        dk = sens[c] - sens[win]
        if dk == 0:
            # Same offset sensitivity: no uniform offset can ever reorder them.
            continue
        best = min(best, abs((weights[win] - weights[c]) / dk))
    return best


def _clade_winner(clades, scores, compat_voting, compat_rule="rexdb"):
    """Pick the winning clade and report whether it is an unambiguous winner.

    Default (score-weighted): sum each clade's per-domain normalized scores and
    take the argmax. The winner is "clear" unless the top two clades tie
    exactly on summed score (vanishingly rare with continuous scores). This
    resolves the TEsorter clade-swap artifact, where sibling clades draw an
    equal number of domain votes and the tie is broken by dict iteration order.

    compat_voting=True restores TEsorter's behaviour: plurality by raw domain
    count. The "clear" rule differs per parser, matching the original code:
      - rexdb: clear when single clade, or top count > 1 AND top count strictly
        exceeds the second count (compared in insertion order).
      - gydb:  clear when single clade, or top count > 1 (no second comparison).

    Returns (max_clade, n_distinct, clear_winner).
    """
    if compat_voting:
        clade_count = Counter(clades)
        max_clade = max(clade_count, key=lambda x: clade_count[x])
        counts = list(clade_count.values())
        if compat_rule == "gydb":
            clear = len(clade_count) == 1 or clade_count[max_clade] > 1
        else:
            clear = len(clade_count) == 1 or (
                clade_count[max_clade] > 1 and len(counts) > 1 and counts[0] > counts[1])
        return max_clade, len(clade_count), clear

    weights = defaultdict(float)
    for c, s in zip(clades, scores):
        weights[c] += s
    max_clade = max(weights, key=lambda x: weights[x])
    ordered = sorted(weights.values(), reverse=True)
    clear = len(weights) == 1 or (len(ordered) > 1 and ordered[0] > ordered[1])
    return max_clade, len(weights), clear


def _demote(display_clade, delta, min_clade_delta):
    """Blank a clade whose vote a small score offset could overturn.

    Applied to the clade string ONLY, after classification is complete.
    Routing it through the existing `clear` flag instead would be wrong: that
    flag also propagates to superfamily and order when a sequence's domains
    span several, so a clade-margin threshold would silently degrade
    Gypsy -> mixture at superfamily and LTR -> mixture at order. Measured on a
    400-sequence GyDB sample, 6 of 42 demotions moved superfamily and 3 moved
    order. A margin on the clade vote is evidence about the clade and nothing
    above it.
    """
    if min_clade_delta <= 0 or delta is None:
        return display_clade
    if delta >= min_clade_delta:
        return display_clade
    if display_clade in ("mixture", "unknown"):
        return display_clade
    return "mixture"


def _classify_rexdb(genes, clades, models, scores, config, compat_voting=False,
                    model_lens=None, min_clade_delta=0.0):
    """REXdb classification logic."""
    max_clade, n_distinct, clear = _clade_winner(clades, scores, compat_voting)
    delta = None if compat_voting else clade_delta_star(clades, scores, model_lens)

    order, superfamily, _, _ = parse_clade_rexdb(
        [m for m, c in zip(models, clades) if c == max_clade][0])

    if clear:
        display_clade = max_clade.split("/")[-1]
    elif n_distinct > 1:
        display_clade = "mixture"
        superfamilies = [parse_clade_rexdb(m)[1] for m in models]
        if len(Counter(superfamilies)) > 1:
            superfamily = "mixture"
            orders = [parse_clade_rexdb(m)[0] for m in models]
            if len(Counter(orders)) > 1:
                order = "mixture"
    else:
        display_clade = max_clade.split("/")[-1]

    # Check completeness
    structures = config["structures"]
    try:
        expected = structures[(order, superfamily)]
        present = [g for g in genes if g in set(expected)]
        complete = "yes" if expected == present else "no"
    except KeyError:
        complete = "unknown"

    # Restrict clade names
    restrict = config.get("clade_restrict")
    if restrict and superfamily not in restrict:
        display_clade = "unknown"
    if display_clade.startswith("Ty"):
        display_clade = "unknown"

    display_clade = _demote(display_clade, delta, min_clade_delta)
    return order, superfamily, display_clade, complete, delta


def _classify_gydb(genes, clades, models, scores, config, compat_voting=False,
                   model_lens=None, min_clade_delta=0.0):
    """GyDB classification logic. Requires clade_map."""
    max_clade, n_distinct, clear = _clade_winner(
        clades, scores, compat_voting, compat_rule="gydb")
    delta = None if compat_voting else clade_delta_star(clades, scores, model_lens)

    # Look up order/superfamily from clade map
    clade_map = config.get("_clade_map", {})
    order, superfamily = clade_map.get(max_clade, ("Unknown", "unknown"))

    if clear:
        display_clade = max_clade
    elif n_distinct > 1:
        display_clade = "mixture"
        superfamilies = [clade_map.get(c, [None, None])[1] for c in clades]
        if len(Counter(superfamilies)) > 1:
            superfamily = "mixture"
            orders = [clade_map.get(c, [None, None])[0] for c in clades]
            if len(Counter(orders)) > 1:
                order = "mixture"
    else:
        display_clade = max_clade

    structures = config["structures"]
    try:
        expected = structures[(order, superfamily)]
        present = [g for g in genes if g in set(expected)]
        complete = "yes" if expected == present else "no"
    except KeyError:
        complete = "unknown"

    display_clade = _demote(display_clade, delta, min_clade_delta)
    return order, superfamily, display_clade, complete, delta


# ---------------------------------------------------------------------------
# Full classification pipeline
# ---------------------------------------------------------------------------

def classify_sequences(hits, config, gydb_clade_map=None, compat_rounding=False,
                       compat_voting=False, min_clade_delta=0.0):
    """Full classification: hmm2best -> filter -> classify per sequence.

    Args:
        hits: dict from deconflict.load_hits()
        config: database config dict
        gydb_clade_map: optional {clade: (order, superfamily)} for GyDB
        compat_voting: if True, decide the clade by raw domain-count plurality
                       (exact TEsorter replication) instead of the default
                       score-weighted vote.

    Returns:
        list of dicts with keys: id, order, superfamily, clade, complete,
        strand, domains
    """
    if config["clade_parser"] == "gydb":
        config = dict(config)
        if gydb_clade_map:
            config["_clade_map"] = gydb_clade_map
        else:
            config["_clade_map"] = load_gydb_clade_map()

    # Step 1: hmm2best
    best_idx = hmm2best(hits, config, compat_rounding=compat_rounding)
    log.info(f"  hmm2best: {len(best_idx)} domain assignments")

    # Step 2: filter
    filtered_idx = apply_filters(hits, best_idx, compat_rounding=compat_rounding)
    log.info(f"  After filter: {len(filtered_idx)} assignments")

    # Step 3: group by base_seq and classify
    seq_domains = defaultdict(list)
    for i in filtered_idx:
        base = hits["base_seq"][i]
        seq_domains[base].append(i)

    results = []
    for base_seq, indices in seq_domains.items():
        # Determine strand
        # Use frame info from target names
        strands = []
        for i in indices:
            target = hits["target"][i]
            if "|fwd" in target:
                strands.append("+")
            elif "|rev" in target:
                strands.append("-")
            else:
                strands.append(".")

        unique_strands = set(strands)
        if len(unique_strands) > 1:
            strand = "?"
        elif len(unique_strands) == 1:
            strand = strands[0]
        else:
            continue

        # Sort by position (env_from)
        sorted_indices = sorted(indices, key=lambda i: hits["env_from"][i])
        if strand == "-":
            sorted_indices.reverse()

        # Extract genes, clades, models
        remap = config["domain_remap"]
        genes = []
        clades = []
        models = []
        scores = []
        model_lens = []
        domain_strs = []

        for i in sorted_indices:
            model = hits["model"][i]
            if config["clade_parser"] == "rexdb":
                _, _, clade, gene = parse_clade_rexdb(model)
            elif config["clade_parser"] == "gydb":
                gene, clade = parse_clade_gydb(model)
            else:
                gene, clade = "SINE", "SINE"

            # Apply remapping for display
            display_gene = remap.get(gene, gene)
            genes.append(display_gene)
            clades.append(clade)
            models.append(model)
            scores.append(float(hits["norm_score"][i]))
            model_lens.append(float(hits["model_len"][i]))
            domain_strs.append(f"{display_gene}|{clade}")

        order, superfamily, max_clade, complete, clade_delta = classify_element(
            genes, clades, models, scores, config, compat_voting=compat_voting,
            model_lens=model_lens, min_clade_delta=min_clade_delta)

        total_norm_score = float(np.sum(hits["norm_score"][sorted_indices]))

        results.append({
            "id": base_seq,
            "order": order,
            "superfamily": superfamily,
            "clade": max_clade,
            "complete": complete,
            "strand": strand,
            "domains": " ".join(domain_strs),
            "score": total_norm_score,
            "clade_delta": clade_delta,
        })

    log.info(f"  Classified: {len(results)} sequences")
    return results


def store_classifications(conn, results, database=None, mode="default",
                          engine=None, stage=0):
    """Store classification results in SQLite.

    The mode column records which search path produced the call
    ('default' or 'hierarchical'), so runs can coexist in a single
    companion database.

    engine/stage carry cascade provenance: which search tool produced the
    evidence behind this call and at which cascade position it was resolved.
    Both are per-call rather than per-run because one database's sequences may
    be resolved by different engines at different stages.
    """
    # Indexes on this table are built by results.finalize_db at the end
    # of the pipeline, not here — this function is called once per
    # database inside the classification loop.
    conn.execute("""
        CREATE TABLE IF NOT EXISTS classifications (
            seq_id      TEXT NOT NULL,
            database    TEXT,
            te_order    TEXT NOT NULL,
            superfamily TEXT NOT NULL,
            clade       TEXT NOT NULL,
            complete    TEXT NOT NULL,
            strand      TEXT NOT NULL,
            domains     TEXT,
            source      TEXT NOT NULL DEFAULT 'hmm',
            mode        TEXT NOT NULL DEFAULT 'default',
            engine      TEXT,
            stage       INTEGER NOT NULL DEFAULT 0
        )
    """)

    rows = [(r["id"], database, r["order"], r["superfamily"], r["clade"],
             r["complete"], r["strand"], r["domains"],
             r.get("blast_source", "hmm") if "blast_source" in r else "hmm",
             mode,
             r.get("engine", engine), r.get("stage", stage))
            for r in results]

    conn.executemany(
        "INSERT INTO classifications (seq_id, database, te_order, superfamily, "
        "clade, complete, strand, domains, source, mode, engine, stage) "
        "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        rows,
    )
    conn.commit()


def export_classification_tsv(results, out_path, include_secondary=False,
                              include_so=False, include_lineage=False,
                              include_clade_delta=False):
    """Export classification results as TSV.

    Default format matches TEsorter cls.tsv (7 columns). If
    include_secondary=True, appends a SecondaryHits column containing
    per-database classifications and evidence scores in descending order.
    If include_so=True, appends the Sequence Ontology term and accession for
    the Order/Superfamily call. If include_lineage=True, appends the sub-clade
    Lineage ('.' when no database resolved one). If include_clade_delta=True,
    appends CladeDelta -- the clade_delta_star offset margin, in bits: how large
    a per-hit score offset the clade call can absorb before the winner changes.
    'inf' means only one clade was in play; '.' means it was not computed.
    All optional columns are
    appended, leaving the positions of the original seven untouched -- EDTA and
    friends read this file positionally as (id, order, superfamily).
    """
    columns = ["#TE", "Order", "Superfamily", "Clade", "Complete",
               "Strand", "Domains"]
    if include_secondary:
        columns.append("SecondaryHits")
    if include_so:
        columns += ["SO_name", "SO_ID"]
    if include_lineage:
        columns.append("Lineage")
    if include_clade_delta:
        columns.append("CladeDelta")
    with open(out_path, "w") as f:
        f.write("\t".join(columns) + "\n")
        for r in results:
            line = [r["id"], r["order"], r["superfamily"], r["clade"],
                    r["complete"], r["strand"], r["domains"]]
            if include_secondary:
                secondary = r.get("secondary") or []
                if secondary:
                    sec_str = ";".join(
                        f"{db}:{o}/{sf}/{cl}={sc:.3f}"
                        for db, o, sf, cl, sc in secondary
                    )
                else:
                    sec_str = "."
                line.append(sec_str)
            if include_so:
                so_name, so_id = so_lookup(r["order"], r["superfamily"])
                line += [so_name, so_id]
            if include_lineage:
                line.append(r.get("lineage") or ".")
            if include_clade_delta:
                d = r.get("clade_delta")
                # "inf" = only one clade in play, nothing to flip; "." = not
                # computed (compat voting, or no model lengths available).
                line.append("." if d is None
                            else ("inf" if d == float("inf") else f"{d:.4f}"))
            f.write("\t".join(line) + "\n")


# A superfamily a database could not resolve. Not a competing claim -- the
# absence of one -- so it does not enter the vote against a database that named
# something. See _scope_aware_superfamily.
UNRESOLVED_SUPERFAMILY = {"unknown", "none", "na", "n/a", ".", "-", ""}


def _scope_aware_superfamily(pool, weighted_winner, scope,
                             unknown_competes=False):
    """Pick the winning (harmonized) superfamily with scope-aware deferral.

    **An unresolved superfamily does not compete with a named one.** REXdb
    carries one catch-all model per non-LTR order, so a LINE element it detects
    comes back ``LINE/unknown`` -- a claim about the order with no claim about
    the superfamily. Scored against the LINE database's ``LINE/L1``, REXdb's
    deeper domain coverage wins the summed-score vote and the element loses its
    superfamily entirely. On RepBase 31.07 that is 8,123 of the 8,264 sequences
    where the reconciled superfamily disagreed with the curated label, and it
    only happens in multi-database runs -- the default. The scope rule below was
    meant to catch exactly this, but `resolves()` returns False for any
    (db, superfamily) missing from clade_scope.tsv, which has rows for gydb and
    rexdb only, so neither side ever qualified and it fell through to score.

    Dropping unresolved entrants is done before the vote, so it applies to the
    plain score fallback too. If every entrant is unresolved, they all stay and
    the answer is unresolved, which is correct. Pass unknown_competes=True to
    restore the previous behaviour.

    Generalizes the paper's two-database superfamily fusion to N databases and
    reduces to it exactly for two. When the databases disagree on superfamily,
    a database's claim is deferred to iff (a) no OTHER participating database
    can resolve that superfamily (it is uncontestable -- only this database can
    express it), and (b) the claiming database can resolve at least one of the
    competing superfamilies (it is competent about the alternatives, so its
    different call is informative). If exactly one superfamily qualifies, it
    wins; otherwise fall back to the plain summed-score vote (``base``).
    """
    if not unknown_competes:
        named = [e for e in pool
                 if str(e[2]["superfamily"]).strip().lower()
                 not in UNRESOLVED_SUPERFAMILY]
        if named:
            pool = named

    sfs = {e[2]["superfamily"] for e in pool}
    base = weighted_winner(pool, "superfamily")
    if len(sfs) == 1:
        return base

    qualifiers = set()
    for db, _r, h in pool:
        s = h["superfamily"]
        others = {e[0] for e in pool if e[0] != db}
        contestable = any(resolves(od, s, scope) for od in others)
        if contestable:
            continue
        competent = any(resolves(db, cs, scope) for cs in sfs if cs != s)
        if competent:
            qualifiers.add(s)

    return qualifiers.pop() if len(qualifiers) == 1 else base


def reconcile_classifications(per_db_results, harmonization=None, scope=None,
                              unknown_competes=False):
    """Reconcile per-database classifications via hierarchical weighted vote.

    For each sequence classified by multiple databases, vote at each taxonomic
    level in sequence (order -> superfamily -> clade), weighted by summed
    normalized domain score per database. At each level, retain only entries
    agreeing with the winning label before voting at the next level. The
    primary database is the highest-scoring entry consistent with the full
    hierarchical winner. All per-database calls are returned as secondary
    hits sorted by evidence strength.

    REXdb and GyDB use different nomenclatures, so a raw string vote can only
    pool evidence at Order. When >= 2 databases classify the same element, the
    vote is therefore run on *harmonized* labels (see clade_harmonize): GyDB's
    ``Pao``/``athila`` and REXdb's ``Bel-Pao``/``Athila`` collapse onto a shared
    taxonomy, so the databases genuinely pool wherever they agree, and the
    primary call is reported with harmonized names.

    The superfamily decision is additionally **scope-aware**: when databases
    disagree on superfamily, the vote defers to the database whose claimed
    superfamily the others cannot resolve at lineage level (a catch-all branch),
    provided that database can resolve the competing superfamilies -- so a
    database does not win a superfamily it can only catch-all detect against one
    that genuinely models it (see _scope_aware_superfamily / clade_scope.tsv).

    Elements classified by a single database are left on their native names
    (strict no-op), so single-db runs are unaffected. The ``secondary`` audit
    list always keeps native names.

    Args:
        per_db_results: dict of {db_name: [result_dict, ...]} where each
            result dict comes from classify_sequences and contains at least
            id, order, superfamily, clade, score.
        harmonization: optional pre-loaded harmonization table (for testing);
            defaults to clade_harmonize.load_harmonization().
        scope: optional pre-loaded scope mask (for testing); defaults to
            clade_harmonize.load_scope().

    Returns:
        list of dicts — one per sequence — with the primary result's fields
        plus 'secondary': list of (db, order, superfamily, clade, score)
        tuples (native names) sorted by score descending across every database
        that classified this sequence.
    """
    if harmonization is None:
        harmonization = load_harmonization()
    if scope is None:
        scope = load_scope()

    by_seq = defaultdict(list)
    for db_name, results in per_db_results.items():
        for r in results:
            by_seq[r["id"]].append((db_name, r))

    def weighted_winner(entries, level):
        # entries are (db, native_result, harmonized_labels); vote on the
        # harmonized label so cross-database calls pool where they agree.
        totals = defaultdict(float)
        for _, r, h in entries:
            totals[h[level]] += r["score"]
        return max(totals.items(), key=lambda kv: kv[1])[0]

    reconciled = []
    for seq_id, entries in by_seq.items():
        multi_db = len({db for db, _ in entries}) >= 2

        # Harmonize labels only when >= 2 databases weigh in; otherwise keep
        # native labels so single-database output is byte-identical to before.
        annotated = []
        for db, r in entries:
            if multi_db:
                o_h, sf_h, c_h = harmonize(
                    db, r["order"], r["superfamily"], r["clade"], harmonization)
                l_h = harmonize_lineage(
                    db, r["order"], r["superfamily"], r["clade"], harmonization)
            else:
                o_h, sf_h, c_h = r["order"], r["superfamily"], r["clade"]
                l_h = ""
            annotated.append((db, r, {"order": o_h, "superfamily": sf_h,
                                      "clade": c_h, "lineage": l_h}))

        pool = annotated
        for level in ("order", "superfamily", "clade"):
            if level == "superfamily" and multi_db:
                winner = _scope_aware_superfamily(pool, weighted_winner, scope,
                                                unknown_competes=unknown_competes)
            else:
                winner = weighted_winner(pool, level)
            pool = [e for e in pool if e[2][level] == winner]

        primary_db, primary_r, primary_h = max(pool, key=lambda e: e[1]["score"])

        # Lineage is deliberately NOT taken from the primary entry. Within the
        # winning clade only some databases resolve a sub-lineage (REXdb names
        # Retand where GyDB can only say Tat), so reading it off the highest
        # scoring entry would silently drop the finer call whenever the coarser
        # database won on score -- the resolution loss this column exists to
        # prevent. Vote among the entries that do resolve one instead.
        resolved = [e for e in pool if e[2]["lineage"]]
        lineage = weighted_winner(resolved, "lineage") if resolved else ""

        secondary = [
            (db, r["order"], r["superfamily"], r["clade"], r["score"])
            for db, r in sorted(entries, key=lambda e: -e[1]["score"])
        ]

        # In multi-db mode report the harmonized consensus labels; the pool
        # filter guarantees primary_h is the hierarchical winner triple.
        reconciled.append({
            **primary_r,
            "order": primary_h["order"],
            "superfamily": primary_h["superfamily"],
            "clade": primary_h["clade"],
            "lineage": lineage,
            "primary_db": primary_db,
            "secondary": secondary,
        })

    return reconciled
