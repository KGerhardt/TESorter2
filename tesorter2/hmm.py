"""
HMM database loading and optimized profile management.

Loads HMM files, auto-detects alphabet from the ALPH field,
builds OptimizedProfiles once at startup, and provides the dict-based
access needed for ad-hoc model subset searches in pass 2.
"""

import gzip
from io import BytesIO

import pyhmmer.easel as easel
import pyhmmer.plan7 as plan7


AMINO_ALPHABET = easel.Alphabet.amino()
DNA_ALPHABET = easel.Alphabet.dna()
RNA_ALPHABET = easel.Alphabet.rna()

_ALPHABET_MAP = {
    "amino": AMINO_ALPHABET,
    "DNA":   DNA_ALPHABET,
    "RNA":   RNA_ALPHABET,
}


def open_hmm_text(hmm_path):
    """Open an HMM file for text reading, transparently through gzip.

    The databases ship compressed: an HMM file is text and gives back ~84% of
    its size, and both readers that matter handle it unaided -- pyhmmer's
    HMMFile detects gzip, and HMMER 3.4's nhmmer binary reads a .gz path
    directly (verified: same model count scanned from either). So there is no
    unpacking step, nothing written into the package directory at first use,
    and no race between concurrent jobs sharing an install.
    """
    if str(hmm_path).endswith(".gz"):
        return gzip.open(hmm_path, "rt", errors="replace")
    return open(hmm_path, errors="replace")


def read_hmm_bytes(hmm_path):
    """Whole HMM file as bytes, decompressing a .gz."""
    if str(hmm_path).endswith(".gz"):
        with gzip.open(hmm_path, "rb") as fh:
            return fh.read()
    with open(hmm_path, "rb") as fh:
        return fh.read()


def peek_alphabet(hmm_path):
    """
    Detect the alphabet of an HMM database by reading the first ALPH field.

    Args:
        hmm_path: path to an HMM file

    Returns:
        easel.Alphabet instance (amino, DNA, or RNA)

    Raises:
        ValueError: if no ALPH line found or unrecognized alphabet
    """
    with open_hmm_text(hmm_path) as f:
        for line in f:
            if line.startswith("ALPH"):
                alph_str = line.split()[1].strip()
                if alph_str in _ALPHABET_MAP:
                    return _ALPHABET_MAP[alph_str]
                raise ValueError(
                    f"Unrecognized HMM alphabet '{alph_str}' in {hmm_path}"
                )
    raise ValueError(f"No ALPH line found in {hmm_path}")


def needs_translation(hmm_path):
    """
    Check whether a database requires amino acid input (i.e. translation).

    Returns True for amino acid HMMs, False for DNA/RNA.
    """
    alphabet = peek_alphabet(hmm_path)
    return alphabet == AMINO_ALPHABET


def load_hmms(hmm_path):
    """
    Load all HMM models from a file.

    Buffers the entire file into memory first for ~20x faster parsing
    vs direct file I/O (especially significant on WSL).

    Args:
        hmm_path: path to an HMM database file (.hmm)

    Returns:
        list of plan7.HMM objects
    """
    return list(plan7.HMMFile(BytesIO(read_hmm_bytes(hmm_path))))


def build_optimized_profiles(hmms, alphabet=None):
    """
    Build OptimizedProfiles from a list of HMMs.

    This is the expensive step that should happen once per worker at startup.
    Each OptimizedProfile has the SSV/MSV filter data pre-computed.

    Args:
        hmms: list of plan7.HMM objects
        alphabet: easel.Alphabet (defaults to amino)

    Returns:
        dict of {model_name: OptimizedProfile}
    """
    if alphabet is None:
        alphabet = AMINO_ALPHABET

    bg = plan7.Background(alphabet)
    profiles = {}

    for hmm in hmms:
        profile = plan7.Profile(hmm.M, alphabet)
        profile.configure(hmm, bg)
        opt = profile.to_optimized()
        profiles[hmm.name] = opt

    return profiles


def load_and_optimize(hmm_path, alphabet=None):
    """
    Load HMMs from file and build optimized profiles in one step.

    Convenience function for worker initialization. Auto-detects
    alphabet if not provided.

    Args:
        hmm_path: path to an HMM database file
        alphabet: easel.Alphabet (auto-detected if None)

    Returns:
        tuple of (list of HMMs, dict of {model_name: OptimizedProfile},
                  easel.Alphabet)
    """
    if alphabet is None:
        alphabet = peek_alphabet(hmm_path)
    hmms = load_hmms(hmm_path)
    profiles = build_optimized_profiles(hmms, alphabet)
    return hmms, profiles, alphabet
