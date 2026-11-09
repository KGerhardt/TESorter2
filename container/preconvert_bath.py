"""Pre-convert the bundled protein HMM databases to BATH format.

Run once at image build time, while site-packages is writable, so the
conversions sit beside the source HMMs and are reused by every run
(the image is read-only under Apptainer).
"""
import glob
import os

from tesorter2.bath_search import resolve_bath_db
from tesorter2.hmm import AMINO_ALPHABET, peek_alphabet
from tesorter2.paths import get_db_dir

for hmm in sorted(glob.glob(os.path.join(get_db_dir(), "*.hmm"))):
    if peek_alphabet(hmm) == AMINO_ALPHABET:
        print("pre-converting", os.path.basename(hmm), "->",
              resolve_bath_db(hmm), flush=True)
