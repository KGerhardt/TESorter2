#!/bin/bash
#SBATCH --job-name=t2_regress
#SBATCH --account=bio250374
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/regress_%j.out
#SBATCH --error=logs/regress_%j.err

# Regression check: does the packaged tesorter2 (src/ -> tesorter2/ package,
# relative imports, renamed CLI) produce byte-identical classifications to the
# pre-refactor code at 59e1af8?
set -euo pipefail

module load python/3.9.5

# Locate the suite rather than hardcoding it: the script's own directory when
# run directly, and under sbatch -- where the script executes from a spool copy
# -- SLURM_SUBMIT_DIR, so submit from _regression/. Each candidate is checked
# for a marker file, because SLURM_SUBMIT_DIR is also set inside interactive
# salloc/srun shells, where it points at wherever that shell was started rather
# than at this script. TESORTER2_DEV overrides everything. (The relative
# --output paths above resolve against the sbatch working directory.)
REG=""
for _c in "${TESORTER2_DEV:+$TESORTER2_DEV/_regression}" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)" \
          "${SLURM_SUBMIT_DIR:-}"; do
    [[ -n "$_c" && -f "$_c/run_regression.sbatch" ]] && { REG="$_c"; break; }
done
[[ -n "$REG" ]] || { echo "cannot locate _regression/; set TESORTER2_DEV" >&2; exit 1; }
DEV="$(dirname "$REG")"
VENV=$REG/venv
mkdir -p "$REG/logs"

# Toolchain (blastp / makeblastdb for BLAST pass-2) and test inputs live outside
# the checkout; the defaults are this machine's, override to run elsewhere.
export PATH="${TESORTER2_TOOLS:-/home/x-dmahechalope/.conda/envs/EDTA-RepGen/bin}:$PATH"
INPUT="${TESORTER2_TEST_LIB:-/anvil/projects/x-bio250374/daniel/bath_vs_hmmer/phase3_data/edta_ltr_100.fa}"

cd "$REG"
rm -rf old_out new_out

echo "=== [1/2] OLD code (pre-refactor, flat src/) ==="
"$VENV/bin/python" "$REG/orig/src/pipeline.py" "$INPUT" \
    -d rexdb -p 8 -o "$REG/old_out" --prefix t

echo "=== [2/2] NEW code (working tree) ==="
# Bundled DBs live under tesorter2/database/ (packaged), not the old top-level
# database/ dir.
#
# The venv holds an installed COPY of the package, so $VENV/bin/tesorter2 would
# test whatever was last pip-installed and report PASS while the working tree
# went unchecked. Run the checkout through PYTHONPATH instead: same interpreter
# and dependencies, code under test is the working tree.
PYTHONPATH="$DEV${PYTHONPATH:+:$PYTHONPATH}" \
"$VENV/bin/python" -m tesorter2.pipeline sequences "$INPUT" \
    -d rexdb -p 8 --db-dir "$DEV/tesorter2/database" -o "$REG/new_out" --prefix t

echo
echo "=== New code must produce everything the old code did ==="
# The package additionally emits the TEsorter-compatible companion files
# (.cls.lib, .cls.pep, .dom.gff3, .dom.tsv, .dom.faa), which the pre-refactor
# src never wrote -- so this is a subset check, not an equality check. Files
# only the OLD code produced are still a failure.
status=0
missing=$(comm -23 <(ls -1 "$REG/old_out") <(ls -1 "$REG/new_out"))
if [[ -n "$missing" ]]; then
    echo "FAIL: outputs produced by old code but not new:"
    echo "$missing"
    status=1
else
    echo "OK: new code produces every old output file"
    added=$(comm -13 <(ls -1 "$REG/old_out") <(ls -1 "$REG/new_out"))
    [[ -n "$added" ]] && { echo "  plus new companion files:"; echo "$added" | sed 's/^/    /'; }
fi

echo
echo "=== DIFF of text outputs actually produced ==="
# .db is SQLite (contains timestamps/rowids) and blast_pass2/ is scratch, so
# compare the deterministic text outputs: classifications and translation.
shopt -s nullglob
for n in "$REG"/new_out/*.tsv "$REG"/new_out/*.aa; do
    f=$(basename "$n"); o="$REG/old_out/$f"
    if [[ ! -f "$o" ]]; then
        # New-only companion output; there is nothing to compare it against.
        echo "NEW ONLY   $f  ($(wc -l < "$n") lines)"; continue
    fi
    # The pre-refactor baseline predates the Sequence-Ontology columns
    # (SO_name/SO_ID, cols 9-10), so compare .cls.tsv on the deterministic
    # original columns 1-8; the .aa translation is compared in full.
    case "$f" in
        *.tsv) ocmp() { cut -f1-8 "$1"; };;
        *)     ocmp() { cat "$1"; };;
    esac
    if diff -q <(ocmp "$o") <(ocmp "$n") >/dev/null; then
        echo "IDENTICAL  $f  ($(wc -l < "$n") lines)"
    else
        echo "DIFFERS    $f"
        diff <(ocmp "$o") <(ocmp "$n") | head -15
        status=1
    fi
done
echo
if [[ $status -eq 0 ]]; then
    echo "RESULT: PASS - refactor is behavior-preserving"
else
    echo "RESULT: FAIL - outputs diverge, investigate above"
fi
exit $status
