#!/bin/bash
#SBATCH --job-name=t2_nhmmer
#SBATCH --account=bio250374
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/nhmmer_%j.out
#SBATCH --error=logs/nhmmer_%j.err
set -euo pipefail
module load python/3.9.5

# See run_regression.sbatch for why the checkout is resolved this way.
REG=""
for _c in "${TESORTER2_DEV:+$TESORTER2_DEV/_regression}" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)" \
          "${SLURM_SUBMIT_DIR:-}"; do
    [[ -n "$_c" && -f "$_c/run_regression.sbatch" ]] && { REG="$_c"; break; }
done
[[ -n "$REG" ]] || { echo "cannot locate _regression/; set TESORTER2_DEV" >&2; exit 1; }
DEV="$(dirname "$REG")"
mkdir -p "$REG/logs"

export PATH="${TESORTER2_TOOLS:-/home/x-dmahechalope/.conda/envs/EDTA-RepGen/bin}:$PATH"
IN="$REG/rice_win.fa"
LIB="${TESORTER2_TEST_LIB:-/anvil/projects/x-bio250374/daniel/bath_vs_hmmer/phase3_data/edta_ltr_100.fa}"
cd "$REG"; rm -rf sine_nh sine_hs

echo "########## SINE db via nhmmer (new default) ##########"
"$REG/venv/bin/python" -m tesorter2.pipeline sequences "$IN" -d sine -p 8 \
    --db-dir "$DEV/tesorter2/database" -o sine_nh --prefix s 2>&1 | grep -E "nhmmer search|Legacy search|hits in|Classified|Loaded"

echo "########## SINE db via hmmsearch (escape hatch) ##########"
"$REG/venv/bin/python" -m tesorter2.pipeline sequences "$IN" -d sine -p 8 --dna-engine hmmsearch \
    --db-dir "$DEV/tesorter2/database" -o sine_hs --prefix s 2>&1 | grep -E "nhmmer search|Legacy search|hits in|Classified|Loaded"

echo
echo "########## Did the AA path stay untouched? (rexdb) ##########"
rm -rf rexdb_check
"$REG/venv/bin/python" -m tesorter2.pipeline sequences "$LIB" -d rexdb -p 8 \
    --db-dir "$DEV/tesorter2/database" -o rexdb_check --prefix t >/dev/null 2>&1
# The baseline predates the SO columns, so compare the original 8 only.
for f in t.cls.tsv t.rexdb.cls.tsv; do
    if diff -q "old_out/$f" <(cut -f1-8 "rexdb_check/$f") >/dev/null; then
        echo "IDENTICAL to pre-refactor baseline: $f"
    else
        echo "REGRESSION in $f"; diff "old_out/$f" "rexdb_check/$f" | head -8
    fi
done

echo
echo "########## strand recorded in the DB? ##########"
"$REG/venv/bin/python" - <<'PY'
import sqlite3
for tag, db in (("nhmmer","sine_nh/s.db"), ("hmmsearch","sine_hs/s.db")):
    c=sqlite3.connect(db)
    rows=c.execute("SELECT strand, COUNT(*) FROM legacy_hits GROUP BY strand ORDER BY 2 DESC").fetchall()
    tot=c.execute("SELECT COUNT(*) FROM legacy_hits").fetchone()[0]
    print(f"  {tag:<10} total={tot:<6} strands={rows}")
PY
