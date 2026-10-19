#!/bin/bash
#SBATCH --job-name=t2_gsine
#SBATCH --account=bio250374
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:40:00
#SBATCH --output=logs/gsine_%j.out
#SBATCH --error=logs/gsine_%j.err

# End-to-end check for genome mode + SINE: genome mode routes DNA databases
# (e.g. sine) through nhmmer to find non-coding elements, alongside the protein
# domains found by BATH (the only protein engine in genome mode). Both engines
# report nucleotide coordinates, so every run emits .dom.fna and never .dom.faa.
# Requires bathsearch/bathconvert for the rexdb runs. Validates:
#   1. mixed rexdb+sine: SINE_element features on BOTH strands, protein domains
#      still typed CDS, one .fna holding every feature and no .faa.
#   2. sine-only: nucleotide-only output, no BATH required.
#   3. rexdb-only: all CDS, .fna only.
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
PY="$REG/venv/bin/python"
DBDIR="$DEV/tesorter2/database"
GENOME="${TESORTER2_TEST_GENOME:-/anvil/projects/x-bio250374/daniel/minirice.fa}"

cd "$REG"; rm -rf gsine_mixed gsine_sine gsine_rexdb
status=0

pass() { echo "  PASS: $1"; }
fail() { echo "  FAIL: $1"; status=1; }

echo "########## 1. mixed genome: rexdb + sine ##########"
"$PY" -m tesorter2.pipeline genome "$GENOME" -d rexdb,sine -p 8 \
    --db-dir "$DBDIR" -o gsine_mixed --prefix g >/dev/null 2>&1
G=gsine_mixed/g.dom.gff3

n_sine=$(grep -c "Classification=SINE" "$G" || true)
n_cds=$(awk -F'\t' '$3=="CDS"' "$G" | wc -l)
n_sine_type=$(awk -F'\t' '$3=="SINE_element"' "$G" | wc -l)
[ "$n_sine" -gt 0 ] && pass "SINE features found ($n_sine)" || fail "no SINE features"
[ "$n_sine_type" -eq "$n_sine" ] && pass "all SINEs typed SINE_element ($n_sine_type)" \
    || fail "SINE type mismatch: $n_sine_type SINE_element vs $n_sine Classification=SINE"
[ "$n_cds" -gt 0 ] && pass "protein domains still typed CDS ($n_cds)" || fail "no CDS features"

# SINE_element rows must carry the SINE SO term, never a CDS-style one.
bad_so=$(awk -F'\t' '$3=="SINE_element" && $9 !~ /Ontology_term=SO:0000206/' "$G" | wc -l)
[ "$bad_so" -eq 0 ] && pass "every SINE_element carries Ontology_term=SO:0000206" \
    || fail "$bad_so SINE_element features missing/incorrect SO term"

strands=$(awk -F'\t' '$3=="SINE_element"{print $7}' "$G" | sort -u | tr -d '\n')
[ "$strands" = "+-" ] && pass "SINEs recovered on BOTH strands" \
    || fail "SINE strands not both +/-: got '$strands'"

# Protein search is BATH-only, so BATH domains and nhmmer elements are both
# nucleotide and share one file. .dom.faa is never written by genome mode.
if [ -s gsine_mixed/g.dom.fna ] && [ ! -e gsine_mixed/g.dom.faa ]; then
    pass "mixed run emits g.dom.fna only (no .faa)"
else
    fail "expected g.dom.fna non-empty and no g.dom.faa"
fi
# Every record must be nucleotide (ACGTN only), BATH domains included.
if grep -v "^>" gsine_mixed/g.dom.fna | grep -qiE "[^ACGTN]"; then
    fail "g.dom.fna contains non-nucleotide characters"
else
    pass "g.dom.fna is nucleotide"
fi
# feature count in GFF must equal records in the single FASTA.
n_feat=$(grep -vc "^#" "$G")
n_seq=$(grep -c "^>" gsine_mixed/g.dom.fna)
[ "$n_feat" -eq "$n_seq" ] && pass "GFF features == FASTA records ($n_feat)" \
    || fail "GFF $n_feat features vs $n_seq FASTA records"
# summary table has a SINE order row.
grep -qE "^SINE\b" gsine_mixed/g.genome.summary.tsv && pass "summary has SINE row" \
    || fail "summary missing SINE row"

echo
echo "########## 2. sine-only genome: nucleotide single-kind ##########"
"$PY" -m tesorter2.pipeline genome "$GENOME" -d sine -p 8 \
    --db-dir "$DBDIR" -o gsine_sine --prefix s >/dev/null 2>&1
if [ -s gsine_sine/s.dom.fna ] && [ ! -e gsine_sine/s.dom.faa ]; then
    pass "sine-only emits .fna only (no .faa)"
else
    fail "sine-only output wrong: expected .fna only"
fi
only_sine=$(awk -F'\t' '$3!="" && $3!="SINE_element" && $0 !~ /^#/' gsine_sine/s.dom.gff3 | wc -l)
[ "$only_sine" -eq 0 ] && pass "sine-only GFF is all SINE_element" \
    || fail "sine-only GFF has $only_sine non-SINE_element features"

echo
echo "########## 3. rexdb-only genome: AA single-kind (no regression) ##########"
"$PY" -m tesorter2.pipeline genome "$GENOME" -d rexdb -p 8 \
    --db-dir "$DBDIR" -o gsine_rexdb --prefix r >/dev/null 2>&1
# BATH reports nucleotide coordinates on the whole chromosome, so even a
# protein-only genome run emits nucleotide sub-sequences.
if [ -s gsine_rexdb/r.dom.fna ] && [ ! -e gsine_rexdb/r.dom.faa ]; then
    pass "rexdb-only emits .fna only (no .faa)"
else
    fail "rexdb-only output wrong: expected .fna only"
fi
# Neither engine needs the genome windowed, so no run writes cut.fa.
[ ! -e gsine_rexdb/cut.fa ] && [ ! -e gsine_sine/cut.fa ] \
    && pass "no run windows the genome (no cut.fa)" \
    || fail "a genome run wrote cut.fa"
only_cds=$(awk -F'\t' '$3!="" && $3!="CDS" && $0 !~ /^#/' gsine_rexdb/r.dom.gff3 | wc -l)
[ "$only_cds" -eq 0 ] && pass "rexdb-only GFF is all CDS" \
    || fail "rexdb-only GFF has $only_cds non-CDS features"

echo
echo "########## 4. reconciliation: SINE annotation independent of co-searched DBs ##########"
# Overlaps resolve within-type, so
# the SINE_element set must be byte-identical whether or not rexdb was co-searched
# (no cross-type eviction).
awk -F'\t' '$3=="SINE_element"{print $1"\t"$4"\t"$5"\t"$7}' gsine_sine/s.dom.gff3 | sort > /tmp/gsine_only.$$
awk -F'\t' '$3=="SINE_element"{print $1"\t"$4"\t"$5"\t"$7}' gsine_mixed/g.dom.gff3 | sort > /tmp/gsine_mix.$$
if diff -q /tmp/gsine_only.$$ /tmp/gsine_mix.$$ >/dev/null; then
    pass "SINE set identical between 'sine' and 'rexdb,sine' ($(wc -l < /tmp/gsine_only.$$) SINEs)"
else
    fail "SINE set differs between sine-only and mixed:"; diff /tmp/gsine_only.$$ /tmp/gsine_mix.$$ | head -6
fi
rm -f /tmp/gsine_only.$$ /tmp/gsine_mix.$$

echo
if [ "$status" -eq 0 ]; then
    echo "RESULT: PASS - genome mode + SINE end-to-end"
else
    echo "RESULT: FAIL - see failures above"
fi
exit $status
