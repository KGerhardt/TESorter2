#!/bin/bash
#SBATCH --job-name=t2_so
#SBATCH --account=bio250374
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:40:00
#SBATCH --output=logs/so_%j.out
#SBATCH --error=logs/so_%j.err
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
IN="${TESORTER2_TEST_LIB:-/anvil/projects/x-bio250374/daniel/bath_vs_hmmer/phase3_data/edta_ltr_100.fa}"
GENOME="${TESORTER2_TEST_GENOME:-/anvil/projects/x-bio250374/daniel/minirice.fa}"
cd "$REG"; rm -rf so_elem so_compat so_genome

echo "=== 1. element mode: SO columns in cls.tsv ==="
"$REG/venv/bin/tesorter2" sequences "$IN" -d rexdb -p 8 -o so_elem --prefix e >/dev/null 2>&1
head -1 so_elem/e.cls.tsv
echo "--- filas de ejemplo (TE, Order, Superfamily, SO_name, SO_ID) ---"
awk -F'\t' 'NR>1{print "  "$2"/"$3" -> "$(NF-1)"  "$NF}' so_elem/e.cls.tsv | sort -u | head -5
echo "--- ¿alguna fila sin SO? ---"
awk -F'\t' 'NR>1 && ($NF=="" || $NF !~ /^SO:[0-9]{7}$/){print "  SIN SO: "$0}' so_elem/e.cls.tsv | head -3
echo "  (vacío = todas tienen SO válido)"

echo
echo "=== 2. compat mode: DEBE seguir siendo 7 columnas exactas ==="
"$REG/venv/bin/tesorter2" sequences "$IN" -d rexdb -p 8 --compat-tesorter-output -o so_compat --prefix c >/dev/null 2>&1
head -1 so_compat/c.cls.tsv
echo "  columnas: $(head -1 so_compat/c.cls.tsv | awk -F'\t' '{print NF}')  (debe ser 7)"

echo
echo "=== 3. genome mode: Ontology_term en el GFF3 ==="
"$REG/venv/bin/tesorter2" genome "$GENOME" -d rexdb -p 8 \
    -o so_genome --prefix g >/dev/null 2>&1
G=$(ls so_genome/*.gff3 2>/dev/null | head -1)
echo "  gff3: $G"
# awk, not `| head -2`: under `set -eo pipefail` head closing the pipe kills
# grep with SIGPIPE and aborts the script before the check below ever runs.
awk '!/^#/ && ++n<=2 {print substr($0, 1, 190)}' "$G"
echo "  features con Ontology_term: $(grep -c "Ontology_term=SO:" "$G" 2>/dev/null || echo 0) / $(grep -vc '^#' "$G" 2>/dev/null || echo 0)"
