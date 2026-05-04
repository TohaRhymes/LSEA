#!/bin/bash
# ============================================================
# 00_create_universe.sh — Universe generation for bioGWAS experiments
# ============================================================
# Creates a universe JSON from variant positions + gene annotations + GMT.
# Run ONCE before any LSEA analysis (scripts 01-03 depend on this).
#
# Original: lsea_test/create_uni.sh
# ============================================================

set -euo pipefail

LSEA_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
ENV_FILE="${LSEA_DIR}/.env"
if [ ! -f "${ENV_FILE}" ]; then echo "ERROR: ${ENV_FILE} not found. Copy .env.example to .env and fill in paths." >&2; exit 1; fi
source "${ENV_FILE}"
[ -f "${CONDA_ACTIVATE:-}" ] && source "${CONDA_ACTIVATE}"

# --- Derived paths ---
TEMPLATE=test10000_path_small_4_path_small_4
SOME_GWAS=${BIOGWAS_DATA_DIR}/3_pathways/in_data/${TEMPLATE}_gwas.tsv
GTF=${GENCODE_GTF}
GMT_SRC=${GMT_C2_KEGG}

VARIANTS=${LSEA_TEST_DIR}/in_data/variants.tsv
GMT=${LSEA_TEST_DIR}/in_data/c2.cp.kegg.v2023.1.Hs.symbols.gmt
BED=${LSEA_TEST_DIR}/in_data/anno.bed
UNIVERSE=${LSEA_TEST_DIR}/in_data/uni.json

INTERVAL=500000

# --- Prepare input data ---
mkdir -p "${LSEA_TEST_DIR}/in_data"

awk '{print $1"\t"$3"\t"$2}' "${SOME_GWAS}" > "${VARIANTS}"
cp "${GMT_SRC}" "${GMT}"
awk -F'\t' 'NR>=6 && $3=="gene"' "${GTF}" | \
    awk '{sub(/^chr/, "", $1); match($0, /gene_name "([^"]+)"/, arr); geneName=arr[1]; print $1"\t"$4"\t"$5"\t"geneName}' \
    > "${BED}"

# --- Run universe generator ---
python3 "${LSEA_DIR}/universe_generator.py" \
    --variants "${VARIANTS}" \
    --variants_colnames chr pos rsid \
    --features "${BED}" "${GMT}" \
    --interval "${INTERVAL}" \
    --out_json "${UNIVERSE}"

echo "Universe created: ${UNIVERSE}"
