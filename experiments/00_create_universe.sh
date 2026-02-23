#!/bin/bash
# ============================================================
# 00_create_universe.sh — Universe generation for bioGWAS experiments
# ============================================================
# Creates a universe JSON from variant positions + gene annotations + GMT.
# Run ONCE before any LSEA analysis (scripts 01-03 depend on this).
#
# Original: lsea_test/create_uni.sh
# Updated:  CLI flags changed to --long_flag format
# ============================================================

set -euo pipefail

# --- Paths (adjust for your server) ---
LSEA_DIR="$(cd "$(dirname "$0")/.." && pwd)"
OLD_DATA_DIR=/media/DATA/gwasim/round2/bioGWAS/tests
OUT_DIR=/media/DATA/gwasim/round2/lsea_test

# Source GWAS file (any one — we only need variant positions)
TEMPLATE=test10000_path_small_4_path_small_4
SOME_GWAS=${OLD_DATA_DIR}/3_pathways/in_data/${TEMPLATE}_gwas.tsv

# Reference gene annotations
GTF=${OLD_DATA_DIR}/data/gencode.v37.annotation.gtf
GMT_SRC=${OLD_DATA_DIR}/data/c2.cp.kegg.v2023.1.Hs.symbols.gmt

# Output files
VARIANTS=${OUT_DIR}/in_data/variants.tsv
GMT=${OUT_DIR}/in_data/c2.cp.kegg.v2023.1.Hs.symbols.gmt
BED=${OUT_DIR}/in_data/anno.bed
UNIVERSE=${OUT_DIR}/in_data/uni.json

INTERVAL=500000

# --- Prepare input data ---
mkdir -p "${OUT_DIR}/in_data"

# Extract variant positions (chr, pos, rsid) from a GWAS summary stats file
awk '{print $1"\t"$3"\t"$2}' "${SOME_GWAS}" > "${VARIANTS}"

# Copy GMT file
cp "${GMT_SRC}" "${GMT}"

# Extract gene annotations from GENCODE GTF -> BED format
# Strips 'chr' prefix from chromosome names for consistency
awk -F'\t' 'NR>=6 && $3=="gene"' "${GTF}" | \
    awk '{sub(/^chr/, "", $1); match($0, /gene_name "([^"]+)"/, arr); geneName=arr[1]; print $1"\t"$4"\t"$5"\t"geneName}' \
    > "${BED}"

# --- Run universe generator ---
python3 "${LSEA_DIR}/universe_generator.py" \
    --variants "${VARIANTS}" \
    --features "${BED}" "${GMT}" \
    --interval "${INTERVAL}" \
    --out_json "${UNIVERSE}"

echo "Universe created: ${UNIVERSE}"
