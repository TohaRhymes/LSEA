#!/bin/bash
# ============================================================
# 04_panukb_universes.sh — Universe creation for Pan-UKB experiments
# ============================================================
# Creates 6 universe JSON files, one per gene set category,
# for the Pan-UK Biobank enrichment analysis.
#
# Gene set categories:
#   c2     — MSigDB C2 (KEGG pathways)
#   gte    — GTEx v8 tissue-specific expression
#   bcm    — Blood cell marker genes
#   go_bp  — GO Biological Process
#   go_cc  — GO Cellular Component
#   go_mf  — GO Molecular Function
#
# Original: panukb_lsea/1.1_prepare_lsea_uni.sh
# Updated:  CLI flags changed to --long_flag format
# ============================================================

set -euo pipefail

# --- Paths (adjust for your server) ---
LSEA_DIR="$(cd "$(dirname "$0")/.." && pwd)"
OLD_DATA_DIR=/media/DATA/gwasim/round2/bioGWAS/tests
PAN_UKB_DIR=/media/DATA/gwasim/round2/panukb
OUT_DIR=/media/DATA/gwasim/round2/panukb_lsea

# Source GWAS file (any Pan-UKB file — we only need variant positions)
SOME_GWAS=${PAN_UKB_DIR}/biomarkers-30600-both_sexes-irnt.tsv.tsv

# Reference gene annotations
GTF=${OLD_DATA_DIR}/data/gencode.v37.annotation.gtf

# Output paths
VARIANTS=${OUT_DIR}/in_data/variants.tsv
BED=${OUT_DIR}/in_data/anno.bed

INTERVAL=500000

# --- GMT files for each gene set category ---
declare -A GMT_FILES
GMT_FILES[c2]=${OLD_DATA_DIR}/data/c2.cp.kegg.v2023.1.Hs.symbols.gmt
GMT_FILES[gte]=/media/DATA/bioinformatics/LSEA/tissues/GTEx8_formatted.gmt
GMT_FILES[bcm]=/media/DATA/bioinformatics/LSEA/tissues/blood_cell_markers.gmt
GMT_FILES[go_bp]=${OLD_DATA_DIR}/data/c5.go.bp.v2024.1.Hs.symbols.gmt
GMT_FILES[go_cc]=${OLD_DATA_DIR}/data/c5.go.cc.v2024.1.Hs.symbols.gmt
GMT_FILES[go_mf]=${OLD_DATA_DIR}/data/c5.go.mf.v2024.1.Hs.symbols.gmt

# --- Prepare input data ---
mkdir -p "${OUT_DIR}/in_data"

# Extract variant positions from Pan-UKB GWAS summary stats
# Format: chr:pos:ref:alt as rsid
awk 'BEGIN {OFS="\t"}
    NR==1 {print "chr", "pos", "rsid"}
    NR > 1 {
        rsid = $1 ":" $2 ":" $3 ":" $4
        print $1, $2, rsid
    }' "${SOME_GWAS}" > "${VARIANTS}"

# Extract gene annotations from GENCODE GTF -> BED format
awk -F'\t' 'NR>=6 && $3=="gene"' "${GTF}" | \
    awk '{sub(/^chr/, "", $1); match($0, /gene_name "([^"]+)"/, arr); geneName=arr[1]; print $1"\t"$4"\t"$5"\t"geneName}' \
    > "${BED}"

# --- Create universe for each gene set category ---
for category in c2 gte bcm go_bp go_cc go_mf; do
    GMT=${GMT_FILES[$category]}
    UNIVERSE=${OUT_DIR}/in_data/uni_${category}.json

    echo "Creating universe for ${category} (GMT: ${GMT})..."

    python3 "${LSEA_DIR}/universe_generator.py" \
        --variants "${VARIANTS}" \
        --features "${BED}" "${GMT}" \
        --interval "${INTERVAL}" \
        --out_json "${UNIVERSE}"

    echo "Universe created: ${UNIVERSE}"
done

echo "All 6 universes created."
