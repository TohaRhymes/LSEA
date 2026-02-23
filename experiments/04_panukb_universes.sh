#!/bin/bash
# ============================================================
# 04_panukb_universes.sh — Universe creation for Pan-UKB experiments
# ============================================================
# Creates 6 universe JSON files, one per gene set category.
# Original: panukb_lsea/1.1_prepare_lsea_uni.sh
# ============================================================

set -euo pipefail

LSEA_DIR="$(cd "$(dirname "$0")/.." && pwd)"
ENV_FILE="${LSEA_DIR}/.env"
if [ ! -f "${ENV_FILE}" ]; then echo "ERROR: ${ENV_FILE} not found." >&2; exit 1; fi
source "${ENV_FILE}"
[ -f "${CONDA_ACTIVATE:-}" ] && source "${CONDA_ACTIVATE}"

# Source GWAS file (any Pan-UKB file — we only need variant positions)
SOME_GWAS=${PANUKB_DATA_DIR}/biomarkers-30600-both_sexes-irnt.tsv.tsv

# Output paths
VARIANTS=${PANUKB_LSEA_DIR}/in_data/variants.tsv
BED=${PANUKB_LSEA_DIR}/in_data/anno.bed

INTERVAL=500000

# --- GMT files for each gene set category ---
declare -A GMT_FILES
GMT_FILES[c2]=${GMT_C2}
GMT_FILES[gte]=${GMT_GTE}
GMT_FILES[bcm]=${GMT_BCM}
GMT_FILES[go_bp]=${GMT_GO_BP}
GMT_FILES[go_cc]=${GMT_GO_CC}
GMT_FILES[go_mf]=${GMT_GO_MF}

# --- Prepare input data ---
mkdir -p "${PANUKB_LSEA_DIR}/in_data"

# Extract variant positions from Pan-UKB GWAS summary stats
awk 'BEGIN {OFS="\t"}
    NR==1 {print "chr", "pos", "rsid"}
    NR > 1 {
        rsid = $1 ":" $2 ":" $3 ":" $4
        print $1, $2, rsid
    }' "${SOME_GWAS}" > "${VARIANTS}"

# Extract gene annotations from GENCODE GTF -> BED format
awk -F'\t' 'NR>=6 && $3=="gene"' "${GENCODE_GTF}" | \
    awk '{sub(/^chr/, "", $1); match($0, /gene_name "([^"]+)"/, arr); geneName=arr[1]; print $1"\t"$4"\t"$5"\t"geneName}' \
    > "${BED}"

# --- Create universe for each gene set category ---
for category in c2 gte bcm go_bp go_cc go_mf; do
    GMT=${GMT_FILES[$category]}
    UNIVERSE=${PANUKB_LSEA_DIR}/in_data/uni_${category}.json

    echo "Creating universe for ${category} (GMT: ${GMT})..."

    python3 "${LSEA_DIR}/universe_generator.py" \
        --variants "${VARIANTS}" \
        --features "${BED}" "${GMT}" \
        --interval "${INTERVAL}" \
        --out_json "${UNIVERSE}"

    echo "Universe created: ${UNIVERSE}"
done

echo "All 6 universes created."
