#!/bin/bash
# ============================================================
# 06_panukb_gwas_on_gwas.sh — GWAS-on-GWAS enrichment (Experiment 2)
# ============================================================
# Steps:
#   2.1: Collect merged_with_line_numbers.bed from Experiment 1 results
#   2.2: Create universe from BED directory (--feature_files_dir)
#   2.3: Run LSEA on each phenotype against the UKB phenotype universe
# Prerequisites: completed Experiment 1 (05_panukb_enrichment.sh)
# Original: panukb_lsea/2.1-2.3_*.sh
# ============================================================

set -euo pipefail

LSEA_DIR="$(cd "$(dirname "$0")/.." && pwd)"
ENV_FILE="${LSEA_DIR}/.env"
if [ ! -f "${ENV_FILE}" ]; then echo "ERROR: ${ENV_FILE} not found." >&2; exit 1; fi
source "${ENV_FILE}"
[ -f "${CONDA_ACTIVATE:-}" ] && source "${CONDA_ACTIVATE}"

PVAL=0.000000007479176476853146  # 0.05/6685228 (Bonferroni)

VARIANTS=${PANUKB_LSEA_DIR}/in_data/variants.tsv
UKB_BED_DIR=${PANUKB_LSEA_DIR}/ukb_universe_bed
UKB_UNIVERSE_DIR=${PANUKB_LSEA_DIR}/ukb_universe
UKB_UNIVERSE=${UKB_UNIVERSE_DIR}/uni_ukb.json

# ============================================================
# Step 2.1: Collect BED files from Experiment 1 results
# ============================================================
echo "Step 2.1: Collecting BED files from Experiment 1 results..."

mkdir -p "${UKB_BED_DIR}"

for dir in "${PANUKB_LSEA_DIR}/lsea_results"/*/; do
    dir="${dir%/}"
    dirname=$(basename "$dir")

    if [[ "$dirname" == "$(basename "$UKB_BED_DIR")" ]]; then
        continue
    fi

    if [[ -f "$dir/merged_with_line_numbers.bed" ]]; then
        # Strip known category suffix (handles multi-word: go_bp, go_cc, go_mf)
        # "pheno_go_bp" -> "pheno.bed", "pheno_c2" -> "pheno.bed"
        new_name="${dirname}"
        for cat_suffix in "_go_bp" "_go_cc" "_go_mf" "_c2" "_gte" "_bcm"; do
            if [[ "${new_name}" == *"${cat_suffix}" ]]; then
                new_name="${new_name%${cat_suffix}}"
                break
            fi
        done
        cp "$dir/merged_with_line_numbers.bed" "${UKB_BED_DIR}/${new_name}.bed"
    else
        echo "File not found in directory: ${dirname}"
    fi
done

echo "Collected $(ls "${UKB_BED_DIR}"/*.bed 2>/dev/null | wc -l) BED files."

# ============================================================
# Step 2.2: Create universe from BED directory
# ============================================================
echo "Step 2.2: Creating UKB phenotype universe..."

mkdir -p "${UKB_UNIVERSE_DIR}"

python3 "${LSEA_DIR}/universe_generator.py" \
    --variants "${VARIANTS}" \
    --variants_colnames chr pos rsid \
    --feature_files_dir "${UKB_BED_DIR}" \
    --interval 500000 \
    --out_json "${UKB_UNIVERSE}"

echo "UKB universe created: ${UKB_UNIVERSE}"

# ============================================================
# Step 2.3: Run LSEA with UKB phenotype universe
# ============================================================
echo "Step 2.3: Running GWAS-on-GWAS enrichment..."

for file in "${PANUKB_DATA_DIR}"/ukb_summstats/*.tsv.tsv; do
    [ -f "${file}" ] || continue

    filename=$(basename "$file")
    TEMPLATE="${filename%.tsv.tsv}"
    GWAS_NORM="${PANUKB_DATA_DIR}/ukb_summstats/${TEMPLATE}.norm.tsv"

    if [ ! -f "${GWAS_NORM}" ]; then
        echo "SKIP (no normalized file): ${TEMPLATE}"
        continue
    fi

    CUR_OUT_DIR="${PANUKB_LSEA_DIR}/lsea_results/${TEMPLATE}_ukb"
    mkdir -p "${CUR_OUT_DIR}"

    echo "Running: ${TEMPLATE} x ukb"

    python3 "${LSEA_DIR}/LSEA_2.4.py" \
        --input "${GWAS_NORM}" \
        --universe "${UKB_UNIVERSE}" \
        --out "${CUR_OUT_DIR}" \
        --plink_dir "${PLINK_DIR}" \
        --bfile "${PANUKB_BFILE}" \
        --column_names chr pos rsid pval \
        --clump_p1 "${PVAL}" \
        --print_all

    echo "FINISHED FOR ${CUR_OUT_DIR} & ${UKB_UNIVERSE}!"
done

echo "GWAS-on-GWAS enrichment complete."
