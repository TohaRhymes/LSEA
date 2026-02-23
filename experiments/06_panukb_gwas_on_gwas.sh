#!/bin/bash
# ============================================================
# 06_panukb_gwas_on_gwas.sh — GWAS-on-GWAS enrichment (Experiment 2)
# ============================================================
# Uses significant loci from other Pan-UKB phenotypes as "gene sets"
# (phenotype-as-gene-set paradigm) for cross-trait enrichment analysis.
#
# Steps:
#   2.1: Collect merged_with_line_numbers.bed from Experiment 1 results
#        into a directory of BED files (one per phenotype).
#        The category suffix (_c2, _gte, etc.) is stripped from dir name.
#   2.2: Create universe from this BED directory using --feature_files_dir
#   2.3: Run LSEA on each phenotype against the UKB phenotype universe
#
# Prerequisites:
#   - Completed Experiment 1 (05_panukb_enrichment.sh)
#   - Pre-normalized GWAS files (*.norm.tsv) in ukb_summstats/
#
# Original: panukb_lsea/2.1_make_universe_bed_dir.sh,
#           panukb_lsea/2.2_prepare_lsea_uni.sh,
#           panukb_lsea/2.3_iter_lsea.sh
# Updated:  CLI flags changed to --long_flag format
# ============================================================

set -euo pipefail

# --- Paths (adjust for your server) ---
LSEA_DIR="$(cd "$(dirname "$0")/.." && pwd)"
PLINK_DIR=/home/achangalidi/tools/plink
OLD_DATA_DIR=/media/DATA/gwasim/round2/bioGWAS/tests
PAN_UKB_DIR=/media/DATA/gwasim/round2/panukb

BFILE=${OLD_DATA_DIR}/data/merged_1000genomes_eur
PVAL=0.000000007479176476853146  # 0.05/6685228 (Bonferroni)

OUT_DIR=/media/DATA/gwasim/round2/panukb_lsea
VARIANTS=${OUT_DIR}/in_data/variants.tsv
UKB_BED_DIR=${OUT_DIR}/ukb_universe_bed
UKB_UNIVERSE_DIR=${OUT_DIR}/ukb_universe
UKB_UNIVERSE=${UKB_UNIVERSE_DIR}/uni_ukb.json

# ============================================================
# Step 2.1: Collect BED files from Experiment 1 results
# ============================================================
echo "Step 2.1: Collecting BED files from Experiment 1 results..."

mkdir -p "${UKB_BED_DIR}"

for dir in "${OUT_DIR}/lsea_results"/*/; do
    dir="${dir%/}"
    dirname=$(basename "$dir")

    # Skip the output directory itself
    if [[ "$dirname" == "$(basename "$UKB_BED_DIR")" ]]; then
        continue
    fi

    if [[ -f "$dir/merged_with_line_numbers.bed" ]]; then
        # Strip category suffix: "pheno_c2" -> "pheno.bed"
        new_name="${dirname%_*}.bed"
        cp "$dir/merged_with_line_numbers.bed" "${UKB_BED_DIR}/${new_name}"
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
    --feature_files_dir "${UKB_BED_DIR}" \
    --interval 500000 \
    --out_json "${UKB_UNIVERSE}"

echo "UKB universe created: ${UKB_UNIVERSE}"

# ============================================================
# Step 2.3: Run LSEA with UKB phenotype universe
# ============================================================
echo "Step 2.3: Running GWAS-on-GWAS enrichment..."

for file in "${PAN_UKB_DIR}"/ukb_summstats/*.tsv.tsv; do
    [ -f "${file}" ] || continue

    filename=$(basename "$file")
    TEMPLATE="${filename%.tsv.tsv}"

    GWAS_NORM="${PAN_UKB_DIR}/ukb_summstats/${TEMPLATE}.norm.tsv"

    if [ ! -f "${GWAS_NORM}" ]; then
        echo "SKIP (no normalized file): ${TEMPLATE}"
        continue
    fi

    CUR_OUT_DIR="${OUT_DIR}/lsea_results/${TEMPLATE}_ukb"
    mkdir -p "${CUR_OUT_DIR}"

    echo "Running: ${TEMPLATE} x ukb"

    python3 "${LSEA_DIR}/LSEA_2.4.py" \
        --input "${GWAS_NORM}" \
        --universe "${UKB_UNIVERSE}" \
        --out "${CUR_OUT_DIR}" \
        --plink_dir "${PLINK_DIR}" \
        --bfile "${BFILE}" \
        --column_names chr pos rsid pval \
        --clump_p1 "${PVAL}" \
        --print_all

    echo "FINISHED FOR ${CUR_OUT_DIR} & ${UKB_UNIVERSE}!"
done

echo "GWAS-on-GWAS enrichment complete."
