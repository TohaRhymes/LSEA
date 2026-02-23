#!/bin/bash
# ============================================================
# 06_panukb_gwas_on_gwas.sh — GWAS-on-GWAS enrichment (Experiment 2)
# ============================================================
# Uses significant loci from other Pan-UKB phenotypes as "gene sets"
# (phenotype-as-gene-set paradigm) for cross-trait enrichment analysis.
#
# Steps:
#   2.1: Collect merged_with_line_numbers.bed from Experiment 1 results
#        into a directory of BED files (one per phenotype)
#   2.2: Create universe from this BED directory using --feature_files_dir
#   2.3: Run LSEA on each phenotype against the UKB phenotype universe
#
# Requires: completed Experiment 1 (05_panukb_enrichment.sh)
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
OUT_DIR=/media/DATA/gwasim/round2/panukb_lsea

BFILE=${OLD_DATA_DIR}/data/merged_1000genomes_eur
PVAL=0.00000000172487933
INTERVAL=500000

VARIANTS=${OUT_DIR}/in_data/variants.tsv
UKB_BED_DIR=${OUT_DIR}/ukb_universe
UKB_UNIVERSE=${UKB_BED_DIR}/uni_ukb.json

# ============================================================
# Step 2.1: Collect BED files from Experiment 1 results
# ============================================================
echo "Step 2.1: Collecting BED files from Experiment 1 results..."

mkdir -p "${UKB_BED_DIR}"

# For each phenotype result directory, copy merged_with_line_numbers.bed
# as a named BED file (phenotype name becomes the "gene set" name)
for result_dir in "${OUT_DIR}/lsea_results"/*; do
    [ -d "${result_dir}" ] || continue
    bed_file="${result_dir}/merged_with_line_numbers.bed"
    if [ -f "${bed_file}" ]; then
        pheno_name=$(basename "${result_dir}")
        cp "${bed_file}" "${UKB_BED_DIR}/${pheno_name}.bed"
    fi
done

echo "Collected $(ls "${UKB_BED_DIR}"/*.bed 2>/dev/null | wc -l) BED files."

# ============================================================
# Step 2.2: Create universe from BED directory
# ============================================================
echo "Step 2.2: Creating UKB phenotype universe..."

python3 "${LSEA_DIR}/universe_generator.py" \
    --variants "${VARIANTS}" \
    --feature_files_dir "${UKB_BED_DIR}" \
    --interval "${INTERVAL}" \
    --out_json "${UKB_UNIVERSE}"

echo "UKB universe created: ${UKB_UNIVERSE}"

# ============================================================
# Step 2.3: Run LSEA with UKB phenotype universe
# ============================================================
echo "Step 2.3: Running GWAS-on-GWAS enrichment..."

GWAS_DIR=${PAN_UKB_DIR}

for gwas_file in "${GWAS_DIR}"/*.tsv; do
    [ -f "${gwas_file}" ] || continue

    pheno=$(basename "${gwas_file}" .tsv)
    RESULT_DIR=${OUT_DIR}/lsea_results/${pheno}_ukb

    echo "Running: ${pheno} x ukb"

    python3 "${LSEA_DIR}/LSEA_2.4.py" \
        --input "${gwas_file}" \
        --universe "${UKB_UNIVERSE}" \
        --out "${RESULT_DIR}" \
        --plink_dir "${PLINK_DIR}" \
        --bfile "${BFILE}" \
        --column_names chr pos rsid pval \
        --clump_p1 "${PVAL}" \
        --print_all

    echo "Done: ${pheno} x ukb"
done

echo "GWAS-on-GWAS enrichment complete."
