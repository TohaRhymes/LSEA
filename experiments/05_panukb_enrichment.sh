#!/bin/bash
# ============================================================
# 05_panukb_enrichment.sh — Pan-UKB enrichment analysis (Experiment 1)
# ============================================================
# Runs LSEA on 139 Pan-UKB phenotypes across 6 gene set categories.
# Each phenotype is tested against each universe independently.
#
# Requires: universes from 04_panukb_universes.sh
# Original: panukb_lsea/1.2_iter_lsea.sh
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

# Gene set categories
CATEGORIES="c2 gte bcm go_bp go_cc go_mf"

# --- Find all Pan-UKB GWAS files ---
# Each file is a phenotype summary statistics TSV
GWAS_DIR=${PAN_UKB_DIR}

echo "Starting Pan-UKB enrichment analysis..."

for gwas_file in "${GWAS_DIR}"/*.tsv; do
    [ -f "${gwas_file}" ] || continue

    # Extract phenotype name from filename (e.g., "biomarkers-30600-both_sexes-irnt.tsv")
    pheno=$(basename "${gwas_file}" .tsv)

    for category in ${CATEGORIES}; do
        UNIVERSE=${OUT_DIR}/in_data/uni_${category}.json
        RESULT_DIR=${OUT_DIR}/lsea_results/${pheno}_${category}

        echo "Running: ${pheno} x ${category}"

        python3 "${LSEA_DIR}/LSEA_2.4.py" \
            --input "${gwas_file}" \
            --universe "${UNIVERSE}" \
            --out "${RESULT_DIR}" \
            --plink_dir "${PLINK_DIR}" \
            --bfile "${BFILE}" \
            --column_names chr pos rsid pval \
            --clump_p1 "${PVAL}" \
            --print_all

        echo "Done: ${pheno} x ${category}"
    done
done

echo "Pan-UKB enrichment analysis complete."
