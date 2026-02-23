#!/bin/bash
# ============================================================
# 05_panukb_enrichment.sh — Pan-UKB enrichment analysis (Experiment 1)
# ============================================================
# Runs LSEA on Pan-UKB phenotypes across 6 gene set categories.
# Each phenotype is tested against each universe independently.
#
# Prerequisites:
#   - Universes from 04_panukb_universes.sh
#   - Pre-normalized GWAS files (*.norm.tsv) in ukb_summstats/
#     Created by preproc_tsv.py (converts neglog10_pval_EUR -> pval,
#     creates rsid=chr:pos:ref:alt, filters to bfile variants)
#
# Original: panukb_lsea/1.2_iter_lsea.sh
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

# Universes for each gene set category
UNIVERSE_C2=${OUT_DIR}/in_data/uni_c2.json
UNIVERSE_GTE=${OUT_DIR}/in_data/uni_gte.json
UNIVERSE_BCM=${OUT_DIR}/in_data/uni_bcm.json
UNIVERSE_GO_BP=${OUT_DIR}/in_data/uni_go_bp.json
UNIVERSE_GO_CC=${OUT_DIR}/in_data/uni_go_cc.json
UNIVERSE_GO_MF=${OUT_DIR}/in_data/uni_go_mf.json

NAMES=("c2" "gte" "bcm" "go_cc" "go_mf" "go_bp")
UNIVERSES=("$UNIVERSE_C2" "$UNIVERSE_GTE" "$UNIVERSE_BCM" "$UNIVERSE_GO_CC" "$UNIVERSE_GO_MF" "$UNIVERSE_GO_BP")

# --- Iterate over all pre-normalized GWAS files ---
for file in "${PAN_UKB_DIR}"/ukb_summstats/*.tsv.tsv; do
    [ -f "${file}" ] || continue

    filename=$(basename "$file")
    TEMPLATE="${filename%.tsv.tsv}"

    GWAS_NORM="${PAN_UKB_DIR}/ukb_summstats/${TEMPLATE}.norm.tsv"

    if [ ! -f "${GWAS_NORM}" ]; then
        echo "SKIP (no normalized file): ${TEMPLATE}"
        continue
    fi

    for i in "${!UNIVERSES[@]}"; do
        universe="${UNIVERSES[i]}"
        name="${NAMES[i]}"

        CUR_OUT_DIR="${OUT_DIR}/lsea_results/${TEMPLATE}_${name}"
        mkdir -p "${CUR_OUT_DIR}"

        echo "Running: ${TEMPLATE} x ${name}"

        python3 "${LSEA_DIR}/LSEA_2.4.py" \
            --input "${GWAS_NORM}" \
            --universe "${universe}" \
            --out "${CUR_OUT_DIR}" \
            --plink_dir "${PLINK_DIR}" \
            --bfile "${BFILE}" \
            --column_names chr pos rsid pval \
            --clump_p1 "${PVAL}"

        echo "FINISHED FOR ${CUR_OUT_DIR} & ${universe}!"
    done
done

echo "Pan-UKB enrichment analysis complete."
