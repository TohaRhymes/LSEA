#!/bin/bash
# ============================================================
# 05_panukb_enrichment.sh — Pan-UKB enrichment analysis (Experiment 1)
# ============================================================
# Runs LSEA on Pan-UKB phenotypes across 6 gene set categories.
# Prerequisites:
#   - Universes from 04_panukb_universes.sh
#   - Pre-normalized GWAS files (*.norm.tsv) in ukb_summstats/
# Original: panukb_lsea/1.2_iter_lsea.sh
# ============================================================

set -euo pipefail

LSEA_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
ENV_FILE="${LSEA_DIR}/.env"
if [ ! -f "${ENV_FILE}" ]; then echo "ERROR: ${ENV_FILE} not found." >&2; exit 1; fi
source "${ENV_FILE}"
[ -f "${CONDA_ACTIVATE:-}" ] && source "${CONDA_ACTIVATE}"

PVAL=0.000000007479176476853146  # 0.05/6685228 (Bonferroni)

NAMES=("c2" "gte" "bcm" "go_cc" "go_mf" "go_bp")
UNIVERSES=(
    "${PANUKB_LSEA_DIR}/in_data/uni_c2.json"
    "${PANUKB_LSEA_DIR}/in_data/uni_gte.json"
    "${PANUKB_LSEA_DIR}/in_data/uni_bcm.json"
    "${PANUKB_LSEA_DIR}/in_data/uni_go_cc.json"
    "${PANUKB_LSEA_DIR}/in_data/uni_go_mf.json"
    "${PANUKB_LSEA_DIR}/in_data/uni_go_bp.json"
)

for file in "${PANUKB_DATA_DIR}"/ukb_summstats/*.tsv.tsv; do
    [ -f "${file}" ] || continue

    filename=$(basename "$file")
    TEMPLATE="${filename%.tsv.tsv}"
    GWAS_NORM="${PANUKB_DATA_DIR}/ukb_summstats/${TEMPLATE}.norm.tsv"

    if [ ! -f "${GWAS_NORM}" ]; then
        echo "SKIP (no normalized file): ${TEMPLATE}"
        continue
    fi

    for i in "${!UNIVERSES[@]}"; do
        universe="${UNIVERSES[i]}"
        name="${NAMES[i]}"

        CUR_OUT_DIR="${PANUKB_LSEA_DIR}/lsea_results/${TEMPLATE}_${name}"
        mkdir -p "${CUR_OUT_DIR}"

        echo "Running: ${TEMPLATE} x ${name}"

        python3 "${LSEA_DIR}/LSEA_2.4.py" \
            --input "${GWAS_NORM}" \
            --universe "${universe}" \
            --out "${CUR_OUT_DIR}" \
            --plink_dir "${PLINK_DIR}" \
            --bfile "${PANUKB_BFILE}" \
            --column_names chr pos rsid pval \
            --clump_p1 "${PVAL}"

        echo "FINISHED FOR ${CUR_OUT_DIR} & ${universe}!"
    done
done

echo "Pan-UKB enrichment analysis complete."
