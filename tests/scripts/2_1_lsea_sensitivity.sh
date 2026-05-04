#!/bin/bash
# ============================================================
# 03_lsea_sensitivity.sh — LSEA sensitivity analysis (varying k)
# ============================================================
# 135 runs: k=1..15, 3 iterations x 3 pathway sizes (big, medium, small).
# Requires: universe from 00_create_universe.sh
# Original: lsea_test/iter_lsea_k.sh
# ============================================================

set -euo pipefail

LSEA_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
ENV_FILE="${LSEA_DIR}/.env"
if [ ! -f "${ENV_FILE}" ]; then echo "ERROR: ${ENV_FILE} not found." >&2; exit 1; fi
source "${ENV_FILE}"
[ -f "${CONDA_ACTIVATE:-}" ] && source "${CONDA_ACTIVATE}"

UNIVERSE=${LSEA_TEST_DIR}/in_data/uni.json
PVAL=0.00000000729730

for k in $(seq 1 15); do
    for i in $(seq 0 2); do
        for path_size in big medium small; do
            echo "STARTED FOR path_size=${path_size} & k=${k} & i=${i}!"

            TEMPLATE="k10000_path_${path_size}_k${k}_${i}_path_${path_size}_k${k}_${i}"
            GWAS=${BIOGWAS_DATA_DIR}/3_pathways/extra_in_data/${TEMPLATE}_gwas.tsv
            OUT=${LSEA_TEST_DIR}/lsea_results/${TEMPLATE}

            python3 "${LSEA_DIR}/LSEA_2.4.py" \
                --input "${GWAS}" \
                --universe "${UNIVERSE}" \
                --out "${OUT}" \
                --plink_dir "${PLINK_DIR}" \
                --bfile "${BIOGWAS_BFILE}" \
                --column_names chr pos rsid pval \
                --clump_p1 "${PVAL}" \
                --print_all

            echo "FINISHED FOR path_size=${path_size} & k=${k} & i=${i}!"
        done
    done
done
