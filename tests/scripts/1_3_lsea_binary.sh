#!/bin/bash
# ============================================================
# 02_lsea_binary.sh — LSEA on binary trait simulations
# ============================================================
# 200 runs: 50 iterations x 4 pathway sizes (small, medium, big, random).
# Requires: universe from 00_create_universe.sh
# Original: lsea_test/iter_lsea_binary.sh
# ============================================================

set -euo pipefail

LSEA_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
ENV_FILE="${LSEA_DIR}/.env"
if [ ! -f "${ENV_FILE}" ]; then echo "ERROR: ${ENV_FILE} not found." >&2; exit 1; fi
source "${ENV_FILE}"
[ -f "${CONDA_ACTIVATE:-}" ] && source "${CONDA_ACTIVATE}"

UNIVERSE=${LSEA_TEST_DIR}/in_data/uni.json
PVAL=0.00000000729730

for i in $(seq 0 49); do
    for path_size in random big medium small; do
        echo "STARTED FOR path_size=${path_size} & i=${i}!"

        TEMPLATE="bin10000_path_${path_size}_${i}_path_${path_size}_${i}"
        GWAS=${BIOGWAS_DATA_DIR}/3_pathways/binary_in_data/${TEMPLATE}_gwas.tsv
        OUT=${LSEA_TEST_DIR}/lsea_results/${TEMPLATE}

        python3 "${LSEA_DIR}/LSEA_2.4.py" \
            --input "${GWAS}" \
            --universe "${UNIVERSE}" \
            --out "${OUT}" \
            --plink_dir "${PLINK_DIR}" \
            --bfile "${BIOGWAS_BFILE}" \
            --column_names chr pos rsid pval \
            --clump_p1 "${PVAL}"

        echo "FINISHED FOR path_size=${path_size} & i=${i}!"
    done
done
