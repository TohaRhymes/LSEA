#!/bin/bash
# ============================================================
# 02_lsea_binary.sh — LSEA on binary trait simulations
# ============================================================
# Runs LSEA enrichment analysis on 200 simulated binary-trait GWAS
# (50 iterations x 4 pathway sizes: small, medium, big, random).
#
# Requires: universe from 00_create_universe.sh
# Original: lsea_test/iter_lsea_binary.sh
# Updated:  CLI flags changed to --long_flag format
# ============================================================

set -euo pipefail

# --- Paths (adjust for your server) ---
LSEA_DIR="$(cd "$(dirname "$0")/.." && pwd)"
PLINK_DIR=/home/achangalidi/tools/plink
OLD_DATA_DIR=/media/DATA/gwasim/round2/bioGWAS/tests
DATA_DIR=/media/DATA/gwasim/round2/lsea_test

BFILE=${OLD_DATA_DIR}/3_pathways/in_data/test10000_filt_sim
UNIVERSE=${DATA_DIR}/in_data/uni.json
PVAL=0.00000000729730

# --- Run LSEA for each pathway size and iteration ---
for i in $(seq 0 49); do
    for path_size in random big medium small; do
        echo "STARTED FOR path_size=${path_size} & i=${i}!"

        TEMPLATE="bin10000_path_${path_size}_${i}_path_${path_size}_${i}"
        GWAS=${OLD_DATA_DIR}/3_pathways/binary_in_data/${TEMPLATE}_gwas.tsv
        OUT=${DATA_DIR}/lsea_results/${TEMPLATE}

        python3 "${LSEA_DIR}/LSEA_2.4.py" \
            --input "${GWAS}" \
            --universe "${UNIVERSE}" \
            --out "${OUT}" \
            --plink_dir "${PLINK_DIR}" \
            --bfile "${BFILE}" \
            --column_names chr pos rsid pval \
            --clump_p1 "${PVAL}"

        echo "FINISHED FOR path_size=${path_size} & i=${i}!"
    done
done
