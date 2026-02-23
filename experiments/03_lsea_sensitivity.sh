#!/bin/bash
# ============================================================
# 03_lsea_sensitivity.sh — LSEA sensitivity analysis (varying k)
# ============================================================
# Runs LSEA on simulations with varying numbers of causal SNPs
# from the target pathway (k=1..15), 3 iterations x 3 pathway sizes.
# Total: 135 runs.
#
# Requires: universe from 00_create_universe.sh
# Original: lsea_test/iter_lsea_k.sh
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

# --- Run LSEA for each k, iteration, and pathway size ---
for k in $(seq 1 15); do
    for i in $(seq 0 2); do
        for path_size in big medium small; do
            echo "STARTED FOR path_size=${path_size} & k=${k} & i=${i}!"

            TEMPLATE="k10000_path_${path_size}_k${k}_${i}_path_${path_size}_k${k}_${i}"
            GWAS=${OLD_DATA_DIR}/3_pathways/extra_in_data/${TEMPLATE}_gwas.tsv
            OUT=${DATA_DIR}/lsea_results/${TEMPLATE}

            python3 "${LSEA_DIR}/LSEA_2.4.py" \
                --input "${GWAS}" \
                --universe "${UNIVERSE}" \
                --out "${OUT}" \
                --plink_dir "${PLINK_DIR}" \
                --bfile "${BFILE}" \
                --column_names chr pos rsid pval \
                --clump_p1 "${PVAL}" \
                --print_all

            echo "FINISHED FOR path_size=${path_size} & k=${k} & i=${i}!"
        done
    done
done
