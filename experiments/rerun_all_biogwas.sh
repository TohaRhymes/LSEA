#!/bin/bash
# ============================================================
# rerun_all_biogwas.sh — Full re-run of ALL bioGWAS experiments
# ============================================================
# Runs all 535 bioGWAS experiments (continuous + binary + sensitivity)
# to a NEW output directory (validation_NEW_full/) so old results
# are NOT overwritten.
#
# Usage:
#   cd /path/to/LSEA
#   nohup bash experiments/rerun_all_biogwas.sh > rerun.log 2>&1 &
#
# Takes several hours. Check progress with: tail -f rerun.log
# ============================================================

set -euo pipefail

LSEA_DIR="$(cd "$(dirname "$0")/.." && pwd)"
echo "LSEA_DIR=${LSEA_DIR}"
echo "Started: $(date)"
echo "Git commit: $(git -C "${LSEA_DIR}" rev-parse HEAD)"
echo ""

# --- Common paths ---
PLINK_DIR=/home/achangalidi/tools/plink
OLD_DATA_DIR=/media/DATA/gwasim/round2/bioGWAS/tests
DATA_DIR=/media/DATA/gwasim/round2/lsea_test
BFILE=${OLD_DATA_DIR}/3_pathways/in_data/test10000_filt_sim
UNIVERSE=${DATA_DIR}/in_data/uni.json
PVAL=0.00000000729730

# Output to a SEPARATE directory — old results untouched
NEW_RESULTS=${DATA_DIR}/validation_NEW_full
mkdir -p "${NEW_RESULTS}"

TOTAL=0
DONE=0

# ============================================================
# Part 1: Continuous traits (200 runs)
# ============================================================
echo "===== Part 1: Continuous traits (200 runs) ====="
for i in $(seq 0 49); do
    for path_size in random big medium small; do
        TEMPLATE="test10000_path_${path_size}_${i}_path_${path_size}_${i}"
        GWAS=${OLD_DATA_DIR}/3_pathways/extra_in_data/${TEMPLATE}_gwas.tsv
        OUT=${NEW_RESULTS}/${TEMPLATE}
        TOTAL=$((TOTAL + 1))

        if [ ! -f "${GWAS}" ]; then
            echo "SKIP (no GWAS file): ${TEMPLATE}"
            continue
        fi

        echo "[$(date +%H:%M:%S)] Running: ${TEMPLATE}"
        python3 "${LSEA_DIR}/LSEA_2.4.py" \
            --input "${GWAS}" \
            --universe "${UNIVERSE}" \
            --out "${OUT}" \
            --plink_dir "${PLINK_DIR}" \
            --bfile "${BFILE}" \
            --column_names chr pos rsid pval \
            --clump_p1 "${PVAL}"
        DONE=$((DONE + 1))
    done
done
echo "Continuous done: ${DONE}/${TOTAL}"
echo ""

# ============================================================
# Part 2: Binary traits (200 runs)
# ============================================================
echo "===== Part 2: Binary traits (200 runs) ====="
for i in $(seq 0 49); do
    for path_size in random big medium small; do
        TEMPLATE="bin10000_path_${path_size}_${i}_path_${path_size}_${i}"
        GWAS=${OLD_DATA_DIR}/3_pathways/binary_in_data/${TEMPLATE}_gwas.tsv
        OUT=${NEW_RESULTS}/${TEMPLATE}
        TOTAL=$((TOTAL + 1))

        if [ ! -f "${GWAS}" ]; then
            echo "SKIP (no GWAS file): ${TEMPLATE}"
            continue
        fi

        echo "[$(date +%H:%M:%S)] Running: ${TEMPLATE}"
        python3 "${LSEA_DIR}/LSEA_2.4.py" \
            --input "${GWAS}" \
            --universe "${UNIVERSE}" \
            --out "${OUT}" \
            --plink_dir "${PLINK_DIR}" \
            --bfile "${BFILE}" \
            --column_names chr pos rsid pval \
            --clump_p1 "${PVAL}"
        DONE=$((DONE + 1))
    done
done
echo "Binary done: ${DONE}/${TOTAL}"
echo ""

# ============================================================
# Part 3: Sensitivity k=1..15 (135 runs)
# ============================================================
echo "===== Part 3: Sensitivity k=1..15 (135 runs) ====="
for k in $(seq 1 15); do
    for i in $(seq 0 2); do
        for path_size in big medium small; do
            TEMPLATE="k10000_path_${path_size}_k${k}_${i}_path_${path_size}_k${k}_${i}"
            GWAS=${OLD_DATA_DIR}/3_pathways/extra_in_data/${TEMPLATE}_gwas.tsv
            OUT=${NEW_RESULTS}/${TEMPLATE}
            TOTAL=$((TOTAL + 1))

            if [ ! -f "${GWAS}" ]; then
                echo "SKIP (no GWAS file): ${TEMPLATE}"
                continue
            fi

            echo "[$(date +%H:%M:%S)] Running: ${TEMPLATE}"
            python3 "${LSEA_DIR}/LSEA_2.4.py" \
                --input "${GWAS}" \
                --universe "${UNIVERSE}" \
                --out "${OUT}" \
                --plink_dir "${PLINK_DIR}" \
                --bfile "${BFILE}" \
                --column_names chr pos rsid pval \
                --clump_p1 "${PVAL}" \
                --print_all
            DONE=$((DONE + 1))
        done
    done
done
echo "Sensitivity done: ${DONE}/${TOTAL}"
echo ""

# ============================================================
# Summary
# ============================================================
echo "=========================================="
echo "Finished: $(date)"
echo "Total: ${DONE}/${TOTAL} experiments completed"
echo "Results in: ${NEW_RESULTS}"
echo ""
echo "To compare with old results:"
echo "  python3 experiments/compare_results.py ${DATA_DIR}/lsea_results ${NEW_RESULTS}"
