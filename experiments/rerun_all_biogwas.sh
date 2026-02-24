#!/bin/bash
# ============================================================
# rerun_all_biogwas.sh — Full re-run of ALL bioGWAS experiments
# ============================================================
# 535 runs (continuous + binary + sensitivity) to validation_NEW_full/.
# Old results are NOT overwritten.
#
# Usage:
#   nohup bash experiments/rerun_all_biogwas.sh > rerun.log 2>&1 &
# ============================================================

set -euo pipefail

LSEA_DIR="$(cd "$(dirname "$0")/.." && pwd)"
ENV_FILE="${LSEA_DIR}/.env"
if [ ! -f "${ENV_FILE}" ]; then echo "ERROR: ${ENV_FILE} not found." >&2; exit 1; fi
source "${ENV_FILE}"
[ -f "${CONDA_ACTIVATE:-}" ] && source "${CONDA_ACTIVATE}"

echo "LSEA_DIR=${LSEA_DIR}"
echo "Started: $(date)"
echo "Git commit: $(git -C "${LSEA_DIR}" rev-parse HEAD)"
echo ""

UNIVERSE=${LSEA_TEST_DIR}/in_data/uni.json
PVAL=0.00000000729730

NEW_RESULTS=${LSEA_TEST_DIR}/validation_NEW_full
mkdir -p "${NEW_RESULTS}"

TOTAL=0
DONE=0

# Part 1: Continuous traits (200 runs)
# Iterations 0-29 are in in_data/, iterations 30-49 are in extra_in_data/
echo "===== Part 1: Continuous traits ====="
for i in $(seq 0 49); do
    for path_size in random big medium small; do
        TEMPLATE="test10000_path_${path_size}_${i}_path_${path_size}_${i}"
        # Try in_data first, then extra_in_data
        GWAS=${BIOGWAS_DATA_DIR}/3_pathways/in_data/${TEMPLATE}_gwas.tsv
        [ ! -f "${GWAS}" ] && GWAS=${BIOGWAS_DATA_DIR}/3_pathways/extra_in_data/${TEMPLATE}_gwas.tsv
        OUT=${NEW_RESULTS}/${TEMPLATE}
        TOTAL=$((TOTAL + 1))

        [ ! -f "${GWAS}" ] && echo "SKIP: ${TEMPLATE}" && continue

        echo "[$(date +%H:%M:%S)] ${TEMPLATE}"
        python3 "${LSEA_DIR}/LSEA_2.4.py" \
            --input "${GWAS}" --universe "${UNIVERSE}" --out "${OUT}" \
            --plink_dir "${PLINK_DIR}" --bfile "${BIOGWAS_BFILE}" \
            --column_names chr pos rsid pval --clump_p1 "${PVAL}"
        DONE=$((DONE + 1))
    done
done

# Part 2: Binary traits (200 runs)
echo "===== Part 2: Binary traits ====="
for i in $(seq 0 49); do
    for path_size in random big medium small; do
        TEMPLATE="bin10000_path_${path_size}_${i}_path_${path_size}_${i}"
        GWAS=${BIOGWAS_DATA_DIR}/3_pathways/binary_in_data/${TEMPLATE}_gwas.tsv
        OUT=${NEW_RESULTS}/${TEMPLATE}
        TOTAL=$((TOTAL + 1))

        [ ! -f "${GWAS}" ] && echo "SKIP: ${TEMPLATE}" && continue

        echo "[$(date +%H:%M:%S)] ${TEMPLATE}"
        python3 "${LSEA_DIR}/LSEA_2.4.py" \
            --input "${GWAS}" --universe "${UNIVERSE}" --out "${OUT}" \
            --plink_dir "${PLINK_DIR}" --bfile "${BIOGWAS_BFILE}" \
            --column_names chr pos rsid pval --clump_p1 "${PVAL}"
        DONE=$((DONE + 1))
    done
done

# Part 3: Sensitivity k=1..15 (135 runs)
echo "===== Part 3: Sensitivity ====="
for k in $(seq 1 15); do
    for i in $(seq 0 2); do
        for path_size in big medium small; do
            TEMPLATE="k10000_path_${path_size}_k${k}_${i}_path_${path_size}_k${k}_${i}"
            GWAS=${BIOGWAS_DATA_DIR}/3_pathways/extra_in_data/${TEMPLATE}_gwas.tsv
            OUT=${NEW_RESULTS}/${TEMPLATE}
            TOTAL=$((TOTAL + 1))

            [ ! -f "${GWAS}" ] && echo "SKIP: ${TEMPLATE}" && continue

            echo "[$(date +%H:%M:%S)] ${TEMPLATE}"
            python3 "${LSEA_DIR}/LSEA_2.4.py" \
                --input "${GWAS}" --universe "${UNIVERSE}" --out "${OUT}" \
                --plink_dir "${PLINK_DIR}" --bfile "${BIOGWAS_BFILE}" \
                --column_names chr pos rsid pval --clump_p1 "${PVAL}" \
                --print_all
            DONE=$((DONE + 1))
        done
    done
done

echo "=========================================="
echo "Finished: $(date)"
echo "Total: ${DONE}/${TOTAL}"
echo "Results in: ${NEW_RESULTS}"
