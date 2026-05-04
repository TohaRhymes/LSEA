#!/bin/bash
# ============================================================
# 08_case_study_gwas_on_gwas.sh — Cross-trait Enrichment for O15
# ============================================================
# GWAS-on-GWAS analysis: O15_HYPTENSPREG (FinnGen R9) tested
# against a UKB-derived universe of 129 phenotypes.
# This is a standalone script — does NOT require re-running
# 5_1_case_study_o15.sh.
# ============================================================

set -euo pipefail

LSEA_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
ENV_FILE="${LSEA_DIR}/.env"
if [ ! -f "${ENV_FILE}" ]; then echo "ERROR: ${ENV_FILE} not found." >&2; exit 1; fi
source "${ENV_FILE}"
[ -f "${CONDA_ACTIVATE:-}" ] && source "${CONDA_ACTIVATE}"

# --- Required .env variables ---
: "${PLINK_DIR:?Set PLINK_DIR in .env}"
: "${PANUKB_BFILE:?Set PANUKB_BFILE in .env}"
: "${CASE_STUDY_DIR:?Set CASE_STUDY_DIR in .env}"
: "${UKB_UNIVERSE:?Set UKB_UNIVERSE in .env}"

WORK_DIR="${CASE_STUDY_DIR}"
BFILE="${PANUKB_BFILE}"
PHENOTYPE="O15_HYPTENSPREG"
CLUMP_P1="5e-8"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"; }

# =============================================================================
# Step 1: Verify input exists
# =============================================================================
log "=== Step 1: Verify inputs ==="

INPUT_FILE="${WORK_DIR}/${PHENOTYPE}.norm.tsv"
if [ ! -f "${INPUT_FILE}" ]; then
    echo "ERROR: ${INPUT_FILE} not found. Run 5_1_case_study_o15.sh first." >&2
    exit 1
fi
log "Input: ${INPUT_FILE}"
log "Universe: ${UKB_UNIVERSE}"

if [ ! -f "${UKB_UNIVERSE}" ]; then
    echo "ERROR: UKB universe not found: ${UKB_UNIVERSE}" >&2
    exit 1
fi

# =============================================================================
# Step 2: Run LSEA with UKB universe
# =============================================================================
log "=== Step 2: LSEA (UKB GWAS-on-GWAS) ==="

OUT_DIR="${WORK_DIR}/lsea_results/O15_ukb"

if [ -d "${OUT_DIR}" ] && ls "${OUT_DIR}"/uni_ukb_result_*.tsv 1>/dev/null 2>&1; then
    log "LSEA UKB: results already exist, skipping"
else
    log "Running LSEA: O15 vs UKB universe (129 phenotypes)..."
    python3 "${LSEA_DIR}/LSEA_2.4.py" \
        --input "${INPUT_FILE}" \
        --universe "${UKB_UNIVERSE}" \
        --out "${OUT_DIR}" \
        --plink_dir "${PLINK_DIR}" \
        --bfile "${BFILE}" \
        --column_names chr pos rsid pval \
        --clump_p1 "${CLUMP_P1}" \
        --print_all
    log "LSEA UKB: done"
fi

# =============================================================================
# Step 3: Compile GWAS-on-GWAS results
# =============================================================================
log "=== Step 3: Compile GWAS-on-GWAS results ==="

python3 "${LSEA_DIR}/tests/preprocessing/compile_case_study.py" \
    --work_dir "${WORK_DIR}" \
    --phenotype "${PHENOTYPE}" \
    --categories c2 gte bcm go_bp go_cc go_mf \
    --gwas_on_gwas_dir "${OUT_DIR}" \
    --manifest "${PANUKB_DATA_DIR:-}/Pan-UK Biobank phenotype manifest - phenotype_manifest.tsv"

log "=== Pipeline complete ==="
log "Results in: ${WORK_DIR}/"
log "  GWAS-on-GWAS: ${OUT_DIR}/"
log "  Summary: ${WORK_DIR}/${PHENOTYPE}_gwas_on_gwas.tsv"
