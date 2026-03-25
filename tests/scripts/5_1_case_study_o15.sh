#!/bin/bash
# ============================================================
# 5_1_case_study_o15.sh — FinnGen R9 Case Study (Experiment V)
# ============================================================
# Full pipeline: download → preprocess → LSEA × 6 → MAGMA × 2 × 6 → PASCAL
# Phenotype: O15_HYPTENSPREG (pregnancy hypertension)
# FinnGen R9: 14,727 cases / 196,143 controls
# ============================================================

set -euo pipefail

LSEA_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
ENV_FILE="${LSEA_DIR}/.env"
if [ ! -f "${ENV_FILE}" ]; then echo "ERROR: ${ENV_FILE} not found." >&2; exit 1; fi
source "${ENV_FILE}"
[ -f "${CONDA_ACTIVATE:-}" ] && source "${CONDA_ACTIVATE}"

# --- Required .env variables ---
: "${PLINK_DIR:?Set PLINK_DIR in .env}"
: "${MAGMA_BIN:?Set MAGMA_BIN in .env}"
: "${PASCAL_DIR:?Set PASCAL_DIR in .env}"
: "${PANUKB_BFILE:?Set PANUKB_BFILE in .env}"
: "${CHAIN_FILE:?Set CHAIN_FILE in .env}"
: "${LIFTOVER_BIN:?Set LIFTOVER_BIN in .env}"
: "${GENE_LOC:?Set GENE_LOC in .env}"
: "${PANUKB_LSEA_DIR:?Set PANUKB_LSEA_DIR in .env}"
: "${CASE_STUDY_DIR:?Set CASE_STUDY_DIR in .env}"
: "${GMT_C2_ALL:?Set GMT_C2_ALL in .env}"
: "${GMT_GTE:?Set GMT_GTE in .env}"
: "${GMT_BCM:?Set GMT_BCM in .env}"
: "${GMT_GO_BP:?Set GMT_GO_BP in .env}"
: "${GMT_GO_CC:?Set GMT_GO_CC in .env}"
: "${GMT_GO_MF:?Set GMT_GO_MF in .env}"

WORK_DIR="${CASE_STUDY_DIR}"
BFILE="${PANUKB_BFILE}"
UNI_DIR="${PANUKB_LSEA_DIR}/in_data"

# FinnGen O15 parameters
FINNGEN_URL="https://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_O15_HYPTENSPREG.gz"
PHENOTYPE="O15_HYPTENSPREG"
N_TOTAL=210870  # 14727 cases + 196143 controls

# LSEA p-value threshold (standard GWAS)
CLUMP_P1="5e-8"

# MAGMA models (summary-stat compatible)
MAGMA_MODELS=("mean" "top")

# Gene set categories
CATEGORIES=("c2" "gte" "bcm" "go_bp" "go_cc" "go_mf")

# GMT file mapping
declare -A GMT_MAP
GMT_MAP[c2]="${GMT_C2_ALL}"
GMT_MAP[gte]="${GMT_GTE}"
GMT_MAP[bcm]="${GMT_BCM}"
GMT_MAP[go_bp]="${GMT_GO_BP}"
GMT_MAP[go_cc]="${GMT_GO_CC}"
GMT_MAP[go_mf]="${GMT_GO_MF}"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"; }

# =============================================================================
# Step 0: Setup
# =============================================================================
log "=== Step 0: Setup ==="
mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}"
mkdir -p lsea_results magma_results pascal_results

# =============================================================================
# Step 1: Download FinnGen sumstats
# =============================================================================
log "=== Step 1: Download ==="

FINNGEN_GZ="${WORK_DIR}/finngen_R9_${PHENOTYPE}.gz"
if [ ! -f "${FINNGEN_GZ}" ]; then
    log "Downloading FinnGen R9 ${PHENOTYPE}..."
    wget -q "${FINNGEN_URL}" -O "${FINNGEN_GZ}"
    log "Downloaded: $(du -h ${FINNGEN_GZ} | cut -f1)"
else
    log "FinnGen file already exists: ${FINNGEN_GZ}"
fi

# =============================================================================
# Step 2: Preprocess (liftover + normalize)
# =============================================================================
log "=== Step 2: Preprocess ==="

OUT_PREFIX="${WORK_DIR}/${PHENOTYPE}"
if [ ! -f "${OUT_PREFIX}.norm.tsv" ]; then
    python3 "${LSEA_DIR}/tests/preprocessing/preproc_finngen.py" \
        "${FINNGEN_GZ}" \
        "${CHAIN_FILE}" \
        "${BFILE}.bim" \
        "${LIFTOVER_BIN}" \
        "${OUT_PREFIX}"
else
    log "Preprocessed file already exists: ${OUT_PREFIX}.norm.tsv"
fi

log "Variant counts:"
log "  LSEA input:  $(( $(wc -l < ${OUT_PREFIX}.norm.tsv) - 1 )) variants"
log "  MAGMA input: $(( $(wc -l < ${OUT_PREFIX}_magma_pval.tsv) - 1 )) variants"
log "  PASCAL input: $(wc -l < ${OUT_PREFIX}_pascal.txt) variants"

# =============================================================================
# Step 3: Run LSEA on 6 categories
# =============================================================================
log "=== Step 3: LSEA (6 categories) ==="

for cat in "${CATEGORIES[@]}"; do
    OUT_DIR="${WORK_DIR}/lsea_results/O15_${cat}"
    UNI_FILE="${UNI_DIR}/uni_${cat}.json"

    if [ -d "${OUT_DIR}" ] && ls "${OUT_DIR}"/annotation_stats_*.tsv 1>/dev/null 2>&1; then
        log "LSEA ${cat}: results already exist, skipping"
        continue
    fi

    log "Running LSEA: ${cat}..."
    python3 "${LSEA_DIR}/LSEA_2.4.py" \
        --input "${OUT_PREFIX}.norm.tsv" \
        --universe "${UNI_FILE}" \
        --out "${OUT_DIR}" \
        --plink_dir "${PLINK_DIR}" \
        --bfile "${BFILE}" \
        --column_names chr pos rsid pval \
        --clump_p1 "${CLUMP_P1}" \
        --print_all

    log "LSEA ${cat}: done"
done

# =============================================================================
# Step 4: Run MAGMA (2 models x 6 categories)
# =============================================================================
log "=== Step 4: MAGMA ==="

# Step 4a: Create MAGMA annotation (once, using 1000G EUR bim)
MAGMA_ANNO="${WORK_DIR}/magma_results/magma_anno_1kg"
if [ ! -f "${MAGMA_ANNO}.genes.annot" ]; then
    log "Creating MAGMA gene annotation..."
    ${MAGMA_BIN} \
        --annotate window=1,0.5 \
        --snp-loc "${BFILE}.bim" \
        --gene-loc "${GENE_LOC}" \
        --out "${MAGMA_ANNO}"
    log "MAGMA annotation done"
else
    log "MAGMA annotation already exists"
fi

# Step 4b: Convert GMT files to SSV format (MAGMA requirement)
for cat in "${CATEGORIES[@]}"; do
    GMT_FILE="${GMT_MAP[$cat]}"
    SSV_FILE="${WORK_DIR}/magma_results/${cat}.ssv"
    if [ ! -f "${SSV_FILE}" ] || [ ! -s "${SSV_FILE}" ]; then
        log "Converting GMT to SSV: ${cat}"
        # GMT format: NAME\tURL\tGENE1\tGENE2\t...
        # SSV format: NAME GENE1 GENE2 ... (skip URL column, filter sets with >=2 genes)
        awk -F'\t' 'NF >= 4 {printf $1; for(i=3;i<=NF;i++) printf " "$i; printf "\n"}' "${GMT_FILE}" > "${SSV_FILE}"
    fi
done

# Step 4c: Gene analysis + gene set analysis for each model x category
for model in "${MAGMA_MODELS[@]}"; do
    # Gene analysis (once per model)
    GENES_OUT="${WORK_DIR}/magma_results/${PHENOTYPE}_genes_${model}"

    if [ ! -f "${GENES_OUT}.genes.raw" ]; then
        log "MAGMA gene analysis: model=${model}..."
        ${MAGMA_BIN} \
            --bfile "${BFILE}" \
            --pval "${OUT_PREFIX}_magma_pval.tsv" use=rsid,pval \
            N=${N_TOTAL} \
            --gene-annot "${MAGMA_ANNO}.genes.annot" \
            --gene-model "snp-wise=${model}" \
            --out "${GENES_OUT}"
        log "MAGMA gene analysis (${model}): done"
    else
        log "MAGMA gene analysis (${model}): results exist, skipping"
    fi

    # Gene set analysis (for each category)
    for cat in "${CATEGORIES[@]}"; do
        SSV_FILE="${WORK_DIR}/magma_results/${cat}.ssv"
        SETS_OUT="${WORK_DIR}/magma_results/${PHENOTYPE}_sets_${model}_${cat}"

        if [ ! -f "${SETS_OUT}.gsa.out" ]; then
            log "MAGMA gene sets: model=${model}, category=${cat}..."
            ${MAGMA_BIN} \
                --gene-results "${GENES_OUT}.genes.raw" \
                --set-annot "${SSV_FILE}" \
                --out "${SETS_OUT}"
            log "MAGMA sets (${model}, ${cat}): done"
        else
            log "MAGMA sets (${model}, ${cat}): results exist, skipping"
        fi
    done
done

# =============================================================================
# Step 5: Run PASCAL
# =============================================================================
log "=== Step 5: PASCAL ==="

PASCAL_INPUT="${OUT_PREFIX}_pascal.txt"
PASCAL_OUT_DIR="${WORK_DIR}/pascal_results"

# PASCAL must be run from its own directory (it uses relative paths to resources/)
cd "${PASCAL_DIR}"

PASCAL_RESULT="${PASCAL_DIR}/output/${PHENOTYPE}_pascal.PathwaySet--msigBIOCARTA_KEGG_REACTOME--sum.txt"
if [ ! -f "${PASCAL_RESULT}" ]; then
    # Copy pascal input to PASCAL dir
    cp "${PASCAL_INPUT}" "${PASCAL_DIR}/${PHENOTYPE}_pascal.txt"

    log "Running PASCAL (failures are non-fatal — PASCAL may fail on sparse inputs)..."
    if ! ./Pascal --runpathway=on --pval="${PHENOTYPE}_pascal.txt" \
        > "${PASCAL_OUT_DIR}/${PHENOTYPE}_pascal.log" 2>&1; then
        log "WARNING: PASCAL exited with non-zero status. Check log: ${PASCAL_OUT_DIR}/${PHENOTYPE}_pascal.log"
    fi

    # Move results back
    if ls "${PASCAL_DIR}/output/${PHENOTYPE}_pascal"* 1>/dev/null 2>&1; then
        cp "${PASCAL_DIR}/output/${PHENOTYPE}_pascal"* "${PASCAL_OUT_DIR}/"
        log "PASCAL: done"
    else
        log "PASCAL: WARNING - no output files found. Check ${PASCAL_OUT_DIR}/${PHENOTYPE}_pascal.log"
    fi

    # Cleanup temp file
    rm -f "${PASCAL_DIR}/${PHENOTYPE}_pascal.txt"
else
    log "PASCAL: results already exist"
    cp "${PASCAL_RESULT}" "${PASCAL_OUT_DIR}/" 2>/dev/null || true
fi

cd "${WORK_DIR}"

# =============================================================================
# Step 6: Compile results
# =============================================================================
log "=== Step 6: Compile results ==="
python3 "${LSEA_DIR}/tests/preprocessing/compile_case_study.py" \
    --work_dir "${WORK_DIR}" \
    --phenotype "${PHENOTYPE}" \
    --categories ${CATEGORIES[@]}

log "=== Pipeline complete ==="
log "Results in: ${WORK_DIR}/"
log "  LSEA:    lsea_results/O15_*/"
log "  MAGMA:   magma_results/${PHENOTYPE}_sets_*"
log "  PASCAL:  pascal_results/"
log "  Summary: ${PHENOTYPE}_comparison_{c2,gte,bcm,...}.tsv + ${PHENOTYPE}_summary_stats.tsv"
