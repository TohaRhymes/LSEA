#!/bin/bash
# ============================================================
# 5_3_case_study_parallel.sh — Multi-phenotype FinnGen Case Study
# ============================================================
# Three FinnGen R9 pregnancy phenotypes in one pipeline:
#   O15_HYPTENSPREG  — pregnancy hypertension (skips steps if data exists)
#   O24_GESTDIAB     — gestational diabetes mellitus
#   O14_PREECLAMPSIA — pre-eclampsia
#
# Sequential phases (each phase parallelised across phenotypes):
#   Phase 1: Preprocess  — parallel (download + liftover + bim filter)
#   Phase 2: PASCAL      — parallel (distinct output prefix per phenotype)
#   Phase 3: MAGMA       — parallel per phenotype (gene → gene-set analysis)
#   Phase 4: LSEA        — parallel (≤MAX_PARALLEL jobs, both p-cutoffs)
#   Phase 5: Compile     — sequential (BH-adjusted q-values for all tools)
#   Phase 6: Timing summary + predictions
#
# Measurement: /usr/bin/time -f "%e\t%M" -o FILE
#   %e = wall time (seconds)   %M = peak RSS (kB)
# ============================================================

set -euo pipefail

LSEA_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
ENV_FILE="${LSEA_DIR}/.env"
if [ ! -f "${ENV_FILE}" ]; then echo "ERROR: ${ENV_FILE} not found." >&2; exit 1; fi
source "${ENV_FILE}"
[ -f "${CONDA_ACTIVATE:-}" ] && source "${CONDA_ACTIVATE}"

: "${PLINK_DIR:?Set PLINK_DIR in .env}"
: "${MAGMA_BIN:?Set MAGMA_BIN in .env}"
: "${PASCAL_DIR:?Set PASCAL_DIR in .env}"
: "${PANUKB_BFILE:?Set PANUKB_BFILE in .env}"
: "${CHAIN_FILE:?Set CHAIN_FILE in .env}"
: "${LIFTOVER_BIN:?Set LIFTOVER_BIN in .env}"
: "${GENE_LOC:?Set GENE_LOC in .env}"
: "${CASE_STUDY_DIR:?Set CASE_STUDY_DIR in .env}"
: "${GMT_C2_ALL:?Set GMT_C2_ALL in .env}"
: "${GMT_GTE:?Set GMT_GTE in .env}"
: "${GMT_BCM:?Set GMT_BCM in .env}"
: "${GMT_GO_BP:?Set GMT_GO_BP in .env}"
: "${GMT_GO_CC:?Set GMT_GO_CC in .env}"
: "${GMT_GO_MF:?Set GMT_GO_MF in .env}"

# Optional — only required for LSEA
: "${PANUKB_LSEA_DIR:?Set PANUKB_LSEA_DIR in .env (dir with pre-built uni_*.json files)}"

# =============================================================================
# Configuration
# =============================================================================
BFILE="${PANUKB_BFILE}"
UNI_DIR="${PANUKB_LSEA_DIR}/in_data"   # pre-built universes (uni_c2.json, etc.)
MAX_PARALLEL=4
CLUMP_P1_VALUES=("5e-8" "1e-5")
MAGMA_MODELS=("mean" "top")
CATEGORIES=("c2" "gte" "bcm" "go_bp" "go_cc" "go_mf")

declare -A GMT_MAP
GMT_MAP[c2]="${GMT_C2_ALL}"
GMT_MAP[gte]="${GMT_GTE}"
GMT_MAP[bcm]="${GMT_BCM}"
GMT_MAP[go_bp]="${GMT_GO_BP}"
GMT_MAP[go_cc]="${GMT_GO_CC}"
GMT_MAP[go_mf]="${GMT_GO_MF}"

# =============================================================================
# Phenotype metadata
# =============================================================================
# N_TOTAL for O24/O14: FinnGen R9 total (342,499) as placeholder.
# Verify exact case counts at https://r9.finngen.fi/pheno/<CODE> and update.
CASE_STUDY_PARENT="$(dirname "${CASE_STUDY_DIR}")"

declare -A FINNGEN_URL N_TOTAL WORK_DIR LSEA_PREFIX

FINNGEN_URL[O15_HYPTENSPREG]="https://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_O15_HYPTENSPREG.gz"
N_TOTAL[O15_HYPTENSPREG]=210870
WORK_DIR[O15_HYPTENSPREG]="${CASE_STUDY_DIR}"          # existing dir, skip if data present
LSEA_PREFIX[O15_HYPTENSPREG]="O15"                     # existing lsea_results/O15_* dirs

FINNGEN_URL[GEST_DIABETES]="https://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_GEST_DIABETES.gz"
N_TOTAL[GEST_DIABETES]=210870   # 13039 cases + 197831 controls
WORK_DIR[GEST_DIABETES]="${CASE_STUDY_PARENT}/case_study_gest_diabetes"
LSEA_PREFIX[GEST_DIABETES]="GEST"

FINNGEN_URL[O15_PREECLAMPS]="https://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_O15_PREECLAMPS.gz"
N_TOTAL[O15_PREECLAMPS]=200929  # 6663 cases + 194266 controls
WORK_DIR[O15_PREECLAMPS]="${CASE_STUDY_PARENT}/case_study_preeclamps"
LSEA_PREFIX[O15_PREECLAMPS]="PREECLAMPS"

PHENOTYPES=("O15_HYPTENSPREG" "GEST_DIABETES" "O15_PREECLAMPS")

# =============================================================================
# Timing log
# =============================================================================
TIMING_LOG="${CASE_STUDY_PARENT}/case_study_timing.tsv"
mkdir -p "${CASE_STUDY_PARENT}"
echo -e "tool\tphenotype\tnote\twall_sec\tpeak_rss_mb" > "${TIMING_LOG}"

log()     { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
log_sep() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ===== $* ====="; }

# Append one timing measurement to TIMING_LOG (file-level atomic via temp+mv)
log_timing() {
    local tool="$1" pheno="$2" note="$3" tfile="$4"
    if [ ! -s "${tfile}" ]; then return; fi
    local wall rss_kb rss_mb
    wall=$(cut -f1 "${tfile}")
    rss_kb=$(cut -f2 "${tfile}")
    rss_mb=$(awk "BEGIN{printf \"%.0f\", ${rss_kb}/1024}")
    echo -e "${tool}\t${pheno}\t${note}\t${wall}\t${rss_mb}" >> "${TIMING_LOG}"
}

# Throttle background jobs to MAX_PARALLEL
wait_for_slot() {
    while [ "$(jobs -rp | wc -l)" -ge "${MAX_PARALLEL}" ]; do
        wait -n 2>/dev/null || true
    done
}

# =============================================================================
# Phase 1: Preprocessing (parallel)
# =============================================================================
log_sep "Phase 1: Preprocessing"

preprocess_phenotype() {
    local pheno="$1"
    local wdir="${WORK_DIR[$pheno]}"
    local url="${FINNGEN_URL[$pheno]}"
    mkdir -p "${wdir}/lsea_results" "${wdir}/magma_results" "${wdir}/pascal_results" "${wdir}/logs"

    local gz="${wdir}/finngen_R9_${pheno}.gz"
    local out_prefix="${wdir}/${pheno}"

    if [ -f "${out_prefix}.norm.tsv" ]; then
        log "${pheno}: preprocessed files exist, skipping"
        return
    fi

    if [ ! -f "${gz}" ]; then
        log "${pheno}: downloading..."
        wget -q "${url}" -O "${gz}"
        log "${pheno}: $(du -h "${gz}" | cut -f1)"
    fi

    log "${pheno}: liftover + bim filter..."
    local tfile; tfile=$(mktemp)
    /usr/bin/time -f "%e\t%M" -o "${tfile}" \
        python3 "${LSEA_DIR}/tests/preprocessing/preproc_finngen.py" \
            "${gz}" "${CHAIN_FILE}" "${BFILE}.bim" "${LIFTOVER_BIN}" "${out_prefix}" \
        > "${wdir}/logs/preproc.log" 2>&1
    log_timing "PREPROC" "${pheno}" "liftover+bim_filter" "${tfile}"
    rm -f "${tfile}"
    log "${pheno}: preprocessing done"
}

for pheno in "${PHENOTYPES[@]}"; do
    preprocess_phenotype "${pheno}" &
done
wait
log "All preprocessing complete"

for pheno in "${PHENOTYPES[@]}"; do
    wdir="${WORK_DIR[$pheno]}"
    norm="${wdir}/${pheno}.norm.tsv"
    mag="${wdir}/${pheno}_magma_pval.tsv"
    pas="${wdir}/${pheno}_pascal.txt"
    log "${pheno}: LSEA=$(( $(wc -l < "${norm}") - 1 )) | MAGMA=$(( $(wc -l < "${mag}") - 1 )) | PASCAL=$(wc -l < "${pas}")"
done

# =============================================================================
# Phase 2: PASCAL (parallel)
# =============================================================================
# Each phenotype has a distinct input-file prefix → outputs don't collide.
# All instances run from ${PASCAL_DIR} using a subshell cd.
log_sep "Phase 2: PASCAL"

run_pascal_one() {
    local pheno="$1"
    local wdir="${WORK_DIR[$pheno]}"
    local pascal_input="${wdir}/${pheno}_pascal.txt"
    local out_dir="${wdir}/pascal_results"
    local logfile="${wdir}/logs/pascal.log"
    local tfile; tfile=$(mktemp)

    # Skip if results already exist
    if ls "${PASCAL_DIR}/output/${pheno}_pascal."* 1>/dev/null 2>&1 || \
       ls "${out_dir}/${pheno}_pascal."* 1>/dev/null 2>&1; then
        log "${pheno}: PASCAL results exist, copying and skipping"
        cp "${PASCAL_DIR}/output/${pheno}_pascal."* "${out_dir}/" 2>/dev/null || true
        rm -f "${tfile}"
        return
    fi

    cp "${pascal_input}" "${PASCAL_DIR}/${pheno}_pascal.txt"
    log "${pheno}: PASCAL starting..."

    (
        cd "${PASCAL_DIR}"
        /usr/bin/time -f "%e\t%M" -o "${tfile}" \
            ./Pascal --runpathway=on --pval="${pheno}_pascal.txt"
    ) > "${logfile}" 2>&1 || true

    log_timing "PASCAL" "${pheno}" "default_genesets" "${tfile}"
    rm -f "${tfile}"

    if ls "${PASCAL_DIR}/output/${pheno}_pascal."* 1>/dev/null 2>&1; then
        cp "${PASCAL_DIR}/output/${pheno}_pascal."* "${out_dir}/"
        rm -f "${PASCAL_DIR}/${pheno}_pascal.txt"
        log "${pheno}: PASCAL done"
    else
        log "${pheno}: WARNING — PASCAL produced no output, check ${logfile}"
        rm -f "${PASCAL_DIR}/${pheno}_pascal.txt"
    fi
}

for pheno in "${PHENOTYPES[@]}"; do
    run_pascal_one "${pheno}" &
done
wait
log "All PASCAL runs complete"

# =============================================================================
# Phase 3: MAGMA (parallel per phenotype)
# =============================================================================
log_sep "Phase 3: MAGMA"

# SSV conversion and annotation (idempotent, done per phenotype)
for pheno in "${PHENOTYPES[@]}"; do
    wdir="${WORK_DIR[$pheno]}"
    for cat in "${CATEGORIES[@]}"; do
        ssv="${wdir}/magma_results/${cat}.ssv"
        if [ ! -s "${ssv}" ]; then
            awk -F'\t' 'NF >= 4 {printf $1; for(i=3;i<=NF;i++) printf " "$i; printf "\n"}' \
                "${GMT_MAP[$cat]}" > "${ssv}"
        fi
    done

    anno="${wdir}/magma_results/magma_anno_1kg"
    if [ ! -f "${anno}.genes.annot" ]; then
        log "${pheno}: MAGMA annotation..."
        "${MAGMA_BIN}" \
            --annotate window=1,0.5 \
            --snp-loc "${BFILE}.bim" \
            --gene-loc "${GENE_LOC}" \
            --out "${anno}" \
            > "${wdir}/logs/magma_anno.log" 2>&1
    fi
done

run_magma_one() {
    local pheno="$1"
    local wdir="${WORK_DIR[$pheno]}"
    local n="${N_TOTAL[$pheno]}"
    local anno="${wdir}/magma_results/magma_anno_1kg"
    local tfile; tfile=$(mktemp)

    for model in "${MAGMA_MODELS[@]}"; do
        local genes_out="${wdir}/magma_results/${pheno}_genes_${model}"

        if [ ! -f "${genes_out}.genes.raw" ]; then
            log "${pheno}: MAGMA gene analysis (${model})..."
            /usr/bin/time -f "%e\t%M" -o "${tfile}" \
                "${MAGMA_BIN}" \
                    --bfile "${BFILE}" \
                    --pval "${wdir}/${pheno}_magma_pval.tsv" use=rsid,pval \
                    N="${n}" \
                    --gene-annot "${anno}.genes.annot" \
                    --gene-model "snp-wise=${model}" \
                    --out "${genes_out}" \
                > "${wdir}/logs/magma_gene_${model}.log" 2>&1
            log_timing "MAGMA_gene" "${pheno}" "${model}" "${tfile}"
        else
            log "${pheno}: MAGMA gene (${model}) exists, skipping"
        fi

        for cat in "${CATEGORIES[@]}"; do
            local sets_out="${wdir}/magma_results/${pheno}_sets_${model}_${cat}"
            if [ ! -f "${sets_out}.gsa.out" ]; then
                log "${pheno}: MAGMA gene-sets (${model}, ${cat})..."
                /usr/bin/time -f "%e\t%M" -o "${tfile}" \
                    "${MAGMA_BIN}" \
                        --gene-results "${genes_out}.genes.raw" \
                        --set-annot "${wdir}/magma_results/${cat}.ssv" \
                        --out "${sets_out}" \
                    > "${wdir}/logs/magma_sets_${model}_${cat}.log" 2>&1
                log_timing "MAGMA_sets" "${pheno}" "${model}_${cat}" "${tfile}"
            fi
        done

        # Free disk: genes.raw (large) and per-gene set results
        rm -f "${genes_out}.genes.raw"
        rm -f "${wdir}/magma_results/${pheno}_sets_${model}_"*.gsa.genes.out
    done

    rm -f "${tfile}"
    log "${pheno}: MAGMA complete"
}

for pheno in "${PHENOTYPES[@]}"; do
    run_magma_one "${pheno}" &
done
wait
log "All MAGMA runs complete"

# =============================================================================
# Phase 4: LSEA (parallel, ≤MAX_PARALLEL jobs, both p-cutoffs)
# =============================================================================
log_sep "Phase 4: LSEA"

run_lsea_one() {
    local pheno="$1" cat="$2" p1="$3"
    local wdir="${WORK_DIR[$pheno]}"
    local prefix="${LSEA_PREFIX[$pheno]}"
    local out_dir="${wdir}/lsea_results/${prefix}_${cat}"
    local uni_file="${UNI_DIR}/uni_${cat}.json"
    local tfile; tfile=$(mktemp)

    if [ ! -f "${uni_file}" ]; then
        log "WARN: ${uni_file} not found — skipping ${pheno} ${cat}"
        rm -f "${tfile}"; return
    fi

    # Normalize p-value to Python float format (5e-8 → 5e-08) for filename matching
    local p1_norm
    p1_norm=$(python3 -c "print(float('${p1}'))")
    if [ -f "${out_dir}/uni_${cat}_result_${p1_norm}.tsv" ]; then
        log "${pheno} ${cat} p=${p1}: results exist, skipping"
        rm -f "${tfile}"; return
    fi

    log "${pheno} ${cat} p=${p1}: LSEA..."
    /usr/bin/time -f "%e\t%M" -o "${tfile}" \
        python3 "${LSEA_DIR}/LSEA_2.4.py" \
            --input "${wdir}/${pheno}.norm.tsv" \
            --universe "${uni_file}" \
            --out "${out_dir}" \
            --plink_dir "${PLINK_DIR}" \
            --bfile "${BFILE}" \
            --column_names chr pos rsid pval \
            --clump_p1 "${p1}" \
            --print_all \
        > "${wdir}/logs/lsea_${cat}_p${p1}.log" 2>&1
    log_timing "LSEA" "${pheno}" "${cat}_p${p1}" "${tfile}"
    rm -f "${tfile}"
    log "${pheno} ${cat} p=${p1}: LSEA done"
}

for pheno in "${PHENOTYPES[@]}"; do
    for cat in "${CATEGORIES[@]}"; do
        for p1 in "${CLUMP_P1_VALUES[@]}"; do
            wait_for_slot
            run_lsea_one "${pheno}" "${cat}" "${p1}" &
        done
    done
done
wait
log "All LSEA runs complete"

# =============================================================================
# Phase 5: Compile results (sequential)
# =============================================================================
log_sep "Phase 5: Compile results"

for pheno in "${PHENOTYPES[@]}"; do
    wdir="${WORK_DIR[$pheno]}"
    prefix="${LSEA_PREFIX[$pheno]}"
    log "${pheno}: compiling..."
    python3 "${LSEA_DIR}/tests/preprocessing/compile_case_study.py" \
        --work_dir "${wdir}" \
        --phenotype "${pheno}" \
        --lsea_prefix "${prefix}" \
        --categories "${CATEGORIES[@]}" \
        > "${wdir}/logs/compile.log" 2>&1
    log "${pheno}: → ${wdir}/${pheno}_comparison_*.tsv"
done

# =============================================================================
# Phase 6: Timing summary + benchmark predictions
# =============================================================================
log_sep "Phase 6: Timing summary"

echo ""
echo "=== TIMING RESULTS ==="
column -t -s $'\t' "${TIMING_LOG}"

echo ""
echo "=== BENCHMARK PREDICTIONS (150 phenotypes × 6 gene sets) ==="
python3 - "${TIMING_LOG}" <<'PYEOF'
import sys, math, statistics

timing_log = sys.argv[1]
data = {}  # key → list of (wall_sec, rss_mb)

with open(timing_log) as f:
    next(f)
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 5:
            continue
        tool, pheno, note, wall, rss = parts
        try:
            w, r = float(wall), float(rss)
        except ValueError:
            continue
        key = tool
        data.setdefault(key, []).append((w, r))

def avg(lst): return sum(lst)/len(lst) if lst else 0
def fmt_days(sec): return f"{sec/86400:.1f} days"

print(f"\n{'Tool':<22} {'Avg wall':<12} {'Peak RSS':<12} {'150 pheno @1T':<17} {'150 pheno @4T':<14} {'Note'}")
print("-" * 90)

N_PHENO, N_CATS, N_P = 150, 6, 2   # 2 p-cutoffs for LSEA

# Preprocessing: one run per phenotype (download excluded — measured separately)
if "PREPROC" in data:
    vals = data["PREPROC"]
    avg_w = avg([v[0] for v in vals])
    max_r = max(v[1] for v in vals)
    t1 = N_PHENO * avg_w
    t4 = t1 / 4
    print(f"{'Preproc (liftover+bim)':<22} {avg_w:>8.0f}s   {max_r:>6.0f} MB   "
          f"{fmt_days(t1):<17} {fmt_days(t4):<14} once per phenotype")

# PASCAL: one run per phenotype covers default gene sets
if "PASCAL" in data:
    vals = data["PASCAL"]
    avg_w = avg([v[0] for v in vals])
    max_r = max(v[1] for v in vals)
    t1 = N_PHENO * avg_w
    t4 = t1 / 4
    print(f"{'PASCAL (per pheno)':<22} {avg_w:>8.0f}s   {max_r:>6.0f} MB   "
          f"{fmt_days(t1):<15} {fmt_days(t4)}")

# MAGMA gene analysis (bottleneck)
if "MAGMA_gene" in data:
    vals = data["MAGMA_gene"]
    avg_w = avg([v[0] for v in vals])
    max_r = max(v[1] for v in vals)
    # 150 pheno × 2 models
    t1 = N_PHENO * 2 * avg_w
    t4 = t1 / 4
    print(f"{'MAGMA gene (2 models)':<22} {avg_w:>8.0f}s   {max_r:>6.0f} MB   "
          f"{fmt_days(t1):<15} {fmt_days(t4)}")

# MAGMA gene-sets
if "MAGMA_sets" in data:
    vals = data["MAGMA_sets"]
    avg_w = avg([v[0] for v in vals])
    max_r = max(v[1] for v in vals)
    t1 = N_PHENO * 2 * N_CATS * avg_w   # 150 × 2 models × 6 cats
    t4 = t1 / 4
    print(f"{'MAGMA sets (2m × 6c)':<22} {avg_w:>8.0f}s   {max_r:>6.0f} MB   "
          f"{fmt_days(t1):<15} {fmt_days(t4)}")

# LSEA
if "LSEA" in data:
    vals = data["LSEA"]
    avg_w = avg([v[0] for v in vals])
    max_r = max(v[1] for v in vals)
    t1 = N_PHENO * N_CATS * N_P * avg_w   # 150 × 6 cats × 2 p-cutoffs
    t4 = t1 / 4
    print(f"{'LSEA (6c × 2 p-cuts)':<22} {avg_w:>8.0f}s   {max_r:>6.0f} MB   "
          f"{fmt_days(t1):<15} {fmt_days(t4)}")

print("\n--- RAM budget @ 4 parallel processes ---")
for key, label in [("PASCAL","PASCAL"), ("MAGMA_gene","MAGMA"), ("LSEA","LSEA")]:
    if key in data:
        max_r = max(v[1] for v in data[key])
        total_gb = 4 * max_r / 1024
        print(f"  {label:<10}: {max_r:.0f} MB/proc × 4 = {total_gb:.1f} GB")

PYEOF

log_sep "Pipeline complete"
log "Results:  ${CASE_STUDY_PARENT}/case_study_*/"
log "Timing:   ${TIMING_LOG}"
