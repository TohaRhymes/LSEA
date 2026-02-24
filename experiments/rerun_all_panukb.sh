#!/bin/bash
# ============================================================
# rerun_all_panukb.sh — Full re-run of ALL Pan-UKB experiments
# ============================================================
# Steps:
#   1. Regenerate all 6 universes (compare with old)
#   2. Run enrichment on all 150 phenotypes x 6 categories (900 runs)
#   3. Run GWAS-on-GWAS enrichment (150 runs)
# Old results are NOT overwritten.
#
# Usage:
#   nohup bash experiments/rerun_all_panukb.sh > rerun_panukb.log 2>&1 &
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

PVAL=0.000000007479176476853146  # 0.05/6685228 (Bonferroni)
INTERVAL=500000

NEW_RESULTS=${PANUKB_LSEA_DIR}/validation_NEW_panukb
mkdir -p "${NEW_RESULTS}"

TOTAL=0
DONE=0

# ============================================================
# Part 1: Regenerate universes
# ============================================================
echo "===== Part 1: Regenerate universes ====="

VARIANTS=${PANUKB_LSEA_DIR}/in_data/variants.tsv
BED=${PANUKB_LSEA_DIR}/in_data/anno.bed

NEW_UNI_DIR=${NEW_RESULTS}/universes
mkdir -p "${NEW_UNI_DIR}"

declare -A GMT_FILES
GMT_FILES[c2]=${GMT_C2}
GMT_FILES[gte]=${GMT_GTE}
GMT_FILES[bcm]=${GMT_BCM}
GMT_FILES[go_bp]=${GMT_GO_BP}
GMT_FILES[go_cc]=${GMT_GO_CC}
GMT_FILES[go_mf]=${GMT_GO_MF}

for category in c2 gte bcm go_bp go_cc go_mf; do
    GMT=${GMT_FILES[$category]}
    OLD_UNIVERSE=${PANUKB_LSEA_DIR}/in_data/uni_${category}.json
    NEW_UNIVERSE=${NEW_UNI_DIR}/uni_${category}.json

    echo "[$(date +%H:%M:%S)] Creating universe: ${category}"

    python3 "${LSEA_DIR}/universe_generator.py" \
        --variants "${VARIANTS}" \
        --variants_colnames chr pos rsid \
        --features "${BED}" "${GMT}" \
        --interval "${INTERVAL}" \
        --out_json "${NEW_UNIVERSE}"

    # Compare with old
    python3 -c "
import json, sys
old = json.load(open('${OLD_UNIVERSE}'))
new = json.load(open('${NEW_UNIVERSE}'))
checks = [
    ('interval', old['interval'] == new['interval']),
    ('universe_intervals_number', old['universe_intervals_number'] == new['universe_intervals_number']),
    ('gene_set_dict count', len(old['gene_set_dict']) == len(new['gene_set_dict'])),
    ('features count', len(old['features']) == len(new['features'])),
    ('interval_counts', old['interval_counts'] == new['interval_counts']),
]
all_ok = all(ok for _, ok in checks)
status = 'IDENTICAL' if all_ok else 'DIFFERENT'
for name, ok in checks:
    print(f'  {\"OK\" if ok else \"DIFF\"}: {name}')
print(f'  Universe ${category}: {status}')
if not all_ok:
    sys.exit(1)
" || echo "WARNING: Universe ${category} differs from old!"

    echo ""
done

# ============================================================
# Part 2: Enrichment (Experiment 1) — 900 runs
# ============================================================
echo "===== Part 2: Enrichment (900 runs) ====="

NAMES=("c2" "gte" "bcm" "go_cc" "go_mf" "go_bp")
UNIVERSES=(
    "${NEW_UNI_DIR}/uni_c2.json"
    "${NEW_UNI_DIR}/uni_gte.json"
    "${NEW_UNI_DIR}/uni_bcm.json"
    "${NEW_UNI_DIR}/uni_go_cc.json"
    "${NEW_UNI_DIR}/uni_go_mf.json"
    "${NEW_UNI_DIR}/uni_go_bp.json"
)

for file in "${PANUKB_DATA_DIR}"/ukb_summstats/*.norm.tsv; do
    [ -f "${file}" ] || continue

    filename=$(basename "$file")
    TEMPLATE="${filename%.norm.tsv}"

    for i in "${!UNIVERSES[@]}"; do
        universe="${UNIVERSES[i]}"
        name="${NAMES[i]}"

        CUR_OUT_DIR="${NEW_RESULTS}/${TEMPLATE}_${name}"
        TOTAL=$((TOTAL + 1))

        echo "[$(date +%H:%M:%S)] ${TEMPLATE} x ${name}"

        python3 "${LSEA_DIR}/LSEA_2.4.py" \
            --input "${file}" \
            --universe "${universe}" \
            --out "${CUR_OUT_DIR}" \
            --plink_dir "${PLINK_DIR}" \
            --bfile "${PANUKB_BFILE}" \
            --column_names chr pos rsid pval \
            --clump_p1 "${PVAL}"
        DONE=$((DONE + 1))
    done
done

echo ""
echo "===== Part 2 done: ${DONE}/${TOTAL} ====="

# ============================================================
# Part 3: GWAS-on-GWAS (Experiment 2) — 150 runs
# ============================================================
echo "===== Part 3: GWAS-on-GWAS ====="

# Step 3.1: Collect BED files from Part 2 results
UKB_BED_DIR=${NEW_RESULTS}/ukb_universe_bed
mkdir -p "${UKB_BED_DIR}"

for dir in "${NEW_RESULTS}"/*/; do
    dir="${dir%/}"
    dirname=$(basename "$dir")
    [[ "$dirname" == "universes" || "$dirname" == "ukb_universe_bed" || "$dirname" == "ukb_universe" ]] && continue

    if [[ -f "$dir/merged_with_line_numbers.bed" ]]; then
        new_name="${dirname%_*}.bed"
        cp "$dir/merged_with_line_numbers.bed" "${UKB_BED_DIR}/${new_name}"
    fi
done

echo "Collected $(ls "${UKB_BED_DIR}"/*.bed 2>/dev/null | wc -l) BED files."

# Step 3.2: Create UKB phenotype universe
UKB_UNIVERSE_DIR=${NEW_RESULTS}/ukb_universe
mkdir -p "${UKB_UNIVERSE_DIR}"
UKB_UNIVERSE=${UKB_UNIVERSE_DIR}/uni_ukb.json

python3 "${LSEA_DIR}/universe_generator.py" \
    --variants "${VARIANTS}" \
    --variants_colnames chr pos rsid \
    --feature_files_dir "${UKB_BED_DIR}" \
    --interval "${INTERVAL}" \
    --out_json "${UKB_UNIVERSE}"

echo "UKB universe created: ${UKB_UNIVERSE}"

# Step 3.3: Run GWAS-on-GWAS enrichment
GWAS_TOTAL=0
GWAS_DONE=0

for file in "${PANUKB_DATA_DIR}"/ukb_summstats/*.norm.tsv; do
    [ -f "${file}" ] || continue

    filename=$(basename "$file")
    TEMPLATE="${filename%.norm.tsv}"

    CUR_OUT_DIR="${NEW_RESULTS}/${TEMPLATE}_ukb"
    GWAS_TOTAL=$((GWAS_TOTAL + 1))

    echo "[$(date +%H:%M:%S)] ${TEMPLATE} x ukb"

    python3 "${LSEA_DIR}/LSEA_2.4.py" \
        --input "${file}" \
        --universe "${UKB_UNIVERSE}" \
        --out "${CUR_OUT_DIR}" \
        --plink_dir "${PLINK_DIR}" \
        --bfile "${PANUKB_BFILE}" \
        --column_names chr pos rsid pval \
        --clump_p1 "${PVAL}" \
        --print_all
    GWAS_DONE=$((GWAS_DONE + 1))
done

echo ""
echo "=========================================="
echo "Finished: $(date)"
echo "Enrichment: ${DONE}/${TOTAL}"
echo "GWAS-on-GWAS: ${GWAS_DONE}/${GWAS_TOTAL}"
echo "Results in: ${NEW_RESULTS}"
