#!/bin/bash
# ============================================================
# validate_single_runs.sh — Run ONE test per experiment type
# ============================================================
# Usage:
#   cd /path/to/LSEA
#   bash experiments/validate_single_runs.sh
# ============================================================

set -euo pipefail

LSEA_DIR="$(cd "$(dirname "$0")/.." && pwd)"
ENV_FILE="${LSEA_DIR}/.env"
if [ ! -f "${ENV_FILE}" ]; then echo "ERROR: ${ENV_FILE} not found." >&2; exit 1; fi
source "${ENV_FILE}"
[ -f "${CONDA_ACTIVATE:-}" ] && source "${CONDA_ACTIVATE}"

echo "LSEA_DIR=${LSEA_DIR}"
echo "Date: $(date)"
echo "Git commit: $(git -C "${LSEA_DIR}" rev-parse HEAD)"
echo ""

UNIVERSE=${LSEA_TEST_DIR}/in_data/uni.json
PVAL=0.00000000729730

VALIDATION_DIR=${LSEA_TEST_DIR}/validation_NEW_single
mkdir -p "${VALIDATION_DIR}"

PASS=0
FAIL=0

compare_stats() {
    local label="$1"
    local old_dir="$2"
    local new_dir="$3"

    echo "--- Comparing: ${label} ---"
    for old_stats in "${old_dir}"/annotation_stats_*.tsv; do
        [ -f "${old_stats}" ] || continue
        fname=$(basename "${old_stats}")
        new_stats="${new_dir}/${fname}"
        if [ ! -f "${new_stats}" ]; then
            echo "  MISSING: ${new_stats}"
            FAIL=$((FAIL + 1))
            continue
        fi
        if diff -q "${old_stats}" "${new_stats}" > /dev/null 2>&1; then
            echo "  OK: ${fname}"
            PASS=$((PASS + 1))
        else
            echo "  DIFF: ${fname}"
            diff "${old_stats}" "${new_stats}" || true
            FAIL=$((FAIL + 1))
        fi
    done
    echo ""
}

# ============================================================
# Test 1: Continuous trait
# ============================================================
echo "========== Test 1: Continuous trait =========="
TEMPLATE="test10000_path_small_0_path_small_0"
GWAS=${BIOGWAS_DATA_DIR}/3_pathways/in_data/${TEMPLATE}_gwas.tsv
OLD_OUT=${LSEA_TEST_DIR}/lsea_results/${TEMPLATE}
NEW_OUT=${VALIDATION_DIR}/${TEMPLATE}

python3 "${LSEA_DIR}/LSEA_2.4.py" \
    --input "${GWAS}" \
    --universe "${UNIVERSE}" \
    --out "${NEW_OUT}" \
    --plink_dir "${PLINK_DIR}" \
    --bfile "${BIOGWAS_BFILE}" \
    --column_names chr pos rsid pval \
    --clump_p1 "${PVAL}"

compare_stats "Continuous (path_small_0)" "${OLD_OUT}" "${NEW_OUT}"

# ============================================================
# Test 2: Binary trait
# ============================================================
echo "========== Test 2: Binary trait =========="
TEMPLATE="bin10000_path_small_0_path_small_0"
GWAS=${BIOGWAS_DATA_DIR}/3_pathways/binary_in_data/${TEMPLATE}_gwas.tsv
OLD_OUT=${LSEA_TEST_DIR}/lsea_results/${TEMPLATE}
NEW_OUT=${VALIDATION_DIR}/${TEMPLATE}

python3 "${LSEA_DIR}/LSEA_2.4.py" \
    --input "${GWAS}" \
    --universe "${UNIVERSE}" \
    --out "${NEW_OUT}" \
    --plink_dir "${PLINK_DIR}" \
    --bfile "${BIOGWAS_BFILE}" \
    --column_names chr pos rsid pval \
    --clump_p1 "${PVAL}"

compare_stats "Binary (path_small_0)" "${OLD_OUT}" "${NEW_OUT}"

# ============================================================
# Test 3: Sensitivity k=15
# ============================================================
echo "========== Test 3: Sensitivity k=15 =========="
TEMPLATE="k10000_path_small_k15_0_path_small_k15_0"
GWAS=${BIOGWAS_DATA_DIR}/3_pathways/extra_in_data/${TEMPLATE}_gwas.tsv
OLD_OUT=${LSEA_TEST_DIR}/lsea_results/${TEMPLATE}
NEW_OUT=${VALIDATION_DIR}/${TEMPLATE}

python3 "${LSEA_DIR}/LSEA_2.4.py" \
    --input "${GWAS}" \
    --universe "${UNIVERSE}" \
    --out "${NEW_OUT}" \
    --plink_dir "${PLINK_DIR}" \
    --bfile "${BIOGWAS_BFILE}" \
    --column_names chr pos rsid pval \
    --clump_p1 "${PVAL}" \
    --print_all

compare_stats "Sensitivity k=15 (path_small_0)" "${OLD_OUT}" "${NEW_OUT}"

# ============================================================
# Test 4: Universe re-generation
# ============================================================
echo "========== Test 4: Universe re-generation =========="
OLD_UNIVERSE=${LSEA_TEST_DIR}/in_data/uni.json
NEW_UNIVERSE=${VALIDATION_DIR}/uni_test.json

VARIANTS=${LSEA_TEST_DIR}/in_data/variants.tsv
BED=${LSEA_TEST_DIR}/in_data/anno.bed
GMT=${LSEA_TEST_DIR}/in_data/c2.cp.kegg.v2023.1.Hs.symbols.gmt

python3 "${LSEA_DIR}/universe_generator.py" \
    --variants "${VARIANTS}" \
    --variants_colnames chr pos rsid \
    --features "${BED}" "${GMT}" \
    --interval 500000 \
    --out_json "${NEW_UNIVERSE}"

python3 -c "
import json, sys
old = json.load(open('${OLD_UNIVERSE}'))
new = json.load(open('${NEW_UNIVERSE}'))
checks = [
    ('interval', old['interval'] == new['interval']),
    ('universe_intervals_number', old['universe_intervals_number'] == new['universe_intervals_number']),
    ('gene_set_dict count', len(old['gene_set_dict']) == len(new['gene_set_dict'])),
    ('features count', len(old['features']) == len(new['features'])),
    ('interval_counts count', len(old['interval_counts']) == len(new['interval_counts'])),
]
all_ok = True
for name, ok in checks:
    status = 'OK' if ok else 'DIFF'
    print(f'  {status}: {name}')
    if not ok:
        all_ok = False
if all_ok:
    print('  Universe: IDENTICAL')
else:
    print('  Universe: DIFFERENT')
    sys.exit(1)
"
echo ""

# ============================================================
# Summary
# ============================================================
echo "=========================================="
echo "Validation summary: PASS=${PASS} FAIL=${FAIL}"
if [ "${FAIL}" -eq 0 ]; then
    echo "ALL TESTS PASSED"
else
    echo "SOME TESTS FAILED"
fi
echo "Results in: ${VALIDATION_DIR}"
