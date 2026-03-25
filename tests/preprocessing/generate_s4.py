#!/usr/bin/env python3
"""
Generate Supplementary Table S4 from O15_HYPTENSPREG comparison tables.
Merges 6 category-level TSVs into one unified table with a 'category' column.
"""

import os
import csv

CASE_STUDY_DIR = os.path.join(os.path.dirname(__file__), "..", "results", "case_study_o15")

CATEGORIES = ["bcm", "gte", "c2", "go_bp", "go_cc", "go_mf"]
CATEGORY_LABELS = {
    "bcm": "Blood Cell Markers",
    "gte": "GTEx Tissues",
    "c2": "C2 (MSigDB Curated)",
    "go_bp": "GO: Biological Process",
    "go_cc": "GO: Cellular Component",
    "go_mf": "GO: Molecular Function",
}

# Unified output columns
OUT_COLUMNS = [
    "category",
    "pathway",
    "lsea_qval",
    "lsea_pval",
    "lsea_loci",
    "magma_mean_pval",
    "magma_top_pval",
    "pascal_chi2_pval",
    "pascal_emp_pval",
]

def main():
    all_rows = []

    for cat in CATEGORIES:
        label = CATEGORY_LABELS[cat]
        fpath = os.path.join(
            CASE_STUDY_DIR,
            f"O15_HYPTENSPREG_comparison_{cat}.tsv",
        )
        if not os.path.exists(fpath):
            print(f"WARNING: missing {fpath}")
            continue

        with open(fpath) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                out_row = {
                    "category": label,
                    "pathway": row.get("pathway", ""),
                    "lsea_qval": row.get("lsea_qval", "NA"),
                    "lsea_pval": row.get("lsea_pval", "NA"),
                    "lsea_loci": row.get("lsea_loci", "NA"),
                    "magma_mean_pval": row.get("magma_mean_pval", "NA"),
                    "magma_top_pval": row.get("magma_top_pval", "NA"),
                    "pascal_chi2_pval": row.get("pascal_chi2_pval", "NA"),
                    "pascal_emp_pval": row.get("pascal_emp_pval", "NA"),
                }
                all_rows.append(out_row)

    # Write unified S4
    out_path = os.path.join(os.path.dirname(__file__), "Supplementary_Table_S4.tsv")
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"Generated: {out_path}")
    print(f"Total rows: {len(all_rows)}")
    for cat in CATEGORIES:
        label = CATEGORY_LABELS[cat]
        count = sum(1 for r in all_rows if r["category"] == label)
        print(f"  {label}: {count} pathways")


if __name__ == "__main__":
    main()
