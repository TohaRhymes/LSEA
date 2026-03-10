#!/usr/bin/env python3
"""
Compile LSEA + MAGMA + PASCAL results for O15_HYPTENSPREG case study.

Outputs:
1. Per-category comparison tables (LSEA q-val vs MAGMA p-vals vs PASCAL p-val)
2. Summary statistics table
3. Highlighted biologically relevant pathways
"""

import argparse
import csv
import os
import sys
from collections import defaultdict


def parse_lsea_results(result_dir, category):
    """Parse LSEA result TSV file for a given category."""
    results = {}

    # Find the result file (pattern: uni_{cat}_result_{p_cutoff}.tsv)
    for fname in os.listdir(result_dir):
        if fname.startswith(f"uni_{category}_result_") and fname.endswith(".tsv"):
            filepath = os.path.join(result_dir, fname)
            with open(filepath, "r") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    gene_set = row["gene_set"]
                    results[gene_set] = {
                        "lsea_pval": float(row["p_value"]),
                        "lsea_qval": float(row["q_value"]),
                        "lsea_loci": int(row["overlapping_loci"]),
                        "lsea_sig": row["significance"],
                    }
            break

    return results


def parse_lsea_stats(result_dir, category):
    """Parse LSEA annotation_stats file, cross-validated against actual results."""
    import re
    stats_file = os.path.join(result_dir, f"annotation_stats_uni_{category}.tsv")
    stats = {}
    if os.path.exists(stats_file):
        with open(stats_file, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                stats[float(row["p_cutoff"])] = {
                    "num_loci": int(row["num_loci"]),
                    "annotated_loci": int(row["annotated_loci"]),
                    "significant_hits": int(row["significant_hits"]),
                    "min_qval": float(row["min_qval"]),
                }

    # Cross-validate against actual result files: annotation_stats can be stale
    # (e.g. annotated_loci=0 while results show overlaps). Recompute from results.
    for fname in os.listdir(result_dir):
        if not (fname.startswith(f"uni_{category}_result_") and fname.endswith(".tsv")):
            continue
        p_cutoff_str = fname.replace(f"uni_{category}_result_", "").replace(".tsv", "")
        try:
            p_cutoff = float(p_cutoff_str)
        except ValueError:
            continue

        all_loci = set()
        sig_hits = 0
        min_q = 1.0
        with open(os.path.join(result_dir, fname)) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                q = float(row["q_value"])
                min_q = min(min_q, q)
                if row["significance"].lower() == "true":
                    sig_hits += 1
                # Count unique annotated loci from description column
                desc = row.get("description", "{}")
                loci_keys = re.findall(r"'(\d+:[0-9]+-[0-9]+)'", desc)
                all_loci.update(loci_keys)

        if p_cutoff not in stats:
            stats[p_cutoff] = {"num_loci": 0, "annotated_loci": 0,
                                "significant_hits": 0, "min_qval": 1.0}
        # Override stale values from annotation_stats
        stats[p_cutoff]["annotated_loci"] = len(all_loci)
        stats[p_cutoff]["significant_hits"] = sig_hits
        stats[p_cutoff]["min_qval"] = min_q

    return stats


def parse_magma_gsa(filepath):
    """Parse MAGMA .gsa.out file."""
    results = {}
    if not os.path.exists(filepath):
        return results

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or line.startswith("VARIABLE"):
                continue
            if not line:
                continue

            # MAGMA gsa.out is fixed-width with FULL_NAME at end
            parts = line.split()
            if len(parts) < 7:
                continue

            # MAGMA gsa.out format:
            #   With FULL_NAME (truncated VARIABLE): VARIABLE TYPE NGENES BETA BETA_STD SE P FULL_NAME
            #   Without FULL_NAME (short names):     VARIABLE TYPE NGENES BETA BETA_STD SE P
            # Detect by checking if last column is numeric (P) or string (FULL_NAME)
            try:
                pval = float(parts[-1])
                # Last column is P (no FULL_NAME) — use VARIABLE as name
                full_name = parts[0]
            except ValueError:
                # Last column is FULL_NAME — P is second-to-last
                full_name = parts[-1]
                try:
                    pval = float(parts[-2])
                except (ValueError, IndexError):
                    continue

            results[full_name] = {"magma_pval": pval}

    return results


def parse_pascal_results(filepath):
    """Parse PASCAL pathway results file."""
    results = {}
    if not os.path.exists(filepath):
        return results

    with open(filepath, "r") as f:
        header = f.readline()  # Name chi2Pvalue empPvalue
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            name = parts[0]
            try:
                chi2_pval = float(parts[1])
                emp_pval = float(parts[2])
            except ValueError:
                continue
            results[name] = {
                "pascal_chi2_pval": chi2_pval,
                "pascal_emp_pval": emp_pval,
            }

    return results


# Biologically relevant keywords for pregnancy hypertension
RELEVANT_KEYWORDS = [
    "hemostasis", "coagulation", "platelet", "thrombosis",
    "immune", "interleukin", "cytokine", "inflammation", "complement",
    "vascular", "angiogenesis", "endotheli", "blood_vessel",
    "hypertens", "preeclampsia", "eclampsia",
    "renin", "angiotensin", "aldosterone",
    "tgf_beta", "vegf", "flt1",
    "placenta", "pregnancy", "trophoblast",
    "oxidative_stress", "reactive_oxygen",
    "nitric_oxide", "no_signaling",
    "blood_pressure", "cardiovascular",
]


def is_relevant(pathway_name):
    """Check if pathway name matches biologically relevant keywords."""
    name_lower = pathway_name.lower()
    return any(kw in name_lower for kw in RELEVANT_KEYWORDS)


def main():
    parser = argparse.ArgumentParser(description="Compile case study results")
    parser.add_argument("--work_dir", required=True)
    parser.add_argument("--phenotype", required=True)
    parser.add_argument("--categories", nargs="+", required=True)
    args = parser.parse_args()

    work_dir = args.work_dir
    phenotype = args.phenotype
    categories = args.categories

    magma_models = ["mean", "top"]

    # --- Parse PASCAL results ---
    pascal_file = None
    pascal_dir = os.path.join(work_dir, "pascal_results")
    for fname in os.listdir(pascal_dir) if os.path.isdir(pascal_dir) else []:
        if fname.endswith("--sum.txt") and phenotype.lower() in fname.lower():
            pascal_file = os.path.join(pascal_dir, fname)
            break
    # Fallback: try common pattern
    if not pascal_file:
        pascal_file = os.path.join(pascal_dir,
            f"{phenotype}_pascal.PathwaySet--msigBIOCARTA_KEGG_REACTOME--sum.txt")

    pascal_results = parse_pascal_results(pascal_file) if pascal_file and os.path.exists(pascal_file) else {}
    print(f"[INFO] PASCAL: {len(pascal_results)} pathways loaded")

    # --- Summary stats ---
    summary_path = os.path.join(work_dir, f"{phenotype}_summary_stats.tsv")
    with open(summary_path, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["category", "p_cutoff", "num_loci", "annotated_loci",
                         "significant_hits", "min_qval"])
        for cat in categories:
            result_dir = os.path.join(work_dir, "lsea_results", f"O15_{cat}")
            stats = parse_lsea_stats(result_dir, cat)
            for p_cutoff, s in sorted(stats.items()):
                writer.writerow([cat, p_cutoff, s["num_loci"], s["annotated_loci"],
                                 s["significant_hits"], s["min_qval"]])

    print(f"[INFO] Summary stats written to: {summary_path}")

    # --- Per-category comparison tables ---
    all_relevant = []

    for cat in categories:
        result_dir = os.path.join(work_dir, "lsea_results", f"O15_{cat}")
        if not os.path.isdir(result_dir):
            print(f"[WARN] LSEA results not found for {cat}, skipping")
            continue

        lsea = parse_lsea_results(result_dir, cat)

        # MAGMA results for this category
        magma = {}
        for model in magma_models:
            gsa_file = os.path.join(work_dir, "magma_results",
                                    f"{phenotype}_sets_{model}_{cat}.gsa.out")
            model_results = parse_magma_gsa(gsa_file)
            for pathway, vals in model_results.items():
                if pathway not in magma:
                    magma[pathway] = {}
                magma[pathway][f"magma_{model}_pval"] = vals["magma_pval"]

        # Merge all pathways
        all_pathways = set(lsea.keys()) | set(magma.keys())

        # Write comparison table
        comp_path = os.path.join(work_dir, f"{phenotype}_comparison_{cat}.tsv")
        with open(comp_path, "w") as f:
            writer = csv.writer(f, delimiter="\t")
            header = ["pathway", "lsea_qval", "lsea_pval", "lsea_loci"]
            for model in magma_models:
                header.append(f"magma_{model}_pval")
            if cat == "c2":  # PASCAL only has C2-like (KEGG/REACTOME) sets
                header.extend(["pascal_chi2_pval", "pascal_emp_pval"])
            header.append("biologically_relevant")
            writer.writerow(header)

            rows = []
            for pathway in sorted(all_pathways):
                row = [pathway]
                l = lsea.get(pathway, {})
                row.append(l.get("lsea_qval", "NA"))
                row.append(l.get("lsea_pval", "NA"))
                row.append(l.get("lsea_loci", "NA"))

                m = magma.get(pathway, {})
                for model in magma_models:
                    row.append(m.get(f"magma_{model}_pval", "NA"))

                if cat == "c2":
                    p = pascal_results.get(pathway, {})
                    row.append(p.get("pascal_chi2_pval", "NA"))
                    row.append(p.get("pascal_emp_pval", "NA"))

                relevant = is_relevant(pathway)
                row.append("*" if relevant else "")
                rows.append(row)

                if relevant and l.get("lsea_qval", 1) < 0.05:
                    all_relevant.append({
                        "category": cat,
                        "pathway": pathway,
                        "lsea_qval": l.get("lsea_qval"),
                        "lsea_loci": l.get("lsea_loci"),
                    })

            # Sort by LSEA q-value
            rows.sort(key=lambda r: float(r[1]) if r[1] != "NA" else 999)
            for row in rows:
                writer.writerow(row)

        print(f"[INFO] {cat}: {len(all_pathways)} pathways → {comp_path}")

    # --- Biologically relevant pathways summary ---
    relevant_path = os.path.join(work_dir, f"{phenotype}_relevant_pathways.tsv")
    with open(relevant_path, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["category", "pathway", "lsea_qval", "lsea_loci"])
        for r in sorted(all_relevant, key=lambda x: x["lsea_qval"]):
            writer.writerow([r["category"], r["pathway"], r["lsea_qval"], r["lsea_loci"]])

    print(f"\n[INFO] === Results Summary ===")
    print(f"  Summary stats:     {summary_path}")
    print(f"  Relevant pathways: {relevant_path} ({len(all_relevant)} pathways)")
    for cat in categories:
        comp_path = os.path.join(work_dir, f"{phenotype}_comparison_{cat}.tsv")
        if os.path.exists(comp_path):
            print(f"  Comparison ({cat}): {comp_path}")

    # --- Print top relevant pathways ---
    if all_relevant:
        print(f"\n[INFO] Top biologically relevant pathways (LSEA q < 0.05):")
        for r in sorted(all_relevant, key=lambda x: x["lsea_qval"])[:20]:
            print(f"  {r['category']:6s}  q={r['lsea_qval']:.2e}  loci={r['lsea_loci']}  {r['pathway']}")


if __name__ == "__main__":
    main()
