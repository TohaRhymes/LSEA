#!/usr/bin/env python3
"""
Compile LSEA + MAGMA + PASCAL results for a FinnGen case study phenotype.

Outputs:
1. Per-category comparison tables — all tools shown with BH-adjusted q-values
2. Summary statistics table
3. Highlighted biologically relevant pathways
4. GWAS-on-GWAS cross-trait enrichment table (when --gwas_on_gwas_dir is given)

P-value adjustment:
- LSEA: q-values already BH-corrected internally (per gene set collection)
- MAGMA: raw p-values → BH-FDR applied here, per collection per model
- PASCAL: raw emp p-values → BH-FDR applied here, per collection
All significance calls use q < 0.05.
"""

import argparse
import csv
import math
import os
import sys
from collections import defaultdict

from statsmodels.stats.multitest import multipletests


def parse_lsea_results(result_dir, category, preferred_p="5e-08"):
    """Parse LSEA result TSV file for a given category.

    Prefers the standard GWAS p-cutoff (5e-08) for the comparison table;
    falls back to any available result file.
    """
    results = {}

    # Collect all result files for this category
    candidates = {}
    for fname in os.listdir(result_dir):
        if fname.startswith(f"uni_{category}_result_") and fname.endswith(".tsv"):
            p_str = fname.replace(f"uni_{category}_result_", "").replace(".tsv", "")
            candidates[p_str] = os.path.join(result_dir, fname)

    if not candidates:
        return results

    # Prefer the standard GWAS cutoff, fall back to any available file
    filepath = candidates.get(preferred_p) or next(iter(candidates.values()))

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


def bh_correct(values_dict, pval_key, qval_key):
    """Apply BH-FDR correction to p-values in a dict-of-dicts.

    values_dict: {name: {pval_key: float, ...}, ...}
    Adds qval_key to each inner dict. Skips NaN/missing.
    Returns values_dict (modified in-place).
    """
    names = list(values_dict.keys())
    pvals = [values_dict[n].get(pval_key, float("nan")) for n in names]
    valid_idx = [i for i, p in enumerate(pvals) if not math.isnan(p)]
    if valid_idx:
        valid_p = [pvals[i] for i in valid_idx]
        _, qvals, _, _ = multipletests(valid_p, method="fdr_bh")
        for rank, i in enumerate(valid_idx):
            values_dict[names[i]][qval_key] = qvals[rank]
    return values_dict


def main():
    parser = argparse.ArgumentParser(description="Compile case study results")
    parser.add_argument("--work_dir", required=True)
    parser.add_argument("--phenotype", required=True)
    parser.add_argument("--lsea_prefix", default=None,
                        help="Prefix for LSEA result dirs (default: same as --phenotype). "
                             "E.g. pass 'O15' to find lsea_results/O15_c2/ for phenotype "
                             "O15_HYPTENSPREG.")
    parser.add_argument("--categories", nargs="+", required=True)
    parser.add_argument("--gwas_on_gwas_dir",
                        help="Directory with UKB GWAS-on-GWAS results (O15_ukb/)")
    parser.add_argument("--manifest",
                        help="Pan-UKB phenotype manifest TSV for phenotype descriptions")
    args = parser.parse_args()

    work_dir = args.work_dir
    phenotype = args.phenotype
    lsea_prefix = args.lsea_prefix or phenotype
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
            result_dir = os.path.join(work_dir, "lsea_results", f"{lsea_prefix}_{cat}")
            stats = parse_lsea_stats(result_dir, cat)
            for p_cutoff, s in sorted(stats.items()):
                writer.writerow([cat, p_cutoff, s["num_loci"], s["annotated_loci"],
                                 s["significant_hits"], s["min_qval"]])

    print(f"[INFO] Summary stats written to: {summary_path}")

    # --- Per-category comparison tables ---
    all_relevant = []

    for cat in categories:
        result_dir = os.path.join(work_dir, "lsea_results", f"{lsea_prefix}_{cat}")
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

        # BH-FDR correction for MAGMA (per model, per category)
        for model in magma_models:
            bh_correct(magma, f"magma_{model}_pval", f"magma_{model}_qval")

        # PASCAL results for this category.
        # New mode: per-category Entrez GMT files (e.g. pascal_results/c2_entrez_result.txt)
        # Old mode: single default run (BIOCARTA/KEGG/REACTOME), compared only against c2.
        pascal_cat = {}
        pascal_cat_file = os.path.join(
            work_dir, "pascal_results", f"{cat}_entrez_result.txt"
        )
        if os.path.exists(pascal_cat_file):
            pascal_cat = parse_pascal_results(pascal_cat_file)
        elif cat == "c2":
            # Fallback: old-style monolithic PASCAL run
            pascal_cat = pascal_results

        # BH-FDR correction for PASCAL emp_pval
        if pascal_cat:
            bh_correct(pascal_cat, "pascal_emp_pval", "pascal_emp_qval")

        has_pascal = bool(pascal_cat)

        # Merge all pathways
        all_pathways = set(lsea.keys()) | set(magma.keys()) | set(pascal_cat.keys())

        # Build column header — all tools show both raw p and BH q-value
        header = ["pathway", "lsea_pval", "lsea_qval", "lsea_loci"]
        for model in magma_models:
            header += [f"magma_{model}_pval", f"magma_{model}_qval"]
        if has_pascal:
            header += ["pascal_emp_pval", "pascal_emp_qval", "pascal_chi2_pval"]
        header.append("biologically_relevant")

        # Write comparison table
        comp_path = os.path.join(work_dir, f"{phenotype}_comparison_{cat}.tsv")
        with open(comp_path, "w") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(header)

            rows = []
            for pathway in sorted(all_pathways):
                row = [pathway]
                l = lsea.get(pathway, {})
                row.append(l.get("lsea_pval", "NA"))
                row.append(l.get("lsea_qval", "NA"))
                row.append(l.get("lsea_loci", "NA"))

                m = magma.get(pathway, {})
                for model in magma_models:
                    row.append(m.get(f"magma_{model}_pval", "NA"))
                    row.append(m.get(f"magma_{model}_qval", "NA"))

                if has_pascal:
                    p = pascal_cat.get(pathway, {})
                    row.append(p.get("pascal_emp_pval", "NA"))
                    row.append(p.get("pascal_emp_qval", "NA"))
                    row.append(p.get("pascal_chi2_pval", "NA"))

                relevant = is_relevant(pathway)
                row.append("*" if relevant else "")
                rows.append(row)

                if relevant and l.get("lsea_qval", 1.0) < 0.05:
                    all_relevant.append({
                        "category": cat,
                        "pathway": pathway,
                        "lsea_qval": l.get("lsea_qval"),
                        "lsea_loci": l.get("lsea_loci"),
                    })

            # Sort by LSEA q-value (col index 2)
            rows.sort(key=lambda r: float(r[2]) if r[2] != "NA" else 999.0)
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

    # --- GWAS-on-GWAS cross-trait enrichment ---
    if args.gwas_on_gwas_dir:
        compile_gwas_on_gwas(work_dir, phenotype, args.gwas_on_gwas_dir,
                             args.manifest)


def load_phenotype_manifest(manifest_path):
    """Load Pan-UKB phenotype manifest and build composite key → description map.

    Pan-UKB phenotype keys have the format:
        {trait_type}-{phenocode}-{pheno_sex}[-{coding}][-{modifier}]
    e.g. 'continuous-135-both_sexes', 'biomarkers-30780-both_sexes-irnt',
         'categorical-20003-both_sexes-1140879802'
    """
    desc_map = {}
    if not manifest_path or not os.path.exists(manifest_path):
        return desc_map
    with open(manifest_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            trait_type = row.get("trait_type", "")
            phenocode = row.get("phenocode", "")
            pheno_sex = row.get("pheno_sex", "")
            coding = row.get("coding", "")
            modifier = row.get("modifier", "")
            desc = row.get("description", "")
            coding_desc = row.get("coding_description", "")

            if not phenocode:
                continue

            # Build composite key matching LSEA gene_set names
            parts = [trait_type, phenocode, pheno_sex]
            if coding:
                parts.append(coding)
            if modifier:
                parts.append(modifier)
            key = "-".join(parts)

            # Build human-readable description
            full_desc = desc
            if coding_desc and coding_desc.strip() and coding_desc.strip() != "NA":
                full_desc = f"{desc}: {coding_desc}"
            desc_map[key] = full_desc
    return desc_map


def compile_gwas_on_gwas(work_dir, phenotype, gwas_dir, manifest_path=None):
    """Parse UKB GWAS-on-GWAS results and produce a ranked table."""
    print(f"\n[INFO] === GWAS-on-GWAS (UKB cross-trait) ===")

    # Find result file
    result_file = None
    if os.path.isdir(gwas_dir):
        for fname in os.listdir(gwas_dir):
            if fname.startswith("uni_ukb_result_") and fname.endswith(".tsv"):
                result_file = os.path.join(gwas_dir, fname)
                break

    if not result_file or not os.path.exists(result_file):
        print(f"[WARN] No GWAS-on-GWAS result file found in {gwas_dir}")
        return

    # Load phenotype descriptions from manifest
    desc_map = load_phenotype_manifest(manifest_path)

    # Parse results
    rows = []
    with open(result_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            pheno_key = row["gene_set"]
            pval = float(row["p_value"])
            qval = float(row["q_value"])
            loci = int(row["overlapping_loci"])
            sig = row["significance"]
            description = desc_map.get(pheno_key, pheno_key)
            rows.append({
                "phenotype_key": pheno_key,
                "description": description,
                "p_value": pval,
                "q_value": qval,
                "overlapping_loci": loci,
                "significant": sig,
            })

    # Sort by p-value
    rows.sort(key=lambda r: r["p_value"])

    # Write output
    out_path = os.path.join(work_dir, f"{phenotype}_gwas_on_gwas.tsv")
    with open(out_path, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["phenotype_key", "description", "p_value", "q_value",
                         "overlapping_loci", "significant"])
        for r in rows:
            writer.writerow([r["phenotype_key"], r["description"],
                             f"{r['p_value']:.6e}", f"{r['q_value']:.6e}",
                             r["overlapping_loci"], r["significant"]])

    print(f"[INFO] GWAS-on-GWAS: {len(rows)} phenotypes → {out_path}")

    # Print top results
    sig_count = sum(1 for r in rows if r["significant"].lower() == "true")
    print(f"[INFO] Significant (q < 0.05): {sig_count} / {len(rows)}")
    print(f"\n[INFO] Top 20 UKB phenotypes by p-value:")
    for r in rows[:20]:
        sig_mark = "*" if r["significant"].lower() == "true" else " "
        print(f"  {sig_mark} p={r['p_value']:.2e}  q={r['q_value']:.2e}  "
              f"loci={r['overlapping_loci']}  {r['description']}")


if __name__ == "__main__":
    main()
