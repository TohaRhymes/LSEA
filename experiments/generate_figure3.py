#!/usr/bin/env python3
"""Generate Figure 3 (BCM + GTE enrichment bars) and supplementary enrichment figures.

Outputs:
- Figure 3: BCM + GTE side-by-side bar plots (main article)
- Supplementary figures: C2, GO:CC, GO:MF, GO:BP individual bar plots
- Combined 3x2 grid of all 6 categories (Supplementary Figure 3)
- Distribution histograms of phenotypes per gene set (Supplementary Figure 4)

Usage:
    python generate_figure3.py --results_dir ./validation_NEW_panukb/ --out_dir ./article_figures/
"""

import argparse
import os
import glob
import pandas as pd
from collections import Counter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


def parse_args():
    p = argparse.ArgumentParser(description="Generate Figure 3 and supplementary enrichment figures")
    p.add_argument("--results_dir", required=True,
                   help="Path to validation results directory (e.g. ./validation_NEW_panukb/)")
    p.add_argument("--out_dir", required=True,
                   help="Output directory for figures")
    p.add_argument("--q_threshold", type=float, default=0.05,
                   help="Q-value significance threshold (default: 0.05)")
    p.add_argument("--top_n", type=int, default=30,
                   help="Maximum gene sets to display per category (default: 30)")
    p.add_argument("--min_phenos", type=int, default=2,
                   help="Min phenotypes to include a gene set in bars (default: 2)")
    return p.parse_args()


# Configuration
PVAL_FILE_PATTERN = "*_result_7.479176476853146e-09.tsv"

PATTERNS = ['*_bcm', '*_gte', '*_c2', '*_go_cc', '*_go_mf', '*_go_bp']
NAMES = ['Blood Cell Markers', 'Genotype-Tissue Expression',
         'C2 collection from MSigDB', 'Gene ontology: cellular component',
         'Gene ontology: molecular function', 'Gene ontology: biological process']
SHORT_NAMES = ["BCM", "GTE", "C2", "GO:CC", "GO:MF", "GO:BP"]
COLORS = ["#C9E4DE", "#DBCDF0", "#F7D9C4", "#FAEDCB", "#F2C6DE", "#C6DEF1"]


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # ---- Load data ----
    print("Loading LSEA results...")
    all_sets = {}
    for pattern in PATTERNS:
        matching_dirs = glob.glob(os.path.join(args.results_dir, pattern))
        sets = {}
        for cur_dir in matching_dirs:
            matching_pvals = glob.glob(os.path.join(cur_dir, PVAL_FILE_PATTERN))
            try:
                data = pd.read_csv(matching_pvals[0], sep='\t')
                data = data[data['q_value'] < args.q_threshold]
                sets[os.path.basename(cur_dir)] = list(data.gene_set)
            except (IndexError, KeyError):
                pass
        all_sets[pattern] = sets
        print(f"  {pattern}: {len(matching_dirs)} dirs, {len(sets)} with significant results")

    # Aggregate counts
    final_dict_result = {}
    for pattern, name in zip(PATTERNS, NAMES):
        res = Counter([item for sublist in all_sets[pattern].values()
                       for item in sublist]).most_common()
        final_dict_result[name] = res
        print(f"  {name}: {len(res)} unique gene sets with significant associations")

    # ---- Global plot settings ----
    plt.rcParams.update({
        'font.size': 14,
        'axes.titlesize': 16,
        'axes.labelsize': 14,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'figure.dpi': 300,
    })

    # ---- Figure 3: BCM + GTE (cell types & tissues) ----
    print("\nGenerating Figure 3 (BCM + GTE)...")
    fig, axes = plt.subplots(1, 2, figsize=(18, 8), constrained_layout=True)

    for ax, (cat_name, color, sname) in zip(axes, [
        ('Blood Cell Markers', COLORS[0], 'BCM'),
        ('Genotype-Tissue Expression', COLORS[1], 'GTE'),
    ]):
        df_all = pd.DataFrame(final_dict_result[cat_name], columns=['geneset', 'hits'])
        df = df_all[df_all['hits'] >= args.min_phenos].head(args.top_n).copy()
        df['geneset'] = df['geneset'].str.replace('_', ' ').str.title()
        sns.barplot(y='geneset', x='hits', data=df, color=color,
                    edgecolor='black', ax=ax)
        ax.set_xlabel('Number of associated phenotypes', fontsize=14)
        ax.set_ylabel('')
        ax.set_title(f'{cat_name}\n'
                     f'(≥{args.min_phenos} phenotypes; top {len(df)} shown)',
                     fontsize=16, fontweight='bold')
        ax.xaxis.set_label_position('top')
        ax.xaxis.tick_top()
        ax.tick_params(axis='y', labelsize=13)

    for fmt in ['pdf', 'png']:
        fig.savefig(os.path.join(args.out_dir, f'Figure3_BCM_GTE.{fmt}'),
                    bbox_inches='tight', dpi=300)
    plt.close(fig)
    print("  Saved Figure3_BCM_GTE.pdf/png")

    # ---- Supplementary: individual plots for C2, GO categories ----
    print("\nGenerating supplementary figures (C2, GO)...")
    for cat_name, color, sname in zip(NAMES[2:], COLORS[2:], SHORT_NAMES[2:]):
        fig, ax = plt.subplots(1, 1, figsize=(12, 10), constrained_layout=True)
        df_all = pd.DataFrame(final_dict_result[cat_name], columns=['geneset', 'hits'])
        df = df_all[df_all['hits'] >= args.min_phenos].head(args.top_n).copy()
        df['geneset'] = df['geneset'].str.replace('_', ' ')
        df['geneset'] = df['geneset'].apply(
            lambda x: x[:60] + '...' if len(x) > 60 else x)
        sns.barplot(y='geneset', x='hits', data=df, color=color,
                    edgecolor='black', ax=ax)
        ax.set_xlabel('Number of associated phenotypes', fontsize=14)
        ax.set_ylabel('')
        ax.set_title(f'{cat_name}\n'
                     f'(≥{args.min_phenos} phenotypes; top {len(df)} shown)',
                     fontsize=16, fontweight='bold')
        ax.xaxis.set_label_position('top')
        ax.xaxis.tick_top()
        ax.tick_params(axis='y', labelsize=11)

        fname_base = f'SuppFig_enrichment_{sname.replace(":", "")}'
        for fmt in ['pdf', 'png']:
            fig.savefig(os.path.join(args.out_dir, f'{fname_base}.{fmt}'),
                        bbox_inches='tight', dpi=300)
        plt.close(fig)
        print(f"  Saved {fname_base}.pdf/png")

    # ---- Combined 3x2 grid (Supplementary Figure 3) ----
    print("\nGenerating combined grid figure...")
    fig, axes = plt.subplots(3, 2, figsize=(28, 18), constrained_layout=True)
    axes_flat = axes.T.flatten()  # Column-first order

    for ax, (cat_name, color) in zip(axes_flat, zip(NAMES, COLORS)):
        df_all = pd.DataFrame(final_dict_result[cat_name],
                              columns=['geneset', 'hits'])
        df = df_all[df_all['hits'] >= args.min_phenos].head(args.top_n).copy()
        df['geneset'] = df['geneset'].str.replace('_', ' ')
        df['geneset'] = df['geneset'].apply(
            lambda x: x[:50] + '...' if len(x) > 50 else x)
        sns.barplot(y='geneset', x='hits', data=df, color=color,
                    edgecolor='black', ax=ax)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_title(f'{cat_name}\n'
                     f'(≥{args.min_phenos} phenotypes; top {len(df)} shown)',
                     fontsize=13, fontweight='bold')
        ax.xaxis.set_label_position('top')
        ax.xaxis.tick_top()
        ax.tick_params(axis='y', labelsize=10)

    for fmt in ['pdf', 'png']:
        fig.savefig(os.path.join(args.out_dir, f'SuppFig_enrichment_all_6categories.{fmt}'),
                    bbox_inches='tight', dpi=300)
    plt.close(fig)
    print("  Saved SuppFig_enrichment_all_6categories.pdf/png")

    # ---- Distribution figure (Supplementary Figure 4) ----
    print("\nGenerating distribution figure...")
    fig, axes = plt.subplots(3, 2, figsize=(16, 14), constrained_layout=True)
    axes_flat = axes.T.flatten()

    for ax, (cat_name, color) in zip(axes_flat, zip(NAMES, COLORS)):
        df = pd.DataFrame(final_dict_result[cat_name], columns=['geneset', 'hits'])
        sns.histplot(df.hits, color=color, edgecolor='black', ax=ax)
        ax.set_xlabel('Number of associated phenotypes', fontsize=12)
        ax.set_ylabel('Number of gene sets', fontsize=12)
        ax.set_title(f'{cat_name}', fontsize=13, fontweight='bold')

    for fmt in ['pdf', 'png']:
        fig.savefig(os.path.join(args.out_dir,
                    f'SuppFig_distribution_phenotypes_per_geneset.{fmt}'),
                    bbox_inches='tight', dpi=300)
    plt.close(fig)
    print("  Saved SuppFig_distribution_phenotypes_per_geneset.pdf/png")

    print(f"\nDone! All figures saved to {args.out_dir}")


if __name__ == '__main__':
    main()
