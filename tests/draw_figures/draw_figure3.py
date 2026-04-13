#!/usr/bin/env python3
"""Generate Figure 3 (BCM + GTE enrichment bars) and supplementary enrichment figures.

Outputs:
- Figure 3: BCM + GTE side-by-side bar plots (main article)
- Supplementary Figure 3: Distribution histograms of phenotypes per gene set
- Supplementary Figure 4: Combined 2x2 grid of C2 + GO categories

Usage:
    python generate_figure3.py --results_dir ./validation_NEW_panukb/ --out_dir ./article_figures/
"""

import argparse
import os
import glob
import math
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
    p.add_argument("--pct_threshold", type=float, default=0.05,
                   help="Fraction of phenotypes for widespread enrichment (default: 0.05)")
    p.add_argument("--max_show_all", type=int, default=40,
                   help="If ≤ this many gene sets qualify, show all (default: 40)")
    return p.parse_args()


# Configuration
PVAL_FILE_PATTERN = "*_result_7.479176476853146e-09.tsv"

PATTERNS = ['*_bcm', '*_gte', '*_c2', '*_go_cc', '*_go_mf', '*_go_bp']
NAMES = ['Blood Cell Markers', 'Genotype-Tissue Expression',
         'C2 collection from MSigDB', 'Gene ontology: cellular component',
         'Gene ontology: molecular function', 'Gene ontology: biological process']
SHORT_NAMES = ["BCM", "GTE", "C2", "GO:CC", "GO:MF", "GO:BP"]
COLORS = ["#A5D4CA", "#C4B0E8", "#F0C5A5", "#F3E0A5", "#E8A8C8", "#A0C8E8"]

# Unified plot style (shared across all figure scripts)
STYLE = {
    'font.size': 13,
    'axes.titlesize': 15,
    'axes.labelsize': 13,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'figure.facecolor': 'white',
    'savefig.facecolor': 'white',
    'axes.facecolor': 'white',
}


def smart_select(df_all, n_total_phenotypes, args):
    """Select gene sets to display using adaptive thresholding.

    For small categories (≤ max_show_all qualifying sets): show all with ≥ min_phenos.
    For large categories: show gene sets enriched in ≥ pct_threshold of phenotypes,
    capped at top_n.
    """
    df_qualifying = df_all[df_all['hits'] >= args.min_phenos]

    if len(df_qualifying) <= args.max_show_all:
        return df_qualifying.copy(), args.min_phenos

    widespread_min = max(args.min_phenos,
                         int(math.ceil(n_total_phenotypes * args.pct_threshold)))
    df = df_all[df_all['hits'] >= widespread_min].head(args.top_n).copy()
    return df, widespread_min


def subtitle(n_shown, threshold, total_qualifying=None):
    """Build subtitle describing the selection criterion."""
    if total_qualifying is not None and total_qualifying == n_shown:
        return f'(all {n_shown} with ≥{threshold} phenotypes)'
    return f'(≥{threshold} phenotypes; top {n_shown} shown)'


def save_fig(fig, out_dir, basename):
    for fmt in ['pdf', 'png']:
        fig.savefig(os.path.join(out_dir, f'{basename}.{fmt}'),
                    bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f"  Saved {basename}.pdf/png")


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    plt.rcParams.update(STYLE)

    # ---- Load data ----
    print("Loading LSEA results...")
    all_sets = {}
    n_phenotypes = 0
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
        n_phenotypes = max(n_phenotypes, len(matching_dirs))
        print(f"  {pattern}: {len(matching_dirs)} dirs, {len(sets)} with significant results")

    print(f"  Total phenotypes: {n_phenotypes}")

    # Aggregate counts
    final_dict_result = {}
    for pattern, name in zip(PATTERNS, NAMES):
        res = Counter([item for sublist in all_sets[pattern].values()
                       for item in sublist]).most_common()
        final_dict_result[name] = res
        print(f"  {name}: {len(res)} unique gene sets with significant associations")

    # ---- Figure 3: BCM + GTE (cell types & tissues) ----
    print("\nGenerating Figure 3 (BCM + GTE)...")
    fig, axes = plt.subplots(1, 2, figsize=(18, 8), constrained_layout=True)

    for ax, (cat_name, color) in zip(axes, [
        ('Blood Cell Markers', COLORS[0]),
        ('Genotype-Tissue Expression', COLORS[1]),
    ]):
        df_all = pd.DataFrame(final_dict_result[cat_name], columns=['geneset', 'hits'])
        df, used_thresh = smart_select(df_all, n_phenotypes, args)
        df['geneset'] = df['geneset'].str.replace('_', ' ').str.title()
        total_q = len(df_all[df_all['hits'] >= args.min_phenos])
        sns.barplot(y='geneset', x='hits', data=df, color=color,
                    edgecolor='#444444', ax=ax)
        ax.set_xlabel('Number of associated phenotypes')
        ax.set_ylabel('')
        ax.set_title(f'{cat_name}\n{subtitle(len(df), used_thresh, total_q)}',
                     fontweight='bold')
        ax.xaxis.set_label_position('top')
        ax.xaxis.tick_top()

    save_fig(fig, args.out_dir, 'Figure3_BCM_GTE')

    # ---- Supplementary Figure 3: Distribution histograms ----
    print("\nGenerating Supplementary Figure 3 (distribution)...")
    fig, axes = plt.subplots(3, 2, figsize=(16, 14), constrained_layout=True)
    axes_flat = axes.T.flatten()

    for ax, (cat_name, color) in zip(axes_flat, zip(NAMES, COLORS)):
        df = pd.DataFrame(final_dict_result[cat_name], columns=['geneset', 'hits'])
        sns.histplot(df.hits, color=color, edgecolor='#444444', ax=ax)
        ax.set_xlabel('Number of associated phenotypes')
        ax.set_ylabel('Number of gene sets')
        ax.set_title(f'{cat_name}', fontweight='bold')

    save_fig(fig, args.out_dir, 'SuppFig3_distribution_phenotypes_per_geneset')

    # ---- Supplementary Figure 4: Combined 2x2 grid (C2 + GO) ----
    print("\nGenerating Supplementary Figure 4 (C2 + GO combined)...")
    fig, axes = plt.subplots(2, 2, figsize=(24, 18), constrained_layout=True)
    axes_flat = axes.flatten()

    for ax, (cat_name, color) in zip(axes_flat, zip(NAMES[2:], COLORS[2:])):
        df_all = pd.DataFrame(final_dict_result[cat_name],
                              columns=['geneset', 'hits'])
        df, used_thresh = smart_select(df_all, n_phenotypes, args)
        total_q = len(df_all[df_all['hits'] >= args.min_phenos])
        df['geneset'] = df['geneset'].str.replace('_', ' ')
        df['geneset'] = df['geneset'].apply(
            lambda x: x[:55] + '...' if len(x) > 55 else x)
        sns.barplot(y='geneset', x='hits', data=df, color=color,
                    edgecolor='#444444', ax=ax)
        ax.set_xlabel('Number of associated phenotypes')
        ax.set_ylabel('')
        ax.set_title(f'{cat_name}\n{subtitle(len(df), used_thresh, total_q)}',
                     fontweight='bold')
        ax.xaxis.set_label_position('top')
        ax.xaxis.tick_top()

    save_fig(fig, args.out_dir, 'SuppFig4_enrichment_C2_GO_combined')

    print(f"\nDone! All figures saved to {args.out_dir}")


if __name__ == '__main__':
    main()
