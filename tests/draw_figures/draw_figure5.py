#!/usr/bin/env python3
"""Generate Figure 5: O15 cross-trait enrichment (GWAS-on-GWAS) bar chart.

Horizontal bar chart of the top UKB phenotypes enriched in O15_HYPTENSPREG loci,
colored by trait type, with significance threshold indicated.

Usage:
    python generate_figure5.py \
        --gwas_on_gwas results/case_study_o15/O15_HYPTENSPREG_gwas_on_gwas.tsv \
        --out_dir ./article_figures/ \
        --top_n 20
"""

import argparse
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

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

# Trait-type palette (consistent with Figure 4)
TRAIT_TYPE_COLORS = {
    'continuous':    '#C9E4DE',
    'biomarkers':    '#F2C6DE',
    'prescriptions': '#DBCDF0',
    'icd10':         '#F7D9C4',
    'phecode':       '#FAEDCB',
    'categorical':   '#C6DEF1',
}
DEFAULT_COLOR = '#E8E8E8'


def parse_args():
    p = argparse.ArgumentParser(
        description="Generate Figure 5: O15 GWAS-on-GWAS cross-trait enrichment")
    p.add_argument("--gwas_on_gwas", required=True,
                   help="Path to O15_HYPTENSPREG_gwas_on_gwas.tsv")
    p.add_argument("--out_dir", required=True)
    p.add_argument("--top_n", type=int, default=20,
                   help="Number of phenotypes to display (default: 20)")
    p.add_argument("--q_threshold", type=float, default=0.05,
                   help="FDR significance threshold (default: 0.05)")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    plt.rcParams.update(STYLE)

    # Load data
    df = pd.read_csv(args.gwas_on_gwas, sep='\t')
    df = df.sort_values('p_value').head(args.top_n).copy()
    df['neg_log10_p'] = -np.log10(df['p_value'])
    df['significant'] = df['q_value'] < args.q_threshold

    # Extract trait type from phenotype_key
    df['trait_type'] = df['phenotype_key'].apply(lambda x: x.split('-')[0])
    df['color'] = df['trait_type'].apply(
        lambda t: TRAIT_TYPE_COLORS.get(t, DEFAULT_COLOR))

    # Clean up descriptions
    df['label'] = df['description'].fillna(df['phenotype_key'])
    df['label'] = df['label'].apply(
        lambda x: x[:50] + '...' if len(str(x)) > 50 else x)

    # Reverse for bottom-to-top ordering (most significant at top)
    df = df.iloc[::-1]

    # Plot
    fig, ax = plt.subplots(figsize=(10, 0.45 * len(df) + 1.5),
                           constrained_layout=True)

    bars = ax.barh(range(len(df)), df['neg_log10_p'],
                   color=df['color'], edgecolor='#444444', linewidth=0.4)

    # Bold labels for significant results
    ytick_labels = []
    for _, row in df.iterrows():
        if row['significant']:
            ytick_labels.append(f"$\\bf{{{row['label']}}}$")
        else:
            ytick_labels.append(row['label'])

    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df['label'])

    # Bold significant labels via font weight
    for i, (label_obj, (_, row)) in enumerate(
            zip(ax.get_yticklabels(), df.iterrows())):
        if row['significant']:
            label_obj.set_fontweight('bold')

    # Annotate with loci count and q-value
    for i, (_, row) in enumerate(df.iterrows()):
        x_pos = row['neg_log10_p'] + 0.05
        txt = f"q={row['q_value']:.3f}, {int(row['overlapping_loci'])} loci"
        ax.text(x_pos, i, txt, va='center', fontsize=9,
                color='#333333')

    # FDR significance threshold line: draw at p-value of last significant hit
    if df['significant'].any():
        # Boundary = p-value of the least significant result with q < threshold
        boundary_p = df.loc[df['significant'], 'p_value'].max()
        ax.axvline(x=-np.log10(boundary_p), color='#CC4444', linestyle='--',
                   linewidth=1.0, alpha=0.6)
        ax.text(-np.log10(boundary_p) + 0.1, len(df) - 0.5,
                f'FDR = {args.q_threshold}', fontsize=9, color='#CC4444',
                va='top')

    ax.set_xlabel('$-\\log_{10}$(p-value)')
    ax.set_title('Cross-trait enrichment: O15 pregnancy hypertension\n'
                 'vs 129 Pan-UKB phenotypes (GWAS-on-GWAS)',
                 fontweight='bold')
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_axisbelow(True)
    ax.grid(axis='x', alpha=0.15, linewidth=0.5)

    # Trait-type legend
    present = set(df['trait_type'])
    elements = [Patch(facecolor=col, edgecolor='#444444', linewidth=0.5,
                      label=cat.capitalize())
                for cat, col in TRAIT_TYPE_COLORS.items()
                if cat in present]
    if elements:
        ax.legend(handles=elements,
                  loc='lower right',
                  fontsize=11,
                  title='Trait type',
                  title_fontsize=12,
                  frameon=True, fancybox=True,
                  framealpha=0.95, edgecolor='#444444')

    for fmt in ['pdf', 'png']:
        fig.savefig(os.path.join(args.out_dir, f'Figure5_O15_gwas_on_gwas.{fmt}'),
                    bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f"Saved Figure5_O15_gwas_on_gwas.pdf/png ({len(df)} phenotypes)")


if __name__ == '__main__':
    main()
