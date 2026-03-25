#!/usr/bin/env python3
"""Generate Supplementary Figures 1 & 2 (sensitivity analysis k-sweep).

Supp Fig 1: 3-panel TPR bar plot (small / medium / big pathway) vs k=1..15.
Supp Fig 2: Scatter of recovered loci from target pathway vs k (all 3 sizes).

Usage:
    python generate_supp_sensitivity.py \
        --results_dir /path/to/validation_NEW_full/ \
        --out_dir ./article_figures/
"""

import argparse
import os
import glob
import re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats


def wilson_ci(p, n, z=1.96):
    """Wilson score 95% CI for a proportion p estimated from n trials."""
    if n == 0:
        return 0.0, 0.0
    center = (p + z**2 / (2 * n)) / (1 + z**2 / n)
    margin = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / (1 + z**2 / n)
    return max(0.0, center - margin), min(1.0, center + margin)


# Target pathway per size
TARGET_PATHWAYS = {
    'small':  'KEGG_STEROID_BIOSYNTHESIS',
    'medium': 'KEGG_PPAR_SIGNALING_PATHWAY',
    'big':    'KEGG_FOCAL_ADHESION',
}

SIZES = ['small', 'medium', 'big']

SIZE_LABELS = {
    'small':  'Small pathway (17 genes)',
    'medium': 'Medium pathway (69 genes)',
    'big':    'Large pathway (199 genes)',
}

# Pastel fill colors (matching Figure 3 palette)
SIZE_COLORS = {
    'small':  '#C9E4DE',
    'medium': '#C6DEF1',
    'big':    '#DBCDF0',
}
SIZE_EDGE = {
    'small':  '#5aaa96',
    'medium': '#5a86c0',
    'big':    '#8a68b8',
}


def parse_args():
    p = argparse.ArgumentParser(
        description="Generate Supp Figs 1 & 2: sensitivity analysis k-sweep")
    p.add_argument("--results_dir", required=True,
                   help="Path to validation_NEW_full/ directory")
    p.add_argument("--out_dir", required=True,
                   help="Output directory for figures")
    return p.parse_args()


def load_sensitivity_data(results_dir):
    """Parse all k-sweep result directories. Returns DataFrame."""
    pattern = re.compile(r'^k10000_path_(\w+?)_k(\d+)_(\d+)_')
    records = []

    for d in glob.glob(os.path.join(results_dir, 'k10000_path_*')):
        m = pattern.match(os.path.basename(d))
        if not m:
            continue
        size, k, it = m.group(1), int(m.group(2)), int(m.group(3))
        if size not in TARGET_PATHWAYS:
            continue

        result_files = glob.glob(os.path.join(d, 'uni_result_*.tsv'))
        if not result_files:
            continue

        try:
            df = pd.read_csv(result_files[0], sep='\t')
            target = TARGET_PATHWAYS[size]
            row = df[df['gene_set'] == target]
            if row.empty:
                records.append({'size': size, 'k': k, 'iter': it,
                                'significant': False, 'loci': 0})
            else:
                row = row.iloc[0]
                if 'significance' in df.columns:
                    sig = str(row['significance']).strip().lower() == 'true'
                else:
                    sig = float(row['q_value']) < 0.05
                loci = int(row['overlapping_loci'])
                records.append({'size': size, 'k': k, 'iter': it,
                                'significant': sig, 'loci': loci})
        except Exception as e:
            print(f"  Warning: could not parse {d}: {e}")

    return pd.DataFrame(records)


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print("Loading sensitivity data...")
    data = load_sensitivity_data(args.results_dir)
    if data.empty:
        print("ERROR: no data loaded — check --results_dir")
        return
    print(f"  Loaded {len(data)} records across "
          f"{data['size'].nunique()} sizes × {data['k'].nunique()} k-values")

    # ---- Unified plot style (shared across all figure scripts) ----
    plt.rcParams.update({
        'font.size':          13,
        'axes.titlesize':     15,
        'axes.labelsize':     13,
        'xtick.labelsize':    11,
        'ytick.labelsize':    11,
        'figure.dpi':         300,
        'savefig.dpi':        300,
        'figure.facecolor':  'white',
        'savefig.facecolor': 'white',
        'axes.facecolor':    'white',
    })

    # ================================================================
    # Supp Fig 1 — TPR vs k, one panel per pathway size
    # ================================================================
    print("\nGenerating Supp Fig 1 (TPR vs k)...")
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5), constrained_layout=True)

    for ax, size in zip(axes, SIZES):
        sub = data[data['size'] == size]
        # Aggregate over 3 iterations using Wilson score 95% CI
        agg = (sub.groupby('k')['significant']
               .agg(['sum', 'count'])
               .reset_index()
               .rename(columns={'sum': 'hits', 'count': 'n'}))
        agg['tpr'] = agg['hits'] / agg['n']
        agg = agg.sort_values('k')

        ks = agg['k'].tolist()
        tpr_vals = agg['tpr'].tolist()
        # Asymmetric Wilson CI bars
        ci_lo, ci_hi = [], []
        for p, n_obs in zip(tpr_vals, agg['n'].tolist()):
            lo, hi = wilson_ci(p, int(n_obs))
            ci_lo.append(p - lo)
            ci_hi.append(hi - p)

        ax.bar(ks, tpr_vals,
               color=SIZE_COLORS[size],
               edgecolor=SIZE_EDGE[size],
               linewidth=0.8,
               width=0.72,
               zorder=2)
        ax.errorbar(ks, tpr_vals,
                    yerr=[ci_lo, ci_hi],
                    fmt='none', color='#333333',
                    linewidth=1.3, capsize=3.5, capthick=1.2,
                    zorder=3)

        ax.set_xlabel('Number of causal SNPs from target pathway (k)',
                      fontsize=12)
        ax.set_ylabel('TPR', fontsize=12)
        ax.set_title(SIZE_LABELS[size], fontsize=14, fontweight='bold')
        ax.set_ylim(0, 1.15)
        ax.set_xlim(0.3, 15.7)
        ax.set_xticks(range(1, 16))
        ax.axhline(1.0, color='#BBBBBB', linewidth=0.7,
                   linestyle='--', zorder=1)
        ax.grid(axis='y', alpha=0.2, linewidth=0.5, zorder=0)
        ax.set_axisbelow(True)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    for fmt in ['pdf', 'png']:
        fig.savefig(os.path.join(args.out_dir, f'SuppFig1_sensitivity_TPR.{fmt}'),
                    bbox_inches='tight', dpi=300)
    plt.close(fig)
    print("  Saved SuppFig1_sensitivity_TPR.pdf/png")

    # ================================================================
    # Supp Fig 2 — Recovered loci vs k, all 3 sizes
    # ================================================================
    print("\nGenerating Supp Fig 2 (Recovered loci vs k)...")
    fig, ax = plt.subplots(figsize=(8, 7), constrained_layout=True)

    # Jitter x slightly so overlapping points are visible
    jitter = {'small': -0.18, 'medium': 0.0, 'big': 0.18}

    for size in SIZES:
        sub = data[data['size'] == size]
        xs = sub['k'].values + jitter[size]
        ys = sub['loci'].values

        ax.scatter(xs, ys,
                   color=SIZE_COLORS[size],
                   edgecolors=SIZE_EDGE[size],
                   linewidths=0.8,
                   s=60, alpha=0.92, zorder=3,
                   label=SIZE_LABELS[size])

        # Regression line (on un-jittered k)
        if len(sub) > 2:
            slope, intercept, r_val, p_val, _ = stats.linregress(
                sub['k'].values, ys)
            xs_line = np.linspace(1, 15, 200)
            ax.plot(xs_line, slope * xs_line + intercept,
                    color=SIZE_EDGE[size], linewidth=2.2, alpha=0.75,
                    zorder=2)

    ax.set_xlabel('Number of causal SNPs from target pathway (k)',
                  fontsize=13)
    ax.set_ylabel('Loci recovered from target pathway',
                  fontsize=13)
    ax.set_title('Sensitivity analysis:\nRecovered loci vs k',
                 fontsize=14, fontweight='bold')
    ax.set_xticks(range(1, 16))
    ax.grid(alpha=0.15, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(fontsize=11, frameon=True, fancybox=True,
              framealpha=0.93, edgecolor='#CCCCCC', loc='upper left')

    for fmt in ['pdf', 'png']:
        fig.savefig(os.path.join(args.out_dir, f'SuppFig2_sensitivity_loci.{fmt}'),
                    bbox_inches='tight', dpi=300)
    plt.close(fig)
    print("  Saved SuppFig2_sensitivity_loci.pdf/png")

    print(f"\nDone! Figures saved to {args.out_dir}")


if __name__ == '__main__':
    main()
