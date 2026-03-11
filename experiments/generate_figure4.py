#!/usr/bin/env python3
"""Generate Figure 4 (multi-panel GWAS-on-GWAS) and Supplementary Figure 9.

Figure 4:
  Panel A — Continuous heatmap of top phenotypes (by number of significant pairs),
             hierarchically clustered, with trait-type color strips and readable names.
  Panel B — Scatter plot of genetic correlation vs LSEA enrichment significance.

Supplementary Figure 9:
  Full 150x129 continuous heatmap with readable phenotype names from the manifest.

Usage:
    python generate_figure4.py \\
        --results_dir ./validation_NEW_panukb/ \\
        --manifest '../panukb/Pan-UK Biobank phenotype manifest - phenotype_manifest.tsv' \\
        --corr_file ../panukb/pairwise/corr_indep.tsv \\
        --out_dir ./article_figures/
"""

import argparse
import os
import glob
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import linkage, leaves_list
from tqdm import tqdm


# ---- Trait-type color palette ----
TRAIT_TYPE_COLORS = {
    'continuous':    '#4ECDC4',
    'biomarkers':    '#FF6B6B',
    'prescriptions': '#95E1D3',
    'icd10':         '#F38181',
    'phecode':       '#FCE38A',
    'categorical':   '#EAFFD0',
}
DEFAULT_TRAIT_COLOR = '#CCCCCC'


def parse_args():
    p = argparse.ArgumentParser(
        description="Generate GWAS-on-GWAS Figure 4 (multi-panel) and Supplementary Figure 9")
    p.add_argument("--results_dir", required=True,
                   help="Path to GWAS-on-GWAS results (e.g. ./validation_NEW_panukb/)")
    p.add_argument("--manifest", required=True,
                   help="Pan-UKB phenotype manifest TSV")
    p.add_argument("--corr_file", required=True,
                   help="Pairwise genetic correlation TSV (corr_indep.tsv)")
    p.add_argument("--out_dir", required=True,
                   help="Output directory for figures")
    p.add_argument("--top_n", type=int, default=40,
                   help="Number of phenotypes for Panel A subset (default: 40)")
    return p.parse_args()


def make_square_matrix(df):
    """Expand a potentially non-square matrix to a square one (union of row/col names)."""
    all_names = sorted(set(df.index) | set(df.columns))
    result = pd.DataFrame(0.0, index=all_names, columns=all_names)
    for col in df.columns:
        for idx in df.index:
            if col in result.columns and idx in result.index:
                result.at[idx, col] = df.at[idx, col]
    return result


def reorder_by_clustering(df):
    """Reorder rows/columns of a square DataFrame by Ward hierarchical clustering."""
    corr = df.corr()
    Z = linkage(corr, method='ward')
    order = leaves_list(Z)
    return df.iloc[order, order]


def load_results(results_dir):
    """Load GWAS-on-GWAS -log10(p) values from result directories."""
    matching_dirs = glob.glob(os.path.join(results_dir, '*_ukb'))
    print(f"  {len(matching_dirs)} phenotype directories found")

    all_sets = {}
    for cur_dir in tqdm(matching_dirs, desc="Loading results"):
        matching_pvals = glob.glob(
            os.path.join(cur_dir, "*_result_7.479176476853146e-09.tsv"))
        try:
            data = pd.read_csv(matching_pvals[0], sep='\t').set_index('gene_set')
            key = os.path.basename(cur_dir).replace('_ukb', '')
            all_sets[key] = data['p_value'].apply(
                lambda x: -np.log10(x) if x > 0 else 0)
        except (IndexError, KeyError):
            pass

    print(f"  {len(all_sets)} phenotypes loaded")
    return all_sets


def load_manifest(manifest_path):
    """Load Pan-UKB manifest; return key→description and key→trait_type mappings."""
    datam = pd.read_csv(manifest_path, sep='\t')
    datam['idx'] = datam.aws_link.apply(
        lambda x: x.replace(
            'https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/', '')
        .replace('.tsv.bgz', ''))

    key2desc = dict(zip(datam['idx'], datam['description']))
    for k in key2desc:
        if str(key2desc[k]) == 'nan':
            key2desc[k] = k

    key2type = dict(zip(datam['idx'], datam['trait_type']))
    return key2desc, key2type


def rename_labels(df, mapping, max_len=50):
    """Rename index/columns using mapping, truncating long names."""
    new_cols = [str(mapping.get(c, c)) for c in df.columns]
    new_idx = [str(mapping.get(i, i)) for i in df.index]
    new_cols = [n[:max_len] + '...' if len(n) > max_len else n for n in new_cols]
    new_idx = [n[:max_len] + '...' if len(n) > max_len else n for n in new_idx]
    df_out = df.copy()
    df_out.columns = new_cols
    df_out.index = new_idx
    return df_out


def make_custom_greens():
    """Custom green colormap: white → green (first bin forced to white)."""
    greens = sns.color_palette("Greens", as_cmap=True)
    new_colors = greens(np.linspace(0, 1, 256))
    new_colors[:1, :] = np.array([1, 1, 1, 1])
    return LinearSegmentedColormap.from_list("custom_greens", new_colors)


def trait_type_color(trait_type):
    return TRAIT_TYPE_COLORS.get(str(trait_type).lower(), DEFAULT_TRAIT_COLOR)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    cmap = make_custom_greens()

    # ---- Load data ----
    print("Loading GWAS-on-GWAS results...")
    all_sets = load_results(args.results_dir)

    print("Loading manifest...")
    key2desc, key2type = load_manifest(args.manifest)

    # Build full square matrix
    pvals_data = pd.DataFrame(all_sets).fillna(0)
    pvals_data = make_square_matrix(pvals_data)
    full_reordered = reorder_by_clustering(pvals_data)

    # Bonferroni threshold for "significant" pair
    threshold = -np.log10(0.05 / (pvals_data.shape[0] ** 2))

    # ==================================================================
    # Supplementary Figure 9: full continuous heatmap with readable names
    # ==================================================================
    print("\nGenerating Supplementary Figure 9 (full continuous heatmap)...")
    supp_named = rename_labels(full_reordered, key2desc, max_len=55)
    fig, ax = plt.subplots(figsize=(36, 38))
    sns.heatmap(supp_named, cmap=cmap,
                xticklabels=supp_named.columns,
                yticklabels=supp_named.index, ax=ax)
    ax.tick_params(axis='both', labelsize=7)
    plt.tight_layout()
    for fmt in ['pdf', 'png']:
        fig.savefig(os.path.join(args.out_dir, f'SuppFig9_heatmap_full.{fmt}'),
                    bbox_inches='tight', dpi=150)
    plt.close(fig)
    print("  Saved SuppFig9_heatmap_full.pdf/png")

    # ==================================================================
    # Select top phenotypes for Panel A
    # ==================================================================
    sig_counts = (pvals_data > threshold).sum(axis=1)
    top_phenos = sig_counts.nlargest(args.top_n).index.tolist()

    if len(top_phenos) < 5:
        print(f"  WARNING: Only {len(top_phenos)} phenotypes above threshold; "
              f"falling back to top-{args.top_n} by total enrichment score.")
        top_phenos = pvals_data.sum(axis=1).nlargest(args.top_n).index.tolist()

    sub_matrix = pvals_data.loc[top_phenos, top_phenos]
    sub_reordered = reorder_by_clustering(sub_matrix)
    print(f"  Selected {len(top_phenos)} phenotypes for Panel A")

    # ==================================================================
    # Figure 4: multi-panel (A = heatmap, B = scatter)
    # ==================================================================
    print("\nGenerating Figure 4 (multi-panel)...")

    fig = plt.figure(figsize=(26, 18))
    gs = gridspec.GridSpec(1, 2, width_ratios=[2.8, 1], wspace=0.35)

    # --- Panel A: heatmap with category color strips -----------------
    gs_left = gridspec.GridSpecFromSubplotSpec(
        2, 2, subplot_spec=gs[0],
        height_ratios=[0.03, 1], width_ratios=[0.03, 1],
        hspace=0.02, wspace=0.02)

    row_colors = [trait_type_color(key2type.get(k, ''))
                  for k in sub_reordered.index]
    col_colors = [trait_type_color(key2type.get(k, ''))
                  for k in sub_reordered.columns]

    # Top color strip (columns)
    ax_top = fig.add_subplot(gs_left[0, 1])
    for i, c in enumerate(col_colors):
        ax_top.add_patch(plt.Rectangle((i, 0), 1, 1, facecolor=c, edgecolor='none'))
    ax_top.set_xlim(0, len(col_colors))
    ax_top.set_ylim(0, 1)
    ax_top.axis('off')

    # Left color strip (rows)
    ax_left = fig.add_subplot(gs_left[1, 0])
    n_rows = len(row_colors)
    for i, c in enumerate(row_colors):
        ax_left.add_patch(plt.Rectangle((0, n_rows - 1 - i), 1, 1,
                                        facecolor=c, edgecolor='none'))
    ax_left.set_xlim(0, 1)
    ax_left.set_ylim(0, n_rows)
    ax_left.axis('off')

    # Main heatmap
    ax_heat = fig.add_subplot(gs_left[1, 1])
    sub_named = rename_labels(sub_reordered, key2desc, max_len=45)
    sns.heatmap(sub_named, cmap=cmap,
                xticklabels=sub_named.columns,
                yticklabels=sub_named.index, ax=ax_heat,
                cbar_kws={'shrink': 0.6, 'label': '$-\\log_{10}$(p-value)'})
    ax_heat.tick_params(axis='x', labelsize=9, rotation=90)
    ax_heat.tick_params(axis='y', labelsize=9)
    ax_heat.set_title('A', fontsize=22, fontweight='bold', loc='left', pad=12)

    # Trait-type legend (only categories present in the subset)
    present_types = set(str(key2type.get(k, '')).lower()
                        for k in sub_reordered.index)
    legend_elements = [Patch(facecolor=col, label=cat.capitalize())
                       for cat, col in TRAIT_TYPE_COLORS.items()
                       if cat in present_types]
    if legend_elements:
        ax_heat.legend(handles=legend_elements, loc='upper left',
                       bbox_to_anchor=(0, -0.02), ncol=3, fontsize=9,
                       title='Trait type', title_fontsize=10,
                       frameon=True, fancybox=True)

    # --- Panel B: scatter (genetic correlation vs LSEA enrichment) ---
    ax_scat = fig.add_subplot(gs[1])

    gen_corr_data = pd.read_csv(args.corr_file, sep='\t')
    gc_list, lsea_list = [], []
    for _, row in gen_corr_data.iterrows():
        ci, cj = row.i_joined, row.j_joined
        if ci == cj:
            continue
        if ci not in pvals_data.index or cj not in pvals_data.columns:
            continue
        val = pvals_data.at[ci, cj]
        gc_list.append(row.entry)
        lsea_list.append(val)

    # y-axis: log10(-log10(p)) to compress the range
    y_vals = [np.log10(x) if x > 0 else 0 for x in lsea_list]
    ax_scat.scatter(gc_list, y_vals,
                    alpha=0.4, s=20, color='#CD853F',
                    edgecolors='#8B6914', linewidths=0.3)
    ax_scat.set_xlabel('Genetic Correlation', fontsize=13)
    ax_scat.set_ylabel('$\\log_{10}(-\\log_{10}$ p-value$)$', fontsize=13)
    ax_scat.set_title('B', fontsize=22, fontweight='bold', loc='left', pad=12)
    ax_scat.tick_params(labelsize=11)

    # Significance threshold line
    if threshold > 0:
        ax_scat.axhline(y=np.log10(threshold), color='red',
                        linestyle='--', linewidth=0.8, alpha=0.6)

    for fmt in ['pdf', 'png']:
        dpi = 150 if fmt == 'pdf' else 300
        fig.savefig(os.path.join(args.out_dir, f'Figure4_gwas_on_gwas.{fmt}'),
                    bbox_inches='tight', dpi=dpi)
    plt.close(fig)
    print(f"  Saved Figure4_gwas_on_gwas.pdf/png ({len(gc_list)} scatter points)")

    print("\nDone!")


if __name__ == '__main__':
    main()
