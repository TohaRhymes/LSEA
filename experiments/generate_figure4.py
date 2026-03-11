#!/usr/bin/env python3
"""Generate Figure 4 (multi-panel GWAS-on-GWAS) and Supplementary Figure 9.

Figure 4:
  Panel A — Continuous heatmap of top phenotypes (by number of significant pairs),
             hierarchically clustered. Tick-label backgrounds are colored by
             trait type (matching the Figure 3 palette).
  Panel B — Square scatter plot of genetic correlation vs LSEA enrichment.

Supplementary Figure 9:
  Full 150x129 continuous heatmap with readable names and trait-type label coloring.

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


# ---- Pastel trait-type palette (consistent with Figure 3 style) ----
TRAIT_TYPE_COLORS = {
    'continuous':    '#C9E4DE',
    'biomarkers':    '#F2C6DE',
    'prescriptions': '#DBCDF0',
    'icd10':         '#F7D9C4',
    'phecode':       '#FAEDCB',
    'categorical':   '#C6DEF1',
}
DEFAULT_TRAIT_COLOR = '#E8E8E8'


def parse_args():
    p = argparse.ArgumentParser(
        description="Generate GWAS-on-GWAS Figure 4 (multi-panel) and Supplementary Figure 9")
    p.add_argument("--results_dir", required=True)
    p.add_argument("--manifest", required=True)
    p.add_argument("--corr_file", required=True)
    p.add_argument("--out_dir", required=True)
    p.add_argument("--top_n", type=int, default=40,
                   help="Number of phenotypes for Panel A subset (default: 40)")
    return p.parse_args()


def make_square_matrix(df):
    all_names = sorted(set(df.index) | set(df.columns))
    result = pd.DataFrame(0.0, index=all_names, columns=all_names)
    for col in df.columns:
        for idx in df.index:
            if col in result.columns and idx in result.index:
                result.at[idx, col] = df.at[idx, col]
    return result


def reorder_by_clustering(df):
    corr = df.corr()
    Z = linkage(corr, method='ward')
    order = leaves_list(Z)
    return df.iloc[order, order]


def load_results(results_dir):
    matching_dirs = glob.glob(os.path.join(results_dir, '*_ukb'))
    print(f"  {len(matching_dirs)} phenotype directories found")
    all_sets = {}
    for cur_dir in tqdm(matching_dirs, desc="Loading results"):
        pvals_files = glob.glob(
            os.path.join(cur_dir, "*_result_7.479176476853146e-09.tsv"))
        try:
            data = pd.read_csv(pvals_files[0], sep='\t').set_index('gene_set')
            key = os.path.basename(cur_dir).replace('_ukb', '')
            all_sets[key] = data['p_value'].apply(
                lambda x: -np.log10(x) if x > 0 else 0)
        except (IndexError, KeyError):
            pass
    print(f"  {len(all_sets)} phenotypes loaded")
    return all_sets


def load_manifest(manifest_path):
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
    """Rename index/columns; also return display→key dict for back-lookup."""
    new_cols, new_idx = [], []
    display2key = {}
    for k in df.columns:
        desc = str(mapping.get(k, k))
        disp = desc[:max_len] + '...' if len(desc) > max_len else desc
        new_cols.append(disp)
        display2key[disp] = k
    for k in df.index:
        desc = str(mapping.get(k, k))
        disp = desc[:max_len] + '...' if len(desc) > max_len else desc
        new_idx.append(disp)
        display2key[disp] = k
    df_out = df.copy()
    df_out.columns = new_cols
    df_out.index = new_idx
    return df_out, display2key


def make_custom_greens():
    greens = sns.color_palette("Greens", as_cmap=True)
    new_colors = greens(np.linspace(0, 1, 256))
    new_colors[:1, :] = np.array([1, 1, 1, 1])
    return LinearSegmentedColormap.from_list("custom_greens", new_colors)


def trait_type_color(trait_type):
    return TRAIT_TYPE_COLORS.get(str(trait_type).lower(), DEFAULT_TRAIT_COLOR)


def color_ticklabels(ax_obj, display2key, key2type, axis='y', bbox_alpha=0.75):
    """Apply colored backgrounds to tick labels based on trait type."""
    get_fn = ax_obj.get_yticklabels if axis == 'y' else ax_obj.get_xticklabels
    for label in get_fn():
        text = label.get_text()
        key = display2key.get(text, text)
        c = trait_type_color(key2type.get(key, ''))
        label.set_bbox({'facecolor': c, 'edgecolor': 'none',
                        'boxstyle': 'round,pad=0.3', 'alpha': bbox_alpha})


def trait_type_legend(ax_obj, key2type, keys_in_plot, fontsize=8.5):
    """Add trait-type legend to axes."""
    present = set(str(key2type.get(k, '')).lower() for k in keys_in_plot)
    elements = [Patch(facecolor=col, edgecolor='#AAAAAA', linewidth=0.5,
                      label=cat.capitalize())
                for cat, col in TRAIT_TYPE_COLORS.items()
                if cat in present]
    if elements:
        ax_obj.legend(handles=elements, loc='upper right',
                      fontsize=fontsize, title='Trait type',
                      title_fontsize=fontsize + 1,
                      frameon=True, fancybox=True,
                      framealpha=0.93, edgecolor='#CCCCCC')


def panel_label(ax_obj, letter, fontsize=20):
    """Place bold panel letter outside top-left corner of axes."""
    ax_obj.text(-0.07, 1.02, letter, transform=ax_obj.transAxes,
                fontsize=fontsize, fontweight='bold', va='bottom', ha='right')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    cmap = make_custom_greens()

    # ---- Load ----
    print("Loading GWAS-on-GWAS results...")
    all_sets = load_results(args.results_dir)
    print("Loading manifest...")
    key2desc, key2type = load_manifest(args.manifest)

    pvals_data = pd.DataFrame(all_sets).fillna(0)
    pvals_data = make_square_matrix(pvals_data)
    full_reordered = reorder_by_clustering(pvals_data)
    threshold = -np.log10(0.05 / (pvals_data.shape[0] ** 2))

    # ==================================================================
    # Supplementary Figure 9: full heatmap, consistent with Figure 4 style
    # ==================================================================
    print("\nGenerating Supplementary Figure 9 (full heatmap)...")
    full_plot = full_reordered.copy()
    np.fill_diagonal(full_plot.values, 0)
    full_keys = list(full_reordered.index)

    flat_full = full_plot.values[full_plot.values > 0]
    vmax_full = np.percentile(flat_full, 95) if len(flat_full) > 0 else None

    supp_named, d2k_full = rename_labels(full_plot, key2desc, max_len=55)

    fig, ax = plt.subplots(figsize=(50, 50))
    sns.heatmap(supp_named, cmap=cmap, square=True,
                vmin=0, vmax=vmax_full,
                xticklabels=True, yticklabels=True, ax=ax,
                linewidths=0.04, linecolor='#F5F5F5',
                cbar_kws={'shrink': 0.22, 'aspect': 30, 'pad': 0.02,
                          'label': '$-\\log_{10}$(p-value)'})
    ax.tick_params(axis='x', labelsize=8, rotation=90)
    ax.tick_params(axis='y', labelsize=8)
    color_ticklabels(ax, d2k_full, key2type, axis='y', bbox_alpha=0.65)
    color_ticklabels(ax, d2k_full, key2type, axis='x', bbox_alpha=0.65)
    trait_type_legend(ax, key2type, full_keys, fontsize=10)

    for fmt in ['pdf', 'png']:
        fig.savefig(os.path.join(args.out_dir, f'SuppFig9_heatmap_full.{fmt}'),
                    bbox_inches='tight', dpi=150)
    plt.close(fig)
    print("  Saved SuppFig9_heatmap_full.pdf/png")

    # ==================================================================
    # Select top phenotypes for Panel A
    # Criterion: phenotypes with the most significant cross-trait pairs
    # (Bonferroni q < 0.05 on pairwise matrix)
    # ==================================================================
    sig_counts = (pvals_data > threshold).sum(axis=1)
    top_phenos = sig_counts.nlargest(args.top_n).index.tolist()
    if len(top_phenos) < 5:
        print(f"  WARNING: falling back to top-{args.top_n} by total score")
        top_phenos = pvals_data.sum(axis=1).nlargest(args.top_n).index.tolist()

    sub_matrix = pvals_data.loc[top_phenos, top_phenos]
    sub_reordered = reorder_by_clustering(sub_matrix)
    print(f"  Selected {len(top_phenos)} phenotypes for Panel A "
          f"(min sig pairs: {int(sig_counts[top_phenos].min())}, "
          f"max: {int(sig_counts[top_phenos].max())})")

    sub_plot = sub_reordered.copy()
    np.fill_diagonal(sub_plot.values, 0)

    flat_sub = sub_plot.values[sub_plot.values > 0]
    vmax_sub = np.percentile(flat_sub, 95) if len(flat_sub) > 0 else None

    sub_named, d2k_sub = rename_labels(sub_plot, key2desc, max_len=45)
    original_keys = list(sub_reordered.index)

    # ==================================================================
    # Figure 4: Panel A (heatmap) + Panel B (scatter)
    # ==================================================================
    print("\nGenerating Figure 4 (multi-panel)...")

    fig = plt.figure(figsize=(24, 13))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.5, 1], wspace=0.35)

    # --- Panel A ---
    ax_heat = fig.add_subplot(gs[0])
    sns.heatmap(sub_named, cmap=cmap, square=True,
                vmin=0, vmax=vmax_sub,
                xticklabels=True, yticklabels=True, ax=ax_heat,
                linewidths=0.15, linecolor='white',
                cbar_kws={'shrink': 0.42, 'aspect': 25, 'pad': 0.02,
                          'label': '$-\\log_{10}$(p-value)'})
    ax_heat.tick_params(axis='x', labelsize=9, rotation=90)
    ax_heat.tick_params(axis='y', labelsize=9)

    # Color label backgrounds
    color_ticklabels(ax_heat, d2k_sub, key2type, axis='y')
    color_ticklabels(ax_heat, d2k_sub, key2type, axis='x')

    # Panel label outside axes
    panel_label(ax_heat, 'A')

    # --- Panel B: square scatter ---
    ax_scat = fig.add_subplot(gs[1])

    gen_corr_data = pd.read_csv(args.corr_file, sep='\t')
    gc_list, lsea_list = [], []
    for _, row in gen_corr_data.iterrows():
        ci, cj = row.i_joined, row.j_joined
        if ci == cj:
            continue
        if ci not in pvals_data.index or cj not in pvals_data.columns:
            continue
        gc_list.append(row.entry)
        lsea_list.append(pvals_data.at[ci, cj])

    y_vals = [np.log10(x) if x > 0 else 0 for x in lsea_list]
    ax_scat.scatter(gc_list, y_vals,
                    alpha=0.35, s=10, color='#D4A574',
                    edgecolors='#A0784C', linewidths=0.2,
                    rasterized=True)
    ax_scat.set_xlabel('Genetic Correlation', fontsize=12)
    ax_scat.set_ylabel('$\\log_{10}(-\\log_{10}$ p-value$)$', fontsize=12)
    ax_scat.tick_params(labelsize=10)
    ax_scat.set_box_aspect(1)   # square axes box
    ax_scat.grid(True, alpha=0.15, linewidth=0.5)
    ax_scat.set_axisbelow(True)
    ax_scat.spines[['top', 'right']].set_visible(False)

    if threshold > 0:
        thr_y = np.log10(threshold)
        ax_scat.axhline(y=thr_y, color='#CC4444',
                        linestyle='--', linewidth=0.8, alpha=0.55,
                        label='Bonferroni threshold')

    panel_label(ax_scat, 'B')

    # --- External trait-type legend (lower-left: empty space below heatmap) ---
    present = set(str(key2type.get(k, '')).lower() for k in original_keys)
    legend_elements = [Patch(facecolor=col, edgecolor='#AAAAAA', linewidth=0.5,
                             label=cat.capitalize())
                       for cat, col in TRAIT_TYPE_COLORS.items()
                       if cat in present]
    if legend_elements:
        fig.legend(handles=legend_elements,
                   loc='lower left',
                   bbox_to_anchor=(0.02, 0.01),
                   ncol=2,
                   fontsize=9,
                   title='Trait type',
                   title_fontsize=10,
                   frameon=True, fancybox=True,
                   framealpha=0.93, edgecolor='#CCCCCC')

    for fmt in ['pdf', 'png']:
        dpi = 150 if fmt == 'pdf' else 300
        fig.savefig(os.path.join(args.out_dir, f'Figure4_gwas_on_gwas.{fmt}'),
                    bbox_inches='tight', dpi=dpi)
    plt.close(fig)
    print(f"  Saved Figure4_gwas_on_gwas.pdf/png ({len(gc_list)} scatter points)")

    print("\nDone!")


if __name__ == '__main__':
    main()
