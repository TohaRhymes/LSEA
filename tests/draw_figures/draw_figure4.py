#!/usr/bin/env python3
"""Generate Figure 4 (multi-panel GWAS-on-GWAS) and Supplementary Figure 5.

Figure 4:
  Panel A — Continuous heatmap of top phenotypes (by number of significant pairs),
             hierarchically clustered. Tick-label backgrounds are colored by
             trait type (matching the Figure 3 palette).
  Panel B — Hexagonal density plot of genetic correlation vs LSEA enrichment
             (FDR-adjusted q-values).
  Panel C — Violin plot of genetic correlations for significant vs
             non-significant enrichment pairs.

Supplementary Figure 5:
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


# ---- Unified plot style (shared across all figure scripts) ----
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
        description="Generate GWAS-on-GWAS Figure 4 (multi-panel) and Supplementary Figure 5")
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
    all_pvals = {}
    all_qvals = {}
    for cur_dir in tqdm(matching_dirs, desc="Loading results"):
        pvals_files = glob.glob(
            os.path.join(cur_dir, "*_result_7.479176476853146e-09.tsv"))
        try:
            data = pd.read_csv(pvals_files[0], sep='\t').set_index('gene_set')
            key = os.path.basename(cur_dir).replace('_ukb', '')
            all_pvals[key] = data['p_value'].apply(
                lambda x: -np.log10(x) if x > 0 else 0)
            all_qvals[key] = data['q_value'].apply(
                lambda x: -np.log10(max(x, 1e-300)) if x >= 0 else 0)
        except (IndexError, KeyError):
            pass
    print(f"  {len(all_pvals)} phenotypes loaded")
    return all_pvals, all_qvals


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
    """Rename index/columns; also return display->key dict for back-lookup."""
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


def add_trait_legend(fig, key2type, keys_in_plot, fontsize=11):
    """Add trait-type legend to figure (lower-left)."""
    present = set(str(key2type.get(k, '')).lower() for k in keys_in_plot)
    elements = [Patch(facecolor=col, edgecolor='#444444', linewidth=0.5,
                      label=cat.capitalize())
                for cat, col in TRAIT_TYPE_COLORS.items()
                if cat in present]
    if elements:
        fig.legend(handles=elements,
                   loc='lower left',
                   bbox_to_anchor=(0.01, 0.01),
                   ncol=1,
                   fontsize=fontsize,
                   title='Trait type',
                   title_fontsize=fontsize + 1,
                   frameon=True, fancybox=True,
                   framealpha=0.95, edgecolor='#444444')


def panel_label(ax_obj, letter, fontsize=20):
    """Place bold panel letter outside top-left corner of axes."""
    ax_obj.text(-0.07, 1.02, letter, transform=ax_obj.transAxes,
                fontsize=fontsize, fontweight='bold', va='bottom', ha='right')


def save_fig(fig, out_dir, basename, dpi=300):
    for fmt in ['pdf', 'png']:
        fig.savefig(os.path.join(out_dir, f'{basename}.{fmt}'),
                    bbox_inches='tight', pad_inches=0.3, dpi=dpi)
    plt.close(fig)
    print(f"  Saved {basename}.pdf/png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    plt.rcParams.update(STYLE)
    cmap = make_custom_greens()

    # ---- Load ----
    print("Loading GWAS-on-GWAS results...")
    all_pvals, all_qvals = load_results(args.results_dir)
    print("Loading manifest...")
    key2desc, key2type = load_manifest(args.manifest)

    pvals_data = pd.DataFrame(all_pvals).fillna(0)
    pvals_data = make_square_matrix(pvals_data)
    qvals_data = pd.DataFrame(all_qvals).fillna(0)
    qvals_data = make_square_matrix(qvals_data)
    full_reordered = reorder_by_clustering(pvals_data)
    # Bonferroni for visualization: 150 queries × 129 targets = 19350 tests
    n_queries = pvals_data.shape[0]
    n_targets = len(all_pvals)  # actual number of loaded phenotypes (columns)
    threshold = -np.log10(0.05 / (n_queries * n_targets))

    # ==================================================================
    # Supplementary Figure 5: full heatmap (style consistent with Fig 4A)
    # ==================================================================
    print("\nGenerating Supplementary Figure 5 (full heatmap)...")
    full_plot = full_reordered.copy()
    np.fill_diagonal(full_plot.values, 0)
    full_keys = list(full_reordered.index)

    flat_full = full_plot.values[full_plot.values > 0]
    vmax_full = np.percentile(flat_full, 95) if len(flat_full) > 0 else None

    supp_named, d2k_full = rename_labels(full_plot, key2desc, max_len=55)

    fig, ax = plt.subplots(figsize=(55, 55))
    sns.heatmap(supp_named, cmap=cmap, square=True,
                vmin=0, vmax=vmax_full,
                xticklabels=True, yticklabels=True, ax=ax,
                linewidths=0.15, linecolor='white',
                cbar_kws={'shrink': 0.35, 'aspect': 30, 'pad': 0.02,
                          'label': '$-\\log_{10}$(p-value)'})
    # Scale fonts for 55x55 figure
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=28)
    cbar.set_label('$-\\log_{10}$(p-value)', fontsize=30)
    ax.tick_params(axis='x', labelsize=26, rotation=90)
    ax.tick_params(axis='y', labelsize=26)
    color_ticklabels(ax, d2k_full, key2type, axis='y', bbox_alpha=0.75)
    color_ticklabels(ax, d2k_full, key2type, axis='x', bbox_alpha=0.75)

    # Legend: vertical, positioned left of the heatmap
    present = set(str(key2type.get(k, '')).lower() for k in full_keys)
    legend_patches = [Patch(facecolor=col, edgecolor='#444444', linewidth=0.5,
                            label=cat.capitalize())
                      for cat, col in TRAIT_TYPE_COLORS.items()
                      if cat in present]
    if legend_patches:
        fig.legend(handles=legend_patches,
                   loc='lower left',
                   bbox_to_anchor=(-0.03, 0.04),
                   ncol=1,
                   fontsize=26,
                   title='Trait type',
                   title_fontsize=28,
                   frameon=True, fancybox=True,
                   framealpha=0.95, edgecolor='#444444')

    save_fig(fig, args.out_dir, 'SuppFig5_heatmap_full', dpi=150)

    # ==================================================================
    # Select top phenotypes for Panel A
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
    # Figure 4: Panel A (heatmap) + Panel B (hexbin) + Panel C (violin)
    # ==================================================================
    print("\nGenerating Figure 4 (multi-panel)...")

    fig = plt.figure(figsize=(24, 14))
    gs_main = gridspec.GridSpec(1, 2, width_ratios=[1.2, 1], wspace=0.25)

    # --- Panel A ---
    ax_heat = fig.add_subplot(gs_main[0])
    sns.heatmap(sub_named, cmap=cmap, square=True,
                vmin=0, vmax=vmax_sub,
                xticklabels=True, yticklabels=True, ax=ax_heat,
                linewidths=0.15, linecolor='white',
                cbar_kws={'shrink': 0.75, 'aspect': 25, 'pad': 0.02,
                          'label': '$-\\log_{10}$(p-value)'})
    ax_heat.set_title('Cross-trait enrichment (top 40 phenotypes)',
                      fontweight='bold', pad=12)
    ax_heat.tick_params(axis='x', labelsize=11, rotation=90)
    ax_heat.tick_params(axis='y', labelsize=11)

    color_ticklabels(ax_heat, d2k_sub, key2type, axis='y')
    color_ticklabels(ax_heat, d2k_sub, key2type, axis='x')
    panel_label(ax_heat, 'A', fontsize=22)

    # --- Right column: Panel B (hexbin) above Panel C (horizontal violin) ---
    from matplotlib.colors import LogNorm
    gs_right = gridspec.GridSpecFromSubplotSpec(
        2, 1, subplot_spec=gs_main[1], height_ratios=[2.8, 1], hspace=0.15)
    ax_hex = fig.add_subplot(gs_right[0])
    ax_viol = fig.add_subplot(gs_right[1], sharex=ax_hex)

    # --- Panel B main: hexbin (genetic correlation vs q-value) ---
    gen_corr_data = pd.read_csv(args.corr_file, sep='\t')
    gc_list, qval_list = [], []
    for _, row in gen_corr_data.iterrows():
        ci, cj = row.i_joined, row.j_joined
        if ci == cj:
            continue
        if ci not in qvals_data.index or cj not in qvals_data.columns:
            continue
        gc_list.append(row.entry)
        qval_list.append(qvals_data.at[ci, cj])

    gc_arr = np.array(gc_list)
    qval_arr = np.array(qval_list)
    # y = log10(-log10(q)) — double-log scale
    y_vals = np.array([np.log10(x) if x > 0 else 0 for x in qval_arr])

    # Pastel-bright hex colormap (white → light lavender → deep rose)
    hex_cmap = LinearSegmentedColormap.from_list('pastel_hex', [
        '#F8F0F8', '#E8C8E0', '#D4A0C8', '#C070A0', '#8C3070', '#4A1040'])
    hb = ax_hex.hexbin(gc_arr, y_vals, gridsize=40, cmap=hex_cmap,
                       mincnt=1, norm=LogNorm(), rasterized=True)
    # Inset colorbar inside the plot area (lower-right, away from data)
    cax = ax_hex.inset_axes([0.85, 0.03, 0.03, 0.35])
    cb = fig.colorbar(hb, cax=cax)
    cb.set_label('Count', fontsize=9)
    cb.ax.tick_params(labelsize=8)
    ax_hex.set_xlabel('Genetic Correlation')
    ax_hex.set_ylabel('$\\log_{10}(-\\log_{10}$ q-value$)$')
    ax_hex.set_title('Genetic correlation vs enrichment',
                      fontweight='bold', pad=10, fontsize=13)
    ax_hex.grid(True, alpha=0.15, linewidth=0.5)
    ax_hex.set_axisbelow(True)
    ax_hex.spines[['top', 'right']].set_visible(False)
    # Zero-correlation reference (visual link to Panel C)
    ax_hex.axvline(x=0, color='#888888', linestyle=':', linewidth=0.6, alpha=0.4)

    # Threshold line at q = 0.05
    thr_q = np.log10(-np.log10(0.05))
    ax_hex.axhline(y=thr_q, color='#CC4444',
                    linestyle='--', linewidth=1.2, alpha=0.7,
                    label='q = 0.05')
    ax_hex.legend(fontsize=10, loc='upper left')
    panel_label(ax_hex, 'B', fontsize=22)

    # --- Panel C: horizontal violin below hex (shared x-axis = genetic correlation) ---
    # Split by enrichment significance: shows GC distribution for sig vs non-sig pairs
    qval_threshold = -np.log10(0.05)  # 1.301 in -log10 space
    sig_mask = qval_arr > qval_threshold
    n_sig = int(sig_mask.sum())
    n_nonsig = len(sig_mask) - n_sig

    sig_label  = f'q < 0.05 (n={n_sig:,})'
    nsig_label = f'q ≥ 0.05 (n={n_nonsig:,})'
    violin_data = pd.DataFrame({
        'gc': gc_arr,
        'Group': np.where(sig_mask, sig_label, nsig_label)
    })
    order = [sig_label, nsig_label]
    sns.violinplot(data=violin_data, x='gc', y='Group',
                   ax=ax_viol, order=order, hue='Group', hue_order=order,
                   palette=['#D45B7A', '#7BA7CC'],
                   inner=None, cut=0, linewidth=0.8, legend=False,
                   orient='h', saturation=0.85, density_norm='width')
    # Overlay real boxplots (thick, visible, clipped to IQR±1.5*IQR)
    for i, grp in enumerate(order):
        subset = violin_data[violin_data['Group'] == grp]['gc']
        bp = ax_viol.boxplot(subset, positions=[i], vert=False, widths=0.22,
                             patch_artist=True, manage_ticks=False,
                             whis=[5, 95],  # 5th-95th percentile (stays within violin)
                             boxprops=dict(facecolor='white', edgecolor='#111111',
                                           linewidth=1.0, alpha=0.9),
                             medianprops=dict(color='#111111', linewidth=2.5),
                             whiskerprops=dict(color='#111111', linewidth=1.0),
                             capprops=dict(color='#111111', linewidth=1.0),
                             flierprops=dict(marker='', markersize=0))
    ax_viol.set_xlabel('Genetic Correlation', fontsize=11)
    ax_viol.set_ylabel('')
    ax_viol.tick_params(axis='y', labelsize=9)
    ax_viol.spines[['top', 'right']].set_visible(False)
    ax_viol.axvline(x=0, color='#888888', linestyle=':', linewidth=0.6, alpha=0.5)
    # Align x-axes: hex keeps its x-label hidden, violin shows it
    ax_hex.tick_params(axis='x', labelbottom=False)
    ax_hex.set_xlabel('')
    # Force identical x-limits
    xlim = ax_hex.get_xlim()
    ax_viol.set_xlim(xlim)
    panel_label(ax_viol, 'C', fontsize=20)

    add_trait_legend(fig, key2type, original_keys, fontsize=11)

    save_fig(fig, args.out_dir, 'Figure4_gwas_on_gwas')

    print(f"\nDone! ({len(gc_list)} points, {n_sig} enriched at q<0.05)")

    # ==================================================================
    # Supplementary heatmaps: trait-type subsets (biomarkers, continuous)
    # ==================================================================
    for trait_filter in ['biomarkers', 'continuous', 'categorical', 'phecode']:
        subset_phenos = [k for k in pvals_data.index
                         if str(key2type.get(k, '')).lower() == trait_filter]
        if len(subset_phenos) < 3:
            print(f"\n  Skipping {trait_filter} heatmap (only {len(subset_phenos)} phenotypes)")
            continue

        print(f"\nGenerating supplementary heatmap: {trait_filter} "
              f"({len(subset_phenos)} phenotypes)...")
        sub_tf = pvals_data.loc[subset_phenos, subset_phenos]
        sub_tf = reorder_by_clustering(sub_tf)
        sub_tf_plot = sub_tf.copy()
        np.fill_diagonal(sub_tf_plot.values, 0)

        flat_tf = sub_tf_plot.values[sub_tf_plot.values > 0]
        vmax_tf = np.percentile(flat_tf, 95) if len(flat_tf) > 0 else None

        max_len = 55 if len(subset_phenos) > 30 else 45
        sub_tf_named, d2k_tf = rename_labels(sub_tf_plot, key2desc, max_len=max_len)

        n = len(subset_phenos)
        figsize = max(10, n * 0.6)
        fig_tf, ax_tf = plt.subplots(figsize=(figsize, figsize))

        # Scale font sizes to matrix dimension
        if n <= 20:
            tick_sz, title_sz, cbar_label_sz, cbar_tick_sz = 12, 16, 14, 11
        elif n <= 40:
            tick_sz, title_sz, cbar_label_sz, cbar_tick_sz = 11, 15, 13, 10
        else:
            tick_sz, title_sz, cbar_label_sz, cbar_tick_sz = 9, 14, 12, 9

        cbar_shrink = min(0.5, 15.0 / figsize)  # shorter bar on big plots

        sns.heatmap(sub_tf_named, cmap=cmap, square=True,
                    vmin=0, vmax=vmax_tf,
                    xticklabels=True, yticklabels=True, ax=ax_tf,
                    linewidths=0.15, linecolor='white',
                    cbar_kws={'shrink': cbar_shrink, 'aspect': 25, 'pad': 0.02,
                              'label': '$-\\log_{10}$(p-value)'})
        cbar = ax_tf.collections[0].colorbar
        cbar.ax.tick_params(labelsize=cbar_tick_sz)
        cbar.set_label('$-\\log_{10}$(p-value)', fontsize=cbar_label_sz)

        ax_tf.set_title(f'Cross-trait enrichment ({trait_filter.capitalize()}, '
                        f'n={n})', fontweight='bold', pad=12,
                        fontsize=title_sz)
        ax_tf.tick_params(axis='x', labelsize=tick_sz, rotation=90)
        ax_tf.tick_params(axis='y', labelsize=tick_sz)

        save_fig(fig_tf, args.out_dir,
                 f'SuppFig_heatmap_{trait_filter}', dpi=200)


if __name__ == '__main__':
    main()
