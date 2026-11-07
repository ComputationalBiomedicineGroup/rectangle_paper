"""Rare cell types in the five real-bulk validation datasets (Figure 2 without Finotello).

A cell type is rare when its median ground-truth fraction across the samples of a dataset is below 5 %. Only cell
types that are scored in Figure 2 are considered, with the same ground truth and label matching
(notebooks/figure2_metrics.py): unmeasured Petitprez values count as 0 there, so they do here, but the scatter
plots draw them as open markers.

Outputs
  output/rare_celltypes_ground_truth.csv    ground-truth range of every scored cell type, with the rare flag
  output/rare_celltypes_metrics.csv         RMSE / correlation per method on the rare cell types
  figures/fig_suppl_rare_celltypes_metrics  RMSE and correlation heatmaps, all methods
  figures/fig_suppl_rare_celltypes_overview RMSE and correlation per dataset, method and rare cell type, in the
                                            layout of fig_suppl1_celltype_method_rmse_correlation
  figures/fig_suppl_rare_celltypes_scatter  estimate vs. ground truth, rare cell types with estimates of every method
  figures/fig_suppl_low_fraction_points     per dataset and method, every (sample, cell type) point whose ground truth
                                            is below the cutoff, coloured by cell type; RMSE / r over those points
  output/low_fraction_points(_metrics)      those points and their per-panel RMSE / r
  figures/fig_suppl_low_fraction_combined   supplementary figure: that scatter grid (A) and the ranked summaries (B-D)

The other methods' per-sample estimates come from the omnideconv runs in data/fig2/real_bulks (unpacked from
real_bulks_27_10_2026.tar.gz, see its README); they reproduce the Figure 2 metrics of every rare cell type exactly.

With --cutoff the median threshold changes (default 0.05); every output carries it as a suffix, e.g. --cutoff 0.1
writes figures/fig_suppl_rare_celltypes_overview_10pct and the default ..._5pct.

Usage: .venv/bin/python scripts/rare_celltypes.py [--cutoff 0.1]
"""
import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm, TwoSlopeNorm
from scipy.stats import pearsonr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'notebooks'))
import figure2_metrics as f2  # noqa: E402

DATASETS = ['petitprez', 'chen', 'hoek', 'morandini', 'dietrich']
RARE_MEDIAN = 0.05
SUFFIX = '_5pct'  # set from the cutoff in main()
METRICS_PATH = ROOT / 'data/fig2/cell_type_level_performance_metrics.csv'
FIGDIR, OUT = ROOT / 'figures', ROOT / 'output'

# Same colours as Figures 2 and 3
METHODS = {
    'rectangle': ('Rectangle', '#000000'), 'secret': ('SECRET', '#A6761D'), 'scdc': ('SCDC', '#B3B3B3'),
    'autogenes': ('AutoGeneS', '#66C2A5'), 'music': ('MuSiC', '#FFD92F'), 'scaden': ('Scaden', '#E5C494'),
    'bayesprism': ('BayesPrism', '#FC8D62'), 'bisque': ('Bisque', '#8DA0CB'), 'dwls': ('DWLS', '#A6D854'),
    'cibersortx': ('CIBERSORTx', '#E78AC3'),
}
DATASET_LABELS = {'petitprez': 'Petitprez', 'chen': 'Chen', 'hoek': 'Hoek', 'morandini': 'Morandini',
                  'dietrich': 'Dietrich'}
# Metric tables and predictions call the same population by different names
CELL_LABELS = {'T cells CD8': 'CD8 T cells', 'T cells CD4 conv': 'CD4 T cells', 'T cells CD4 Total': 'CD4 T cells',
               'T cells CD4': 'CD4 T cells'}
INK, MUTED = '#222222', '#888888'
# Cell-type colours of fig_suppl1_celltype_method_rmse_correlation (tab20 over its sorted cell types)
CELL_COLORS = {'B cells': '#1f77b4', 'CD4 T cells': '#aec7e8', 'CD8 T cells': '#ff7f0e', 'Lymphocytes': '#ffbb78',
               'Macrophages': '#2ca02c', 'Monocytes': '#98df8a', 'NK cells': '#d62728', 'T cells': '#ff9896',
               'Tregs': '#9467bd', 'mDC': '#c5b0d5'}
RECTANGLE_BAND, RECTANGLE_LINE, RECTANGLE_TEXT = '#fff3c4', '#d4a017', '#6b4e00'


def scored_truth():
    """Ground truth of the cell types Figure 2 scores (those Rectangle estimates), per dataset."""
    truth = {}
    for dataset in DATASETS:
        t, p = f2.load_truth_and_predictions(ROOT, dataset, f2.RECTANGLE_PREDICTIONS[dataset])
        truth[dataset] = t.loc[sorted(set(t.index) & set(p.index))]
    return truth


def ground_truth_table(truth):
    rows = []
    for dataset, t in truth.items():
        raw = pd.read_csv(ROOT / f2.TRUTH[dataset][0], index_col=0)
        for cell_type, values in t.iterrows():
            unmeasured = 0
            if f2.TRUTH[dataset][1] is not None:  # labelled by position; NaN = not measured
                unmeasured = int(raw.iloc[f2.TRUTH[dataset][1].index(cell_type)].isna().sum())
            rows.append(dict(dataset=dataset, cell_type=CELL_LABELS.get(cell_type, cell_type), truth_label=cell_type,
                             n_samples=len(values), n_unmeasured=unmeasured, min=values.min(), median=values.median(),
                             mean=values.mean(), max=values.max(), share_below_5pct=(values < 0.05).mean(),
                             rare=values.median() < RARE_MEDIAN))
    return pd.DataFrame(rows)


def rare_metrics(rare):
    """RMSE / correlation of every method on the rare cell types: the Figure 2 table plus SECRET."""
    metrics = pd.read_csv(METRICS_PATH)
    secret, _ = f2.current_method_metrics(ROOT, 'secret', f2.SECRET_PREDICTIONS)
    metrics = pd.concat([metrics, secret], ignore_index=True)
    keys = set(zip(rare['dataset'], rare['cell_type']))
    metrics = metrics[[(d, c) in keys for d, c in zip(metrics['bulk_ds'], metrics['cell_type'])]]
    metrics = metrics.merge(rare[['dataset', 'cell_type', 'mean', 'median']].rename(columns={'dataset': 'bulk_ds'}),
                            on=['bulk_ds', 'cell_type'])
    metrics['rmse_over_mean_abundance'] = metrics['RMSE'] / metrics['mean']
    missing = {m for m in METHODS} - set(metrics['method'])
    if missing:
        raise ValueError(f'No rare-cell metrics for {sorted(missing)}')
    return metrics


def column_label(row):
    return f"{DATASET_LABELS[row.dataset]}\n{row.cell_type}\nmedian {row.median * 100:.1f}%"


def method_order(metrics):
    """All methods, Rectangle included, by their mean RMSE rank on the rare cell types."""
    ranks = metrics.pivot_table(index='method', columns=['bulk_ds', 'cell_type'], values='RMSE').rank().mean(axis=1)
    return ranks.sort_values(kind='stable').index.tolist()


def plot_metrics(metrics, rare):
    columns = list(zip(rare['dataset'], rare['cell_type']))
    labels = [column_label(row) for row in rare.itertuples()]
    order = method_order(metrics)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), layout='constrained', facecolor='white')
    panels = [('RMSE', 'RMSE (fraction)', 'Blues', LogNorm(0.004, 0.5)),
              ('cor', 'Pearson r', 'RdBu_r', TwoSlopeNorm(0, -1, 1))]
    for ax, (value, title, cmap, norm) in zip(axes, panels):
        table = metrics.pivot_table(index='method', columns=['bulk_ds', 'cell_type'], values=value)
        table = table.reindex(index=order, columns=pd.MultiIndex.from_tuples(columns))
        image = ax.imshow(table.to_numpy(dtype=float), cmap=cmap, norm=norm, aspect='auto')
        for i, method in enumerate(order):
            for j, column in enumerate(columns):
                v = table.iat[i, j]
                if not np.isfinite(v):
                    ax.text(j, i, 'n/a', ha='center', va='center', fontsize=8, color=MUTED)
                    continue
                text = f'{v:.3f}' if value == 'RMSE' else f'{v:.2f}'
                rgba = image.cmap(image.norm(v))
                light = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2] > 0.55
                ax.text(j, i, text, ha='center', va='center', fontsize=8.5, color=INK if light else 'white',
                        fontweight='bold' if method == 'rectangle' else 'normal')
        ax.set_xticks(range(len(columns)), labels, fontsize=8.5)
        ax.set_yticks(range(len(order)), [METHODS[m][0] for m in order], fontsize=9.5)
        ax.get_yticklabels()[order.index('rectangle')].set_fontweight('bold')
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_title(title, loc='left', fontsize=11, color=INK)
        colorbar = fig.colorbar(image, ax=ax, shrink=0.8)
        if value == 'RMSE':
            ticks = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]
            colorbar.set_ticks(ticks, labels=[f'{t:g}' for t in ticks])
            colorbar.ax.minorticks_off()
        colorbar.ax.tick_params(labelsize=8)
    save(fig, 'fig_suppl_rare_celltypes_metrics')


def plot_overview(metrics, table):
    """All five datasets x all methods x the rare cell types, RMSE (A) and correlation (B)."""
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FixedLocator, NullLocator

    panels = [('RMSE', 'RMSE (log scale)', 'A'), ('cor', 'Correlation', 'B')]
    any_undefined = False
    fig = plt.figure(figsize=(13, 13), facecolor='white', layout='constrained')
    for subfig, (value, ylabel, letter) in zip(fig.subfigures(2, 1, hspace=0.04), panels):
        subfig.suptitle(letter, x=0.0, ha='left', fontsize=24, fontweight='bold')
        axes = subfig.subplots(2, 3).ravel()
        for ax, dataset in zip(axes, DATASETS):
            ax.set_title(DATASET_LABELS[dataset], fontsize=12, fontweight='bold', color=INK)
            data = metrics[metrics['bulk_ds'] == dataset]
            if data.empty:  # no rare cell type: name the least abundant scored one
                lowest = table[table['dataset'] == dataset].sort_values('median').iloc[0]
                ax.text(0.5, 0.5, f'no rare cell type\n(lowest median: {lowest.cell_type}, {lowest["median"] * 100:.1f}%)',
                        ha='center', va='center', fontsize=10, color=MUTED, transform=ax.transAxes)
                ax.set_xticks([])
                ax.set_yticks([])
                continue
            # All methods, Rectangle included, best first: lowest mean RMSE (A) or highest mean correlation (B)
            # on this dataset's rare cell types
            order = (data.groupby('method')[value].mean()
                     .sort_values(ascending=value == 'RMSE', kind='stable').index.tolist())
            rect_x = order.index('rectangle')
            ax.axvspan(rect_x - 0.5, rect_x + 0.5, color=RECTANGLE_BAND, zorder=0)
            ax.axvline(rect_x, color=RECTANGLE_LINE, linewidth=1.5, zorder=1)
            cell_types = sorted(data['cell_type'].unique())
            for k, (cell_type, group) in enumerate(data.groupby('cell_type')):
                # side by side within a method when a dataset has several rare cell types
                offset = (k - (len(cell_types) - 1) / 2) * 0.18
                x = np.array([order.index(m) + offset for m in group['method']])
                y = group[value].to_numpy(dtype=float)
                defined = np.isfinite(y)
                ax.scatter(x[defined], y[defined], s=60, color=CELL_COLORS[cell_type], edgecolor='white',
                           linewidth=0.8, zorder=3)
                if (~defined).any():  # correlation undefined (constant estimate): an x on the axis floor
                    any_undefined = True
                    ax.scatter(x[~defined], np.full((~defined).sum(), -1.0), s=55, marker='x',
                               color=CELL_COLORS[cell_type], linewidth=1.8, zorder=3, clip_on=False)
            ax.set_xticks(range(len(order)), [METHODS[m][0] for m in order], rotation=90, fontsize=10)
            ax.get_xticklabels()[rect_x].set_fontweight('bold')
            ax.get_xticklabels()[rect_x].set_color(RECTANGLE_TEXT)
            ax.set_xlim(-0.5, len(order) - 0.5)
            if value == 'RMSE':
                ax.set_yscale('log')
                ax.set_ylim(0.003, 0.6)
                ticks = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]
                ax.yaxis.set_major_locator(FixedLocator(ticks))
                ax.yaxis.set_major_formatter(lambda v, _: f'{v:g}')
                ax.yaxis.set_minor_locator(NullLocator())
            else:
                ax.set_ylim(-1.05, 1.05)
                ax.axhline(0, color='#cccccc', linestyle='--', linewidth=0.8, zorder=1)
            ax.grid(axis='y', color='#eeeeee', linewidth=0.6, zorder=0)
            ax.tick_params(axis='y', labelsize=10)
        legend_ax = axes[-1]
        legend_ax.axis('off')
        present = [c for c in CELL_COLORS if c in set(metrics['cell_type'])]
        handles = [Line2D([0], [0], marker='o', linestyle='', markersize=9, markerfacecolor=CELL_COLORS[c],
                          markeredgecolor='white', label=c) for c in present]
        if value == 'cor' and any_undefined:
            handles.append(Line2D([0], [0], marker='x', linestyle='', markersize=9, color=MUTED, markeredgewidth=1.8,
                                  label='r undefined (constant estimate)'))
        legend_ax.legend(handles=handles, loc='center', frameon=False, fontsize=11, title=f'Cell types with median < {RARE_MEDIAN:.0%}',
                         title_fontsize=11)
        subfig.supylabel(ylabel, fontsize=12)
    save(fig, 'fig_suppl_rare_celltypes_overview')


def predictions_for(method, dataset):
    """Per-sample predictions (cell types x samples) or None when they are not available."""
    if method == 'rectangle':
        return f2.load_truth_and_predictions(ROOT, dataset, f2.RECTANGLE_PREDICTIONS[dataset])[1]
    if method == 'secret':
        return f2.load_truth_and_predictions(ROOT, dataset, f2.SECRET_PREDICTIONS[dataset])[1]
    return f2.omnideconv_predictions(ROOT, method, dataset)


def has_all_estimates(row):
    for method in METHODS:
        predictions = predictions_for(method, row.dataset)
        if predictions is None or row.truth_label not in predictions.index:
            return False
    return True


def plot_scatter(truth, rare, metrics):
    """Only rare cell types with per-sample estimates of every method; the metrics figure covers all of them."""
    rare = rare[[has_all_estimates(row) for row in rare.itertuples()]].reset_index(drop=True)
    print('Scatter plots for:', ', '.join(f'{row.dataset} {row.cell_type}' for row in rare.itertuples()))
    # Rectangle first, the others alphabetically; the same columns in every row
    methods = ['rectangle'] + sorted((m for m in METHODS if m != 'rectangle'), key=lambda m: METHODS[m][0].lower())
    any_unmeasured = False
    fig, axes = plt.subplots(len(rare), len(methods), figsize=(1.75 * len(methods), 1.85 * len(rare) + 0.4),
                             layout='constrained', facecolor='white', squeeze=False)
    for i, row in enumerate(rare.itertuples()):
        t = truth[row.dataset].loc[row.truth_label]
        raw = pd.read_csv(ROOT / f2.TRUTH[row.dataset][0], index_col=0)
        unmeasured = pd.Series(False, index=t.index)
        if f2.TRUTH[row.dataset][1] is not None:
            unmeasured = raw.iloc[f2.TRUTH[row.dataset][1].index(row.truth_label)].isna().reindex(t.index)
        any_unmeasured |= bool(unmeasured.any())
        for j, method in enumerate(methods):
            ax = axes[i, j]
            label, color = METHODS[method]
            predictions = predictions_for(method, row.dataset)
            if i == 0:  # method names above the first row, the metrics are the panel titles
                ax.annotate(label, (0.5, 1.22), xycoords='axes fraction', ha='center', va='bottom', fontsize=10,
                            color=INK, fontweight='bold' if method == 'rectangle' else 'normal')
            if j == 0:
                ax.set_ylabel(f'{DATASET_LABELS[row.dataset]}\n{row.cell_type}\nestimate', fontsize=8.5)
            if predictions is None or row.truth_label not in predictions.index:
                ax.text(0.5, 0.5, 'per-sample\nestimates\nnot available', ha='center', va='center', fontsize=7.5,
                        color=MUTED, transform=ax.transAxes)
                ax.set_xticks([])
                ax.set_yticks([])
                for spine in ax.spines.values():
                    spine.set_color('#dddddd')
                continue
            e = predictions.loc[row.truth_label].reindex(t.index).astype(float)
            # x: the ground-truth range, the same for every method of a row; y: each method's own range
            x_limit = t.max() * 1.08
            y_limit = max(t.max(), e.max()) * 1.08
            ax.plot([0, y_limit], [0, y_limit], color=MUTED, linestyle='--', linewidth=0.8, zorder=1)
            measured = ~unmeasured.to_numpy()
            ax.scatter(t[measured], e[measured], s=16, color=color, edgecolor=INK, linewidth=0.4, zorder=2)
            if (~measured).any():
                ax.scatter(t[~measured], e[~measured], s=16, facecolor='white', edgecolor=INK, linewidth=0.6, zorder=2)
            rmse = float(np.sqrt(np.mean((e - t) ** 2)))
            r = float(pearsonr(t, e).statistic) if np.ptp(t) and np.ptp(e) else np.nan
            ax.set_title(f'RMSE {rmse:.3f}  r {r:.2f}', fontsize=7.5, color=INK, pad=3)
            ax.set_xlim(0, x_limit)
            ax.set_ylim(0, y_limit)
            ax.tick_params(labelsize=6.5, length=2)
            ax.xaxis.set_major_formatter(lambda v, _, d=0 if x_limit >= 0.2 else 1: f'{v * 100:.{d}f}%')
            ax.yaxis.set_major_formatter(lambda v, _, d=0 if y_limit >= 0.2 else 1: f'{v * 100:.{d}f}%')
            ax.locator_params(nbins=3)
            ax.spines[['top', 'right']].set_visible(False)
            if i == len(rare) - 1:
                ax.set_xlabel('ground truth', fontsize=8)
    note = 'Dashed: identity.'
    if any_unmeasured:
        note = 'Open markers: not measured in the ground truth, scored as 0.  ' + note
    fig.text(0.005, -0.012, note, fontsize=8, color=MUTED)
    save(fig, 'fig_suppl_rare_celltypes_scatter')


def unmeasured_mask(dataset, truth):
    """True where the ground truth was not measured (scored as 0, as in Figure 2), cell types x samples."""
    mask = pd.DataFrame(False, index=truth.index, columns=truth.columns)
    labels = f2.TRUTH[dataset][1]
    if labels is not None:  # labelled by position; NaN = not measured
        raw = pd.read_csv(ROOT / f2.TRUTH[dataset][0], index_col=0)
        raw.index = labels
        mask = raw.reindex(index=truth.index, columns=truth.columns).isna()
    return mask


def low_fraction_points(truth):
    """All (sample, cell type) pairs with ground truth below the cutoff, with every method's estimate."""
    rows = []
    for dataset, t in truth.items():
        unmeasured = unmeasured_mask(dataset, t)
        for method in METHODS:
            predictions = predictions_for(method, dataset)
            if predictions is None:
                raise FileNotFoundError(f'No {method} estimates for {dataset}')
            for cell_type in t.index:
                if cell_type not in predictions.index:
                    continue
                e = predictions.loc[cell_type].reindex(t.columns).astype(float)
                for sample, value in t.loc[cell_type].items():
                    if value < RARE_MEDIAN:
                        rows.append(dict(dataset=dataset, method=method, sample=sample,
                                         cell_type=CELL_LABELS.get(cell_type, cell_type), truth=value,
                                         estimate=e[sample], unmeasured=bool(unmeasured.at[cell_type, sample])))
    return pd.DataFrame(rows)


def plot_low_fraction_points(points):
    """One panel per dataset and method: every point with ground truth below the cutoff, coloured by cell type;
    RMSE and r over those points."""
    from matplotlib.lines import Line2D

    methods = ['rectangle'] + sorted((m for m in METHODS if m != 'rectangle'), key=lambda m: METHODS[m][0].lower())
    stats = []
    fig, axes = plt.subplots(len(DATASETS), len(methods), figsize=(1.75 * len(methods), 1.85 * len(DATASETS) + 0.8),
                             layout='constrained', facecolor='white', squeeze=False)
    for i, dataset in enumerate(DATASETS):
        for j, method in enumerate(methods):
            ax = axes[i, j]
            if i == 0:
                ax.annotate(METHODS[method][0], (0.5, 1.25), xycoords='axes fraction', ha='center', va='bottom',
                            fontsize=12, color=INK, fontweight='bold' if method == 'rectangle' else 'normal')
            data = points[(points['dataset'] == dataset) & (points['method'] == method)]
            if j == 0:  # n is the same for every method of a dataset, so it goes into the row label
                ax.set_ylabel(f'{DATASET_LABELS[dataset]} (n={len(data)})\nestimate', fontsize=10.5)
            t, e = data['truth'].to_numpy(), data['estimate'].to_numpy()
            rmse = float(np.sqrt(np.mean((e - t) ** 2))) if len(data) else np.nan
            r = float(pearsonr(t, e).statistic) if len(data) > 1 and np.ptp(t) and np.ptp(e) else np.nan
            above = float(np.mean(e > RARE_MEDIAN)) if len(e) else np.nan
            stats.append(dict(dataset=dataset, method=method, n_points=len(data),
                              n_cell_types=data['cell_type'].nunique(), RMSE=rmse, cor=r, share_above_cutoff=above))
            y_limit = max(RARE_MEDIAN, float(e.max()) if len(e) else 0) * 1.08
            ax.plot([0, y_limit], [0, y_limit], color=MUTED, linestyle='--', linewidth=0.8, zorder=1)
            # every point has ground truth below the cutoff, so an estimate above it is too high
            ax.axhline(RARE_MEDIAN, color=MUTED, linestyle='--', linewidth=0.8, zorder=1)  # styled like the identity
            for cell_type, group in data.groupby('cell_type'):
                measured = ~group['unmeasured']
                ax.scatter(group.loc[measured, 'truth'], group.loc[measured, 'estimate'], s=12,
                           color=CELL_COLORS[cell_type], edgecolor=INK, linewidth=0.3, zorder=2)
                if (~measured).any():
                    ax.scatter(group.loc[~measured, 'truth'], group.loc[~measured, 'estimate'], s=12,
                               facecolor='white', edgecolor=CELL_COLORS[cell_type], linewidth=0.8, zorder=2)
            ax.set_title(f'RMSE {rmse:.3f}  r {r:.2f}', fontsize=9, color=INK, pad=3)
            ax.set_xlim(0, RARE_MEDIAN * 1.04)
            ax.set_ylim(0, y_limit)
            ax.tick_params(labelsize=8.5, length=2)
            ax.xaxis.set_major_formatter(lambda v, _: f'{v * 100:.0f}%')
            ax.yaxis.set_major_formatter(lambda v, _, d=0 if y_limit >= 0.2 else 1: f'{v * 100:.{d}f}%')
            ax.locator_params(nbins=3)
            ax.spines[['top', 'right']].set_visible(False)
            if i == len(DATASETS) - 1:
                ax.set_xlabel('ground truth', fontsize=10.5)
    present = [c for c in CELL_COLORS if c in set(points['cell_type'])]
    handles = [Line2D([0], [0], marker='o', linestyle='', markersize=8, markerfacecolor=CELL_COLORS[c],
                      markeredgecolor=INK, markeredgewidth=0.3, label=c) for c in present]
    if points['unmeasured'].any():
        handles.append(Line2D([0], [0], marker='o', linestyle='', markersize=8, markerfacecolor='white',
                              markeredgecolor=MUTED, label='not measured (scored as 0)'))
    handles.append(Line2D([0], [0], color=MUTED, linestyle='--', linewidth=0.8,
                          label=f'identity; {RARE_MEDIAN:.0%} (estimates above are too high)'))
    fig.legend(handles=handles, loc='outside lower center', ncol=6, frameon=False, fontsize=11,
               title=f'Points with ground truth < {RARE_MEDIAN:.0%}', title_fontsize=11)
    save(fig, 'fig_suppl_low_fraction_points')
    return pd.DataFrame(stats)


def low_fraction_celltype_metrics(points, point_stats):
    """RMSE / r per dataset, method and cell type on the low-fraction points, plus the pooled 'all' rows."""
    rows = []
    for (dataset, method, cell_type), g in points.groupby(['dataset', 'method', 'cell_type']):
        t, e = g['truth'].to_numpy(), g['estimate'].to_numpy()
        r = float(pearsonr(t, e).statistic) if len(g) > 1 and np.ptp(t) and np.ptp(e) else np.nan
        rows.append(dict(dataset=dataset, method=method, cell_type=cell_type, n_points=len(g),
                         RMSE=float(np.sqrt(np.mean((e - t) ** 2))), cor=r))
    pooled = point_stats.assign(cell_type='all')[['dataset', 'method', 'cell_type', 'n_points', 'RMSE', 'cor']]
    return pd.concat([pd.DataFrame(rows), pooled], ignore_index=True)


def plot_low_fraction_celltypes(cell_metrics):
    """Layout of fig_suppl1_celltype_method_rmse_correlation on the low-fraction points: per dataset, methods sorted
    by the pooled metric (triangle), dots per cell type."""
    from matplotlib.lines import Line2D

    panels = [('RMSE', 'RMSE', 'A'), ('cor', 'Correlation', 'B')]
    rmse_max = float(cell_metrics['RMSE'].max())
    fig = plt.figure(figsize=(13, 13), facecolor='white', layout='constrained')
    any_undefined = False
    for subfig, (value, ylabel, letter) in zip(fig.subfigures(2, 1, hspace=0.04), panels):
        subfig.suptitle(letter, x=0.0, ha='left', fontsize=24, fontweight='bold')
        axes = subfig.subplots(2, 3).ravel()
        for ax, dataset in zip(axes, DATASETS):
            ax.set_title(DATASET_LABELS[dataset], fontsize=12, fontweight='bold', color=INK)
            data = cell_metrics[cell_metrics['dataset'] == dataset]
            pooled = data[data['cell_type'] == 'all'].set_index('method')[value]
            order = pooled.sort_values(ascending=value == 'RMSE', kind='stable').index.tolist()
            rect_x = order.index('rectangle')
            ax.axvspan(rect_x - 0.5, rect_x + 0.5, color=RECTANGLE_BAND, zorder=0)
            ax.axvline(rect_x, color=RECTANGLE_LINE, linewidth=1.5, zorder=1)
            for cell_type, group in data[data['cell_type'] != 'all'].groupby('cell_type'):
                x = np.array([order.index(m) for m in group['method']])
                y = group[value].to_numpy(dtype=float)
                defined = np.isfinite(y)
                ax.scatter(x[defined], y[defined], s=55, color=CELL_COLORS[cell_type], edgecolor='white',
                           linewidth=0.8, zorder=3)
                if (~defined).any():  # r undefined: an x on the axis floor, as 'absent cell' in suppl1
                    any_undefined = True
                    ax.scatter(x[~defined], np.full((~defined).sum(), -1.0), s=45, marker='x',
                               color=CELL_COLORS[cell_type], linewidth=1.6, zorder=3, clip_on=False)
            ax.scatter([order.index(m) for m in pooled.index], pooled.to_numpy(), s=110, marker='^',
                       facecolor='none', edgecolor=INK, linewidth=1.6, zorder=4)
            ax.set_xticks(range(len(order)), [METHODS[m][0] for m in order], rotation=90, fontsize=10)
            ax.get_xticklabels()[rect_x].set_fontweight('bold')
            ax.get_xticklabels()[rect_x].set_color(RECTANGLE_TEXT)
            ax.set_xlim(-0.5, len(order) - 0.5)
            if value == 'RMSE':
                ax.set_ylim(0, rmse_max * 1.05)
            else:
                ax.set_ylim(-1.05, 1.05)
                ax.axhline(0, color='#cccccc', linestyle='--', linewidth=0.8, zorder=1)
            ax.grid(axis='y', color='#eeeeee', linewidth=0.6, zorder=0)
            ax.tick_params(axis='y', labelsize=10)
        legend_ax = axes[-1]
        legend_ax.axis('off')
        present = [c for c in CELL_COLORS if c in set(cell_metrics['cell_type'])]
        handles = [Line2D([0], [0], marker='o', linestyle='', markersize=9, markerfacecolor=CELL_COLORS[c],
                          markeredgecolor='white', label=c) for c in present]
        handles.append(Line2D([0], [0], marker='^', linestyle='', markersize=10, markerfacecolor='none',
                              markeredgecolor=INK, markeredgewidth=1.6, label='all'))
        if value == 'cor' and any_undefined:
            handles.append(Line2D([0], [0], marker='x', linestyle='', markersize=9, color=MUTED,
                                  markeredgewidth=1.6, label='r undefined'))
        legend_ax.legend(handles=handles, loc='center', frameon=False, fontsize=11, ncol=2,
                         title=f'Points with ground truth < {RARE_MEDIAN:.0%}', title_fontsize=11)
        subfig.supylabel(ylabel, fontsize=12)
    save(fig, 'fig_suppl_low_fraction_celltypes')


def plot_low_fraction_summary(point_stats, letters='ABC', path=None):
    """Methods ranked on the low-fraction points: RMSE stacked over datasets (as Figure 2A), and RMSE and
    correlation per dataset (dots) with their mean over datasets (bar), each ranked by that mean.

    `letters` label the three panels; with `path` the figure is written there (PDF) instead of to figures/."""
    from matplotlib.lines import Line2D

    # Set2 over the six Figure 2 datasets, as in Figure 2A (Finotello takes the first colour)
    colors = {'petitprez': '#fc8d62', 'chen': '#8da0cb', 'hoek': '#e78ac3', 'morandini': '#a6d854',
              'dietrich': '#ffd92f'}
    rmse = point_stats.pivot(index='method', columns='dataset', values='RMSE')[DATASETS]
    cor = point_stats.pivot(index='method', columns='dataset', values='cor')[DATASETS]
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2), layout='constrained', facecolor='white')
    panels = [(axes[0], rmse.sum(axis=1).sort_values(), 'stacked', 'Total RMSE (sum over datasets)', None),
              (axes[1], rmse.mean(axis=1).sort_values(), 'rmse', 'RMSE (dot: dataset, bar: mean)', rmse),
              (axes[2], cor.mean(axis=1).sort_values(ascending=False), 'cor', 'Correlation (dot: dataset, bar: mean)',
               cor)]
    for (ax, order, kind, xlabel, table), letter in zip(panels, letters):
        methods = order.index.tolist()
        y = np.arange(len(methods))
        rect_y = methods.index('rectangle')
        ax.axhspan(rect_y - 0.5, rect_y + 0.5, color=RECTANGLE_BAND, zorder=0)
        if kind == 'stacked':
            left = np.zeros(len(methods))
            for dataset in DATASETS:
                values = rmse.loc[methods, dataset].to_numpy()
                ax.barh(y, values, left=left, height=0.7, color=colors[dataset], edgecolor='white', linewidth=0.8,
                        zorder=2)
                left += values
            for yi, total in zip(y, left):
                ax.text(total, yi, f' {total:.2f}', va='center', fontsize=8.5, color=INK)
            ax.set_xlim(0, left.max() * 1.12)
        else:
            ax.barh(y, order.to_numpy(), height=0.12, color=INK, zorder=2)
            ax.scatter(order.to_numpy(), y, marker='|', s=220, color=INK, linewidth=2, zorder=3)
            for dataset in DATASETS:
                ax.scatter(table.loc[methods, dataset], y, s=45, color=colors[dataset], edgecolor='white',
                           linewidth=0.6, zorder=4)
            if kind == 'cor':
                ax.axvline(0, color='#cccccc', linestyle='--', linewidth=0.8, zorder=1)
                label_x = 1.08
                ax.set_xlim(-1.05, 1.25)
                ax.set_xticks([-1, -0.5, 0, 0.5, 1])
            else:
                label_x = float(np.nanmax(table.to_numpy())) * 1.04
                ax.set_xlim(0, label_x * 1.12)
            for yi, mean in zip(y, order.to_numpy()):
                ax.text(label_x, yi, f'{mean:.3f}' if kind == 'rmse' else f'{mean:.2f}', va='center', ha='left',
                        fontsize=8.5, color=INK)
        ax.set_yticks(y, [METHODS[m][0] for m in methods], fontsize=10)
        ax.get_yticklabels()[rect_y].set_fontweight('bold')
        ax.get_yticklabels()[rect_y].set_color(RECTANGLE_TEXT)
        ax.invert_yaxis()  # best on top
        ax.set_xlabel(xlabel, fontsize=10)
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(axis='x', color='#eeeeee', linewidth=0.6, zorder=0)
        ax.set_title(letter, loc='left', fontsize=18, fontweight='bold', x=-0.18)
    handles = [Line2D([0], [0], marker='s', linestyle='', markersize=10, markerfacecolor=colors[d],
                      markeredgecolor='white', label=DATASET_LABELS[d]) for d in DATASETS]
    fig.legend(handles=handles, loc='outside lower center', ncol=len(DATASETS), frameon=False, fontsize=10,
               title=f'Dataset (points with ground truth < {RARE_MEDIAN:.0%})', title_fontsize=10)
    if path is None:
        save(fig, 'fig_suppl_low_fraction_summary')
    else:
        fig.savefig(path, bbox_inches='tight', facecolor='white')
        plt.close(fig)


def combine_low_fraction_figure(point_stats):
    """Supplementary figure: the low-fraction scatter grid (A) above the three ranked summaries (B-D).

    The saved panel PDFs are embedded as vectors, as in the Figure 2 assembly (figure_2.ipynb)."""
    import tempfile

    import pymupdf

    scatter_path = FIGDIR / f'fig_suppl_low_fraction_points{SUFFIX}.pdf'
    out = FIGDIR / f'fig_suppl_low_fraction_combined{SUFFIX}'
    with tempfile.TemporaryDirectory() as tmp:
        summary_path = Path(tmp) / 'summary.pdf'
        plot_low_fraction_summary(point_stats, letters='BCD', path=summary_path)  # A is the scatter grid
        with pymupdf.open(scatter_path) as scatter, pymupdf.open(summary_path) as summary, pymupdf.open() as combined:
            page_width, margin, gap, label = 1000, 12, 14, 26
            width = page_width - 2 * margin
            heights = [width * doc[0].rect.height / doc[0].rect.width for doc in (scatter, summary)]
            page = combined.new_page(width=page_width, height=margin + label + heights[0] + gap + heights[1] + margin)
            page.draw_rect(page.rect, color=None, fill=(1, 1, 1), overlay=False)
            top = margin + label
            page.show_pdf_page(pymupdf.Rect(margin, top, margin + width, top + heights[0]), scatter, 0)
            page.insert_text((margin, top - 6), 'A', fontsize=24, fontname='hebo')
            bottom = top + heights[0] + gap
            page.show_pdf_page(pymupdf.Rect(margin, bottom, margin + width, bottom + heights[1]), summary, 0)
            combined.save(out.with_suffix('.pdf'), garbage=4, deflate=True)
            page.get_pixmap(dpi=300, alpha=False).save(out.with_suffix('.png'))


def save(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(FIGDIR / f'{name}{SUFFIX}.{ext}', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)


def main():
    global RARE_MEDIAN, SUFFIX
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--cutoff', type=float, default=RARE_MEDIAN, help='median ground-truth fraction below which '
                        'a cell type counts as rare')
    args = parser.parse_args()
    RARE_MEDIAN = args.cutoff
    SUFFIX = f'_{RARE_MEDIAN * 100:g}pct'

    truth = scored_truth()
    table = ground_truth_table(truth)
    table.to_csv(OUT / f'rare_celltypes_ground_truth{SUFFIX}.csv', index=False)
    rare = table[table['rare']].reset_index(drop=True)
    print(table.round(3).to_string(index=False))
    metrics = rare_metrics(rare)
    metrics.to_csv(OUT / f'rare_celltypes_metrics{SUFFIX}.csv', index=False)
    plot_metrics(metrics, rare)
    plot_overview(metrics, table)
    plot_scatter(truth, rare, metrics)
    points = low_fraction_points(truth)
    points.to_csv(OUT / f'low_fraction_points{SUFFIX}.csv', index=False)
    point_stats = plot_low_fraction_points(points)
    point_stats.to_csv(OUT / f'low_fraction_points_metrics{SUFFIX}.csv', index=False)
    plot_low_fraction_summary(point_stats)
    combine_low_fraction_figure(point_stats)
    cell_metrics = low_fraction_celltype_metrics(points, point_stats)
    cell_metrics.to_csv(OUT / f'low_fraction_celltype_metrics{SUFFIX}.csv', index=False)
    plot_low_fraction_celltypes(cell_metrics)
    print(point_stats.pivot(index='method', columns='dataset', values='RMSE')[DATASETS].round(3).to_string())
    print(f'{len(rare)} rare cell types; figures in {FIGDIR}')


if __name__ == '__main__':
    main()
