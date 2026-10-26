import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
MPL_CACHE_DIR = Path("/tmp") / "rectangle_paper_matplotlib_cache"
MPL_CACHE_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(Path("/tmp") / "rectangle_paper_xdg_cache"))

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch, Rectangle


DATA_DIR = ROOT / "data" / "fig0"
FIGURE_DIR = ROOT / "figures"
DATA_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

RAW_TABLE_PATH = DATA_DIR / "funky_heatmap_raw_metrics.csv"
SCORE_TABLE_PATH = DATA_DIR / "funky_heatmap_scores.csv"
CATEGORY_SCORE_TABLE_PATH = DATA_DIR / "funky_heatmap_category_scores.csv"
LONG_TABLE_PATH = DATA_DIR / "funky_heatmap_long_metrics.csv"
FIGURE_PATH = FIGURE_DIR / "fig6.pdf"
FIGURE_PNG_PATH = FIGURE_PATH.with_suffix(".png")

METHOD_LABELS = {
    "autogenes": "AutoGeneS",
    "bayesprism": "BayesPrism",
    "bisque": "Bisque",
    "cibersortx": "CIBERSORTx",
    "dwls": "DWLS",
    "music": "MuSiC",
    "musics": "MuSiC",
    "rectangle": "Rectangle",
    "scaden": "Scaden",
    "scdc": "SCDC",
    "secret": "SECRET",
}

METHOD_ORDER = [
    "Rectangle",
    "SCDC",
    "AutoGeneS",
    "MuSiC",
    "Scaden",
    "BayesPrism",
    "Bisque",
    "DWLS",
    "CIBERSORTx",
    "SECRET",
]

DATASET_LABELS = {
    "chen": "Chen",
    "hoek": "Hoek",
    "morandini": "Morandini",
    "petitprez": "Petitprez",
    "dietrich": "Dietrich",
}

# As Figure 2 (figure_2.ipynb), which leaves out Finotello
FIG2_DATASETS = ["petitprez", "chen", "hoek", "morandini", "dietrich"]

METRIC_GROUPS = {
    "Accuracy": [
        *(f"accuracy_rmse_{dataset}" for dataset in FIG2_DATASETS),
        *(f"accuracy_corr_{dataset}" for dataset in FIG2_DATASETS),
    ],
    "Efficiency": ["runtime_hao_20k_s", "memory_hao_20k_mib"],
    "Spillover": ["spillover_pct_hao", "spillover_pct_hoek", "spillover_pct_linsley"],
    "Resolution": ["resolution_rank_wu", "resolution_rank_lambrechts"],
    "Unknown": ["unknown_rmse_all", "unknown_corr_all"],
}

CATEGORY_SCORE_COLUMNS = {
    group: f"score_{group.lower()}" for group in METRIC_GROUPS
}

LOWER_IS_BETTER_PATTERNS = ("rmse", "rank", "runtime", "memory", "spillover")
UNKNOWN_CONTENT_FRACTIONS = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9]


def _canonical_method(value):
    key = str(value).strip().replace(" ", "").replace("-", "").lower()
    return METHOD_LABELS.get(key, str(value).strip())


def _empty_raw_table():
    methods = pd.Index(METHOD_ORDER, name="method")
    columns = [column for group_columns in METRIC_GROUPS.values() for column in group_columns]
    return pd.DataFrame(index=methods, columns=columns, dtype=float)


def _figure2_cell_type_metrics():
    # The per-cell-type scores Figure 2 plots (figure_2.ipynb, suppl. RMSE/correlation figure): the saved table,
    # with the other methods' Morandini scores from results_main_morandini_ncbi and SECRET scored from its saved
    # predictions, as Rectangle
    sys.path.insert(0, str(ROOT / "notebooks"))
    import figure2_metrics

    metrics = pd.read_csv(ROOT / "data" / "fig2" / "cell_type_level_performance_metrics.csv")
    metrics["method"] = metrics["method"].str.lower()
    metrics["bulk_ds"] = metrics["bulk_ds"].str.lower()
    metrics = metrics[~((metrics["bulk_ds"] == "morandini") & (metrics["method"] != "rectangle"))]
    secret_cells, _ = figure2_metrics.current_method_metrics(ROOT, "secret", figure2_metrics.SECRET_PREDICTIONS)
    return pd.concat([metrics, figure2_metrics.morandini_comparator_metrics(ROOT), secret_cells], ignore_index=True)


def _add_accuracy_metrics(raw):
    metrics = _figure2_cell_type_metrics()
    metrics["method"] = metrics["method"].map(_canonical_method)
    metrics["bulk_ds"] = metrics["bulk_ds"].astype(str).str.lower()
    metrics = metrics[
        metrics["method"].isin(raw.index)
        & metrics["bulk_ds"].isin(FIG2_DATASETS)
        & (metrics["cell_type"].astype(str).str.lower() != "all")
    ].copy()

    grouped = metrics.groupby(["method", "bulk_ds"], observed=True)
    rmse_sum = grouped["RMSE"].sum(min_count=1)
    corr_mean = grouped["cor"].mean()

    for method in raw.index:
        for dataset in FIG2_DATASETS:
            key = (method, dataset)
            if key in rmse_sum.index:
                raw.loc[method, f"accuracy_rmse_{dataset}"] = rmse_sum.loc[key]
            if key in corr_mean.index:
                raw.loc[method, f"accuracy_corr_{dataset}"] = corr_mean.loc[key]


def _add_efficiency_metrics(raw):
    # Figure 2B/C values, written by figure_2.ipynb
    benchmark = pd.read_csv(ROOT / "data" / "fig2" / "runtimes" / "benchmark_20k_summary.csv")
    benchmark["method"] = benchmark["Method"].map(_canonical_method)
    benchmark = benchmark.set_index("method")
    for method in raw.index.intersection(benchmark.index):
        raw.loc[method, "runtime_hao_20k_s"] = benchmark.loc[method, "Time_Mean"]
        raw.loc[method, "memory_hao_20k_mib"] = benchmark.loc[method, "Memory_Mean"]


def _add_spillover_metrics(raw):
    # As the Figure 3A plot (figure_3.ipynb): other methods from spillover_total.csv, Rectangle and SECRET as
    # computed in figure_3.ipynb
    spillover = pd.concat([
        pd.read_csv(ROOT / "data" / "fig3" / "spillover_total.csv"),
        pd.read_csv(ROOT / "data" / "fig3" / "spillover_rectangle_secret.csv"),
    ])
    spillover["method"] = spillover["Method"].map(_canonical_method)
    spillover = spillover.drop_duplicates("method", keep="last").set_index("method")
    for method in raw.index.intersection(spillover.index):
        raw.loc[method, "spillover_pct_hao"] = spillover.loc[method, "Hao"] * 100
        raw.loc[method, "spillover_pct_hoek"] = spillover.loc[method, "Hoek"] * 100
        raw.loc[method, "spillover_pct_linsley"] = spillover.loc[method, "Linsley"] * 100


def _resolution_mean_rank(path, methods):
    rmse = pd.read_csv(path, index_col=0)
    rmse.index = rmse.index.map(_canonical_method)
    rmse = rmse[rmse.index.isin(methods)].apply(pd.to_numeric, errors="coerce")
    # Rank the heatmap's methods on each fine cell type (1 = lowest RMSE), then average over the cell types. A method
    # without an estimate for a type is ranked on the others only (MuSiC and CIBERSORTx miss one type each)
    return rmse.rank(axis=0, method="average").mean(axis=1)


def _add_resolution_metrics(raw):
    # Fine-cell-type RMSE ranks, as in the Figure 3 bump charts (figure_3.ipynb); every cell type counts equally
    resolution_sources = {
        "resolution_rank_wu": ROOT / "data" / "fig3" / "wu_fine_method_rmse_all.csv",
        "resolution_rank_lambrechts": ROOT / "data" / "fig3" / "lambrechts_fine_method_rmse_all.csv",
    }
    for column, path in resolution_sources.items():
        values = _resolution_mean_rank(path, raw.index)
        for method in raw.index.intersection(values.index):
            raw.loc[method, column] = values.loc[method]


def _rmse(values):
    true_values = pd.to_numeric(values["true_value"], errors="coerce")
    predicted_values = pd.to_numeric(values["predicted_value"], errors="coerce")
    valid = true_values.notna() & predicted_values.notna()
    if not valid.any():
        return np.nan
    return float(np.sqrt(np.mean((predicted_values[valid] - true_values[valid]) ** 2)))


def _correlation(values):
    true_values = pd.to_numeric(values["true_value"], errors="coerce")
    predicted_values = pd.to_numeric(values["predicted_value"], errors="coerce")
    valid = true_values.notna() & predicted_values.notna()
    if valid.sum() < 2:
        return np.nan
    true_values = true_values[valid]
    predicted_values = predicted_values[valid]
    if np.isclose(true_values.std(ddof=0), 0) or np.isclose(predicted_values.std(ddof=0), 0):
        return np.nan
    return float(np.corrcoef(true_values, predicted_values)[0, 1])


def _unknown_summary_all(unknown, value_column, value_function):
    group_columns = ["method", "replicate", "dataset", "unknown_fraction"]
    replicate_values = (
        unknown.groupby(group_columns, observed=True)
        .apply(value_function, include_groups=False)
        .rename(value_column)
        .reset_index()
    )
    return (
        replicate_values.groupby(["method", "dataset", "unknown_fraction"], observed=True)[value_column]
        .mean()
        .reset_index()
    )


def _add_unknown_metrics(raw):
    # Figure 4C (figure_4.ipynb): the run with the reference split by cell, scored on the four known cell types
    # (Rectangle's and SECRET's 'Unknown' rows left out)
    unknown = pd.read_csv(ROOT / "data" / "fig4" / "unknown_content_cellsplit" / "unknown_content_results.csv")
    unknown = unknown[unknown["unknown_fraction"].isin(UNKNOWN_CONTENT_FRACTIONS) & (unknown["celltype"] != "Unknown")].copy()
    unknown["method"] = unknown["method"].map(_canonical_method)
    unknown = unknown[unknown["method"].isin(raw.index)]

    # Match the Figure 4C "all-cells" panels: compute all-cell metrics within
    # each method/replicate/dataset/unknown-fraction group, average replicates,
    # then average those plotted all-cell summary values for the heatmap.
    rmse_summary_all = _unknown_summary_all(unknown, "rmse", _rmse)
    correlation_summary_all = _unknown_summary_all(unknown, "corr", _correlation)
    unknown_summary = pd.DataFrame(index=raw.index)
    unknown_summary["unknown_rmse_all"] = rmse_summary_all.groupby("method", observed=True)["rmse"].mean()
    unknown_summary["unknown_corr_all"] = correlation_summary_all.groupby("method", observed=True)["corr"].mean()

    for method in raw.index.intersection(unknown_summary.index):
        raw.loc[method, "unknown_rmse_all"] = unknown_summary.loc[method, "unknown_rmse_all"]
        raw.loc[method, "unknown_corr_all"] = unknown_summary.loc[method, "unknown_corr_all"]


def build_raw_table():
    raw = _empty_raw_table()
    _add_accuracy_metrics(raw)
    _add_efficiency_metrics(raw)
    _add_spillover_metrics(raw)
    _add_resolution_metrics(raw)
    _add_unknown_metrics(raw)
    return raw


def normalize_metrics(raw):
    scores = raw.copy()
    for column in scores.columns:
        values = pd.to_numeric(scores[column], errors="coerce")
        valid = values.dropna()
        if valid.empty:
            scores[column] = np.nan
            continue

        min_value = valid.min()
        max_value = valid.max()
        if np.isclose(max_value, min_value):
            scaled = pd.Series(np.nan, index=values.index, dtype=float)
            scaled.loc[valid.index] = 1.0
        elif any(pattern in column for pattern in LOWER_IS_BETTER_PATTERNS):
            scaled = (max_value - values) / (max_value - min_value)
        else:
            scaled = (values - min_value) / (max_value - min_value)
        scores[column] = scaled
    return scores


def compute_category_scores(metric_scores):
    category_scores = pd.DataFrame(index=metric_scores.index)
    for group, columns in METRIC_GROUPS.items():
        category_column = CATEGORY_SCORE_COLUMNS[group]
        # Missing metrics count as zero contribution while preserving equal
        # category weight across methods.
        category_scores[category_column] = metric_scores[columns].fillna(0).sum(axis=1) / len(columns)

    category_scores["overall_score"] = category_scores[list(CATEGORY_SCORE_COLUMNS.values())].sum(axis=1)
    category_scores["overall_score_scaled"] = category_scores["overall_score"] / len(METRIC_GROUPS)
    return category_scores


def build_long_table(raw, metric_scores):
    rows = []
    for group, columns in METRIC_GROUPS.items():
        for column in columns:
            for method in raw.index:
                rows.append(
                    {
                        "method": method,
                        "metric_group": group,
                        "metric": column,
                        "raw_value": raw.loc[method, column],
                        "score": metric_scores.loc[method, column],
                        "direction": "lower_is_better"
                        if any(pattern in column for pattern in LOWER_IS_BETTER_PATTERNS)
                        else "higher_is_better",
                    }
                )
    return pd.DataFrame(rows)


def _sort_by_overall(raw, metric_scores, category_scores):
    sorted_methods = (
        category_scores.sort_values(
            ["overall_score", "overall_score_scaled"],
            ascending=[False, False],
        )
        .index.tolist()
    )
    return (
        raw.reindex(sorted_methods),
        metric_scores.reindex(sorted_methods),
        category_scores.reindex(sorted_methods),
    )


def _metric_display_name(metric):
    replacements = {
        "accuracy_rmse_": "RMSE ",
        "accuracy_corr_": "Correlation ",
        "runtime_hao_20k_s": "Runtime Hao 20k",
        "memory_hao_20k_mib": "Memory Hao 20k",
        "spillover_pct_": "Spillover ",
        "resolution_rank_": "RMSE rank ",  # mean over the fine cell types
        "unknown_rmse_all": "All-cells RMSE",
        "unknown_corr_all": "All-cells correlation",
    }
    label = metric
    for old, new in replacements.items():
        label = label.replace(old, new)
    for key, value in DATASET_LABELS.items():
        label = label.replace(key, value)
    return (
        label.replace("_", " ")
        .title()
        .replace("Rmse", "RMSE")
        .replace("Cumulative Correlation", "Cumulative correlation")
    )


def _rank_scaled(values):
    valid = pd.to_numeric(values, errors="coerce").dropna()
    ranks = pd.Series(np.nan, index=values.index, dtype=float)
    if valid.empty:
        return ranks
    if len(valid) == 1:
        ranks.loc[valid.index] = 1.0
        return ranks
    raw_ranks = valid.rank(method="first", ascending=False)
    ranks.loc[valid.index] = 1.0 - (raw_ranks - 1.0) / (len(valid) - 1.0)
    return ranks


def _group_palette():
    return {
        "Overall": "#53B8C5",
        "Accuracy": "#8D98C8",
        "Efficiency": "#94B6CE",
        "Spillover": "#EBA0A2",
        "Resolution": "#DD6AA8",
        "Unknown": "#E4A247",
    }


def _score_cmaps():
    return {
        "Accuracy": LinearSegmentedColormap.from_list("fig0_accuracy", ["#f5f3fb", "#b5bce0", "#5964ad"]),
        "Efficiency": LinearSegmentedColormap.from_list("fig0_efficiency", ["#eff7f8", "#b4d5df", "#3f8faa"]),
        "Spillover": LinearSegmentedColormap.from_list("fig0_spillover", ["#fff1ec", "#f1aaa3", "#c74755"]),
        "Resolution": LinearSegmentedColormap.from_list("fig0_resolution", ["#fff0f6", "#e894bf", "#b22e80"]),
        "Unknown": LinearSegmentedColormap.from_list("fig0_unknown", ["#fff4df", "#e7b66f", "#aa601a"]),
        "Overall": LinearSegmentedColormap.from_list("fig0_overall", ["#eef8f9", "#8dd3da", "#198998"]),
    }


def _draw_header(ax, left, right, label, color, fontsize=9):
    gap = 0.08
    width = max(right - left - 2 * gap, 0.35)
    patch = FancyBboxPatch(
        (left + gap, -1.35),
        width,
        0.56,
        boxstyle="round,pad=0.02,rounding_size=0.16",
        linewidth=0,
        facecolor=color,
        clip_on=False,
        zorder=6,
    )
    ax.add_patch(patch)
    ax.text(
        (left + right) / 2,
        -1.07,
        label,
        ha="center",
        va="center",
        fontsize=fontsize,
        color="white",
        zorder=7,
    )


def _draw_score_legend(ax, score_cmaps):
    x0 = 1.03
    y0 = 0.72
    scores = [0.0, 0.25, 0.5, 0.75, 1.0]
    rank_count = len(METHOD_ORDER)
    groups = [*METRIC_GROUPS.keys(), "Overall"]
    group_abbreviations = {
        "Accuracy": "Acc",
        "Efficiency": "Eff",
        "Spillover": "Spi",
        "Resolution": "Res",
        "Unknown": "Unk",
        "Overall": "All",
    }

    ax.text(
        x0,
        y0 + rank_count * 0.014 + 0.028,
        "Rank color",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
        color="#111111",
        clip_on=False,
    )
    rect_width = 0.017
    rect_height = 0.014
    column_gap = 0.0015
    for group_idx, group in enumerate(groups):
        x = x0 + group_idx * (rect_width + column_gap)
        for rank_idx in range(rank_count):
            rank_value = rank_idx / (rank_count - 1)
            ax.add_patch(
                Rectangle(
                    (x, y0 + rank_idx * rect_height),
                    rect_width,
                    rect_height,
                    transform=ax.transAxes,
                    facecolor=score_cmaps[group](rank_value),
                    edgecolor="#222222",
                    linewidth=0.25,
                    clip_on=False,
                    zorder=8,
                )
            )
        ax.text(
            x + rect_width / 2,
            y0 - 0.014,
            group_abbreviations[group],
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=5.8,
            color="#111111",
            clip_on=False,
        )
    arrow_x = x0 + len(groups) * (rect_width + column_gap) + 0.01
    ax.annotate(
        "",
        xy=(arrow_x, y0 + rank_count * rect_height),
        xytext=(arrow_x, y0),
        xycoords=ax.transAxes,
        textcoords=ax.transAxes,
        arrowprops={"arrowstyle": "->", "linewidth": 0.8, "color": "#111111"},
        annotation_clip=False,
    )
    ax.text(
        arrow_x + 0.008,
        y0 + rank_count * rect_height,
        "Rank 1",
        transform=ax.transAxes,
        ha="left",
        va="center",
        fontsize=7.2,
        clip_on=False,
    )
    ax.text(arrow_x + 0.008, y0, f"Rank {rank_count}", transform=ax.transAxes, ha="left", va="center", fontsize=7.2, clip_on=False)

    size_y = 0.51
    ax.text(
        x0,
        size_y + 0.10,
        "Score size",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
        color="#111111",
        clip_on=False,
    )
    for score_idx, score in enumerate(scores):
        x = x0 + 0.018 + score_idx * 0.022
        ax.scatter(
            x,
            size_y + 0.035,
            s=8 + score * 265,
            facecolors="white",
            edgecolors="#777777",
            linewidths=0.75,
            transform=ax.transAxes,
            clip_on=False,
            zorder=8,
        )
        ax.text(
            x,
            size_y - 0.02,
            f"{score:g}",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=6.8,
            clip_on=False,
        )

    na_y = size_y - 0.09
    ax.text(x0, na_y, "NA", transform=ax.transAxes, ha="left", va="center", fontsize=7.2, color="#111111", clip_on=False)
    ax.scatter(
        x0 + 0.035,
        na_y,
        marker="x",
        s=24,
        color="#b8b8b8",
        linewidths=0.9,
        transform=ax.transAxes,
        clip_on=False,
        zorder=8,
    )


def plot_funky_heatmap(raw, metric_scores, category_scores, save_path=FIGURE_PATH):
    methods = raw.index.tolist()
    column_specs = []
    for group, columns in METRIC_GROUPS.items():
        for column in columns:
            column_specs.append({"kind": "metric", "group": group, "column": column})
        column_specs.append({"kind": "category", "group": group, "column": CATEGORY_SCORE_COLUMNS[group]})
    column_specs.append({"kind": "overall", "group": "Overall", "column": "overall_score_scaled"})
    column_widths = [1.0 if spec["kind"] != "overall" else 3.0 for spec in column_specs]
    column_lefts = np.r_[0, np.cumsum(column_widths[:-1])]
    column_centers = column_lefts + np.array(column_widths) / 2
    total_width = float(np.sum(column_widths))

    fig_width = max(14.4, 0.37 * len(column_specs) + 1.9)
    fig_height = max(5.7, 0.43 * len(methods) + 2.1)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.set_facecolor("white")
    fig.patch.set_facecolor("white")

    score_cmaps = _score_cmaps()
    group_palette = _group_palette()
    method_column_left = -2.75
    row_right = total_width + 0.8
    metric_ranks = metric_scores.apply(_rank_scaled)
    category_rank_columns = [*CATEGORY_SCORE_COLUMNS.values(), "overall_score_scaled"]
    category_ranks = category_scores[category_rank_columns].apply(_rank_scaled)

    y_positions = np.arange(len(methods))
    category_overall_indices = [idx for idx, spec in enumerate(column_specs) if spec["kind"] == "category"]

    for y_idx, method in enumerate(methods):
        row_color = "#f2f2f2" if y_idx % 2 == 0 else "#ffffff"
        ax.add_patch(
            Rectangle(
                (method_column_left, y_idx - 0.5),
                row_right - method_column_left,
                1.0,
                facecolor=row_color,
                edgecolor="none",
                zorder=0,
            )
        )
        ax.text(
            method_column_left + 0.15,
            y_idx,
            method,
            ha="left",
            va="center",
            fontsize=10,
            fontweight="bold" if method == "Rectangle" else "normal",
            color="#111111",
            zorder=6,
        )

    for category_idx in category_overall_indices:
        group_color = group_palette[column_specs[category_idx]["group"]]
        overall_left = column_centers[category_idx] - 0.44
        overall_right = column_centers[category_idx] + 0.44
        ax.add_patch(
            Rectangle(
                (overall_left, -0.5),
                overall_right - overall_left,
                len(methods),
                facecolor=group_color,
                edgecolor="none",
                alpha=0.16,
                zorder=0.5,
            )
        )
        ax.plot([overall_left, overall_left], [-0.5, len(methods) - 0.5], color=group_color, linewidth=0.7, alpha=0.45, zorder=2)
        ax.plot([overall_right, overall_right], [-0.5, len(methods) - 0.5], color=group_color, linewidth=0.7, alpha=0.45, zorder=2)

    for y_idx, method in enumerate(methods):
        for x_idx, spec in enumerate(column_specs):
            x_center = column_centers[x_idx]
            kind = spec["kind"]
            column = spec["column"]
            cmap = score_cmaps[spec["group"]]
            if kind == "metric":
                score = metric_scores.loc[method, column]
                if pd.isna(score):
                    ax.scatter(x_center, y_idx, marker="x", s=24, color="#b8b8b8", linewidths=0.9, zorder=5)
                    continue
                rank = metric_ranks.loc[method, column]
                size = 8 + float(score) * 265
                ax.scatter(
                    x_center,
                    y_idx,
                    s=size,
                    c=[cmap(float(rank))],
                    edgecolors="#202020",
                    linewidths=0.45,
                    zorder=5,
                )
            else:
                score = category_scores.loc[method, column]
                if pd.isna(score):
                    ax.scatter(x_center, y_idx, marker="x", s=30, color="#b8b8b8", linewidths=1.0, zorder=5)
                    continue
                if kind == "overall":
                    cmap = score_cmaps["Overall"]
                    rank = category_ranks.loc[method, "overall_score_scaled"]
                    bar_left = column_lefts[x_idx] + 0.28
                    bar_max_width = column_widths[x_idx] - 0.78
                    bar_width = bar_max_width * float(score)
                    ax.add_patch(
                        Rectangle(
                            (bar_left, y_idx - 0.33),
                            bar_width,
                            0.66,
                            facecolor=cmap(float(rank)),
                            edgecolor="#222222",
                            linewidth=0.55,
                            zorder=4,
                        )
                    )
                    ax.text(
                        column_lefts[x_idx] + column_widths[x_idx] + 0.10,
                        y_idx,
                        f"{category_scores.loc[method, 'overall_score']:.1f}",
                        ha="right",
                        va="center",
                        fontsize=8.0,
                        fontweight="bold",
                        color="#111111",
                        zorder=6,
                    )
                else:
                    rank = category_ranks.loc[method, column]
                    bar_max_width = 0.72
                    bar_left = x_center - bar_max_width / 2
                    ax.add_patch(
                        Rectangle(
                            (bar_left, y_idx - 0.23),
                            bar_max_width,
                            0.46,
                            facecolor="#e4e4e4",
                            edgecolor="none",
                            zorder=3,
                        )
                    )
                    ax.add_patch(
                        Rectangle(
                            (bar_left, y_idx - 0.23),
                            bar_max_width * float(score),
                            0.46,
                            facecolor=cmap(float(rank)),
                            edgecolor="#222222",
                            linewidth=0.5,
                            zorder=4,
                        )
                    )
                    ax.add_patch(
                        Rectangle(
                            (bar_left, y_idx - 0.23),
                            bar_max_width,
                            0.46,
                            facecolor="none",
                            edgecolor="#cfcfcf",
                            linewidth=0.4,
                            zorder=5,
                        )
                    )

    group_start = 0
    for group, columns in METRIC_GROUPS.items():
        group_end = group_start + len(columns)
        _draw_header(
            ax,
            column_lefts[group_start],
            column_lefts[group_end] + column_widths[group_end],
            group,
            group_palette[group],
        )
        group_start = group_end + 1
    _draw_header(
        ax,
        column_lefts[-1],
        column_lefts[-1] + column_widths[-1],
        "Overall",
        group_palette["Overall"],
        fontsize=8.5,
    )

    ax.set_xlim(method_column_left, row_right)
    ax.set_ylim(len(methods) - 0.5, -1.75)
    ax.set_xticks(column_centers)
    x_labels = []
    for spec in column_specs:
        if spec["kind"] == "metric":
            x_labels.append(_metric_display_name(spec["column"]))
        elif spec["kind"] == "category":
            x_labels.append(f"Overall\n{spec['group'].lower()}")
        else:
            x_labels.append("Overall\nscore")
    ax.set_xticklabels(x_labels, rotation=90, ha="center", fontsize=8)
    for label, spec in zip(ax.get_xticklabels(), column_specs):
        if spec["kind"] in {"category", "overall"}:
            label.set_fontweight("bold")
    ax.set_yticks([])

    ax.tick_params(axis="both", which="both", length=0)
    for y_boundary in np.arange(-0.5, len(methods), 1):
        ax.axhline(y_boundary, color="#ffffff", linewidth=0.7, zorder=1)
    for spine in ax.spines.values():
        spine.set_visible(False)

    _draw_score_legend(ax, score_cmaps)

    plt.tight_layout()
    fig.savefig(save_path, bbox_inches="tight")
    fig.savefig(Path(save_path).with_suffix(".png"), dpi=300, bbox_inches="tight")
    return fig, ax


def build_all():
    raw_metrics = build_raw_table()
    metric_scores = normalize_metrics(raw_metrics)
    category_scores = compute_category_scores(metric_scores)
    raw_metrics, metric_scores, category_scores = _sort_by_overall(raw_metrics, metric_scores, category_scores)

    score_metrics = pd.concat([category_scores, metric_scores], axis=1)
    long_metrics = build_long_table(raw_metrics, metric_scores)

    raw_metrics.to_csv(RAW_TABLE_PATH)
    score_metrics.to_csv(SCORE_TABLE_PATH)
    category_scores.to_csv(CATEGORY_SCORE_TABLE_PATH)
    long_metrics.to_csv(LONG_TABLE_PATH, index=False)
    fig, ax = plot_funky_heatmap(raw_metrics, metric_scores, category_scores)

    return raw_metrics, score_metrics, category_scores, long_metrics, fig, ax


if __name__ == "__main__":
    raw_metrics, score_metrics, category_scores, long_metrics, _, _ = build_all()
    print(f"Wrote {RAW_TABLE_PATH.relative_to(ROOT)}")
    print(f"Wrote {SCORE_TABLE_PATH.relative_to(ROOT)}")
    print(f"Wrote {CATEGORY_SCORE_TABLE_PATH.relative_to(ROOT)}")
    print(f"Wrote {LONG_TABLE_PATH.relative_to(ROOT)}")
    print(f"Wrote {FIGURE_PATH.relative_to(ROOT)}")
    print(f"Wrote {FIGURE_PNG_PATH.relative_to(ROOT)}")
    print(f"Raw table shape: {raw_metrics.shape}")
    print(f"Score table shape: {score_metrics.shape}")
    print(f"Category score table shape: {category_scores.shape}")
