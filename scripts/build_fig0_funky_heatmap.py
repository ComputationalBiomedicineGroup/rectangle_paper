import os
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
FIGURE_PATH = FIGURE_DIR / "fig0_funky_heatmap_prototype.pdf"

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
]

DATASET_LABELS = {
    "chen": "Chen",
    "finotello": "Finotello",
    "hoek": "Hoek",
    "morandini": "Morandini",
    "petitprez": "Petitprez",
    "dietrich": "Dietrich",
}

FIG2_DATASETS = ["chen", "finotello", "hoek", "morandini", "petitprez", "dietrich"]

METRIC_GROUPS = {
    "Accuracy": [
        *(f"accuracy_rmse_{dataset}" for dataset in FIG2_DATASETS),
        *(f"accuracy_corr_{dataset}" for dataset in FIG2_DATASETS),
    ],
    "Efficiency": ["runtime_hao_20k_s", "memory_hao_20k_mib"],
    "Spillover": ["spillover_pct_hao", "spillover_pct_hoek", "spillover_pct_linsley"],
    "Resolution": ["resolution_rmse_wu", "resolution_rmse_lambrechts"],
    "Unknown": ["unknown_rmse_all", "unknown_corr_all"],
}

CATEGORY_SCORE_COLUMNS = {
    group: f"score_{group.lower()}" for group in METRIC_GROUPS
}

LOWER_IS_BETTER_PATTERNS = ("rmse", "runtime", "memory", "spillover")


def _canonical_method(value):
    key = str(value).strip().replace(" ", "").replace("-", "").lower()
    return METHOD_LABELS.get(key, str(value).strip())


def _empty_raw_table():
    methods = pd.Index(METHOD_ORDER, name="method")
    columns = [column for group_columns in METRIC_GROUPS.values() for column in group_columns]
    return pd.DataFrame(index=methods, columns=columns, dtype=float)


def _add_accuracy_metrics(raw):
    metrics = pd.read_csv(ROOT / "data" / "fig2" / "cell_type_level_performance_metrics.csv")
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
    benchmark = pd.read_csv(ROOT / "output" / "benchmark_table_20k_ref.csv")
    benchmark["method"] = benchmark["Method"].map(_canonical_method)
    benchmark = benchmark.set_index("method")
    for method in raw.index.intersection(benchmark.index):
        raw.loc[method, "runtime_hao_20k_s"] = benchmark.loc[method, "Runtime mean (s)"]
        raw.loc[method, "memory_hao_20k_mib"] = benchmark.loc[method, "Memory mean (MiB)"]


def _add_spillover_metrics(raw):
    spillover = pd.read_csv(ROOT / "data" / "fig3" / "spillover_total.csv")
    spillover["method"] = spillover["Method"].map(_canonical_method)
    spillover = spillover.set_index("method")
    for method in raw.index.intersection(spillover.index):
        raw.loc[method, "spillover_pct_hao"] = spillover.loc[method, "Hao"] * 100
        raw.loc[method, "spillover_pct_hoek"] = spillover.loc[method, "Hoek"] * 100
        raw.loc[method, "spillover_pct_linsley"] = spillover.loc[method, "Linsley"] * 100


def _resolution_sum(path):
    df = pd.read_csv(path)
    df = df.rename(columns={df.columns[0]: "method"})
    df["method"] = df["method"].map(_canonical_method)
    numeric = df.drop(columns=["method"]).apply(pd.to_numeric, errors="coerce")
    values = numeric.sum(axis=1, min_count=1)
    return pd.Series(values.to_numpy(), index=df["method"])


def _add_resolution_metrics(raw):
    resolution_sources = {
        "resolution_rmse_wu": ROOT / "data" / "fig3" / "wu_fine_method_rmse_all.csv",
        "resolution_rmse_lambrechts": ROOT / "data" / "fig3" / "lambrechts_fine_method_rmse_all.csv",
    }
    for column, path in resolution_sources.items():
        values = _resolution_sum(path)
        for method in raw.index.intersection(values.index):
            raw.loc[method, column] = values.loc[method]


def _add_unknown_metrics(raw):
    unknown = pd.read_csv(ROOT / "data" / "fig4" / "unknown_content_results.csv")
    unknown = unknown[unknown["unknown_fraction"].isin([0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9])].copy()
    unknown["method"] = unknown["method"].map(_canonical_method)
    unknown = unknown[unknown["method"].isin(raw.index)]

    rows = []
    for method, method_df in unknown.groupby("method", observed=True):
        true_values = method_df["true_value"].astype(float)
        predicted_values = method_df["predicted_value"].astype(float)
        rows.append(
            {
                "method": method,
                "unknown_rmse_all": float(np.sqrt(np.mean((predicted_values - true_values) ** 2))),
                "unknown_corr_all": float(np.corrcoef(true_values, predicted_values)[0, 1]),
            }
        )
    unknown_summary = pd.DataFrame(rows).set_index("method")
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
        "accuracy_corr_": "Corr ",
        "runtime_hao_20k_s": "Runtime Hao 20k",
        "memory_hao_20k_mib": "Memory Hao 20k",
        "spillover_pct_": "Spillover ",
        "resolution_rmse_": "Fine RMSE ",
        "unknown_rmse_all": "Unknown RMSE",
        "unknown_corr_all": "Unknown corr",
    }
    label = metric
    for old, new in replacements.items():
        label = label.replace(old, new)
    for key, value in DATASET_LABELS.items():
        label = label.replace(key, value)
    return label.replace("_", " ").title().replace("Rmse", "RMSE").replace("Corr", "corr")


def _group_palette():
    return {
        "Overall": "#53B8C5",
        "Accuracy": "#8D98C8",
        "Efficiency": "#94B6CE",
        "Spillover": "#EBA0A2",
        "Resolution": "#DD6AA8",
        "Unknown": "#E99BC5",
    }


def _score_cmap():
    return LinearSegmentedColormap.from_list(
        "fig0_score",
        ["#eef1c7", "#b8d5cf", "#47b9c7", "#33499d", "#4b197f"],
    )


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

    cmap = _score_cmap()
    group_palette = _group_palette()

    y_positions = np.arange(len(methods))

    for y_idx, method in enumerate(methods):
        row_color = "#f2f2f2" if y_idx % 2 == 0 else "#ffffff"
        ax.add_patch(
            Rectangle(
                (-0.5, y_idx - 0.5),
                total_width,
                1.0,
                facecolor=row_color,
                edgecolor="none",
                zorder=0,
            )
        )

    for y_idx, method in enumerate(methods):
        for x_idx, spec in enumerate(column_specs):
            x_center = column_centers[x_idx]
            kind = spec["kind"]
            column = spec["column"]
            if kind == "metric":
                score = metric_scores.loc[method, column]
                if pd.isna(score):
                    ax.scatter(x_center, y_idx, marker="x", s=24, color="#b8b8b8", linewidths=0.9, zorder=5)
                    continue
                size = 8 + float(score) * 265
                ax.scatter(
                    x_center,
                    y_idx,
                    s=size,
                    c=[cmap(float(score))],
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
                    bar_left = column_lefts[x_idx] + 0.28
                    bar_max_width = column_widths[x_idx] - 0.78
                    bar_width = bar_max_width * float(score)
                    ax.add_patch(
                        Rectangle(
                            (bar_left, y_idx - 0.33),
                            bar_width,
                            0.66,
                            facecolor=cmap(float(score)),
                            edgecolor="#222222",
                            linewidth=0.55,
                            zorder=4,
                        )
                    )
                    ax.text(
                        column_lefts[x_idx] + column_widths[x_idx] - 0.08,
                        y_idx,
                        f"{category_scores.loc[method, 'overall_score']:.1f}",
                        ha="right",
                        va="center",
                        fontsize=8.0,
                        color="#111111",
                        zorder=6,
                    )
                else:
                    ax.add_patch(
                        Rectangle(
                            (x_center - 0.36, y_idx - 0.36),
                            0.72,
                            0.72,
                            facecolor=cmap(float(score)),
                            edgecolor="#222222",
                            linewidth=0.55,
                            zorder=4,
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
        ax.axvline(column_lefts[group_end + 1], color="#222222", linewidth=0.9, alpha=0.7, zorder=3)
        group_start = group_end + 1
    _draw_header(
        ax,
        column_lefts[-1],
        column_lefts[-1] + column_widths[-1],
        "Overall",
        group_palette["Overall"],
        fontsize=8.5,
    )

    ax.set_xlim(-0.5, total_width + 0.45)
    ax.set_ylim(len(methods) - 0.5, -1.75)
    ax.set_xticks(column_centers)
    x_labels = []
    for spec in column_specs:
        if spec["kind"] == "metric":
            x_labels.append(_metric_display_name(spec["column"]))
        elif spec["kind"] == "category":
            x_labels.append("Weighted\ncategory")
        else:
            x_labels.append("Overall\nscore")
    ax.set_xticklabels(x_labels, rotation=90, ha="center", fontsize=8)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(methods, fontsize=10)

    for label in ax.get_yticklabels():
        if label.get_text() == "Rectangle":
            label.set_fontweight("bold")
            label.set_color("#6b4f00")

    ax.tick_params(axis="both", length=0)
    ax.set_xticks(np.r_[column_lefts - 0.5, total_width - 0.5], minor=True)
    ax.set_yticks(np.arange(-0.5, len(methods), 1), minor=True)
    ax.grid(which="minor", color="#ffffff", linewidth=0.7, zorder=1)
    for spine in ax.spines.values():
        spine.set_visible(False)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 1))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.022, pad=0.03)
    cbar.set_label("Score", rotation=90)

    legend_scores = [0.25, 0.5, 0.75, 1.0]
    handles = [
        ax.scatter([], [], s=8 + score * 265, color=cmap(score), edgecolors="#202020", linewidths=0.45, label=f"{score:.2g}")
        for score in legend_scores
    ]
    handles.append(ax.scatter([], [], marker="x", s=24, color="#b8b8b8", linewidths=0.9, label="NA"))
    ax.legend(handles=handles, title="Metric score", frameon=False, bbox_to_anchor=(1.14, 0.52), loc="center left")

    ax.set_title("Funky heatmap prototype: weighted deconvolution performance", loc="left", fontsize=13, fontweight="bold", pad=24)
    plt.tight_layout()
    fig.savefig(save_path, bbox_inches="tight")
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
    print(f"Raw table shape: {raw_metrics.shape}")
    print(f"Score table shape: {score_metrics.shape}")
    print(f"Category score table shape: {category_scores.shape}")
