"""Supplementary figure: runtime and peak memory of every method against the number of bulk samples.

Same HaoSub reference (14,744 cells) for all three bulk datasets, run with scripts/wu20k_benchmark.py:
* Hoek: 8 bulks, output/hoek_benchmark/results.csv
* Dietrich: 59 bulks, output/dietrich_benchmark/results.csv
* Morandini: 159 bulks, output/morandini_benchmark/results.csv
Mean +- SD over the successful runs (three, or one when the first run took longer than 45 min).
Failed runs (6 h time limit) are drawn as crosses on the limit. DWLS ran natively (--native-r): its MAST step needs
~45-48 GB, more than the 30.3 GiB Docker VM.
Run from the repository root:
    .venv/bin/python scripts/suppl_method_scalability.py
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, NullLocator

ROOT = Path(__file__).resolve().parents[1]
FIGURE_PATH = ROOT / "figures" / "fig_suppl_method_scalability.pdf"
TABLE_PATH = ROOT / "data" / "fig2" / "runtimes" / "method_scalability_summary.csv"

# Figure 2 colours
METHOD_COLORS = {
    "Rectangle": "#000000",
    "SCDC": "#B3B3B3",
    "AutoGeneS": "#66C2A5",
    "MuSiC": "#FFD92F",
    "Scaden": "#E5C494",
    "BayesPrism": "#FC8D62",
    "Bisque": "#8DA0CB",
    "DWLS": "#A6D854",
    "CIBERSORTx": "#E78AC3",
    "SECRET": "#A6761D",
}
METHOD_LABELS = {method.lower(): method for method in METHOD_COLORS}
DATASETS = {
    "Hoek": "Hoek\n8 bulks",
    "Dietrich": "Dietrich\n59 bulks",
    "Morandini": "Morandini\n159 bulks",
}
TIME_LIMIT_S = 6 * 3600
INK = "#222222"


def load_summary():
    rows = []
    for dataset in DATASETS:
        runs = pd.read_csv(ROOT / "output" / f"{dataset.lower()}_benchmark" / "results.csv")
        runs["method"] = runs["method"].map(METHOD_LABELS)
        ok = runs[runs["status"] == "ok"].groupby("method")
        summary = pd.DataFrame({
            "runtime_mean": ok["runtime_s"].mean(),
            "runtime_sd": ok["runtime_s"].std(),
            "memory_mean": ok["peak_memory_mib"].mean() / 1024,
            "memory_sd": ok["peak_memory_mib"].std() / 1024,
            "n_runs": ok.size(),
        })
        failed = runs[~runs["method"].isin(summary.index)].groupby("method")["status"].first()
        summary = summary.reindex(summary.index.union(failed.index))
        summary["status"] = failed.reindex(summary.index).fillna("ok")
        rows.append(summary.rename_axis("method").reset_index().assign(dataset=dataset))
    summary = pd.concat(rows, ignore_index=True)
    assert summary["method"].notna().all() and set(summary["method"]) == set(METHOD_COLORS)
    return summary


def plot(summary):
    datasets = list(DATASETS)
    x = {dataset: i for i, dataset in enumerate(datasets)}
    # Legend in the order of the Morandini runtime; failed methods last
    morandini = summary[summary["dataset"] == "Morandini"].set_index("method")
    order = morandini.assign(_t=morandini["runtime_mean"].fillna(1e9)).sort_values("_t").index.tolist()

    style = {
        "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
        "axes.linewidth": 0.7, "axes.edgecolor": INK, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5, "pdf.fonttype": 42, "ps.fonttype": 42,
    }
    with mpl.rc_context(style):
        fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), layout="constrained")
        panels = [
            (axes[0], "runtime", "Runtime (log scale)", TIME_LIMIT_S, "6 h limit", "timeout"),
            (axes[1], "memory", "Peak memory (GiB)", None, None, None),
        ]
        for ax, key, ylabel, limit, limit_label, failure in panels:
            if key == "runtime":  # no run hit a memory limit
                ax.axhline(limit, color="#999999", linewidth=0.7, linestyle=(0, (3, 2)), zorder=1)
                # Between the first two datasets, where no failure cross sits on the line
                ax.text(0.5, limit, limit_label, ha="center", va="bottom", fontsize=6.5, color="#666666")
            for method in order:
                data = summary[summary["method"] == method].set_index("dataset").reindex(datasets)
                color = METHOD_COLORS[method]
                edge = INK if method in {"MuSiC", "SCDC", "Scaden"} else color  # light colours need an outline
                width = 1.8 if method == "Rectangle" else 1.1
                xs = [x[d] for d in datasets]
                ax.plot(xs, data[f"{key}_mean"], color=color, linewidth=width, zorder=3 if method == "Rectangle" else 2)
                # Filled: mean of repeated runs; open: a single run (no SD). The SD bars are drawn in ink over the
                # markers, so they stay visible where they are shorter than a marker
                single = data["n_runs"] == 1
                ax.plot(xs, data[f"{key}_mean"].where(~single), "o", markersize=4, color=color,
                        markeredgecolor=edge, markeredgewidth=0.5, zorder=4)
                ax.plot(xs, data[f"{key}_mean"].where(single), "o", markersize=4, color="white",
                        markeredgecolor=color if edge == color else INK, markeredgewidth=1.0, zorder=4)
                ax.errorbar(xs, data[f"{key}_mean"], yerr=data[f"{key}_sd"], fmt="none", ecolor=INK, elinewidth=0.7,
                            capsize=2.5, capthick=0.7, zorder=6)
                hit = data["status"].astype(str).str.contains(failure) if failure else data["status"].isna()
                for dataset in data.index[hit]:
                    ax.plot(x[dataset], limit, marker="X", markersize=7, color=color, markeredgecolor=INK,
                            markeredgewidth=0.5, zorder=5, clip_on=False)
            ax.set_xticks(range(len(datasets)))
            ax.set_xticklabels(DATASETS.values())
            ax.set_xlim(-0.3, len(datasets) - 0.5)
            ax.set_ylabel(ylabel)
            ax.spines[["top", "right"]].set_visible(False)

        runtime_ax, memory_ax = axes
        runtime_ax.set_yscale("log")
        ticks = {10: "10 s", 60: "1 min", 600: "10 min", 3600: "1 h", TIME_LIMIT_S: "6 h"}
        runtime_ax.yaxis.set_major_locator(FixedLocator(list(ticks)))
        runtime_ax.yaxis.set_minor_locator(NullLocator())
        runtime_ax.set_yticklabels(ticks.values())
        runtime_ax.set_ylim(8, TIME_LIMIT_S * 1.6)
        memory_ax.set_ylim(bottom=0)

        handles = [Line2D([0], [0], color=METHOD_COLORS[m], marker="o", markersize=4, linewidth=1.8 if m == "Rectangle" else 1.1,
                          markeredgecolor=INK if m in {"MuSiC", "SCDC", "Scaden"} else METHOD_COLORS[m], markeredgewidth=0.5,
                          label=m) for m in order]
        handles.append(Line2D([0], [0], linestyle="none", marker="o", markersize=4, color="#999999",
                              markeredgecolor=INK, markeredgewidth=0.5, label="Mean ± SD, 3 runs"))
        handles.append(Line2D([0], [0], linestyle="none", marker="o", markersize=4, color="white",
                              markeredgecolor=INK, markeredgewidth=1.0, label="Single run"))
        handles.append(Line2D([0], [0], linestyle="none", marker="X", markersize=7, color="#cccccc",
                              markeredgecolor=INK, markeredgewidth=0.5, label="Failed run"))
        legend = fig.legend(handles=handles, loc="outside right center", frameon=False, fontsize=7, handlelength=1.6)
        for text in legend.get_texts():
            if text.get_text() == "Rectangle":
                text.set_fontweight("bold")

    fig.savefig(FIGURE_PATH, bbox_inches="tight")
    fig.savefig(FIGURE_PATH.with_suffix(".png"), dpi=300, bbox_inches="tight")
    return fig


if __name__ == "__main__":
    summary = load_summary()
    summary.to_csv(TABLE_PATH, index=False)
    plot(summary)
    print(summary.sort_values(["method", "dataset"]).round(2).to_string(index=False))
    print(f"Wrote {FIGURE_PATH.relative_to(ROOT)} (+ .png) and {TABLE_PATH.relative_to(ROOT)}")
