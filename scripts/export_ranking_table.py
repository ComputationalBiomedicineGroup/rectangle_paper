"""
Export expanded cumulative ranking table for Figure 5c methods comparison.
Columns: Method | RMSE per dataset | RMSE cumul. + rank | R per dataset | R cumul. + rank | Mean rank
"""
import re
from pathlib import Path
import pandas as pd

base_dir = Path("data/fig5")
datasets = ["Eng2019_cortex_svz", "Eng2019_ob", "Wang2018_visp"]

proper_method_names = {
    "spotlight": "SPOTlight", "music": "MuSiC", "cell2location": "Cell2location",
    "rctd": "RCTD", "stereoscope": "Stereoscope", "spatialdwls": "SpatialDWLS",
    "destvi": "DestVI", "nnls": "NNLS", "dstg": "DSTG", "seurat": "Seurat",
    "tangram": "Tangram", "stride": "STRIDE", "rectangle": "Rectangle",
}

# ── 1. Load all metrics ──────────────────────────────────────────────────────
rows = []
for dataset in datasets:
    dataset_dir = base_dir / dataset
    if dataset == "Wang2018_visp":
        pattern = re.compile(r"metrics_(.+)_Wang2018_visp_.*$")
        glob_pattern = "metrics_*_Wang2018_visp_*"
    else:
        pattern = re.compile(rf"metrics_(.+)_{re.escape(dataset)}_fov\d+$")
        glob_pattern = f"metrics_*_{dataset}_fov*"
    for path in sorted(dataset_dir.glob(glob_pattern)):
        m = pattern.match(path.name)
        if not m:
            continue
        method = m.group(1)
        df = pd.read_csv(path, sep=r"\s+", engine="python")
        df.columns = [c.strip('"') for c in df.columns]
        rows.append({
            "Set": dataset,
            "Method": method,
            "Correlation": float(df.loc[0, "corr"]),
            "RMSE": float(df.loc[0, "RMSE"]),
        })

all_meths = (
    pd.DataFrame(rows)
    .groupby(["Set", "Method"], as_index=False)
    .agg({"Correlation": "mean", "RMSE": "mean"})
)

# ── 2. Per-dataset rankings ──────────────────────────────────────────────────
method_rankings = all_meths.copy()
method_rankings["Method"] = method_rankings["Method"].map(proper_method_names).fillna(method_rankings["Method"])
method_rankings["rank_rmse"] = method_rankings.groupby("Set")["RMSE"].rank(method="min", ascending=True)
method_rankings["rank_correlation"] = method_rankings.groupby("Set")["Correlation"].rank(method="min", ascending=False)
method_rankings["mean_rank"] = (method_rankings["rank_rmse"] + method_rankings["rank_correlation"]) / 2

# ── 3. Cumulative metrics ────────────────────────────────────────────────────
cumul = method_rankings.groupby("Method", as_index=False).agg(
    cumulative_rmse=("RMSE", "sum"),
    cumulative_correlation=("Correlation", "sum"),
    datasets=("Set", "nunique"),
)
cumul["rank_cumulative_rmse"] = cumul["cumulative_rmse"].rank(method="min", ascending=True)
cumul["rank_cumulative_correlation"] = cumul["cumulative_correlation"].rank(method="min", ascending=False)
cumul["mean_rank"] = (cumul["rank_cumulative_rmse"] + cumul["rank_cumulative_correlation"]) / 2
cumul = cumul.sort_values(
    ["mean_rank", "rank_cumulative_rmse", "rank_cumulative_correlation", "Method"]
).reset_index(drop=True)

# ── 4. Pivot per-dataset values wide ────────────────────────────────────────
dataset_short = {
    "Eng2019_cortex_svz": "Eng2019_cortex_svz",
    "Eng2019_ob": "Eng2019_ob",
    "Wang2018_visp": "Wang2018_visp",
}
_mr = method_rankings.copy()
_mr["Set_short"] = _mr["Set"].map(dataset_short)

rmse_wide = _mr.pivot(index="Method", columns="Set_short", values="RMSE")
rmse_wide.columns = [f"RMSE ({c})" for c in rmse_wide.columns]

corr_wide = _mr.pivot(index="Method", columns="Set_short", values="Correlation")
corr_wide.columns = [f"R ({c})" for c in corr_wide.columns]

_cumul_idx = cumul.set_index("Method")[[
    "cumulative_rmse", "cumulative_correlation",
    "rank_cumulative_rmse", "rank_cumulative_correlation", "mean_rank",
]].rename(columns={
    "cumulative_rmse": "RMSE (cumul.)",
    "cumulative_correlation": "R (cumul.)",
    "rank_cumulative_rmse": "Rank RMSE",
    "rank_cumulative_correlation": "Rank R",
    "mean_rank": "Mean rank",
})

# ── 5. Assemble & order columns ──────────────────────────────────────────────
expanded = (
    rmse_wide
    .join(corr_wide)
    .join(_cumul_idx)
    .sort_values("Mean rank")
    .reset_index()
)

col_order = (
    ["Method"]
    + [f"RMSE ({s})" for s in ["cortex", "ob", "visp"]]
    + ["RMSE (cumul.)", "Rank RMSE"]
    + [f"R ({s})" for s in ["cortex", "ob", "visp"]]
    + ["R (cumul.)", "Rank R", "Mean rank"]
)
expanded = expanded[col_order]

float_cols = [c for c in expanded.columns if c != "Method"]
expanded[float_cols] = expanded[float_cols].round(4)

# ── 6. Print & export ────────────────────────────────────────────────────────
print(expanded.to_string(index=False))

out_path = Path("data/fig5/cumulative_ranking_table.csv")
out_path.parent.mkdir(parents=True, exist_ok=True)
expanded.to_csv(out_path, index=False)
print(f"\nSaved → {out_path}  ({len(expanded)} methods × {len(expanded.columns)} columns)")

