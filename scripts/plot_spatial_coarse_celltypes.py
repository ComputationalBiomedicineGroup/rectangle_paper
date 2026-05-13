#!/usr/bin/env python3
"""Create per-celltype coarse spatial plots from Rectangle estimations.

This script groups fine cell types into coarse labels by splitting each column name
at the first underscore (e.g. ``Ext_L23`` -> ``Ext``), then generates one spatial
PDF per coarse cell type using the same visual style as the existing Figure 5b plot.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import spatialdata as sd
from anndata import AnnData
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Colormap, Normalize
from spatialdata.models import TableModel

_AX_BG_COLOR = "white"


def transparent_sequential_cmap(
    color: str,
    n: int = 256,
    bg: str = "white",
    name: str | None = None,
) -> mpl.colors.ListedColormap:
    """Create a cmap from `bg` to `color` where value 0 is fully transparent."""
    rgb_target = np.array(mpl.colors.to_rgb(color), dtype=float)
    if bg == "white":
        bg_color = np.ones(3, dtype=float)
    elif bg == "black":
        bg_color = np.zeros(3, dtype=float)
    else:
        raise ValueError(f"bg must be 'white' or 'black', got {bg!r}")

    t = np.linspace(0.0, 1.0, n)
    rgb = (1.0 - t)[:, None] * bg_color[None, :] + t[:, None] * rgb_target[None, :]
    alpha = t
    rgba = np.concatenate([rgb, alpha[:, None]], axis=1)
    return mpl.colors.ListedColormap(rgba, name or f"{color}_transparent_seq")


def add_celltype_colorbar(
    ax: plt.Axes,
    values: np.ndarray,
    cmap: Colormap,
    title: str,
    *,
    bg: str = "white",
) -> None:
    """Add a compact horizontal colorbar with min/max ticks."""
    vals = np.asarray(values, float)
    vals = vals[~np.isnan(vals)]
    if vals.size == 0:
        vmin, vmax = 0.0, 1.0
    else:
        vmin, vmax = float(vals.min()), float(vals.max())
        if vmin == vmax:
            vmin -= 0.5
            vmax += 0.5

    base = cmap(np.linspace(0.0, 1.0, 256))
    base[:, 3] = 1.0
    cmap_opaque = mpl.colors.ListedColormap(base, name=f"{cmap.name}_opaque")

    ax.set_facecolor(bg)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)

    bar_ax = ax.inset_axes([0.05, 0.1, 0.9, 0.6])
    cbar = ColorbarBase(
        bar_ax,
        cmap=cmap_opaque,
        norm=Normalize(vmin=vmin, vmax=vmax),
        orientation="horizontal",
        extend="both",
        extendfrac=0.05,
    )
    cbar.set_ticks([vmin, vmax])
    cbar.set_ticklabels([f"{vmin:.2f}", f"{vmax:.2f}"])
    cbar.ax.tick_params(colors="black", labelsize=8, length=2, pad=2)

    ax.text(0.5, 0.82, title, ha="center", va="bottom", color="black", fontsize=9, transform=ax.transAxes)


def coarse_grain_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate columns by prefix before the first underscore."""
    coarse_names = [col.split("_", 1)[0] for col in df.columns]
    coarse_df = df.T.groupby(coarse_names).sum().T
    return coarse_df


def add_results_table_to_spatialdata(
    sdata,
    estimations: pd.DataFrame,
    *,
    table_name: str = "rectangle_results_coarse",
    region_name: str = "ST8059048",
    source_table: str = "table",
) -> None:
    """Add coarse estimations as a SpatialData table aligned by spot id."""
    table_data = AnnData(estimations)
    adata = TableModel.parse(table_data)

    spot_id_df = adata.obs.merge(
        sdata.tables[source_table].obs,
        left_index=True,
        right_index=True,
        how="left",
    )
    good_idx = spot_id_df.index[~spot_id_df["spot_id"].isna()]
    adata = adata[good_idx, :].copy()

    adata.uns["spatialdata_attrs"] = {
        "region": region_name,
        "region_key": "region",
        "instance_key": "spot_id",
    }
    adata.obs["region"] = pd.Categorical([region_name] * len(adata))
    adata.obs["spot_id"] = spot_id_df.loc[good_idx, "spot_id"].astype(int).values

    adata.X = adata.to_df().fillna(0).values
    sdata.tables[table_name] = adata


def plot_individual_coarse_celltypes(
    sdata,
    table_name: str,
    coarse_celltypes: list[str],
    output_dir: Path,
) -> None:
    """Generate one PDF per coarse cell type with a single overlay layer."""
    output_dir.mkdir(parents=True, exist_ok=True)
    palette = mpl.cm.get_cmap("tab20", len(coarse_celltypes))

    table_df = sdata.tables[table_name].to_df()
    sdata[table_name].obs["is_bad"] = (table_df.sum(axis=1) == 0).astype(str)

    for idx, celltype in enumerate(coarse_celltypes):
        color = mpl.colors.to_hex(palette(idx))
        cmap = transparent_sequential_cmap(color, n=256, bg=_AX_BG_COLOR, name=f"{celltype}_transparent_seq")

        fig = plt.figure(figsize=(4, 6), facecolor="white")
        gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[5, 1], hspace=0.2)
        top_ax = fig.add_subplot(gs[0])
        bottom_ax = fig.add_subplot(gs[1])

        (
            sdata.pl.render_shapes(
                color="is_bad",
                groups=["True", "False"],
                palette=["gray", "white"],
                shape="visium_hex",
                table_name=table_name,
            ).pl.show("ST8059048", ax=top_ax, legend_loc=None, colorbar=False)
        )

        (
            sdata.pl.render_shapes(
                color=celltype,
                cmap=cmap,
                shape="visium_hex",
                table_name=table_name,
            ).pl.show("ST8059048", ax=top_ax, colorbar=False)
        )

        top_ax.set_facecolor("white")
        top_ax.set_title(celltype, fontsize=11)
        top_ax.set_xlabel("")
        top_ax.set_ylabel("")

        add_celltype_colorbar(bottom_ax, table_df[celltype].to_numpy(), cmap, celltype, bg="white")

        output_path = output_dir / f"5b_spatial_median_coarse_{celltype}.pdf"
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sdata-path",
        type=Path,
        default=Path("data/ST8059048.zarr"),
        help="Path to ST8059048 SpatialData zarr",
    )
    parser.add_argument(
        "--estimations-csv",
        type=Path,
        default=Path("data/fig5/rectangle_ST8059048_estimations_median.csv"),
        help="Rectangle estimations CSV",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("figures/coarse_celltypes"),
        help="Output folder for generated PDFs",
    )
    parser.add_argument(
        "--min-spot-count",
        type=float,
        default=500,
        help="Drop spots with total bulk count lower than this threshold",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Only print coarse cell types and exit",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    sdata = sd.read_zarr(args.sdata_path)
    data_table = sdata["table"]
    bulks = data_table.to_df()
    counts_per_row = bulks.sum(axis=1)

    estimations = pd.read_csv(args.estimations_csv, index_col=0)
    low_count_spots = counts_per_row[counts_per_row < args.min_spot_count].index
    estimations = estimations.drop(index=low_count_spots, errors="ignore")

    coarse_estimations = coarse_grain_columns(estimations)
    coarse_celltypes = coarse_estimations.columns.tolist()
    print(f"Coarse cell types ({len(coarse_celltypes)}): {coarse_celltypes}")

    if args.dry_run:
        return

    add_results_table_to_spatialdata(sdata, coarse_estimations)
    plot_individual_coarse_celltypes(
        sdata=sdata,
        table_name="rectangle_results_coarse",
        coarse_celltypes=coarse_celltypes,
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()

