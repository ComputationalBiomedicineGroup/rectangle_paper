# Rectangle Paper Utilities

## Coarse Spatial Plots

Generate per-celltype spatial plots where fine labels are collapsed to the prefix before `_` (for example, `Ext_L23` -> `Ext`):

```bash
python scripts/plot_spatial_coarse_celltypes.py
```

Dry run (prints detected coarse cell types only):

```bash
python scripts/plot_spatial_coarse_celltypes.py --dry-run
```

Outputs are written to `figures/coarse_celltypes/` as one PDF per coarse cell type.

