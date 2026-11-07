"""
Combine lambrechts_fine_method_rmse_all.csv and wu_fine_method_rmse_all.csv
into:
  1. output/rmse_combined_long.csv        – long format (Dataset | Method | Cell type | RMSE)
  2. output/rmse_combined_mean_summary.csv – mean RMSE per method × dataset (wide, paper-ready)
"""

import os
import pandas as pd
import numpy as np

DATA   = os.path.join(os.path.dirname(__file__), '..', 'data', 'fig3')
OUT    = os.path.join(os.path.dirname(__file__), '..', 'output')
os.makedirs(OUT, exist_ok=True)

METHOD_ORDER = ['Rectangle', 'SCDC', 'AutoGeneS', 'MuSiC', 'Scaden',
                'BayesPrism', 'Bisque', 'DWLS', 'CIBERSORTx']

# ── load ─────────────────────────────────────────────────────────────────────
lam = pd.read_csv(f'{DATA}/lambrechts_fine_method_rmse_all.csv', index_col=0)
wu  = pd.read_csv(f'{DATA}/wu_fine_method_rmse_all.csv',         index_col=0)

# ── 1. Long format ───────────────────────────────────────────────────────────
def to_long(df, dataset_label):
    long = (
        df
        .reset_index()
        .rename(columns={'index': 'Method'})
        .melt(id_vars='Method', var_name='Cell type', value_name='RMSE')
    )
    long.insert(0, 'Dataset', dataset_label)
    return long

lam_long = to_long(lam, 'Lambrechts')
wu_long  = to_long(wu,  'Wu')

combined_long = pd.concat([lam_long, wu_long], ignore_index=True)
combined_long['Method'] = pd.Categorical(
    combined_long['Method'], categories=METHOD_ORDER, ordered=True
)
combined_long = combined_long.sort_values(['Dataset', 'Method', 'Cell type']).reset_index(drop=True)
combined_long['RMSE'] = combined_long['RMSE'].round(6)

combined_long.to_csv(f'{OUT}/rmse_combined_long.csv', index=False)
print(f"Long format saved  → {OUT}/rmse_combined_long.csv  ({len(combined_long)} rows)")

# ── 2. Mean RMSE summary (paper table) ───────────────────────────────────────
def mean_rmse(df, dataset_label):
    s = df.mean(axis=1).rename('Mean RMSE')
    s.index.name = 'Method'
    out = s.reset_index()
    out.insert(0, 'Dataset', dataset_label)
    return out

summary = pd.concat([mean_rmse(lam, 'Lambrechts'), mean_rmse(wu, 'Wu')], ignore_index=True)

# Pivot to wide: rows = Method, columns = Dataset
pivot = (
    summary
    .pivot(index='Method', columns='Dataset', values='Mean RMSE')
    .rename_axis(None, axis=1)
    .reset_index()
)
# reorder methods
pivot['_ord'] = pivot['Method'].map({m: i for i, m in enumerate(METHOD_ORDER)})
pivot = pivot.sort_values('_ord').drop(columns='_ord').reset_index(drop=True)
pivot['Lambrechts'] = pivot['Lambrechts'].round(6)
pivot['Wu']         = pivot['Wu'].round(6)

pivot.to_csv(f'{OUT}/rmse_combined_mean_summary.csv', index=False)

# ── print ─────────────────────────────────────────────────────────────────────
sep = '-' * 60
print()
print('Mean RMSE per method (fine cell types)')
print(sep)
print(pivot.to_string(index=False))
print()
print(f"Summary saved      → {OUT}/rmse_combined_mean_summary.csv")

