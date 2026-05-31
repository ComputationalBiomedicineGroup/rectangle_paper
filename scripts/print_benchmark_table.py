"""
Print runtime & memory benchmark tables for Figure 3C.

- Table A: Full reference  (fig3c_scdc_rectangle_comparison_fail_no_log)
- Table B: Subsampled 20k  (fig3c_scdc_rectangle_comparison_20k_no_log)
"""

import os
import pandas as pd
import numpy as np

DATA = os.path.join(os.path.dirname(__file__), '..', 'data', 'fig3')


def fix_mem(mean_val, sd_val=None):
    """If memory was stored as machine words, convert to MiB."""
    if mean_val > 1e6:
        factor = 8 / (1024 * 1024)
        return mean_val * factor, (sd_val * factor if sd_val is not None else None)
    return mean_val, sd_val


# ---------------------------------------------------------------------------
# TABLE A  --  Full reference  (up to 147 k cells)
# ---------------------------------------------------------------------------
rect_mem_s = pd.read_csv(f'{DATA}/rectangle_memory_summary_qs.csv')
rect_run_s = pd.read_csv(f'{DATA}/rectangle_runtime_summary_qs.csv')
ag_mem_s   = pd.read_csv(f'{DATA}/autogenes_memory_summary.csv')
ag_run_s   = pd.read_csv(f'{DATA}/autogenes_runtime_summary.csv')
scdc       = pd.read_csv(f'{DATA}/scdc.csv')
music      = pd.read_csv(f'{DATA}/music.csv')
scaden     = pd.read_csv(f'{DATA}/scaden.csv')
dwls       = pd.read_csv(f'{DATA}/dwls_benchmark_results.csv')
bisque     = pd.read_csv(f'{DATA}/bisque_benchmark_results.csv')
bp         = pd.read_csv(f'{DATA}/bayesprism_benchmark_results.csv')
cibx_rt    = pd.read_csv(f'{DATA}/cibersortx_runtime_80k.csv')
cibx_mem   = pd.read_csv(f'{DATA}/cibersortx_memory_80k.csv')

n_full = int(rect_mem_s['n_cells'].max())
rect_mem_row = rect_mem_s[rect_mem_s['n_cells'] == n_full].iloc[0]
rect_run_row = rect_run_s[rect_run_s['n_cells'] == n_full].iloc[0]
ag_mem_row   = ag_mem_s[ag_mem_s['n_cells'] == n_full].iloc[0]
ag_run_row   = ag_run_s[ag_run_s['n_cells'] == n_full].iloc[0]

scdc_mem, scdc_mem_sd = fix_mem(scdc['Peak_Mem_MiB'].mean(), scdc['Peak_Mem_MiB'].std())
dwls_mem, dwls_mem_sd = fix_mem(dwls['Peak_Mem_MiB'].mean(), dwls['Peak_Mem_MiB'].std())

rows_full = [
    dict(Method='Rectangle',  Ref_size=f'{n_full:,}',
         Runtime_mean_s=round(rect_run_row['time_mean'], 1),
         Runtime_SD_s=round(rect_run_row['time_sd'], 1),
         Memory_mean_MiB=round(rect_mem_row['memory_mean'], 0),
         Memory_SD_MiB=round(rect_mem_row['memory_sd'], 0),
         Note=''),
    dict(Method='SCDC',       Ref_size=f'{n_full:,}',
         Runtime_mean_s=round(scdc['Time_s'].mean(), 1),
         Runtime_SD_s=round(scdc['Time_s'].std(), 1),
         Memory_mean_MiB=round(scdc_mem, 0),
         Memory_SD_MiB=round(scdc_mem_sd, 0),
         Note=''),
    dict(Method='Autogenes',  Ref_size=f'{n_full:,}',
         Runtime_mean_s=round(ag_run_row['time_mean'], 1),
         Runtime_SD_s=round(ag_run_row['time_sd'], 1),
         Memory_mean_MiB=round(ag_mem_row['memory_mean'], 0),
         Memory_SD_MiB=round(ag_mem_row['memory_sd'], 0),
         Note=''),
    dict(Method='MuSiC',      Ref_size=f'{n_full:,}',
         Runtime_mean_s=round(music['Time_s'].mean(), 1),
         Runtime_SD_s=round(music['Time_s'].std(), 1),
         Memory_mean_MiB=round(music['Peak_Mem_MiB'].mean(), 0),
         Memory_SD_MiB=round(music['Peak_Mem_MiB'].std(), 0),
         Note=''),
    dict(Method='Scaden',     Ref_size='<=90k',
         Runtime_mean_s=round(scaden['Time_s'].mean(), 1),
         Runtime_SD_s=round(scaden['Time_s'].std(), 1),
         Memory_mean_MiB=round(scaden['Peak_Mem_MiB'].mean(), 0),
         Memory_SD_MiB=round(scaden['Peak_Mem_MiB'].std(), 0),
         Note='No solution for >90k cells'),
    dict(Method='BayesPrism', Ref_size='<=50k',
         Runtime_mean_s=round(bp['Time_s'].mean(), 1),
         Runtime_SD_s=round(bp['Time_s'].std(), 1),
         Memory_mean_MiB=round(bp['Peak_Mem_MiB'].mean(), 0),
         Memory_SD_MiB=round(bp['Peak_Mem_MiB'].std(), 0),
         Note='No solution for >50k cells'),
    dict(Method='Bisque',     Ref_size='<=35k',
         Runtime_mean_s=round(bisque['Time_s'].mean(), 1),
         Runtime_SD_s=round(bisque['Time_s'].std(), 1),
         Memory_mean_MiB=round(bisque['Peak_Mem_MiB'].mean(), 0),
         Memory_SD_MiB=round(bisque['Peak_Mem_MiB'].std(), 0),
         Note='No solution for >35k cells'),
    dict(Method='DWLS',       Ref_size='<=20k',
         Runtime_mean_s=round(dwls['Time_s'].mean(), 1),
         Runtime_SD_s=round(dwls['Time_s'].std(), 1),
         Memory_mean_MiB=round(dwls_mem, 0),
         Memory_SD_MiB=round(dwls_mem_sd, 0),
         Note='No solution for >20k cells'),
    dict(Method='CIBERSORTx', Ref_size='<=80k',
         Runtime_mean_s=round(cibx_rt['time'].mean(), 1),
         Runtime_SD_s=round(cibx_rt['time'].std(), 1),
         Memory_mean_MiB=round(cibx_mem['memory_mib'].mean(), 0),
         Memory_SD_MiB=round(cibx_mem['memory_mib'].std(), 0),
         Note='No solution for >80k cells (cloud-only tool)'),
]

df_full = pd.DataFrame(rows_full)
df_full.columns = ['Method', 'Ref size', 'Runtime mean (s)', 'Runtime SD (s)',
                   'Memory mean (MiB)', 'Memory SD (MiB)', 'Note']

# ---------------------------------------------------------------------------
# TABLE B  --  Subsampled reference  (20 000 cells)
# ---------------------------------------------------------------------------
rect_mem_s2 = pd.read_csv(f'{DATA}/rectangle_memory_summary_qs.csv')
rect_run_s2 = pd.read_csv(f'{DATA}/rectangle_runtime_summary_qs.csv')
ag_mem_s2   = pd.read_csv(f'{DATA}/autogenes_memory_summary_20k.csv')
ag_run_s2   = pd.read_csv(f'{DATA}/autogenes_runtime_summary_20k.csv')
scdc2       = pd.read_csv(f'{DATA}/scdc_benchmark_results_20k.csv')
music2      = pd.read_csv(f'{DATA}/music_benchmark_results_20k.csv')
scaden2     = pd.read_csv(f'{DATA}/scaden_20k.csv')
dwls2       = pd.read_csv(f'{DATA}/dwls_benchmark_results.csv')
bisque2     = pd.read_csv(f'{DATA}/bisque_benchmark_results_20k.csv')
bp2         = pd.read_csv(f'{DATA}/bayesprism_benchmark_results_20k.csv')
cibx_rt2    = pd.read_csv(f'{DATA}/cibersortx_runtime_20k.csv')
cibx_mem2   = pd.read_csv(f'{DATA}/cibersortx_memory_20k.csv')

rect_mem_row2 = rect_mem_s2[rect_mem_s2['n_cells'] == 20000].iloc[0]
rect_run_row2 = rect_run_s2[rect_run_s2['n_cells'] == 20000].iloc[0]
ag_mem_row2   = ag_mem_s2[ag_mem_s2['n_cells'] == 20000].iloc[0]
ag_run_row2   = ag_run_s2[ag_run_s2['n_cells'] == 20000].iloc[0]

scdc2_mem, scdc2_mem_sd = fix_mem(scdc2['Peak_Mem_MiB'].mean(), scdc2['Peak_Mem_MiB'].std())
dwls2_mem, dwls2_mem_sd = fix_mem(dwls2['Peak_Mem_MiB'].mean(), dwls2['Peak_Mem_MiB'].std())

rows_20k = [
    dict(Method='Rectangle',  Runtime_mean_s=round(rect_run_row2['time_mean'], 1),
         Runtime_SD_s=round(rect_run_row2['time_sd'], 1),
         Memory_mean_MiB=round(rect_mem_row2['memory_mean'], 0),
         Memory_SD_MiB=round(rect_mem_row2['memory_sd'], 0)),
    dict(Method='SCDC',       Runtime_mean_s=round(scdc2['Time_s'].mean(), 1),
         Runtime_SD_s=round(scdc2['Time_s'].std(), 1),
         Memory_mean_MiB=round(scdc2_mem, 0),
         Memory_SD_MiB=round(scdc2_mem_sd, 0)),
    dict(Method='Autogenes',  Runtime_mean_s=round(ag_run_row2['time_mean'], 1),
         Runtime_SD_s=round(ag_run_row2['time_sd'], 1),
         Memory_mean_MiB=round(ag_mem_row2['memory_mean'], 0),
         Memory_SD_MiB=round(ag_mem_row2['memory_sd'], 0)),
    dict(Method='MuSiC',      Runtime_mean_s=round(music2['Time_s'].mean(), 1),
         Runtime_SD_s=round(music2['Time_s'].std(), 1),
         Memory_mean_MiB=round(music2['Peak_Mem_MiB'].mean(), 0),
         Memory_SD_MiB=round(music2['Peak_Mem_MiB'].std(), 0)),
    dict(Method='Scaden',     Runtime_mean_s=round(scaden2['Time_s'].mean(), 1),
         Runtime_SD_s=round(scaden2['Time_s'].std(), 1),
         Memory_mean_MiB=round(scaden2['Peak_Mem_MiB'].mean(), 0),
         Memory_SD_MiB=round(scaden2['Peak_Mem_MiB'].std(), 0)),
    dict(Method='BayesPrism', Runtime_mean_s=round(bp2['Time_s'].mean(), 1),
         Runtime_SD_s=round(bp2['Time_s'].std(), 1),
         Memory_mean_MiB=round(bp2['Peak_Mem_MiB'].mean(), 0),
         Memory_SD_MiB=round(bp2['Peak_Mem_MiB'].std(), 0)),
    dict(Method='Bisque',     Runtime_mean_s=round(bisque2['Time_s'].mean(), 1),
         Runtime_SD_s=round(bisque2['Time_s'].std(), 1),
         Memory_mean_MiB=round(bisque2['Peak_Mem_MiB'].mean(), 0),
         Memory_SD_MiB=round(bisque2['Peak_Mem_MiB'].std(), 0)),
    dict(Method='DWLS',       Runtime_mean_s=round(dwls2['Time_s'].mean(), 1),
         Runtime_SD_s=round(dwls2['Time_s'].std(), 1),
         Memory_mean_MiB=round(dwls2_mem, 0),
         Memory_SD_MiB=round(dwls2_mem_sd, 0)),
    dict(Method='CIBERSORTx', Runtime_mean_s=round(cibx_rt2['time'].mean(), 1),
         Runtime_SD_s=round(cibx_rt2['time'].std(), 1),
         Memory_mean_MiB=round(cibx_mem2['memory_mib'].mean(), 0),
         Memory_SD_MiB=round(cibx_mem2['memory_mib'].std(), 0)),
]

df_20k = pd.DataFrame(rows_20k)
df_20k.insert(1, 'Ref size', '20,000')
df_20k.columns = ['Method', 'Ref size', 'Runtime mean (s)', 'Runtime SD (s)',
                  'Memory mean (MiB)', 'Memory SD (MiB)']

# ---------------------------------------------------------------------------
# Combined pivot table  (one row per method, Full + 20k side by side)
# ---------------------------------------------------------------------------
method_order = ['Rectangle', 'SCDC', 'Autogenes', 'MuSiC', 'Scaden',
                'BayesPrism', 'Bisque', 'DWLS', 'CIBERSORTx']

# Build lookup dicts keyed by Method
full_lookup = df_full.set_index('Method')
k20_lookup  = df_20k.set_index('Method')

def fmt(val, unit=''):
    if val == 'N/A' or (isinstance(val, float) and np.isnan(val)):
        return 'N/A'
    return f'{val:,.1f}{unit}' if isinstance(val, float) else str(val)

def fmt_sd(val):
    if val == 'N/A' or (isinstance(val, float) and np.isnan(val)):
        return ''
    return f' ± {val:,.1f}'

rows_pivot = []
for m in method_order:
    fr = full_lookup.loc[m]
    kr = k20_lookup.loc[m]
    note = str(fr.get('Note', '')).strip()
    rows_pivot.append({
        'Method':                         m,
        'Full ref size':                  fr['Ref size'],
        'Full Runtime mean ± SD (s)':     fmt(fr['Runtime mean (s)']) + fmt_sd(fr['Runtime SD (s)']),
        'Full Memory mean ± SD (MiB)':    fmt(fr['Memory mean (MiB)']) + fmt_sd(fr['Memory SD (MiB)']),
        '20k ref size':                   kr['Ref size'],
        '20k Runtime mean ± SD (s)':      fmt(kr['Runtime mean (s)']) + fmt_sd(kr['Runtime SD (s)']),
        '20k Memory mean ± SD (MiB)':     fmt(kr['Memory mean (MiB)']) + fmt_sd(kr['Memory SD (MiB)']),
        'Note':                           note,
    })

df_pivot = pd.DataFrame(rows_pivot)

sep = '-' * 150
print()
print('COMBINED TABLE  --  Runtime & Memory for all methods, full and 20k-subsampled reference')
print(sep)
print(df_pivot.to_string(index=False))
print()

# Save
out_dir = os.path.join(os.path.dirname(__file__), '..', 'output')
os.makedirs(out_dir, exist_ok=True)
df_pivot.to_csv(f'{out_dir}/benchmark_table_combined.csv', index=False)
df_full.to_csv(f'{out_dir}/benchmark_table_full_ref.csv', index=False)
df_20k.to_csv(f'{out_dir}/benchmark_table_20k_ref.csv', index=False)
print(f'CSVs saved to {os.path.abspath(out_dir)}/')

