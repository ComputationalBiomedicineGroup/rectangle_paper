"""Evaluate Figure 2's current saved predictions, without rerunning deconvolution."""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr


def evaluate_pairs(truth, predictions, method, dataset):
    """Align by both labels before calculating pooled and per-cell metrics."""
    if not truth.index.is_unique or not truth.columns.is_unique:
        raise ValueError(f'{dataset}: duplicate ground-truth labels')
    if not predictions.index.is_unique or not predictions.columns.is_unique:
        raise ValueError(f'{dataset}/{method}: duplicate prediction labels')
    cells = sorted(set(truth.index) & set(predictions.index))
    samples = sorted(set(truth.columns) & set(predictions.columns))
    if not cells or not samples:
        raise ValueError(f'{dataset}/{method}: no matching cell types or samples')
    actual = truth.loc[cells, samples].to_numpy(dtype=float)
    estimated = predictions.loc[cells, samples].to_numpy(dtype=float)
    if not np.isfinite(actual).all() or not np.isfinite(estimated).all():
        raise ValueError(f'{dataset}/{method}: non-finite evaluation values')
    errors = estimated - actual
    per_cell_rmse = np.sqrt(np.mean(errors ** 2, axis=1))
    overall_rmse = float(np.sqrt(np.mean(errors ** 2)))
    if not np.isclose(overall_rmse, np.sqrt(np.mean(per_cell_rmse ** 2))):
        raise ValueError(f'{dataset}/{method}: pooled and per-cell RMSE disagree')
    rename = {'T cells CD4 Total': 'CD4 T cells', 'T cells CD4 conv': 'CD4 T cells',
              'T cells CD4': 'CD4 T cells', 'T cells CD8': 'CD8 T cells',
              'Neutrophilss': 'Neutrophils'}
    rows = []
    for index, cell in enumerate(cells):
        cor = (float(pearsonr(actual[index], estimated[index]).statistic)
               if np.ptp(actual[index]) and np.ptp(estimated[index]) else np.nan)
        rows.append({'cell_type': rename.get(cell, cell), 'method': method,
                     'bulk_ds': dataset, 'cor': cor, 'RMSE': per_cell_rmse[index]})
    overall = {'cell_type': 'all', 'method': method, 'bulk_ds': dataset,
               'cor': float(pearsonr(actual.ravel(), estimated.ravel()).statistic),
               'RMSE': overall_rmse, 'n_samples': len(samples), 'n_cell_types': len(cells),
               'n_truth_samples': len(truth.columns), 'n_prediction_samples': len(predictions.columns)}
    return rows, overall


RECTANGLE_PREDICTIONS = {
    'finotello': 'output/hao_estimations_fino.csv',
    'chen': 'output/muris_estimations_chen.csv',
    'petitprez': 'output/muris_estimations_peti.csv',
    'hoek': 'output/hao_estimations_hoek.csv',
    'morandini': 'output/hao_estimations_morandini.csv',
    'dietrich': 'output/hao_estimations_dietrich.csv',
}

# Written by scripts/run_secret_fig2.R (HaoSub reference for the human, Tabula Muris for the mouse bulks)
SECRET_PREDICTIONS = {dataset: f'output/secret_fig2/secret_estimations_{dataset}.csv'
                      for dataset in RECTANGLE_PREDICTIONS}
# SECRET without the unknown component (fractions sum to 1)
SECRET_NO_UNKNOWN_PREDICTIONS = {dataset: f'output/secret_fig2/secret_estimations_{dataset}_no_unknown.csv'
                                 for dataset in RECTANGLE_PREDICTIONS}


def current_rectangle_metrics(project_root):
    return current_method_metrics(project_root, 'rectangle', RECTANGLE_PREDICTIONS)


def current_method_metrics(project_root, method, prediction_paths):
    """Score one method's saved predictions (samples x cell types) the way Rectangle's are scored.

    Datasets missing from prediction_paths are skipped."""
    rows, overall = [], []
    for dataset in TRUTH:
        prediction_path = prediction_paths.get(dataset)
        if prediction_path is None:
            continue
        truth, predictions = load_truth_and_predictions(project_root, dataset, prediction_path)
        cell_rows, total = evaluate_pairs(truth, predictions, method, dataset)
        rows.extend(cell_rows)
        total['prediction_file'] = prediction_path
        overall.append(total)
    return pd.DataFrame(rows), pd.DataFrame(overall)


# Ground truth per dataset: path and, where the file has no usable row names, the cell-type labels
TRUTH = {
    'finotello': ('data/finotello_facs.csv',
                  ['NK cells', 'B cells', 'Tregs', 'mDC', 'Monocytes', 'Neutrophilss',
                   'T cells CD8', 'T cells CD4 conv', 'Other']),
    'chen': ('data/chen_facs.csv', ['Monocytes', 'B cells', 'T cells CD8', 'T cells CD4 Total']),
    'petitprez': ('data/petitprez_facs.csv',
                  ['T cells', 'T cells CD8', 'B cells', 'Tregs', 'NK cells', 'B cells ger',
                   'Basophils', 'Mast', 'Macrophages', 'Eosinophils', 'Neutrophils',
                   'Monocytes', 'B derived', 'Mono Macro', 'NK_or_T', 'Granulocytes']),
    'hoek': ('data/hoek_facs.csv', ['T cells', 'Monocytes', 'B cells', 'mDC', 'NK cells']),
    'morandini': ('data/fig2/morandini/morandini_facs.csv', None),
    'dietrich': ('data/fig2/deconvMe/dietrich_facs.csv', None),
}


def load_truth(project_root, dataset):
    """Ground truth (cell types x samples) with the labels the predictions are matched on."""
    truth_path, cell_labels = TRUTH[dataset]
    truth = pd.read_csv(Path(project_root) / truth_path, index_col=0)
    if cell_labels is not None:
        truth.index = cell_labels
        truth = truth.fillna(0)
    if dataset in ('petitprez', 'dietrich'):
        truth = truth.drop(index='T cells', errors='ignore')
    return truth


def load_truth_and_predictions(project_root, dataset, prediction_path):
    """Ground truth and one method's saved predictions (samples x cell types in the file), both as
    cell types x samples, with the combined cell types each ground truth reports."""
    truth = load_truth(project_root, dataset)
    predictions = pd.read_csv(Path(project_root) / prediction_path, index_col=0)
    return truth, combine_predictions(dataset, predictions)


def combine_predictions(dataset, predictions):
    """Adds the combined cell types a dataset's ground truth reports to predictions (samples x cell types) and
    returns them as cell types x samples."""
    predictions = predictions.copy()
    if dataset == 'chen':
        predictions['T cells CD4 Total'] = predictions['T cells CD4 conv'] + predictions['Tregs']
    if dataset in ('hoek', 'morandini'):
        predictions['T cells'] = predictions[['T cells CD4 conv', 'T cells CD8', 'Tregs']].sum(axis=1)
    if dataset == 'morandini':
        predictions['T cells CD4'] = predictions['T cells CD4 conv']
        predictions['Lymphocytes'] = predictions[['T cells', 'NK cells', 'ILC', 'B cells', 'Plasma cells']].sum(axis=1)
    return predictions.T


def dietrich_comparator_predictions(project_root):
    """The other methods' Dietrich predictions (deconvMe runs) as {method: cell types x samples}."""
    predictions = {}
    for path in sorted((Path(project_root) / 'data/fig2/deconvMe').glob('*/deconvolution.csv')):
        method = path.parent.name.split('_')[0].lower()
        if method == 'rectangle':
            continue
        frame = pd.read_csv(path, index_col=0)
        frame.index = frame.index.astype(str).str.replace('^X', '', regex=True).str.replace('.', '-', regex=False)
        predictions[method] = frame.T
    return predictions


# omnideconv estimates of the eight comparators on the real bulks (unpacked from data/fig2/real_bulks_27_10_2026.tar.gz):
# <reference>/<method>_<reference>_<sc norm>_<dataset>_<bulk norm>_ct0_rep0/deconvolution.rds
REAL_BULKS = 'data/fig2/real_bulks'
OMNIDECONV_METHODS = ['autogenes', 'bayesprism', 'bisque', 'cibersortx', 'dwls', 'music', 'scaden', 'scdc']


def omnideconv_predictions(project_root, method, dataset):
    """A comparator's estimates as cell types x samples (with the combined cell types of the ground truth),
    or None when the run is missing."""
    import pyreadr

    paths = list((Path(project_root) / REAL_BULKS).glob(f'*/{method}_*_{dataset}_*_ct0_rep0/deconvolution.rds'))
    if len(paths) != 1:
        return None
    estimates = next(iter(pyreadr.read_r(paths[0]).values()))  # samples x cell types
    estimates.index = estimates.index.astype(str)
    if not set(estimates.index) & set(load_truth(project_root, dataset).columns):
        # CIBERSORTx returns sample names changed by R's make.names (95-CO350-EOT -> X95.CO350.EOT)
        estimates.index = estimates.index.str.replace('^X', '', regex=True).str.replace('.', '-', regex=False)
    return combine_predictions(dataset, estimates)


def omnideconv_comparator_metrics(project_root, dataset):
    """Scores the comparators from their estimates exactly like Rectangle (the omnideconv metrics of Chen compared
    CD4 conv only with the CD4 total ground truth, Rectangle CD4 conv + Tregs)."""
    truth = load_truth(project_root, dataset)
    rows, overall = [], []
    for method in OMNIDECONV_METHODS:
        predictions = omnideconv_predictions(project_root, method, dataset)
        if predictions is None:
            raise FileNotFoundError(f'No {method} estimates for {dataset} in {REAL_BULKS}')
        cell_rows, total = evaluate_pairs(truth, predictions, method, dataset)
        rows.extend(cell_rows)
        overall.append(total)
    return pd.DataFrame(rows), pd.DataFrame(overall)


def current_dietrich_comparator_metrics(project_root):
    truth = load_truth(project_root, 'dietrich')
    rows, overall = [], []
    for method, predictions in dietrich_comparator_predictions(project_root).items():
        cell_rows, total = evaluate_pairs(truth, predictions, method, 'dietrich')
        rows.extend(cell_rows)
        overall.append(total)
    return pd.DataFrame(rows), pd.DataFrame(overall)


def morandini_comparator_metrics(project_root):
    results = Path(project_root) / 'data/fig2/results_main_morandini_ncbi'
    frames = []
    for path in sorted(results.glob('*/cor_cell_type.csv')):
        cor = pd.read_csv(path)[['cell_type', 'cor']]
        rmse = pd.read_csv(path.parent / 'rmse_cell_type.csv')[['cell_type', 'RMSE']]
        frame = cor.merge(rmse, on='cell_type', validate='one_to_one')
        frame['method'] = path.parent.name.split('_')[0].lower()
        frame['bulk_ds'] = 'morandini'
        frame['cell_type'] = frame['cell_type'].replace({
            'T cells CD4': 'CD4 T cells', 'T cells CD8': 'CD8 T cells',
            'Neutrophilss': 'Neutrophils',
        })
        frames.append(frame[['cell_type', 'method', 'bulk_ds', 'cor', 'RMSE']])
    if not frames:
        raise FileNotFoundError(f'No Morandini comparator metrics found in {results}')
    return pd.concat(frames, ignore_index=True)
