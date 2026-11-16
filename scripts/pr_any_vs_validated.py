#!/usr/bin/env python3
"""Does read rejection help? PR comparison of the validated and the non-rejected table of ONE run.

Both tables come from the same reads and the same k-mers (kmer_block_detect writes kmer_hits and kmer_hits_any
in one pass; block_coverage_depth.py --kmer_hits_any turns them into two tables with the same column names), so
a difference is the effect of read rejection alone.

  pr_any_vs_validated.py --validated coverage_depth_all.tsv --any coverage_depth_any_all.tsv --truth truth.tsv

truth: tab-separated, columns strain, sample, truth_bin (1 = strain present); optional truth_label.
For every metric: average precision and recall at fixed precision for both tables, and a bootstrap over groups
(default column genome_list = donor; falls back to strain) for the difference validated - any.
A zi metric is used as "coverage set to 0 when zi_frac > cut" (NaN zi keeps the coverage), as in the SC2 analysis.
Scores that are NaN count as 0 (no evidence); rows are never dropped.
"""
import argparse
import sys

import numpy as np
import pandas as pd


def pr_points(y, s):
    """precision and recall at every distinct score threshold (calls = score >= threshold), highest first"""
    order = np.argsort(-s, kind='mergesort')
    y, s = y[order], s[order]
    last = np.r_[np.flatnonzero(np.diff(s) != 0), len(s) - 1]      # last row of each tied score
    tp = np.cumsum(y)[last]
    n_called = last + 1
    return tp / n_called, tp / max(int(y.sum()), 1)


def average_precision_score(y, s):
    prec, rec = pr_points(y, s)
    return float(np.sum(np.diff(np.r_[0.0, rec]) * prec))


def recall_at(y, s, p):
    prec, rec = pr_points(y, s)
    ok = prec >= p
    return float(rec[ok].max()) if ok.any() else 0.0


def cv_scores(df, feats, y, groups, n_splits=6, seed=0):
    """donor-held-out predictions of a small XGBoost model: every row is scored by a model that never saw its
    group. Returns None if xgboost is not installed."""
    try:
        import xgboost as xgb
    except ImportError:
        return None
    X = df[feats].astype(float).replace([np.inf, -np.inf], np.nan).to_numpy()
    uniq = np.unique(groups)
    fold_of = {g: i % min(n_splits, len(uniq)) for i, g in enumerate(np.random.default_rng(seed).permutation(uniq))}
    fold = np.array([fold_of[g] for g in groups])
    p = np.zeros(len(df))
    for f in np.unique(fold):
        tr, te = fold != f, fold == f
        if y[tr].sum() == 0 or y[tr].sum() == tr.sum():
            continue
        m = xgb.XGBClassifier(n_estimators=300, max_depth=3, learning_rate=0.05, subsample=0.8,
                              colsample_bytree=0.8, random_state=seed, eval_metric='logloss')
        m.fit(X[tr], y[tr])
        p[te] = m.predict_proba(X[te])[:, 1]
    return p


def scores(df, metric, zi_cut):
    if metric.endswith('_zi_pruned'):
        base = metric[:-len('_zi_pruned')]
        cov, zi = df[f'{base}_coverage'], df[f'{base}_zi_frac']
        return cov.where(~(zi > zi_cut), 0.0).fillna(0.0).to_numpy(float)
    return df[metric].fillna(0.0).to_numpy(float)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--validated', required=True, help='coverage_depth_all.tsv (with read rejection)')
    ap.add_argument('--any', required=True, help='coverage_depth_any_all.tsv (without read rejection)')
    ap.add_argument('--truth', required=True, help='strain, sample, truth_bin [, truth_label]')
    ap.add_argument('--metrics', default='kmer_all_coverage,kmer_all_zi_pruned,kmer_snp_coverage,kmer_snp_zi_pruned',
                    help='comma-separated; <subset>_zi_pruned = coverage zeroed when zi_frac > --zi_cut')
    ap.add_argument('--zi_cut', type=float, default=0.6)
    ap.add_argument('--model_features', default='kmer_all_coverage,kmer_all_zi_frac',
                    help='features of the coverage + zero-inflation hit-call model (XGBoost, cross-validated with '
                         'whole groups held out, same folds for both tables); empty = skip. Needs xgboost.')
    ap.add_argument('--group', default='genome_list', help='bootstrap unit (donor); falls back to strain')
    ap.add_argument('--bootstrap', type=int, default=1000)
    ap.add_argument('--precisions', default='1.0,0.99,0.95,0.90')
    ap.add_argument('--output', help='write the summary table here (tsv)')
    a = ap.parse_args()

    key = ['strain', 'sample']
    truth = pd.read_csv(a.truth, sep='\t', dtype={'sample': str})
    val = pd.read_csv(a.validated, sep='\t', dtype={'sample': str}, low_memory=False)
    anyt = pd.read_csv(a.any, sep='\t', dtype={'sample': str}, low_memory=False)
    if set(map(tuple, val[key].to_numpy())) != set(map(tuple, anyt[key].to_numpy())):
        sys.exit('ERROR: the two tables do not hold the same strain x sample rows')
    val = val.merge(truth, on=key, how='inner').sort_values(key).reset_index(drop=True)
    anyt = anyt.merge(truth, on=key, how='inner').sort_values(key).reset_index(drop=True)
    if len(val) == 0:
        sys.exit('ERROR: no strain x sample row of the tables is in the truth table')
    y = val['truth_bin'].to_numpy(int)
    grp_col = a.group if a.group in val.columns else ('genome_list' if 'genome_list' in truth.columns else 'strain')
    groups = val[grp_col].astype(str).to_numpy()
    uniq = np.unique(groups)
    idx_by = {g: np.flatnonzero(groups == g) for g in uniq}
    precs = [float(p) for p in a.precisions.split(',')]
    print(f'{len(val)} rows, {int(y.sum())} positives, bootstrap over {len(uniq)} groups of "{grp_col}"')

    rng = np.random.default_rng(1)
    rows = []
    todo = [(m, None) for m in a.metrics.split(',') if m]
    feats = [f for f in a.model_features.split(',') if f]
    if feats:
        if any(f not in val.columns or f not in anyt.columns for f in feats):
            print(f'skipping the model: a feature of {feats} is missing in one of the tables')
        elif len(uniq) < 3:
            print('skipping the model: fewer than 3 groups to hold out')
        else:
            pv, pa = cv_scores(val, feats, y, groups), cv_scores(anyt, feats, y, groups)
            if pv is None:
                print('skipping the model: xgboost is not installed')
            else:
                todo.append(('xgb_cv(' + '+'.join(feats) + ')', (pv, pa)))
    for metric, pre in todo:
        if pre is not None:
            sv, sa = pre
        else:
            need = ([metric[:-10] + '_coverage', metric[:-10] + '_zi_frac'] if metric.endswith('_zi_pruned') else [metric])
            if any(c not in val.columns or c not in anyt.columns for c in need):
                print(f'skipping {metric}: column missing in one of the tables')
                continue
            sv, sa = scores(val, metric, a.zi_cut), scores(anyt, metric, a.zi_cut)
        stats = {'AP': lambda yy, ss: average_precision_score(yy, ss)}
        stats.update({f'R@P{int(round(p * 100))}': (lambda yy, ss, p=p: recall_at(yy, ss, p)) for p in precs})
        boots = {n: [] for n in stats}
        for _ in range(a.bootstrap):
            ii = np.concatenate([idx_by[g] for g in rng.choice(uniq, len(uniq))])
            if y[ii].sum() == 0 or y[ii].sum() == len(ii):
                continue
            for n, f in stats.items():
                boots[n].append(f(y[ii], sv[ii]) - f(y[ii], sa[ii]))
        for n, f in stats.items():
            b = np.array(boots[n])
            rows.append(dict(metric=metric, stat=n, any=f(y, sa), validated=f(y, sv), diff=f(y, sv) - f(y, sa),
                             ci_lo=np.quantile(b, 0.025) if len(b) else np.nan,
                             ci_hi=np.quantile(b, 0.975) if len(b) else np.nan))
    out = pd.DataFrame(rows)
    pd.set_option('display.width', 200)
    print(out.round(4).to_string(index=False))
    if a.output:
        out.to_csv(a.output, sep='\t', index=False)


if __name__ == '__main__':
    main()
