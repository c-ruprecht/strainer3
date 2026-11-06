import pandas as pd
import os
import argparse
import glob
import json
import numpy as np
import xgboost as xgb

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, '..', 'data')


def concat_files(paths, label):
    """read and concat tsv files that carry a 'strain' column"""
    li_dfs = []
    for df_path in paths:
        df = pd.read_csv(df_path, sep='\t')
        if 'strain' in df.columns:
            li_dfs.append(df)
        else:
            print(f"Add 'strain' to {df_path} before concatenating")
    if not li_dfs:
        raise ValueError(f"Found {label} files, but none have a 'strain' column")
    return pd.concat(li_dfs, ignore_index=True)


def load_model(prefix):
    """load the xgb model (<prefix>.json) and its metadata (<prefix>_meta.json) exported by the notebook"""
    with open(f'{prefix}_meta.json') as fh:
        meta = json.load(fh)
    model = xgb.XGBClassifier()
    model.load_model(f'{prefix}.json')
    return model, meta


def xgb_hit_call(df, model, meta, cutoff_rule=None):
    """add xgb_prob (P(present)) and xgb_call (0/1) to df using the features stored in meta"""
    features = meta['features']
    missing = [f for f in features if f not in df.columns]
    if missing:
        raise KeyError(f'model needs columns {missing}, not found in coverage table')

    if cutoff_rule:
        if cutoff_rule not in meta['cutoffs']:
            raise KeyError(f"cutoff rule '{cutoff_rule}' not in model, choose from {list(meta['cutoffs'])}")
        threshold = meta['cutoffs'][cutoff_rule]['threshold']
    else:
        cutoff_rule, threshold = meta['cutoff_rule'], meta['threshold']

    # same preparation as in training: inf -> +/-1e12, NaN kept (xgboost routes missing values)
    X = df[features].replace([np.inf, -np.inf], [1e12, -1e12]).astype(float)
    df = df.copy()
    df['xgb_prob'] = model.predict_proba(X)[:, 1]
    df['xgb_call'] = (df['xgb_prob'] >= threshold).astype(int)
    print(f"xgb hit call: features {features}, rule '{cutoff_rule}', threshold {threshold:.3f} "
          f"-> {df['xgb_call'].sum()} of {len(df)} strain x sample calls")
    return df


def main():
    # to do, basename + strain name addition
    parser = argparse.ArgumentParser(description='Using coverage and depth to make hit calls')
    parser.add_argument('--coverage_file', help='a coverage_depth.tsv file')
    parser.add_argument('--hist_file', help='a zi_hist.tsv.gz file (coverage_depth_locations.py --export_hist)')
    parser.add_argument('--coverage_dir', help = "Directory in which coverage files are contained, will discover 'coverage_depth.tsv' and 'zi_hist.tsv.gz' files in subdirectories and concat them")
    parser.add_argument('--model', help = 'xgb model prefix: <prefix>.json + <prefix>_meta.json',
                        default=os.path.join(DATA_DIR, 'xgb_hitcall_2param'))
    parser.add_argument('--cutoff_rule', default='precision >= 0.95',
                        help = "probability cutoff stored in the model meta, e.g. 'max F1' or 'precision >= 0.99' "
                        "(default: 'precision >= 0.95')")
    parser.add_argument('--output_dir', help='directory where output is saved')
    parser.add_argument('--verbose',action='store_true', help = 'incluces all possible metrics from coverage_depth into hit call table')
    args = parser.parse_args()

    if not (args.coverage_file or args.coverage_dir):
        parser.error('give --coverage_file or --coverage_dir')

    os.makedirs(args.output_dir, exist_ok=True)

    df_cov = None
    df_hist = None

    if args.coverage_file:
        df_cov = pd.read_csv(args.coverage_file, sep='\t')

    if args.hist_file:
        df_hist = pd.read_csv(args.hist_file, sep='\t')

    if args.coverage_dir:
        coverage_dir = args.coverage_dir
        print(f"searching {coverage_dir} (exists: {os.path.isdir(coverage_dir)})")

        glob_cov = glob.glob(os.path.join(coverage_dir, '**', '*coverage_depth.tsv'), recursive=True)
        print(glob_cov)
        if not glob_cov:
            raise FileNotFoundError(f"No *coverage_depth.tsv found under {coverage_dir}")
        df_cov = concat_files(glob_cov, 'coverage')

        # zero-inflation histograms, only present if coverage_depth_locations.py ran with --export_hist
        glob_hist = glob.glob(os.path.join(coverage_dir, '**', '*zi_hist.tsv.gz'), recursive=True)
        print(glob_hist)
        if glob_hist:
            df_hist = concat_files(glob_hist, 'histogram')

            # fitting settings per run: histograms can be refit with any settings, but the exported metrics
            # in coverage_depth.tsv are only comparable if all runs used the same ones
            glob_params = glob.glob(os.path.join(coverage_dir, '**', 'zi_params.json'), recursive=True)
            params = {}
            for p in glob_params:
                with open(p) as fh:
                    params[p] = json.dumps(json.load(fh), sort_keys=True)
            if len(set(params.values())) > 1:
                print('WARNING: runs used different zero-inflation settings:')
                for p, v in params.items():
                    print(f'  {p}: {v}')

    if df_hist is not None:
        hist_path = os.path.join(args.output_dir, 'zi_hist.tsv.gz')
        df_hist.to_csv(hist_path, sep='\t', index=None)
        print(f"Saved {df_hist[['strain', 'sample']].drop_duplicates().shape[0]} strain x sample histograms to \n {hist_path}")


    df_cov['depth_ratio'] = df_cov['depth_rare_kmers']/df_cov['depth_nonrare_kmers']

    # make the hit call with the xgb model
    model, meta = load_model(args.model)
    df_cov = xgb_hit_call(df_cov, model, meta, args.cutoff_rule)

    # non verbose: the model's input columns, its probability and the call
    df_hits = df_cov[['strain', 'sample'] + meta['features'] + ['xgb_prob', 'xgb_call']].copy()
    if args.verbose:
        df_hits = df_cov.copy()
    print(df_hits)
    df_hits.sort_values('xgb_prob', ascending= False).to_csv(os.path.join(args.output_dir, 'hit_calls.tsv'), sep = '\t', index = None)


if __name__ == "__main__":
        main()