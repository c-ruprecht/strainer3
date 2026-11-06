import pandas as pd
import os
import argparse
import glob
import json
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


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


def main():
    # to do, basename + strain name addition
    parser = argparse.ArgumentParser(description='Using coverage and depth to make hit calls')
    parser.add_argument('--coverage_file', help='a coverage_depth.tsv file')
    parser.add_argument('--hist_file', help='a zi_hist.tsv.gz file (coverage_depth_locations.py --export_hist)')
    parser.add_argument('--coverage_dir', help = "Directory in which coverage files are contained, will discover 'coverage_depth.tsv' and 'zi_hist.tsv.gz' files in subdirectories and concat them")
    parser.add_argument('--cutoffs', help = 'Cutoff file path from strainer3', default=os.path.join(SCRIPT_DIR,'..', 'data', 'cutoffs.tsv'))
    parser.add_argument('--output_dir', help='directory where output is saved')
    parser.add_argument('--verbose',action='store_true', help = 'incluces all possible metrics from coverage_depth into hit call table')
    args = parser.parse_args()

    if not (args.coverage_file or args.coverage_dir):
        parser.error('give --coverage_file or --coverage_dir')

    os.makedirs(args.output_dir, exist_ok=True)
    # concat hit files for different strains
    df_cutoffs = pd.read_csv(args.cutoffs, sep = '\t')
    print(df_cutoffs)

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


    # get rare cutoffs
    df_cutoffs_s = df_cutoffs[['rare kmers', 'PR']].sort_values('rare kmers')
    df_cov = df_cov.sort_values('coverage_rare_kmers')
    df_cov = pd.merge_asof(df_cov, df_cutoffs_s,
                           left_on='coverage_rare_kmers', right_on='rare kmers',
                           direction='backward')

    df_cov['hit_call_pr'] = df_cov['PR'].astype(str) + ' (' + df_cov['rare kmers'].astype(str) + ')'
    df_cov.loc[df_cov['PR'].isna(), 'hit_call_pr'] = np.nan
    df_cov = df_cov.drop(columns=['PR', 'rare kmers'])

    df_cov['depth_ratio'] = df_cov['depth_rare_kmers']/df_cov['depth_nonrare_kmers']

    # get all metrics neccessary for xgb hit calling
    df_hits = df_cov[[ 'strain','sample', 'hit_call_pr', 'coverage_rare_kmers', 'rare_count']].copy()

    # make the hit call

    print(df_hits)
    if args.verbose:
        df_hits = df_cov.copy()
        print(df_hits)
    df_hits.sort_values('coverage_rare_kmers', ascending= False).to_csv(os.path.join(args.output_dir, 'hit_calls.tsv'), sep = '\t', index = None)


if __name__ == "__main__":
        main()