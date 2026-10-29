import pandas as pd
import os
import argparse
import glob
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))



def main():
    # to do, basename + strain name addition
    parser = argparse.ArgumentParser(description='Using coverage and depth to make hit calls')
    parser.add_argument('--coverage_file', help='a global_counts.tsv.zst file')
    parser.add_argument('--coverage_dir', help = "Directory in which doverage files are contained, will discover 'coverage_depth.tsv' files in subdirectories and concat them")
    parser.add_argument('--cutoffs', help = 'Cutoff file path from strainer3', default=os.path.join(SCRIPT_DIR,'..', 'data', 'cutoffs.tsv'))
    parser.add_argument('--output_dir', help='directory where output is saved')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    # concat hit files for different strains
    df_cutoffs = pd.read_csv(args.cutoffs, sep = '\t')
    print(df_cutoffs)

    df_cov = None

    if args.coverage_file:
        df_cov = pd.read_csv(args.coverage_file, sep='\t')

    if args.coverage_dir:
        coverage_dir = args.coverage_dir
        print(f"searching {coverage_dir} (exists: {os.path.isdir(coverage_dir)})")

        glob_cov = glob.glob(os.path.join(coverage_dir, '**', '*coverage_depth.tsv'), recursive=True)
        print(glob_cov)
        if not glob_cov:
            raise FileNotFoundError(f"No *coverage_depth.tsv found under {coverage_dir}")

        li_dfs = []
        for df_path in glob_cov:
            df = pd.read_csv(df_path, sep='\t')
            if 'strain' in df.columns:
                li_dfs.append(df)
            else:
                print(f"Add 'strain' to {df_path} before concatenating")

        if not li_dfs:
            raise ValueError("Found coverage files, but none have a 'strain' column")
        df_cov = pd.concat(li_dfs, ignore_index=True)

    
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

    df_hits = df_cov[[ 'strain','sample', 'hit_call_pr', 'coverage_rare_kmers',  'pair_cov_geo', 'depth_ratio', 'rare_count']].copy()
    print(df_hits)
    df_hits.sort_values('coverage_rare_kmers', ascending= False).to_csv(os.path.join(args.output_dir, 'hit_calls.tsv'), sep = '\t', index = None)


    #concat only if strain column exists
if __name__ == "__main__":
        main()