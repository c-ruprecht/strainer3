import pandas as pd
import os
import argparse

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))



def main():
    # to do, basename + strain name addition
    parser = argparse.ArgumentParser(description='Using coverage and depth to make hit calls')
    parser.add_argument('--coverage', help='a global_counts.tsv.zst file', required=True)
    parser.add_argument('--cutoffs', help = 'Cutoff file path from strainer3',
                        default=os.path.join(SCRIPT_DIR,'..', 'data', 'cutoffs.tsv'))
    parser.add_argument('--output_dir', help='directory where output is saved')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    # concat hit files for different strains
    df_cutoffs = pd.read_csv(args.cutoffs, sep = '\t')
    print(df_cutoffs)

    # call hit with PR values
    df_cov = pd.read_csv(args.coverage, sep = '\t')
    print(df_cov)

if __name__ == "__main__":
        main()