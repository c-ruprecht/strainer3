import pandas as pd
import os
import argparse
import numpy as np

def get_sample_coverage(df_hits, df_counts):

    df_counts = df_counts.drop(columns=['reference_count', 'drug_count']).set_index('#kmer').sum(axis=1)
    df_counts = df_counts.reset_index(name= 'count')
    df = df_hits.set_index('#kmer')

    # get full kmer coverages based on columns
    df_cov_depth = pd.DataFrame({'all_kmers_hit': (df > 0).sum(axis=0),
                                'all_kmers_total': len(df),
                                'coverage_all_kmers': round((df > 0).sum(axis=0)/len(df),4),
                                'depth_all_kmers': round(df.mean(axis=0),2),
                            })


    df_counts_uni = df_counts.loc[df_counts['count']==0].copy()
    rare_count = '==0'
    if len(df_counts_uni) < 1000:
        df_counts_uni = df_counts.loc[df_counts['count'] <= 1].copy()
        rare_count = '<=1'
    df_counts_uni['rare_count'] = rare_count

    df_unique = df_hits.loc[df_hits['#kmer'].isin(df_counts_uni['#kmer'])].set_index('#kmer')

    df_cov_depth_uni = pd.DataFrame({'rare_kmers_hit': (df_unique > 0).sum(axis=0),
                                'rare_kmers_total': len(df_unique),
                                'coverage_rare_kmers': round((df_unique > 0).sum(axis=0)/len(df_unique),4),
                                'depth_rare_kmers': round(df_unique.mean(axis=0),2),
                                'rare_count': rare_count,
                            })

    df_nonunique = df_hits.loc[df_hits['#kmer'].isin(df_counts_uni['#kmer'])==False].set_index('#kmer')
    
    df_cov_depth_nonuni = pd.DataFrame({'nonrare_kmers_hit': (df_nonunique > 0).sum(axis=0),
                                'nonrare_kmers_total': len(df_nonunique),
                                'coverage_nonrare_kmers': round((df_nonunique > 0).sum(axis=0)/len(df_nonunique),4),
                                'depth_nonrare_kmers': round(df_nonunique.mean(axis=0),2),
                            })
    df_cov_depth = pd.concat([df_cov_depth, df_cov_depth_uni, df_cov_depth_nonuni], axis = 1)    

    return df_cov_depth


def main():
    # to do, basename + strain name addition
    parser = argparse.ArgumentParser(description='Calcualting coverage and depth with kmer locations')
    parser.add_argument('--counts', help='a global_counts.tsv.zst file', required=True)
    parser.add_argument('--hits', help='a target_strain .kmer_hits.tsv.gz', required = True)
    parser.add_argument('--output_dir', help='directory where output is saved')
    parser.add_argument('--strain', help = 'add a name for the target strain you are tracking')
   
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    hits_path = args.hits
    counts_path = args.counts
    outpath = args.output_dir


    df_hits = pd.read_csv(hits_path, sep = '\t')
    total_evaluated = df_hits.loc[df_hits['#kmer']=='total_evaluated'].drop(columns='#kmer').T.reset_index()#.to_dict()
    total_evaluated.columns = ['sample', 'total_evaluated']
    print(total_evaluated)

    df_hits = df_hits.loc[df_hits['#kmer'] !='total_evaluated'].copy()
    df_counts = pd.read_csv(counts_path, sep = '\t')

    # get sample coverage for all kemrs and unique/rare kmers
    df_cov_depth = get_sample_coverage(df_hits=df_hits, df_counts= df_counts)

    df_cov_depth = df_cov_depth.reset_index(names = 'sample')
    
    df_cov_depth = df_cov_depth.merge(total_evaluated, on='sample', how='left')
    if args.strain:
        df_cov_depth['strain'] = args.strain 

    df_cov_depth.set_index(['strain', 'sample']).to_csv(os.path.join(outpath,"coverage_depth.tsv"), sep = '\t')
    print(f'Saved coverage dataframe to \n {outpath}')
    print(df_cov_depth.head(5))



if __name__ == '__main__':
    main()