import pandas as pd
import os
import argparse
import numpy as np
import ahocorasick
from Bio import SeqIO
from Bio.Seq import Seq
import gzip

def get_sample_coverage(df_hits, df_counts, rare_t = 1000):

    df_counts = df_counts.drop(columns=['reference_count', 'drug_count']).set_index('#kmer').sum(axis=1)
    df_counts = df_counts.reset_index(name= 'count')
    df = df_hits.set_index('#kmer')

    # get full kmer coverages based on columns
    df_cov_depth = pd.DataFrame({'all_kmers_hit': (df > 0).sum(axis=0),
                                'all_kmers_total': len(df),
                                'coverage_all_kmers': round((df > 0).sum(axis=0)/len(df),4),
                                'depth_all_kmers': round(df.mean(axis=0),2),
                            })


    max_count = int(df_counts['count'].max())
    for t in range(0, max_count + 1):
        df_counts_uni = df_counts.loc[df_counts['count'] <= t].copy()
        rare_count = '==0' if t == 0 else f'<={t}'
        if len(df_counts_uni) >= rare_t:
            break
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

    return df_cov_depth, df_counts_uni

def load_genome(genome_path):
    opener = gzip.open if genome_path.endswith('.gz') else open
    records = {}
    with opener(genome_path, 'rt') as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            records[record.id] = record.seq
    return records

def build_mapped_kmers_ahocorasick(records, kmers, terminal_dist):
    # Build Aho-Corasick automaton with forward and reverse complement kmers
    # Important: only pass single-count kmers; only the first hit per kmer is kept
    kmers = list(kmers)
    A = ahocorasick.Automaton()
    for kmer in kmers:
        A.add_word(kmer, (kmer, False))
        A.add_word(str(Seq(kmer).reverse_complement()), (kmer, True))
    A.make_automaton()

    found = set()
    rows = []
    for record_id, seq in records.items():
        for pos, (kmer, is_rc) in A.iter(str(seq)):
            if kmer not in found:
                rows.append((record_id, kmer, pos - len(kmer) + 1, is_rc))
                found.add(kmer)

    df = pd.DataFrame(rows, columns=['contig_id', '#kmer', 'kmer_position', 'reverse_complement'])

    if len(df) < len(kmers):
        print(f'WARNING: not all kmers found in genome ({len(df)} / {len(kmers)})')
    elif len(df) > len(kmers):
        print('WARNING: kmers found more than once')
    else:
        print(f'  {len(df)} kmers mapped (all unique)')

    dict_len = {cid: len(seq) for cid, seq in records.items()}
    df['contig_length'] = df['contig_id'].map(dict_len)
    df['terminal_kmer'] = (
        (df['kmer_position'] < terminal_dist) |
        ((df['contig_length'] - df['kmer_position']) < terminal_dist)
    )
    df['label'] = df['terminal_kmer'].map({True: 'terminal', False: 'internal'})

    n_terminal = int(df['terminal_kmer'].sum())
    print(f'  Terminal kmers: {n_terminal} / {len(df)}')

    return df, dict_len

# get kmer blocks
def add_block_ids(df, min_block_size=30):
    df = df.sort_values(['contig_id', 'kmer_position']).copy()
    pos_break = df.groupby('contig_id')['kmer_position'].diff().ne(1)
    df['block_id'] = (pos_break).cumsum()
    df['block_size'] = df.groupby('block_id')['block_id'].transform('size')
    df['block_id'] = pd.factorize(df['block_id'])[0] + 1

    return df

# get block pair coverage
def call_blocks(df_counts, df_hits):
    df_counts = df_counts[['#kmer','block_id']].copy()
    samples = [c for c in df_hits.columns if c != '#kmer']
    block_total = df_counts.groupby('block_id').size()

    # inner merge: kmers without a hits row are zeros, accounted for via block_total
    df = pd.merge(df_hits, df_counts[['#kmer', 'block_id']], on='#kmer', how='inner')

    block_hit = ((df[samples] > 0).groupby(df['block_id']).sum()
                 .reindex(block_total.index, fill_value=0))
    block_sum = (df[samples].groupby(df['block_id']).sum()
                 .reindex(block_total.index, fill_value=0))
    del df

    block_coverage = block_hit.div(block_total, axis=0).round(4)
    block_depth = block_sum.div(block_total, axis=0).round(2)

    df_gpd = pd.concat({'block_hit': block_hit,
                        'block_coverage': block_coverage,
                        'block_depth': block_depth}, axis=1)
    df_gpd = df_gpd.stack(level=-1, future_stack=True).rename_axis(['block_id', 'sample']).reset_index()
    df_gpd['sample'] = df_gpd['sample'].astype('category')
    df_gpd['block_total'] = df_gpd['block_id'].map(block_total)

    return df_gpd[['block_id',  'sample', 'block_hit', 'block_total', 'block_coverage', 'block_depth', ]]

def get_block_pair_coverage(df_called_blocks):
    # blocks x samples coverage matrix
    df_cov = (df_called_blocks
              .pivot(index='block_id', columns='sample', values='block_coverage')
              .fillna(0))
    df_cov.columns = df_cov.columns.astype(str)
    C = df_cov.to_numpy(dtype=np.float64)                   # B x S
    B = C.shape[0]
    P = B * (B - 1) // 2                                    # all block pairs

    # sum over all pairs a<b of sqrt(C[a]*C[b])
    tot = np.sqrt(C).sum(axis=0)
    pair_sum = (tot * tot - C.sum(axis=0)) / 2.0

    return pd.DataFrame({
        'sample':        df_cov.columns,
        'blocks_total':  B,
        'pairs_total':   P,
        'pair_cov_geo':  pair_sum / P if P else 0.0,
    })


def main():
    # to do, basename + strain name addition
    parser = argparse.ArgumentParser(description='Calcualting coverage and depth with kmer locations')
    parser.add_argument('--counts', help='a global_counts.tsv.zst file', required=True)
    parser.add_argument('--hits', help='a target_strain .kmer_hits.tsv.gz', required = True)
    parser.add_argument('--output_dir', help='directory where output is saved')
    parser.add_argument('--strain', help = 'add a name for the target strain you are tracking')
    parser.add_argument('--reference', '-r', help = 'Path to reference file of a strain. Can be used to map kmers to locations')
    parser.add_argument('--terminal_dist', type = int, default = 300, help = 'distance from end of contig. Kmers within that range will be dropped')
    parser.add_argument('--rare_min_count', type = int, default = 1000, help = 'Minimal rare kmers required, if less count 0 kmers are present after scrub, count 1 kmers are included')

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

    # get hit locations and strip terminal kmers
    if args.reference:
        df_locs, _ = build_mapped_kmers_ahocorasick(records= load_genome(args.reference),
                                                    terminal_dist= args.terminal_dist,
                                                    kmers = df_hits['#kmer'].to_list())
        print(f"Removing {len(df_locs.loc[df_locs['label']=='terminal'])} terminal kmers")
        df_internal = df_locs.loc[df_locs['label']=='internal']
        df_hits = df_hits.loc[df_hits['#kmer'].isin(df_internal['#kmer'])].copy()
        df_loc_counts = pd.merge(df_locs, df_counts, on = '#kmer', how = 'inner')
        # add block ids don't drop anything
        df_loc_counts = add_block_ids(df_loc_counts)

    # get sample coverage for all kemrs and unique/rare kmers
    df_cov_depth, df_cnt_uni = get_sample_coverage(df_hits=df_hits, df_counts= df_counts, rare_t = args.rare_min_count)
    
    if args.reference:
        df_loc_counts = pd.merge(df_loc_counts, df_cnt_uni, on = '#kmer', how = 'left')
        df_loc_counts.set_index('#kmer').to_csv(os.path.join(outpath,"kmer_locations.tsv"), sep = '\t')
        # create block pair coverage
        df_rare = df_loc_counts.dropna(subset = 'rare_count')
        df_blk_call = call_blocks(df_counts=df_rare, df_hits=df_hits)
        df_blk_call.set_index(['sample','block_id']).to_csv(os.path.join(outpath,"block_calls.tsv"), sep = '\t')
        df_blk_cov = get_block_pair_coverage(df_blk_call)

    df_cov_depth = df_cov_depth.reset_index(names = 'sample')
    df_cov_depth = pd.merge(df_cov_depth, df_blk_cov, on = 'sample', how = 'left')
    df_cov_depth = df_cov_depth.merge(total_evaluated, on='sample', how='left')
    # get relative abundance


    if args.strain:
        df_cov_depth['strain'] = args.strain 
    df_cov_depth.set_index(['strain', 'sample']).to_csv(os.path.join(outpath,"coverage_depth.tsv"), sep = '\t')
    print(f'Saved coverage dataframe to \n {outpath}')
    print(df_cov_depth.head(5))



if __name__ == '__main__':
    main()