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

# get kmer blocks: a new block starts at a new contig, a position gap, or a change in count
def add_block_ids(df, count_cols = ['pangenome_count', 'metagenome_count', 'drug_count']):
    df = df.sort_values(['contig_id', 'kmer_position']).copy()
    df['count'] = df[count_cols].sum(axis=1)
    new_contig = df['contig_id'].ne(df['contig_id'].shift())
    pos_break = df.groupby('contig_id')['kmer_position'].diff().ne(1)
    cnt_break = df['count'].ne(df['count'].shift())
    df['block_id'] = (new_contig | pos_break | cnt_break).cumsum()      # starts at 1
    df['block_size'] = df.groupby('block_id')['block_id'].transform('size')
    return df

#  get per-block coverage/depth; carries the (single) count value of each block
def call_blocks(df_blocks, df_hits):
    df_blocks = df_blocks[['#kmer', 'block_id', 'count']].copy()
    samples = [c for c in df_hits.columns if c != '#kmer']
    block_total = df_blocks.groupby('block_id').size()
 
    n_counts = df_blocks.groupby('block_id')['count'].nunique()
    assert (n_counts == 1).all(), 'blocks must have exactly one count value'
    block_count = df_blocks.groupby('block_id')['count'].first()
 
    # inner merge: kmers without a hits row are zeros, accounted for via block_total
    df = pd.merge(df_hits, df_blocks[['#kmer', 'block_id']], on='#kmer', how='inner')
 
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
    df_gpd['block_count'] = df_gpd['block_id'].map(block_count)
 
    return df_gpd[['block_id', 'sample', 'block_count', 'block_hit', 'block_total',
                   'block_coverage', 'block_depth']]

def get_block_pair_coverage(df_called_blocks):
    # blocks x samples coverage and depth matrices
    df_cov = (df_called_blocks
              .pivot(index='block_id', columns='sample', values='block_coverage')
              .fillna(0))
    df_cov.columns = df_cov.columns.astype(str)
    C = df_cov.to_numpy(dtype=np.float64)                   # B x S
    D = (df_called_blocks
         .pivot(index='block_id', columns='sample', values='block_depth')
         .fillna(0)
         .reindex(df_cov.index)
         .to_numpy(dtype=np.float64))                       # B x S, same order as C
    B = C.shape[0]
    P = B * (B - 1) // 2                                    # all block pairs

    # sum over all pairs a<b of sqrt(C[a]*C[b])
    tot = np.sqrt(C).sum(axis=0)
    pair_sum = (tot * tot - C.sum(axis=0)) / 2.0

    # sum over all pairs a<b of min(D[a], D[b]): k-th smallest (0-based) is the min of B-1-k pairs
    w = (B - 1 - np.arange(B))[:, None]
    depth_sum = (np.sort(D, axis=0) * w).sum(axis=0)

    return pd.DataFrame({
        'sample':         df_cov.columns,
        'blocks_total':   B,
        'pairs_total':    P,
        'pair_cov_geo':   np.round(pair_sum / P, 4) if P else 0.0,
        'pair_depth_min': np.round(depth_sum / P, 2) if P else 0.0,
    })

# create pairs and triplets with weights
def e2(v):                                                   # sum over pairs a<b of v[a]*v[b], per column
    return (v.sum(axis=0) ** 2 - (v * v).sum(axis=0)) / 2.0
 
 
def wmin_sum(V, W, k):
    """sum over all pairs (k=2) or triplets (k=3) of prod(weights) * min(values)"""
    idx = np.argsort(V, axis=0, kind='stable')
    Vs = np.take_along_axis(V, idx, 0)
    Ws = np.take_along_axis(W, idx, 0)
    S1 = Ws.sum(axis=0) - np.cumsum(Ws, axis=0)              # weight of the elements sorted after each one
    if k == 2:
        later = S1
    else:
        S2 = (Ws ** 2).sum(axis=0) - np.cumsum(Ws ** 2, axis=0)
        later = (S1 ** 2 - S2) / 2.0
    return (Vs * Ws * later).sum(axis=0)
 
 
def count_weight(count, alpha=1.0):
    """count = number of other genomes carrying the kmer; rarer -> heavier. count 0 -> weight 1"""
    return (1.0 + np.asarray(count, dtype=np.float64)) ** (-alpha)
 
 
def get_block_weighted_coverage(df_called_blocks, alpha=1.0, size_weight=False):
    """
    Weighted mean over all block pairs, weight of a pair = w_a * w_b,
    w = (1 + block_count)^-alpha, optionally multiplied by block_total (kmers per block).
    """
    def piv(col):
        return df_called_blocks.pivot(index='block_id', columns='sample', values=col).fillna(0)

    df_cov = piv('block_coverage')
    df_cov.columns = df_cov.columns.astype(str)
    C = df_cov.to_numpy(dtype=np.float64)                    # B x S coverage
    D = piv('block_depth').reindex(df_cov.index).to_numpy(dtype=np.float64)
    B, S = C.shape

    # one count / size per block, identical across samples
    blk = (df_called_blocks.groupby('block_id', observed=True)[['block_count', 'block_total']]
           .first().reindex(df_cov.index))
    w = count_weight(blk['block_count'].to_numpy(), alpha)
    if size_weight:
        w = w * blk['block_total'].to_numpy(dtype=np.float64)
    W = np.repeat(w[:, None], S, axis=1)                     # B x S (needed for take_along_axis)

    def ratio(num, den):
        num = np.asarray(num, dtype=np.float64)
        den = np.broadcast_to(np.asarray(den, dtype=np.float64), num.shape)
        return np.divide(num, den, out=np.zeros_like(num), where=den > 0)

    den2 = e2(W)                                             # sum over pairs of w_a*w_b

    # sum over pairs of w_a*w_b*(D_a + D_b) = sum_a w_a*D_a*(W_total - w_a)
    pair_sum = (W * D * (W.sum(axis=0) - W)).sum(axis=0)
    min_sum = wmin_sum(D, W, 2)
    max_sum = pair_sum - min_sum                             # min + max = D_a + D_b

    return pd.DataFrame({
        'sample':              df_cov.columns,
        'blocks_total_cw':     B,
        'pair_cov_geo_cw':     np.round(ratio(e2(W * np.sqrt(C)), den2), 4),
        'pair_depth_min_cw':   np.round(ratio(min_sum, den2), 2),
        'pair_depth_mean_cw':  np.round(ratio(pair_sum / 2.0, den2), 2),
        'pair_depth_max_cw':   np.round(ratio(max_sum, den2), 2),
        'count_weight_alpha':  alpha,
    })

def depth_weight(D, baseline, beta=1.0, mode='over', pseudo=1.0):
    """
    per block x sample weight from depth discrepancy to the rare-kmer baseline depth.
    r = (D + pseudo) / (baseline + pseudo)
    mode 'over': penalise only depth above baseline (shared regions inflated by relatives)
    mode 'both': penalise both directions (careful: zero-depth blocks then drop out of coverage)
    """
    r = (D + pseudo) / (np.asarray(baseline, dtype=np.float64)[None, :] + pseudo)
    dev = np.maximum(r, 1.0) if mode == 'over' else np.maximum(r, 1.0 / r)
    return dev ** (-beta)


def get_block_depth_weighted_coverage(df_called_blocks, baseline, beta=1.0, alpha=0.0,
                                      mode='over', pseudo=1.0, size_weight=False):
    """
    Weighted pair statistics over all blocks; block weight = count weight * depth weight.
    baseline: pd.Series sample -> baseline depth (e.g. depth_rare_kmers)
    alpha = 0 -> depth weighting only.
    """
    def piv(col):
        return df_called_blocks.pivot(index='block_id', columns='sample', values=col).fillna(0)

    df_cov = piv('block_coverage')
    df_cov.columns = df_cov.columns.astype(str)
    C = df_cov.to_numpy(dtype=np.float64)                    # B x S coverage
    D = piv('block_depth').reindex(df_cov.index).to_numpy(dtype=np.float64)
    base = baseline.reindex(df_cov.columns).fillna(0).to_numpy(dtype=np.float64)

    blk = (df_called_blocks.groupby('block_id', observed=True)[['block_count', 'block_total']]
           .first().reindex(df_cov.index))
    w = count_weight(blk['block_count'].to_numpy(), alpha)
    if size_weight:
        w = w * blk['block_total'].to_numpy(dtype=np.float64)
    W = w[:, None] * depth_weight(D, base, beta=beta, mode=mode, pseudo=pseudo)   # B x S, per sample

    def ratio(num, den):
        return np.divide(num, den, out=np.zeros_like(num, dtype=np.float64), where=den > 0)

    den2 = e2(W)
    pair_sum = (W * D * (W.sum(axis=0) - W)).sum(axis=0)
    min_sum = wmin_sum(D, W, 2)

    return pd.DataFrame({
        'sample':             df_cov.columns,
        'depth_baseline_dw':  base,
        'eff_blocks_dw':      np.round(ratio(W.sum(axis=0) ** 2, (W ** 2).sum(axis=0)), 1),
        'pair_cov_geo_dw':    np.round(ratio(e2(W * np.sqrt(C)), den2), 4),
        'pair_depth_min_dw':  np.round(ratio(min_sum, den2), 2),
        'pair_depth_mean_dw': np.round(ratio(pair_sum / 2.0, den2), 2),
        'pair_depth_max_dw':  np.round(ratio(pair_sum - min_sum, den2), 2),
        'depth_weight_beta':  beta,
        'count_weight_alpha': alpha,
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
    parser.add_argument('--weight_alpha', type=float, default=0 ,help='block weight = (1 + count)^-alpha; higher favours rare blocks more')
    parser.add_argument('--size_weight', action='store_true', help='also multiply block weight by kmers per block')
    parser.add_argument('--depth_beta', type=float, default= 0,
                                        help='depth weight = max(r,1)^-beta, r = (block_depth+pseudo)/(rare_depth+pseudo)')
    parser.add_argument('--depth_mode', choices=['over', 'both'], default='over',
                                        help="'over' penalises only depth above the rare baseline")
    parser.add_argument('--depth_pseudo', type=float, default=1.0, help='pseudocount for the depth ratio')
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
    total_evaluated['total_strain_kmers'] = len(df_counts)
    
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
        print(df_loc_counts)
        
    
    # get sample coverage for all kemrs and unique/rare kmers
    df_cov_depth, df_cnt_uni = get_sample_coverage(df_hits=df_hits, df_counts= df_counts, rare_t = args.rare_min_count)
    df_cov_depth = df_cov_depth.reset_index(names='sample')

    if args.reference:
        df_loc_counts = pd.merge(df_loc_counts, df_cnt_uni[['#kmer', 'rare_count']], on='#kmer', how='left')
        df_loc_counts.set_index('#kmer').to_csv(os.path.join(outpath, "kmer_locations.tsv"), sep='\t')

        ### block calls for all kmers (rare blocks are a subset: blocks never mix rare and non-rare)
        df_blk_call_all = call_blocks(df_blocks=df_loc_counts, df_hits=df_hits)
        rare_blocks = df_loc_counts.loc[df_loc_counts['rare_count'].notna(), 'block_id'].unique()
        df_blk_call_all['rare_block'] = df_blk_call_all['block_id'].isin(rare_blocks)
        df_blk_call_all.set_index(['sample', 'block_id']).to_csv(os.path.join(outpath, "block_calls.tsv"), sep='\t')

        df_blk_call_rare = df_blk_call_all.loc[df_blk_call_all['rare_block']]

        df_blk_cov_w_all = get_block_weighted_coverage(df_blk_call_all, alpha=args.weight_alpha,
                                                       size_weight=args.size_weight)
        df_blk_cov_w_rare = get_block_weighted_coverage(df_blk_call_rare, alpha=args.weight_alpha,
                                                        size_weight=args.size_weight)
        # depth based weights
        baseline = df_cov_depth.set_index('sample')['depth_rare_kmers']
        df_blk_cov_dw = get_block_depth_weighted_coverage(df_blk_call_all, baseline,
                                                          beta=args.depth_beta, alpha=args.weight_alpha,
                                                          mode=args.depth_mode, pseudo=args.depth_pseudo,
                                                          size_weight=args.size_weight)

        df_cov_depth = (df_cov_depth
                        .merge(df_blk_cov_w_all.set_index('sample').add_prefix('all_').reset_index(),
                               on='sample', how='left')
                        .merge(df_blk_cov_w_rare.set_index('sample').add_prefix('rare_').reset_index(),
                               on='sample', how='left')
                        .merge(df_blk_cov_dw, on='sample', how='left'))
        
    # get relative abundance
    df_cov_depth = df_cov_depth.merge(total_evaluated, on='sample', how='left')
    df_cov_depth['relative_abundance'] =round((df_cov_depth['depth_rare_kmers'] * df_cov_depth['total_strain_kmers'] / total_evaluated['total_evaluated']) * 100,4)

    if args.strain:
        df_cov_depth['strain'] = args.strain 
    df_cov_depth.set_index(['strain', 'sample']).to_csv(os.path.join(outpath,"coverage_depth.tsv"), sep = '\t')
    print(f'Saved coverage dataframe to \n {outpath}')
    print(df_cov_depth.head(5))



if __name__ == '__main__':
    main()