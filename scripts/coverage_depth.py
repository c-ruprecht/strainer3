import pandas as pd
import os
import argparse
import numpy as np
import ahocorasick
from Bio import SeqIO
from Bio.Seq import Seq
import gzip
#plotting
import plotly.express as px
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from pycirclize import Circos
from scipy.stats import poisson
import re
from matplotlib.colors import LogNorm, Normalize
from matplotlib.cm import ScalarMappable
from scipy.optimize import brentq



def robust_depth(x, f_bad=0.02, q=1 - 1e-4, n_iter=20):
    """Depth estimate robust to a high-depth k-mer subset and an assembly zero floor.
    Returns lam, coverage on kept k-mers, fraction of k-mers trimmed as high-depth, threshold."""
    x = np.asarray(x)
    n = len(x)
    n1, n2 = (x == 1).sum(), (x == 2).sum()
    lam = 2 * n2 / n1 if n1 > 0 else np.median(x)      # low-count start: ignores zeros + tail
    lam = max(lam, 1e-6)

    for _ in range(n_iter):
        thr = poisson.ppf(q, lam)
        keep = x <= thr
        p0 = (x[keep] == 0).mean() if keep.any() else 1.0
        p0_adj = np.clip((p0 - f_bad) / (1 - f_bad), 1 / n, 1)   # remove assembly floor
        lam_new = -np.log(p0_adj)
        if abs(lam_new - lam) < 1e-6 * max(lam, 1):
            lam = lam_new
            break
        lam = max(lam_new, 1e-6)

    thr = poisson.ppf(q, lam)
    keep = x <= thr
    return dict(depth_robust=lam,
                cov_kept=(x[keep] > 0).mean() if keep.any() else 0.0,
                frac_trimmed=1 - keep.mean(),
                thr=thr,
                lam_mean=x.mean())



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

    max_t = 1
    for t in range(0, max_t + 1):
        df_counts_uni = df_counts.loc[df_counts['count'] <= t].copy()
        if len(df_counts_uni) >= rare_t:
            break
    rare_count = '==0' if t == 0 else f'<={t}'
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

    return df_cov_depth, df_counts_uni, df_unique, df_nonunique

def ztp_lambda(m_pos):
    """MLE of λ for a zero-truncated Poisson with mean m_pos."""
    if not np.isfinite(m_pos) or m_pos <= 1:
        return 0.0            # all positives are 1s → λ → 0
    return brentq(lambda l: l / (1 - np.exp(-l)) - m_pos, 1e-9, m_pos)

def zero_inflate_poisson(df_hits):
    df = df_hits.set_index('#kmer')
    hits = (df > 0).sum(axis=0)
    total = len(df)
    
    out = pd.DataFrame({'obs_0_cnts': total - hits ,
                        'coverage': (df > 0).sum(axis=0)/len((df > 0)),
                        'depth_pos': df[df > 0].mean(axis=0),
                        })
    
    out['lambda_ztp'] = out['depth_pos'].apply(ztp_lambda)
    p0 = np.exp(-out['lambda_ztp'])

    # Zeros expected among the k-mers that are actually present
    exp0 = (hits * p0 / (1 - p0)).replace([np.inf, -np.inf], np.nan)
    out['exp_0_cnts'] = exp0.round().astype('Int64')    # Coverage expected if every k-mer were present
    out['exp_coverage'] = 1 - p0

    out['excess_0'] = out['obs_0_cnts'] - out['exp_0_cnts']     # structurally absent k-mers
    out['frac_absent'] = (out['excess_0'] / total).clip(lower=0) 
    out['zi_index'] = out['excess_0'] / out['obs_0_cnts']              # share of zeros that are structural
    out['zero_ratio'] = out['obs_0_cnts'] / out['exp_0_cnts']          # obs/exp zeros; 1 = pure Poisson
    out = out.reset_index(names=['sample'])
    return out

def zero_inflate_stats(obs_0, total, lam, exp0):
    """shared output columns for the zero-inflation estimators (same definitions as zero_inflate_poisson)"""
    excess_0 = obs_0 - exp0
    return {'lambda': lam,
            'exp_0_cnts': exp0,
            'exp_coverage': 1 - np.exp(-lam) if np.isfinite(lam) else np.nan,
            'excess_0': excess_0,
            'frac_absent': np.clip(excess_0 / total, 0, 1) if np.isfinite(exp0) else np.nan,   # structurally absent kmers
            'zi_index': excess_0 / obs_0 if obs_0 > 0 else 0.0,                               # share of zeros that are structural
            'zero_ratio': obs_0 / exp0 if exp0 > 0 else np.inf}                               # obs/exp zeros; 1 = pure Poisson

def lowcount_0_inflation(x, min_f1=5):
    """
    Low-count estimator: depth only from kmers seen once (f1) and twice (f2), ignores zeros and the high-depth tail.
    Poisson: f2/f1 = lam/2 -> expected zeros f0 = f1/lam = f1^2/(2 f2)  (Chao1, bias corrected: f1(f1-1)/(2(f2+1)))
    """
    f = np.bincount(np.asarray(x, dtype=np.int64), minlength=3)
    total, obs_0, f1, f2 = f.sum(), f[0], f[1], f[2]
    if f1 < min_f1:
        # few singletons: deep coverage (no zeros expected) or no signal at all
        deep = (total - obs_0) >= min_f1
        out = zero_inflate_stats(obs_0, total, np.nan, 0.0 if deep else np.nan)
        return {**out, 'f1': f1, 'f2': f2, 'reliable': deep}
    exp0 = f1 * (f1 - 1) / (2.0 * (f2 + 1))
    lam = 2.0 * (f2 + 1) / f1
    return {**zero_inflate_stats(obs_0, total, lam, exp0), 'f1': f1, 'f2': f2, 'reliable': True}

def trunc_poisson_lambda(f, thr):
    """MLE of lam for a Poisson truncated to 1..thr, from the count histogram f"""
    k = np.arange(len(f))
    keep = (k >= 1) & (k <= thr)
    n_kept = f[keep].sum()
    if n_kept == 0:
        return np.nan, 0
    m = (k[keep] * f[keep]).sum() / n_kept
    if m <= 1:                                     # all positives are 1s -> lam -> 0
        return 0.0, n_kept
    # E[X | 1 <= X <= thr] = lam * P(X <= thr-1) / (P(X <= thr) - P(0))
    g = lambda l: l * poisson.cdf(thr - 1, l) / (poisson.cdf(thr, l) - np.exp(-l)) - m
    hi = max(2 * m, 10.0)
    while g(hi) < 0 and hi < 1e7:
        hi *= 2
    return (brentq(g, 1e-9, hi) if g(hi) >= 0 else np.nan), n_kept

def robust_0_inflation(x, q=0.99, min_pos=10, n_iter=30):
    """
    Robust ZTP: zero- and right-truncated Poisson fitted on the Poisson body of the positives.
    start lam = 2 f2/f1, thr = Poisson q-quantile, kmers > thr are high-depth (present, but not used for lam),
    refit lam on 1..thr, iterate until thr is stable. Zeros are never used to estimate lam.
    """
    f = np.bincount(np.asarray(x, dtype=np.int64), minlength=3)
    total, obs_0 = f.sum(), f[0]
    pos = total - obs_0
    if pos < min_pos:
        return {**zero_inflate_stats(obs_0, total, np.nan, np.nan), 'frac_trimmed': np.nan, 'thr': np.nan, 'reliable': False}

    k = np.arange(len(f))
    lam = 2.0 * f[2] / f[1] if f[1] > 0 and f[2] > 0 else np.median(np.repeat(k[1:], f[1:]))
    lam = max(lam, 1e-3)
    thr_prev = None
    for _ in range(n_iter):
        thr = int(max(poisson.ppf(q, lam), 2))
        if thr == thr_prev:
            break
        lam_new, _ = trunc_poisson_lambda(f, thr)
        if not np.isfinite(lam_new):
            break
        lam, thr_prev = lam_new, thr

    thr = int(max(poisson.ppf(q, lam), 2)) if lam > 0 else 2
    lam_fit, n_kept = trunc_poisson_lambda(f, thr)
    lam = lam_fit if np.isfinite(lam_fit) else lam
    frac_trimmed = f[thr + 1:].sum() / pos
    if lam <= 0:
        return {**zero_inflate_stats(obs_0, total, 0.0, np.nan), 'frac_trimmed': frac_trimmed, 'thr': thr, 'reliable': False}

    # present kmers in the body = kept / P(1 <= X <= thr); zeros expected among them
    exp0 = n_kept / (poisson.cdf(thr, lam) - np.exp(-lam)) * np.exp(-lam)
    return {**zero_inflate_stats(obs_0, total, lam, exp0), 'frac_trimmed': frac_trimmed, 'thr': thr, 'reliable': True}

def zero_inflate_lowcount(df_hits, min_f1=5):
    df = df_hits.set_index('#kmer')
    out = pd.DataFrame.from_dict({s: lowcount_0_inflation(df[s].to_numpy(), min_f1=min_f1) for s in df.columns}, orient='index')
    out = out.reset_index(names=['sample'])
    return out

def zero_inflate_robust(df_hits, q=0.99, min_pos=10):
    df = df_hits.set_index('#kmer')
    out = pd.DataFrame.from_dict({s: robust_0_inflation(df[s].to_numpy(), q=q, min_pos=min_pos) for s in df.columns}, orient='index')
    out = out.reset_index(names=['sample'])
    return out



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
# get kmer blocks: a new block starts at a new contig, a position gap, or (if count_cols is given) a change in count
def add_block_ids(df, count_cols = ['pangenome_count', 'metagenome_count', 'drug_count']):
    df = df.sort_values(['contig_id', 'kmer_position']).copy()

    new_contig = df['contig_id'].ne(df['contig_id'].shift())
    pos_break  = df.groupby('contig_id')['kmer_position'].diff().ne(1)
    brk = new_contig | pos_break

    if count_cols:
        df['count'] = df[count_cols].sum(axis=1)
        brk = brk | df['count'].ne(df['count'].shift())

    df['block_id'] = brk.cumsum()                                         # starts at 1
    df['block_size'] = df.groupby('block_id')['block_id'].transform('size')
    return df

#  get per-block coverage/depth; carries the (single) count value of each block
def call_blocks(df_blocks, df_hits, counts='count'):
    samples = [c for c in df_hits.columns if c != '#kmer']
    cols = ['#kmer', 'block_id'] + ([counts] if counts else [])
    df_blocks = df_blocks[cols].copy()
    block_total = df_blocks.groupby('block_id').size()

    block_count = None
    if counts:
        n_counts = df_blocks.groupby('block_id')[counts].nunique()
        assert (n_counts == 1).all(), 'blocks must have exactly one count value'
        block_count = df_blocks.groupby('block_id')[counts].first()

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

    out_cols = ['block_id', 'sample', 'block_hit', 'block_total', 'block_coverage', 'block_depth']
    if block_count is not None:
        df_gpd['block_count'] = df_gpd['block_id'].map(block_count)
        out_cols.insert(2, 'block_count')
    return df_gpd[out_cols]

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
 
def get_weighted_coverage(df_hits, df_counts, alpha = 1):

    df_counts['count_weight'] = df_counts['count']
    return df

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

### plots
def plot_map(df_loc, outdir, prefix, df_hits=None, samples=None, ymax_min=20, s=3,
             quality='both', hit_cmap='viridis', hit_log=True, hit_q=0.995):
    """
    df_loc : kmer locations with counts (output of add_block_ids + rare_count merge)
    df_hits: optional '#kmer' x samples hit table -> one extra map per sample with a hits track
    samples: optional subset of sample columns to plot (default: all)
    """
    df_loc = df_loc.copy()
    df_loc['count'] = (df_loc[['pangenome_count', 'metagenome_count', 'drug_count']]
                       .fillna(0).sum(axis=1).astype(int))
    df_loc['rare_count'] = df_loc['rare_count'].fillna('not_rare')

    # ---------- shared layout: one sector, contigs concatenated ----------
    contigs = (df_loc.drop_duplicates('contig_id')[['contig_id', 'contig_length']]
                     .sort_values('contig_length', ascending=False))
    contigs['offset'] = contigs['contig_length'].cumsum() - contigs['contig_length']
    genome_len = int(contigs['contig_length'].sum())

    df = df_loc.merge(contigs[['contig_id', 'offset']], on='contig_id')
    df['gpos'] = df['kmer_position'] + df['offset']
    df = df.sort_values('gpos').reset_index(drop=True)

    # positional blocks (no count break), as extents on the concatenated genome axis
    blk = add_block_ids(df[['contig_id', 'kmer_position', 'gpos']], count_cols=None)
    blocks = blk.groupby('block_id')['gpos'].agg(start='min', end='max').sort_index()
    blocks['end'] = blocks['end'] + 1                         # a block of n kmers spans n positions

    ymax = max(ymax_min, int(df['count'].max()))

    def base_circos(with_hits=False):
        r_ctg, r_cnt, r_blk, r_hit = (((90, 94), (68, 86), (62, 66), (36, 58)) if with_hits
                                      else ((78, 82), (55, 75), None, None))
        circos = Circos({'genome': genome_len}, start=0, end=360)
        sector = circos.sectors[0]
        track_ctg = sector.add_track(r_ctg)
        for i, (st, ln) in enumerate(zip(contigs['offset'], contigs['contig_length'])):
            track_ctg.rect(st, st + ln, fc='grey' if i % 2 else 'lightgrey', ec='none')
        track_cnt = sector.add_track(r_cnt, r_pad_ratio=0.1)
        track_cnt.axis()
        track_hit = None
        if with_hits:
            track_blk = sector.add_track(r_blk)
            for i, (st, en) in enumerate(zip(blocks['start'], blocks['end'])):
                track_blk.rect(st, en, fc='grey' if i % 2 else 'lightgrey', ec='none')
            track_hit = sector.add_track(r_hit, r_pad_ratio=0.1)
            track_hit.axis()
        return circos, track_cnt, track_hit

    def draw_groups(track, groups, colors):
        for key in groups:                       # later groups are drawn on top
            d = df[df['_group'] == key]
            if len(d):
                track.scatter(d['gpos'].to_numpy(), d['count'].to_numpy(),
                              vmin=0, vmax=ymax, s=s, color=colors[key])

    def save(circos, handles, title, name, cbar=None):
        fig = circos.plotfig()
        circos.ax.set_title(prefix if cbar is None else f'{prefix} | {name.split(".", 1)[-1]}',
                            fontsize=11, pad=12)
        circos.ax.legend(handles=handles, title=title,
                         loc='upper left' if cbar else 'center left',
                         bbox_to_anchor=(1.02, 0.95) if cbar else (1.02, 0.5),
                         frameon=False, fontsize=8)
        if cbar is not None:
            norm, cmap, label = cbar
            cax = circos.ax.inset_axes([1.05, 0.05, 0.025, 0.35])
            cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax)
            cb.set_label(label, fontsize=8)
            cb.ax.tick_params(labelsize=7)
        base = os.path.join(outdir, f"{prefix}.{name}")
        fmts = ['svg', 'png'] if quality == 'both' else [quality]
        for fmt in fmts:
            kw = dict(format='svg') if fmt == 'svg' else dict(dpi=300)
            fig.savefig(f'{base}.{fmt}', bbox_inches='tight', **kw)
            print(f"saved {base}.{fmt}")
        plt.close(fig)

    # ---------- 1) one distinct colour per count ----------
    palette = px.colors.qualitative.Dark24
    count_colors = {c: palette[c % len(palette)] for c in range(ymax + 1)}
    present = sorted(int(c) for c in df['count'].unique())

    df['_group'] = df['count'].astype(int)
    circos, track, _ = base_circos()
    draw_groups(track, present, count_colors)
    save(circos, [Patch(facecolor=count_colors[c], label=str(c)) for c in present],
         'k-mer count', 'by_count')

    # ---------- 2) rare vs not_rare ----------
    colors2 = {'rare': 'blue', 'not_rare': 'red'}
    df['_group'] = np.where(df['rare_count'] == 'not_rare', 'not_rare', 'rare')
    circos, track, _ = base_circos()
    draw_groups(track, ['rare', 'not_rare'], colors2)
    save(circos, [Patch(facecolor=colors2[k], label=k.replace('_', ' '))
                  for k in ['not_rare', 'rare']],
         'k-mer class', 'rare')

    # ---------- 3) rare vs not_rare vs not_rare terminal ----------
    colors3 = {'rare': 'blue', 'not_rare': 'red', 'not_rare_terminal': 'orange'}
    groups3 = ['rare', 'not_rare', 'not_rare_terminal']
    df['_group'] = np.select(
        [(df['rare_count'] == 'not_rare') & (df['terminal_kmer'] == True),
         df['rare_count'] == 'not_rare'],
        ['not_rare_terminal', 'not_rare'],
        default='rare',
    )
    handles3 = [Patch(facecolor=colors3[k], label=k.replace('_', ' ')) for k in groups3[::-1]]
    circos, track, _ = base_circos()
    draw_groups(track, groups3, colors3)
    save(circos, handles3, 'k-mer class', 'rare_terminal')

    # ---------- 4) per sample: count track (by_count) + block track + hits track (colour = hit count) ----------
    if df_hits is not None:
        hit_cols = [c for c in df_hits.columns if c != '#kmer']
        samples = hit_cols if samples is None else [c for c in samples if c in hit_cols]
        # kmers without a hits row are zero hits
        dh = (df[['#kmer', 'gpos']]
              .merge(df_hits[['#kmer'] + samples], on='#kmer', how='left')
              .fillna({c: 0 for c in samples}))
        pos_all = dh['gpos'].to_numpy()
        cmap = plt.get_cmap(hit_cmap)

        df['_group'] = df['count'].astype(int)               # by_count colouring for the count track
        handles_cnt = [Patch(facecolor=count_colors[c], label=str(c)) for c in present]

        for smp in samples:
            h = dh[smp].to_numpy(dtype=np.float64)
            nz = h > 0
            hmax = max(2.0, float(np.quantile(h[nz], hit_q))) if nz.any() else 2.0
            if hit_log:
                norm = LogNorm(vmin=1, vmax=hmax)
                yfun = lambda v: np.log10(1 + np.minimum(v, hmax))
            else:
                norm = Normalize(vmin=0, vmax=hmax)
                yfun = lambda v: np.minimum(v, hmax)
            y_max = float(yfun(np.array(hmax)))

            circos, track_cnt, track_hit = base_circos(with_hits=True)
            draw_groups(track_cnt, present, count_colors)

            if (~nz).any():                                  # zero-hit kmers: grey baseline
                track_hit.scatter(pos_all[~nz], np.zeros((~nz).sum()),
                                  vmin=0, vmax=y_max, s=s, color='lightgrey')
            if nz.any():
                order = np.argsort(h[nz], kind='stable')     # high hits drawn on top
                hv = h[nz][order]
                track_hit.scatter(pos_all[nz][order], yfun(hv),
                                  vmin=0, vmax=y_max, s=s,
                                  c=np.minimum(hv, hmax), cmap=cmap, norm=norm)

            safe = re.sub(r'[^\w.-]', '_', str(smp))
            save(circos, handles_cnt + [Patch(facecolor='lightgrey', label='0 hits')],
                 'k-mer count', f'hits.{safe}',
                 cbar=(norm, cmap, 'hit count' + (' (log)' if hit_log else '')))

    return df.drop(columns='_group')

def plot_hit_distribution(df_hits, outdir, prefix, samples=None, quality='both'):
    """One histogram of k-mer hit counts per sample."""
    df = df_hits.set_index('#kmer') if '#kmer' in df_hits.columns else df_hits
    samples = list(df.columns) if samples is None else [s for s in samples if s in df.columns]

    paths = []
    for smp in samples:
        fig = px.histogram(df[smp], width=800, template='simple_white',
                           title=f'{prefix} | {smp}')
        safe = re.sub(r'[^\w.-]', '_', str(smp))
        base = os.path.join(outdir, f'{prefix}.hit_distribution.{safe}')
        for fmt in (['svg', 'png'] if quality == 'both' else [quality]):
            path = f'{base}.{fmt}'
            fig.write_image(path, scale=2 if fmt == 'png' else None)
            print(f'saved {path}')
            paths.append(path)
    return paths





def main():
    # to do, basename + strain name addition
    parser = argparse.ArgumentParser(description='Calcualting coverage and depth with kmer locations')
    parser.add_argument('--counts', help='a global_counts.tsv.zst file', required=True)
    parser.add_argument('--hits', help='a target_strain .kmer_hits.tsv.gz', required = True)
    parser.add_argument('--output_dir', help='directory where output is saved')
    parser.add_argument('--strain', help = 'add a name for the target strain you are tracking')
    parser.add_argument('--reference', '-r', help = 'Path to reference file of a strain. Can be used to map kmers to locations')
    parser.add_argument('--terminal_dist', type = int, default = 0, help = 'distance from end of contig. Kmers within that range will be dropped')
    parser.add_argument('--rare_min_count', type = int, default = 1000, help = 'Minimal rare kmers required, if less count 0 kmers are present after scrub, count 1 kmers are included')
    parser.add_argument('--figures', choices=['png','svg', 'both','none'], default='none', help = 'Change to svg, png or both if you want figures to be saved.')
    
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    hits_path = args.hits
    counts_path = args.counts
    outpath = args.output_dir


    df_hits = pd.read_csv(hits_path, sep = '\t')
    df_counts = pd.read_csv(counts_path, sep = '\t')
    
    # get total evaluated kmers for relative abundance evaluation 
    total_evaluated = df_hits.loc[df_hits['#kmer']=='total_evaluated'].drop(columns='#kmer').T.reset_index()#.to_dict()
    total_evaluated.columns = ['sample', 'total_evaluated']
    total_evaluated['total_strain_kmers'] = len(df_counts)
    print(total_evaluated)
    df_hits = df_hits.loc[df_hits['#kmer'] !='total_evaluated'].copy()

    
    
    # get hit locations and strip terminal kmers
    if args.reference:
        df_locs, _ = build_mapped_kmers_ahocorasick(records= load_genome(args.reference),
                                                    terminal_dist= args.terminal_dist,
                                                    kmers = df_hits['#kmer'].to_list())
        
        print(f"Removing {len(df_locs.loc[df_locs['label']=='terminal'])} terminal kmers")
        
        df_internal = df_locs.loc[df_locs['label']=='internal']
        df_hits = df_hits.loc[df_hits['#kmer'].isin(df_internal['#kmer'])].copy()
        df_loc_counts = pd.merge(df_internal, df_counts, on='#kmer', how='inner')
        
        # add block ids don't drop anything, blocks are only split by position
        df_loc_counts = add_block_ids(df_loc_counts, count_cols= None)
        
    
    # get sample coverage for all kemrs and unique/rare kmers
    df_cov_depth, df_cnt_uni, df_hits_unique, df_hits_shared = get_sample_coverage(df_hits=df_hits, df_counts= df_counts, rare_t = args.rare_min_count)
    df_cov_depth = df_cov_depth.reset_index(names='sample')

    if args.reference:
        df_loc_counts = pd.merge(df_loc_counts, df_cnt_uni[['#kmer', 'rare_count']], on='#kmer', how='left')
        df_loc_counts.set_index('#kmer').to_csv(os.path.join(outpath, "kmer_locations.tsv"), sep='\t')
        df_blocks = call_blocks(df_hits= df_hits, df_blocks=df_loc_counts, counts = None)
        df_blocks.to_csv(os.path.join(outpath, "blocks.tsv"), sep='\t')

    # get theoretical 0 counts using poisson of non 0s
    df_zpt_all = zero_inflate_poisson(df_hits).set_index('sample').add_prefix('zptall_').reset_index()
    df_zpt_rare = zero_inflate_poisson(df_hits_unique.reset_index()).set_index('sample').add_prefix('zptrare_').reset_index()
    df_zpt_share = zero_inflate_poisson(df_hits_shared.reset_index()).set_index('sample').add_prefix('zptshare_').reset_index()

    # low-count estimator (f1, f2) -> robust to high-depth tail
    df_zlc_all = zero_inflate_lowcount(df_hits).set_index('sample').add_prefix('zlcall_').reset_index()
    df_zlc_rare = zero_inflate_lowcount(df_hits_unique.reset_index()).set_index('sample').add_prefix('zlcrare_').reset_index()
    df_zlc_share = zero_inflate_lowcount(df_hits_shared.reset_index()).set_index('sample').add_prefix('zlcshare_').reset_index()

    # robust zero truncated poisson, high depth kmers trimmed
    df_zrob_all = zero_inflate_robust(df_hits).set_index('sample').add_prefix('zroball_').reset_index()
    df_zrob_rare = zero_inflate_robust(df_hits_unique.reset_index()).set_index('sample').add_prefix('zrobrare_').reset_index()
    df_zrob_share = zero_inflate_robust(df_hits_shared.reset_index()).set_index('sample').add_prefix('zrobshare_').reset_index()

    df_blocks_gpd = df_blocks.groupby('sample').agg(**{'block_count': ('block_id', 'count'),
                                                       '0_cov_blocks': ('block_coverage', lambda x: sum(x==0))})
 

    df_cov_depth = (df_cov_depth.merge(df_zpt_all, on='sample', how='left')
                    .merge(df_zpt_rare, on='sample', how='left')
                    .merge(df_zpt_share, on='sample', how='left')
                    .merge(df_zlc_all, on='sample', how='left')
                    .merge(df_zlc_rare, on='sample', how='left')
                    .merge(df_zlc_share, on='sample', how='left')
                    .merge(df_zrob_all, on='sample', how='left')
                    .merge(df_zrob_rare, on='sample', how='left')
                    .merge(df_zrob_share, on='sample', how='left')
                    .merge(df_blocks_gpd, on='sample', how='left'))

    # get relative abundance
    df_cov_depth = df_cov_depth.merge(total_evaluated, on='sample', how='left')
    df_cov_depth['relative_abundance'] =round((df_cov_depth['depth_rare_kmers'] * df_cov_depth['total_strain_kmers'] / df_cov_depth['total_evaluated']) * 100,4)
    df_cov_depth['relative_abundance_zroball'] = round((df_cov_depth['zroball_lambda'] * df_cov_depth['total_strain_kmers'] / df_cov_depth['total_evaluated']) * 100, 4)
    df_cov_depth['relative_abundance_zrobrare'] = round((df_cov_depth['zrobrare_lambda'] * df_cov_depth['total_strain_kmers'] / df_cov_depth['total_evaluated']) * 100, 4)
    
    # Exports and figures
    if args.strain:
        df_cov_depth['strain'] = args.strain 
    df_cov_depth.set_index(['strain', 'sample']).round(4).to_csv(os.path.join(outpath,"coverage_depth.tsv"), sep = '\t')
    print(f'Saved coverage dataframe to \n {outpath}')
    print(df_cov_depth.head(5))

    if args.figures != 'none':
        if not args.reference:
            plot_hit_distribution(df_hits, args.output_dir, basename, quality=args.figures)
        basename = args.strain if args.strain else 'contig_map'
        plot_map(df_loc_counts, args.output_dir, basename,
                 df_hits=df_hits, quality=args.figures)
        plot_hit_distribution(df_hits, args.output_dir, basename, quality=args.figures)


        


if __name__ == '__main__':
    main()