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
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from pycirclize import Circos
from scipy.stats import poisson
from scipy.special import gammaln
import re
import json
from matplotlib.colors import LogNorm, Normalize
from matplotlib.cm import ScalarMappable
from scipy.optimize import brentq



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


### zero inflation: mixture of zero-truncated Poissons, target depth = lowest real peak
# fitting settings; all exposed on the command line (--zi_*) and written to zi_params.json
ZI_PARAMS = dict(k_max=3,          # max number of Poisson components (chosen by BIC)
                 min_weight=0.05,  # min share of the positives for a real peak
                 noise_lam=1.0,    # peaks below this λ are low-count noise if deeper peaks hold the majority
                 dip=0.25,         # min density dip between two maxima to call them separate genomes
                 min_pos=10)       # min positive kmers to fit at all

def ztp_lambda(m_pos):
    """MLE of λ for a zero-truncated Poisson with mean m_pos (mean of the positive counts)"""
    if not np.isfinite(m_pos) or m_pos <= 1 + 1e-9:
        return 1e-6                                            # all positives are 1s -> λ -> 0
    return brentq(lambda l: l / (1 - np.exp(-l)) - m_pos, 1e-9, m_pos)

def ztp_mixture(x, K=2, n_iter=300, tol=1e-8):
    """EM for a mixture of K zero-truncated Poissons on the positive counts (histogram based). Zeros are never used."""
    f = np.bincount(np.asarray(x, dtype=np.int64)); k = np.arange(len(f))
    kk, w = k[1:][f[1:] > 0], f[1:][f[1:] > 0].astype(float)
    n = w.sum()
    lam = np.quantile(np.repeat(kk, w.astype(int)), (np.arange(K) + 0.5) / K) + 0.1   # starts spread over quantiles
    pi = np.full(K, 1 / K); ll_old = -np.inf
    with np.errstate(divide='ignore', invalid='ignore'):
        for _ in range(n_iter):
            logp = (kk[None, :] * np.log(lam[:, None]) - lam[:, None] - gammaln(kk[None, :] + 1)
                    - np.log1p(-np.exp(-lam[:, None])) + np.log(pi[:, None]))
            mx = logp.max(0); lse = mx + np.log(np.exp(logp - mx).sum(0))
            ll = (w * lse).sum()
            resp = np.exp(logp - lse)
            nk = (resp * w).sum(1)
            pi = nk / n
            lam = np.array([ztp_lambda((resp[j] * w * kk).sum() / nk[j]) if nk[j] > 0 else 1e-6 for j in range(K)])
            if ll - ll_old < tol * abs(ll):
                break
            ll_old = ll
    order = np.argsort(lam)
    return dict(lam=lam[order], pi=pi[order], n_pos=nk[order], bic=-2 * ll + (2 * K - 1) * np.log(n))

def mode_groups(lam, npos, real, dip=0.25):
    """
    Group the real mixture components into modes of the fitted density (k >= 1).
    Two neighbouring maxima are separate modes only if the density dips between them by >= dip of the lower maximum;
    otherwise they are one hump: overdispersed coverage of one genome, not two genomes.
    Returns the mode index of each real component (0 = lowest mode).
    """
    lr, nr = lam[real], npos[real]
    kmax = int(lr.max() + 5 * np.sqrt(max(lr.max(), 1)) + 5)
    k = np.arange(1, kmax + 1)
    dens = (nr[:, None] * poisson.pmf(k[None, :], lr[:, None]) / (1 - np.exp(-lr[:, None]))).sum(0)
    peaks = [i for i in range(len(k)) if dens[i] >= dens[max(i - 1, 0)] and dens[i] >= dens[min(i + 1, len(k) - 1)]
             and (i == 0 or dens[i] > dens[i - 1])]
    valleys = []                                             # count positions separating real modes
    cur = peaks[0] if peaks else 0
    for nxt in peaks[1:]:
        v = cur + int(np.argmin(dens[cur:nxt + 1]))
        if dens[v] <= (1 - dip) * min(dens[cur], dens[nxt]):
            valleys.append(k[v]); cur = nxt
        elif dens[nxt] > dens[cur]:
            cur = nxt
    return np.searchsorted(np.array(valleys, float), lr)

def robust_0_inflation(x, k_max=3, min_weight=0.05, noise_lam=1.0, dip=0.25, min_pos=10):
    """
    Robust zero inflation from a mixture of zero-truncated Poissons fitted on the positive counts (K = 1..k_max by BIC).
    Target depth = lowest mode: every target kmer carries at least the target depth, kmers shared with a relative
    sit higher, multi-copy kmers higher still.
      real peak  : weight >= min_weight; a peak with λ < noise_lam is low-count noise only if peaks >= 5 x noise_lam
                   hold the majority of the positives
      target     : real peaks forming the lowest mode of the fitted density. Peaks without a dip between them are one
                   hump (overdispersed coverage of one genome) and are merged.
    zi_frac  : (total - present kmers) / total         -> zero inflation fraction, share of the strain not in the sample
    zi_index : (obs_0 - exp_0) / obs_0                  -> zero inflation index, share of the zeros that are structural
    """
    x = np.asarray(x, dtype=np.int64)
    total, obs_0 = len(x), int((x == 0).sum())
    pos = total - obs_0
    nan = {'lambda': np.nan, 'lambda_present': np.nan, 'n_modes': np.nan, 'peaks': np.nan, 'target_weight': np.nan,
           'mode_ratio': np.nan,
           'zi_frac': np.nan, 'zi_index': np.nan, 'exp_0_cnts': np.nan}
    if pos < min_pos:
        return nan

    best = min((ztp_mixture(x, K) for K in range(1, k_max + 1)), key=lambda r: r['bic'])
    lam, pi, npos = best['lam'], best['pi'], best['n_pos']
    deep = pi[lam >= 5 * noise_lam].sum()                                     # positives held by clearly deeper peaks
    real = (pi >= min_weight) & ~((lam < noise_lam) & (deep >= 0.5))
    if not real.any():
        real = pi == pi.max()
    mode = np.full(len(lam), -1)
    mode[real] = mode_groups(lam, npos, real, dip=dip)
    target = mode == 0
    lam_t = np.average(lam[target], weights=npos[target])

    p0 = np.exp(-lam[real])
    present = (npos[real] / (1 - p0)).sum()                                   # real peaks incl. their zeros
    exp_0 = (npos[real] * p0 / (1 - p0)).sum()                                # zeros expected among present kmers
    excess_0 = obs_0 - exp_0
    return {'lambda': lam_t,                                                            # depth of the target (lowest) peak
            'lambda_present': np.average(lam[real], weights=npos[real]),                # depth over all present kmers
            'n_modes': int(mode.max() + 1),                                             # separate depth modes (genomes)
            'peaks': ';'.join(f'{l:.2f}:{w:.3f}:{m}' for l, w, m in zip(lam[real], pi[real], mode[real])),  # λ:weight:mode
            'target_weight': pi[target].sum(),                                          # share of positives in the target mode
            'mode_ratio': (np.average(lam[mode == 1], weights=npos[mode == 1]) / lam_t) if (mode == 1).any() else 1.0,
            'zi_frac': float(np.clip((total - present) / total, 0, 1)),
            'zi_index': excess_0 / obs_0 if obs_0 > 0 else 0.0,
            'exp_0_cnts': exp_0,
            # kept for plotting
            '_lams': lam, '_pis': pi, '_npos': npos, '_real': real, '_target': target}

def zero_inflate_robust(df_hits, **kw):
    df = df_hits.set_index('#kmer')
    out = pd.DataFrame.from_dict({s: robust_0_inflation(df[s].to_numpy(), **kw) for s in df.columns}, orient='index')
    out = out.drop(columns=[c for c in out.columns if c.startswith('_') or c == 'exp_0_cnts'])
    out = out.reset_index(names=['sample'])
    return out


def hit_histograms(subsets, strain):
    """
    kmer hit count histogram per sample and kmer subset: everything the zero-inflation fit uses.
    long format: strain, sample, subset, subset_kmers, count, n_kmers (only counts that occur)
    """
    rows = []
    for sub, df in subsets.items():
        df = df.set_index('#kmer') if '#kmer' in df.columns else df
        for smp in df.columns:
            f = np.bincount(df[smp].to_numpy(dtype=np.int64))
            k = np.nonzero(f)[0]
            rows.append(pd.DataFrame({'strain': strain, 'sample': smp, 'subset': sub, 'subset_kmers': len(df),
                                      'count': k, 'n_kmers': f[k]}))
    return pd.concat(rows, ignore_index=True)

def refit_from_histograms(df_hist, **params):
    """
    rerun the zero-inflation fit on exported histograms (e.g. to optimise ZI_PARAMS against a truth table)
    returns one row per strain x sample x subset
    """
    p = {**ZI_PARAMS, **params}
    rows = []
    for (strain, smp, sub), d in df_hist.groupby(['strain', 'sample', 'subset'], sort=False):
        x = np.repeat(d['count'].to_numpy(), d['n_kmers'].to_numpy())
        r = robust_0_inflation(x, **p)
        rows.append({'strain': strain, 'sample': smp, 'subset': sub,
                     **{k: v for k, v in r.items() if not k.startswith('_')}})
    return pd.DataFrame(rows)


def window_evenness(df_loc, df_hits, window=200):
    """
    Spatial patchiness along the genome. Mapped kmers are cut into windows of `window` consecutive kmers per contig.
    A present strain is covered evenly, a relative leaves whole regions uncovered (zeros cluster in windows).
    zero_window_frac     : fraction of windows with zero coverage
    window_cov_dispersion: variance of window coverage / variance expected if kmers were hit independently (1 = even)
    window_depth_cv      : coefficient of variation of window depth
    """
    samples = [c for c in df_hits.columns if c != '#kmer']
    loc = df_loc[['#kmer', 'contig_id', 'kmer_position']].sort_values(['contig_id', 'kmer_position']).copy()
    loc['window'] = loc.groupby('contig_id').cumcount() // window
    loc['window'] = loc['contig_id'].astype(str) + '_' + loc['window'].astype(str)
    df = loc.merge(df_hits, on='#kmer', how='left').fillna({s: 0 for s in samples})
    n_win = df.groupby('window').size()
    df = df.loc[df['window'].map(n_win) >= window // 2]          # drop short contig-end windows
    g = df.groupby('window')
    n = g.size().to_numpy(float)
    cov = (df[samples] > 0).groupby(df['window']).mean()
    dep = df[samples].groupby(df['window']).mean()
    rows = {}
    for smp in samples:
        c, d = cov[smp].to_numpy(float), dep[smp].to_numpy(float)
        p = np.average(c, weights=n)
        var_exp = np.average(p * (1 - p) / n, weights=n)
        var_obs = np.average((c - p) ** 2, weights=n)
        if var_exp > 0:
            disp = var_obs / var_exp
        else:                                                    # coverage 0 or 1 everywhere: perfectly even
            disp = 1.0 if var_obs < 1e-12 else np.nan
        rows[smp] = {'zero_window_frac': (c == 0).mean() if len(c) else np.nan,
                     'window_cov_dispersion': disp,
                     'window_depth_cv': d.std() / d.mean() if d.mean() > 0 else np.nan}
    out = pd.DataFrame.from_dict(rows, orient='index').reset_index(names=['sample'])
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

# colours of the zero-inflation histogram
C_TARGET, C_OTHER, C_EXP0, C_EXC0, C_FIT, C_ALL = '#2a78d6', '#b4b2a9', '#86b6ef', '#e34948', '#0b0b0b', '#52514e'

def zi_traces(x, r, showlegend=True):
    """
    one zero-inflation histogram
    blue  : counts belonging mostly to the target peak (used for λ), grey: other peaks (relative / high-depth / noise)
    zeros : expected from the real peaks (light blue) / excess = zero inflation (red)
    lines : target Poisson (black), all real peaks (grey dashed)
    """
    x = np.asarray(x, dtype=np.int64)
    lam, pi, npos, real, target = r['_lams'], r['_pis'], r['_npos'], r['_real'], r['_target']
    big = real & (pi >= 0.10)                                # small far peaks (high-depth) go to the last bar
    top = lam[big].max() if big.any() else lam[real].max()
    kmax = int(top + 4 * np.sqrt(max(top, 1)) + 2)
    n_over = int((x > kmax).sum())
    f = np.bincount(np.minimum(x, kmax), minlength=kmax + 1)  # counts above kmax lumped in the last bar
    k = np.arange(kmax + 1)

    # expected kmers per peak at each count; at 0 the zeros expected among the present kmers of that peak
    comp = np.array([npos[j] * poisson.pmf(k, lam[j]) / (1 - np.exp(-lam[j])) for j in range(len(lam))])
    comp[:, 0] = npos * np.exp(-lam) / (1 - np.exp(-lam))
    tgt = comp[target].sum(0)
    is_tgt = tgt > 0.5 * comp.sum(0)
    if n_over:
        is_tgt[kmax] = False
    obs_0 = f[0]
    exp_0 = min(r['exp_0_cnts'], obs_0)
    leg = dict(showlegend=showlegend)
    return [go.Bar(x=k, y=np.where(is_tgt & (k > 0), f, 0), marker_color=C_TARGET,
                   name='target peak (used for λ)', legendgroup='target', **leg),
            go.Bar(x=k, y=np.where(~is_tgt & (k > 0), f, 0), marker_color=C_OTHER,
                   name='other peaks (relative / high-depth / noise)', legendgroup='other',
                   customdata=[f'+{n_over} kmers > {kmax}' if (i == kmax and n_over) else '' for i in k],
                   hovertemplate='%{x}: %{y} %{customdata}<extra>other</extra>', **leg),
            go.Bar(x=[0], y=[exp_0], marker_color=C_EXP0, name='0: expected (Poisson)', legendgroup='exp0', **leg),
            go.Bar(x=[0], y=[obs_0 - exp_0], marker_color=C_EXC0, name='0: excess (zero inflation)',
                   legendgroup='exc0', **leg),
            go.Scatter(x=k, y=tgt, mode='lines+markers', line=dict(color=C_FIT, width=2), marker=dict(size=4),
                       name='target Poisson', legendgroup='fit', **leg),
            go.Scatter(x=k, y=comp[real].sum(0), mode='lines', line=dict(color=C_ALL, width=1.5, dash='dash'),
                       name='all peaks', legendgroup='all', **leg)]

def zi_label(r):
    if not np.isfinite(r['lambda']):
        return 'too few positive kmers'
    n = r['n_modes']
    peaks = ', '.join(f'{float(l):.1f} ({float(w):.0%})' for l, w, _ in (p.split(':') for p in r['peaks'].split(';')))
    return (f"<b>target depth λ = {r['lambda']:.2f}</b>  (all present kmers: {r['lambda_present']:.2f})<br>"
            f"{n} depth mode{'s' if n > 1 else ''} | Poisson peaks λ (weight): {peaks}<br>"
            f"zi_frac = {r['zi_frac']:.3f}   zi_index = {r['zi_index']:.2f}")

def plot_zi_distribution(subsets, outdir, prefix, samples=None, quality='both', zi_params=None):
    """
    One figure per sample: zero-inflation histogram with the mixture fit for each kmer subset.
    subsets: dict name -> '#kmer' x samples hit table, e.g. {'all': df_hits, 'rare': ..., 'share': ...}
    """
    subsets = {k: (v.set_index('#kmer') if '#kmer' in v.columns else v) for k, v in subsets.items() if len(v)}
    cols = list(next(iter(subsets.values())).columns)
    samples = cols if samples is None else [s for s in samples if s in cols]

    paths = []
    for smp in samples:
        fits = {name: robust_0_inflation(df[smp].to_numpy(), **(zi_params or ZI_PARAMS)) for name, df in subsets.items()}
        fig = make_subplots(rows=1, cols=len(fits), horizontal_spacing=0.06,
                            subplot_titles=[f'{name} ({len(subsets[name])} kmers)<br><sup>{zi_label(r)}</sup>'
                                            for name, r in fits.items()])
        first = True
        for i, (name, r) in enumerate(fits.items()):
            if not np.isfinite(r['lambda']):
                continue
            for t in zi_traces(subsets[name][smp].to_numpy(), r, showlegend=first):
                fig.add_trace(t, row=1, col=i + 1)
            first = False
        n = len(fits)
        fig.update_layout(barmode='stack', bargap=0.1, template='simple_white', width=max(520 * n, 760),
                          height=520 if n > 1 else 620, title=f'{prefix} | {smp}', margin=dict(t=150, b=40),
                          legend=dict(orientation='h', yanchor='top', y=-0.2 if n > 1 else -0.16, x=0))
        fig.update_xaxes(title_text='kmer hit count')
        fig.update_yaxes(title_text='kmers', col=1)
        fig.update_annotations(font_size=11)
        safe = re.sub(r'[^\w.-]', '_', str(smp))
        base = os.path.join(outdir, f'{prefix}.zi_distribution.{safe}')
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
    parser.add_argument('--export_hist', action='store_true', help = 'export kmer hit count histograms per sample and kmer subset (input for parameter optimisation)')
    parser.add_argument('--zi_k_max', type = int, default = ZI_PARAMS['k_max'], help = 'max number of Poisson components in the mixture')
    parser.add_argument('--zi_min_weight', type = float, default = ZI_PARAMS['min_weight'], help = 'min share of positive kmers for a real peak')
    parser.add_argument('--zi_noise_lam', type = float, default = ZI_PARAMS['noise_lam'], help = 'peaks below this depth are low-count noise if deeper peaks hold the majority')
    parser.add_argument('--zi_dip', type = float, default = ZI_PARAMS['dip'], help = 'min density dip between two maxima to call them separate genomes')
    parser.add_argument('--zi_min_pos', type = int, default = ZI_PARAMS['min_pos'], help = 'min positive kmers to fit')

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    hits_path = args.hits
    counts_path = args.counts
    outpath = args.output_dir
    basename = args.strain if args.strain else 'contig_map'
    zi_params = dict(k_max=args.zi_k_max, min_weight=args.zi_min_weight, noise_lam=args.zi_noise_lam,
                     dip=args.zi_dip, min_pos=args.zi_min_pos)
    with open(os.path.join(outpath, 'zi_params.json'), 'w') as fh:
        json.dump(zi_params, fh, indent=2)


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

    # robust zero inflation: mixture of zero truncated poissons, depth = lowest real peak
    # a subset is only fitted if it has kmers and is not identical to all kmers (e.g. every kmer rare -> no share subset)
    zi_subsets = {'all': df_hits, 'rare': df_hits_unique.reset_index(), 'share': df_hits_shared.reset_index()}
    zi_subsets = {k: v for k, v in zi_subsets.items() if len(v) and (k == 'all' or len(v) < len(df_hits))}
    for sub, df_sub in zi_subsets.items():
        df_zi = zero_inflate_robust(df_sub, **zi_params).set_index('sample').add_prefix(f'{sub}_').reset_index()
        df_cov_depth = df_cov_depth.merge(df_zi, on='sample', how='left')

    # blocks and spatial evenness along the genome (needs --reference for kmer positions)
    if args.reference:
        df_blocks_gpd = df_blocks.groupby('sample', observed=True).agg(**{'block_count': ('block_id', 'count'),
                                                                          '0_cov_blocks': ('block_coverage', lambda x: sum(x==0))}).reset_index()
        df_cov_depth = df_cov_depth.merge(df_blocks_gpd, on='sample', how='left')
        df_cov_depth = df_cov_depth.merge(window_evenness(df_loc_counts, df_hits), on='sample', how='left')

    # strain GC content: a genus-level fingerprint (taxon proxy) -> for stratifying / checking, not as a model feature
    df_cov_depth['strain_gc'] = round(np.mean([(k.upper().count('G') + k.upper().count('C')) / len(k) for k in df_hits['#kmer']]), 4)

    # get relative abundance
    df_cov_depth = df_cov_depth.merge(total_evaluated, on='sample', how='left')

    # relative signature: rare kmers missing while shared kmers are present
    df_cov_depth['percentage_rare'] = df_cov_depth['rare_kmers_total'] / df_cov_depth['all_kmers_total']
    if 'rare' in zi_subsets and 'share' in zi_subsets:
        df_cov_depth['zi_frac_rare_minus_share'] = df_cov_depth['rare_zi_frac'] - df_cov_depth['share_zi_frac']
        # depth of strain-specific vs shared kmers among present kmers: a relative inflates shared depth (ratio < 1)
        df_cov_depth['depth_rare_over_share'] = df_cov_depth['rare_lambda'] / df_cov_depth['share_lambda']

    # relative abundance from the target peak depth (no zeros, no relative / high-depth peaks)
    for sub in zi_subsets:
        if sub == 'share':
            continue
        df_cov_depth[f'relative_abundance_{sub}'] = round((df_cov_depth[f'{sub}_lambda'] * df_cov_depth['total_strain_kmers'] / df_cov_depth['total_evaluated']) * 100, 4)
    df_cov_depth['relative_abundance_present'] = round((df_cov_depth['all_lambda_present'] * df_cov_depth['total_strain_kmers'] / df_cov_depth['total_evaluated']) * 100, 4)


    # Exports and figures
    df_cov_depth['strain'] = args.strain if args.strain else basename
    df_cov_depth.set_index(['strain', 'sample']).round(4).to_csv(os.path.join(outpath,"coverage_depth.tsv"), sep = '\t')
    print(f'Saved coverage dataframe to \n {outpath}')
    print(df_cov_depth.head(5))

    # histograms for parameter optimisation on real data (refit_from_histograms)
    if args.export_hist:
        hist_path = os.path.join(outpath, f'{basename}.zi_hist.tsv.gz')
        hit_histograms(zi_subsets, basename).to_csv(hist_path, sep='\t', index=False)
        print(f'Saved hit histograms to \n {hist_path}')

    if args.figures != 'none':
        plot_zi_distribution(zi_subsets, args.output_dir, basename, quality=args.figures, zi_params=zi_params)
        if args.reference:
            plot_map(df_loc_counts, args.output_dir, basename,
                     df_hits=df_hits, quality=args.figures)





if __name__ == '__main__':
    main()