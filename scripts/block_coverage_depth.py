#!/usr/bin/env python
"""
block_coverage_depth.py

Feature table for hit calling from kmer_block_detect output (the block-aware counterpart of coverage_depth.py).
One row per strain x sample in block_coverage_depth.tsv, so tables from many runs can be concatenated and compared to a
truth table (PR curves).

Inputs (one strain, all samples of a kmer_block_detect run):
  --kmer_hits     <strain>.kmer_hits.tsv.gz      validated fragments per k-mer (+ total_evaluated row)
  --block_hits    <strain>.block_hits.tsv.gz     frags / starts per block and sample
  --blocks        <strain>.blocks.tsv.gz         kmer_scrub_filter.py block table (n_kmers, ref_copies, ...)
  --detect_stats  <strain>.detect_stats.tsv      reads, fragments, accepted / rejected groups
  --kmer_loci     <strain>.kmer_loci.tsv.gz      optional: k-mer -> block, positions (k-mer subsets, window evenness)
  --reference     strain fasta                   optional: genome length for relative abundance
                                                 (else --genome_length, else the last block end, with a warning)

Feature groups (column prefixes):
  reads / QC     reads, read_len, total_bases, frags_with_hits, reject_frac, ambiguous_frac
  blk_           block breadth: blocks hit, fraction hit (all / SNP blocks), frags, starts
  snp_starts_    read starts per SNP block (n_kmers == --snp_kmers, single copy): equal exposure -> Poisson
                   _ztp_lambda    zero-truncated Poisson (positives only)
                   _lambda        robust zero- and right-truncated Poisson (same fit as coverage_depth.py)
                   _zi_frac / _zi_index / _frac_trimmed / _thr
                   _ztnb_mu / _ztnb_size   zero-truncated negative binomial (overdispersion-aware)
                   _zero_lambda (+ _lo / _hi)  inverse Poisson from the fraction of zero blocks
                   _zero_vs_robust             zero-based / robust lambda (1 = plain Poisson)
  best_          one count per placed read, for the block it covers most completely (ties -> closest to its 5' end);
                 exposure per block from the block layout + read lengths (--read_lengths), so counts are Poisson:
                   _rate / _zip_rate / _zip_pi / _zero_rate / _dispersion   all single-copy blocks, exposure-weighted
                   _iso_*   isolated SNP blocks (identical exposure): robust zero inflation, inverse-Poisson lambda
  rep_           repeat families (multi-copy blocks): reads on any copy count once to the family (placed + ambiguous),
                 exposure = sum over copies -> rep_rate, rep_zip_rate / _pi, rep_rate_ratio (vs single-copy rate:
                 copy-number change of mobile elements in the sample strain)
  allunits_      single-copy blocks + repeat families together (exposure-weighted Poisson / ZIP / zero-based)
  diagnostics    best_obs_zero_frac / best_exp_zero_frac / diag_best_zero_excess (empty blocks vs Poisson),
                 best_zero_run_max / best_zero_run_p (do empty blocks cluster along the genome?),
                 best_diag_vs_starts_rate, best_diag_iso_vs_all_rate, best_diag_bigblock_vs_snp_rate (~1 if the
                 exposure model holds), diag_depth_*_vs_best_zip (agreement of the depth estimates),
                 diag_frac_reads_scanned, diag_frac_fragments_hit, strain_filter_* (kmer_scrub_filter.py stats:
                 dropped k-mers per reason, run / segment / block / gap size quantiles, repeat share)
  rate_          start rate per k-mer position by block size class; rate_ratio_long_vs_snp (read rejection QC)
  kmer_<subset>_ coverage, mean depth and robust zero inflation of validated k-mer hits
                 subsets: all, snp (k-mers in SNP blocks), long (k-mers in blocks > --long_kmers)
                 (snp / long need --kmer_loci)
  window evenness (needs --kmer_loci): zero_window_frac, window_cov_dispersion, window_depth_cv
  depth_<est> / relative_abundance_<est>   fold coverage and % of sequenced bases for every estimate
  detection_limit_rel_abund                 95% upper bound (%) when no SNP block has a start
"""

import argparse
import gzip
import json
import os
import re

import numpy as np
import pandas as pd
from scipy.optimize import brentq, minimize
from scipy.stats import nbinom, poisson


### robust zero inflation: zero- and right-truncated Poisson (identical to coverage_depth.py)
ZI_PARAMS = dict(q=0.95, min_pos=10)
COUNT_COLS = ['pangenome_count', 'metagenome_count', 'drug_count']


def zero_inflate_stats(obs_0, total, lam, exp0):
    excess_0 = obs_0 - exp0
    return {'lambda': lam,
            'zi_frac': np.clip(excess_0 / total, 0, 1) if np.isfinite(exp0) else np.nan,
            'zi_index': excess_0 / obs_0 if obs_0 > 0 else 0.0,
            'exp_0_cnts': exp0}


def trunc_poisson_lambda(f, thr):
    k = np.arange(len(f))
    keep = (k >= 1) & (k <= thr)
    n_kept = f[keep].sum()
    if n_kept == 0:
        return np.nan, 0
    m = (k[keep] * f[keep]).sum() / n_kept
    if m <= 1 + 1e-6:
        return 0.0, n_kept
    g = lambda l: l * poisson.cdf(thr - 1, l) / (poisson.cdf(thr, l) - np.exp(-l)) - m
    lo, hi = 1e-9, max(2 * m, 10.0)
    while g(hi) < 0 and hi < 1e7:
        hi *= 2
    if not (g(lo) < 0 <= g(hi)):
        return np.nan, n_kept
    return brentq(g, lo, hi), n_kept


def robust_0_inflation(x, q=0.95, min_pos=10, n_iter=30):
    f = np.bincount(np.asarray(x, dtype=np.int64), minlength=3)
    total, obs_0 = f.sum(), f[0]
    pos = total - obs_0
    nan = {**zero_inflate_stats(obs_0, total, np.nan, np.nan), 'zi_index': np.nan, 'frac_trimmed': np.nan, 'thr': np.nan}
    if pos < min_pos:
        return nan
    k = np.arange(len(f))

    def fit_from(lam):
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
        return lam, f[1:thr + 1].sum() / pos

    pos_med = np.median(np.repeat(k[1:], f[1:]))
    lam, share = fit_from(2.0 * f[2] / f[1]) if f[1] > 0 and f[2] > 0 else fit_from(pos_med)
    if share < 0.5:
        lam_m, share_m = fit_from(pos_med)
        if share_m > share:
            lam = lam_m

    thr = int(max(poisson.ppf(q, lam), 2)) if lam > 0 else 2
    lam_fit, n_kept = trunc_poisson_lambda(f, thr)
    lam = lam_fit if np.isfinite(lam_fit) else lam
    frac_trimmed = f[thr + 1:].sum() / pos
    if not np.isfinite(lam) or lam <= 0:
        return {**nan, 'frac_trimmed': frac_trimmed, 'thr': thr}
    exp0 = n_kept / (poisson.cdf(thr, lam) - np.exp(-lam)) * np.exp(-lam)
    return {**zero_inflate_stats(obs_0, total, lam, exp0), 'frac_trimmed': frac_trimmed, 'thr': thr}


def robust_0_inflation_safe(x, label='', **kw):
    try:
        return robust_0_inflation(x, **kw)
    except Exception as e:
        f = np.bincount(np.asarray(x, dtype=np.int64))
        k = np.nonzero(f)[0]
        print(f'WARNING: zero-inflation fit failed for {label}: {type(e).__name__}: {e}')
        print(f'  histogram count:n {dict(zip(k.tolist(), f[k].tolist()))}')
        return robust_0_inflation(np.zeros(0, np.int64), **kw)


### simple estimators for equal-exposure counts (read starts per SNP block)
def ztp_lambda(x, n_iter=200):
    """zero-truncated Poisson MLE: mean of positives = lam / (1 - exp(-lam))"""
    x = np.asarray(x)
    pos = x[x > 0]
    if len(pos) == 0:
        return 0.0
    m = pos.mean()
    if m <= 1:
        return np.nan                     # all positives are 1s: not identifiable
    lam = m
    for _ in range(n_iter):
        lam = m * (1 - np.exp(-lam))
    return lam


def ztnb_fit(x, min_pos=10):
    """zero-truncated negative binomial: returns mean mu and size (dispersion; large = Poisson)"""
    x = np.asarray(x)
    pos = x[x > 0]
    if len(pos) < min_pos or pos.mean() <= 1:
        return np.nan, np.nan

    def nll(params):
        mu, size = np.exp(params)
        p = size / (size + mu)
        return -(nbinom.logpmf(pos, size, p) - np.log1p(-nbinom.pmf(0, size, p))).sum()

    fit = minimize(nll, x0=[np.log(pos.mean()), 0.0], method='Nelder-Mead',
                   options=dict(xatol=1e-6, fatol=1e-8, maxiter=2000))
    mu, size = np.exp(fit.x)
    return (mu, size) if fit.success else (np.nan, np.nan)


def lambda_from_zeros(x, z=1.96):
    """inverse Poisson: lam = -ln(fraction of zeros); 95% CI from the binomial error; upper bound if all zero"""
    x = np.asarray(x)
    n = len(x)
    if n == 0:
        return np.nan, np.nan, np.nan
    p0 = (x == 0).mean()
    if p0 == 1:
        return 0.0, 0.0, -np.log(0.05) / n
    if p0 == 0:
        return np.inf, np.nan, np.inf
    se = np.sqrt(p0 * (1 - p0) / n)
    return -np.log(p0), -np.log(min(p0 + z * se, 1.0)), -np.log(max(p0 - z * se, 1e-12))


def rate_from_zeros(starts, n_kmers, n_iter=100):
    """unequal exposure: solve sum_b exp(-r * n_b) = observed zero blocks"""
    starts = np.asarray(starts)
    n_kmers = np.asarray(n_kmers, float)
    obs = (starts == 0).sum()
    if obs == len(starts):
        return 0.0
    if obs == 0:
        return np.inf
    lo, hi = 0.0, 1.0
    while np.exp(-hi * n_kmers).sum() > obs:
        hi *= 2
    for _ in range(n_iter):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if np.exp(-mid * n_kmers).sum() > obs else (lo, mid)
    return (lo + hi) / 2


### exposure of the 'best' count (one block per read)
# Every placed read is counted once: single-copy blocks first, then the block it covers most completely (largest
# share of the block's k-mers), ties -> the block closest to the read's 5' end. That choice depends only on where the read starts,
# its length and strand, so per-block counts are independent Poisson(rate x exposure_b); exposure_b (start
# positions) is computed by sweeping read starts over the block layout for every read length bin.
def _assign(gs, ge, nk, lo_ok, hi_ok, w, strand, uniq):
    """
    For reads whose k-mer starts span [p, p + w] (w = L - k) on one strand: which block does each start
    position p give the read to? Returns the number of start positions per block (len(gs) array).
    gs / ge: first / last k-mer start of each block on a global axis; lo_ok / hi_ok: valid p per block's contig.
    """
    n_rows = (ge - gs + 1) + w                                   # p from gs - w to ge
    b_idx = np.repeat(np.arange(len(gs)), n_rows)
    off = np.arange(n_rows.sum()) - np.repeat(np.cumsum(n_rows) - n_rows, n_rows)
    p = np.repeat(gs - w, n_rows) + off
    keep = (p >= lo_ok[b_idx]) & (p <= hi_ok[b_idx])
    p, b_idx = p[keep], b_idx[keep]
    s_b, e_b = gs[b_idx], ge[b_idx]
    overlap = np.minimum(e_b, p + w) - np.maximum(s_b, p) + 1
    frac = overlap / nk[b_idx]
    dist = (np.maximum(s_b, p) - p) if strand > 0 else ((p + w) - np.minimum(e_b, p + w))
    # by p, then single-copy blocks first, then most complete, then closest to the 5' end (as kmer_block_detect)
    order = np.lexsort((dist, -frac, -uniq[b_idx], p))
    p, b_idx = p[order], b_idx[order]
    first = np.ones(len(p), bool)
    first[1:] = p[1:] != p[:-1]
    return np.bincount(b_idx[first], minlength=len(gs))


def exposure_by_length(df_blocks, contig_len, lengths, k):
    """
    exposure (start positions) per block for each read length in `lengths`, averaged over both strands.
    df_blocks: block_id, contig_id, start, stop (0-based, stop = end of the last k-mer), n_kmers
    returns DataFrame index block_id, one column per length
    """
    b = df_blocks.sort_values(['contig_id', 'start']).reset_index(drop=True)
    contigs = sorted(b['contig_id'].unique())
    clen = {c: int(contig_len.get(c, b.loc[b['contig_id'] == c, 'stop'].max() + 1)) for c in contigs}
    pad = 10 * max(max(lengths), 1000)
    offset, o = {}, 0
    for c in contigs:
        offset[c] = o
        o += clen[c] + pad
    off = b['contig_id'].map(offset).to_numpy(np.int64)
    gs = b['start'].to_numpy(np.int64) + off
    ge = (b['stop'].to_numpy(np.int64) - k + 1) + off
    nk = (ge - gs + 1).astype(float)
    cl = b['contig_id'].map(clen).to_numpy(np.int64)
    uniq = (b['ref_copies'].to_numpy() == 1).astype(int) if 'ref_copies' in b.columns else np.ones(len(b), int)
    out = {}
    for L in lengths:
        w = int(L) - k
        if w < 0:
            continue
        lo_ok, hi_ok = off, off + cl - k - w                     # read must fit on its contig
        e = (_assign(gs, ge, nk, lo_ok, hi_ok, w, +1, uniq) + _assign(gs, ge, nk, lo_ok, hi_ok, w, -1, uniq)) / 2
        out[int(L)] = e
    return pd.DataFrame(out, index=b['block_id'].to_numpy())


def length_bins(df_len, k, width=10):
    """read-length histogram (sample, length, n) -> bin representative lengths and per-sample weights"""
    d = df_len.loc[df_len['length'] >= k].copy()
    d['bin'] = (d['length'] - k) // width
    rep = d.groupby('bin').apply(lambda x: int(round(np.average(x['length'], weights=x['n']))),
                                 include_groups=False)
    d['rep'] = d['bin'].map(rep)
    w = d.groupby(['sample', 'rep'])['n'].sum().unstack(fill_value=0)
    w = w.div(w.sum(axis=1), axis=0)
    mean_len = d.groupby('sample').apply(lambda x: np.average(x['length'], weights=x['n']), include_groups=False)
    return sorted(rep.unique()), w, mean_len


def sample_exposure(exp_by_len, weights):
    """per-sample exposure per block: sum over length bins of weight x exposure(L). returns blocks x samples"""
    cols = [c for c in weights.columns if c in exp_by_len.columns]
    return exp_by_len[cols].to_numpy() @ weights[cols].T.to_numpy() if cols else None


def zip_rate(y, e, n_iter=500, tol=1e-10, min_pos=10):
    """zero-inflated Poisson with exposure: y_b ~ pi*0 + (1-pi)*Poisson(r * e_b). EM; returns r, pi"""
    y = np.asarray(y, float); e = np.asarray(e, float)
    if y.sum() == 0:
        return 0.0, np.nan
    if (y > 0).sum() < min_pos:                    # too few positive units: pi not identifiable
        return y.sum() / e.sum(), np.nan
    zero = y == 0
    r = y.sum() / e.sum()
    # start pi at the zero excess (pi = 0 would be a fixed point of the EM and never move)
    pi = float(np.clip(zero.mean() - np.exp(-r * e).mean(), 0.05, 0.95))
    for _ in range(n_iter):
        p0 = np.exp(-r * e)
        z = np.where(zero, pi / (pi + (1 - pi) * p0), 0.0)
        pi_new = z.mean()
        r_new = y.sum() / ((1 - z) * e).sum()
        done = abs(pi_new - pi) < tol and abs(r_new - r) < tol * max(r, 1e-12)
        r, pi = r_new, pi_new
        if done:
            break
    return r, pi



### spatial evenness (identical to coverage_depth.py)
def _scaled(df, col, scale):
    """df[col] * scale, or None if the column does not exist (estimator not available for this strain)"""
    return df[col] * scale if col in df.columns else None


def window_evenness(df_loc, df_hits, window=200):
    samples = [c for c in df_hits.columns if c != '#kmer']
    loc = df_loc[['#kmer', 'contig_id', 'kmer_position']].sort_values(['contig_id', 'kmer_position']).copy()
    loc['window'] = loc.groupby('contig_id').cumcount() // window
    loc['window'] = loc['contig_id'].astype(str) + '_' + loc['window'].astype(str)
    df = loc.merge(df_hits, on='#kmer', how='left').fillna({s: 0 for s in samples})
    n_win = df.groupby('window').size()
    df = df.loc[df['window'].map(n_win) >= window // 2]
    if df.empty:
        # few rare k-mers spread over many contigs: no window reaches window // 2 k-mers -> no evenness metrics
        print(f"NOTE: no window with >= {window // 2} k-mers; window_* metrics left empty")
        return pd.DataFrame({'sample': samples, 'zero_window_frac': np.nan,
                             'window_cov_dispersion': np.nan, 'window_depth_cv': np.nan})
    n = df.groupby('window').size().to_numpy(float)
    cov = (df[samples] > 0).groupby(df['window']).mean()
    dep = df[samples].groupby(df['window']).mean()
    rows = {}
    for smp in samples:
        c, d = cov[smp].to_numpy(float), dep[smp].to_numpy(float)
        p = np.average(c, weights=n)
        var_exp = np.average(p * (1 - p) / n, weights=n)
        var_obs = np.average((c - p) ** 2, weights=n)
        disp = var_obs / var_exp if var_exp > 0 else (1.0 if var_obs < 1e-12 else np.nan)
        rows[smp] = {'zero_window_frac': (c == 0).mean() if len(c) else np.nan,
                     'window_cov_dispersion': disp,
                     'window_depth_cv': d.std() / d.mean() if d.mean() > 0 else np.nan}
    return pd.DataFrame.from_dict(rows, orient='index').reset_index(names=['sample'])


### plots
# colours of the zero-inflation histogram (same as coverage_depth.py)
C_BODY, C_TRUNC, C_EXP0, C_EXC0, C_FIT = '#2a78d6', '#b4b2a9', '#86b6ef', '#e34948', '#0b0b0b'


def _fig_formats(quality):
    return ['svg', 'png'] if quality == 'both' else [quality]


def zi_traces(x, r, showlegend=True):
    """one zero-inflation histogram: Poisson body (blue), trimmed tail (grey), expected / excess zeros, fitted Poisson"""
    import plotly.graph_objects as go
    x = np.asarray(x, dtype=np.int64)
    lam, thr = r['lambda'], int(r['thr'])
    kmax = int(max(thr + 3 * np.sqrt(max(lam, 1)) + 2, thr + 3))
    n_over = int((x > kmax).sum())
    f = np.bincount(np.minimum(x, kmax), minlength=kmax + 1)
    k = np.arange(kmax + 1)
    body = (k >= 1) & (k <= thr)
    obs_0 = f[0]
    exp_0 = min(r['exp_0_cnts'], obs_0)
    body_present = f[body].sum() / (poisson.cdf(thr, lam) - np.exp(-lam))
    leg = dict(showlegend=showlegend)
    return [go.Bar(x=k, y=np.where(body, f, 0), marker_color=C_BODY, name='Poisson body 1..thr (used for λ)',
                   legendgroup='body', **leg),
            go.Bar(x=k, y=np.where(k > thr, f, 0), marker_color=C_TRUNC, name='high counts > thr (present, trimmed)',
                   legendgroup='trunc',
                   customdata=[f'+{n_over} > {kmax}' if (i == kmax and n_over) else '' for i in k],
                   hovertemplate='%{x}: %{y} %{customdata}<extra>trimmed</extra>', **leg),
            go.Bar(x=[0], y=[exp_0], marker_color=C_EXP0, name='0: expected (Poisson)', legendgroup='exp0', **leg),
            go.Bar(x=[0], y=[obs_0 - exp_0], marker_color=C_EXC0, name='0: excess (zero inflation)',
                   legendgroup='exc0', **leg),
            go.Scatter(x=k, y=body_present * poisson.pmf(k, lam), mode='lines+markers',
                       line=dict(color=C_FIT, width=2), marker=dict(size=4), name='fitted Poisson',
                       legendgroup='fit', **leg)]


def zi_label(r):
    if not np.isfinite(r['lambda']):
        return 'λ not identifiable (too few positives or all positives are 1s)'
    return (f"<b>λ = {r['lambda']:.2f}</b>  (body 1..{int(r['thr'])}, {r['frac_trimmed']:.1%} of positives trimmed)<br>"
            f"zi_frac = {r['zi_frac']:.3f}   zi_index = {r['zi_index']:.2f}")


def plot_zi_distribution(subsets, outdir, prefix, samples, quality, zi_params):
    """
    one figure per sample, one panel per count subset (x = count, y = number of blocks / k-mers)
    subsets: name -> DataFrame (rows = blocks or k-mers, columns = samples)
    """
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    paths = []
    for smp in samples:
        fits = {name: robust_0_inflation_safe(df[smp].to_numpy(), label=f'{prefix} {name} {smp}', **zi_params)
                for name, df in subsets.items() if smp in df.columns}
        fig = make_subplots(rows=1, cols=len(fits), horizontal_spacing=0.06,
                            subplot_titles=[f'{name} ({len(subsets[name])})<br><sup>{zi_label(r)}</sup>'
                                            for name, r in fits.items()])
        first = True
        for i, (name, r) in enumerate(fits.items()):
            if not np.isfinite(r['lambda']):
                continue
            for t in zi_traces(subsets[name][smp].to_numpy(), r, showlegend=first):
                fig.add_trace(t, row=1, col=i + 1)
            first = False
        n = len(fits)
        fig.update_layout(barmode='stack', bargap=0.1, template='simple_white', width=max(480 * n, 760),
                          height=520, title=f'{prefix} | {smp}', margin=dict(t=150, b=40),
                          legend=dict(orientation='h', yanchor='top', y=-0.2, x=0))
        fig.update_xaxes(title_text='count (read starts per SNP block / fragments per k-mer)')
        fig.update_yaxes(title_text='blocks / k-mers', col=1)
        fig.update_annotations(font_size=11)
        safe = re.sub(r'[^\w.-]', '_', str(smp))
        base = os.path.join(outdir, f'{prefix}.zi_distribution.{safe}')
        for fmt in _fig_formats(quality):
            path = f'{base}.{fmt}'
            fig.write_image(path, scale=2 if fmt == 'png' else None)
            print(f'saved {path}')
            paths.append(path)
    return paths


def plot_block_map(df_blocks, df_bl, contig_len, outdir, prefix, samples, quality,
                   hit_cmap='viridis', hit_q=0.995, min_lw=0.25):
    """
    Circular genome maps with every block drawn as a bar of its true genomic width (start..stop).
    Bars narrower than a pixel are kept visible with a thin edge of the same colour (min_lw).

    strain map  ({prefix}.block_map.<fmt>):
        contigs | blocks coloured by scrub count (pangenome + metagenome + drug) | repeat blocks (ref_copies > 1)
    per sample ({prefix}.block_map.<sample>.<fmt>):
        contigs | blocks by scrub count | fragments per block (bar height = log10(1 + frags), colour = frags)
                | read starts per block (bar height) | blocks without any fragment (red)
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.cm import ScalarMappable
    from matplotlib.patches import Patch
    from pycirclize import Circos
    import plotly.express as px

    contig_len = dict(contig_len)
    missing = set(df_blocks['contig_id']) - set(contig_len)
    if missing:
        print(f'WARNING: {len(missing)} block contigs not in the reference headers; using their last block end as length')
        ends = df_blocks.groupby('contig_id')['stop'].max() + 1
        contig_len.update({c: int(ends[c]) for c in missing})
    contigs = (pd.DataFrame({'contig_id': list(contig_len), 'contig_length': list(contig_len.values())})
               .sort_values('contig_length', ascending=False))
    contigs['offset'] = contigs['contig_length'].cumsum() - contigs['contig_length']
    genome_len = int(contigs['contig_length'].sum())

    blk = df_blocks.merge(contigs[['contig_id', 'offset']], on='contig_id', how='left')
    blk['g0'] = blk['start'] + blk['offset']
    blk['g1'] = blk['stop'] + 1 + blk['offset']
    blk['count'] = blk[COUNT_COLS].fillna(0).sum(axis=1).astype(int)
    blk = blk.sort_values('g0').reset_index(drop=True)

    palette = px.colors.qualitative.Dark24
    present = sorted(blk['count'].unique())
    count_colors = {c: palette[i % len(palette)] for i, c in enumerate(present)}

    def rect(track, g0, g1, color, frac=1.0):
        r0, r1 = track.r_lim
        track.rect(g0, g1, r_lim=(r0, r0 + (r1 - r0) * frac), fc=color, ec=color, lw=min_lw)

    def base(radii):
        circos = Circos({'genome': genome_len}, start=0, end=360)
        sector = circos.sectors[0]
        t_ctg = sector.add_track(radii['ctg'])
        for i, (st, ln) in enumerate(zip(contigs['offset'], contigs['contig_length'])):
            t_ctg.rect(st, st + ln, fc='grey' if i % 2 else 'lightgrey', ec='none')
        t_cnt = sector.add_track(radii['cnt'])
        for c, d in blk.groupby('count'):
            for g0, g1 in zip(d['g0'], d['g1']):
                rect(t_cnt, g0, g1, count_colors[c])
        return circos, sector

    def save(circos, handles, title, name, cbar=None):
        fig = circos.plotfig()
        circos.ax.set_title(f'{prefix} | {name}', fontsize=11, pad=12)
        circos.ax.legend(handles=handles, title=title, loc='upper left', bbox_to_anchor=(1.02, 0.95),
                         frameon=False, fontsize=8, title_fontsize=9)
        if cbar is not None:
            norm, cmap, label = cbar
            cax = circos.ax.inset_axes([1.05, 0.05, 0.025, 0.35])
            cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax)
            cb.set_label(label, fontsize=8)
            cb.ax.tick_params(labelsize=7)
        out = os.path.join(outdir, f'{prefix}.block_map{"" if name == "strain" else "." + name}')
        for fmt in _fig_formats(quality):
            kw = dict(format='svg') if fmt == 'svg' else dict(dpi=300)
            fig.savefig(f'{out}.{fmt}', bbox_inches='tight', **kw)
            print(f'saved {out}.{fmt}')
        plt.close(fig)

    count_handles = [Patch(facecolor=count_colors[c], label=f'scrub count {c}') for c in present]

    # ── strain map ──
    circos, sector = base({'ctg': (88, 92), 'cnt': (76, 84)})
    t_rep = sector.add_track((68, 72))
    for g0, g1 in zip(blk.loc[blk['ref_copies'] > 1, 'g0'], blk.loc[blk['ref_copies'] > 1, 'g1']):
        rect(t_rep, g0, g1, 'darkorange')
    save(circos, count_handles + [Patch(facecolor='darkorange', label='repeat block (ref_copies > 1)')],
         'blocks', 'strain')

    # ── per sample ──
    cmap = plt.get_cmap(hit_cmap)
    for smp in samples:
        d = blk.merge(df_bl.loc[df_bl['sample'] == smp, ['block_id', 'frags', 'starts']], on='block_id', how='left')
        d[['frags', 'starts']] = d[['frags', 'starts']].fillna(0)
        nz = d['frags'] > 0
        fmax = max(2.0, float(np.quantile(d.loc[nz, 'frags'], hit_q))) if nz.any() else 2.0
        smax = max(1.0, float(np.quantile(d.loc[d['starts'] > 0, 'starts'], hit_q))) if (d['starts'] > 0).any() else 1.0
        norm = LogNorm(vmin=1, vmax=fmax)

        circos, sector = base({'ctg': (92, 95), 'cnt': (86, 90)})
        t_frg = sector.add_track((62, 84))
        t_str = sector.add_track((44, 60))
        t_zero = sector.add_track((38, 42))
        t_frg.axis(fc='none', ec='lightgrey', lw=0.3)
        t_str.axis(fc='none', ec='lightgrey', lw=0.3)
        for g0, g1, f in zip(d.loc[nz, 'g0'], d.loc[nz, 'g1'], d.loc[nz, 'frags']):
            frac = np.log10(1 + min(f, fmax)) / np.log10(1 + fmax)
            rect(t_frg, g0, g1, cmap(norm(min(f, fmax))), frac=max(frac, 0.03))
        for g0, g1, s_ in zip(d.loc[d['starts'] > 0, 'g0'], d.loc[d['starts'] > 0, 'g1'], d.loc[d['starts'] > 0, 'starts']):
            rect(t_str, g0, g1, 'steelblue', frac=max(min(s_, smax) / smax, 0.03))
        for g0, g1 in zip(d.loc[~nz, 'g0'], d.loc[~nz, 'g1']):
            rect(t_zero, g0, g1, '#e34948')

        handles = count_handles + [Patch(facecolor='steelblue', label=f'read starts (bar height, max {smax:.0f})'),
                                   Patch(facecolor='#e34948', label=f'0 fragments ({int((~nz).sum())} blocks)')]
        safe = re.sub(r'[^\w.-]', '_', str(smp))
        save(circos, handles, 'blocks', safe, cbar=(norm, cmap, 'fragments per block (log)'))


### diagnostics helpers
def _max_run(flags):
    best = cur = 0
    for v in flags:
        cur = cur + 1 if v else 0
        best = max(best, cur)
    return best


def _max_run_rows(z):
    """longest run of True per row of a 2-D boolean array"""
    z = np.asarray(z, bool)
    pad = np.zeros((z.shape[0], 1), bool)
    zz = np.concatenate([pad, z, pad], axis=1).astype(np.int8)
    d = np.diff(zz, axis=1)
    out = np.zeros(z.shape[0], int)
    for i in range(z.shape[0]):
        st, en = np.flatnonzero(d[i] == 1), np.flatnonzero(d[i] == -1)
        out[i] = (en - st).max() if len(st) else 0
    return out


def zero_run_test(zero_flags, p0, n_sim=200, seed=0):
    """longest run of consecutive empty blocks along the genome, and the share of simulations with a run at
    least as long when every block is empty with its own Poisson probability p0 = exp(-rate * exposure)
    (small p = empty blocks cluster beyond what their exposures explain: a region missing from the strain).
    Blocks with low exposure (back-to-back chunks) are empty more often, so a plain permutation would mislead."""
    z = np.asarray(zero_flags, bool)
    if len(z) == 0 or not z.any() or z.all():
        return _max_run(z), np.nan
    obs = _max_run(z)
    rng = np.random.default_rng(seed)
    sims = _max_run_rows(rng.random((n_sim, len(z))) < np.asarray(p0, float)[None, :])
    return obs, (int((sims >= obs).sum()) + 1) / (n_sim + 1)


### helpers
def contig_lengths(path):
    opener = gzip.open if path.endswith('.gz') else open
    lens, name = {}, None
    with opener(path, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                name = line[1:].split()[0]
                lens[name] = 0
            elif name is not None:
                lens[name] += len(line.strip())
    return lens


def genome_length(path):
    opener = gzip.open if path.endswith('.gz') else open
    n = 0
    with opener(path, 'rt') as fh:
        for line in fh:
            if not line.startswith('>'):
                n += len(line.strip())
    return n


def hist_rows(strain, sample, subset, x):
    f = np.bincount(np.asarray(x, dtype=np.int64))
    k = np.nonzero(f)[0]
    return pd.DataFrame({'strain': strain, 'sample': sample, 'subset': subset, 'subset_n': len(x),
                         'count': k, 'n': f[k]})


def kmer_subset_features(df_hits, subset, zi_params, label):
    """coverage, mean depth and robust zero inflation per sample for one k-mer subset"""
    df = df_hits.set_index('#kmer')
    rows = {}
    for smp in df.columns:
        x = df[smp].to_numpy()
        r = robust_0_inflation_safe(x, label=f'{label} kmer_{subset} {smp}', **zi_params)
        rows[smp] = {'n': len(x),
                     'coverage': (x > 0).mean() if len(x) else np.nan,
                     'depth_mean': x.mean() if len(x) else np.nan,
                     'lambda': r['lambda'], 'zi_frac': r['zi_frac'], 'zi_index': r['zi_index'],
                     'frac_trimmed': r['frac_trimmed'],
                     'zero_lambda': lambda_from_zeros(x)[0]}
    return (pd.DataFrame.from_dict(rows, orient='index').add_prefix(f'kmer_{subset}_')
            .reset_index(names=['sample']))


def main():
    parser = argparse.ArgumentParser(description='Hit-calling features from kmer_block_detect output (one strain, all samples)')
    parser.add_argument('--kmer_hits', required=True, help='<strain>.kmer_hits.tsv.gz from kmer_block_detect')
    parser.add_argument('--block_hits', required=True, help='<strain>.block_hits.tsv.gz from kmer_block_detect')
    parser.add_argument('--blocks', required=True, help='<strain>.blocks.tsv.gz from kmer_scrub_filter.py --reference')
    parser.add_argument('--detect_stats', required=True, help='<strain>.detect_stats.tsv from kmer_block_detect')
    parser.add_argument('--read_lengths', help='<strain>.read_lengths.tsv.gz from kmer_block_detect (exposure of the '
                        "'best' count); default: next to --block_hits")
    parser.add_argument('--family_hits', help='<strain>.family_hits.tsv.gz from kmer_block_detect (repeat families); '
                        'default: next to --block_hits')
    parser.add_argument('--selection', help='<strain>.selection.json from kmer_scrub_filter.py (filter diagnostics, '
                        'copied into strain_filter_* columns); default: next to --blocks')
    parser.add_argument('--kmer_loci', help='<strain>.kmer_loci.tsv.gz: enables snp/long k-mer subsets and window evenness')
    parser.add_argument('--reference', '-r', help='strain fasta, for the genome length')
    parser.add_argument('--genome_length', type=int, help='genome length in bp if no --reference is given')
    parser.add_argument('--strain', required=True, help='strain name (index of the output table)')
    parser.add_argument('--output_dir', required=True)
    parser.add_argument('--snp_kmers', type=int, default=31, help='block size of a single-SNP block (= k); equal-exposure set for the start model')
    parser.add_argument('--long_kmers', type=int, default=120, help='blocks above this many k-mers are "long" (rejection QC, kmer_long subset)')
    parser.add_argument('--zi_q', type=float, default=ZI_PARAMS['q'])
    parser.add_argument('--zi_min_pos', type=int, default=ZI_PARAMS['min_pos'])
    parser.add_argument('--figures', choices=['png', 'svg', 'both', 'none'], default='none',
                        help='save figures: zero-inflation histograms per sample and circular block maps')
    parser.add_argument('--min_iso_blocks', type=int, default=30,
                        help='minimum isolated SNP blocks for the best_iso_* fits (default 30)')
    parser.add_argument('--export_blocks', action='store_true',
                        help='export the per-block table (position, exposure, frags, starts, best per sample) '
                             'as <strain>.blocks_long.tsv.gz, for debugging coverage along the genome')
    parser.add_argument('--export_hist', action='store_true',
                        help='export count histograms (SNP-block starts, k-mer subsets) for refitting / parameter optimisation')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    strain = args.strain
    zi_params = dict(q=args.zi_q, min_pos=args.zi_min_pos)
    with open(os.path.join(args.output_dir, 'zi_params.json'), 'w') as fh:
        json.dump({**zi_params, 'snp_kmers': args.snp_kmers, 'long_kmers': args.long_kmers}, fh, indent=2)

    df_kmer_hits = pd.read_csv(args.kmer_hits, sep='\t')
    df_block_hits = pd.read_csv(args.block_hits, sep='\t')
    df_blocks = pd.read_csv(args.blocks, sep='\t')
    df_stats = pd.read_csv(args.detect_stats, sep='\t')
    df_block_hits['sample'] = df_block_hits['sample'].astype(str)
    df_stats['sample'] = df_stats['sample'].astype(str)

    k = len(df_kmer_hits.loc[df_kmer_hits['#kmer'] != 'total_evaluated', '#kmer'].iloc[0])
    snp_n = args.snp_kmers

    contig_len = None
    if args.reference:
        contig_len = contig_lengths(args.reference)
        genome_len = int(sum(contig_len.values()))
    elif args.genome_length:
        genome_len = args.genome_length
    else:
        genome_len = int((df_blocks.groupby('contig_id')['stop'].max() + 1).sum())
        print(f'WARNING: no --reference / --genome_length; using the last block end per contig ({genome_len} bp), '
              f'relative abundances are slightly too low')

    # total evaluated row
    total_evaluated = df_kmer_hits.loc[df_kmer_hits['#kmer'] == 'total_evaluated'].drop(columns='#kmer').T.reset_index()
    total_evaluated.columns = ['sample', 'total_evaluated']
    total_evaluated['sample'] = total_evaluated['sample'].astype(str)
    df_kmer_hits = df_kmer_hits.loc[df_kmer_hits['#kmer'] != 'total_evaluated'].copy()
    df_kmer_hits.columns = ['#kmer'] + [str(c) for c in df_kmer_hits.columns[1:]]
    samples = [c for c in df_kmer_hits.columns if c != '#kmer']

    # ── reads / QC ──
    df_cov = df_stats.copy()
    df_cov['read_len'] = df_cov['kmer_positions'] / df_cov['reads'] + k - 1
    df_cov['total_bases'] = df_cov['kmer_positions'] + df_cov['reads'] * (k - 1)
    groups = df_cov['groups_accepted'] + df_cov['groups_rejected']
    df_cov['reject_frac'] = np.where(groups > 0, df_cov['groups_rejected'] / groups.where(groups > 0, 1), np.nan)
    pieces = df_cov['pieces_placed'] + df_cov['pieces_ambiguous']
    df_cov['ambiguous_frac'] = np.where(pieces > 0, df_cov['pieces_ambiguous'] / pieces.where(pieces > 0, 1), np.nan)
    df_cov = df_cov.merge(total_evaluated, on='sample', how='left')

    # ── blocks ──
    df_bl = df_block_hits.merge(df_blocks, on='block_id', how='left')
    single = df_bl['ref_copies'] == 1
    df_bl['is_snp'] = single & (df_bl['n_kmers'] == snp_n)
    df_bl['is_long'] = single & (df_bl['n_kmers'] > args.long_kmers)

    df_blk = (
        df_bl.groupby('sample')
        .agg(**{'blk_n': ('block_id', 'size'),
                'blk_hit': ('frags', lambda x: (x > 0).sum()),
                'blk_frac_hit': ('frags', lambda x: (x > 0).mean()),
                'blk_frags_total': ('frags', 'sum'),
                'blk_frags_mean': ('frags', 'mean'),
                'blk_starts_total': ('starts', 'sum')})
        .reset_index()
    )
    df_blk_snp = (
        df_bl.loc[df_bl['is_snp']].groupby('sample')
        .agg(**{'blk_snp_n': ('block_id', 'size'),
                'blk_snp_frac_hit': ('frags', lambda x: (x > 0).mean()),
                'blk_snp_frags_mean': ('frags', 'mean')})
        .reset_index()
    )
    df_cov = df_cov.merge(df_blk, on='sample', how='left').merge(df_blk_snp, on='sample', how='left')

    # start rates by block size class (single copy)
    df_rate = (
        df_bl.loc[single]
        .assign(size_class=lambda d: np.select([d['is_snp'], d['is_long']], ['snp', 'long'], 'mid'))
        .groupby(['sample', 'size_class'])
        .agg(starts=('starts', 'sum'), kmers=('n_kmers', 'sum'))
        .assign(rate=lambda d: d['starts'] / d['kmers'])['rate']
        .unstack()
        .reindex(columns=['snp', 'mid', 'long'])
        .add_prefix('rate_')
        .reset_index()
    )
    df_rate['rate_ratio_long_vs_snp'] = df_rate['rate_long'] / df_rate['rate_snp'].replace(0, np.nan)
    df_rate_all0 = (
        df_bl.loc[single].groupby('sample')
        .apply(lambda d: rate_from_zeros(d['starts'], d['n_kmers']), include_groups=False)
        .rename('rate_zero_allblocks').reset_index()
    )
    df_cov = df_cov.merge(df_rate, on='sample', how='left').merge(df_rate_all0, on='sample', how='left')

    # ── read starts on SNP blocks: equal exposure -> Poisson family fits ──
    hists = []
    rows = {}
    for smp, d in df_bl.loc[df_bl['is_snp']].groupby('sample'):
        x = d['starts'].to_numpy()
        r = robust_0_inflation_safe(x, label=f'{strain} snp_starts {smp}', **zi_params)
        lam0, lam0_lo, lam0_hi = lambda_from_zeros(x)
        mu_nb, size_nb = ztnb_fit(x, min_pos=zi_params['min_pos'])
        rows[smp] = {'n': len(x), 'mean': x.mean(), 'obs_zero': (x == 0).mean(),
                     'ztp_lambda': ztp_lambda(x),
                     'lambda': r['lambda'], 'zi_frac': r['zi_frac'], 'zi_index': r['zi_index'],
                     'frac_trimmed': r['frac_trimmed'], 'thr': r['thr'],
                     'ztnb_mu': mu_nb, 'ztnb_size': size_nb,
                     'zero_lambda': lam0, 'zero_lambda_lo': lam0_lo, 'zero_lambda_hi': lam0_hi}
        if args.export_hist:
            hists.append(hist_rows(strain, smp, 'snp_starts', x))
    df_snp = pd.DataFrame.from_dict(rows, orient='index').add_prefix('snp_starts_').reset_index(names=['sample'])
    if len(df_snp):
        df_snp['snp_starts_zero_vs_robust'] = df_snp['snp_starts_zero_lambda'] / df_snp['snp_starts_lambda']
    df_cov = df_cov.merge(df_snp, on='sample', how='left')


    # ── 'best': one block per read (most completely covered, ties -> closest to the 5' end) ──
    # exposure from the block layout and the read lengths -> Poisson with ~(n_kmers + L - k) start positions
    # per isolated block, i.e. far fewer empty blocks at low depth than read starts
    len_path = args.read_lengths or args.block_hits.replace('.block_hits.tsv.gz', '.read_lengths.tsv.gz')
    best_scale = None
    if 'best' in df_bl.columns and os.path.exists(len_path):
        df_len = pd.read_csv(len_path, sep='\t')
        df_len['sample'] = df_len['sample'].astype(str)
        reps, wts, mean_len = length_bins(df_len, k)
        clen_b = contig_len if contig_len is not None else (df_blocks.groupby('contig_id')['stop'].max() + 1).to_dict()
        E_len = exposure_by_length(df_blocks, clen_b, reps, k)
        nk = df_blocks.set_index('block_id').loc[E_len.index, 'n_kmers'].to_numpy()
        iso = np.all([np.isclose(E_len[L].to_numpy(), nk + (L - k)) for L in E_len.columns], axis=0)
        iso_ids = set(E_len.index[iso & (nk == snp_n)])
        if len(iso_ids) < args.min_iso_blocks:      # e.g. near-whole-genome rare sets: blocks are back to back
            print(f"NOTE: only {len(iso_ids)} isolated SNP blocks (< --min_iso_blocks {args.min_iso_blocks}); "
                  f"best_iso_* left empty")
            iso_ids = set()
        E_s = pd.DataFrame(sample_exposure(E_len, wts), index=E_len.index, columns=wts.index)
        if args.export_blocks:
            blong = df_bl[['block_id', 'sample', 'frags', 'starts', 'best', 'contig_id', 'start', 'stop', 'n_kmers',
                           'ref_copies', 'family_id', 'segment_id']].copy()
            blong['exposure'] = [E_s.at[b, s_] if s_ in E_s.columns else np.nan
                                 for b, s_ in zip(blong['block_id'], blong['sample'])]
            blong['isolated'] = blong['block_id'].isin(iso_ids).astype(int)
            bpath = os.path.join(args.output_dir, f'{strain}.blocks_long.tsv.gz')
            blong.sort_values(['sample', 'contig_id', 'start']).to_csv(bpath, sep='\t', index=False)
            print(f'Saved per-block table to\n {bpath}')
        rows = {}
        for smp, d in df_bl.loc[single].groupby('sample'):
            if smp not in E_s.columns:
                continue
            y = d['best'].to_numpy()
            e = E_s[smp].reindex(d['block_id']).to_numpy()
            r_naive = y.sum() / e.sum()
            r_zip, pi_zip = zip_rate(y, e, min_pos=zi_params['min_pos'])
            r_zero = rate_from_zeros(y, e)
            disp = (((y - r_naive * e) ** 2 / (r_naive * e)).sum() / max(len(y) - 1, 1)) if r_naive > 0 else np.nan
            di = d.loc[d['block_id'].isin(iso_ids)]
            x = di['best'].to_numpy()
            e_iso = float(E_s.loc[list(di['block_id']), smp].mean()) if len(di) else np.nan
            rr = robust_0_inflation_safe(x, label=f'{strain} best_iso {smp}', **zi_params) if len(x) else \
                robust_0_inflation(np.zeros(0, np.int64), **zi_params)
            lam0, lam0_lo, lam0_hi = lambda_from_zeros(x) if len(x) else (np.nan,) * 3
            rows[smp] = {'total': int(y.sum()), 'n': len(y), 'frac_hit': (y > 0).mean(),
                         'exposure_mean': e.mean(), 'rate': r_naive, 'zip_rate': r_zip, 'zip_pi': pi_zip,
                         'zero_rate': r_zero, 'dispersion': disp,
                         'iso_n': len(x), 'iso_exposure': e_iso, 'iso_obs_zero': (x == 0).mean() if len(x) else np.nan,
                         'iso_lambda': rr['lambda'], 'iso_zi_frac': rr['zi_frac'], 'iso_zi_index': rr['zi_index'],
                         'iso_frac_trimmed': rr['frac_trimmed'],
                         'iso_zero_lambda': lam0, 'iso_zero_lambda_lo': lam0_lo, 'iso_zero_lambda_hi': lam0_hi,
                         'rate_upper95': -np.log(0.05) / e.sum(),
                         'mean_read_len': mean_len.get(smp, np.nan)}
            # diagnostics: exposure model, zero excess, clustering of empty blocks, size dependence
            order = np.lexsort((d['start'].to_numpy(), d['contig_id'].astype(str).to_numpy()))
            p0 = np.exp(-r_naive * e[order]) if r_naive > 0 else np.ones(len(e))
            run_max, run_p = zero_run_test(y[order] == 0, p0)
            r_starts = d['starts'].sum() / d['n_kmers'].sum()
            snp_m = (d['n_kmers'] == snp_n).to_numpy()
            r_snp = y[snp_m].sum() / e[snp_m].sum() if snp_m.any() and e[snp_m].sum() > 0 else np.nan
            r_big = y[~snp_m].sum() / e[~snp_m].sum() if (~snp_m).any() and e[~snp_m].sum() > 0 else np.nan
            rows[smp].update({
                'obs_zero_frac': (y == 0).mean(),
                'exp_zero_frac': np.exp(-r_naive * e).mean() if r_naive > 0 else np.nan,
                'zero_run_max': run_max, 'zero_run_p': run_p,
                'diag_vs_starts_rate': r_naive / r_starts if r_starts > 0 else np.nan,
                'diag_iso_vs_all_rate': (rr['lambda'] / e_iso) / r_naive if r_naive > 0 and e_iso else np.nan,
                'diag_bigblock_vs_snp_rate': r_big / r_snp if r_snp and r_snp > 0 else np.nan,
            })
            if args.export_hist and len(x):
                hists.append(hist_rows(strain, smp, 'best_iso', x))
        df_best = pd.DataFrame.from_dict(rows, orient='index').add_prefix('best_').reset_index(names=['sample'])
        if len(df_best):
            # no empty block at all -> the zero-based lambda is unbounded; leave the ratio undefined
            df_best['best_iso_zero_vs_robust'] = (df_best['best_iso_zero_lambda'] / df_best['best_iso_lambda']
                                                  ).replace([np.inf, -np.inf], np.nan)
            df_cov = df_cov.merge(df_best, on='sample', how='left')

        # ── repeat families: every read on any copy counts once to the family, exposure = sum over its copies ──
        # placed reads: their best block is a multi-copy block; reads inside a repeat that fit several copies
        # equally: best_ambiguous from family_hits. Copy number is taken from the reference, so a strain whose
        # mobile elements gained or lost copies shows up as rep_rate_ratio != 1 (feature, not an error).
        fam_path = args.family_hits or args.block_hits.replace('.block_hits.tsv.gz', '.family_hits.tsv.gz')
        multi = df_bl['ref_copies'] > 1
        if multi.any():
            df_fam = (pd.read_csv(fam_path, sep='\t') if os.path.exists(fam_path)
                      else pd.DataFrame(columns=['family_id', 'sample', 'best_ambiguous']))
            df_fam['sample'] = df_fam['sample'].astype(str)
            if 'best_ambiguous' not in df_fam.columns:
                df_fam['best_ambiguous'] = 0
            rep_rows, all_rows = {}, {}
            for smp, d in df_bl.loc[multi].groupby('sample'):
                if smp not in E_s.columns:
                    continue
                d = d.assign(E=E_s[smp].reindex(d['block_id']).to_numpy())
                fam = d.groupby('family_id').agg(placed=('best', 'sum'), E=('E', 'sum'),
                                                 copies=('ref_copies', 'max'))
                amb = df_fam.loc[df_fam['sample'] == smp].set_index('family_id')['best_ambiguous']
                fam['amb'] = amb.reindex(fam.index).fillna(0).to_numpy()
                fam['count'] = fam['placed'] + fam['amb']
                y, e = fam['count'].to_numpy(), fam['E'].to_numpy()
                r_rep = y.sum() / e.sum() if e.sum() > 0 else np.nan
                r_zip_rep, pi_rep = zip_rate(y, e, min_pos=zi_params['min_pos'])
                rep_rows[smp] = {'n_families': len(fam), 'copies_mean': fam['copies'].mean(),
                                 'count': int(y.sum()), 'ambiguous': int(fam['amb'].sum()),
                                 'exposure': e.sum(), 'frac_hit': (y > 0).mean(),
                                 'rate': r_rep, 'zip_rate': r_zip_rep, 'zip_pi': pi_rep}
                # combined: single-copy blocks + repeat families as one set of Poisson units
                ds_ = df_bl.loc[single & (df_bl['sample'] == smp)]
                y_all = np.concatenate([ds_['best'].to_numpy(), y])
                e_all = np.concatenate([E_s[smp].reindex(ds_['block_id']).to_numpy(), e])
                r_zip_all, pi_all = zip_rate(y_all, e_all, min_pos=zi_params['min_pos'])
                all_rows[smp] = {'n_units': len(y_all), 'frac_hit': (y_all > 0).mean(),
                                 'rate': y_all.sum() / e_all.sum(), 'zip_rate': r_zip_all, 'zip_pi': pi_all,
                                 'zero_rate': rate_from_zeros(y_all, e_all),
                                 'rate_upper95': -np.log(0.05) / e_all.sum()}
            if rep_rows:
                df_rep = pd.DataFrame.from_dict(rep_rows, orient='index').add_prefix('rep_').reset_index(names=['sample'])
                df_all = pd.DataFrame.from_dict(all_rows, orient='index').add_prefix('allunits_').reset_index(names=['sample'])
                df_cov = df_cov.merge(df_rep, on='sample', how='left').merge(df_all, on='sample', how='left')
                df_cov['rep_rate_ratio'] = df_cov['rep_rate'] / df_cov['best_rate'].replace(0, np.nan)
        best_wide = df_bl.loc[df_bl['block_id'].isin(iso_ids)].pivot(index='block_id', columns='sample', values='best')
    else:
        best_wide = None
        print(f"NOTE: no 'best' column or no read-length file ({len_path}); skipping the best-block features")

    # ── k-mer subsets: classic coverage / depth / robust zero inflation on validated hits ──
    kmer_subsets = {'all': df_kmer_hits}
    df_loc = None
    if args.kmer_loci:
        df_loc = pd.read_csv(args.kmer_loci, sep='\t')
        df_loc = df_loc.loc[df_loc['ref_copies'] == 1].merge(df_blocks[['block_id', 'n_kmers']], on='block_id', how='left')
        snp_kmers = df_loc.loc[df_loc['n_kmers'] == snp_n, '#kmer']
        long_kmers = df_loc.loc[df_loc['n_kmers'] > args.long_kmers, '#kmer']
        kmer_subsets['snp'] = df_kmer_hits.loc[df_kmer_hits['#kmer'].isin(snp_kmers)]
        kmer_subsets['long'] = df_kmer_hits.loc[df_kmer_hits['#kmer'].isin(long_kmers)]
    for sub, df_sub in kmer_subsets.items():
        if len(df_sub) == 0:
            continue
        df_cov = df_cov.merge(kmer_subset_features(df_sub, sub, zi_params, strain), on='sample', how='left')
        if args.export_hist:
            for smp in samples:
                hists.append(hist_rows(strain, smp, f'kmer_{sub}', df_sub[smp].to_numpy()))

    if df_loc is not None:
        df_cov = df_cov.merge(window_evenness(df_loc, df_kmer_hits), on='sample', how='left')

    # ── depth (fold coverage) and relative abundance (% of sequenced bases) for every estimate ──
    # an estimator whose column is missing (e.g. no SNP-size blocks -> no snp_starts_*) gives None and is skipped
    L = df_cov['read_len']
    starts_scale = L / snp_n                   # starts per SNP block -> start rate per position -> read depth
    # fragments per k-mer -> read depth: a read covers L - k + 1 k-mers, and kmer_hits counts FRAGMENTS
    # (both mates of a pair count once), so scale by reads per fragment (2 for paired, 1 for single-end)
    kmer_scale = L / (L - k + 1) * df_cov['reads'] / df_cov['fragments']
    depth = {
        'starts_robust': _scaled(df_cov, 'snp_starts_lambda', starts_scale),
        'starts_ztp': _scaled(df_cov, 'snp_starts_ztp_lambda', starts_scale),
        'starts_ztnb': _scaled(df_cov, 'snp_starts_ztnb_mu', starts_scale),
        'starts_naive': _scaled(df_cov, 'snp_starts_mean', starts_scale),
        'starts_zero': _scaled(df_cov, 'snp_starts_zero_lambda', starts_scale),
        'starts_zero_lo': _scaled(df_cov, 'snp_starts_zero_lambda_lo', starts_scale),
        'starts_zero_hi': _scaled(df_cov, 'snp_starts_zero_lambda_hi', starts_scale),
        'starts_zero_allblocks': df_cov['rate_zero_allblocks'] * L,
    }
    if 'best_rate' in df_cov:
        Lb = df_cov['best_mean_read_len']             # rate (reads per start position) x mean read length = depth
        depth.update({
            'best_naive': df_cov['best_rate'] * Lb,
            'best_zip': df_cov['best_zip_rate'] * Lb,
            'best_zero': df_cov['best_zero_rate'] * Lb,
            'best_robust': df_cov['best_iso_lambda'] / df_cov['best_iso_exposure'] * Lb,
            'best_iso_zero': df_cov['best_iso_zero_lambda'] / df_cov['best_iso_exposure'] * Lb,
            'best_iso_zero_lo': df_cov['best_iso_zero_lambda_lo'] / df_cov['best_iso_exposure'] * Lb,
            'best_iso_zero_hi': df_cov['best_iso_zero_lambda_hi'] / df_cov['best_iso_exposure'] * Lb,
        })
        if 'rep_rate' in df_cov:
            depth.update({
                'rep': df_cov['rep_rate'] * Lb,
                'rep_zip': df_cov['rep_zip_rate'] * Lb,
                'allunits': df_cov['allunits_rate'] * Lb,
                'allunits_zip': df_cov['allunits_zip_rate'] * Lb,
                'allunits_zero': df_cov['allunits_zero_rate'] * Lb,
            })
    for sub in kmer_subsets:
        if f'kmer_{sub}_lambda' in df_cov:
            depth[f'kmer_{sub}_robust'] = df_cov[f'kmer_{sub}_lambda'] * kmer_scale
            depth[f'kmer_{sub}_naive'] = df_cov[f'kmer_{sub}_depth_mean'] * kmer_scale
            depth[f'kmer_{sub}_zero'] = df_cov[f'kmer_{sub}_zero_lambda'] * kmer_scale
    for est, dep in depth.items():
        if dep is None:
            continue
        df_cov[f'depth_{est}'] = dep
        df_cov[f'relative_abundance_{est}'] = dep * genome_len / df_cov['total_bases'] * 100

    # zero-inflated depth is the depth of the PRESENT part of the genome; scale by (1 - pi) for the share of
    # sequenced bases (a strain missing 20% of the reference blocks contributes 20% fewer bases)
    for est, pi_col in (('best_zip', 'best_zip_pi'), ('allunits_zip', 'allunits_zip_pi')):
        if f'relative_abundance_{est}' in df_cov and pi_col in df_cov:
            df_cov[f'relative_abundance_{est}_present'] = df_cov[f'relative_abundance_{est}'] * (1 - df_cov[pi_col].fillna(0))

    # detection limit: 95% upper bound if nothing were seen (reported for every sample)
    if 'snp_starts_n' in df_cov:
        df_cov['detection_limit_rel_abund_starts'] = (-np.log(0.05) / df_cov['snp_starts_n']) * starts_scale \
            * genome_len / df_cov['total_bases'] * 100
    if 'allunits_rate_upper95' in df_cov:          # repeats included: lowest limit
        df_cov['detection_limit_rel_abund'] = df_cov['allunits_rate_upper95'] * df_cov['best_mean_read_len'] \
            * genome_len / df_cov['total_bases'] * 100
    elif 'best_rate_upper95' in df_cov:
        df_cov['detection_limit_rel_abund'] = df_cov['best_rate_upper95'] * df_cov['best_mean_read_len'] \
            * genome_len / df_cov['total_bases'] * 100
    elif 'detection_limit_rel_abund_starts' in df_cov:
        df_cov['detection_limit_rel_abund'] = df_cov['detection_limit_rel_abund_starts']

    # ── diagnostics (debugging after a run) ──
    if 'best_mean_read_len' in df_cov and len_path and os.path.exists(len_path):
        n_pieces = df_len.groupby('sample')['n'].sum()
        df_cov['diag_frac_reads_scanned'] = df_cov['sample'].map(n_pieces) / df_cov['reads']   # pieces >= k per read
    df_cov['diag_frac_fragments_hit'] = df_cov['frags_with_hits'] / df_cov['fragments']
    if 'depth_best_zip' in df_cov:
        for other in ('kmer_all_robust', 'starts_robust', 'allunits_zip'):
            if f'depth_{other}' in df_cov:
                df_cov[f'diag_depth_{other}_vs_best_zip'] = df_cov[f'depth_{other}'] / df_cov['depth_best_zip'].replace(0, np.nan)
    if 'best_obs_zero_frac' in df_cov:
        df_cov['diag_best_zero_excess'] = df_cov['best_obs_zero_frac'] - df_cov['best_exp_zero_frac']

    # filter diagnostics of this strain (kmer_scrub_filter.py selection.json) -> strain_filter_* columns
    sel_path = args.selection or args.blocks.replace('.blocks.tsv.gz', '.selection.json')
    if os.path.exists(sel_path):
        with open(sel_path) as fh:
            sel = json.load(fh)
        for key, v in sel.items():
            if isinstance(v, (int, float, str, bool)) or v is None:
                name = key if key.startswith('filter_') else f'filter_{key}'
                df_cov[f'strain_{name}'] = v

    # strain-level descriptors (constant per strain, for stratifying)
    df_cov['strain'] = strain
    df_cov['genome_length'] = genome_len
    df_cov['strain_n_kmers'] = len(df_kmer_hits)
    df_cov['strain_n_blocks'] = len(df_blocks)
    df_cov['strain_n_snp_blocks'] = int((df_blocks['n_kmers'].eq(snp_n) & df_blocks['ref_copies'].eq(1)).sum())
    rep_b = df_blocks.loc[df_blocks['ref_copies'] > 1]
    df_cov['strain_n_repeat_families'] = int(rep_b['family_id'].nunique())
    df_cov['strain_repeat_kmer_frac'] = round(rep_b['n_kmers'].sum() / df_blocks['n_kmers'].sum(), 4)
    df_cov['strain_rare_count'] = int((df_blocks['pangenome_count'] + df_blocks['metagenome_count']).max())
    df_cov['strain_gc'] = round(np.mean([(s.count('G') + s.count('C')) / len(s) for s in df_kmer_hits['#kmer']]), 4)

    num = df_cov.select_dtypes(include=[np.number]).columns
    df_cov[num] = df_cov[num].replace([np.inf, -np.inf], np.nan)
    out = df_cov.set_index(['strain', 'sample']).round(6)
    out_path = os.path.join(args.output_dir, 'block_coverage_depth.tsv')
    out.to_csv(out_path, sep='\t')
    print(f'Saved {out.shape[0]} strain x sample rows, {out.shape[1]} features to\n {out_path}')

    if args.export_hist and hists:
        hist_path = os.path.join(args.output_dir, f'{strain}.block_hist.tsv.gz')
        pd.concat(hists, ignore_index=True).to_csv(hist_path, sep='\t', index=False)
        print(f'Saved count histograms to\n {hist_path}')

    if args.figures != 'none':
        zi_subsets = {}
        snp_wide = (df_bl.loc[df_bl['is_snp']].pivot(index='block_id', columns='sample', values='starts'))
        if best_wide is not None and len(best_wide):
            zi_subsets['best_iso'] = best_wide
        if len(snp_wide):
            zi_subsets['snp_starts'] = snp_wide
        for sub, df_sub in kmer_subsets.items():
            if len(df_sub):
                zi_subsets[f'kmer_{sub}'] = df_sub.set_index('#kmer')
        plot_zi_distribution(zi_subsets, args.output_dir, strain, samples, args.figures, zi_params)
        if contig_len is None:
            contig_len = (df_blocks.groupby('contig_id')['stop'].max() + 1).to_dict()
        plot_block_map(df_blocks, df_bl, contig_len, args.output_dir, strain, samples, args.figures)

    show = ['blk_frac_hit', 'best_total', 'best_frac_hit', 'best_dispersion', 'best_zip_pi', 'best_iso_zi_frac',
            'best_iso_zero_vs_robust', 'rep_count', 'rep_rate_ratio', 'allunits_zip_pi', 'kmer_all_coverage', 'reject_frac',
            'relative_abundance_best_robust', 'relative_abundance_best_zip', 'relative_abundance_allunits_zip',
            'relative_abundance_starts_robust', 'detection_limit_rel_abund']
    print(out[[c for c in show if c in out.columns]].round(4).T)


if __name__ == '__main__':
    main()