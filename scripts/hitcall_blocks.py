import gzip
import numpy as np
import pandas as pd
import plotly.express as px

ref_path = '/metrica/projects/faith/strainer/20261110-readhash/dev/D6331_Faecalibacterium_prausnitzii.fasta.gz'
kmer_hits_path = '/metrica/projects/faith/strainer/20261110-readhash/dev_read_locations/strain_detect/D6331_Faecalibacterium_prausnitzii.kmer_hits.tsv.gz'
block_hits_path = '/metrica/projects/faith/strainer/20261110-readhash/dev_read_locations/strain_detect/D6331_Faecalibacterium_prausnitzii.block_hits.tsv.gz'
blocks_path = '/metrica/projects/faith/strainer/20261110-readhash/dev_read_locations/filter/D6331_Faecalibacterium_prausnitzii.blocks.tsv.gz'
detect_stats_path = '/metrica/projects/faith/strainer/20261110-readhash/dev_read_locations/strain_detect/D6331_Faecalibacterium_prausnitzii.detect_stats.tsv'
outpath = '/metrica/projects/faith/strainer/20261110-readhash/dev/output-test'

basename = 'test'
k = 31


def ztp_lambda(mean_pos, n_iter=200):
    # zero-truncated Poisson: solve mean_pos = lam / (1 - exp(-lam)) for lam
    # uses only the non-zero counts, so it is not pulled down by excess zeros
    if not np.isfinite(mean_pos) or mean_pos <= 1:
        return 0.0
    lam = mean_pos
    for _ in range(n_iter):
        lam = mean_pos * (1 - np.exp(-lam))
    return lam


def zero_inflation(counts):
    # robust zero inflation: lambda from the positives, pi from the zero excess
    counts = np.asarray(counts)
    n = len(counts)
    pos = counts[counts > 0]
    obs_zero = (counts == 0).mean()
    lam = ztp_lambda(pos.mean()) if len(pos) else 0.0
    exp_zero = np.exp(-lam)
    pi = (obs_zero - exp_zero) / (1 - exp_zero) if lam > 0 else np.nan
    return {'n': n, 'lambda': lam, 'obs_zero': obs_zero, 'exp_zero': exp_zero,
            'pi': np.clip(pi, 0, 1) if lam > 0 else np.nan}


def genome_length(path):
    opener = gzip.open if path.endswith('.gz') else open
    n = 0
    with opener(path, 'rt') as fh:
        for line in fh:
            if not line.startswith('>'):
                n += len(line.strip())
    return n


def lambda_from_zeros(counts, z=1.96):
    # inverse Poisson: P(0) = exp(-lambda)  ->  lambda = -ln(fraction of zeros)
    # assumes NO zero inflation; returns lambda and a 95% CI from the binomial error on the zero fraction
    counts = np.asarray(counts)
    n = len(counts)
    p0 = (counts == 0).mean()
    if p0 == 1:
        return 0.0, 0.0, -np.log(0.05) / n   # 95% upper bound when nothing was seen: exp(-lam * n) = 0.05
    if p0 == 0:
        return np.inf, np.nan, np.inf
    se = np.sqrt(p0 * (1 - p0) / n)
    lo_p, hi_p = min(p0 + z * se, 1.0), max(p0 - z * se, 1e-12)
    return -np.log(p0), -np.log(lo_p), -np.log(hi_p)


def rate_from_zeros(starts, n_kmers, n_iter=100):
    # same with unequal exposure: solve sum_b exp(-r * n_b) = observed zero blocks (bisection)
    starts = np.asarray(starts); n_kmers = np.asarray(n_kmers, float)
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
        if np.exp(-mid * n_kmers).sum() > obs:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


df_kmer_hits = pd.read_csv(kmer_hits_path, sep='\t')
df_block_hits = pd.read_csv(block_hits_path, sep='\t')
df_blocks = pd.read_csv(blocks_path, sep='\t')
df_stats = pd.read_csv(detect_stats_path, sep='\t')
genome_len = genome_length(ref_path)

# total evaluated row
total_evaluated = df_kmer_hits.loc[df_kmer_hits['#kmer'] == 'total_evaluated'].drop(columns='#kmer').T.reset_index()
total_evaluated.columns = ['sample', 'total_evaluated']
df_kmer_hits = df_kmer_hits.loc[df_kmer_hits['#kmer'] != 'total_evaluated'].copy()

# long formats
df_kmer_long = df_kmer_hits.melt(id_vars='#kmer', var_name='sample', value_name='hit_count')
df_block_long = df_block_hits.merge(df_blocks, on='block_id', how='left')

# ── 1. starts: Poisson on single-SNP blocks (identical exposure: 31 k-mers, single copy) ──
# the start rate drops with block size (a read must be error-free across its whole overlap
# with a segment), so only same-size blocks share one Poisson rate
df_snp = df_block_long.loc[(df_block_long['n_kmers'] == k) & (df_block_long['ref_copies'] == 1)]

df_starts = (
    df_snp
    .groupby('sample')['starts']
    .apply(lambda x: pd.Series(zero_inflation(x)))
    .unstack()
    .add_prefix('snp_starts_')
    .reset_index()
)
df_starts['snp_start_rate'] = df_starts['snp_starts_lambda'] / k          # starts per k-mer position
df_starts['snp_starts_exp_zero_blocks'] = df_starts['snp_starts_exp_zero'] * df_starts['snp_starts_n']
df_starts['snp_starts_obs_zero_blocks'] = df_starts['snp_starts_obs_zero'] * df_starts['snp_starts_n']

# ── 2. block level summary over all blocks ──
df_block_summary = (
    df_block_long
    .groupby('sample')
    .agg(**{'blocks': ('block_id', 'size'),
            'blocks_hit': ('frags', lambda x: (x > 0).sum()),
            'frac_blocks_hit': ('frags', lambda x: (x > 0).mean()),
            'block_frags_mean': ('frags', 'mean'),
            'starts_total': ('starts', 'sum'),
            'kmer_positions': ('n_kmers', 'sum')})
    .reset_index()
)
# QC: start rate in long blocks vs SNP blocks (low = many reads rejected, e.g. untrimmed adapters/errors)
df_long_blocks = df_block_long.loc[(df_block_long['n_kmers'] > 120) & (df_block_long['ref_copies'] == 1)]
df_qc = (
    df_long_blocks.groupby('sample')
    .agg(**{'long_starts': ('starts', 'sum'), 'long_kmers': ('n_kmers', 'sum')})
    .reset_index()
)
df_qc['long_block_start_rate'] = df_qc['long_starts'] / df_qc['long_kmers']

# ── 3. classic k-mer level: coverage of rare k-mers + robust zero inflation ──
df_kmer_summary = (
    df_kmer_long
    .groupby('sample')['hit_count']
    .apply(lambda x: pd.Series(zero_inflation(x)))
    .unstack()
    .add_prefix('kmer_')
    .reset_index()
)
df_kmer_summary['coverage_rare_kmers'] = 1 - df_kmer_summary['kmer_obs_zero']
df_kmer_summary['expected_coverage_rare_kmers'] = 1 - df_kmer_summary['kmer_exp_zero']
df_kmer_summary['kmer_depth_mean'] = df_kmer_long.groupby('sample')['hit_count'].mean().values

# ── combine ──
df_hitcall = (
    df_block_summary
    .merge(df_starts, on='sample', how='left')
    .merge(df_qc[['sample', 'long_block_start_rate']], on='sample', how='left')
    .merge(df_kmer_summary, on='sample', how='left')
    .merge(total_evaluated, on='sample', how='left')
)
df_hitcall['rate_ratio_long_vs_snp'] = df_hitcall['long_block_start_rate'] / df_hitcall['snp_start_rate']

# ── 4. depth and relative abundance from the different means ──
# sample size: reads and mean read length (from the k-mer positions scanned per read)
df_stats['read_len'] = df_stats['kmer_positions'] / df_stats['reads'] + k - 1
df_stats['total_bases'] = df_stats['kmer_positions'] + df_stats['reads'] * (k - 1)
df_hitcall = df_hitcall.merge(df_stats[['sample', 'reads', 'read_len', 'total_bases']], on='sample', how='left')

# naive means (zeros included) next to the robust zero-truncated lambdas
df_naive = (
    df_snp.groupby('sample')
    .agg(**{'snp_starts_mean': ('starts', 'mean')})
    .reset_index()
)
df_hitcall = df_hitcall.merge(df_naive, on='sample', how='left')

L = df_hitcall['read_len']
# read starts: start rate per position x read length = read depth (fold coverage)
df_hitcall['depth_starts_robust'] = df_hitcall['snp_starts_lambda'] / k * L
df_hitcall['depth_starts_naive'] = df_hitcall['snp_starts_mean'] / k * L
# k-mer hits: fragments per k-mer -> read depth (a read of length L covers L - k + 1 k-mers;
# kmer_hits counts fragments, so x reads per fragment: 2 for paired, 1 for single-end)
df_hitcall = df_hitcall.merge(df_stats[['sample', 'fragments']], on='sample', how='left')
kmer_scale = L / (L - k + 1) * df_hitcall['reads'] / df_hitcall['fragments']
df_hitcall['depth_kmer_robust'] = df_hitcall['kmer_lambda'] * kmer_scale
df_hitcall['depth_kmer_naive'] = df_hitcall['kmer_depth_mean'] * kmer_scale

# inverse: lambda from the fraction of ZERO blocks / k-mers (assumes no zero inflation)
df_zero_snp = (
    df_snp.groupby('sample')['starts']
    .apply(lambda x: pd.Series(lambda_from_zeros(x), index=['lam', 'lam_lo', 'lam_hi']))
    .unstack()
    .add_prefix('snp_starts_zero_')
    .reset_index()
)
df_zero_kmer = (
    df_kmer_long.groupby('sample')['hit_count']
    .apply(lambda x: pd.Series(lambda_from_zeros(x), index=['lam', 'lam_lo', 'lam_hi']))
    .unstack()
    .add_prefix('kmer_zero_')
    .reset_index()
)
# all single-copy blocks, each with its own exposure n_kmers (mixes block sizes -> biased low, see rate table)
df_zero_all = (
    df_block_long.loc[df_block_long['ref_copies'] == 1]
    .groupby('sample')
    .apply(lambda d: rate_from_zeros(d['starts'], d['n_kmers']), include_groups=False)
    .rename('all_blocks_zero_rate')
    .reset_index()
)
df_hitcall = (
    df_hitcall
    .merge(df_zero_snp, on='sample', how='left')
    .merge(df_zero_kmer, on='sample', how='left')
    .merge(df_zero_all, on='sample', how='left')
)

df_hitcall['depth_starts_zero'] = df_hitcall['snp_starts_zero_lam'] / k * L
df_hitcall['depth_starts_zero_lo'] = df_hitcall['snp_starts_zero_lam_lo'] / k * L
df_hitcall['depth_starts_zero_hi'] = df_hitcall['snp_starts_zero_lam_hi'] / k * L
df_hitcall['depth_starts_zero_allblocks'] = df_hitcall['all_blocks_zero_rate'] * L
df_hitcall['depth_kmer_zero'] = df_hitcall['kmer_zero_lam'] * kmer_scale

# zero-based vs positive-based lambda: ~1 = plain Poisson, < 1 = excess zeros (technical or missing blocks)
df_hitcall['snp_zero_vs_robust'] = df_hitcall['snp_starts_zero_lam'] / df_hitcall['snp_starts_lambda']

# relative (sequence) abundance: bases from this genome / all sequenced bases
for m in ['starts_robust', 'starts_naive', 'starts_zero', 'starts_zero_lo', 'starts_zero_hi',
          'starts_zero_allblocks', 'kmer_robust', 'kmer_naive', 'kmer_zero']:
    df_hitcall[f'rel_abund_{m}'] = df_hitcall[f'depth_{m}'] * genome_len / df_hitcall['total_bases']

cols = ['sample', 'frac_blocks_hit', 'block_frags_mean', 'starts_total',
        'snp_starts_n', 'snp_start_rate', 'snp_starts_obs_zero_blocks', 'snp_starts_exp_zero_blocks', 'snp_starts_pi',
        'coverage_rare_kmers', 'expected_coverage_rare_kmers', 'kmer_lambda', 'kmer_pi',
        'rate_ratio_long_vs_snp']
display(df_hitcall[cols].round(4))

df_abund = (
    df_hitcall
    .melt(id_vars='sample',
          value_vars=[c for c in df_hitcall.columns if c.startswith('rel_abund_')],
          var_name='estimate', value_name='rel_abund')
    .assign(estimate=lambda x: x['estimate'].str.replace('rel_abund_', ''))
)
df_abund['depth'] = df_abund.apply(lambda r: df_hitcall.set_index('sample').loc[r['sample'], f"depth_{r['estimate']}"], axis=1)
display(df_abund.pivot(index='estimate', columns='sample', values=['depth', 'rel_abund']).round(4))
display(df_hitcall[['sample', 'snp_starts_lambda', 'snp_starts_zero_lam', 'snp_zero_vs_robust']].round(4))

# ── plot: observed vs expected starts per SNP block (Poisson with the robust lambda) ──
df_plot = (
    df_snp.groupby(['sample', 'starts']).size().rename('observed').reset_index()
    .merge(df_starts[['sample', 'snp_starts_lambda', 'snp_starts_n', 'snp_starts_pi']], on='sample')
)
from scipy.stats import poisson
df_plot['expected'] = (1 - df_plot['snp_starts_pi'].fillna(0)) * df_plot['snp_starts_n'] * \
    poisson.pmf(df_plot['starts'], df_plot['snp_starts_lambda'])
df_plot.loc[df_plot['starts'] == 0, 'expected'] += df_plot['snp_starts_pi'].fillna(0) * df_plot['snp_starts_n']

fig = px.bar(df_plot.melt(id_vars=['sample', 'starts'], value_vars=['observed', 'expected'],
                          var_name='type', value_name='blocks'),
             x='starts', y='blocks', color='type', barmode='group', facet_col='sample',
             title='Read starts per SNP block (31 k-mers): observed vs zero-inflated Poisson')
fig.show()
