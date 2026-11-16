import sys, numpy as np, pandas as pd
O, S = sys.argv[1], sys.argv[2]
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 200)
def zi(counts):
    """zero-inflation fraction of a count vector (zero-truncated Poisson mean -> lambda -> expected zeros)"""
    c = np.asarray(counts, float); pos = c[c > 0]
    if len(pos) < 20: return np.nan
    m = pos.mean(); lam = max(m - 1, 1e-3)
    for _ in range(100): lam = m * (1 - np.exp(-lam))
    if lam < 1e-6: return np.nan
    return max(0.0, 1 - (len(pos) / len(c)) / (1 - np.exp(-lam)))
na = pd.read_csv(f'{O}/det_nonaware/{S}.kmer_hits.tsv.gz', sep='\t'); na = na[na['#kmer'] != 'total_evaluated'].set_index('#kmer')
bk = pd.read_csv(f'{O}/det_block/{S}.kmer_hits.tsv.gz', sep='\t'); bk = bk[bk['#kmer'] != 'total_evaluated'].set_index('#kmer')
pr = pd.read_csv(f'{O}/det_pair/{S}.kmer_hits.tsv.gz', sep='\t'); pr = pr[pr['#kmer'] != 'total_evaluated'].set_index('#kmer')
bh = pd.read_csv(f'{O}/det_pair/{S}.block_hits.tsv.gz', sep='\t')
st = pd.read_csv(f'{O}/det_pair/{S}.detect_stats.tsv', sep='\t').set_index('sample')
rows = []
for smp in pr.columns:
    k, d, lib = smp.split('_'); b = bh[bh['sample'] == smp]
    sup = (b.frags > 0).sum(); ref = ((b.frags == 0) & (b.refused > 0)).sum()
    rows.append(dict(kind=k, lib=lib, depth=int(d[1:]), cov_nonaware=(na[smp] > 0).mean(), cov_block=(bk[smp] > 0).mean(), cov_pair=(pr[smp] > 0).mean(),
                     zi_nonaware=zi(na[smp]), zi_pair=zi(pr[smp]), blocks_refused=ref / max(sup + ref, 1),
                     frags_contra=st.loc[smp, 'frags_contradicted'] / max(st.loc[smp, 'frags_hit'], 1)))
d = pd.DataFrame(rows)
names = {'true': 'true strain', 'd999': '99.9 diffuse', 'd995': '99.5 diffuse', 't999': '99.9 tracts', 't995': '99.5 tracts', 't990': '99.0 tracts'}
d['kind'] = pd.Categorical(d.kind.map(names), list(names.values()))
d.to_csv(f'{O}/depth_summary.tsv', sep='\t', index=False)
for lib in ('long', 'short'):
    for v in ('cov_nonaware', 'cov_block', 'cov_pair', 'zi_nonaware', 'zi_pair', 'blocks_refused', 'frags_contra'):
        print(f'== {lib} insert: {v} by depth =='); print(d[d.lib == lib].pivot_table(index='kind', columns='depth', values=v, observed=True).round(3).to_string())
