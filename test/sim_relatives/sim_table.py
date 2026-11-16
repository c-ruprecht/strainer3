import sys, numpy as np, pandas as pd
O, S = sys.argv[1], sys.argv[2]
modes = ['nonaware', 'block', 'veto', 'veto_ctx', 'pairs_veto', 'pairs_veto_ctx']
res = {}
for m in modes:
    h = pd.read_csv(f'{O}/det_{m}/{S}.kmer_hits.tsv.gz', sep='\t')
    h = h[h['#kmer'] != 'total_evaluated'].set_index('#kmer')
    res[m] = pd.DataFrame({'cov': (h > 0).mean(), 'hits': h.sum()})
cov = pd.DataFrame({m: res[m]['cov'] for m in modes}); hits = pd.DataFrame({m: res[m]['hits'] for m in modes})
pd.set_option('display.width', 250)
def lab(i):
    p = i.split('_'); 
    return ('true', 0.0, p[1], p[2]) if p[0] == 'true' else (p[1], float(p[2]), p[3], p[4])
idx = pd.MultiIndex.from_tuples([lab(i) for i in cov.index], names=['kind', 'div', 'depth', 'lib'])
cov.index = idx; hits.index = idx
print(f'== {S}: k-mer coverage (fraction of rare k-mers hit) ==')
print(cov.sort_index().round(4).to_string())
# relative to the true strain at the same depth and library, same method
print('\n== relative coverage / true-strain coverage at the same depth (3x) and library ==')
rows = []
for (kind, div, depth, lib), r in cov.iterrows():
    if kind == 'true': continue
    t = cov.loc[('true', 0.0, depth, lib)]
    rows.append(dict(model=kind, ani=100 * (1 - div), lib=lib, **(r / t).round(3).to_dict()))
out = pd.DataFrame(rows).sort_values(['model', 'lib', 'ani'], ascending=[True, True, False])
print(out.to_string(index=False))
print('\n== true strain: hits kept relative to non-aware ==')
t = hits.xs('true', level='kind'); print((t.div(t['nonaware'], axis=0)).round(3).to_string())
out.to_csv(f'{O}/sim_ratio.tsv', sep='\t', index=False)
# reads-level contradiction rate
for m in ('block', 'veto_ctx', 'pairs_veto_ctx'):
    st = pd.read_csv(f'{O}/det_{m}/{S}.detect_stats.tsv', sep='\t').set_index('sample')
    st.index = idx
    c = (st.reads_contradicted / st.reads_hit).rename('contradicted_reads_frac')
    g = (st.groups_rejected / (st.groups_accepted + st.groups_rejected)).rename('reject_frac')
    pc = (st.pairs_discordant / st.pairs_placed.where(st.pairs_placed > 0)).rename('pairs_discordant_frac')
    print(f'\n== {m}: read-level contradiction, group reject_frac ==')
    print(pd.concat([c, g, st.pairs_placed, pc], axis=1).sort_index().round(3).to_string())
