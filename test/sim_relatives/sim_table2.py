import sys, numpy as np, pandas as pd
O, S = sys.argv[1], sys.argv[2]
modes = ['nonaware', 'block', 'pairs_veto_ctx', 'raw_indep', 'raw_pair_noctx', 'raw_pair', 'raw_pair_q20']
pd.set_option('display.width', 250)
def lab(i):
    p = i.split('_')
    return ('true', 100.0, p[1], p[2]) if p[0] == 'true' else (p[1], round(100 * (1 - float(p[2])), 2), p[3], p[4])
cov, hits, con = {}, {}, {}
for m in modes:
    h = pd.read_csv(f'{O}/det_{m}/{S}.kmer_hits.tsv.gz', sep='\t'); h = h[h['#kmer'] != 'total_evaluated'].set_index('#kmer')
    cov[m] = (h > 0).mean(); hits[m] = h.sum()
    try:
        st = pd.read_csv(f'{O}/det_{m}/{S}.detect_stats.tsv', sep='\t').set_index('sample')
        con[m] = st.reads_contradicted / st.reads_hit
    except Exception: pass
cov, hits, con = pd.DataFrame(cov), pd.DataFrame(hits), pd.DataFrame(con)
for d in (cov, hits, con): d.index = pd.MultiIndex.from_tuples([lab(i) for i in d.index], names=['kind', 'ani', 'depth', 'lib'])
rows = []
for (kind, ani, depth, lib), r in cov.iterrows():
    if kind == 'true': continue
    rows.append(dict(model=kind, lib=lib, ani=ani, **(r / cov.loc[('true', 100.0, depth, lib)]).round(3).to_dict()))
print(f'== {S}: relative coverage / true coverage, 3x ==')
print(pd.DataFrame(rows).sort_values(['model', 'lib', 'ani'], ascending=[True, True, False]).to_string(index=False))
t = hits.xs('true', level='kind')
print('\n== true strain: hits relative to non-aware =='); print(t.div(t['nonaware'], axis=0).round(3).to_string())
print('\n== contradicted reads fraction =='); print(con.sort_index(ascending=[True, False, True, True]).round(3).to_string())
