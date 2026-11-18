"""Hit-dense query set for the race test: k-mers sampled from the test genomes
themselves (so every column gets many concurrent hits), forward and reverse
complement, plus random absent k-mers. Usage: make_dense_queries.py out.tsv.gz k genome.fa[.gz] ..."""
import gzip, random, sys
random.seed(11)
out, k, files = sys.argv[1], int(sys.argv[2]), sys.argv[3:]
comp = str.maketrans('ACGT', 'TGCA')
rows = []
for fn in files:
    op = gzip.open if fn.endswith('.gz') else open
    seq = ''.join(l.strip() for l in op(fn, 'rt') if not l.startswith('>')).upper()
    for _ in range(40000):
        p = random.randint(0, len(seq) - k); km = seq[p:p+k]
        if set(km) <= set('ACGT'):
            rows.append(km if random.random() < 0.7 else km.translate(comp)[::-1])
rows += [''.join(random.choice('ACGT') for _ in range(k)) for _ in range(20000)]
random.shuffle(rows)
with gzip.open(out, 'wt') as f:
    f.write('#kmer\n' + '\n'.join(rows) + '\n')
