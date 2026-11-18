"""Adversarial reads/reference for the differential tests (k=31): IUPAC/U/lowercase/N
import random, gzip, os, sys
random.seed(3)
out = sys.argv[1]; os.makedirs(out, exist_ok=True)
K = 31
comp = str.maketrans('ACGT', 'TGCA')
rc = lambda s: s.translate(comp)[::-1]
def noisy(n, exotic=0.004):
    s = []
    for _ in range(n):
        r = random.random()
        if   r < 0.002: s.append('N' * random.randint(1, 20))
        elif r < 0.002 + exotic: s.append(random.choice('RYKMSWBDHVU'))
        elif r < 0.012: s.append(random.choice('acgt'))
        else: s.append(random.choice('ACGT'))
    return ''.join(s)
ref = [noisy(random.randint(8000, 30000)) for _ in range(3)] + ['ACGT' * 10 + 'GAATTC' * 20]
with gzip.open(f'{out}/ref.fna.gz', 'wt') as f:
    for i, c in enumerate(ref):
        f.write(f'>r{i}\n' + '\n'.join(c[j:j+70] for j in range(0, len(c), 70)) + '\n')
up = [c.upper() for c in ref]
# informative k-mers: reference windows (incl. IUPAC ones), revcomps, plus broken lines
inf = []
for _ in range(1500):
    c = random.choice(up); p = random.randint(0, len(c) - K); w = c[p:p+K]
    if 'N' in w: continue
    inf.append(w if random.random() < 0.6 or not set(w) <= set('ACGT') else rc(w))
inf += ['ACGT' * 5, 'A' * 40, ''.join(random.choice('ACGT') for _ in range(K))]   # wrong length / absent
inf += [random.choice(up)[100:100+K].lower()]
with open(f'{out}/inf.txt', 'w') as f:
    f.write('#kmer\n' + '\n'.join(inf) + '\n')
def read_from_ref():
    c = random.choice(up); L = random.choice([150, 150, 150, 100, 60, 31, 30, 20, 0 if random.random() < .1 else 150])
    if L == 0: return random.choice('ACGT') * 5
    p = random.randint(0, max(0, len(c) - L)); s = list(c[p:p+L])
    for i in range(len(s)):
        r = random.random()
        if   r < 0.003: s[i] = 'N'
        elif r < 0.005: s[i] = random.choice('RYUKM.-E')
        elif r < 0.02:  s[i] = s[i].lower()
        elif r < 0.03:  s[i] = random.choice('ACGT')
    s = ''.join(s)
    return rc(s.upper()) if random.random() < 0.5 and set(s) <= set('ACGT') else s
def fq(name, seqs):
    with gzip.open(name, 'wt') as f:
        for i, s in enumerate(seqs): f.write(f'@q{i}\n{s}\n+\n{"F" * len(s)}\n')
n = 3000
r1 = [read_from_ref() for _ in range(n)]; r2 = [read_from_ref() for _ in range(n)]
fq(f'{out}/r1.fq.gz', r1); fq(f'{out}/r2.fq.gz', r2)
fq(f'{out}/inter.fq.gz', [x for p in zip(r1, r2) for x in p])
fq(f'{out}/se.fq.gz', [read_from_ref() for _ in range(n)])
with open(f'{out}/batch.txt', 'w') as f:
    f.write(f'PE\t{out}/r1.fq.gz\t{out}/r2.fq.gz\nPEI\t{out}/inter.fq.gz\nSE\t{out}/se.fq.gz\nSE\t{out}/r2.fq.gz\n')
# background: an IUPAC-rich genome copy + reads
bgg = [noisy(len(c), exotic=0.02) if i % 2 else c for i, c in enumerate(ref)]
with open(f'{out}/bg.fasta', 'w') as f:
    for i, c in enumerate(bgg): f.write(f'>b{i}\n{c}\n')
with open(f'{out}/bg.txt', 'w') as f: f.write(f'{out}/bg.fasta\n{out}/r1.fq.gz\n{out}/se.fq.gz\n')
