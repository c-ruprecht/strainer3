"""Adversarial inputs for differential testing against v0.
Genomes: multi-line FASTA, lowercase runs, N runs, IUPAC codes (R,Y,K,M,S,W,B,D,H,V),
'U' bases, CRLF line endings, contigs shorter than k, empty lines, plain + gz files.
Query sets per k: sampled kmers forward + revcomp (duplicates of the same canonical
kmer), random absent kmers, extra columns before '#kmer', comment lines."""
import random, gzip, os, sys
random.seed(7)
out = sys.argv[1]; os.makedirs(out, exist_ok=True)
comp = str.maketrans('ACGT', 'TGCA')
rc = lambda s: s.translate(comp)[::-1]
def genome(n):
    s = []
    for _ in range(n):
        r = random.random()
        if   r < 0.003: s.append('N' * random.randint(1, 40))
        elif r < 0.006: s.append(random.choice('RYKMSWBDHV'))
        elif r < 0.008: s.append('U')
        elif r < 0.02:  s.append(''.join(random.choice('acgt') for _ in range(random.randint(5, 60))))
        else:           s.append(random.choice('ACGT'))
    return ''.join(s)
# low-complexity + palindrome-rich contig (ACGT / GAATTC repeats) for even-k palindromes
genomes = []
for gi in range(6):
    contigs = [genome(random.randint(2000, 20000)) for _ in range(random.randint(1, 5))]
    contigs.append(random.choice('ACGT') * 5)                # shorter than k
    contigs.append('GAATTCACGTACGT' * 50 + 'AAAAAAAAAA' * 10 + 'TTTTTTTTTT' * 10)
    genomes.append(contigs)
    name = f"{out}/g{gi}." + ("fna.gz" if gi % 2 else "fasta")
    lines = []
    for ci, c in enumerate(contigs):
        lines.append(f">c{ci} desc")
        w = random.choice([60, 80, 1000])
        lines += [c[i:i+w] for i in range(0, len(c), w)] + ([""] if ci % 2 else [])
    eol = "\r\n" if gi == 3 else "\n"
    data = eol.join(lines) + eol
    (gzip.open(name, 'wt') if name.endswith('.gz') else open(name, 'w')).write(data)
with open(f"{out}/targets.txt", 'w') as f:
    for gi in range(4): f.write(f"{out}/g{gi}." + ("fna.gz" if gi % 2 else "fasta") + "\n")
with open(f"{out}/background.txt", 'w') as f:
    for gi in range(4, 6): f.write(f"{out}/g{gi}." + ("fna.gz" if gi % 2 else "fasta") + "\n")
allseq = ''.join(c.upper() for g in genomes for c in g)
for k in (4, 15, 21, 31, 32, 33):
    rows = []
    for _ in range(3000):
        p = random.randint(0, len(allseq) - k); km = allseq[p:p+k]
        if not set(km) <= set('ACGT'): continue
        rows.append(km)
        if random.random() < 0.3: rows.append(rc(km))       # same canonical kmer again
    rows += [''.join(random.choice('ACGT') for _ in range(k)) for _ in range(500)]
    random.shuffle(rows)
    with gzip.open(f"{out}/q_k{k}.tsv.gz", 'wt') as f:
        f.write("id\t#kmer\tscore\n#comment line\n")
        for i, r in enumerate(rows): f.write(f"x{i}\t{r}\t{random.random():.3f}\n")
# query sets that must take the legacy string path: lowercase / N inside query kmers
for tag, mut in (("lower", str.lower), ("withN", lambda s: s[:5] + 'N' + s[6:])):
    rows = []
    for _ in range(2000):
        p = random.randint(0, len(allseq) - 21); km = allseq[p:p+21]
        if set(km) <= set('ACGT'): rows.append(mut(km) if random.random() < 0.2 else km)
    with open(f"{out}/q_{tag}.tsv", 'w') as f:
        f.write("#kmer\n" + "\n".join(rows) + "\n")
