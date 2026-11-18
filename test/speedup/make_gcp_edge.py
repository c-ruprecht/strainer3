"""80 related synthetic genomes (5 ancestors, point mutations, indels, exotic bases)
so k-mers are shared by 1..80 strains and N >= 64 triggers the bitset path."""
import random, gzip, os, sys
random.seed(5)
out = sys.argv[1]; os.makedirs(out + '/dir', exist_ok=True)
anc = [''.join(random.choice('ACGT') for _ in range(40000)) for _ in range(5)]
def mutate(s, rate):
    s = list(s)
    for i in range(len(s)):
        r = random.random()
        if   r < rate:          s[i] = random.choice('ACGT')
        elif r < rate * 1.05:   s[i] = random.choice('NRYKMSWU')
        elif r < rate * 1.2:    s[i] = s[i].lower()
    return ''.join(s)
paths = []
for g in range(80):
    a = anc[g % 5] if g < 70 else anc[0][:20000] + anc[1][:20000]
    seq = mutate(a, random.choice([0.0005, 0.002, 0.01, 0.03]))
    cuts = sorted(random.sample(range(1000, len(seq) - 1000), 3))
    contigs = [seq[:cuts[0]], seq[cuts[0]:cuts[1]], 'ACG', seq[cuts[1]:]]
    ext = ['.fna.gz', '.fa', '.fasta', '.fna'][g % 4]
    d = out + '/dir' if g >= 74 else out
    p = f'{d}/g{g:02d}{ext}'
    eol = '\r\n' if g == 7 else '\n'
    txt = eol.join(f'>c{i}{eol}' + eol.join(c[j:j+80] for j in range(0, len(c), 80)) for i, c in enumerate(contigs)) + eol
    (gzip.open(p, 'wt') if p.endswith('.gz') else open(p, 'w')).write(txt)
    if g < 74: paths.append(p)
with open(out + '/list.txt', 'w') as f:
    f.write('# genome list with a comment, a blank line, a missing file and a directory\n\n')
    for p in paths[:40]: f.write(p + '\n')
    f.write(out + '/does_not_exist.fna\n')
    f.write(f'  {paths[40]}  \n')
    for p in paths[41:]: f.write(p + '\n')
    f.write(out + '/dir\n')
