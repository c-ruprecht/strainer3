import random, gzip, os, sys
random.seed(9)
out = sys.argv[1]; ref = sys.argv[2]; os.makedirs(out, exist_ok=True)
seq = []
for l in gzip.open(ref, 'rt'):
    if l.startswith('>'): seq.append('')
    else: seq[-1] += l.strip()
def mut(s, rate):
    s = list(s)
    for i in range(len(s)):
        r = random.random()
        if r < rate: s[i] = random.choice('ACGT')
        elif r < rate * 1.1: s[i] = random.choice('NRYKMU')
        elif r < rate * 1.3: s[i] = s[i].lower()
    return ''.join(s)
paths = []
for g in range(8):
    p = f'{out}/m{g}.' + ('fna.gz' if g % 2 else 'fasta')
    txt = ''.join(f'>c{i}\n{mut(c, random.choice([0.001, 0.01, 0.05]))}\nACG\n' for i, c in enumerate(seq))
    (gzip.open(p, 'wt') if p.endswith('.gz') else open(p, 'w')).write(txt)
    paths.append(p)
open(f'{out}/A.txt', 'w').write('\n'.join(paths[:5]) + '\n\n  \n' + paths[0] + '\n')            # duplicate -> skipped
open(f'{out}/C.txt', 'w').write(ref + '\n' + '\n'.join(paths[5:]) + '\n' + f'{out}/missing.fasta\n')  # ref skipped, missing -> warning
