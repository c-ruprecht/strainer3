#!/usr/bin/env python3
"""Simulate a relative of a genome and paired-end reads from it.
  mutate: diffuse SNPs (uniform) or tract model (clonal frame + recombined tracts), 5 % of events are 1-bp indels
  reads : paired-end, per-base qualities that decay along the read (R2 worse), substitution errors drawn from the
          qualities, adapter read-through when the insert is shorter than the read
"""
import argparse, gzip, sys
import numpy as np

COMP = bytes.maketrans(b'ACGTN', b'TGCAN')
AD1 = b'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
AD2 = b'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'


def read_fasta(path):
    recs, name, buf = [], None, []
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None: recs.append((name, ''.join(buf).upper()))
                name, buf = line[1:].split()[0], []
            else: buf.append(line.strip())
    if name is not None: recs.append((name, ''.join(buf).upper()))
    return recs


def mutate(recs, div, model, rng, tract_len=5000, tract_div=0.03, frame_div=0.00005):
    """returns mutated records and the realised divergence"""
    out, nmut, ntot = [], 0, 0
    for name, seq in recs:
        n = len(seq)
        if model == 'diffuse':
            p = np.full(n, div)
        else:   # tracts: fraction f of the genome at tract_div, the rest at frame_div, mean = div
            f = min(1.0, max(0.0, (div - frame_div) / (tract_div - frame_div)))
            p = np.full(n, frame_div)
            n_tr = int(round(f * n / tract_len))
            for s in rng.integers(0, max(n - tract_len, 1), n_tr):
                p[s:s + tract_len] = tract_div
        hit = np.flatnonzero(rng.random(n) < p)
        b = bytearray(seq.encode())
        dels = set()
        ins = {}
        for i in hit:
            r = rng.random()
            if r < 0.95:
                cur = chr(b[i]); alt = [x for x in 'ACGT' if x != cur]
                b[i] = ord(alt[rng.integers(0, len(alt))])
            elif r < 0.975: dels.add(int(i))
            else: ins[int(i)] = 'ACGT'[rng.integers(0, 4)]
        if dels or ins:
            parts = []
            for i in range(n):
                if i in dels: continue
                parts.append(chr(b[i]))
                if i in ins: parts.append(ins[i])
            new = ''.join(parts)
        else: new = b.decode()
        out.append((name, new)); nmut += len(hit); ntot += n
    return out, nmut / max(ntot, 1)


def qual_profile(L, mate, rng, n):
    """n x L matrix of phred scores"""
    end = 33.0 if mate == 0 else 30.0
    mean = np.linspace(37.5, end, L)[None, :]
    q = mean + rng.normal(0, 3.0, (n, L)) + rng.normal(0, 2.0, (n, 1))
    bad = rng.random(n) < (0.03 if mate == 0 else 0.06)        # reads with a poor tail
    start = rng.integers(L // 3, L, n)
    tail = np.arange(L)[None, :] >= start[:, None]
    q = np.where(bad[:, None] & tail, q - 18, q)
    return np.clip(np.round(q), 2, 40).astype(np.int16)


def simulate(recs, depth, read_len, ins_mean, ins_sd, rng, prefix, out1, out2):
    total = sum(len(s) for _, s in recs)
    n_pairs = int(depth * total / (2 * read_len))
    lens = np.array([len(s) for _, s in recs], float)
    which = rng.choice(len(recs), n_pairs, p=lens / lens.sum())
    ins = np.clip(np.round(rng.normal(ins_mean, ins_sd, n_pairs)), 35, None).astype(int)
    Q = [qual_profile(read_len, 0, rng, n_pairs), qual_profile(read_len, 1, rng, n_pairs)]
    bases = np.frombuffer(b'ACGT', dtype=np.uint8)
    seqb = [s.encode() for _, s in recs]
    for i in range(n_pairs):
        s = seqb[which[i]]
        L = min(ins[i], len(s))
        st = rng.integers(0, len(s) - L + 1)
        frag = s[st:st + L]
        if rng.random() < 0.5: frag = frag.translate(COMP)[::-1]
        mates = (frag, frag.translate(COMP)[::-1])
        for m, (out, ad) in enumerate(((out1, AD1), (out2, AD2))):
            r = mates[m][:read_len]
            if len(r) < read_len:      # adapter read-through, then random bases
                r = (r + ad + bytes(bases[rng.integers(0, 4, read_len)]))[:read_len]
            a = np.frombuffer(r, dtype=np.uint8).copy()
            q = Q[m][i]
            err = rng.random(read_len) < 10.0 ** (-q / 10.0)
            if err.any():
                for j in np.flatnonzero(err):
                    alt = bases[bases != a[j]]
                    a[j] = alt[rng.integers(0, len(alt))]
            out.write(b'@%s_%d/%d\n%s\n+\n%s\n' % (prefix.encode(), i, m + 1, a.tobytes(), (q + 33).astype(np.uint8).tobytes()))
    return n_pairs


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--genome', required=True); ap.add_argument('--div', type=float, default=0.0)
    ap.add_argument('--model', choices=['diffuse', 'tracts'], default='diffuse')
    ap.add_argument('--depth', type=float, required=True); ap.add_argument('--read_len', type=int, default=150)
    ap.add_argument('--ins_mean', type=float, required=True); ap.add_argument('--ins_sd', type=float, required=True)
    ap.add_argument('--seed', type=int, default=1); ap.add_argument('--out', required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    recs = read_fasta(a.genome)
    real = 0.0
    if a.div > 0: recs, real = mutate(recs, a.div, a.model, rng)
    with gzip.open(a.out + '_R1.fq.gz', 'wb', compresslevel=1) as o1, gzip.open(a.out + '_R2.fq.gz', 'wb', compresslevel=1) as o2:
        n = simulate(recs, a.depth, a.read_len, a.ins_mean, a.ins_sd, rng, 'r', o1, o2)
    print(f'{a.out}\tpairs={n}\trealised_div={real:.5f}', file=sys.stderr)
