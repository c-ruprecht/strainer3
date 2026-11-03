#!/usr/bin/env python3
"""
Synthetic FASTQ cases for kmer_metagenome_db --pair-dedup, plus an
independent reference count for each case.

  make_cases.py <outdir> [k]

Writes, per case <c>:
  <c>_R1.fq / <c>_R2.fq      mates (same record names, R2 = mate strand)
  <c>_I.fq                   the same pairs interleaved
  <c>.dedup.tsv              expected counts with pair dedup  (kmer \t count)
  <c>.nodedup.tsv            expected counts without          (kmer \t count)
Rows are sorted by canonical k-mer, as kmdb_dump prints them.
Canonical = max(kmer, revcomp) as strings, which equals max() on the
2-bit A0C1G2T3 encoding the C code uses.
"""
import os
import random
import sys
from collections import Counter

COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def kmers(seq, k):
    """canonical k-mers of every window without N"""
    out = []
    for i in range(len(seq) - k + 1):
        w = seq[i:i + k]
        if "N" in w:
            continue
        out.append(max(w, rc(w)))
    return out


def rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def write_case(outdir, name, pairs, k):
    with open(f"{outdir}/{name}_R1.fq", "w") as f1, \
         open(f"{outdir}/{name}_R2.fq", "w") as f2, \
         open(f"{outdir}/{name}_I.fq", "w") as fi:
        for i, (a, b) in enumerate(pairs):
            r1 = f"@frag{i} 1:N:0:ACGT\n{a}\n+\n{'I' * len(a)}\n"
            r2 = f"@frag{i} 2:N:0:ACGT\n{b}\n+\n{'I' * len(b)}\n"
            f1.write(r1); f2.write(r2); fi.write(r1 + r2)
    ded, nod = Counter(), Counter()
    for a, b in pairs:
        ka, kb = kmers(a, k), kmers(b, k)
        nod.update(ka); nod.update(kb)
        ded.update(set(ka) | set(kb))
    for tag, cnt in (("dedup", ded), ("nodedup", nod)):
        with open(f"{outdir}/{name}.{tag}.tsv", "w") as f:
            for km in sorted(cnt):
                f.write(f"{km}\t{cnt[km]}\n")


def main():
    outdir = sys.argv[1]
    k = int(sys.argv[2]) if len(sys.argv) > 2 else 31
    os.makedirs(outdir, exist_ok=True)
    rng = random.Random(42)
    L = 100

    # 1. full overlap: R2 is the exact reverse complement of R1
    a = rand_seq(rng, L)
    write_case(outdir, "full", [(a, rc(a))], k)

    # 2. no overlap: mates from unrelated sequence
    write_case(outdir, "none", [(rand_seq(rng, L), rand_seq(rng, L))], k)

    # 3. partial overlap: 150 bp insert, mates share frag[50:100]
    frag = rand_seq(rng, 150)
    write_case(outdir, "partial", [(frag[:100], rc(frag[50:]))], k)

    # 4. internal repeat of one k-mer inside R1, unrelated R2
    x = rand_seq(rng, k)
    r1 = rand_seq(rng, 20) + x + rand_seq(rng, 10) + x
    r1 += rand_seq(rng, L - len(r1))
    write_case(outdir, "repeat", [(r1, rand_seq(rng, L))], k)

    # 5. bulk: many pairs from a few small "genomes", random inserts
    #    (overlapping and not), occasional N, short mates (< k) and k-mers
    #    shared across fragments
    genomes = [rand_seq(rng, 3000) for _ in range(3)]
    pairs = []
    for _ in range(3000):
        g = rng.choice(genomes)
        ins = rng.randint(40, 260)
        s = rng.randint(0, len(g) - ins)
        f = g[s:s + ins]
        if rng.random() < 0.5:
            f = rc(f)
        la = min(rng.choice([L, L, L, 25, 70]), len(f))
        lb = min(rng.choice([L, L, L, 25, 70]), len(f))
        m1, m2 = list(f[:la]), list(rc(f[len(f) - lb:]))
        if rng.random() < 0.1:
            m1[rng.randrange(len(m1))] = "N"
        pairs.append(("".join(m1), "".join(m2)))
    write_case(outdir, "bulk", pairs, k)

    # 6. mates out of sync (names) and different read counts: error cases
    with open(f"{outdir}/bad_R1.fq", "w") as f1, open(f"{outdir}/bad_R2.fq", "w") as f2:
        for i in range(3):
            s = rand_seq(rng, L)
            f1.write(f"@frag{i}/1\n{s}\n+\n{'I' * L}\n")
            f2.write(f"@frag{i + (i == 2)}/2\n{rc(s)}\n+\n{'I' * L}\n")
    with open(f"{outdir}/short_R2.fq", "w") as f2:
        s = rand_seq(rng, L)
        f2.write(f"@frag0 2:N:0:ACGT\n{s}\n+\n{'I' * L}\n")
    with open(f"{outdir}/odd_I.fq", "w") as fi:
        for i in range(3):
            s = rand_seq(rng, L)
            fi.write(f"@frag{i // 2}\n{s}\n+\n{'I' * L}\n")


if __name__ == "__main__":
    main()
