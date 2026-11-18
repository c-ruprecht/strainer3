# C tool speed-ups — tests and benchmarks

Four tools were sped up; every change keeps the outputs identical to the
original (git `c18dbc4`) and is checked by a script in this directory.
The legacy `strain_detect` (strainer2 read-level detector, not built by the
Makefile, superseded by kmer_strain_detect) was removed.

| tool | test script | measured speed-up (2-core test box) |
|---|---|---|
| kmer_strain_detect_genome | `run_tests.sh` | 22–23x |
| kmer_strain_detect (FASTQ + .kmdb targets) | `run_tests_kmer_strain_detect.sh` | 6.8x (`-t 1`), 5.8x (`-t 2`) |
| kmer_scrub_count (+ kmer_scrub_streaming.c) | `run_tests_scrub_count.sh` | 1.4–1.8x |
| genome_compare_presence (new `-t`) | `run_tests_presence.sh` | 1.6x at `-t 2` |

Every script builds the original from git history next to the new source,
runs differential tests (real test data plus adversarial synthetic input),
a race test at `-t 1..64` on few cores, and a ThreadSanitizer build. Each
prints `ALL TESTS PASSED` / `TESTS FAILED` and exits 0 / 1.

---

## kmer_strain_detect_genome

`bash test/speedup/run_tests.sh` checks the rewritten `src/kmer_strain_detect_genome.c`
against the original implementation (taken from git commit `c18dbc4`, plus the one-byte
`-G` overflow fix so it can run). Every check compares the decompressed output TSV.

| test | what it covers |
|---|---|
| identity | test data: 18 D6331 targets + 52 F. prausnitzii pangenome background, `-t 1` / `-t 2` |
| edge | synthetic genomes with lowercase, N runs, IUPAC codes, `U`, CRLF, contigs shorter than k, palindromes; k = 4, 15, 21, 31, 32, 33; query files with lowercase/N (legacy path) |
| race | hit-dense scenario (~400k hits over 20 columns, same genome listed 4x), `-t 1..64` on few cores (forces preemption) x3 repeats |
| tsan | ThreadSanitizer build, 8 threads, must report 0 data races |

## Steps and timings

Scenario: `D6331_Faecalibacterium_prausnitzii.scrubbed_kmers.gz` (142,973 k-mers, k=31, from
`example.sh` steps 1–2) vs 18 targets + 52 background genomes (210 Mbp). 2-core Xeon
(Sapphire Rapids), zlib-ng, `hyperfine` means. Every step: identical TSV and identical .gz
bytes vs v0, 0 edge differences, 0 race mismatches, 0 TSan reports.

| step | change | `-t 1` | `-t 2` |
|---|---|---:|---:|
| v0 | original (+ `-G` overflow fix) | 49.7 s | 36.2 s |
| v1 | thread-local position counter (false sharing on `total_kmer_counts`); join worker threads | 50.4 s | 23.9 s |
| v2 | 2-bit rolling k-mers, integer open-addressing table | 6.69 s | 4.12 s |
| v3 | IUPAC/U fallback per window instead of per read | 5.74 s | 3.44 s |
| v4 | private cache-aligned count column per genome; 16-byte {key,row} slots | 5.43 s | 3.28 s |
| v5 | L2-resident bit prefilter in front of the table | 3.10 s | 2.15 s |
| v6 | row formatting into one buffer + `gzwrite` (same bytes) | 2.85 s | 1.79 s |
| v7 | batched, prefetched table lookups | 2.35 s | 1.56 s |
| v8 | opt-in `-z <level>` (default 9 = unchanged output) | **2.28 s** | **1.55 s** |

Overall: 21.8x (`-t 1`), 23.4x (`-t 2`). In v0 the 2-thread run burned 70 s of CPU for
49 s of work — that is the false sharing v1 removes.

Large k-mer set (2,864,528 k-mers = full global count table, `-t 2`): v0 88.8 s → v8 38.5 s
at `-z 9` (identical bytes), 11.5 s at `-z 6`, 6.4 s at `-z 1` (identical TSV, .gz 13% / 76%
larger). Peak RSS 1.25 → 1.44 GB.

## Where the time goes now (v8)

- Output: single-stream level-9 deflate is serial — 0.58 s of 1.55 s on the scrubbed set,
  34 of 38 s on the 2.86M set. Only a lower `-z` (or parallel gzip, which also changes the
  bytes) helps.
- Scan: ~1/3 gzip inflate + FASTA parsing (0.53 s single-threaded for all 70 files),
  ~2/3 the k-mer loop, whose hottest instruction is the prefilter bit load (L2 latency).

## Bugs found in the original

- `-G` crashes (`*** buffer overflow detected ***`) on glibc with fortify (Ubuntu default):
  `malloc(strlen(basename) + 9)` has no room for the NUL after `"backmeta_"`.
- Thread pool: workers were never joined; the pool mutex/struct were destroyed after a
  10 ms sleep while workers could still be exiting (3 ThreadSanitizer reports). Count
  updates themselves were race-free.

---

## kmer_strain_detect

`bash test/speedup/run_tests_kmer_strain_detect.sh`  (run `make` in src first)

Changes: query k-mers get a row each and every target a private,
cache-line-aligned count column (was one per-k-mer vector with a column per
thread: false sharing). FASTQ targets with k <= 31: lookups go to a 2-bit
index (bit prefilter, {key,row} slots, batched prefetched lookups) instead of
strcmp + string hash per window; window selection, orientation and pair dedup
are unchanged (k > 31 keeps the string path). .kmdb targets: the same sorted
block merge, now writing the target's column; with fewer targets than threads
a database's blocks are split over the spare threads (a key lives in one
block, so rows are disjoint). Output rows formatted into one buffer per row
(same gz bytes). Pool workers are joined (the original freed the pool after a
10 ms sleep: 6 ThreadSanitizer reports). New opt-in `-z <level>` (default 9,
unchanged bytes).

Output identity: exact for every -t (decompressed TSV, sorted stderr
statistics, exit code); FASTQ and .kmdb targets built from the same reads
give identical columns including total_evaluated.

| 2 PE metagenomes (500k pairs each) + their 2 .kmdb, -G background, 143k k-mers | `-t 1` | `-t 2` |
|---|---:|---:|
| original | 50.7 s | 28.0 s |
| new | 7.5 s | 4.8 s |

FASTQ-only (the 2 PE metagenomes, `-t 2`): 27.8 s -> 3.9 s.
.kmdb targets were already a sorted merge: ~3 s -> 2.7 s; what is left there
is the level-9 gzip of the output table (one .kmdb, 573k k-mers: 9.0 s at the
default level, 0.94 s with `-z 1`, identical TSV) and the block decoding.

---

## kmer_scrub_count

`bash test/speedup/run_tests_scrub_count.sh`  (run `make` in src first, for kmer_metagenome_db)

The FASTA/FASTQ scan was already 2-bit. Changes: {key,row} slots with a bit
prefilter and batched prefetched lookups; per-worker private scratch counts by
row (were a column in every k-mer's malloc'd vector, reached through two
dependent loads); the "reference count > 0" flag precomputed per row; global
counts accumulated per list in a dense array and folded into the hash column
after the join (unsigned adds, identical); the sorted reference k-mers and the
index built once instead of once per -A/-B/-C list; global-count and presence
rows formatted into a 4 MB buffer handed to zstd (same byte stream, same level).

Output identity: `-t 1`: all three files byte-identical. More threads: the
global counts byte-identical; summary and presence compared with scrub_ids
mapped to sample names (`canon_ksc.py`) because the original assigns scrub_ids
in sample-completion order.

| example.sh step 1 (52 genomes, 1 .kmdb, 18 drug genomes) | `-t 1` | `-t 2` |
|---|---:|---:|
| original, no -P | 25.8 s | 16.9 s |
| new, no -P | 14.0 s | 12.4 s |
| original, -P 0 | 29.8 s | 25.5 s |
| new, -P 0 | 18.9 s | 18.3 s |

Now: serial work dominates. zstd level 9 is single-threaded: 3.2 s for the
global counts (115 MB) and 3.7 s for the presence file (138 MB) — as fast as
`zstd -9` itself. Multithreaded zstd would halve that but produces different
.zst bytes (same content), so it was not switched on. Then the reference hash
build (1.2 s, string hash) and, with -P, the single presence writer thread.

---

## genome_compare_presence

`bash test/speedup/run_tests_presence.sh`

Was fully serial. New `-t` (default 4; the output is identical for any value):
the default sorted mode builds the index in 64 hash shards — genomes read in
parallel batches, each thread inserts the batch in strain_id order into the
shards it owns, so every id list and per-genome count is unchanged. `-U` (hash
order) keeps one table filled in the original order (readers parse ahead, one
thread inserts with prefetch), so even its order is identical. The sort keeps
keys inline and runs as a parallel bucket sort; the sweep runs in parallel
parts written in order, with private intersection matrices and atomic OR for
shared bitset words.

| 25 F. prausnitzii genomes, 23.5M distinct 31-mers | time |
|---|---:|
| original | 37.9 s |
| new `-t 1` | 33.3 s |
| new `-t 2` | 23.5 s |
| original `-U` | 25.8 s |
| new `-t 2 -U` | 19.9 s |

Now: memory-bound. ~7–11 s of each run is kernel time zeroing the pages of the
growing 2 GB hash; the rest is random access into it (index build and sweep).
These parts scale with cores; on 2 cores the gain is capped.
