#!/bin/bash
# Differential + race-safety tests for src/kmer_strain_detect_genome.
#
# The reference is the original string-hash implementation taken from git
# history (commit $REF), with the one-byte -G overflow fix applied, so it can
# run with background genomes. Every test compares the decompressed TSV of the
# new binary against the reference.
#
#   1. identity   test data (D6331 targets + F. prausnitzii pangenome background), -t 1/2
#   2. edge       synthetic adversarial genomes (lowercase, N runs, IUPAC, U, CRLF,
#                 short contigs, palindromes) x k in {4,15,21,31,32,33} + legacy-path queries
#   3. race       hit-dense scenario, -t 1..64 (oversubscribed) x3 repeats
#   4. tsan       ThreadSanitizer build, 8 threads; must report 0 data races
#
# usage: bash run_tests.sh        (from anywhere; needs gcc, zlib, python3, git)
set -u
here=$(cd "$(dirname "$0")" && pwd); repo=$(cd "$here/../.." && pwd)
src=$repo/src; T=$repo/test; W=${WORK:-$here/work}; REF=${REF:-c18dbc4}
mkdir -p "$W"; fail=0
libs="$src/BIO_sequence.c $src/BIO_hash.c $src/genome_compare.c $src/up2bit.c"
cc() { gcc -std=gnu17 -O2 -I"$src" "$@" $libs -lz -lm -lpthread 2>/dev/null; }

echo "== build"
git -C "$repo" show "$REF:src/kmer_strain_detect_genome.c" \
  | sed 's/malloc(strlen(bg\[i\].basename) + 9)/malloc(strlen(bg[i].basename) + 10)/' > "$W/ref.c" || exit 1
cc -o "$W/ref" "$W/ref.c" || { echo "reference build failed"; exit 1; }
cc -o "$W/new" "$src/kmer_strain_detect_genome.c" || { echo "new build failed"; exit 1; }

cmp_run() {   # cmp_run <label> <ref-args> -- <new-args...>
  "$W/ref" $REFARGS -t 1 -o "$W/r.gz" 2>/dev/null
  "$W/new" $REFARGS "$@" -o "$W/n.gz" 2>/dev/null
  if cmp -s <(zcat "$W/r.gz") <(zcat "$W/n.gz"); then return 0; else return 1; fi
}

echo "== 1. identity on test data"
ls "$T"/D6331_genomes/*.fasta.gz > "$W/targets.txt"; ls "$T"/scrub_pangenome/*.fna.gz > "$W/background.txt"
python3 "$here/make_dense_queries.py" "$W/q_test.tsv.gz" 31 "$T/D6331_genomes/D6331_Faecalibacterium_prausnitzii.fasta.gz" "$(head -1 "$W/background.txt")"
REFARGS="-k $W/q_test.tsv.gz -B $W/targets.txt -G $W/background.txt"
"$W/ref" $REFARGS -t 1 -o "$W/ref_test.gz" 2>/dev/null; zcat "$W/ref_test.gz" > "$W/ref_test.tsv"
for t in 1 2; do
  "$W/new" $REFARGS -t $t -o "$W/n.gz" 2>/dev/null
  zcat "$W/n.gz" | cmp -s - "$W/ref_test.tsv" && echo "  t=$t identical" || { echo "  t=$t DIFFERS"; fail=1; }
done

echo "== 2. edge cases"
python3 "$here/make_edge.py" "$W/edge"
n=0; d=0
for q in "$W"/edge/q_*; do for t in 1 3; do
  REFARGS="-k $q -B $W/edge/targets.txt -G $W/edge/background.txt"; n=$((n+1))
  cmp_run -t $t || { d=$((d+1)); echo "  DIFF $(basename "$q") t=$t"; }
done; done
echo "  $n runs, $d differences"; [ $d -eq 0 ] || fail=1

echo "== 3. race (hit-dense, oversubscribed threads)"
F=$T/D6331_genomes/D6331_Faecalibacterium_prausnitzii.fasta.gz
( for i in 1 2 3 4; do echo "$F"; done; sed -n 1,8p "$W/background.txt" ) > "$W/race_t.txt"; sed -n 9,16p "$W/background.txt" > "$W/race_b.txt"
REFARGS="-k $W/q_test.tsv.gz -B $W/race_t.txt -G $W/race_b.txt"
"$W/ref" $REFARGS -t 1 -o "$W/race_ref.gz" 2>/dev/null; zcat "$W/race_ref.gz" > "$W/race_ref.tsv"
n=0; d=0
for t in 1 2 3 4 8 16 32 64; do for r in 1 2 3; do
  "$W/new" $REFARGS -t $t -o "$W/n.gz" 2>/dev/null; n=$((n+1))
  zcat "$W/n.gz" | cmp -s - "$W/race_ref.tsv" || { d=$((d+1)); echo "  MISMATCH t=$t rep=$r"; }
done; done
echo "  $n runs, $d mismatches"; [ $d -eq 0 ] || fail=1

echo "== 4. ThreadSanitizer"
if gcc -std=gnu17 -O1 -g -fsanitize=thread -I"$src" -o "$W/new_tsan" "$src/kmer_strain_detect_genome.c" $libs -lz -lm -lpthread 2>/dev/null; then
  ( echo "$F"; echo "$F"; sed -n 1,2p "$W/background.txt" ) > "$W/ts_t.txt"; sed -n 3,4p "$W/background.txt" > "$W/ts_b.txt"
  TSAN_OPTIONS="halt_on_error=0" "$W/new_tsan" -k "$W/q_test.tsv.gz" -B "$W/ts_t.txt" -G "$W/ts_b.txt" -t 8 -o "$W/ts.gz" 2> "$W/tsan.log"
  nr=$(grep -c "WARNING: ThreadSanitizer" "$W/tsan.log")
  echo "  reports: $nr"; [ "$nr" -eq 0 ] || { fail=1; grep -A6 "WARNING: ThreadSanitizer" "$W/tsan.log" | head -20; }
else echo "  skipped (compiler has no -fsanitize=thread)"; fi

[ $fail -eq 0 ] && echo "ALL TESTS PASSED" || echo "TESTS FAILED"
exit $fail
