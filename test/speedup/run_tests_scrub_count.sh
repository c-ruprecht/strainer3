#!/bin/bash
# kmer_scrub_count: new src/kmer_scrub_count.c + kmer_scrub_streaming.c vs the
# originals (git $REF). -t 1: every output byte-identical. More threads: the
# global counts byte-identical; summary/presence compared with scrub_ids mapped
# to sample names (canon_ksc.py), since the original assigns ids in completion order.
set -u; source "$(dirname "$0")/common.sh"
printf 'int main(void){return 0;}' | gcc -x c - -lzstd -o /dev/null 2>/dev/null || { echo "needs libzstd"; exit 1; }
kl="$src/kmdb.c $libs -lz -lzstd -lm -lpthread"
echo "== build"
gitsrc kmer_scrub_count.c > "$W/ksc_ref.c"; gitsrc kmer_scrub_streaming.c > "$W/kss_ref.c"
gcc -std=gnu17 -O2 -I"$src" -DHAVE_ZSTD -o "$W/ksc_ref" "$W/ksc_ref.c" "$W/kss_ref.c" $kl 2>/dev/null &&
gcc -std=gnu17 -O2 -I"$src" -DHAVE_ZSTD -o "$W/ksc_new" "$src/kmer_scrub_count.c" "$src/kmer_scrub_streaming.c" $kl 2>/dev/null || { echo "build failed"; exit 1; }
cd "$TD"
R=D6331_genomes/D6331_Faecalibacterium_prausnitzii.fasta.gz
python3 "$here/make_sd_edge.py" "$W/sdedge" >/dev/null; ER=$W/sdedge/ref.fna.gz
python3 "$here/make_ksc_edge.py" "$W/kedge" $ER >/dev/null; E=$W/kedge
# a .kmdb sample (needs src/kmer_metagenome_db: run make in src first)
[ -x "$src/kmer_metagenome_db" ] && "$src/kmer_metagenome_db" -m sample -O $E -n readsdb -i $W/sdedge/r1.fq.gz,$W/sdedge/r2.fq.gz >/dev/null 2>&1
[ -f $E/readsdb.kmdb ] || echo "  (note: src/kmer_metagenome_db not built: edge set runs without a .kmdb)"
printf "$W/sdedge/r1.fq.gz\n$W/sdedge/se.fq.gz\n" > $E/B.txt; [ -f $E/readsdb.kmdb ] && echo $E/readsdb.kmdb >> $E/B.txt
files() { ls $1.global_counts.tsv.zst $1.summary.tsv $1.presence.tsv.zst 2>/dev/null; }
exact() { local lab=$1; shift; rm -f $W/r.* $W/n.*
  "$W/ksc_ref" "$@" -t 1 -o $W/r 2>/dev/null; "$W/ksc_new" "$@" -t 1 -o $W/n 2>/dev/null
  local ok=1; for f in $(files $W/r); do cmp -s $f ${f/\/r./\/n.} || ok=0; done
  [ "$(files $W/r | wc -l)" = "$(files $W/n | wc -l)" ] || ok=0
  [ $ok = 1 ] && echo "  ok      $lab" || { echo "  DIFFERS $lab"; fail=1; }; }
multi() { local lab=$1 t=$2; shift 2; rm -f $W/r.* $W/n.*
  "$W/ksc_ref" "$@" -t 1 -o $W/r 2>/dev/null; "$W/ksc_new" "$@" -t $t -o $W/n 2>/dev/null
  if cmp -s $W/r.global_counts.tsv.zst $W/n.global_counts.tsv.zst && cmp -s <(python3 "$here/canon_ksc.py" $W/r) <(python3 "$here/canon_ksc.py" $W/n)
  then echo "  ok      $lab (-t $t)"; else echo "  DIFFERS $lab (-t $t)"; fail=1; fi; }
L="-A lists/genomes_to_scrub.txt -B lists/metagenomes_to_scrub.txt -C lists/drug_scrub.txt"
echo "== identity, -t 1, every byte"
exact "example.sh step 1"                  -r $R $L -T 0.995
exact "  + presence -P 0"                  -r $R $L -T 0.995 -P 0
exact "  + -P 5 and coverage gate -T 0.045" -r $R $L -T 0.045 -P 5
exact "single files (-A genome, -B .kmdb)" -r $R -A scrub_pangenome/MGYG000006888.fna.gz -B scrub_metagenomes/ATCC-20-staggerd.kmdb -P 0
exact "edge: IUPAC reference, reads, kmdb, duplicate/missing/reference-skip" -r $ER -A $E/A.txt -B $E/B.txt -C $E/C.txt -P 0
exact "edge: -P 2 -T 0.5"                  -r $ER -A $E/A.txt -B $E/B.txt -C $E/C.txt -P 2 -T 0.5
echo "== multi-threaded"
multi "example.sh step 1 -P 0" 2 -r $R $L -T 0.995 -P 0
multi "edge -P 0"              4 -r $ER -A $E/A.txt -B $E/B.txt -C $E/C.txt -P 0
echo "== race: -t 1..64 (example.sh step 1, -P 0)"
"$W/ksc_ref" -r $R $L -T 0.995 -P 0 -t 1 -o $W/race_ref 2>/dev/null; python3 "$here/canon_ksc.py" $W/race_ref > $W/race_ref.canon; n=0; d=0
for t in 1 2 3 4 8 16 32 64; do
  "$W/ksc_new" -r $R $L -T 0.995 -P 0 -t $t -o $W/race_new 2>/dev/null; n=$((n+1))
  { cmp -s $W/race_ref.global_counts.tsv.zst $W/race_new.global_counts.tsv.zst && cmp -s <(python3 "$here/canon_ksc.py" $W/race_new) $W/race_ref.canon; } || { d=$((d+1)); echo "  MISMATCH -t $t"; }
done; echo "  $n runs, $d mismatches"; [ $d -eq 0 ] || fail=1
echo "== ThreadSanitizer"
if gcc -std=gnu17 -O1 -g -fsanitize=thread -I"$src" -DHAVE_ZSTD -o "$W/ksc_tsan" "$src/kmer_scrub_count.c" "$src/kmer_scrub_streaming.c" $kl 2>/dev/null; then
  TSAN_OPTIONS="halt_on_error=0" "$W/ksc_tsan" -r $ER -A $E/A.txt -B $E/B.txt -C $E/C.txt -P 0 -t 8 -o $W/ts 2> $W/tsan.log
  nr=$(grep -c "WARNING: ThreadSanitizer" $W/tsan.log); echo "  reports: $nr"; [ "$nr" -eq 0 ] || fail=1
else echo "  skipped (no -fsanitize=thread)"; fi
report
