#!/bin/bash
# kmer_strain_detect: new src/kmer_strain_detect.c vs the original (git $REF).
# Output is deterministic for any -t (one column per target), so every check
# is exact: decompressed TSV, stderr statistics (sorted) and exit code.
# Also checks that FASTQ and .kmdb targets built from the same reads agree.
# Needs libzstd and src/kmer_metagenome_db (run make in src first).
set -u; source "$(dirname "$0")/common.sh"
kl="$src/kmdb.c $libs -lz -lzstd -lm -lpthread"
[ -x "$src/kmer_metagenome_db" ] || { echo "build src/kmer_metagenome_db first (make)"; exit 1; }
echo "== build"
gitsrc kmer_strain_detect.c > "$W/ksd_ref.c"
gcc -std=gnu17 -O2 -I"$src" -DHAVE_ZSTD -o "$W/ksd_ref" "$W/ksd_ref.c" $kl 2>/dev/null &&
gcc -std=gnu17 -O2 -I"$src" -DHAVE_ZSTD -o "$W/ksd_new" "$src/kmer_strain_detect.c" $kl 2>/dev/null || { echo "build failed"; exit 1; }
TM=$TD/target_metagenomes; M=$W/mg; mkdir -p $M
zcat $TM/D6331-231101-C4_R1_001.fastq.downsampled.gz | head -400000 | gzip -1 > $M/pos_R1.fq.gz
zcat $TM/D6331-231101-C4_R2_001.fastq.downsampled.gz | head -400000 | gzip -1 > $M/pos_R2.fq.gz
zcat $TM/ATCC-10-even_S21_L001_R1_001.fastq.downsampled.gz | head -200000 | gzip -1 > $M/neg_R1.fq.gz
zcat $TM/ATCC-10-even_S21_L001_R2_001.fastq.downsampled.gz | head -200000 | gzip -1 > $M/neg_R2.fq.gz
paste <(zcat $M/pos_R1.fq.gz | paste - - - -) <(zcat $M/pos_R2.fq.gz | paste - - - -) | tr '\t' '\n' | gzip -1 > $M/pos_inter.fq.gz
"$src/kmer_metagenome_db" -m sample -O $M -n pos -i $M/pos_R1.fq.gz,$M/pos_R2.fq.gz >/dev/null 2>&1
"$src/kmer_metagenome_db" -m sample -O $M -n neg -i $M/neg_R1.fq.gz,$M/neg_R2.fq.gz >/dev/null 2>&1
printf "PE\t$M/pos_R1.fq.gz\t$M/pos_R2.fq.gz\nPE\t$M/neg_R1.fq.gz\t$M/neg_R2.fq.gz\n$M/pos.kmdb\n$M/neg.kmdb\n" > $W/targets.txt
printf "PEI\t$M/pos_inter.fq.gz\nSE\t$M/neg_R1.fq.gz\nPE\t$M/pos_R1.fq.gz\t$M/pos_R2.fq.gz\n" > $W/modes.txt
printf "SE\t$M/pos_R2.fq.gz\n$M/neg.kmdb\n" > $W/bg.txt
echo $M/pos.kmdb > $W/one.txt
( for i in 1 2 3; do cat $W/modes.txt; done; echo $M/pos.kmdb; echo $M/neg.kmdb ) > $W/race.txt
REF=$TD/D6331_genomes/D6331_Faecalibacterium_prausnitzii.fasta.gz
SEQ="zcat $REF | grep -v '>' | tr -d '\n'"
(echo -e "#kmer\tinfo"; eval $SEQ | fold -w 31 | awk 'NR%3==0 && length($0)==31 {print $0"\tx"}') > $W/k_dense.tsv
(echo "# scrub stats"; eval $SEQ | fold -w 21 | awk 'NR%7==0 && length($0)==21') | gzip -1 > $W/k21.gz
eval $SEQ | fold -w 33 | awk 'NR%20==0 && length($0)==33' > $W/k33.txt
: > $W/k_empty.txt
python3 "$here/make_sd_edge.py" "$W/edge" >/dev/null; E=$W/edge
paste <(zcat $E/r1.fq.gz | paste - - - -) <(zcat $E/r2.fq.gz | paste - - - -) | awk -F'\t' '{print $1"/1\n"$2"\n+\n"$4"\n"$5"/2\n"$6"\n+\n"$8}' | gzip -1 > $E/inter_ok.fq.gz
"$src/kmer_metagenome_db" -m sample -O $E -n readsdb -i $E/r1.fq.gz,$E/r2.fq.gz >/dev/null 2>&1
printf "PE\t$E/r1.fq.gz\t$E/r2.fq.gz\nPEI\t$E/inter_ok.fq.gz\nSE\t$E/se.fq.gz\n$E/readsdb.kmdb\n" > $W/edge_targets.txt
(echo "#kmer"; zcat $E/ref.fna.gz | grep -v ">" | tr -d '\n' | fold -w 31 | awk 'NR%2==0') > $W/k_edge.txt
same() { local lab=$1 t=$2; shift 2; rm -f $W/r.gz $W/n.gz
  "$W/ksd_ref" "$@" -t 1 -o $W/r.gz 2> $W/r.err; local rr=$?
  "$W/ksd_new" "$@" -t $t -o $W/n.gz 2> $W/n.err; local nr=$?
  local f='cat'; [ $rr != 0 ] && f='grep ^error'   # before an error exit, which stats lines print depends on timing
  if [ $rr = $nr ] && cmp -s <(zcat $W/r.gz 2>/dev/null) <(zcat $W/n.gz 2>/dev/null) &&
     cmp -s <(grep -v "^output" $W/r.err | $f | sort) <(grep -v "^output" $W/n.err | $f | sort)
  then echo "  ok      $lab (-t $t)"; else echo "  DIFFERS $lab (-t $t)"; fail=1; fi; }
echo "== identity (TSV, stderr statistics, exit code)"
for t in 1 2; do same "2 PE fastq + 2 kmdb, -G background" $t -k $W/k_dense.tsv -B $W/targets.txt -G $W/bg.txt; done
same "PEI/SE/PE"                     3 -k $W/k_dense.tsv -B $W/modes.txt -G $W/bg.txt
same "  --no-pair-dedup"             3 -k $W/k_dense.tsv -B $W/modes.txt --no-pair-dedup
same "headerless k=21 + kmdb (k mismatch: error)" 2 -k $W/k21.gz -B $W/targets.txt
same "headerless k=21, fastq"        3 -k $W/k21.gz -B $W/modes.txt
same "k=33 (string path)"            3 -k $W/k33.txt -B $W/modes.txt
same "empty k-mer file"              2 -k $W/k_empty.txt -B $W/modes.txt
same "edge: N/lowercase/IUPAC, invalid queries, kmdb" 3 -k $W/k_edge.txt -B $W/edge_targets.txt
same "edge --no-pair-dedup"          3 -k $W/k_edge.txt -B $W/edge_targets.txt --no-pair-dedup
echo "== FASTQ and .kmdb of the same reads agree"
"$W/ksd_new" -k $W/k_dense.tsv -B $W/targets.txt -t 2 -o $W/fk.gz 2>/dev/null
d=$(zcat $W/fk.gz | awk -F'\t' 'NR>1 && ($2!=$4 || $3!=$5) {n++} END{print n+0}')
[ "$d" = 0 ] && echo "  ok      fastq == kmdb columns (incl. total_evaluated)" || { echo "  DIFFERS fastq vs kmdb: $d rows"; fail=1; }
echo "== race: 11 targets, -t 1..64 x2; one .kmdb with its blocks split, -t 2..64"
"$W/ksd_ref" -k $W/k_dense.tsv -B $W/race.txt -G $W/bg.txt -t 1 -o $W/race_ref.gz 2>/dev/null; zcat $W/race_ref.gz > $W/race_ref.tsv
"$W/ksd_ref" -k $W/k_dense.tsv -B $W/one.txt -t 1 -o $W/one_ref.gz 2>/dev/null; zcat $W/one_ref.gz > $W/one_ref.tsv
n=0; d=0
for t in 1 2 3 4 8 16 32 64; do for r in 1 2; do
  "$W/ksd_new" -k $W/k_dense.tsv -B $W/race.txt -G $W/bg.txt -t $t -o $W/n.gz 2>/dev/null; n=$((n+1))
  zcat $W/n.gz | cmp -s - $W/race_ref.tsv || { d=$((d+1)); echo "  MISMATCH -t $t rep $r"; }
done
  "$W/ksd_new" -k $W/k_dense.tsv -B $W/one.txt -t $t -z 1 -o $W/n.gz 2>/dev/null; n=$((n+1))
  zcat $W/n.gz | cmp -s - $W/one_ref.tsv || { d=$((d+1)); echo "  MISMATCH single kmdb -t $t"; }
done; echo "  $n runs, $d mismatches"; [ $d -eq 0 ] || fail=1
echo "== ThreadSanitizer"
if gcc -std=gnu17 -O1 -g -fsanitize=thread -I"$src" -DHAVE_ZSTD -o "$W/ksd_tsan" "$src/kmer_strain_detect.c" $kl 2>/dev/null; then
  : > $W/tsan.log
  TSAN_OPTIONS="halt_on_error=0" "$W/ksd_tsan" -k $W/k_edge.txt -B $W/edge_targets.txt -t 8 -o $W/ts.gz 2>> $W/tsan.log
  TSAN_OPTIONS="halt_on_error=0" "$W/ksd_tsan" -k $W/k_dense.tsv -B $W/one.txt -t 8 -z 1 -o $W/ts.gz 2>> $W/tsan.log
  nr=$(grep -c "WARNING: ThreadSanitizer" $W/tsan.log); echo "  reports: $nr"; [ "$nr" -eq 0 ] || fail=1
else echo "  skipped (no -fsanitize=thread)"; fi
report
