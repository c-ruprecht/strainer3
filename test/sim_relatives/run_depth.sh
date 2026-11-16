#!/usr/bin/env bash
# depth series: true strain and relatives at 1-30x, two insert sizes; non-aware counting vs pair mode
set -euo pipefail
S=$1; LOCI=$2; CTX=$3; KM=$4; G=$5; O=$6
SC=$(cd "$(dirname "$0")" && pwd); SRC=${SRC:-$SC/../../src}; BIN=${BIN:-$SRC/kmer_block_detect}
mkdir -p $O/reads $O/fastp; : > $O/samples.tsv; printf 'sample_name\ttype\tfile1\n' > $O/targets.tsv
seed=900
for LIB in long:350:60 short:120:35; do L=${LIB%%:*}; R=${LIB#*:}; IM=${R%%:*}; IS=${R#*:}
 for X in true:0:diffuse d999:0.001:diffuse d995:0.005:diffuse t999:0.001:tracts t995:0.005:tracts t990:0.01:tracts; do K=${X%%:*}; R2=${X#*:}; DIV=${R2%%:*}; MODEL=${R2#*:}
  for D in 1 3 10 30; do N=${K}_d${D}_$L; seed=$((seed+1))
   python3 -I $SC/simreads.py --genome $G --div $DIV --model $MODEL --depth $D --ins_mean $IM --ins_sd $IS --seed $seed --out $O/reads/$N 2>> $O/sim.log
   fastp -w 2 -i $O/reads/${N}_R1.fq.gz -I $O/reads/${N}_R2.fq.gz --dedup -o $O/fastp/$N.R1.fq.gz -O $O/fastp/$N.R2.fq.gz --json $O/fastp/$N.json --html /dev/null 2>/dev/null
   cat $O/fastp/$N.R1.fq.gz $O/fastp/$N.R2.fq.gz > $O/fastp/$N.both.fq.gz
   printf '%s\t%s\t%s\n' $N $O/fastp/$N.R1.fq.gz $O/fastp/$N.R2.fq.gz >> $O/samples.tsv
   printf '%s\tSE\t%s\n' $N $O/fastp/$N.both.fq.gz >> $O/targets.tsv
   rm -f $O/reads/${N}_R?.fq.gz
  done; done; done
mkdir -p $O/det_nonaware $O/det_pair $O/det_block
$SRC/kmer_strain_detect -k $KM -B $O/targets.tsv -O $O/det_nonaware -n $S -t 2 > $O/det_nonaware/log 2>&1
$BIN -k $LOCI -s $O/samples.tsv -O $O/det_block -n $S -t 2 2> $O/det_block/log
$BIN -k $LOCI -s $O/samples.tsv -O $O/det_pair -n $S -t 2 -P -V -Q 20 -c $CTX 2> $O/det_pair/log
echo done
