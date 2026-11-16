#!/usr/bin/env bash
# usage: run_sim_raw.sh <strain> <loci> <context_loci> <genome.fasta.gz> <outdir>   (same seeds as run_sim.sh)
# Pair mode on raw pairs: same simulated samples as run_sim.sh, but fastp --dedup only (no merging), then
# kmer_block_detect -P in four variants. Summarise with sim_table2.py <outdir> <strain> (needs run_sim.sh output too).
set -euo pipefail
S=$1; LOCI=$2; CTX=$3; G=$4; O=$5
SC=$(cd "$(dirname "$0")" && pwd); BIN=${BIN:-$SC/../../src/kmer_block_detect}
mkdir -p $O/reads $O/fastp; : > $O/samples_raw.tsv
seed=100
mk () { local N=$1; seed=$((seed+1))
  python3 -I $SC/simreads.py --genome $G --div $2 --model $3 --depth $4 --ins_mean $5 --ins_sd $6 --seed $seed --out $O/reads/$N 2>> $O/sim_raw.log
  fastp -w 2 -i $O/reads/${N}_R1.fq.gz -I $O/reads/${N}_R2.fq.gz --dedup -o $O/fastp/$N.dd.R1.fq.gz -O $O/fastp/$N.dd.R2.fq.gz --json $O/fastp/$N.dd.json --html /dev/null 2> /dev/null
  printf '%s\t%s\t%s\n' $N $O/fastp/$N.dd.R1.fq.gz $O/fastp/$N.dd.R2.fq.gz >> $O/samples_raw.tsv
  rm -f $O/reads/${N}_R1.fq.gz $O/reads/${N}_R2.fq.gz; }
for LIB in short:120:35 long:350:60; do L=${LIB%%:*}; R=${LIB#*:}; IM=${R%%:*}; IS=${R#*:}
  for D in 1 3; do mk true_d${D}_$L 0 diffuse $D $IM $IS; done
  for MODEL in diffuse tracts; do for DIV in 0.0005 0.001 0.003 0.005 0.01 0.02; do mk rel_${MODEL}_${DIV}_d3_$L $DIV $MODEL 3 $IM $IS; done; done
done
run () { local D=$1; shift; mkdir -p $O/det_$D; $BIN -k $LOCI -s $O/samples_raw.tsv -O $O/det_$D -n $S -t 2 -m 1 "$@" 2> $O/det_$D/log; }
run raw_indep -V -c $CTX
run raw_pair -V -P -c $CTX
run raw_pair_q20 -V -P -Q 20 -c $CTX
run raw_pair_noctx -V -P
echo done $S
