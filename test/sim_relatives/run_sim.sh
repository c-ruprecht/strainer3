#!/usr/bin/env bash
# usage: run_sim.sh <strain> <loci.tsv.gz> <context_loci.tsv.gz> <scrubbed_kmers.gz> <genome.fasta.gz> <outdir>
# Simulates the true strain and in silico relatives (diffuse SNPs / recombination tracts, ANI 99.95-98) in a
# short-insert and a long-insert library, runs fastp and compares non-aware counting with kmer_block_detect
# (default, -V, -V -c, each on single-end and pair-preserving input). Summarise with sim_table.py <outdir> <strain>.
# The fastp --dedup step needs about 6 GB; run one simulation at a time on a small machine.
set -euo pipefail
S=$1; LOCI=$2; CTX=$3; KM=$4; G=$5; O=$6
HERE=$(cd "$(dirname "$0")" && pwd); SRC=${SRC:-$HERE/../../src}; SC=$HERE
mkdir -p $O/reads $O/fastp
: > $O/samples_se.tsv; : > $O/samples_pp.tsv; printf 'sample_name\ttype\tfile1\n' > $O/targets_se.tsv
FP="--detect_adapter_for_pe --trim_poly_g --cut_front --cut_front_window_size 4 --cut_front_mean_quality 20 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 --qualified_quality_phred 20 --unqualified_percent_limit 10 --n_base_limit 2 --length_required 50 --low_complexity_filter --complexity_threshold 30 --merge --correction --overlap_len_require 30 --overlap_diff_limit 3 --overlap_diff_percent_limit 10 --dedup"
seed=100
mk () {  # name div model depth ins_mean ins_sd
  local N=$1; seed=$((seed+1))
  python3 -I $SC/simreads.py --genome $G --div $2 --model $3 --depth $4 --ins_mean $5 --ins_sd $6 --seed $seed --out $O/reads/$N 2>> $O/sim.log
  fastp -w 2 -i $O/reads/${N}_R1.fq.gz -I $O/reads/${N}_R2.fq.gz $FP --include_unmerged --merged_out $O/fastp/$N.se.fq.gz --json $O/fastp/$N.se.json --html /dev/null 2> /dev/null
  fastp -w 2 -i $O/reads/${N}_R1.fq.gz -I $O/reads/${N}_R2.fq.gz $FP --merged_out $O/fastp/$N.m.fq.gz --out1 $O/fastp/$N.R1.fq.gz --out2 $O/fastp/$N.R2.fq.gz --unpaired1 $O/fastp/$N.u.fq.gz --unpaired2 $O/fastp/$N.u.fq.gz --json $O/fastp/$N.pp.json --html /dev/null 2> /dev/null
  printf '%s\t%s\n' $N $O/fastp/$N.se.fq.gz >> $O/samples_se.tsv
  printf '%s\tSE\t%s\n' $N $O/fastp/$N.se.fq.gz >> $O/targets_se.tsv
  printf '%s\t%s\n%s\t%s\t%s\n%s\t%s\n' $N $O/fastp/$N.m.fq.gz $N $O/fastp/$N.R1.fq.gz $O/fastp/$N.R2.fq.gz $N $O/fastp/$N.u.fq.gz >> $O/samples_pp.tsv
  rm -f $O/reads/${N}_R1.fq.gz $O/reads/${N}_R2.fq.gz
}
for LIB in short:120:35 long:350:60; do L=${LIB%%:*}; R=${LIB#*:}; IM=${R%%:*}; IS=${R#*:}
  for D in 1 3; do mk true_d${D}_$L 0 diffuse $D $IM $IS; done
  for MODEL in diffuse tracts; do for DIV in 0.0005 0.001 0.003 0.005 0.01 0.02; do
    mk rel_${MODEL}_${DIV}_d3_$L $DIV $MODEL 3 $IM $IS
  done; done
done
run () { local D=$1 SH=$2; shift 2; mkdir -p $O/det_$D; $SRC/kmer_block_detect -k $LOCI -s $SH -O $O/det_$D -n $S -t 2 -m 1 "$@" 2> $O/det_$D/log; }
mkdir -p $O/det_nonaware; $SRC/kmer_strain_detect -k $KM -B $O/targets_se.tsv -O $O/det_nonaware -n $S -t 2 > $O/det_nonaware/log 2>&1
run block $O/samples_se.tsv
run veto $O/samples_se.tsv -V
run veto_ctx $O/samples_se.tsv -V -c $CTX
run pairs_veto $O/samples_pp.tsv -V
run pairs_veto_ctx $O/samples_pp.tsv -V -c $CTX
echo done $S
