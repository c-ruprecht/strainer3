#!/usr/bin/env bash
# Example strainer run with block-aware detection straight from FASTQ:
#   - kmer_scrub_filter.py --reference places the rare k-mers on the genome (all copies)
#     and annotates segments / blocks / repeat families
#   - kmer_block_detect scans the FASTQ files directly (no kmdb). A read only counts for
#     a segment if it hits every segment k-mer inside its footprint on one diagonal, so
#     reads with mismatches (errors, relatives) are rejected and mates are counted once.
#   - fastp (strict): adapter / poly-G / quality trimming, overlap merge with error correction, dedup.
#     Merged and unmerged reads go into ONE file per sample (--include_unmerged), scanned as single-end:
#     reads with adapter read-through or bad tails would otherwise be rejected by the completeness check.
#   - block_coverage_depth.py turns the block / k-mer counts into hit-calling features
#     (Poisson on read starts per SNP block, robust zero inflation, relative abundance)
# Run from strainer3/test after building the tools (cd ../src && make).
set -e

OUT=output_blocks
STRAIN=D6331_Faecalibacterium_prausnitzii
REF=D6331_genomes/${STRAIN}.fasta.gz
THREADS=12
RUN_FASTP=1          # 0 = scan the raw read pairs (for comparing rejection rates)

echo "Example strainer run for ${STRAIN} (block-aware detection from FASTQ)"
echo 'STEP 0: Removing old test data if present'
rm -rf "$OUT"
mkdir -p "$OUT"

printf '%*s' 20 '' | tr ' ' '='; printf ' STEP1 '; printf '%*s\n' 20 '' | tr ' ' '='

echo 'STEP1: kmer_scrub_count (counting the frequency of target genome kmers in unrelated genomes and metagenomes to learn their frequency)'
../src/kmer_scrub_count \
-r ${REF} \
-A lists/genomes_to_scrub.txt \
-B lists/metagenomes_to_scrub.txt \
-C lists/drug_scrub.txt \
-t ${THREADS} \
-T 0.995 \
-O "$OUT"/scrub

zstdcat "$OUT"/scrub/${STRAIN}.global_counts.tsv.zst | head -n 5
echo 'STEP1: complete'

printf '%*s' 20 '' | tr ' ' '='; printf ' STEP2 '; printf '%*s\n' 20 '' | tr ' ' '='

echo 'STEP2: keeping kmers never seen outside the strain (scrub count 0; count 1 added if < 1000 kmers)'
echo '       and placing them on the genome'
echo '       (all copies up to --max_ref_copies; runs >= 31 kmers cut into equal 31-61 kmer segments)'
conda run -n strainer3_env python ../scripts/kmer_scrub_filter.py \
-s "$OUT"/scrub/${STRAIN}.global_counts.tsv.zst \
-O "$OUT"/filter \
-n ${STRAIN} \
--selection count \
--min_kmers 1000 \
--reference ${REF} \
--max_ref_copies 5 \
--min_segment_kmers 31 \
--max_segment_kmers 31 \
--chunk_mode balanced

zcat "$OUT"/filter/${STRAIN}.kmer_loci.tsv.gz | head -n 3
zcat "$OUT"/filter/${STRAIN}.blocks.tsv.gz | head -n 3
cat "$OUT"/filter/${STRAIN}.selection.json
echo 'STEP2: complete'

printf '%*s' 20 '' | tr ' ' '='; printf ' STEP3 '; printf '%*s\n' 20 '' | tr ' ' '='

echo ''
echo 'STEP3: Read preprocessing with fastp and sample sheet for the target metagenomes'
# raw reads: sample name <TAB> R1 <TAB> R2   (POSIX sh compatible: no arrays)
printf '%s\t%s\t%s\n' \
  ATCC-10-even target_metagenomes/ATCC-10-even_S21_L001_R1_001.fastq.downsampled.gz target_metagenomes/ATCC-10-even_S21_L001_R2_001.fastq.downsampled.gz \
  D6331        target_metagenomes/D6331-231101-C4_R1_001.fastq.downsampled.gz        target_metagenomes/D6331-231101-C4_R2_001.fastq.downsampled.gz \
  > "$OUT"/raw_samples.tsv
echo 'sample 1 atcc - negative'
echo 'sample 2 d6331 - positive'

mkdir -p "$OUT"/fastp
: > "$OUT"/samples.tsv
TAB=$(printf '\t')
# read the sheet on fd 3 so conda run / fastp cannot consume the loop's stdin
while IFS="$TAB" read -r NAME R1 R2 <&3; do
  [ -z "$NAME" ] && continue
  if [ "$RUN_FASTP" -eq 1 ]; then
    # strict settings for the block completeness check (a read must match every rare k-mer it overlaps):
    #   --detect_adapter_for_pe        adapter read-through is the main cause of rejected reads
    #   --trim_poly_g                  NovaSeq / NextSeq G-runs
    #   --cut_front / --cut_right      sliding-window trim (4 bp, mean Q20); cut_right drops the rest of the read
    #   -q 20 -u 10                    drop reads with > 10% bases below Q20
    #   -n 2                           at most 2 N per read (kmer_block_detect splits reads at N)
    #   -l 50                          short-insert libraries: 50 bp still carries 20 k-mers; the completeness check guards the rest
    #   -y -Y 30                       low-complexity filter
    #   --merge --correction           merge overlapping mates; mismatches in the overlap corrected by quality
    #   --overlap_len_require 30 --overlap_diff_limit 3 --overlap_diff_percent_limit 10   strict merge
    #   --include_unmerged             merged + unmerged reads in one file (single-end scan)
    #   --dedup                        PCR duplicates inflate counts and overdispersion
    conda run -n strainer3_env fastp \
      -i "$R1" -I "$R2" \
      --detect_adapter_for_pe \
      --trim_poly_g \
      --cut_front --cut_front_window_size 4 --cut_front_mean_quality 20 \
      --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 \
      --qualified_quality_phred 20 --unqualified_percent_limit 10 \
      --n_base_limit 2 \
      --length_required 50 \
      --low_complexity_filter --complexity_threshold 30 \
      --merge --correction \
      --overlap_len_require 30 --overlap_diff_limit 3 --overlap_diff_percent_limit 10 \
      --include_unmerged \
      --merged_out "$OUT"/fastp/${NAME}.fastp.fq.gz \
      --dedup \
      --thread $(( THREADS < 16 ? THREADS : 16 )) \
      --json "$OUT"/fastp/${NAME}.fastp.json \
      --html "$OUT"/fastp/${NAME}.fastp.html \
      2> "$OUT"/fastp/${NAME}.fastp.log
    printf '%s\t%s\n' "$NAME" "$OUT/fastp/${NAME}.fastp.fq.gz" >> "$OUT"/samples.tsv
    conda run -n strainer3_env python -c "
import json
j = json.load(open('$OUT/fastp/${NAME}.fastp.json'))
b, a = j['summary']['before_filtering'], j['summary']['after_filtering']
print('{}: reads {:,} -> {:,}, mean length {} -> {}, duplication {:.1%}'.format(
    '${NAME}', b.get('total_reads', 0), a.get('total_reads', 0),
    b.get('read1_mean_length', 'NA'), a.get('read1_mean_length', 'NA'), j.get('duplication', {}).get('rate', 0)))
"
  else
    printf '%s\t%s\t%s\n' "$NAME" "$R1" "$R2" >> "$OUT"/samples.tsv
  fi
done 3< "$OUT"/raw_samples.tsv
cat "$OUT"/samples.tsv
echo 'STEP3 - DONE'

printf '%*s' 20 '' | tr ' ' '='; printf ' STEP4 '; printf '%*s\n' 20 '' | tr ' ' '='

echo ''
echo 'STEP4: Block-aware strain detection directly from FASTQ'
echo '       (fastp output is scanned single-end: every merged / unmerged read is its own fragment)'
../src/kmer_block_detect \
-k "$OUT"/filter/${STRAIN}.kmer_loci.tsv.gz \
-s "$OUT"/samples.tsv \
-O "$OUT"/strain_detect \
-n ${STRAIN} \
-t ${THREADS} \
-m 1

echo 'kmer hits (fragments per kmer and sample, validated reads only; last row total_evaluated):'
zcat "$OUT"/strain_detect/${STRAIN}.kmer_hits.tsv.gz | head -n 4
zcat "$OUT"/strain_detect/${STRAIN}.kmer_hits.tsv.gz | tail -n 1
echo 'block hits (frags = fragments placed on the block, starts = reads starting in it):'
zcat "$OUT"/strain_detect/${STRAIN}.block_hits.tsv.gz | head -n 4
echo 'read / group statistics:'
cat "$OUT"/strain_detect/${STRAIN}.detect_stats.tsv
echo 'STEP4 - DONE'

printf '%*s' 20 '' | tr ' ' '='; printf ' STEP5 '; printf '%*s\n' 20 '' | tr ' ' '='
echo 'STEP5: Block features for hit calling (read starts per SNP block, zero inflation, relative abundance)'
conda run -n strainer3_env python ../scripts/block_coverage_depth.py \
--kmer_hits "$OUT"/strain_detect/${STRAIN}.kmer_hits.tsv.gz \
--block_hits "$OUT"/strain_detect/${STRAIN}.block_hits.tsv.gz \
--detect_stats "$OUT"/strain_detect/${STRAIN}.detect_stats.tsv \
--blocks "$OUT"/filter/${STRAIN}.blocks.tsv.gz \
--kmer_loci "$OUT"/filter/${STRAIN}.kmer_loci.tsv.gz \
--reference ${REF} \
--strain ${STRAIN} \
--output_dir "$OUT"/block_coverage \
--export_hist

echo 'STEP5b: Classic coverage / depth on the validated kmer hits (input for the current hit_call model)'
# kmer_hits.tsv.gz keeps the kmer_strain_detect layout, so coverage_depth.py runs unchanged
conda run -n strainer3_env python ../scripts/coverage_depth.py \
--counts "$OUT"/scrub/${STRAIN}.global_counts.tsv.zst \
--hits "$OUT"/strain_detect/${STRAIN}.kmer_hits.tsv.gz \
--strain ${STRAIN} \
--output_dir "$OUT"/coverage \
--reference ${REF} \
--export_hist \
--figures 'png'
echo 'STEP5 - DONE'

printf '%*s' 20 '' | tr ' ' '='; printf ' STEP6 '; printf '%*s\n' 20 '' | tr ' ' '='
echo 'STEP6: Call hits'
# no sylph profiles in this example -> only the 2-feature model (the 5-feature model needs --sylph_file)
# NOTE: the model was trained on unvalidated kmer counts; retrain on block_coverage_depth.tsv features
#       (strain_challenge_2 + truth table) before trusting calls
conda run -n strainer3_env python ../scripts/hit_call.py \
--coverage_dir "$OUT" \
--model ../data/xgb_hitcall_2param \
--output_dir "$OUT"

echo 'block-based summary (relative abundance in % of sequenced bases):'
conda run -n strainer3_env python -c "
import pandas as pd
d = pd.read_csv('$OUT/block_coverage/block_coverage_depth.tsv', sep='\t')
cols = ['strain', 'sample', 'blk_frac_hit', 'blk_starts_total', 'snp_starts_lambda', 'snp_starts_zi_frac',
        'snp_starts_zero_vs_robust', 'reject_frac', 'rate_ratio_long_vs_snp',
        'relative_abundance_starts_robust', 'relative_abundance_kmer_snp_robust', 'detection_limit_rel_abund']
print(d[[c for c in cols if c in d.columns]].round(4).T.to_string(header=False))
"

printf '%*s' 20 '' | tr ' ' '='; printf ' FINISHED '; printf '%*s\n' 20 '' | tr ' ' '='

cat <<'EOF'
     _             _                     _____
 ___| |_ _ __ __ _(_)_ __   ___ _ __    |___ /
/ __| __| '__/ _` | | '_ \ / _ \ '__|     |_ \
\__ \ |_| | | (_| | | | | |  __/ |       ___) |
|___/\__|_|  \__,_|_|_| |_|\___|_|      |____/
EOF
printf '%*s' 20 '' | tr ' ' '='; printf ' ENJOY '; printf '%*s\n' 20 '' | tr ' ' '='
