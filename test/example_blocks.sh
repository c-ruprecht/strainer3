#!/usr/bin/env bash
# Example strainer run with read-aware (pair mode) detection straight from FASTQ:
#   - kmer_scrub_filter.py --reference places the rare k-mers on the genome (all copies), annotates
#     segments / blocks / repeat families, and with --context_bp also writes the single-copy reference
#     k-mers within 300 bp of a rare k-mer (context: used to check reads, never counted)
#   - fastp: duplicate removal only. No merging and no extra trimming; fastp still trims the adapters of
#     read pairs by overlap. The reads stay R1 / R2.
#   - kmer_block_detect -P -V -Q 20 -c scans the pairs directly (no kmdb):
#       -c  a read must match every rare AND context k-mer it covers on one diagonal
#       -P  the two mates are one fragment: overlapping mates are booked once, a break in one mate is
#           forgiven if the other mate reads the same reference k-mers cleanly, mates must face each other
#       -Q  a break next to a base below Q20 is a sequencing error, not a difference to the strain
#       -V  a fragment that contradicts the strain anywhere is dropped as a whole
#     It writes TWO k-mer tables from the same reads:
#       kmer_hits.tsv.gz       fragments whose reads passed the checks (read-aware)
#       kmer_hits_any.tsv.gz   every fragment that holds the k-mer (what a non-aware detector counts)
#   - block_coverage_depth.py turns both into hit-calling features: block_coverage_depth.tsv (read-aware)
#     and block_coverage_depth_any.tsv (no rejection), same rows and the same k-mer column names
# Run from strainer3/test after building the tools (cd ../src && make).
# By default only the first 100,000 read pairs of each sample are used, so the example finishes in a few
# minutes; MAX_PAIRS=0 bash example_blocks.sh runs all reads, FIGURES=png adds the coverage figures.
# RUN is the prefix for python / fastp; set RUN="" if they are already on the PATH.
set -e

OUT=output_blocks
STRAIN=D6331_Faecalibacterium_prausnitzii
REF=D6331_genomes/${STRAIN}.fasta.gz
THREADS=12
CONTEXT_BP=300       # reference context around the rare k-mers (about the insert size)
ERROR_QUALITY=20     # kmer_block_detect -Q
RUN_FASTP=1          # 0 = scan the raw read pairs
MAX_PAIRS=0 #${MAX_PAIRS-100000}   # read pairs per sample used by the example (0 = all 500,000); small = fast
FIGURES=${FIGURES-none}         # png draws the classic coverage figures (slow)
RUN=${RUN-conda run -n strainer3_env}

echo "Example strainer run for ${STRAIN} (read-aware detection from FASTQ pairs)"
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

echo 'STEP2: keeping kmers never seen outside the strain and placing them on the genome'
echo '       (scrub count 0; count 1 added if fewer than max(1000, 0.1 % of the strain kmers) are left after'
echo '        block filtering; all copies up to --max_ref_copies; runs >= 31 kmers cut into 31-61 kmer segments)'
echo "       plus the context kmers within ${CONTEXT_BP} bp of a rare kmer"
$RUN python ../scripts/kmer_scrub_filter.py \
-s "$OUT"/scrub/${STRAIN}.global_counts.tsv.zst \
-O "$OUT"/filter \
-n ${STRAIN} \
--selection count \
--min_kmers 1000 \
--min_kmers_frac 0.001 \
--reference ${REF} \
--max_ref_copies 5 \
--min_segment_kmers 31 \
--max_segment_kmers 31 \
--chunk_mode balanced \
--context_bp ${CONTEXT_BP}

zcat "$OUT"/filter/${STRAIN}.kmer_loci.tsv.gz | head -n 3
zcat "$OUT"/filter/${STRAIN}.blocks.tsv.gz | head -n 3
zcat "$OUT"/filter/${STRAIN}.context_loci.tsv.gz | head -n 3
cat "$OUT"/filter/${STRAIN}.selection.json
echo 'STEP2: complete'

printf '%*s' 20 '' | tr ' ' '='; printf ' STEP3 '; printf '%*s\n' 20 '' | tr ' ' '='

echo ''
echo 'STEP3: Read preprocessing with fastp (duplicate removal only) and sample sheet for the target metagenomes'
# raw reads: sample name <TAB> R1 <TAB> R2   (POSIX sh compatible: no arrays)
# A sample sequenced in several runs gets one line per run with the SAME name: kmer_block_detect pools them.
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
  if [ "$MAX_PAIRS" -gt 0 ]; then
    # the first MAX_PAIRS pairs of each file: enough to see every step work
    for M in 1 2; do
      SRC=$R1; [ $M -eq 2 ] && SRC=$R2
      gzip -dc "$SRC" | head -n $(( MAX_PAIRS * 4 )) | gzip -1 > "$OUT"/fastp/${NAME}.first.R${M}.fq.gz || true
    done
    R1="$OUT"/fastp/${NAME}.first.R1.fq.gz; R2="$OUT"/fastp/${NAME}.first.R2.fq.gz
  fi
  if [ "$RUN_FASTP" -eq 1 ]; then
    # --dedup only: PCR duplicates inflate counts and overdispersion.
    # Nothing is merged or quality-trimmed here; the detector deals with overlap (-P) and base quality (-Q).
    # fastp's defaults still apply: adapter trimming of pairs by overlap, and its basic quality / length filter.
    $RUN fastp \
      -i "$R1" -I "$R2" \
      -o "$OUT"/fastp/${NAME}.fastp.R1.fq.gz -O "$OUT"/fastp/${NAME}.fastp.R2.fq.gz \
      --dedup \
      --thread $(( THREADS < 16 ? THREADS : 16 )) \
      --json "$OUT"/fastp/${NAME}.fastp.json \
      --html "$OUT"/fastp/${NAME}.fastp.html \
      2> "$OUT"/fastp/${NAME}.fastp.log
    printf '%s\t%s\t%s\n' "$NAME" "$OUT/fastp/${NAME}.fastp.R1.fq.gz" "$OUT/fastp/${NAME}.fastp.R2.fq.gz" >> "$OUT"/samples.tsv
    $RUN python -c "
import json
j = json.load(open('$OUT/fastp/${NAME}.fastp.json'))
b, a = j['summary']['before_filtering'], j['summary']['after_filtering']
print('{}: reads {:,} -> {:,}, mean length {} -> {}, duplication {:.1%}, insert peak {}'.format(
    '${NAME}', b.get('total_reads', 0), a.get('total_reads', 0),
    b.get('read1_mean_length', 'NA'), a.get('read1_mean_length', 'NA'), j.get('duplication', {}).get('rate', 0),
    j.get('insert_size', {}).get('peak', 'NA')))
"
  else
    printf '%s\t%s\t%s\n' "$NAME" "$R1" "$R2" >> "$OUT"/samples.tsv
  fi
done 3< "$OUT"/raw_samples.tsv
cat "$OUT"/samples.tsv
echo 'STEP3 - DONE'

printf '%*s' 20 '' | tr ' ' '='; printf ' STEP4 '; printf '%*s\n' 20 '' | tr ' ' '='

echo ''
echo 'STEP4: Read-aware strain detection directly from the FASTQ pairs'
../src/kmer_block_detect \
-k "$OUT"/filter/${STRAIN}.kmer_loci.tsv.gz \
-c "$OUT"/filter/${STRAIN}.context_loci.tsv.gz \
-s "$OUT"/samples.tsv \
-O "$OUT"/strain_detect \
-n ${STRAIN} \
-t ${THREADS} \
-m 1 \
-P -V -Q ${ERROR_QUALITY}

echo 'kmer hits, read-aware (fragments per kmer and sample whose reads passed the checks; last row total_evaluated):'
zcat "$OUT"/strain_detect/${STRAIN}.kmer_hits.tsv.gz | head -n 4
zcat "$OUT"/strain_detect/${STRAIN}.kmer_hits.tsv.gz | tail -n 1
echo 'kmer hits, no rejection (every fragment that holds the kmer):'
zcat "$OUT"/strain_detect/${STRAIN}.kmer_hits_any.tsv.gz | head -n 4
echo 'block hits (frags = fragments placed on the block, starts = reads starting in it,'
echo '            refused = fragments that hit the block but contradict the strain):'
zcat "$OUT"/strain_detect/${STRAIN}.block_hits.tsv.gz | head -n 4
echo 'read / group / pair statistics:'
cat "$OUT"/strain_detect/${STRAIN}.detect_stats.tsv
echo 'STEP4 - DONE'

printf '%*s' 20 '' | tr ' ' '='; printf ' STEP5 '; printf '%*s\n' 20 '' | tr ' ' '='
echo 'STEP5: Block features for hit calling, read-aware and without rejection'
$RUN python ../scripts/block_coverage_depth.py \
--kmer_hits "$OUT"/strain_detect/${STRAIN}.kmer_hits.tsv.gz \
--kmer_hits_any "$OUT"/strain_detect/${STRAIN}.kmer_hits_any.tsv.gz \
--block_hits "$OUT"/strain_detect/${STRAIN}.block_hits.tsv.gz \
--detect_stats "$OUT"/strain_detect/${STRAIN}.detect_stats.tsv \
--blocks "$OUT"/filter/${STRAIN}.blocks.tsv.gz \
--kmer_loci "$OUT"/filter/${STRAIN}.kmer_loci.tsv.gz \
--reference ${REF} \
--strain ${STRAIN} \
--output_dir "$OUT"/block_coverage \
--export_hist

echo 'STEP5b: Classic coverage / depth for the current hit_call model, once per k-mer table'
# both tables keep the kmer_strain_detect layout, so coverage_depth.py runs unchanged on either
for KIND in validated any; do
  HITS="$OUT"/strain_detect/${STRAIN}.kmer_hits.tsv.gz
  [ "$KIND" = any ] && HITS="$OUT"/strain_detect/${STRAIN}.kmer_hits_any.tsv.gz
  $RUN python ../scripts/coverage_depth.py \
  --counts "$OUT"/scrub/${STRAIN}.global_counts.tsv.zst \
  --hits "$HITS" \
  --strain ${STRAIN} \
  --output_dir "$OUT"/coverage_${KIND} \
  --reference ${REF} \
  --export_hist \
  --figures "$FIGURES"
done
echo 'STEP5 - DONE'

printf '%*s' 20 '' | tr ' ' '='; printf ' STEP6 '; printf '%*s\n' 20 '' | tr ' ' '='
echo 'STEP6: Call hits on both tables (same model, same reads: the difference is the effect of read rejection)'
# no sylph profiles in this example -> only the 2-feature model (the 5-feature model needs --sylph_file)
# NOTE: the model was trained on unvalidated kmer counts; retrain before trusting calls on the read-aware table
for KIND in validated any; do
  $RUN python ../scripts/hit_call.py \
  --coverage_dir "$OUT"/coverage_${KIND} \
  --model ../data/xgb_hitcall_2param \
  --output_dir "$OUT"/hit_calls_${KIND}
done

echo 'read-aware against no rejection (k-mer coverage, zero inflation, relative abundance in % of sequenced bases):'
$RUN python -c "
import pandas as pd
d = pd.read_csv('$OUT/block_coverage/block_coverage_depth.tsv', sep='\t')
cols = ['strain', 'sample', 'any_kmer_all_coverage', 'kmer_all_coverage', 'rej_kmer_all_coverage',
        'any_kmer_all_zi_frac', 'kmer_all_zi_frac', 'rej_kmer_all_hits_frac', 'frags_hit', 'frags_contradicted',
        'pairs_overlap', 'groups_rescued', 'relative_abundance_any_kmer_all_naive',
        'relative_abundance_kmer_all_naive', 'relative_abundance_best_naive', 'detection_limit_rel_abund']
print(d[[c for c in cols if c in d.columns]].round(4).T.to_string(header=False))
for kind in ('validated', 'any'):
    h = pd.read_csv('$OUT/hit_calls_%s/hit_calls.tsv' % kind, sep='\t')
    print('hit calls, %s:' % kind)
    print(h[[c for c in h.columns if c in ('strain', 'sample') or c.endswith(('_prob', '_call'))]].to_string(index=False))
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
