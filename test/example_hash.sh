#!/usr/bin/env bash
# Example strainer run with read hashes (kmer_metagenome_db / kmer_strain_detect --read-hash):
# same as example.sh, plus per k-mer min read hashes so coverage_depth.py can count the reads behind the hits.
# Run from strainer3/test after building the tools (cd ../src && make).
set -e

OUT=output_hash
STRAIN=D6331_Faecalibacterium_prausnitzii

echo "Example strainer run for ${STRAIN} (pair dedup + read hashes)"
echo 'STEP 0: Removing old test data if present'
rm -rf "$OUT"

printf '%*s' 20 '' | tr ' ' '='; printf ' STEP1 '; printf '%*s\n' 20 '' | tr ' ' '='

echo 'STEP1: kmer_scrub_count (counting the frequency of target genome kmers in unrelated genomes and metagenomes to learn their frequency)'
../src/kmer_scrub_count \
-r D6331_genomes/${STRAIN}.fasta.gz \
-A lists/genomes_to_scrub.txt \
-B lists/metagenomes_to_scrub.txt \
-C lists/drug_scrub.txt \
-t 12 \
-T 0.995 \
-O "$OUT"/scrub

zstdcat "$OUT"/scrub/${STRAIN}.global_counts.tsv.zst | head -n 5
echo 'STEP1: complete'

printf '%*s' 20 '' | tr ' ' '='; printf ' STEP2 '; printf '%*s\n' 20 '' | tr ' ' '='

echo 'STEP2: keeping the 1% most rare kmers in the target genome'
conda run -n strainer3_env python ../scripts/kmer_scrub_filter.py \
-s "$OUT"/scrub/${STRAIN}.global_counts.tsv.zst \
-O "$OUT"/filter \
-n ${STRAIN} \
-m 0.01 \
--keep_equal_rare \
--reference D6331_genomes/${STRAIN}.fasta.gz \
--min_segment_kmers 31
echo 'STEP2: complete'

printf '%*s' 20 '' | tr ' ' '='; printf ' STEP3 '; printf '%*s\n' 20 '' | tr ' ' '='

echo ''
echo 'STEP3: Create kmer db of target metagenomes with read hashes'
echo '       (pair dedup: mates joined with a comma; --read-hash: min read-pair hash per k-mer, kmdb version 2)'
echo 'sample 1 atcc - negative'
../src/kmer_metagenome_db -m sample --read-hash \
-O "$OUT"/target_samples_kmdb \
-n ATCC-10-even \
-i target_metagenomes/ATCC-10-even_S21_L001_R1_001.fastq.downsampled.gz,target_metagenomes/ATCC-10-even_S21_L001_R2_001.fastq.downsampled.gz

echo 'sample 2 d6331 - positive'
../src/kmer_metagenome_db -m sample --read-hash \
-O "$OUT"/target_samples_kmdb \
-n D6331 \
-i target_metagenomes/D6331-231101-C4_R1_001.fastq.downsampled.gz,target_metagenomes/D6331-231101-C4_R2_001.fastq.downsampled.gz

echo 'STEP3 - DONE'

printf '%*s' 20 '' | tr ' ' '='; printf ' STEP4 '; printf '%*s\n' 20 '' | tr ' ' '='

echo ''
echo 'STEP4: Detect strains in target metagenomes (with read hashes)'
echo 'Creating target list for strain detect'
for f in "$OUT"/target_samples_kmdb/*.kmdb; do echo "$f"; done > "$OUT"/targets.txt

../src/kmer_strain_detect --read-hash \
-k "$OUT"/filter/${STRAIN}.scrubbed_kmers.gz \
-B "$OUT"/targets.txt \
-O "$OUT"/strain_detect \
-n ${STRAIN}

echo 'kmer hits (count per sample):'
zcat "$OUT"/strain_detect/${STRAIN}.kmer_hits.tsv.gz | head -n 4
echo 'read hashes (min read-pair hash per kmer and sample, 0 = not seen; last row n_fragments):'
zcat "$OUT"/strain_detect/${STRAIN}.kmer_readhash.tsv.gz | head -n 4
zcat "$OUT"/strain_detect/${STRAIN}.kmer_readhash.tsv.gz | tail -n 1
echo 'STEP4 - DONE'

printf '%*s' 20 '' | tr ' ' '='; printf ' STEP5 '; printf '%*s\n' 20 '' | tr ' ' '='
echo 'STEP5: Calculating coverage and depth for rare kmers (+ read support from the read hashes)'
# the read-hash file next to --hits is picked up automatically (or give --readhash <file>)
conda run -n strainer3_env python ../scripts/coverage_depth.py \
--counts "$OUT"/scrub/${STRAIN}.global_counts.tsv.zst \
--hits "$OUT"/strain_detect/${STRAIN}.kmer_hits.tsv.gz \
--strain ${STRAIN} \
--output_dir "$OUT"/coverage \
--reference D6331_genomes/${STRAIN}.fasta.gz \
--export_hist \
--export_kmers \
--figures 'png'

echo 'read support per sample (reads behind the rare-kmer hits):'
conda run -n strainer3_env python -c "
import pandas as pd
d = pd.read_csv('$OUT/coverage/coverage_depth.tsv', sep='\t')
cols = ['sample', 'coverage_rare_kmers', 'rare_kmers_hit', 'rare_reads', 'rare_reads_eff',
        'rare_kmers_per_read', 'rare_top_read_share', 'rare_multi_frac', 'n_fragments']
print(d[[c for c in cols if c in d.columns]].to_string(index=False))
"
echo 'STEP5 - DONE'

printf '%*s' 20 '' | tr ' ' '='; printf ' STEP6 '; printf '%*s\n' 20 '' | tr ' ' '='
echo 'STEP6: Call hits'
# no sylph profiles in this example -> only the 2-feature model (the 5-feature model needs --sylph_file)
conda run -n strainer3_env python ../scripts/hit_call.py \
--coverage_dir "$OUT" \
--model ../data/xgb_hitcall_2param \
--output_dir "$OUT"

printf '%*s' 20 '' | tr ' ' '='; printf ' FINISHED '; printf '%*s\n' 20 '' | tr ' ' '='

cat <<'EOF'
     _             _                     _____
 ___| |_ _ __ __ _(_)_ __   ___ _ __    |___ /
/ __| __| '__/ _` | | '_ \ / _ \ '__|     |_ \
\__ \ |_| | | (_| | | | | |  __/ |       ___) |
|___/\__|_|  \__,_|_|_| |_|\___|_|      |____/
EOF
printf '%*s' 20 '' | tr ' ' '='; printf ' ENJOY '; printf '%*s\n' 20 '' | tr ' ' '='
