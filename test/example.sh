
echo 'Example strainer run for D6331_Faecalibacterium_prausnitzii'
echo 'STEP 0: Removing old test data if present'
rm -rf output

printf '%*s' 20 '' | tr ' ' '='; printf ' STEP1 '; printf '%*s\n' 20 '' | tr ' ' '='

echo 'STEP1: kmer_scrub_count (counting the frequency of target genome kmers in unrelated genomes and metagenomes to learn their frequency'
../src/kmer_scrub_count \
-r D6331_genomes/D6331_Faecalibacterium_prausnitzii.fasta.gz \
-A lists/genomes_to_scrub.txt \
-B lists/metagenomes_to_scrub.txt \
-C lists/drug_scrub.txt \
-t 4 \
-T 0.995 \
-O output/scrub 

zstdcat output/scrub/D6331_Faecalibacterium_prausnitzii.global_counts.tsv.zst | head -n 5
echo 'STEP1: complete'

printf '%*s' 20 '' | tr ' ' '='; printf ' STEP2 '; printf '%*s\n' 20 '' | tr ' ' '='

echo 'STEP2: keeping the 1% most rare kmers in the target genome'
conda run -n strainer3_env python ../scripts/kmer_scrub_filter.py \
-s output/scrub/D6331_Faecalibacterium_prausnitzii.global_counts.tsv.zst \
-O output/filter \
-n D6331_Faecalibacterium_prausnitzii
-m 0.01 
echo 'STEP2: complete'

printf '%*s' 20 '' | tr ' ' '='; printf ' STEP3 '; printf '%*s\n' 20 '' | tr ' ' '='

echo ''
echo 'STEP3: Create kmer db of target metagenomes'
echo 'sample 1 atcc - negative'
../src/kmer_metagenome_db -m sample \
-O output/target_samples_kmdb \
-n ATCC-10-even \
target_metagenomes/ATCC-10-even_S21_L001_R1_001.fastq.downsampled.gz target_metagenomes/ATCC-10-even_S21_L001_R2_001.fastq.downsampled.gz

echo 'sample 2 d6331 - positive'
../src/kmer_metagenome_db -m sample \
-O output/target_samples_kmdb \
-n D6331 \
target_metagenomes/D6331-231101-C4_R1_001.fastq.downsampled.gz target_metagenomes/D6331-231101-C4_R2_001.fastq.downsampled.gz

echo 'STEP3 - DONE'

printf '%*s' 20 '' | tr ' ' '='; printf ' STEP4 '; printf '%*s\n' 20 '' | tr ' ' '='

echo ''
echo 'STEP4: Detect strains in target metagenomes'
echo 'Creating target list for strain detect'
for f in output/target_samples_kmdb/*.kmdb; do echo "$f"; done > output/targets.txt

../src/kmer_strain_detect \
-k output/filter/D6331_Faecalibacterium_prausnitzii.scrubbed_kmers.gz \
-B output/targets.txt \
-O output/strain_detect \
-n D6331_Faecalibacterium_prausnitzii


printf '%*s' 20 '' | tr ' ' '='; printf ' STEP5 '; printf '%*s\n' 20 '' | tr ' ' '='
echo 'STEP5: Calculating coverage and depth for rare kmers'
conda run -n strainer3_env python ../scripts/coverage_depth.py \
--counts output/scrub/D6331_Faecalibacterium_prausnitzii.global_counts.tsv.zst \
--hits output/strain_detect/D6331_Faecalibacterium_prausnitzii.kmer_hits.tsv.gz \
--strain D6331_Faecalibacterium_prausnitzii \
--output_dir output/coverage

printf '%*s' 20 '' | tr ' ' '='; printf ' STEP6 '; printf '%*s\n' 20 '' | tr ' ' '='
echo 'STEP6: Call hits'
conda run -n strainer3_env python ../scripts/hit_call.py \
--coverage_dir output \
--output_dir output 

printf '%*s' 20 '' | tr ' ' '='; printf ' FINISHED '; printf '%*s\n' 20 '' | tr ' ' '='

cat <<'EOF'
     _             _                     _____
 ___| |_ _ __ __ _(_)_ __   ___ _ __    |___ /
/ __| __| '__/ _` | | '_ \ / _ \ '__|     |_ \
\__ \ |_| | | (_| | | | | |  __/ |       ___) |
|___/\__|_|  \__,_|_|_| |_|\___|_|      |____/
EOF
printf '%*s' 20 '' | tr ' ' '='; printf ' ENJOY '; printf '%*s\n' 20 '' | tr ' ' '='
