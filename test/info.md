# select pangenome for 
dest=/metrica/codebase/strainer3/test/scrub_pangenome; mkdir -p "$dest"
awk -F'\t' 'NR>1 && $23 ~ /Faecalibacterium prausnitzii/ {print $15}' /metrica/databases/genomes/scrubdb-derep-hq-comp99-cont05/representatives_t96.tsv \
  | xargs -I{} cp -L {} "$dest"/
ls -1 "$dest"/*.fna.gz > "$dest/../genomes.list"; wc -l "$dest/../genomes.list"   

# small scrub metagenome db build in this case build with c 1 since it is already downsampled

/metrica/codebase/strainer3/src/kmer_metagenome_db -T scrub \
-A /metrica/codebase/strainer3/test/genomes_to_scrub.txt \
-i /metrica/databases/genomes/strainer_test_data/target_metagenomes/ATCC-20-staggerd_S22_L001_R1_001.fastq.downsampled.gz,/metrica/databases/genomes/strainer_test_data/target_metagenomes/ATCC-20-staggerd_S22_L001_R2_001.fastq.downsampled.gz \
-o /metrica/codebase/strainer3/test/scrub_metagenomes/ATCC-20-staggerd.kmdb \
-c 2

conda run -n strainer3_env python ../scripts/coverage_depth.py \
--counts output/scrub/D6331_Faecalibacterium_prausnitzii.global_counts.tsv.zst \
--hits output/strain_detect/D6331_Faecalibacterium_prausnitzii.kmer_hits.tsv.gz \
--strain D6331_Faecalibacterium_prausnitzii \
--output_dir output/coverage \
--reference D6331_genomes/D6331_Faecalibacterium_prausnitzii.fasta.gz