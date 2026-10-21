# Strainer 3
## Change notes from strainer 2
- strainer3 is multithreaded in key steps to speed up the process of scrubbing and strain detection
- ships with a exhaustive pangenome and metagenome scrubbing database for human gut microbiome strain tracking. If used for tracking bacterial strains in other environments a custom metagenome scrub database likely needs to be built.
- target metagenomes are converted to bit packed kmer files before running strain detect to speed up (runs strains x target metagenomes)

## 

# How to set up
## intallation
git clone
cd src 
make

create strainer3_env
## Run test
- install strainer3_env first
cd test
bash example.sh

## Pulling scrub database

# Example workflows

# Functions
## Strainer3
### kmer_scrub_count
Output creates 2 files:
1. *.global_counts.tsv.zst
2. *.summary.tsv

If a -P flag is given a presence file will be created containing a list of scrub_ids from the summary file per kmer.
To reduce file sizes a presence list cap is recommended (-P 20), for a full presence list add -P 0,

###

## little helpers

### genome_compare

### genome_compare_presence
- inverted index