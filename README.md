# Strainer 3
## Change notes from strainer 2
- strainer3 is multithreaded in key steps to speed up the process of scrubbing and strain detection
- ships with a exhaustive pangenome and metagenome scrubbing database for human gut microbiome strain tracking. If used for tracking bacterial strains in other environments a custom metagenome scrub database likely needs to be built.
- target metagenomes are converted to kmer files before running strain detect to speed up (runs strains x target metagenomes)


# Set up
## intallation
````bash
git clone https://github.com/c-ruprecht/strainer3.git
cd strainer3/src 
make
cd ../
conda env create -f strainer3_env.yml
````


## Run tests
There are two tests in strainer3 based on commercially available in vitro cell mix sequencing data (downsampled).
The examples track strains in Zymo D6331 in a metagenome of D6331 and a metagenome of ATCC 10 even mix. There are 3 overlaping species between both mixes. The first test is a simple bash script of a tracking a single strain. Secondly, a more complex snakemake runs all target strains in D6331 in both samples.

| Species | ATCC 10 Strain Even Mix | Zymo D6331 |
|---|---|---|
| *Bifidobacterium adolescentis* | ATCC 15703 | LMG 10502 |
| *Enterococcus faecalis* | ATCC 47077| IP101412 AER FAA 2 |
| *Escherichia coli* | ATCC 700926 | JM109, B-3008, B-2207, B-766, B-1109 |


### bash script
````bash
cd test
bash example.sh
````

### snakemake
Output is created in test directory output-snake. Hit calls for all strains are in output-snake/hit_calls.tsv. The hit calls are based on the coverage_rare_kmers.  

````bash
conda activate strainer3_env
cd test
snakemake -s strainer3.smk --cores 12 --quiet rules
````

## Pulling scrub database

The scrub database contains ~160k dereplicated genomes of high quality (checkm2 completeness >99%, contamination <0.5%) and ~500 metagenomes of individual patients. A sourmash .rocksdb was created from the genomes which can be used to subset the genome database for your target strains. 
The metagenomes where built into kmdb files, scrubbing their kmers using the genome database. In each metagenome kmers are kept with count >= 2.
- full size: 
- genomes: 157296, 198.2 GiB
- metagenomes: 511, 153.8 GiB
- sourmash db: 1, 

# Understanding results:
## hit calls
- explain hit_call_pr
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