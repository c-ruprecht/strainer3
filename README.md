# Strainer 3
## Change notes from strainer 2
- strainer3 is multithreaded in key steps to speed up the process of scrubbing and strain detection
- ships with a exhaustive pangenome and metagenome scrubbing database for human gut microbiome strain tracking. If used for tracking bacterial strains in other environments a custom metagenome scrub database likely needs to be built.
- target metagenomes are converted to kmer files before running strain detect to speed up (runs strains x target metagenomes)


# Set up
## intallation
git clone https://github.com/c-ruprecht/strainer3.git
cd strainer3/src 
make

cd ../
conda env create -f strainer3_env.yml

## Run test
There are two tests in strainer3 based on commercially available in vitro cell mix sequencing data (downsampled).
The examples track strains in Zymo D6331 in a metagenome of D6331 and a metagenome of ATCC 10 even mix. There are 3 overlaping species between both mixes. The first test is a simple bash script of a tracking a single strain. Secondly, a more complex snakemake runs all target strains in D6331 in both samples.

| Species | ATCC 10 Strain Even Mix (10% each) | Zymo D6331 strain(s) |
|---|---|---|
| *Bifidobacterium adolescentis* | [ATCC 15703](https://www.atcc.org/products/15703) | LMG 10502 |
| *Enterococcus faecalis* | [ATCC 47077](https://www.atcc.org/products/47077) | IP101412 AER FAA 2 |
| *Escherichia coli* | [ATCC 700926](https://www.atcc.org/products/700926) | JM109, B-3008, B-2207, B-766, B-1109 |


### bash script
cd test
bash example.sh
### snakemake
conda activate strainer3_env
cd test
snakemake -s strainer3.smk --cores 12 --quiet rules

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