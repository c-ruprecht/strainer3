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
- full size: ~400 GiB
- genomes: 157296
- metagenomes: 511
- sourmash db: 1

To download the full database, unpack and clean up the download directory. During the building process you require roughly 2x the full database size (~800 GiB) of free disk space.

````bash
conda activate strainer3_env
python scripts/scrubdb_hf.py pull \
  --repo         cruprecht/strainer3-scrub-db \
  --revision     v2026.10 \
  --download-dir /path/to/download/dir \
  --dest         /path/to/final/db \
  --cleanup
````
# Running strainer3:
- there is an example workflow that runs a full snakemake for all genomes in a set of genome lists
## Preparation
### dereplicate strains you want to track and create one or more genome list files
- you can use tools like skani, or also genome_compare_presence to dereplicate strains
- for incomplete genomes with high contamination skani is recommended 
- a genome list is a simple .txt file with one genome name per file

### Create sample sheets
- multiple samples can be merged into one target kmdb
- simple tsv separated file with columns sample_name and seq_files. seq_files should contain a comma separated list

### Modify config files:
#### Strainer specific config:
- look at example_workflows/config_strainer3.yaml
- Need to modify:
  - scrub_db_path -> folder where the scrub database was built
  - output directory -> where files will be saved
  - target_metagenomes -> path to the sample sheet file
  - strainer_path -> git repo of strainer3
  - genome_lists: one list entry per genome list to track
  - locations: parent folder of where genomes and metagenomes can be found

#### Config for cluster/local run
- snakemake requires a general configfile to bennamed config.yaml and it is called by its parent directory.
- example_workflows/profiles/lsf/config.yaml
-   adapt to lsf settings on your cluster, project_name, log locations etc.
- example_workflows/profiles/env_lsf.yaml you can modify rule specific allocations e.g.: 
  - kmer_scrub: {threads: 45, mem_mb: 32000,  runtime: 720, lsf_queue: premium}

## Run strainer
````bash
conda activate strainer3_env
snakemake -s example_workflows/strainer3_containment \
--configfiles example_workflows/config_strainer3.yaml example_workflows/profiles/env_lsf.yaml \
--profile example_workflows/profiles/lsf
````


# Understanding results:
## output folder structure
## Output folder
| Wildcard | Meaning |
|---|---|
| `{target}` | target metagenome, from `sample_name` in the sample sheet |
| `{glist}` | genome list, the file name of the list without extension (`drug_scrub.txt` becomes `drug_scrub`) |
| `{strain}` | strain genome, the file name without sequence extension (`X.fna.gz` becomes `X`) |

A strain that appears in several genome lists is processed once per list, because each list is used as the drug scrub for its own strains.

```
{output_dir}/
├── hit_calls.tsv                                   # hit calls for all genome lists combined
├── targets.txt                                     # list of target metagenome databases
├── target_samples_kmdb/
│   └── {target}.kmdb                               # k-mer database per target metagenome
├── drug_lists/
│   └── {glist}.drug_scrub.txt                      # genomes of each genome list, used as drug scrub
├── sourmash/                                       # only if containment is not "all"
│   ├── sketches/{strain}.sig
│   ├── all_vs_db.csv                               # strains vs scrub database
│   ├── hits/{strain}.csv
│   └── {strain}.genomes_to_scrub.txt               # scrub genomes passing the containment cutoff
├── genome_lists/
│   └── {glist}/
│       ├── hit_calls.tsv                           # hit calls for this genome list
│       └── {strain}/
│           ├── scrub/{strain}.global_counts.tsv.zst
│           ├── filter/{strain}.scrubbed_kmers.gz
│           ├── strain_detect/{strain}.kmer_hits.tsv.gz
│           └── coverage/coverage_depth.tsv
└── 0-logs/                                         # one log per rule and job
```

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