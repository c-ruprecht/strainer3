import os
import pandas as pd

OUTPUTDIR = 'output-snake'
DRUG = 'D6331'
STRAINDIR = f"{OUTPUTDIR}/{DRUG}"
DRUGLIST = 'lists/drug_scrub.txt'
PANSCRUB = 'lists/genomes_to_scrub.txt'
METASCRUB = 'lists/metagenomes_to_scrub.txt'

# target metagenomes: sample name -> (R1, R2)
dict_targets = {'ATCC-10-even': ('target_metagenomes/ATCC-10-even_S21_L001_R1_001.fastq.downsampled.gz',
                                 'target_metagenomes/ATCC-10-even_S21_L001_R2_001.fastq.downsampled.gz'),
                'D6331': ('target_metagenomes/D6331-231101-C4_R1_001.fastq.downsampled.gz',
                          'target_metagenomes/D6331-231101-C4_R2_001.fastq.downsampled.gz'),
                }
SAMPLES = list(dict_targets)

df_target_strains = pd.read_csv(DRUGLIST, header=None, names=['strain_path'])

df_target_strains['strain'] = (df_target_strains['strain_path']
                               .str.replace(r'^.*/', '', regex=True)
                               .str.replace(r'\.fasta\.gz$', '', regex=True))
dict_strains = df_target_strains.set_index('strain')['strain_path'].to_dict()
STRAINS = list(dict_strains)


rule all:
    input:
        OUTPUTDIR + '/hit_calls.tsv'


rule kmer_scrub:
    input:
        strain_path = lambda wc: dict_strains[wc.strain],
        gscrub = PANSCRUB,
        mscrub = METASCRUB,
        dscrub = DRUGLIST
    output:
        STRAINDIR + '/{strain}/scrub/{strain}.global_counts.tsv.zst'
    params:
        outdir = lambda wc: f"{STRAINDIR}/{wc.strain}/scrub"
    threads: 4
    shell:
        '''
        ../src/kmer_scrub_count \
        -r {input.strain_path} \
        -A {input.gscrub} \
        -B {input.mscrub} \
        -C {input.dscrub} \
        -t {threads} \
        -T 0.995 \
        -O {params.outdir}
        '''


rule kmer_filter:
    input:
        STRAINDIR + '/{strain}/scrub/{strain}.global_counts.tsv.zst'
    output:
        STRAINDIR + '/{strain}/filter/{strain}.scrubbed_kmers.gz'
    params:
        outdir = lambda wc: f"{STRAINDIR}/{wc.strain}/filter"
    threads: 1
    shell:
        """
        conda run -n strainer3_env python ../scripts/kmer_scrub_filter.py \
        -s {input} \
        -O {params.outdir} \
        -n {wildcards.strain} \
        -m 0.01
        """


# shared by all strains: one kmdb per target metagenome
rule target_kmdb:
    input:
        r1 = lambda wc: dict_targets[wc.sample][0],
        r2 = lambda wc: dict_targets[wc.sample][1]
    output:
        OUTPUTDIR + '/target_samples_kmdb/{sample}.kmdb'
    params:
        outdir = OUTPUTDIR + '/target_samples_kmdb'
    threads: 4
    shell:
        '''
        ../src/kmer_metagenome_db -m sample \
        -O {params.outdir} \
        -n {wildcards.sample} \
        {input.r1} {input.r2}
        '''


rule target_list:
    input:
        expand(OUTPUTDIR + '/target_samples_kmdb/{sample}.kmdb', sample=SAMPLES)
    output:
        OUTPUTDIR + '/targets.txt'
    shell:
        '''
        printf '%s\\n' {input} > {output}
        '''


rule strain_detect:
    input:
        kmers = STRAINDIR + '/{strain}/filter/{strain}.scrubbed_kmers.gz',
        targets = OUTPUTDIR + '/targets.txt',
        kmdbs = expand(OUTPUTDIR + '/target_samples_kmdb/{sample}.kmdb', sample=SAMPLES)
    output:
        STRAINDIR + '/{strain}/strain_detect/{strain}.kmer_hits.tsv.gz'
    params:
        outdir = lambda wc: f"{STRAINDIR}/{wc.strain}/strain_detect"
    threads: 4
    shell:
        '''
        ../src/kmer_strain_detect \
        -k {input.kmers} \
        -B {input.targets} \
        -O {params.outdir} \
        -n {wildcards.strain} \
        -t {threads}
        '''


rule coverage_depth:
    input:
        counts = STRAINDIR + '/{strain}/scrub/{strain}.global_counts.tsv.zst',
        hits = STRAINDIR + '/{strain}/strain_detect/{strain}.kmer_hits.tsv.gz'
    output:
        STRAINDIR + '/{strain}/coverage/coverage_depth.tsv'
    params:
        outdir = lambda wc: f"{STRAINDIR}/{wc.strain}/coverage"
    threads: 1
    shell:
        """
        conda run -n strainer3_env python ../scripts/coverage_depth.py \
        --counts {input.counts} \
        --hits {input.hits} \
        --strain {wildcards.strain} \
        --output_dir {params.outdir}
        """


rule hit_call:
    input:
        expand(STRAINDIR + '/{strain}/coverage/coverage_depth.tsv', strain=STRAINS)
    output:
        OUTPUTDIR + '/hit_calls.tsv'
    params:
        covdir = STRAINDIR,
        outdir = OUTPUTDIR
    threads: 1
    shell:
        """
        conda run -n strainer3_env python ../scripts/hit_call.py \
        --coverage_dir {params.covdir} \
        --output_dir {params.outdir}
        """