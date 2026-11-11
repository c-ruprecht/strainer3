import os
import sys
import json
import pandas as pd

# Small test of the block workflow (run from strainer3/test after `cd ../src && make`):
#   kmer_scrub_count -> kmer_scrub_filter.py --reference (count selection, balanced 31-kmer segments)
#   fastp (strict, merged + unmerged in one file) -> kmer_block_detect (FASTQ, all samples)
#   -> block_coverage_depth.py -> coverage_depth_all.tsv

OUTPUTDIR = 'output-snake-blocks'
DRUG = 'D6331'
STRAINDIR = f"{OUTPUTDIR}/{DRUG}"
DRUGLIST = 'lists/drug_scrub.txt'
PANSCRUB = 'lists/genomes_to_scrub.txt'
METASCRUB = 'lists/metagenomes_to_scrub.txt'
ENV = 'conda run -n strainer3_env'

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

# subset for quick tests:  --config n_strains=3   (first 3 of the list)
#                     or:  --config strains=D6331_Faecalibacterium_prausnitzii,D6331_Bacteroides_fragilis
if config.get('strains'):
    _want = [s.strip() for s in str(config['strains']).split(',') if s.strip()]
    _missing = [s for s in _want if s not in dict_strains]
    if _missing:
        raise ValueError(f"unknown strains {_missing}; available: {STRAINS}")
    STRAINS = _want
elif config.get('n_strains'):
    STRAINS = STRAINS[: int(config['n_strains'])]
print(f"[test] {len(STRAINS)} strain(s): {', '.join(STRAINS)}", file=sys.stderr)

wildcard_constraints:
    strain = '[^/]+',
    sample = '[^/]+',


onsuccess:
    rows = []
    for dirpath, _, files in os.walk(OUTPUTDIR + '/0-benchmarks'):
        for f in files:
            if f.endswith('.tsv'):
                df = pd.read_csv(os.path.join(dirpath, f), sep='\t')
                df.insert(0, 'job', f[:-4])
                df.insert(0, 'rule', os.path.basename(dirpath))
                rows.append(df)
    if rows:
        pd.concat(rows, ignore_index=True).to_csv(OUTPUTDIR + '/benchmarks_all.tsv', sep='\t', index=False)
        print(f"[benchmarks] {len(rows)} jobs -> {OUTPUTDIR}/benchmarks_all.tsv", file=sys.stderr)


rule all:
    input:
        OUTPUTDIR + '/coverage_depth_all.tsv',
        OUTPUTDIR + '/block_hist_all.tsv.gz',
        OUTPUTDIR + '/blocks_long_all.tsv.gz'


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


# rare k-mers placed on the genome: scrub count 0 (count 1 added if < 1000 k-mers), all copies up to 5,
# runs >= 31 k-mers cut into equal 31..61 k-mer segments
rule kmer_filter:
    input:
        counts = STRAINDIR + '/{strain}/scrub/{strain}.global_counts.tsv.zst',
        ref = lambda wc: dict_strains[wc.strain]
    output:
        loci = STRAINDIR + '/{strain}/filter/{strain}.kmer_loci.tsv.gz',
        blocks = STRAINDIR + '/{strain}/filter/{strain}.blocks.tsv.gz',
        kmers = STRAINDIR + '/{strain}/filter/{strain}.scrubbed_kmers.gz',
        selection = STRAINDIR + '/{strain}/filter/{strain}.selection.json'
    params:
        outdir = lambda wc: f"{STRAINDIR}/{wc.strain}/filter"
    benchmark:
        OUTPUTDIR + '/0-benchmarks/kmer_filter/{strain}.tsv'
    threads: 1
    shell:
        """
        {ENV} python ../scripts/kmer_scrub_filter.py \
        -s {input.counts} \
        -O {params.outdir} \
        -n {wildcards.strain} \
        --selection count \
        --min_kmers 1000 \
        --reference {input.ref} \
        --max_ref_copies 5 \
        --min_segment_kmers 31 \
        --max_segment_kmers 31 \
        --chunk_mode balanced
        """


# shared by all strains: strict fastp per target, merged + unmerged reads in one file (scanned single-end)
rule fastp:
    input:
        r1 = lambda wc: dict_targets[wc.sample][0],
        r2 = lambda wc: dict_targets[wc.sample][1]
    output:
        reads = OUTPUTDIR + '/fastp/{sample}.fastp.fq.gz',
        json = OUTPUTDIR + '/fastp/{sample}.fastp.json',
        html = OUTPUTDIR + '/fastp/{sample}.fastp.html'
    log:
        OUTPUTDIR + '/fastp/{sample}.fastp.log'
    benchmark:
        OUTPUTDIR + '/0-benchmarks/fastp/{sample}.tsv'
    threads: 4
    shell:
        '''
        {ENV} fastp -w {threads} \
        -i {input.r1} -I {input.r2} \
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
        --merged_out {output.reads} \
        --dedup \
        --json {output.json} --html {output.html} > {log} 2>&1
        '''


rule sample_sheet:
    input:
        expand(OUTPUTDIR + '/fastp/{sample}.fastp.fq.gz', sample=SAMPLES)
    output:
        OUTPUTDIR + '/samples.tsv'
    run:
        with open(output[0], 'w') as fh:
            for smp, path in zip(SAMPLES, input):
                fh.write(f"{smp}\t{os.path.abspath(path)}\n")


rule block_detect:
    input:
        loci = STRAINDIR + '/{strain}/filter/{strain}.kmer_loci.tsv.gz',
        samples = OUTPUTDIR + '/samples.tsv',
        reads = expand(OUTPUTDIR + '/fastp/{sample}.fastp.fq.gz', sample=SAMPLES)
    output:
        kmer_hits = STRAINDIR + '/{strain}/block_detect/{strain}.kmer_hits.tsv.gz',
        block_hits = STRAINDIR + '/{strain}/block_detect/{strain}.block_hits.tsv.gz',
        family = STRAINDIR + '/{strain}/block_detect/{strain}.family_hits.tsv.gz',
        stats = STRAINDIR + '/{strain}/block_detect/{strain}.detect_stats.tsv',
        lengths = STRAINDIR + '/{strain}/block_detect/{strain}.read_lengths.tsv.gz'
    params:
        outdir = lambda wc: f"{STRAINDIR}/{wc.strain}/block_detect"
    benchmark:
        OUTPUTDIR + '/0-benchmarks/block_detect/{strain}.tsv'
    threads: 2
    shell:
        '''
        ../src/kmer_block_detect \
        -k {input.loci} \
        -s {input.samples} \
        -O {params.outdir} \
        -n {wildcards.strain} \
        -t {threads} \
        -m 1
        '''


rule block_coverage_depth:
    input:
        kmer_hits = STRAINDIR + '/{strain}/block_detect/{strain}.kmer_hits.tsv.gz',
        block_hits = STRAINDIR + '/{strain}/block_detect/{strain}.block_hits.tsv.gz',
        stats = STRAINDIR + '/{strain}/block_detect/{strain}.detect_stats.tsv',
        lengths = STRAINDIR + '/{strain}/block_detect/{strain}.read_lengths.tsv.gz',
        family = STRAINDIR + '/{strain}/block_detect/{strain}.family_hits.tsv.gz',
        blocks = STRAINDIR + '/{strain}/filter/{strain}.blocks.tsv.gz',
        loci = STRAINDIR + '/{strain}/filter/{strain}.kmer_loci.tsv.gz',
        selection = STRAINDIR + '/{strain}/filter/{strain}.selection.json',
        ref = lambda wc: dict_strains[wc.strain]
    output:
        cov = STRAINDIR + '/{strain}/block_coverage/block_coverage_depth.tsv',
        hist = STRAINDIR + '/{strain}/block_coverage/{strain}.block_hist.tsv.gz',
        blocks_long = STRAINDIR + '/{strain}/block_coverage/{strain}.blocks_long.tsv.gz'
    params:
        outdir = lambda wc: f"{STRAINDIR}/{wc.strain}/block_coverage"
    benchmark:
        OUTPUTDIR + '/0-benchmarks/block_coverage_depth/{strain}.tsv'
    threads: 1
    shell:
        """
        {ENV} python ../scripts/block_coverage_depth.py \
        --kmer_hits {input.kmer_hits} \
        --block_hits {input.block_hits} \
        --detect_stats {input.stats} \
        --read_lengths {input.lengths} \
        --family_hits {input.family} \
        --blocks {input.blocks} \
        --kmer_loci {input.loci} \
        --selection {input.selection} \
        --reference {input.ref} \
        --strain {wildcards.strain} \
        --output_dir {params.outdir} \
        --export_hist \
        --export_blocks \
        --figures png
        """


rule coverage_depth_concat:
    input:
        cov = expand(STRAINDIR + '/{strain}/block_coverage/block_coverage_depth.tsv', strain=STRAINS),
        hist = expand(STRAINDIR + '/{strain}/block_coverage/{strain}.block_hist.tsv.gz', strain=STRAINS),
        blocks_long = expand(STRAINDIR + '/{strain}/block_coverage/{strain}.blocks_long.tsv.gz', strain=STRAINS),
        fastp = expand(OUTPUTDIR + '/fastp/{sample}.fastp.json', sample=SAMPLES)
    output:
        cov = OUTPUTDIR + '/coverage_depth_all.tsv',
        hist = OUTPUTDIR + '/block_hist_all.tsv.gz',
        blocks_long = OUTPUTDIR + '/blocks_long_all.tsv.gz'
    run:
        pd.concat([pd.read_csv(p, sep='\t') for p in input.hist], ignore_index=True) \
          .to_csv(output.hist, sep='\t', index=False)
        pd.concat([pd.read_csv(p, sep='\t').assign(strain=s_) for s_, p in zip(STRAINS, input.blocks_long)],
                  ignore_index=True).to_csv(output.blocks_long, sep='\t', index=False)
        tab = pd.concat([pd.read_csv(p, sep='\t') for p in input.cov], ignore_index=True)
        rows = []
        for smp, p in zip(SAMPLES, input.fastp):
            j = json.load(open(p))
            b, a = j['summary']['before_filtering'], j['summary']['after_filtering']
            rows.append(dict(sample=smp, fastp_reads_in=b.get('total_reads'), fastp_reads_out=a.get('total_reads'),
                             fastp_too_short=j.get('filtering_result', {}).get('too_short_reads'),
                             fastp_insert_peak=j.get('insert_size', {}).get('peak'),
                             fastp_duplication=j.get('duplication', {}).get('rate')))
        tab['sample'] = tab['sample'].astype(str)
        tab = tab.merge(pd.DataFrame(rows), on='sample', how='left')
        tab.to_csv(output.cov, sep='\t', index=False)
        cols = ['strain', 'sample', 'best_total', 'best_frac_hit', 'best_zip_pi', 'best_zero_run_p',
                'relative_abundance_best_zip_present', 'detection_limit_rel_abund', 'reject_frac']
        print(tab[[c for c in cols if c in tab.columns]].round(4).to_string(index=False))