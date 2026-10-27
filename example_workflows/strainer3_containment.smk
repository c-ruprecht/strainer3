import os
import re
import sys
import pandas as pd

# =============================================================================
# strainer3 workflow
#
#   [containment != "all"]  sourmash sketch (per strain) -> manysearch vs scrub DB (once)
#                           -> per strain: genomes_to_scrub list (DB genomes with containment >= cutoff)
#   [containment == "all"]  skips sourmash, scrubs against the whole genomes_to_scrub_full.txt
#   -> kmer_scrub_count -> kmer_scrub_filter.py -> kmer_strain_detect
#   -> coverage_depth.py (reference, terminal distance, rare kmer count) -> hit_call.py per genome list
#   -> hit_calls.tsv (all genome lists, with containment / genome_list columns)
#
# Machine-specific settings (module loading, conda call, threads / memory / queue per rule)
# come from the config, so the same Snakefile runs on the cluster and locally:
#   shell_prefix: "ml anaconda3/latest &&"                  # "" locally
#   env_run:      "conda run --no-capture-output -n strainer3_env"   # "" if already in the env
#   resources:                                              # optional per-rule overrides
#     kmer_scrub: {threads: 28, mem_mb: 16000, runtime: 720, lsf_queue: premium}
#
# Output layout
#   {output_dir}/targets.txt
#   {output_dir}/hit_calls.tsv
#   {output_dir}/target_samples_kmdb/{target}.kmdb
#   {output_dir}/drug_lists/{glist}.drug_scrub.txt
#   {output_dir}/sourmash/{sketches/{strain}.sig, all_vs_db.csv, hits/{strain}.csv, {strain}.genomes_to_scrub.txt}
#   {output_dir}/{glist}/hit_calls.tsv
#   {output_dir}/{glist}/{strain}/{scrub,filter,strain_detect,coverage}/
#   {output_dir}/0-logs/...
# =============================================================================

strainer_dir = config["strainer_path"]
src_dir      = config.get("src_dir", os.path.join(strainer_dir, "src"))
scripts_dir  = config.get("scripts_dir", os.path.join(strainer_dir, "scripts"))
output_dir   = config["output_dir"]
scrub_db     = config["scrub_db_path"]

# scrub DB as restored by scrubdb_hf.py
ROCKSDB     = os.path.join(scrub_db, config.get("rocksdb", "high_qual-comp99-cont05-maxreps2000.rocksdb"))
GENOME_LIST = config.get("genome_scrub_list", os.path.join(scrub_db, "genomes_to_scrub_full.txt"))
META_LIST   = config.get("metagenome_scrub_list", os.path.join(scrub_db, "metagenomes_to_scrub.txt"))

# input tables (new key names, old ones still accepted)
SAMPLE_SHEET = config.get("sample_sheet", config.get("target_metagenomes"))
GENOME_LISTS = config.get("genome_list", config.get("genome_lists"))
if isinstance(GENOME_LISTS, str):
    GENOME_LISTS = [GENOME_LISTS]
if not SAMPLE_SHEET or not GENOME_LISTS:
    raise ValueError("config needs 'sample_sheet' (or 'target_metagenomes') and "
                     "'genome_list' (or 'genome_lists')")

# sourmash containment cutoff for the pangenome scrub list; "all" = whole genome database, no prefilter
CONT     = str(config.get("containment", 0.0001))
CONT_ALL = CONT == "all"
CONT_COL = config.get("containment_column", "containment")   # manysearch column the cutoff is applied to

DRUG_SCRUB     = bool(config.get("drug_scrub", True))
SCRUB_T        = config.get("scrub_t", 0.999)                                  # kmer_scrub_count -T
FILTER_MIN     = config.get("filter_fraction", config.get("filter_min", 0.01)) # kmer_scrub_filter.py -m
TERMINAL_DIST  = config.get("terminal_dist", 500)                              # coverage_depth.py --terminal_dist
RARE_MIN_COUNT = config.get("rare_kmer_count", 1000)                           # coverage_depth.py --rare_min_count


# ---------- execution environment (cluster vs local) ----------
SHELL_PREFIX = config.get("shell_prefix", "")    # e.g. "ml anaconda3/latest &&"
ENV_RUN      = config.get("env_run", "")         # e.g. "conda run --no-capture-output -n strainer3_env"
PY       = f"{SHELL_PREFIX} {ENV_RUN} python -u".strip()
SOURMASH = f"{SHELL_PREFIX} {ENV_RUN} sourmash".strip()

# per-rule resources: defaults below (the cluster values), overridden by config["resources"][rule][key]
DEFAULT_RES = {
    "target_kmdb":         dict(threads=10, mem_mb=120000, runtime=120, lsf_queue="premium"),
    "sourmash_sketch":     dict(threads=1,  mem_mb=4000,   runtime=30,  lsf_queue="premium"),
    "sourmash_search_all": dict(threads=20, mem_mb=64000,  runtime=480, lsf_queue="premium"),
    "sourmash_split_hits": dict(threads=1,  mem_mb=32000,  runtime=60,  lsf_queue="premium"),
    "kmer_scrub":          dict(threads=45, mem_mb=32000,  runtime=720, lsf_queue="premium"),
    "kmer_filter":         dict(threads=1,  mem_mb=16000,  runtime=60,  lsf_queue="premium"),
    "strain_detect":       dict(threads=20, mem_mb=32000,  runtime=120, lsf_queue="premium"),
    "coverage_depth":      dict(threads=1,  mem_mb=16000,  runtime=60,  lsf_queue="premium"),
    "hit_call":            dict(threads=1,  mem_mb=8000,   runtime=30,  lsf_queue="premium"),
}

_unknown = set(config.get("resources", {})) - set(DEFAULT_RES)
if _unknown:
    raise ValueError(f"config 'resources' has unknown rule(s): {sorted(_unknown)}; "
                     f"valid: {sorted(DEFAULT_RES)}")

def R(rule_name, key):
    return config.get("resources", {}).get(rule_name, {}).get(key, DEFAULT_RES[rule_name][key])


# ---------- helpers ----------
_INDEX = {}   # location -> {basename: [paths]}, built once per location (one os.walk instead of a glob per file)

def _index(loc):
    loc = os.path.abspath(loc)
    if loc not in _INDEX:
        idx = {}
        for dirpath, _, filenames in os.walk(loc, followlinks=True):
            for name in filenames:
                idx.setdefault(name, []).append(os.path.join(dirpath, name))
        _INDEX[loc] = idx
        print(f"[index] {loc}: {sum(map(len, idx.values()))} files", file=sys.stderr)
    return _INDEX[loc]

def resolve(fname, locations, exts=("", ".fna", ".fna.gz", ".fa", ".fa.gz", ".fasta", ".fasta.gz")):
    if os.path.isabs(fname) or os.path.exists(fname):
        return fname
    for loc in locations:
        idx = _index(loc)
        for ext in exts:
            target = fname + ext
            hits = [h for h in idx.get(os.path.basename(target), [])
                    if h.endswith(os.sep + target)]
            if len(hits) > 1:
                print(f"WARNING: {fname}{ext} matched {len(hits)} files, using {hits[0]}", file=sys.stderr)
            if hits:
                return hits[0]
    raise FileNotFoundError(
        f"{fname!r} not found under {[os.path.abspath(l) for l in locations]} "
        f"(cwd={os.getcwd()}, tried extensions {exts})"
    )

def strip_seq_ext(path):
    """GCF_000123.4_genomic.fna.gz -> GCF_000123.4_genomic (same rule as the C derive_basename)."""
    name = os.path.basename(path)
    for exts in ((".gz", ".zst", ".bz2", ".xz"), (".fa", ".fna", ".fasta", ".ffn", ".fas", ".fsa")):
        for e in exts:
            if name.lower().endswith(e):
                name = name[: -len(e)]
                break
    return name

def sd(*parts):
    """path inside the per genome-list / strain directory"""
    return os.path.join(output_dir, "genome_lists", "{glist}", "{strain}", *parts)


# ---------- wildcards: targets ----------
df_targets = pd.read_csv(SAMPLE_SHEET, sep="\t")
dict_targets = {
    row.sample_name: [resolve(f.strip(), config["metagenome_locations"])
                      for f in row.seq_files.split(",") if f.strip()]
    for row in df_targets.itertuples()
}

# ---------- wildcards: genome lists / strains ----------
# strain name = file name without sequence extensions, i.e. what the strainer3 binaries name their outputs
dfs = []
for glist_path in GENOME_LISTS:
    df_g = pd.read_csv(glist_path, header=None, names=["genome"], sep="\t")
    df_g["genome"] = df_g["genome"].str.strip()
    df_g = df_g[df_g["genome"] != ""]
    df_g["genome_list"] = os.path.splitext(os.path.basename(glist_path))[0]
    df_g["genome_path"] = df_g["genome"].apply(lambda g: resolve(g, config["genome_locations"]))
    df_g["strain"] = df_g["genome_path"].apply(strip_seq_ext)
    dfs.append(df_g)

df_glist = pd.concat(dfs, ignore_index=True)
_clash = df_glist.groupby("strain")["genome_path"].nunique()
if (_clash > 1).any():
    raise ValueError(f"strain names map to different files: {list(_clash[_clash > 1].index)[:5]}")

dict_glist    = df_glist.set_index(["genome_list", "strain"])["genome_path"].to_dict()
GLIST_PATHS   = df_glist.groupby("genome_list")["genome_path"].apply(list).to_dict()
GLIST_STRAINS = df_glist.groupby("genome_list")["strain"].apply(list).to_dict()
STRAIN_PATHS  = df_glist.drop_duplicates("strain").set_index("strain")["genome_path"].to_dict()
ALL_STRAINS   = sorted(STRAIN_PATHS)

wildcard_constraints:
    target = "|".join(re.escape(t) for t in dict_targets),
    glist  = "|".join(re.escape(g) for g in GLIST_PATHS),
    strain = "|".join(re.escape(s) for s in ALL_STRAINS),

# trivial file-writing rules run on the submit node, not as LSF jobs
localrules: all, pangenome_scrub_list, drug_scrub_list, target_list, hit_calls_concat


# ---------- 1) target k-mer databases (shared by all strains) ----------
rule target_kmdb:
    input:
        lambda wc: dict_targets[wc.target]
    output:
        os.path.join(output_dir, "target_samples_kmdb", "{target}.kmdb")
    params:
        bin    = os.path.join(src_dir, "kmer_metagenome_db"),
        outdir = lambda wc, output: os.path.dirname(output[0]),
    threads: R("target_kmdb", "threads")
    resources:
        lsf_queue = R("target_kmdb", "lsf_queue"),
        mem_mb    = R("target_kmdb", "mem_mb"),
        runtime   = R("target_kmdb", "runtime"),
    log:
        os.path.join(output_dir, "0-logs", "target_kmdb", "{target}.log")
    shell:
        """
        {params.bin} -m sample \
            -O {params.outdir} \
            -n {wildcards.target} \
            -t {threads} \
            {input} > {log} 2>&1
        """


rule target_list:
    input:
        expand(rules.target_kmdb.output[0], target=dict_targets)
    output:
        os.path.join(output_dir, "targets.txt")
    run:
        with open(output[0], "w") as fh:
            fh.write("".join(p + "\n" for p in input))


# ---------- 2a) sourmash prefilter for the pangenome scrub (skipped when containment is "all") ----------
rule sourmash_sketch:
    # one sketch per strain, shared across genome lists
    input:
        lambda wc: STRAIN_PATHS[wc.strain]
    output:
        os.path.join(output_dir, "sourmash", "sketches", "{strain}.sig")
    params:
        sourmash = SOURMASH,
    threads: R("sourmash_sketch", "threads")
    resources:
        lsf_queue = R("sourmash_sketch", "lsf_queue"),
        mem_mb    = R("sourmash_sketch", "mem_mb"),
        runtime   = R("sourmash_sketch", "runtime"),
    log:
        os.path.join(output_dir, "0-logs", "sourmash", "sketch", "{strain}.log")
    shell:
        """
        {params.sourmash} sketch dna {input} \
            -p k=31,scaled=1000 \
            --name {wildcards.strain} \
            -o {output} > {log} 2>&1
        """


rule sourmash_search_all:
    # one manysearch of all strains against the DB (DB opened once)
    input:
        sketches = expand(rules.sourmash_sketch.output[0], strain=ALL_STRAINS),
    output:
        csv      = os.path.join(output_dir, "sourmash", "all_vs_db.csv"),
        pathlist = temp(os.path.join(output_dir, "sourmash", "query_sketches.txt")),
    params:
        sourmash    = SOURMASH,
        db          = ROCKSDB,
        containment = 0 if CONT_ALL else float(CONT),
    threads: R("sourmash_search_all", "threads")
    resources:
        lsf_queue = R("sourmash_search_all", "lsf_queue"),
        mem_mb    = R("sourmash_search_all", "mem_mb"),
        runtime   = R("sourmash_search_all", "runtime"),
    log:
        os.path.join(output_dir, "0-logs", "sourmash", "search_all.log")
    run:
        with open(output.pathlist, "w") as fh:
            fh.write("\n".join(input.sketches) + "\n")
        shell(
            "{params.sourmash} scripts manysearch {output.pathlist} {params.db} "
            "-k 31 -s 1000 "
            "--threshold {params.containment} "
            "-c {threads} "
            "-o {output.csv} > {log} 2>&1"
        )


rule sourmash_split_hits:
    # the (large) search table is read once and split into one csv per strain
    input:
        csv = rules.sourmash_search_all.output.csv,
    output:
        expand(os.path.join(output_dir, "sourmash", "hits", "{strain}.csv"), strain=ALL_STRAINS),
    threads: R("sourmash_split_hits", "threads")
    resources:
        lsf_queue = R("sourmash_split_hits", "lsf_queue"),
        mem_mb    = R("sourmash_split_hits", "mem_mb"),
        runtime   = R("sourmash_split_hits", "runtime"),
    run:
        df = pd.read_csv(input.csv)
        if CONT_COL not in df.columns:
            raise ValueError(f"column '{CONT_COL}' not in {input.csv}; columns: {list(df.columns)}")
        df["_q"] = df["query_name"].astype(str).str.split().str[0]
        groups = dict(tuple(df.groupby("_q")))
        for s, out in zip(ALL_STRAINS, output):
            groups.get(s, df.iloc[0:0]).drop(columns="_q").to_csv(out, index=False)


rule pangenome_scrub_list:
    # per strain: the scrub DB genomes with containment >= cutoff (shared across genome lists)
    input:
        hits      = os.path.join(output_dir, "sourmash", "hits", "{strain}.csv"),
        pangenome = GENOME_LIST,
    output:
        os.path.join(output_dir, "sourmash", "{strain}.genomes_to_scrub.txt"),
    run:
        df = pd.read_csv(input.hits)
        df = df[df[CONT_COL] >= float(CONT)]

        # map sourmash match_name -> path in the scrub DB genome list
        with open(input.pangenome) as fh:
            by_name = {strip_seq_ext(p): p for p in (l.strip() for l in fh) if p}

        keep, missing = set(), []
        for m in df["match_name"].astype(str):
            key = strip_seq_ext(m.split()[0])
            if key in by_name:
                keep.add(by_name[key])
            else:
                missing.append(m)
        if missing:
            print(f"WARNING {wildcards.strain}/cont{CONT}: {len(missing)} sourmash matches "
                  f"not found in {input.pangenome}, e.g. {missing[:3]}", file=sys.stderr)
        print(f"{wildcards.strain}: cont {CONT} -> {len(keep)} scrub genomes", file=sys.stderr)
        with open(output[0], "w") as fh:
            fh.write("".join(p + "\n" for p in sorted(keep)))


# ---------- 2b) drug scrub list (all genomes of the list) ----------
rule drug_scrub_list:
    input:
        lambda wc: GLIST_PATHS[wc.glist]
    output:
        os.path.join(output_dir, "drug_lists", "{glist}.drug_scrub.txt")
    run:
        with open(output[0], "w") as fh:
            fh.write("".join(p + "\n" for p in input))


# ---------- 3) scrub counts ----------
rule kmer_scrub:
    input:
        ref    = lambda wc: dict_glist[(wc.glist, wc.strain)],
        # "all" -> the full DB genome list, otherwise the sourmash-filtered list
        gscrub = lambda wc: GENOME_LIST if CONT_ALL else
                 rules.pangenome_scrub_list.output[0].format(strain=wc.strain),
        mscrub = META_LIST,
        dscrub = rules.drug_scrub_list.output[0] if DRUG_SCRUB else [],
    output:
        sd("scrub", "{strain}.global_counts.tsv.zst"),
    params:
        bin     = os.path.join(src_dir, "kmer_scrub_count"),
        outdir  = lambda wc, output: os.path.dirname(output[0]),
        drug    = lambda wc, input: f"-C {input.dscrub}" if DRUG_SCRUB else "",
        scrub_t = SCRUB_T,
    threads: R("kmer_scrub", "threads")
    resources:
        lsf_queue = R("kmer_scrub", "lsf_queue"),
        mem_mb    = R("kmer_scrub", "mem_mb"),
        runtime   = R("kmer_scrub", "runtime"),
    log:
        os.path.join(output_dir, "0-logs", "kmer_scrub", "{glist}", "{strain}.log")
    shell:
        """
        {params.bin} \
            -r {input.ref} \
            -A {input.gscrub} \
            -B {input.mscrub} \
            {params.drug} \
            -t {threads} \
            -T {params.scrub_t} \
            -O {params.outdir} > {log} 2>&1
        """


# ---------- 4) filter scrubbed k-mers ----------
rule kmer_filter:
    input:
        rules.kmer_scrub.output[0],
    output:
        sd("filter", "{strain}.scrubbed_kmers.gz"),
    params:
        py     = PY,
        script = os.path.join(scripts_dir, "kmer_scrub_filter.py"),
        outdir = lambda wc, output: os.path.dirname(output[0]),
        min_f  = FILTER_MIN,
    threads: R("kmer_filter", "threads")
    resources:
        lsf_queue = R("kmer_filter", "lsf_queue"),
        mem_mb    = R("kmer_filter", "mem_mb"),
        runtime   = R("kmer_filter", "runtime"),
    log:
        os.path.join(output_dir, "0-logs", "kmer_filter", "{glist}", "{strain}.log")
    shell:
        """
        {params.py} {params.script} \
            -s {input} \
            -O {params.outdir} \
            -n {wildcards.strain} \
            -m {params.min_f} > {log} 2>&1
        """


# ---------- 5) strain detection against the target kmdbs ----------
rule strain_detect:
    input:
        kmers   = rules.kmer_filter.output[0],
        targets = rules.target_list.output[0],
        kmdbs   = expand(rules.target_kmdb.output[0], target=dict_targets),
    output:
        sd("strain_detect", "{strain}.kmer_hits.tsv.gz"),
    params:
        bin    = os.path.join(src_dir, "kmer_strain_detect"),
        outdir = lambda wc, output: os.path.dirname(output[0]),
    threads: R("strain_detect", "threads")
    resources:
        lsf_queue = R("strain_detect", "lsf_queue"),
        mem_mb    = R("strain_detect", "mem_mb"),
        runtime   = R("strain_detect", "runtime"),
    log:
        os.path.join(output_dir, "0-logs", "strain_detect", "{glist}", "{strain}.log")
    shell:
        """
        {params.bin} \
            -k {input.kmers} \
            -B {input.targets} \
            -O {params.outdir} \
            -n {wildcards.strain} \
            -t {threads} > {log} 2>&1
        """


# ---------- 6) coverage depth (reference genome, terminal distance, rare kmer count) ----------
rule coverage_depth:
    input:
        counts = rules.kmer_scrub.output[0],
        hits   = rules.strain_detect.output[0],
        ref    = lambda wc: dict_glist[(wc.glist, wc.strain)],
    output:
        sd("coverage", "coverage_depth.tsv"),
    params:
        py     = PY,
        script = os.path.join(scripts_dir, "coverage_depth.py"),
        outdir = lambda wc, output: os.path.dirname(output[0]),
        td     = TERMINAL_DIST,
        rmc    = RARE_MIN_COUNT,
    threads: R("coverage_depth", "threads")
    resources:
        lsf_queue = R("coverage_depth", "lsf_queue"),
        mem_mb    = R("coverage_depth", "mem_mb"),
        runtime   = R("coverage_depth", "runtime"),
    log:
        os.path.join(output_dir, "0-logs", "coverage_depth", "{glist}", "{strain}.log")
    shell:
        """
        {params.py} {params.script} \
            --counts {input.counts} \
            --hits {input.hits} \
            --reference {input.ref} \
            --terminal_dist {params.td} \
            --rare_min_count {params.rmc} \
            --strain {wildcards.strain} \
            --output_dir {params.outdir} > {log} 2>&1
        """


# ---------- 7) hit calls per genome list ----------
rule hit_call:
    input:
        lambda wc: expand(rules.coverage_depth.output[0], glist=wc.glist,
                          strain=GLIST_STRAINS[wc.glist]),
    output:
        os.path.join(output_dir, "genome_lists", "{glist}", "hit_calls.tsv"),
    params:
        py     = PY,
        script = os.path.join(scripts_dir, "hit_call.py"),
        covdir = lambda wc, output: os.path.dirname(output[0]),
    threads: R("hit_call", "threads")
    resources:
        lsf_queue = R("hit_call", "lsf_queue"),
        mem_mb    = R("hit_call", "mem_mb"),
        runtime   = R("hit_call", "runtime"),
    log:
        os.path.join(output_dir, "0-logs", "hit_call", "{glist}.log")
    shell:
        """
        {params.py} {params.script} \
            --coverage_dir {params.covdir} \
            --output_dir {params.covdir} > {log} 2>&1
        """


HIT_KEYS = list(GLIST_PATHS)

rule hit_calls_concat:
    # all genome-list hit calls -> one table with leading containment / genome_list columns
    input:
        [rules.hit_call.output[0].format(glist=g) for g in HIT_KEYS]
    output:
        os.path.join(output_dir, "hit_calls.tsv")
    run:
        dfs = []
        for g, p in zip(HIT_KEYS, input):
            df = pd.read_csv(p, sep="\t")
            df.insert(0, "genome_list", g)
            df.insert(0, "containment", CONT)    # cutoff as given, or "all"
            dfs.append(df)
        pd.concat(dfs, ignore_index=True).to_csv(output[0], sep="\t", index=False)


# ---------- targets (last, so it can reference the rules above) ----------
rule all:
    default_target: True
    input:
        rules.hit_calls_concat.output,