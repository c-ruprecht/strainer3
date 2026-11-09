import pandas as pd
import os
import argparse
import glob
import json
import numpy as np
import xgboost as xgb

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, '..', 'data')
# default: both models -> columns xgb2feat_prob / xgb2feat_call and xgb5feat_prob / xgb5feat_call
DEFAULT_MODELS = [os.path.join(DATA_DIR, 'xgb_hitcall_2param'), os.path.join(DATA_DIR, 'hitcaller_xgb_5feat_sylph')]

# sylph profile columns, carried over as sylph_<column> (as in the training table)
SYLPH_COLS = ['Taxonomic_abundance', 'Sequence_abundance', 'Adjusted_ANI', 'Eff_cov', 'ANI_5-95_percentile',
              'Eff_lambda', 'Lambda_5-95_percentile', 'Median_cov', 'Mean_cov_geq1', 'Containment_ind',
              'Naive_ANI', 'kmers_reassigned', 'Contig_name']
def concat_files(paths, label):
    """read and concat tsv files that carry a 'strain' column"""
    li_dfs = []
    for df_path in paths:
        df = pd.read_csv(df_path, sep='\t')
        if 'strain' in df.columns:
            li_dfs.append(df)
        else:
            print(f"Add 'strain' to {df_path} before concatenating")
    if not li_dfs:
        raise ValueError(f"Found {label} files, but none have a 'strain' column")
    return pd.concat(li_dfs, ignore_index=True)


def strip_seq_ext(path):
    """GCF_000123.4_genomic.fna.gz -> GCF_000123.4_genomic (same rule as the Snakefile / C derive_basename)"""
    name = os.path.basename(str(path))
    for exts in ((".gz", ".zst", ".bz2", ".xz"), (".fa", ".fna", ".fasta", ".ffn", ".fas", ".fsa")):
        for e in exts:
            if name.lower().endswith(e):
                name = name[: -len(e)]
                break
    return name


def reads_to_target(path):
    """sylph Sample_file -> target name. The workflow sketches the seqtk output {target}_R1.fastq.gz (PE)
    or {target}.fastq.gz (SE); a .kmdb suffix is stripped too, as in the coverage table"""
    name = os.path.basename(str(path))
    for suf in ('.gz', '.fastq', '.fq', '.kmdb'):
        if name.endswith(suf):
            name = name[: -len(suf)]
    for suf in ('_R1', '_R2'):
        if name.endswith(suf):
            name = name[: -len(suf)]
    return name


def read_sylph(paths):
    """sylph profile tsv(s) -> one row per strain x sample with sylph_<col> columns.
    strain / sample are taken from 'strain' / 'sample' columns if present, otherwise derived from
    Genome_file and Sample_file with the workflow's naming rules (see strip_seq_ext / reads_to_target)."""
    dfs = []
    for p in paths:
        try:
            df = pd.read_csv(p, sep='\t')
        except pd.errors.EmptyDataError:   # sylph profile writes nothing when no genome passes its cutoffs
            df = pd.DataFrame()
        if df.empty:
            print(f'WARNING: {p} is empty (sylph found none of the genomes in any sample)')
            continue
        if 'strain' not in df.columns:
            df['strain'] = df['Genome_file'].map(strip_seq_ext)
        if 'sample' not in df.columns:
            df['sample'] = df['Sample_file'].map(reads_to_target)
        dfs.append(df)
    cols = ['strain', 'sample'] + SYLPH_COLS
    if not dfs:
        return pd.DataFrame(columns=cols).rename(columns={c: f'sylph_{c}' for c in SYLPH_COLS})
    df = pd.concat(dfs, ignore_index=True)
    df = df[['strain', 'sample'] + [c for c in SYLPH_COLS if c in df.columns]]
    dup = df.duplicated(['strain', 'sample'])
    if dup.any():
        print(f'WARNING: {dup.sum()} duplicate strain x sample rows in sylph output, keeping the first')
        df = df[~dup]
    return df.rename(columns={c: f'sylph_{c}' for c in SYLPH_COLS})


def add_sylph(df_cov, df_sylph):
    """left-join sylph onto the coverage table. Genomes sylph did not report stay NaN, which is how the model
    was trained (NaN = not detected by sylph; do NOT fill with 0). Hits to background DBs (strains not in
    the coverage table) drop out in the join."""
    df_sylph = df_sylph.copy()
    for c in df_sylph.columns:
        if c.startswith('sylph_') and c not in ('sylph_Contig_name', 'sylph_Containment_ind',
                                                 'sylph_ANI_5-95_percentile', 'sylph_Lambda_5-95_percentile'):
            df_sylph[c] = pd.to_numeric(df_sylph[c], errors='coerce')   # Eff_lambda can be 'HIGH' -> NaN
    df_cov = df_cov.drop(columns=[c for c in df_cov.columns if c.startswith('sylph_')])
    key_sample = df_cov['sample'].map(reads_to_target)                  # tolerate '<target>.kmdb' sample names
    out = df_cov.assign(_sample=key_sample).merge(df_sylph.rename(columns={'sample': '_sample'}),
                                                  on=['strain', '_sample'], how='left', validate='many_to_one')
    out = out.drop(columns='_sample')

    if df_sylph.empty:
        print('WARNING: no sylph rows; all sylph features NaN (treated as "not detected by sylph")')
        return out
    # sanity checks on the join keys: a naming mismatch silently turns every sylph value into NaN
    cov_samples, sy_samples = set(key_sample), set(df_sylph['sample'])
    cov_strains, sy_strains = set(df_cov['strain']), set(df_sylph['strain'])
    print(f"sylph: {len(df_sylph)} rows | samples matched {len(cov_samples & sy_samples)}/{len(cov_samples)} | "
          f"strains matched {len(cov_strains & sy_strains)}/{len(cov_strains)} | "
          f"strain x sample with a sylph hit {out['sylph_Naive_ANI'].notna().sum()}/{len(out)}")
    if not cov_samples & sy_samples or not cov_strains & sy_strains:
        raise ValueError('no sylph strain/sample names match the coverage table, e.g. coverage '
                         f"{sorted(cov_samples)[:3]} / {sorted(cov_strains)[:2]} vs sylph "
                         f"{sorted(sy_samples)[:3]} / {sorted(sy_strains)[:2]}. Add 'strain' and 'sample' "
                         'columns to the sylph tsv that match the coverage table.')
    if cov_samples - sy_samples:
        print(f'note: samples without any sylph hit to these strains: {sorted(cov_samples - sy_samples)[:10]}')
    return out


def model_label(meta):
    """column prefix of a model: xgb<n features>feat, e.g. xgb2feat / xgb5feat"""
    return f"xgb{len(meta['features'])}feat"


def parse_model_arg(arg):
    """'<prefix>' or '<label>=<prefix>' -> (label or None, prefix)"""
    if '=' in arg and not os.path.exists(arg + '.json'):
        label, prefix = arg.split('=', 1)
        return label, prefix
    return None, arg


def load_model(prefix):
    """load the xgb model (<prefix>.json) and its metadata (<prefix>_meta.json) exported by the notebook"""
    with open(f'{prefix}_meta.json') as fh:
        meta = json.load(fh)
    model = xgb.XGBClassifier()
    model.load_model(f'{prefix}.json')
    return model, meta


def check_ranges(df, meta):
    """warn when inputs fall outside what the model saw in training (e.g. samples deeper than ~8M reads)"""
    for f, (lo, hi) in meta.get('feature_ranges', {}).items():
        x = df[f].replace([np.inf, -np.inf], np.nan)
        n_lo, n_hi = int((x < lo).sum()), int((x > hi).sum())
        if n_lo or n_hi:
            print(f'WARNING: {f}: {n_lo} rows below / {n_hi} rows above the training range [{lo:.4g}, {hi:.4g}]'
                  + (' -> downsample reads into the training range' if f == 'total_evaluated' and n_hi else ''))


def xgb_hit_call(df, model, meta, cutoff_rule=None, label='xgb'):
    """add <label>_prob (P(present)) and <label>_call (0/1) to df using the features stored in meta"""
    features = meta['features']
    missing = [f for f in features if f not in df.columns]
    if missing:
        hint = ' (sylph features: give --sylph_file)' \
            if any(f.startswith('sylph_') for f in missing) else ''
        raise KeyError(f'model needs columns {missing}, not found in coverage table{hint}')

    if cutoff_rule:
        if cutoff_rule not in meta.get('cutoffs', {}):
            raise KeyError(f"cutoff rule '{cutoff_rule}' not in model, choose from {list(meta['cutoffs'])}")
        threshold = meta['cutoffs'][cutoff_rule]['threshold']
    else:
        cutoff_rule, threshold = meta['cutoff_rule'], meta['threshold']

    check_ranges(df, meta)
    # same preparation as in training: inf -> +/-1e12, NaN kept (xgboost routes missing values)
    X = df[features].replace([np.inf, -np.inf], [1e12, -1e12]).astype(float)
    df = df.copy()
    prob, call = f'{label}_prob', f'{label}_call'
    df[prob] = model.predict_proba(X)[:, 1]
    df[call] = (df[prob] >= threshold).astype(int)
    print(f"{label} hit call: model {meta.get('model_name', '?')}, features {features}, rule '{cutoff_rule}', "
          f"threshold {threshold:.3f} -> {df[call].sum()} of {len(df)} strain x sample calls")
    return df


def main():
    # to do, basename + strain name addition
    parser = argparse.ArgumentParser(description='Using coverage and depth to make hit calls')
    parser.add_argument('--coverage_file', help='a coverage_depth.tsv file')
    parser.add_argument('--hist_file', help='a zi_hist.tsv.gz file (coverage_depth_locations.py --export_hist)')
    parser.add_argument('--sylph_file', nargs='+',
                        help='sylph profile tsv(s) of the same samples (same depth) against the same genomes, '
                        'e.g. ds{ds}/sylph/{glist}.sylph_profile.tsv; needed if the model uses sylph_* features')
    parser.add_argument('--coverage_dir', help = "Directory in which coverage files are contained, will discover "
                        "'coverage_depth.tsv' and 'zi_hist.tsv.gz' files in subdirectories and concat them")
    parser.add_argument('--model', nargs='+', default=DEFAULT_MODELS,
                        help = 'one or more xgb model prefixes (<prefix>.json + <prefix>_meta.json). Each model '
                        'writes <label>_prob and <label>_call; label = xgb<n features>feat (e.g. xgb2feat, xgb5feat), '
                        'or give it explicitly as <label>=<prefix>. Default: the 2- and 5-feature models in data/')
    parser.add_argument('--cutoff_rule', nargs='+', default=['precision >= 0.95'],
                        help = "probability cutoff stored in the model meta, e.g. 'max F1' or 'precision >= 0.99' "
                        "(default: 'precision >= 0.95'). One rule for all models, or one per --model in the same order")
    parser.add_argument('--output_dir', help='directory where output is saved')
    parser.add_argument('--verbose',action='store_true', help = 'incluces all possible metrics from coverage_depth into hit call table')
    args = parser.parse_args()

    if not (args.coverage_file or args.coverage_dir):
        parser.error('give --coverage_file or --coverage_dir')

    os.makedirs(args.output_dir, exist_ok=True)

    # load the models first so a missing input fails before the slow file discovery
    if len(args.cutoff_rule) not in (1, len(args.model)):
        parser.error(f'give one --cutoff_rule for all models or one per model ({len(args.model)})')
    rules = args.cutoff_rule * len(args.model) if len(args.cutoff_rule) == 1 else args.cutoff_rule
    models = []   # (label, prefix, model, meta, cutoff rule)
    for arg, rule in zip(args.model, rules):
        label, prefix = parse_model_arg(arg)
        model, meta = load_model(prefix)
        models.append((label or model_label(meta), prefix, model, meta, rule))
    labels = [m[0] for m in models]
    if len(set(labels)) != len(labels):
        parser.error(f'models share a column label {labels}; name them with --model <label>=<prefix>')
    sylph_models = [p for _, p, _, meta, _ in models if any(f.startswith('sylph_') for f in meta['features'])]

    df_cov = None
    df_hist = None
    sylph_paths = list(args.sylph_file or [])

    if args.coverage_file:
        df_cov = pd.read_csv(args.coverage_file, sep='\t')

    if args.hist_file:
        df_hist = pd.read_csv(args.hist_file, sep='\t')

    if args.coverage_dir:
        coverage_dir = args.coverage_dir
        print(f"searching {coverage_dir} (exists: {os.path.isdir(coverage_dir)})")

        glob_cov = glob.glob(os.path.join(coverage_dir, '**', '*coverage_depth.tsv'), recursive=True)
        print(glob_cov)
        if not glob_cov:
            raise FileNotFoundError(f"No *coverage_depth.tsv found under {coverage_dir}")
        df_cov = concat_files(glob_cov, 'coverage')

        # zero-inflation histograms, only present if coverage_depth_locations.py ran with --export_hist
        glob_hist = glob.glob(os.path.join(coverage_dir, '**', '*zi_hist.tsv.gz'), recursive=True)
        print(glob_hist)
        if glob_hist:
            df_hist = concat_files(glob_hist, 'histogram')

            # fitting settings per run: histograms can be refit with any settings, but the exported metrics
            # in coverage_depth.tsv are only comparable if all runs used the same ones
            glob_params = glob.glob(os.path.join(coverage_dir, '**', 'zi_params.json'), recursive=True)
            params = {}
            for p in glob_params:
                with open(p) as fh:
                    params[p] = json.dumps(json.load(fh), sort_keys=True)
            if len(set(params.values())) > 1:
                print('WARNING: runs used different zero-inflation settings:')
                for p, v in params.items():
                    print(f'  {p}: {v}')

    if df_hist is not None:
        hist_path = os.path.join(args.output_dir, 'zi_hist.tsv.gz')
        df_hist.to_csv(hist_path, sep='\t', index=None)
        print(f"Saved {df_hist[['strain', 'sample']].drop_duplicates().shape[0]} strain x sample histograms to \n {hist_path}")

    df_cov['depth_ratio'] = df_cov['depth_rare_kmers']/df_cov['depth_nonrare_kmers']

    if sylph_paths:
        df_cov = add_sylph(df_cov, read_sylph(sylph_paths))
    elif sylph_models:
        parser.error(f"model(s) {sylph_models} use sylph features, give --sylph_file "
                     f"(or only --model {DATA_DIR}/xgb_hitcall_2param for the strainer-only model)")

    # make the hit calls, one <label>_prob / <label>_call pair per model
    for label, prefix, model, meta, rule in models:
        df_cov = xgb_hit_call(df_cov, model, meta, rule, label=label)
    out_cols = [c for label in labels for c in (f'{label}_prob', f'{label}_call')]
    if len(models) > 1:
        calls = df_cov[[f'{l}_call' for l in labels]]
        print('call agreement (rows called by each combination of models):')
        print(calls.value_counts().rename('n').to_frame().to_string())

    # non verbose: the models' input columns, their probabilities and calls
    features = list(dict.fromkeys(f for *_, meta, _ in models for f in meta['features']))
    df_hits = df_cov[['strain', 'sample'] + features + out_cols].copy()
    if args.verbose:
        df_hits = df_cov.copy()
    print(df_hits)
    # sorted by the last model given (default: the 5-feature model)
    df_hits.sort_values(f'{labels[-1]}_prob', ascending= False).to_csv(os.path.join(args.output_dir, 'hit_calls.tsv'), sep = '\t', index = None)


if __name__ == "__main__":
        main()