#!/usr/bin/env python
"""
kmer_scrub_filter.py

Takes the output of kmer_scrub_count(_individual) (every k-mer of a strain with its
counts in the reference, pangenomes, metagenomes and drug strains) and selects the
most informative (rare) k-mers.

Block mode (--reference)
------------------------
When a reference genome is given, the selected k-mers are placed on it with
Aho-Corasick. EVERY occurrence is kept, so k-mers present in up to
--max_ref_copies copies of the reference are retained, and each copy gets its own
locus. Loci are grouped into:

  segment_id  run of consecutive k-mer positions on one contig (diff == 1).
              Reads are validated against segments (diagonal / completeness check).
  block_id    sub-run of a segment with identical ref_copies, pangenome_count,
              metagenome_count and drug_count. Statistics are reported per block.
  family_id   segments that share at least one k-mer (copies of a repeat) are
              merged into one family. Unique segments are their own family.
              Reads that cannot be placed on one copy are counted per family.

Outputs (in --output_dir, prefixed with --basename):
  .scrubbed_kmers.gz    '#kmer' only, one row per unique k-mer (legacy format)
  .kmer_loci.tsv.gz     one row per locus (block mode only):
                        #kmer contig_id kmer_position ref_rc ref_copies
                        segment_id block_id family_id
                        ref_rc = True if the #kmer string, as written, occurs on the
                        reverse strand of the reference at this locus.
                        kmer_position is the 0-based start on the forward strand.
  .blocks.tsv.gz        one row per block (block mode only):
                        block_id segment_id family_id contig_id start stop
                        length_bp n_kmers ref_copies pangenome_count
                        metagenome_count drug_count
                        start/stop are 0-based inclusive; stop is the end of the
                        last k-mer.
"""

import argparse
import sys
import gzip
import io
import json
import subprocess
import os

import re

import numpy as np
import pandas as pd

_ACGT = re.compile(r'[ACGT]+')

COUNT_COLS = ['pangenome_count', 'metagenome_count', 'drug_count']
_RC = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')


def revcomp(s):
    return s.translate(_RC)[::-1]


def open_counts(path):
    """Open a kmer_scrub_count(_individual) table as text.
    Handles .zst (global_counts.tsv.zst), .gz, and plain text."""
    if path.endswith('.zst'):
        try:
            import zstandard
            fh = open(path, 'rb')
            return io.TextIOWrapper(zstandard.ZstdDecompressor().stream_reader(fh), encoding='utf-8')
        except ImportError:
            proc = subprocess.Popen(['zstd', '-dc', path], stdout=subprocess.PIPE, text=True)
            return proc.stdout
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path)


# ── original hash-based scrubbing (unchanged logic) ──────────────────────────

def scrub_max_kmers(min_frac_to_keep, kmer_hash, total_kmers):
    min_kmer_count = -1
    fraction_kept = -1.0
    kmers_to_scrub = {}
    total_kmers = float(total_kmers)

    # keep making the threshold more stringent until the fraction kept is reached
    while fraction_kept < min_frac_to_keep:
        min_kmer_count += 1
        kmer_hits = 0
        for key, val in kmer_hash.items():
            if val > min_kmer_count:
                kmer_hits += 1
        fraction_kept = 1 - (kmer_hits / total_kmers)
        sys.stderr.write("kept " + str(fraction_kept) + " with threshold " + str(min_kmer_count) + "\n")

    for key, val in kmer_hash.items():
        if val > min_kmer_count:
            kmers_to_scrub[key] = val

    sys.stderr.write("threshold was " + str(min_kmer_count) + " left with " + str(len(kmers_to_scrub))
                     + " out of " + str(total_kmers) + " that will be scrubbed\n")
    return kmers_to_scrub


def drug_scrub(drug_genome_hash, strain_hash, all_kmers):
    strain_hash_scrubbed = strain_hash
    for key in drug_genome_hash:
        if key in strain_hash_scrubbed:
            del strain_hash_scrubbed[key]
    return strain_hash_scrubbed


def independent_scrub(min_fraction, pangenome_hash, metagenome_hash, strain_hash, all_kmers):
    strain_hash_scrubbed = strain_hash
    kmers_to_scrub_pan = scrub_max_kmers(min_fraction, pangenome_hash, all_kmers)
    kmers_to_scrub_meta = scrub_max_kmers(min_fraction, metagenome_hash, all_kmers)

    for key in kmers_to_scrub_pan:
        if key in strain_hash_scrubbed:
            del strain_hash_scrubbed[key]
    for key in kmers_to_scrub_meta:
        if key in strain_hash_scrubbed:
            del strain_hash_scrubbed[key]
    return strain_hash_scrubbed


def joint_scrub(min_fraction, pangenome_hash, metagenome_hash, strain_hash, all_kmers, num_drug_scrubbed):
    # NOTE: mutates pangenome_hash, metagenome_hash and strain_hash values (to fractions)

    metagenome_sum = sum(metagenome_hash.values())
    for key in metagenome_hash:
        metagenome_hash[key] /= float(metagenome_sum)

    pangenome_sum = sum(pangenome_hash.values())
    for key in pangenome_hash:
        pangenome_hash[key] /= float(pangenome_sum)

    # combine the two hashes so we remove the kmers in order of frequency
    for key in strain_hash:
        strain_hash[key] = 0
        if key in metagenome_hash and metagenome_hash[key] > strain_hash[key]:
            strain_hash[key] = metagenome_hash[key]
        if key in pangenome_hash and pangenome_hash[key] > strain_hash[key]:
            strain_hash[key] = pangenome_hash[key]

    pairs = sorted(strain_hash.items(), key=lambda x: x[1], reverse=True)

    # remove them from most to least frequent
    num_scrubbed = float(num_drug_scrubbed)
    for pair in pairs:
        if (1 - ((num_scrubbed + 1) / all_kmers)) > min_fraction:
            num_scrubbed += 1.0
            if pair[0] in strain_hash:
                del strain_hash[pair[0]]

    return strain_hash


# ── percentile selection ─────────────────────────────────────────────────────

def get_lowest_percentile(df, percentile=0.01, drug_scrub='count_hard', max_ref_copies=5):
    """Select k-mers in the rarest percentile of pangenome and metagenome counts.

    The reference is no longer filtered by percentile: any k-mer present in
    1..max_ref_copies copies of the reference is kept (use max_ref_copies=1 for the
    old single-copy behaviour)."""
    df = df.copy()
    if 'drug_count' not in df.columns or df['drug_count'].isna().all():
        print('#No drug counts found, continuing without drugscrub')
        df['drug_count'] = 0

    ref_ok = df['reference_count'].between(1, max_ref_copies)
    rare = (
        (df['pangenome_count'] <= df['pangenome_count'].quantile(percentile)) &
        (df['metagenome_count'] <= df['metagenome_count'].quantile(percentile))
    )

    if drug_scrub == 'percentile':
        drug_ok = df['drug_count'] <= df['drug_count'].quantile(percentile)
    elif drug_scrub == 'count_hard':
        drug_ok = df['drug_count'] == 0
    else:
        raise ValueError(f'unknown drug_scrub mode: {drug_scrub}')

    return df[ref_ok & rare & drug_ok].copy()


# ── reference mapping and block annotation ───────────────────────────────────

def load_genome(genome_path):
    from Bio import SeqIO
    opener = gzip.open if genome_path.endswith('.gz') else open
    records = {}
    with opener(genome_path, 'rt') as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            records[record.id] = str(record.seq).upper()
    return records


def map_kmers_all_loci(records, kmers):
    """Find EVERY occurrence of each k-mer (either strand) in the reference.

    Returns one row per locus: contig_id, #kmer, kmer_position (0-based, forward
    strand), ref_rc (True if the k-mer as given occurs on the reverse strand)."""
    import ahocorasick

    # deduplicate by canonical form so a k-mer and its reverse complement in the
    # input cannot overwrite each other in the automaton
    seen = set()
    uniq = []
    for kmer in kmers:
        rc = revcomp(kmer)
        canon = max(kmer, rc)
        if canon in seen:
            continue
        seen.add(canon)
        uniq.append(kmer)
    if len(uniq) < len(kmers):
        print(f'#WARNING: {len(kmers) - len(uniq)} input k-mers were reverse-complement duplicates')

    A = ahocorasick.Automaton()
    for kmer in uniq:
        rc = revcomp(kmer)
        if rc != kmer:                       # palindromes: add once, as forward
            A.add_word(rc, (kmer, True))
        A.add_word(kmer, (kmer, False))
    A.make_automaton()

    rows = []
    for cid, seq in records.items():
        for end, (kmer, is_rc) in A.iter(seq):
            rows.append((cid, kmer, end - len(kmer) + 1, is_rc))

    return pd.DataFrame(rows, columns=['contig_id', '#kmer', 'kmer_position', 'ref_rc'])


def _assign_segments(df, max_kmers=0, mode='fixed'):
    """df sorted by contig_id, kmer_position. Returns 1-based segment ids.
    A segment is a run of consecutive k-mer positions on one contig. With max_kmers > 0, longer runs are cut:
      fixed    : consecutive chunks of max_kmers k-mers, the last chunk keeps the remainder (which can be
                 dropped by the minimum segment size)
      balanced : floor(n / max_kmers) chunks of (almost) equal size, i.e. between max_kmers and
                 2 * max_kmers - 1 k-mers each: no remainder, nothing of a run >= max_kmers is lost
    so every segment has (about) the same exposure and the same error-free span a read needs."""
    new_contig = df['contig_id'].ne(df['contig_id'].shift())
    pos_break = df['kmer_position'].diff().ne(1)
    brk = new_contig | pos_break
    if max_kmers and max_kmers > 0:
        run = brk.cumsum()
        rank = brk.astype(int).groupby(run).cumcount()
        if mode == 'balanced':
            n = rank.groupby(run).transform('size')
            m = np.maximum(n // max_kmers, 1)
            chunk = (rank * m) // n
            brk = brk | (chunk.ne(chunk.shift()) & ~brk)
        else:
            brk = brk | (rank % max_kmers == 0)
    return brk.cumsum()


def _repeat_families(df):
    """Union segments that share a k-mer. Returns dict segment_id -> family_id
    (family_id = smallest segment_id in the family)."""
    parent = {s: s for s in df['segment_id'].unique()}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    multi = df.loc[df['ref_copies'] > 1, ['#kmer', 'segment_id']]
    for segs in multi.groupby('#kmer')['segment_id'].unique():
        a = find(segs[0])
        for s in segs[1:]:
            b = find(s)
            if a != b:
                lo, hi = min(a, b), max(a, b)
                parent[hi] = lo
                a = lo
    return {s: find(s) for s in parent}


def _quantiles(x, prefix):
    x = pd.Series(x)
    if len(x) == 0:
        return {}
    q = x.quantile([0.05, 0.25, 0.5, 0.75, 0.95])
    return {f'{prefix}_n': int(len(x)), f'{prefix}_min': int(x.min()), f'{prefix}_q05': float(q[0.05]),
            f'{prefix}_q25': float(q[0.25]), f'{prefix}_median': float(q[0.5]), f'{prefix}_q75': float(q[0.75]),
            f'{prefix}_q95': float(q[0.95]), f'{prefix}_max': int(x.max())}


def build_blocks(df_sel, records, k, terminal_dist=0, max_ref_copies=5, min_segment_kmers=1,
                 max_segment_kmers=0, chunk_mode='fixed'):
    """Place selected k-mers on the reference and annotate segments, blocks and
    repeat families. Whenever a locus has to be dropped, the WHOLE k-mer is dropped
    (all copies), so a read can never be placed on the wrong copy because the
    true copy was removed."""
    stats = {}
    df_loci = map_kmers_all_loci(records, df_sel['#kmer'].tolist())

    n_in = df_sel['#kmer'].nunique()
    n_found = df_loci['#kmer'].nunique()
    stats['input_kmers'] = n_in
    stats['not_found'] = n_in - n_found
    if n_found < n_in:
        print(f'#WARNING: {n_in - n_found} / {n_in} selected k-mers not found in the reference')

    # observed copy number in the reference
    df_loci['ref_copies'] = df_loci.groupby('#kmer')['#kmer'].transform('size')
    ref_expected = df_sel.set_index('#kmer')['reference_count']
    mismatch = df_loci.drop_duplicates('#kmer')
    mismatch = mismatch[mismatch['ref_copies'].to_numpy()
                        != ref_expected.reindex(mismatch['#kmer']).to_numpy()]
    stats['copy_mismatch'] = len(mismatch)
    if len(mismatch):
        print(f'#WARNING: {len(mismatch)} k-mers have a different copy number in the reference '
              f'than reference_count (using the observed number)')

    def drop(mask_kmers, label):
        # whole k-mers are dropped (all copies); loci removed only because another copy failed are counted
        nonlocal df_loci
        bad = set(df_loci.loc[mask_kmers, '#kmer'])
        stats[label] = stats.get(label, 0) + len(bad)
        if bad:
            gone = df_loci['#kmer'].isin(bad)
            extra = int(gone.sum() - mask_kmers.sum())
            if extra > 0:
                stats[f'{label}_other_copy_loci'] = stats.get(f'{label}_other_copy_loci', 0) + extra
            df_loci = df_loci[~gone]

    # too many copies
    drop(df_loci['ref_copies'] > max_ref_copies, 'too_many_copies')

    # tandem / overlapping copies on the same contig (low complexity)
    df_loci = df_loci.sort_values(['#kmer', 'contig_id', 'kmer_position'])
    same = df_loci['#kmer'].eq(df_loci['#kmer'].shift()) & df_loci['contig_id'].eq(df_loci['contig_id'].shift())
    close = same & (df_loci['kmer_position'] - df_loci['kmer_position'].shift() < k)
    drop(close, 'tandem')

    # terminal k-mers (any copy terminal -> drop the k-mer)
    if terminal_dist > 0:
        clen = df_loci['contig_id'].map({c: len(s) for c, s in records.items()})
        terminal = (df_loci['kmer_position'] < terminal_dist) | \
                   ((clen - (df_loci['kmer_position'] + k)) < terminal_dist)
        drop(terminal, 'terminal')

    stats['kmers_before_segments'] = int(df_loci['#kmer'].nunique())
    df_loci = df_loci.sort_values(['contig_id', 'kmer_position']).reset_index(drop=True)
    stats.update(_quantiles(df_loci.groupby(_assign_segments(df_loci)).size(), 'runs_before'))

    # counts per k-mer
    df_loci = df_loci.merge(df_sel[['#kmer'] + COUNT_COLS], on='#kmer', how='left')
    df_loci[COUNT_COLS] = df_loci[COUNT_COLS].fillna(0).astype(np.int64)

    # segments (optionally cut into chunks of max_segment_kmers); iteratively drop k-mers that sit in
    # too-short segments / chunks. Remainder chunks below min_segment_kmers are dropped the same way.
    while True:
        df_loci = df_loci.sort_values(['contig_id', 'kmer_position']).reset_index(drop=True)
        df_loci['segment_id'] = _assign_segments(df_loci, max_segment_kmers, chunk_mode)
        if min_segment_kmers <= 1:
            break
        seg_size = df_loci.groupby('segment_id')['segment_id'].transform('size')
        small = seg_size < min_segment_kmers
        if not small.any():
            break
        # short runs vs remainders of long runs cut into chunks (fixed mode), counted separately
        run = _assign_segments(df_loci)
        run_size = run.map(run.value_counts())
        remainder = small & (run_size >= min_segment_kmers)
        if remainder.any():
            drop(remainder, 'chunk_remainder')
            continue
        drop(small, 'short_segment')

    # blocks: break inside a segment whenever any count column or copy number changes
    break_cols = ['ref_copies'] + COUNT_COLS
    brk = df_loci['segment_id'].ne(df_loci['segment_id'].shift())
    for c in break_cols:
        brk |= df_loci[c].ne(df_loci[c].shift())
    df_loci['block_id'] = brk.cumsum()

    # repeat families
    fam = _repeat_families(df_loci)
    df_loci['family_id'] = df_loci['segment_id'].map(fam)
    # renumber families 1..F in genome order
    fam_order = {f: i + 1 for i, f in enumerate(pd.unique(df_loci['family_id']))}
    df_loci['family_id'] = df_loci['family_id'].map(fam_order)

    df_blocks = (
        df_loci.groupby('block_id', sort=True)
        .agg(segment_id=('segment_id', 'first'),
             family_id=('family_id', 'first'),
             contig_id=('contig_id', 'first'),
             start=('kmer_position', 'min'),
             stop=('kmer_position', 'max'),
             n_kmers=('#kmer', 'size'),
             ref_copies=('ref_copies', 'first'),
             pangenome_count=('pangenome_count', 'first'),
             metagenome_count=('metagenome_count', 'first'),
             drug_count=('drug_count', 'first'))
        .reset_index()
    )
    df_blocks['stop'] = df_blocks['stop'] + k - 1
    df_blocks['length_bp'] = df_blocks['stop'] - df_blocks['start'] + 1
    df_blocks = df_blocks[['block_id', 'segment_id', 'family_id', 'contig_id', 'start', 'stop',
                           'length_bp', 'n_kmers', 'ref_copies'] + COUNT_COLS]

    df_loci = df_loci[['#kmer', 'contig_id', 'kmer_position', 'ref_rc', 'ref_copies',
                       'segment_id', 'block_id', 'family_id']]

    stats['kept_kmers'] = df_loci['#kmer'].nunique()
    stats['loci'] = len(df_loci)
    stats['multi_copy_kmers'] = int((df_loci.drop_duplicates('#kmer')['ref_copies'] > 1).sum())
    stats['contiguous_runs'] = int(_assign_segments(df_loci).nunique())
    stats['segments'] = df_loci['segment_id'].nunique()
    stats['blocks'] = len(df_blocks)
    stats['families'] = df_loci['family_id'].nunique()
    rep = df_blocks['ref_copies'] > 1
    stats['repeat_families'] = int(df_blocks.loc[rep, 'family_id'].nunique())
    stats['repeat_blocks'] = int(rep.sum())
    stats['repeat_loci'] = int(df_blocks.loc[rep, 'n_kmers'].sum())
    stats['repeat_loci_frac'] = round(stats['repeat_loci'] / max(len(df_loci), 1), 5)
    stats['blocks_snp_size'] = int(((df_blocks['n_kmers'] == k) & ~rep).sum())
    stats['contigs_with_blocks'] = int(df_blocks['contig_id'].nunique())
    stats['kmer_retention'] = round(stats['kept_kmers'] / max(stats['input_kmers'], 1), 5)
    stats.update(_quantiles(df_loci.groupby(_assign_segments(df_loci)).size(), 'runs'))
    stats.update(_quantiles(df_loci.groupby('segment_id').size(), 'segments'))
    stats.update(_quantiles(df_blocks['n_kmers'], 'blocks'))
    stats.update(_quantiles(df_blocks.loc[rep, 'n_kmers'], 'repeat_blocks'))
    # spacing between consecutive blocks on a contig, in k-mer start positions (0 = back to back): the
    # exposure of the 'best' count depends on it (isolated blocks need a gap >= read length - k)
    gaps = (df_blocks.sort_values(['contig_id', 'start']).groupby('contig_id')
            .apply(lambda d: (d['start'].to_numpy()[1:] - (d['stop'].to_numpy()[:-1] - k + 1) - 1),
                   include_groups=False))
    gaps = np.concatenate([g for g in gaps if len(g)]) if len(gaps) else np.array([])
    stats.update(_quantiles(gaps, 'block_gap_kmers'))
    return df_loci, df_blocks, stats


def select_by_count(df, min_kmers=1000, max_count=1, max_ref_copies=5):
    """Keep k-mers never seen outside the strain: scrub count (pangenome + metagenome) == 0, drug count == 0,
    1..max_ref_copies copies in the reference. If that leaves fewer than min_kmers, relax the scrub count
    step by step (<= 1, <= 2, ... up to max_count). Returns the selection and the count threshold used."""
    df = df.copy()
    df['scrub_count'] = df['pangenome_count'] + df['metagenome_count']
    base = df['reference_count'].between(1, max_ref_copies) & (df['drug_count'] == 0)
    for t in range(0, max_count + 1):
        sel = df[base & (df['scrub_count'] <= t)]
        print(f"#count selection: scrub count <= {t}: {len(sel)} kmers")
        if len(sel) >= min_kmers:
            break
    if len(sel) < min_kmers:
        print(f"#WARNING: only {len(sel)} kmers with scrub count <= {max_count} (wanted {min_kmers}); "
              f"keeping them all")
    return sel.drop(columns='scrub_count').copy(), t


# ── main ─────────────────────────────────────────────────────────────────────

def read_counts_df(path):
    df = pd.read_csv(open_counts(path), sep='\t')
    df = df.rename(columns={df.columns[0]: '#kmer'})
    if 'drug_count' not in df.columns:
        df['drug_count'] = 0
    return df


def drop_non_acgt(df):
    """k-mers with ambiguity codes (N, Y, R, ... in the reference) can never match a read: remove them BEFORE
    selection, so they neither count towards --min_kmers / the percentiles nor end up in segments or blocks"""
    ok = df['#kmer'].str.upper().str.fullmatch('[ACGT]+')
    n_bad = int((~ok).sum())
    if n_bad:
        print(f"#removed {n_bad} kmers with non-ACGT characters (ambiguity codes in the reference) before selection")
    df = df.loc[ok].copy()
    df['#kmer'] = df['#kmer'].str.upper()
    return df, n_bad


def read_counts_list(paths):
    """several count files for the same strain (e.g. scrub run in batches): sum the scrub counts,
    reference counts must be identical"""
    df = None
    for path in paths:
        d = read_counts_df(path).set_index('#kmer')
        if df is None:
            df = d
            continue
        if not d.index.equals(df.index) or not (d['reference_count'] == df['reference_count']).all():
            sys.exit(f"error: {path} does not have the same kmers / reference counts as the first file")
        for c in COUNT_COLS:
            df[c] = df[c] + d[c]
    return df.reset_index()


def main():
    parser = argparse.ArgumentParser(description="\
    DESCRIPTION: Script will take as input the output from the program kmer_scrub_count \n\
    that lists every kmer in a strain of interest and the frequency of that kmer in a \n\
    set of pangenomes and a set of metagenomes.\n\
    With this program we then filter the set of kmers for the strain to select those that are most informative.\n\
    With --reference, the selected kmers are placed on the genome (all copies) and annotated with\n\
    segment, block and repeat-family ids for block-aware strain detection.\n\
    ")

    parser.add_argument("--scrub_count_file", "-s", help="Input file with the kmer counts to the pangenome and metagenomes")
    parser.add_argument("--scrub_count_list", "-l", help="A text file containing a list of input files containing kmer counts")
    parser.add_argument("--min_fraction", "-m", help="The minimum fraction of kmers to keep (so don't scrub everything); "
                        "default 0.04; range (0.0-1.0)", default=0.04, type=float)
    parser.add_argument("--independent", "-i", action='store_true',
                        help="Independently scrub metagenome and pangenome (without this flag, it scrubs them "
                             "together so there is definitely -m min_fraction left)")
    parser.add_argument("--output_dir", "-O", help="Output directory for scrubbed kmer file", required=True)
    parser.add_argument("--basename", '-n', help="Rename filename")
    parser.add_argument("--selection", choices=['count', 'percentile', 'scrub'], default='count',
                        help="How rare kmers are chosen. count (default): scrub count (pangenome + metagenome) == 0 "
                             "and no drug hits; if fewer than --min_kmers, add scrub count 1 (up to --max_rare_count). "
                             "percentile: rarest --min_fraction percentile (old --keep_equal_rare). "
                             "scrub: iterative joint/independent scrub down to --min_fraction (old default).")
    parser.add_argument("--min_kmers", type=int, default=1000,
                        help="count selection: minimum number of kmers before scrub count 1 kmers are added (default 1000)")
    parser.add_argument("--max_rare_count", type=int, default=1,
                        help="count selection: highest scrub count that may be added to reach --min_kmers (default 1)")
    parser.add_argument("--keep_equal_rare", action='store_true',
                        help='alias for --selection percentile (kept for old scripts)')

    # block mode
    parser.add_argument("--reference", "-r", help="Reference genome (fasta[.gz]). Enables block mode: writes "
                        "kmer_loci and blocks tables in addition to the kmer list.")
    parser.add_argument("--max_ref_copies", type=int, default=5,
                        help="Keep kmers present in 1..N copies of the reference (default 5; 1 = single-copy only)")
    parser.add_argument("--terminal_dist", type=int, default=0,
                        help="Drop kmers with any copy within this many bp of a contig end (default 0 = off)")
    parser.add_argument("--min_segment_kmers", type=int, default=1,
                        help="Drop kmers in segments with fewer than N consecutive kmers (default 1 = off). "
                             "A single private SNP gives a segment of k kmers (31 for k=31).")
    parser.add_argument("--max_segment_kmers", type=int, default=0,
                        help="Cut longer runs of consecutive k-mers into segments of N k-mers (default 0 = off). "
                             "Each chunk is validated separately by kmer_block_detect, so one read error only "
                             "rejects the chunk it falls in, and all segments get the same exposure. "
                             "Use k (31) to match single-SNP segments. The remainder chunk is kept as its own "
                             "segment and dropped if it is shorter than --min_segment_kmers/--min_segment_bp.")
    parser.add_argument("--chunk_mode", choices=['fixed', 'balanced'], default='fixed',
                        help="How --max_segment_kmers cuts long runs. fixed: chunks of N, the remainder is its own "
                             "segment (dropped if below the minimum). balanced: floor(n/N) equal chunks of N..2N-1 "
                             "k-mers, nothing of a run >= N is lost (also protects repeat copies, because a dropped "
                             "k-mer is removed at all its copies).")
    parser.add_argument("--min_segment_bp", type=int, default=0,
                        help="Same filter in bp: drop segments shorter than N bp (n_kmers + k - 1). "
                             "A single private SNP spans 2k-1 bp (61 for k=31). Default 0 = off. "
                             "If both are given, the stricter one is used.")

    args = parser.parse_args()
    if args.min_fraction < 0.0 or args.min_fraction > 1.0:
        sys.exit("error: --min_fraction (-m) must be between 0.0 and 1.0 (" + str(args.min_fraction) + ")")
    if not args.scrub_count_file and not args.scrub_count_list:
        sys.exit("error: one of --scrub_count_file or --scrub_count_list must be provided.")
    if args.scrub_count_file and args.scrub_count_list:
        sys.exit("error: can provide only one of either --scrub_count_file or --scrub_count_list.")
    if args.keep_equal_rare:
        args.selection = 'percentile'
    if args.selection == 'percentile' and args.scrub_count_list:
        sys.exit("error: --selection percentile currently supports only --scrub_count_file.")
    if args.max_ref_copies < 1:
        sys.exit("error: --max_ref_copies must be >= 1")

    selection_info = {'selection': args.selection}
    if args.selection == 'count':
        if args.scrub_count_file:
            df_scrub = read_counts_df(args.scrub_count_file)
        else:
            with open(args.scrub_count_list) as fh:
                df_scrub = read_counts_list([l.strip() for l in fh if l.strip()])
        df_scrub, n_non_acgt = drop_non_acgt(df_scrub)
        all_kmers = len(df_scrub)
        df_sel, t = select_by_count(df_scrub, min_kmers=args.min_kmers, max_count=args.max_rare_count,
                                    max_ref_copies=args.max_ref_copies)
        selection_info.update(rare_count=t, min_kmers=args.min_kmers, selected_kmers=len(df_sel))
        del df_scrub

    elif args.selection == 'percentile':
        df_scrub = read_counts_df(args.scrub_count_file)
        df_scrub, n_non_acgt = drop_non_acgt(df_scrub)
        all_kmers = len(df_scrub)
        df_sel = get_lowest_percentile(df_scrub, percentile=args.min_fraction,
                                       drug_scrub='count_hard', max_ref_copies=args.max_ref_copies)
        selection_info.update(min_fraction=args.min_fraction, selected_kmers=len(df_sel))
        del df_scrub

    else:
        scrub_count_files = []
        if args.scrub_count_file:
            scrub_count_files.append(args.scrub_count_file)
        else:
            with open(args.scrub_count_list) as file:
                for line in file:
                    if line.strip():
                        scrub_count_files.append(line.rstrip())

        n_non_acgt = 0         # kmers with ambiguity codes, skipped before selection
        strain_hash = {}       # n times kmer detected in self genome
        metagenome_hash = {}   # n times kmer detected in any scrub metagenome
        pangenome_hash = {}    # n times kmer detected in any scrub genome
        drug_genome_hash = {}  # n times kmer detected in non-self consortia genomes
        drug_filter = 0
        all_kmers = 0
        previous_strain_hash = None

        for i, file in enumerate(scrub_count_files):
            # strain_hash and all_kmers must be identical across files and are not summed
            if i > 0:
                previous_strain_hash = strain_hash
            strain_hash = {}
            all_kmers = 0

            with open_counts(file) as reader:
                for line in reader:
                    if line.startswith('#'):
                        continue
                    content = line.rstrip('\n').split('\t')
                    key = content[0].upper()
                    if not _ACGT.fullmatch(key):
                        if i == 0:
                            n_non_acgt += 1
                        continue
                    all_kmers += 1
                    ref_c, pan_c, meta_c = int(content[1]), int(content[2]), int(content[3])
                    strain_hash[key] = ref_c
                    if pan_c > 0:
                        pangenome_hash[key] = pangenome_hash.get(key, 0) + pan_c
                    if meta_c > 0:
                        metagenome_hash[key] = metagenome_hash.get(key, 0) + meta_c
                    if len(content) == 5:
                        drug_filter = 1
                        drug_c = int(content[4])
                        if drug_c > 0:
                            drug_genome_hash[key] = drug_genome_hash.get(key, 0) + drug_c

            if i > 0 and strain_hash != previous_strain_hash:
                sys.exit("error: input files do not have identical hash and strain hash values.")

        print("#total kmers in strain:" + str(all_kmers) + "," + str(len(strain_hash)) + " pangenome: "
              + str(len(pangenome_hash)) + " metagenome: " + str(len(metagenome_hash)))

        # raw counts, captured before the scrub functions mutate the hashes
        ref_counts = dict(strain_hash)
        raw_pan = dict(pangenome_hash)
        raw_meta = dict(metagenome_hash)
        raw_drug = dict(drug_genome_hash)

        drug_scrubbed = 0
        if drug_filter:
            print("#total kmers cross drug:" + str(len(drug_genome_hash)))
            strain_hash = drug_scrub(drug_genome_hash, strain_hash, all_kmers)
            fraction_remaining = float(len(strain_hash)) / float(all_kmers)
            drug_scrubbed = all_kmers - len(strain_hash)
            print("#fraction kmers remaining drug post scrub:" + str(fraction_remaining))
            print("#drug_scrubbed kmers:" + str(drug_scrubbed))
            if fraction_remaining < args.min_fraction * 2:
                raise Exception('ERROR: too few kmers remain after drug scrub. Are your drug strains too similar?')

        if args.independent:
            strain_hash_scrubbed = independent_scrub(args.min_fraction, pangenome_hash, metagenome_hash,
                                                     strain_hash, all_kmers)
        else:
            strain_hash_scrubbed = joint_scrub(args.min_fraction, pangenome_hash, metagenome_hash,
                                               strain_hash, all_kmers, drug_scrubbed)

        sel = list(strain_hash_scrubbed)
        df_sel = pd.DataFrame({
            '#kmer': sel,
            'reference_count': [ref_counts[k] for k in sel],
            'pangenome_count': [raw_pan.get(k, 0) for k in sel],
            'metagenome_count': [raw_meta.get(k, 0) for k in sel],
            'drug_count': [raw_drug.get(k, 0) for k in sel],
        })
        n_before = len(df_sel)
        df_sel = df_sel[df_sel['reference_count'].between(1, args.max_ref_copies)].copy()
        print(f"#removed {n_before - len(df_sel)} kmers with reference_count > {args.max_ref_copies}")
        del ref_counts, raw_pan, raw_meta, raw_drug

    selection_info['non_acgt_kmers'] = n_non_acgt          # removed before selection
    if len(df_sel) == 0:
        sys.exit("error: no kmers selected")

    os.makedirs(args.output_dir, exist_ok=True)
    prefix = f"{args.basename}." if args.basename else ""

    if args.reference:
        if len(df_sel) == 0:
            sys.exit("error: no kmers selected, nothing to place on the reference")
        k = len(df_sel['#kmer'].iloc[0])
        min_seg = max(args.min_segment_kmers, args.min_segment_bp - k + 1, 1)
        if 0 < args.max_segment_kmers < min_seg:
            sys.exit(f"error: --max_segment_kmers ({args.max_segment_kmers}) is smaller than the minimum "
                     f"segment size ({min_seg} kmers): every chunk would be dropped")
        print(f"#min segment size: {min_seg} kmers ({min_seg + k - 1} bp)"
              + (f", segments cut into {args.chunk_mode} chunks of {args.max_segment_kmers} kmers"
                 if args.max_segment_kmers > 0 else ""))
        records = load_genome(args.reference)
        df_loci, df_blocks, stats = build_blocks(
            df_sel, records, k,
            terminal_dist=args.terminal_dist,
            max_ref_copies=args.max_ref_copies,
            min_segment_kmers=min_seg,
            max_segment_kmers=args.max_segment_kmers,
            chunk_mode=args.chunk_mode,
        )
        for key, val in stats.items():
            print(f"#block_mode {key}: {val}")
        if stats['kept_kmers'] == 0:
            sys.exit("error: no kmers left after block filtering (check --min_segment_kmers/--min_segment_bp)")

        df_loci.to_csv(os.path.join(args.output_dir, prefix + "kmer_loci.tsv.gz"),
                       sep='\t', index=False, compression='gzip')
        df_blocks.to_csv(os.path.join(args.output_dir, prefix + "blocks.tsv.gz"),
                         sep='\t', index=False, compression='gzip')
        kept = pd.unique(df_loci['#kmer'])
    else:
        kept = df_sel['#kmer'].to_numpy()

    print("#post scrub kmers " + str(len(kept)) + " out of " + str(all_kmers))
    if args.selection == 'scrub':
        selection_info.update(min_fraction=args.min_fraction)
    selection_info.update(total_strain_kmers=int(all_kmers), kept_kmers=int(len(kept)))
    if args.reference:
        selection_info.update(max_ref_copies=args.max_ref_copies, min_segment_kmers=min_seg,
                              max_segment_kmers=args.max_segment_kmers, chunk_mode=args.chunk_mode,
                              **{f'filter_{key}': (int(v) if isinstance(v, (np.integer,)) else
                                                   float(v) if isinstance(v, np.floating) else v)
                                 for key, v in stats.items()})
    if selection_info.get('kept_kmers', 0) < args.min_kmers and args.selection == 'count':
        print(f"#WARNING: {len(kept)} kmers left after block filtering (< --min_kmers {args.min_kmers})")
    with open(os.path.join(args.output_dir, prefix + "selection.json"), 'w') as fh:
        json.dump(selection_info, fh, indent=2)

    with gzip.open(os.path.join(args.output_dir, prefix + "scrubbed_kmers.gz"), 'wt') as out:
        out.write("#kmer\n")
        for key in kept:
            out.write(key + "\n")


if __name__ == "__main__":
    main()
