"""Canonical form of kmer_scrub_count outputs for multi-threaded runs, where the
original assigns scrub_ids in sample-completion order: summary rows without the
id (sorted), presence rows with ids replaced by type:sample_id (sorted)."""
import sys, io, subprocess, hashlib
pre = sys.argv[1]
name = {}
rows = []
for i, line in enumerate(open(pre + '.summary.tsv')):
    f = line.rstrip('\n').split('\t')
    if i == 0: rows.append('\t'.join(f[1:])); continue
    name[f[0]] = f[1] + ':' + f[2]
    rows.append('\t'.join(f[1:]))
out = [rows[0]] + sorted(rows[1:])
try:
    p = subprocess.run(['zstdcat', pre + '.presence.tsv.zst'], capture_output=True, check=True, text=True).stdout
    h = hashlib.sha256()
    pres = []
    for line in p.splitlines()[1:]:
        k, ids = line.split('\t')
        pres.append(k + '\t' + ','.join(sorted(name[x] for x in ids.split(','))))
    pres.sort()
    for l in pres: h.update(l.encode() + b'\n')
    out.append('presence ' + str(len(pres)) + ' rows sha256 ' + h.hexdigest())
except (FileNotFoundError, subprocess.CalledProcessError):
    out.append('presence none')
print('\n'.join(out))
