#!/usr/bin/env python3
"""
scrubdb_hf.py - package, upload and restore the scrub database via a Hugging Face dataset repo.

Subcommands
  pack     Stage everything into sharded archives + a manifest (this folder IS the repo layout)
  upload   Create the dataset repo if needed and upload the staging folder (resumable)
  pull     Download from HF, then run `restore`
  restore  Verify checksums and unpack into --dest

Repo / staging layout
  MANIFEST.json                               layout used by `restore`
  README.md                                   dataset card (written once, edit freely)
  manifest/checksums.sha256                   sha256sum-compatible, paths relative to repo root
  manifest/<section>.tsv                      name, size, shard (genomes, metagenomes)
  genomes/genomes.NNN.tar                     all files of one genome folder, sharded
  metagenomes/metagenomes.NNN.tar.zst         all *.kmdb entries of one folder, sharded
  rocksdb/<name>.tar.zst.part-NNN             sourmash RocksDB, compressed and split

Restored layout (--dest)
  genomes/  metagenomes/  <rocksdb name>/
  genomes_to_scrub_full.txt    absolute path of every entry in genomes/, one per line
  metagenomes_to_scrub.txt     absolute path of every entry in metagenomes/, one per line

The RocksDB is archived by its contents (not its directory name), so --rocksdb-name only decides
the directory name it gets on restore.
"""
import argparse
import datetime
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

FORMAT_VERSION = 4
SECTIONS = ("genomes", "metagenomes")
LIST_FILES = {"genomes": "genomes_to_scrub_full.txt", "metagenomes": "metagenomes_to_scrub.txt"}
ROCKSDB_EXCLUDES = ("LOCK", "LOG", "LOG.old.*")
COMPRESSED_EXTS = (".gz", ".bz2", ".xz", ".zst")

CODECS = {
    "none": {"ext": ".tar", "comp": "cat", "decomp": "cat", "level": 0},
    "zstd": {"ext": ".tar.zst", "comp": "zstd -q -T{threads} -{level}", "decomp": "zstd -dc -q", "level": 3},
    "pigz": {"ext": ".tar.gz", "comp": "pigz -p {threads} -{level}", "decomp": "pigz -dc", "level": 6},
    "gzip": {"ext": ".tar.gz", "comp": "gzip -{level}", "decomp": "gzip -dc", "level": 6},
}

q = shlex.quote


# ----------------------------------------------------------------------------- helpers

def log(msg):
    print(f"[{datetime.datetime.now():%H:%M:%S}] {msg}", file=sys.stderr, flush=True)


def die(msg):
    log(f"ERROR: {msg}")
    sys.exit(1)


def run(cmd):
    """Run a shell pipeline; fail if any stage fails."""
    subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True)


def sha256(path, bufsize=16 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(bufsize), b""):
            h.update(chunk)
    return h.hexdigest()


def pick_codec(requested, fallback_order=("zstd", "pigz", "gzip")):
    order = [requested] if requested else list(fallback_order)
    for name in order:
        if shutil.which(CODECS[name]["comp"].split()[0]):
            return name
    die(f"no compressor found (tried {order}); `conda install -c conda-forge zstd`")


def entry_size(p):
    """Size of a file, or total size of all files under a directory (symlinks followed)."""
    p = Path(p)
    if p.is_dir():
        return sum(f.stat().st_size for f in p.rglob("*") if f.is_file())
    return p.stat().st_size


def write_tsv(path, header, rows):
    with open(path, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")


def read_tsv(path):
    with open(path) as f:
        f.readline()
        return [line.rstrip("\n").split("\t") for line in f if line.strip()]


def make_shards(sizes, target_bytes):
    shards, cur, cur_size = [], [], 0
    for n in sorted(sizes):
        if cur and cur_size + sizes[n] > target_bytes:
            shards.append(cur)
            cur, cur_size = [], 0
        cur.append(n)
        cur_size += sizes[n]
    if cur:
        shards.append(cur)
    return shards


# ----------------------------------------------------------------------------- pack

def list_genomes(d):
    """All files (incl. symlinks to files) directly in d; subdirectories are skipped."""
    d = Path(d)
    if not d.is_dir():
        die(f"missing directory: {d}")
    names = []
    for e in os.scandir(d):
        if e.is_dir(follow_symlinks=True):
            log(f"  skipping subdirectory {e.path}")
            continue
        if not os.path.exists(e.path):
            die(f"broken symlink: {e.path}")
        names.append(e.name)
    if not names:
        die(f"no files in {d}")
    return names


def list_metagenomes(d, pattern):
    """All entries in d matching pattern (files or directories)."""
    d = Path(d)
    if not d.is_dir():
        die(f"missing directory: {d}")
    names = sorted(p.name for p in d.glob(pattern))
    if not names:
        die(f"no '{pattern}' entries in {d}")
    broken = [n for n in names if not (d / n).exists()]
    if broken:
        die(f"{len(broken)} broken symlinks in {d}, e.g. {broken[:3]}")
    return names


def pack_shard(src_dir, names, out, staging, codec, threads):
    """tar (dereferencing symlinks) + optional compression of one shard. Resumable via a state file."""
    rel = out.relative_to(staging)
    state = staging / ".state" / (str(rel) + ".json")
    list_hash = hashlib.sha256(("\0".join(names) + "\0" + str(src_dir)).encode()).hexdigest()
    if out.exists() and state.exists():
        st = json.loads(state.read_text())
        if st.get("list_hash") == list_hash:
            return str(rel), st["sha256"], "reused"

    state.parent.mkdir(parents=True, exist_ok=True)
    filelist = Path(str(out) + ".filelist")
    tmp = Path(str(out) + ".tmp")
    filelist.write_bytes(b"\0".join(n.encode() for n in names) + b"\0")
    comp = codec["comp"].format(threads=threads, level=codec["level"])
    run(f"tar -c -h -f - -C {q(str(src_dir))} --null -T {q(str(filelist))} | {comp} > {q(str(tmp))}")
    tmp.rename(out)
    filelist.unlink()
    digest = sha256(out)
    state.write_text(json.dumps({"list_hash": list_hash, "sha256": digest}))
    return str(rel), digest, "built"


def pack_section(section, src_dir, names, staging, codec_name, shard_bytes, threads, jobs):
    codec = CODECS[codec_name]
    sizes = {n: entry_size(src_dir / n) for n in names}
    shards = make_shards(sizes, shard_bytes)
    log(f"  {section}: {len(names)} entries, {sum(sizes.values()) / 1024 ** 3:.1f} GiB, "
        f"{len(shards)} shard(s), codec={codec_name}")

    outdir = staging / section
    outdir.mkdir(parents=True, exist_ok=True)
    paths = [outdir / f"{section}.{i:03d}{codec['ext']}" for i in range(len(shards))]
    for old in outdir.iterdir():          # drop shards left over from an earlier, different layout
        if old not in paths:
            old.unlink()

    threads_per_job = max(1, threads // jobs)
    checksums = {}
    with ThreadPoolExecutor(max_workers=jobs) as ex:
        futs = [ex.submit(pack_shard, src_dir, n, o, staging, codec, threads_per_job)
                for n, o in zip(shards, paths)]
        for k, f in enumerate(as_completed(futs), 1):
            rel, digest, status = f.result()
            checksums[rel] = digest
            log(f"    [{k}/{len(paths)}] {rel} ({status})")

    rows = [(n, sizes[n], paths[i].name) for i, ns in enumerate(shards) for n in ns]
    write_tsv(staging / "manifest" / f"{section}.tsv", ["name", "size", "shard"], rows)
    info = {"codec": codec_name, "count": len(names),
            "shards": [str(p.relative_to(staging)) for p in paths]}
    return info, checksums


def pack_rocksdb(rdb, name, staging, codec_name, threads, part_bytes):
    """Archive the *contents* of the RocksDB dir, so the restored dir name is free to choose."""
    codec = CODECS[codec_name]
    outdir = staging / "rocksdb"
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = outdir / f"{name}{codec['ext']}.part-"
    state = staging / ".state" / "rocksdb.json"
    if state.exists():
        st = json.loads(state.read_text())
        if (st.get("source") == str(rdb) and st.get("name") == name and st.get("codec") == codec_name
                and all((staging / p).exists() for p in st["parts"])):
            log("  rocksdb: reusing existing parts")
            return st
    for old in outdir.iterdir():
        old.unlink()
    excl = " ".join(f"--exclude={q(x)}" for x in ROCKSDB_EXCLUDES)
    comp = codec["comp"].format(threads=threads, level=codec["level"])
    run(f"tar -c -f - -C {q(str(rdb))} {excl} . | {comp} "
        f"| split -a 3 -d -b {part_bytes} - {q(str(prefix))}")
    parts = sorted(outdir.glob(f"{name}{codec['ext']}.part-*"))
    st = {
        "source": str(rdb),     # kept in the local .state/ folder only (never uploaded)
        "name": name,
        "codec": codec_name,
        "parts": [str(p.relative_to(staging)) for p in parts],
        "sha256": {str(p.relative_to(staging)): sha256(p) for p in parts},
    }
    state.parent.mkdir(parents=True, exist_ok=True)
    state.write_text(json.dumps(st))
    log(f"  rocksdb: {len(parts)} part(s), restored as '{name}/'")
    return st


def cmd_pack(a):
    staging = Path(a.staging).resolve()
    shard_bytes = int(a.shard_gb * 1024 ** 3)
    srcs = {
        "genomes": Path(a.genomes_dir).resolve() if a.genomes_dir else None,
        "metagenomes": Path(a.metagenomes_dir).resolve() if a.metagenomes_dir else None,
    }
    rdb = Path(a.rocksdb).resolve() if a.rocksdb else None
    if not any([*srcs.values(), rdb]):
        die("nothing to pack: give at least one of --genomes-dir, --metagenomes-dir, --rocksdb")
    for src in filter(None, [*srcs.values(), rdb]):
        if staging == src or src in staging.parents or staging in src.parents:
            die(f"--staging {staging} overlaps the source {src}; use a separate folder (e.g. scratch)")
    if rdb and not (rdb / "CURRENT").exists():
        die(f"{rdb} does not look like a RocksDB directory (no CURRENT file)")
    rdb_name = a.rocksdb_name or (rdb.name if rdb else None)

    names = {}
    if srcs["genomes"]:
        names["genomes"] = list_genomes(srcs["genomes"])
    if srcs["metagenomes"]:
        names["metagenomes"] = list_metagenomes(srcs["metagenomes"], a.metagenome_pattern)

    codecs = {}
    if "genomes" in names:     # genomes are normally already gzipped -> plain tar
        g = names["genomes"]
        n_comp = sum(n.endswith(COMPRESSED_EXTS) for n in g)
        codecs["genomes"] = pick_codec(a.genome_codec) if a.genome_codec else (
            "none" if n_comp > len(g) / 2 else pick_codec(None))
    if "metagenomes" in names:
        codecs["metagenomes"] = pick_codec(a.metagenome_codec)
    rdb_codec = pick_codec(a.rocksdb_codec) if rdb else None

    if a.dry_run:
        for s, ns in names.items():
            log(f"  {s}: {len(ns)} entries from {srcs[s]}, codec {codecs[s]}, e.g. {ns[:3]}")
        if rdb:
            log(f"  rocksdb: {rdb} -> '{rdb_name}' (codec {rdb_codec})")
        log("Dry run - nothing written.")
        return

    (staging / "manifest").mkdir(parents=True, exist_ok=True)
    checksums, sections = {}, {}
    for s, ns in names.items():
        sections[s], cs = pack_section(s, srcs[s], ns, staging, codecs[s], shard_bytes, a.threads, a.jobs)
        checksums.update(cs)

    rdb_info = None
    if rdb:
        log(f"Packing RocksDB {rdb} with {rdb_codec} (make sure nothing is writing to it)")
        st = pack_rocksdb(rdb, rdb_name, staging, rdb_codec, a.threads, int(a.part_gb * 1024 ** 3))
        checksums.update(st["sha256"])
        rdb_info = {"name": st["name"], "codec": st["codec"], "parts": st["parts"]}

    with open(staging / "manifest" / "checksums.sha256", "w") as f:
        for rel in sorted(checksums):
            f.write(f"{checksums[rel]}  {rel}\n")

    manifest = {
        "format_version": FORMAT_VERSION,
        "created": datetime.datetime.now().isoformat(timespec="seconds"),
        "sections": sections,
        "rocksdb": rdb_info,
    }
    (staging / "MANIFEST.json").write_text(json.dumps(manifest, indent=2))

    readme = staging / "README.md"
    if not readme.exists():
        readme.write_text(README_TEMPLATE)
    log(f"Staging complete: {staging}")


README_TEMPLATE = """---
license: other
tags: [metagenomics, bacterial-genomes, strain-detection]
---

# Scrub database

- `genomes/`: reference genomes (sharded tar archives of .fna.gz files)
- `metagenomes/`: k-mer databases (*.kmdb), sharded archives
- `rocksdb/`: sourmash RocksDB index (split archive)

Restore with `python scrubdb_hf.py pull --repo <user>/<repo> --download-dir <dl> --dest <dir>`.
See `MANIFEST.json` for the machine-readable layout.
"""


# ----------------------------------------------------------------------------- upload

def cmd_upload(a):
    from huggingface_hub import HfApi

    staging = Path(a.staging).resolve()
    if not (staging / "MANIFEST.json").exists():
        die(f"{staging}/MANIFEST.json not found - run `pack` first")
    api = HfApi()
    log(f"Logged in as: {api.whoami()['name']}")
    url = api.create_repo(a.repo, repo_type="dataset", private=not a.public, exist_ok=True)
    log(f"Repo: {url}")
    api.upload_large_folder(
        repo_id=a.repo,
        folder_path=staging,
        repo_type="dataset",
        ignore_patterns=[".state/*", "*.tmp", "*.filelist"],
        num_workers=a.workers,
    )
    if a.tag:
        api.create_tag(a.repo, tag=a.tag, repo_type="dataset", exist_ok=True)
        log(f"Tagged revision as {a.tag}")
    log("Upload finished.")


# ----------------------------------------------------------------------------- pull / restore

def wanted(a):
    """Which parts the user asked for."""
    return {
        "genomes": not a.no_genomes,
        "metagenomes": not a.no_metagenomes,
        "rocksdb": not a.no_rocksdb,
    }


def cmd_pull(a):
    from huggingface_hub import snapshot_download

    w = wanted(a)
    allow = ["MANIFEST.json", "README.md", "manifest/*"] + [f"{k}/*" for k, v in w.items() if v]
    log(f"Downloading {a.repo}@{a.revision or 'main'} -> {a.download_dir}")
    snapshot_download(
        repo_id=a.repo,
        repo_type="dataset",
        revision=a.revision,
        local_dir=a.download_dir,
        allow_patterns=allow,
        max_workers=a.workers,
    )
    cmd_restore(a)


def cmd_restore(a):
    dl = Path(a.download_dir).resolve()
    dest = Path(a.dest).resolve()
    man = json.loads((dl / "MANIFEST.json").read_text())
    if man["format_version"] != FORMAT_VERSION:
        die(f"unsupported format_version {man['format_version']} (script expects {FORMAT_VERSION})")
    w = wanted(a)
    sections = {s: info for s, info in man["sections"].items() if w[s]}
    rdb = man["rocksdb"] if w["rocksdb"] else None
    rdb_name = a.rocksdb_name or (rdb["name"] if rdb else None)

    marks = dest / ".restore_state"          # tied to the destination, not the download dir
    marks.mkdir(parents=True, exist_ok=True)

    def mark_of(rel):
        return marks / (rel.replace("/", "__") + ".ok")

    for cname in {i["codec"] for i in sections.values()} | ({rdb["codec"]} if rdb else set()):
        tool = CODECS[cname]["decomp"].split()[0]
        if not shutil.which(tool):
            die(f"`{tool}` is required to unpack this repo")

    present = [p for i in sections.values() for p in i["shards"]]
    present += rdb["parts"] if rdb else []
    missing = [p for p in present if not (dl / p).exists()]
    if missing:
        die(f"{len(missing)} file(s) missing from {dl}, e.g. {missing[:3]}")

    # --- checksums (skipping archives already unpacked in a previous run)
    if not a.skip_verify and present:
        sums = {}
        for line in (dl / "manifest" / "checksums.sha256").read_text().splitlines():
            h, rel = line.split(None, 1)
            sums[rel.strip()] = h
        todo = [p for i in sections.values() for p in i["shards"] if not mark_of(p).exists()]
        if rdb and not (marks / "rocksdb.ok").exists():
            todo += rdb["parts"]
        log(f"Verifying sha256 of {len(todo)} file(s)")
        with ThreadPoolExecutor(max_workers=a.jobs) as ex:
            for rel, h in zip(todo, ex.map(lambda p: sha256(dl / p), todo)):
                if h != sums.get(rel):
                    die(f"checksum mismatch: {rel}")
        log("  all checksums OK")

    # --- sharded sections
    for s, info in sections.items():
        decomp = CODECS[info["codec"]]["decomp"]
        target = dest / s
        target.mkdir(parents=True, exist_ok=True)

        def unpack(rel, decomp=decomp, target=target):
            if mark_of(rel).exists():
                return rel, "already unpacked"
            run(f"{decomp} {q(str(dl / rel))} | tar -x --no-same-owner -f - -C {q(str(target))}")
            mark_of(rel).touch()
            return rel, "unpacked"

        log(f"Unpacking {len(info['shards'])} {s} shard(s) into {target}")
        with ThreadPoolExecutor(max_workers=a.jobs) as ex:
            futs = [ex.submit(unpack, r) for r in info["shards"]]
            for k, f in enumerate(as_completed(futs), 1):
                rel, status = f.result()
                log(f"  [{k}/{len(info['shards'])}] {rel} ({status})")

    # --- rocksdb
    if rdb:
        rtarget = (Path(a.rocksdb_dest).resolve() if a.rocksdb_dest else dest) / rdb_name
        mark = marks / "rocksdb.ok"
        if mark.exists() and rtarget.exists():
            log(f"RocksDB already unpacked at {rtarget}")
        else:
            if rtarget.exists():
                if not a.force:
                    die(f"{rtarget} exists; use --force to overwrite")
                shutil.rmtree(rtarget)
            rtarget.mkdir(parents=True)
            parts = " ".join(q(str(dl / p)) for p in rdb["parts"])
            log(f"Unpacking RocksDB -> {rtarget}")
            run(f"cat {parts} | {CODECS[rdb['codec']]['decomp']} | tar -x --no-same-owner -f - -C {q(str(rtarget))}")
            mark.touch()

    # --- verify sizes of unpacked entries
    errors = 0
    for s in sections:
        for name, size, _ in read_tsv(dl / "manifest" / f"{s}.tsv"):
            p = dest / s / name
            if not p.exists() or entry_size(p) != int(size):
                errors += 1
                if errors <= 5:
                    log(f"  BAD: {p}")
    if errors:
        die(f"verification failed: {errors} problem(s)")
    log("Verified: " + ", ".join(f"{s}={i['count']}" for s, i in sections.items())
        + (f", rocksdb={rdb_name}" if rdb else ""))

    # --- path lists for the scrub step (absolute paths of the restored entries)
    for s in sections:
        out = dest / LIST_FILES[s]
        names = sorted(r[0] for r in read_tsv(dl / "manifest" / f"{s}.tsv"))
        out.write_text("".join(f"{dest / s / n}\n" for n in names))
        log(f"Wrote {out.name} ({len(names)} paths)")

    if a.cleanup:
        for p in present:
            (dl / p).unlink(missing_ok=True)
        log("Removed downloaded archives (--cleanup)")
    log(f"Done: {dest}")


# ----------------------------------------------------------------------------- CLI

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("pack", help="stage everything into archives + manifest")
    p.add_argument("--genomes-dir", help="folder whose files go into genomes/")
    p.add_argument("--metagenomes-dir", help="folder whose *.kmdb entries go into metagenomes/")
    p.add_argument("--metagenome-pattern", default="*.kmdb", help="glob for --metagenomes-dir (default *.kmdb)")
    p.add_argument("--rocksdb", help="sourmash RocksDB directory to include")
    p.add_argument("--rocksdb-name", help="directory name the RocksDB gets on restore (default: its current name)")
    p.add_argument("--staging", required=True, help="output folder (= repo layout)")
    p.add_argument("--shard-gb", type=float, default=10, help="GiB per shard (default 10)")
    p.add_argument("--part-gb", type=float, default=10, help="GiB per RocksDB archive part (default 10)")
    p.add_argument("--genome-codec", choices=list(CODECS),
                   help="default: none if genomes are already compressed, else zstd")
    p.add_argument("--metagenome-codec", choices=list(CODECS), help="default: zstd, else pigz, else gzip")
    p.add_argument("--rocksdb-codec", choices=list(CODECS), help="default: zstd, else pigz, else gzip")
    p.add_argument("--threads", type=int, default=os.cpu_count() or 4)
    p.add_argument("--jobs", type=int, default=4, help="shards built in parallel")
    p.add_argument("--dry-run", action="store_true")
    p.set_defaults(func=cmd_pack)

    p = sub.add_parser("upload", help="upload staging folder to a HF dataset repo")
    p.add_argument("--staging", required=True)
    p.add_argument("--repo", required=True, help="e.g. yourname/scrubdb-derep-hq")
    p.add_argument("--public", action="store_true", help="create repo as public (default private)")
    p.add_argument("--tag", help="tag this upload, e.g. v2026.10")
    p.add_argument("--workers", type=int, default=8)
    p.set_defaults(func=cmd_upload)

    def add_restore_args(p):
        p.add_argument("--download-dir", required=True, help="where archives are / will be downloaded")
        p.add_argument("--dest", required=True, help="directory to restore into")
        p.add_argument("--rocksdb-dest", help="parent dir for the RocksDB (default: --dest)")
        p.add_argument("--rocksdb-name", help="override the RocksDB directory name on restore")
        p.add_argument("--jobs", type=int, default=4)
        p.add_argument("--skip-verify", action="store_true", help="skip sha256 of downloaded files")
        p.add_argument("--cleanup", action="store_true", help="delete archives after successful restore")
        p.add_argument("--force", action="store_true", help="overwrite existing RocksDB")
        p.add_argument("--no-genomes", action="store_true")
        p.add_argument("--no-metagenomes", action="store_true")
        p.add_argument("--no-rocksdb", action="store_true")

    p = sub.add_parser("pull", help="download from HF and restore")
    p.add_argument("--repo", required=True)
    p.add_argument("--revision", help="branch, tag or commit (default main)")
    p.add_argument("--workers", type=int, default=8)
    add_restore_args(p)
    p.set_defaults(func=cmd_pull)

    p = sub.add_parser("restore", help="restore from an already-downloaded folder")
    add_restore_args(p)
    p.set_defaults(func=cmd_restore)

    a = ap.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()