#!/usr/bin/env bash
# Tests for kmer_metagenome_db --pair-dedup (sample mode).
#   bash test/pair_dedup/run_tests.sh            (from the repo root or anywhere)
# Optional: BASE=/path/to/old/kmer_metagenome_db also checks that single-end
# and --no-pair-dedup output is byte-for-byte the k-mer counts of that binary.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
SRC="$HERE/../../src"
BIN="${BIN:-$SRC/kmer_metagenome_db}"
OUT="${OUT:-$HERE/out}"
rm -rf "$OUT"; mkdir -p "$OUT"
CC="${CC:-gcc}"
$CC -O2 -I"$SRC" -o "$OUT/kmdb_dump" "$HERE/kmdb_dump.c" "$SRC/kmdb.c" -lzstd || exit 1
DUMP="$OUT/kmdb_dump"
python3 "$HERE/make_cases.py" "$OUT/cases" 31 || exit 1
C="$OUT/cases"

pass=0; fail=0
ok()   { echo "  ok    $1"; pass=$((pass+1)); }
bad()  { echo "  FAIL  $1"; fail=$((fail+1)); }
same() { if cmp -s "$2" "$3"; then ok "$1"; else bad "$1"; diff "$2" "$3" | head -5; fi; }

db() {  # db <name> <args...> -> $OUT/<name>.tsv (dump), stderr in .log
	local n=$1; shift
	"$BIN" -m sample -o "$OUT/$n.kmdb" "$@" 2> "$OUT/$n.log" || { echo "  (build $n failed)"; cat "$OUT/$n.log"; return 1; }
	"$DUMP" "$OUT/$n.kmdb" > "$OUT/$n.tsv"
}

echo "== expected counts (independent Python reference)"
for c in full none partial repeat bulk; do
	db "$c.pe"     -t 3 -i "$C/${c}_R1.fq,$C/${c}_R2.fq"
	db "$c.pe_off" -t 3 --no-pair-dedup -i "$C/${c}_R1.fq,$C/${c}_R2.fq"
	same "$c: PE dedup == reference"            "$OUT/$c.pe.tsv"     "$C/$c.dedup.tsv"
	same "$c: PE --no-pair-dedup == reference"  "$OUT/$c.pe_off.tsv" "$C/$c.nodedup.tsv"
done

echo "== specific properties"
awk '$2!=1{e=1} END{exit e}' "$OUT/full.pe.tsv"     && ok "full overlap: every k-mer 1 with dedup"    || bad "full overlap dedup"
awk '$2!=2{e=1} END{exit e}' "$OUT/full.pe_off.tsv" && ok "full overlap: every k-mer 2 without"       || bad "full overlap no-dedup"
same "no overlap: dedup == no-dedup" "$OUT/none.pe.tsv" "$OUT/none.pe_off.tsv"
n_diff=$(paste "$OUT/partial.pe.tsv" "$OUT/partial.pe_off.tsv" | awk '$2!=$4' | wc -l)
[ "$n_diff" -eq 20 ] && ok "partial overlap: exactly the 20 overlap k-mers differ (50 bp - 31 + 1)" \
                     || bad "partial overlap: $n_diff k-mers differ, expected 20"
r_max=$(sort -k2,2n "$OUT/repeat.pe.tsv" | tail -1 | cut -f2)
r_off=$(sort -k2,2n "$OUT/repeat.pe_off.tsv" | tail -1 | cut -f2)
[ "$r_max" = 1 ] && [ "$r_off" = 2 ] && ok "internal repeat: 1 with dedup (2 without)" || bad "internal repeat ($r_max / $r_off)"

echo "== single-end unchanged"
db se_pos    -t 3 "$C/bulk_R1.fq" "$C/bulk_R2.fq"
db se_i      -t 3 -i "$C/bulk_R1.fq" -i "$C/bulk_R2.fq"
same "SE (two positional files) == no-dedup reference" "$OUT/se_pos.tsv" "$C/bulk.nodedup.tsv"
same "SE (two -i) == no-dedup reference"               "$OUT/se_i.tsv"   "$C/bulk.nodedup.tsv"
grep -q "WARNING: pair dedup is on but no read pairs" "$OUT/se_pos.log" && ok "SE with 2 files warns" || bad "SE warning missing"
if [ -n "${BASE:-}" ]; then
	"$BASE" -m sample -t 3 -o "$OUT/base_se.kmdb" "$C/bulk_R1.fq" "$C/bulk_R2.fq" 2>/dev/null
	"$BASE" -m sample -t 3 -o "$OUT/base_comma.kmdb" -i "$C/bulk_R1.fq,$C/bulk_R2.fq" 2>/dev/null
	"$DUMP" "$OUT/base_se.kmdb" > "$OUT/base_se.tsv"; "$DUMP" "$OUT/base_comma.kmdb" > "$OUT/base_comma.tsv"
	same "SE == old binary"                     "$OUT/se_pos.tsv"      "$OUT/base_se.tsv"
	same "PE --no-pair-dedup == old binary"     "$OUT/bulk.pe_off.tsv" "$OUT/base_comma.tsv"
	for key in n_reads kmer_positions valid_kmer_positions; do
		a=$("$DUMP" -m $key "$OUT/bulk.pe.kmdb"); b=$("$DUMP" -m $key "$OUT/base_se.kmdb")
		[ "$a" = "$b" ] && ok "metadata $key unchanged by dedup ($a)" || bad "metadata $key $a vs old $b"
	done
fi

echo "== PEI == PE, threads, low-memory passes"
db bulk.pei  -t 3 --interleaved -i "$C/bulk_I.fq"
same "PEI == PE" "$OUT/bulk.pei.tsv" "$OUT/bulk.pe.tsv"
db bulk.t1  -t 1 -i "$C/bulk_R1.fq,$C/bulk_R2.fq"
db bulk.t8  -t 8 -i "$C/bulk_R1.fq,$C/bulk_R2.fq"
same "1 thread == 3 threads" "$OUT/bulk.t1.tsv" "$OUT/bulk.pe.tsv"
same "8 threads == 3 threads" "$OUT/bulk.t8.tsv" "$OUT/bulk.pe.tsv"
db bulk.M3  -t 3 -M 3 -i "$C/bulk_R1.fq,$C/bulk_R2.fq"
same "-M 3 passes == single pass" "$OUT/bulk.M3.tsv" "$OUT/bulk.pe.tsv"
db bulk.R   -t 2 -R 1 -i "$C/bulk_R1.fq,$C/bulk_R2.fq"
same "-R (estimate pass + presized) == single pass" "$OUT/bulk.R.tsv" "$OUT/bulk.pe.tsv"
db bulk.lanes -t 3 -i "$C/full_R1.fq,$C/full_R2.fq" -i "$C/bulk_R1.fq,$C/bulk_R2.fq"
python3 - "$OUT/bulk.lanes.tsv" "$C/full.dedup.tsv" "$C/bulk.dedup.tsv" <<'EOF' && ok "two -i pairs = sum of both" || bad "two -i pairs"
import sys
from collections import Counter
def rd(p): return Counter({l.split()[0]: int(l.split()[1]) for l in open(p)})
sys.exit(0 if rd(sys.argv[1]) == rd(sys.argv[2]) + rd(sys.argv[3]) else 1)
EOF
grep -q "^pair dedup ON:.*increments skipped" "$OUT/bulk.pe.log" && ok "stderr reports skipped increments" || bad "stderr stats line"
"$BIN" -m sample -o "$OUT/st.kmdb" --stats "$OUT/st.tsv" -i "$C/bulk_R1.fq,$C/bulk_R2.fq" 2>/dev/null
python3 - "$OUT/st.tsv" "$C/bulk.nodedup.tsv" "$C/bulk.dedup.tsv" <<'EOF' && ok "--stats: skipped = no-dedup sum - dedup sum" || bad "--stats numbers"
import sys
row = dict(zip(*[l.rstrip("\n").split("\t") for l in open(sys.argv[1])]))
s = lambda p: sum(int(l.split()[1]) for l in open(p))
nod, ded = s(sys.argv[2]), s(sys.argv[3])
ok = int(row["dedup_skipped"]) == nod - ded and int(row["valid_kmer_positions"]) == nod \
     and int(row["count_sum"]) == ded
print(f"    hits {nod}, skipped {row['dedup_skipped']} ({float(row['dedup_frac_all']):.2%})")
sys.exit(0 if ok else 1)
EOF

echo "== total_evaluated = sum of all counts (kmer_metagenome_db)"
for n in bulk.pe bulk.pe_off se_pos bulk.M3; do
	te=$("$DUMP" -m total_evaluated "$OUT/$n.kmdb")
	cs=$(awk '{s+=$2} END{print s}' "$OUT/$n.tsv")
	[ "$te" = "$cs" ] && ok "$n: total_evaluated == sum of counts ($te)" || bad "$n: total_evaluated $te vs sum $cs"
done

echo "== kmer_strain_detect"
SD="${BIN_SD:-$SRC/kmer_strain_detect}"
# query = every k-mer of the bulk sample, so the count columns hold the whole
# spectrum and total_evaluated must equal each column's sum
cut -f1 "$C/bulk.nodedup.tsv" > "$OUT/q.txt"
printf 'PE\t%s\t%s\n'  "$C/bulk_R1.fq" "$C/bulk_R2.fq"  > "$OUT/t_pe.txt"
printf 'PEI\t%s\n'     "$C/bulk_I.fq"                   > "$OUT/t_pei.txt"
printf 'SE\t%s\n'      "$C/bulk_R1.fq"                  > "$OUT/t_se.txt"
printf '%s\n' "$OUT/bulk.pe.kmdb"                       > "$OUT/t_kon.txt"
printf '%s\n' "$OUT/bulk.pe_off.kmdb"                   > "$OUT/t_koff.txt"
sd() {  # sd <name> <targets> [opts] -> $OUT/sd_<name>.tsv (#kmer + count), .tot
	local n=$1 t=$2; shift 2
	"$SD" -k "$OUT/q.txt" -B "$t" -t 2 -o "$OUT/sd_$n.tsv.gz" "$@" 2> "$OUT/sd_$n.log" \
		|| { echo "  (strain_detect $n failed)"; cat "$OUT/sd_$n.log"; return 1; }
	zcat "$OUT/sd_$n.tsv.gz" | awk 'NR>1 && $1!="total_evaluated"' > "$OUT/sd_$n.tsv"
	zcat "$OUT/sd_$n.tsv.gz" | awk '$1=="total_evaluated"{print $2}' > "$OUT/sd_$n.tot"
}
sd pe   "$OUT/t_pe.txt"
sd pei  "$OUT/t_pei.txt"
sd peoff "$OUT/t_pe.txt" --no-pair-dedup
sd kon  "$OUT/t_kon.txt"
sd koff "$OUT/t_koff.txt"
sd se   "$OUT/t_se.txt"
same "FASTQ PE dedup == reference"             "$OUT/sd_pe.tsv"    "$C/bulk.dedup.tsv"
same "FASTQ PE --no-pair-dedup == reference"   "$OUT/sd_peoff.tsv" "$C/bulk.nodedup.tsv"
same "FASTQ PEI == FASTQ PE"                   "$OUT/sd_pei.tsv"   "$OUT/sd_pe.tsv"
same "FASTQ PE dedup == .kmdb dedup"           "$OUT/sd_pe.tsv"    "$OUT/sd_kon.tsv"
same "FASTQ PE no-dedup == .kmdb no-dedup"     "$OUT/sd_peoff.tsv" "$OUT/sd_koff.tsv"
same "total_evaluated FASTQ == .kmdb (dedup)"     "$OUT/sd_pe.tot"    "$OUT/sd_kon.tot"
same "total_evaluated FASTQ == .kmdb (no dedup)"  "$OUT/sd_peoff.tot" "$OUT/sd_koff.tot"
for n in pe pei peoff kon koff se; do
	cs=$(awk '{s+=$2} END{print s}' "$OUT/sd_$n.tsv"); te=$(cat "$OUT/sd_$n.tot")
	[ "$cs" = "$te" ] && ok "strain_detect $n: total_evaluated == sum of all counts ($te)" \
	                  || bad "strain_detect $n: total_evaluated $te vs sum $cs"
done
grep -q "skipped .* query k-mer hits" "$OUT/sd_pe.log" && ok "strain_detect reports skipped hits" || bad "strain_detect stats line"
printf 'PE\t%s\t%s\n' "$C/bad_R1.fq" "$C/bad_R2.fq" > "$OUT/t_bad.txt"
if "$SD" -k "$OUT/q.txt" -B "$OUT/t_bad.txt" -o "$OUT/sd_bad.tsv.gz" 2> "$OUT/sd_bad.log"; then bad "strain_detect: mates out of sync (exit 0)"
else grep -q "mates out of sync" "$OUT/sd_bad.log" && ok "strain_detect: mates out of sync" || bad "strain_detect: wrong error"; fi
if [ -n "${BASE_SD:-}" ]; then
	"$BASE_SD" -k "$OUT/q.txt" -B "$OUT/t_se.txt" -o "$OUT/sd_base_se.tsv.gz" 2>/dev/null
	zcat "$OUT/sd_base_se.tsv.gz" | awk 'NR>1 && $1!="total_evaluated"' > "$OUT/sd_base_se.tsv"
	same "strain_detect SE counts == old binary" "$OUT/sd_se.tsv" "$OUT/sd_base_se.tsv"
fi
if [ -n "${BASE:-}" ]; then   # database from before the total_evaluated key
	printf '%s\n' "$OUT/base_comma.kmdb" > "$OUT/t_kold.txt"
	sd kold "$OUT/t_kold.txt"
	same "old .kmdb (no total_evaluated key): falls back to its sum of counts" "$OUT/sd_kold.tot" "$OUT/sd_koff.tot"
fi

echo "== errors"
err() {  # err <label> <expected stderr regex> <args...>
	local l=$1 re=$2; shift 2
	if "$BIN" -o "$OUT/err.kmdb" "$@" 2> "$OUT/err.log"; then bad "$l (exit 0)";
	elif grep -qE -e "$re" "$OUT/err.log"; then ok "$l"; else bad "$l"; cat "$OUT/err.log"; fi
	[ -e "$OUT/err.kmdb" ] && bad "$l left a database" && rm -f "$OUT/err.kmdb"
}
err "scrub mode rejects --pair-dedup"  "only for -m sample" -m scrub --pair-dedup -i "$C/full_R1.fq,$C/full_R2.fq"
err "mate names out of sync"           "mates out of sync"  -m sample -i "$C/bad_R1.fq,$C/bad_R2.fq"
err "different read counts"            "different numbers of reads" -m sample -i "$C/bulk_R1.fq,$C/short_R2.fq"
err "odd interleaved file"             "odd number of records" -m sample --interleaved -i "$C/odd_I.fq"
err "three files in one -i"            "a comma joins the two mates" -m sample -i "$C/full_R1.fq,$C/full_R2.fq,$C/none_R1.fq"
err "--interleaved with a comma pair"  "--interleaved takes one file" -m sample --interleaved -i "$C/full_R1.fq,$C/full_R2.fq"
"$BIN" -m scrub -c 1 -o "$OUT/scrub.kmdb" -i "$C/full_R1.fq,$C/full_R2.fq" 2>/dev/null \
	&& ok "scrub mode still runs (dedup off)" || bad "scrub mode run"

echo
echo "passed $pass, failed $fail"
[ "$fail" -eq 0 ]
