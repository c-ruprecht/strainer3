#!/bin/bash
# genome_compare_presence: new src/genome_compare_presence.c vs the original
# (git $REF): all four outputs must be byte-identical for every option and -t.
set -u; source "$(dirname "$0")/common.sh"
zl="-lz -lm -lpthread"; zd=""; if printf 'int main(void){return 0;}' | gcc -x c - -lzstd -o /dev/null 2>/dev/null; then zd="-DHAVE_ZSTD"; zl="$zl -lzstd"; fi
echo "== build"
gitsrc genome_compare_presence.c > "$W/gcp_ref.c"
gcc -std=gnu17 -O2 -I"$src" $zd -o "$W/gcp_ref" "$W/gcp_ref.c" $zl 2>/dev/null && gcc -std=gnu17 -O2 -I"$src" $zd -o "$W/gcp_new" "$src/genome_compare_presence.c" $zl 2>/dev/null || { echo "build failed"; exit 1; }
ls $TD/scrub_pangenome/*.fna.gz | head -12 > $W/g12.txt
python3 "$here/make_gcp_edge.py" "$W/gedge" >/dev/null; E=$W/gedge/list.txt
same() { local lab=$1 t=$2; shift 2; rm -f $W/r.* $W/n.*
  "$W/gcp_ref" "$@" -o $W/r 2>/dev/null; "$W/gcp_new" "$@" -t $t -o $W/n 2>/dev/null
  local ok=1; for f in $W/r.manifest.tsv $W/r.identity.tsv $W/r.pairs.tsv $W/r.presence.tsv*; do cmp -s $f ${f/\/r./\/n.} || ok=0; done
  [ $ok = 1 ] && echo "  ok      $lab (-t $t)" || { echo "  DIFFERS $lab (-t $t)"; fail=1; }; }
echo "== 12 F. prausnitzii genomes"
for t in 1 2 4; do same "default" $t -A $W/g12.txt; done
same "-U hash order" 3 -A $W/g12.txt -U
same "-z 3"          3 -A $W/g12.txt -z 3
same "-x 0 bitset"   3 -A $W/g12.txt -x 0
same "-F"            3 -A $W/g12.txt -F
same "-m 2 -P 6"     3 -A $W/g12.txt -m 2 -P 6
echo "== 80 synthetic genomes (exotic bases, CRLF, short contigs, dir/list input)"
for a in "" "-U" "-x 5" "-x 0" "-k 4" "-k 32" "-k 33" "-k 51 -U" "-k 64" "-k 1" "-z 5 -m 3"; do same "options: ${a:-default}" 3 -A $E $a; done
same "directory input" 3 -A $W/gedge/dir
echo "== race: -t 1..64 (-x 5: direct + bitset paths)"
"$W/gcp_ref" -A $E -x 5 -o $W/race_ref 2>/dev/null; n=0; d=0
for t in 1 2 3 4 8 16 32 64; do
  "$W/gcp_new" -A $E -x 5 -t $t -o $W/race_new 2>/dev/null; n=$((n+1))
  for f in manifest.tsv identity.tsv pairs.tsv presence.tsv; do cmp -s $W/race_ref.$f $W/race_new.$f || { d=$((d+1)); echo "  MISMATCH -t $t $f"; }; done
done; echo "  $n runs, $d mismatches"; [ $d -eq 0 ] || fail=1
echo "== ThreadSanitizer"
if gcc -std=gnu17 -O1 -g -fsanitize=thread -I"$src" $zd -o "$W/gcp_tsan" "$src/genome_compare_presence.c" $zl 2>/dev/null; then
  for a in "-x 5" "-U"; do TSAN_OPTIONS="halt_on_error=0" "$W/gcp_tsan" -A $E $a -t 8 -o $W/ts 2> $W/tsan.log
    nr=$(grep -c "WARNING: ThreadSanitizer" $W/tsan.log); echo "  $a reports: $nr"; [ "$nr" -eq 0 ] || fail=1; done
else echo "  skipped (no -fsanitize=thread)"; fi
report
