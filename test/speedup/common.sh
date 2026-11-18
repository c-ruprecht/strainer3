# shared setup for the run_tests_*.sh scripts (sourced)
here=$(cd "$(dirname "$0")" && pwd); repo=$(cd "$here/../.." && pwd)
src=$repo/src; TD=$repo/test; W=${WORK:-$here/work}; REF=${REF:-c18dbc4}
mkdir -p "$W"; fail=0
libs="$src/BIO_sequence.c $src/BIO_hash.c $src/genome_compare.c $src/up2bit.c"
gitsrc() { git -C "$repo" show "$REF:src/$1"; }
report() { [ $fail -eq 0 ] && echo "ALL TESTS PASSED" || echo "TESTS FAILED"; exit $fail; }
