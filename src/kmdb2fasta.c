/*
	kmdb2fasta — write the k-mers of a .kmdb as unitigs (<db>.kmdb.fasta).

	Overlapping k-mers are glued into unitigs, the maximal non-branching
	paths of the de Bruijn graph of the k-mer set (both strands, canonical
	k-mers as in kmdb: max(forward, reverse complement)). The FASTA spells
	exactly the database's k-mer set:

	  - every k-mer of the database occurs in exactly one unitig, once;
	  - every k-mer window of the FASTA is a database k-mer.

	So counting the FASTA (e.g. kmer_scrub_count -B <db>.kmdb.fasta) gives
	the same k-mers as the database, each once. Counts stored in sample
	databases are NOT written: this is a presence-only export.

	A unitig of L bases holds L - k + 1 k-mers, so long unitigs cost about
	1 byte (plain) or ~0.25 bytes (gzip/zstd) per k-mer instead of ~6.

	Usage: kmdb2fasta [-o out.fasta] [-z] [-l level] [-w width] [-L] [-V] in.kmdb
	  -o  output (default <in>.fasta, or <in>.fasta.gz with -z; "-" = stdout).
	      An -o name ending in .gz turns on gzip by itself.
	  -z  gzip the output (level 6: about 0.4 bytes per k-mer instead of ~1.5;
	      readers using zlib/kseq, e.g. kmer_scrub_count -B, read it directly)
	  -l  gzip level 1..9 (default 6; implies -z)
	  -w  wrap sequence lines at this width (default 0 = one line)
	  -L  long headers: ">u<id> len=<bp> kmers=<n>" (default ">u<id>"; headers
	      are a real share of the file when unitigs are short)
	  -V  verify: re-read the written FASTA (plain or .gz) and check that its
	      k-mers are exactly the database's, each once (not with -o -)

	Memory: 8 bytes per k-mer (sorted keys) + 1 bit (visited) + a 2^24-entry
	prefix index (128 MB) for lookups. Single-threaded.
*/

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <getopt.h>
#include <errno.h>
#include <time.h>
#include <zlib.h>
#include "kmdb.h"

static void *xmalloc(size_t n)
{
	void *p = malloc(n ? n : 1);
	if (!p) { fprintf(stderr, "out of memory (%zu bytes)\n", n); exit(EXIT_FAILURE); }
	return p;
}

static double now_s(void)
{
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* ── k-mer set: sorted canonical keys + prefix index ─────────────────── */

typedef struct {
	uint64_t *keys;
	size_t    n;
	int       k;
	uint64_t  mask;
	int       pshift;     /* key >> pshift = prefix bucket */
	size_t   *pidx;       /* bucket b holds keys[pidx[b] .. pidx[b+1]) */
	uint8_t  *vis;        /* visited bitset */
} kset_t;

static void build_index(kset_t *s)
{
	int pbits = 2 * s->k < 24 ? 2 * s->k : 24;
	s->pshift = 2 * s->k - pbits;
	size_t nb = (size_t)1 << pbits;
	s->pidx = xmalloc((nb + 1) * sizeof(size_t));
	size_t j = 0;
	for (size_t b = 0; b < nb; b++) {
		s->pidx[b] = j;
		while (j < s->n && (s->keys[j] >> s->pshift) == b) j++;
	}
	s->pidx[nb] = s->n;
}

/* index of canonical key x, or -1 */
static inline int64_t find(const kset_t *s, uint64_t x)
{
	size_t b = (size_t)(x >> s->pshift);
	size_t lo = s->pidx[b], hi = s->pidx[b + 1];
	while (lo < hi) {
		size_t mid = lo + (hi - lo) / 2;
		uint64_t v = s->keys[mid];
		if (v < x) lo = mid + 1;
		else if (v > x) hi = mid;
		else return (int64_t)mid;
	}
	return -1;
}

static inline int64_t lookup(const kset_t *s, uint64_t oriented)
{
	return find(s, kmdb_canonical(oriented, s->k));
}

static inline int  is_vis(const kset_t *s, size_t i) { return s->vis[i >> 3] >> (i & 7) & 1; }
static inline void set_vis(kset_t *s, size_t i)      { s->vis[i >> 3] |= (uint8_t)(1u << (i & 7)); }

/* Successors of oriented k-mer u: count, and the (last) one found. */
static inline int successors(const kset_t *s, uint64_t u, uint64_t *v, int64_t *iv)
{
	int n = 0;
	for (uint64_t b = 0; b < 4; b++) {
		uint64_t w = ((u << 2) | b) & s->mask;
		int64_t i = lookup(s, w);
		if (i >= 0) { n++; *v = w; *iv = i; }
	}
	return n;
}

static inline int n_predecessors(const kset_t *s, uint64_t v)
{
	int n = 0;
	for (uint64_t b = 0; b < 4; b++) {
		uint64_t w = (b << (2 * (s->k - 1))) | (v >> 2);
		if (lookup(s, w) >= 0) n++;
	}
	return n;
}

/* ── growable byte buffer ────────────────────────────────────────────── */

typedef struct { uint8_t *b; size_t n, cap; } buf_t;

static inline void buf_push(buf_t *x, uint8_t c)
{
	if (x->n == x->cap) {
		x->cap = x->cap ? x->cap * 2 : 1024;
		x->b = realloc(x->b, x->cap);
		if (!x->b) { fprintf(stderr, "out of memory\n"); exit(EXIT_FAILURE); }
	}
	x->b[x->n++] = c;
}

/* Walk right from oriented k-mer u while the path is unique in both
   directions and unvisited; append the new last base (0..3) of each step. */
static void extend(kset_t *s, uint64_t u, buf_t *out)
{
	for (;;) {
		uint64_t v;
		int64_t iv;
		if (successors(s, u, &v, &iv) != 1) return;
		if (n_predecessors(s, v) != 1) return;   /* v has another way in */
		if (is_vis(s, (size_t)iv)) return;       /* cycle / already used */
		set_vis(s, (size_t)iv);
		buf_push(out, (uint8_t)(v & 3));
		u = v;
	}
}

/* ── output ──────────────────────────────────────────────────────────── */

static const char NT[4] = { 'A', 'C', 'G', 'T' };

/* Output: plain FILE*, or gzip via deflate with a full header (original
   file name + modification time, like `gzip` writes). gzopen() leaves both
   out, which desktop archive managers reject. */
typedef struct {
	FILE     *fp;
	int       gz;
	z_stream  zs;
	gz_header hdr;
	char     *name;
	uint8_t  *zbuf;
} out_t;

enum { ZCHUNK = 1 << 20 };

static void out_die(const char *what)
{
	fprintf(stderr, "error: %s\n", what);
	exit(EXIT_FAILURE);
}

static void out_open(out_t *o, FILE *fp, int gz, int level, const char *path)
{
	memset(o, 0, sizeof *o);
	o->fp = fp;
	o->gz = gz;
	setvbuf(fp, NULL, _IOFBF, 8u << 20);
	if (!gz) return;
	if (deflateInit2(&o->zs, level, Z_DEFLATED, 15 + 16, 8, Z_DEFAULT_STRATEGY) != Z_OK)
		out_die("deflateInit2 failed");
	/* stored name = output file name without directory and ".gz" */
	const char *base = strrchr(path, '/');
	base = base ? base + 1 : path;
	size_t bl = strlen(base);
	if (bl > 3 && strcmp(base + bl - 3, ".gz") == 0) bl -= 3;
	o->name = xmalloc(bl + 1);
	memcpy(o->name, base, bl);
	o->name[bl] = '\0';
	o->hdr.text = 0;
	o->hdr.time = (uLong)time(NULL);
	o->hdr.os   = 3;                               /* Unix */
	o->hdr.name = (Bytef *)o->name;
	if (deflateSetHeader(&o->zs, &o->hdr) != Z_OK) out_die("deflateSetHeader failed");
	o->zbuf = xmalloc(ZCHUNK);
}

static void out_deflate(out_t *o, int flush)
{
	do {
		o->zs.next_out  = o->zbuf;
		o->zs.avail_out = ZCHUNK;
		int rc = deflate(&o->zs, flush);
		if (rc == Z_STREAM_ERROR) out_die("deflate failed");
		size_t have = ZCHUNK - o->zs.avail_out;
		if (have && fwrite(o->zbuf, 1, have, o->fp) != have) out_die("write failed (disk full?)");
	} while (o->zs.avail_out == 0);
}

static void gz_put(out_t *o, const char *p, size_t n)
{
	if (!o->gz) {
		if (n && fwrite(p, 1, n, o->fp) != n) out_die("write failed (disk full?)");
		return;
	}
	while (n > 0) {
		uInt chunk = n > (1u << 30) ? (1u << 30) : (uInt)n;
		o->zs.next_in  = (Bytef *)(uintptr_t)p;
		o->zs.avail_in = chunk;
		out_deflate(o, Z_NO_FLUSH);
		p += chunk; n -= chunk;
	}
}

static int out_close(out_t *o, int is_stdout)
{
	int rc = 0;
	if (o->gz) {
		o->zs.next_in = NULL;
		o->zs.avail_in = 0;
		out_deflate(o, Z_FINISH);
		deflateEnd(&o->zs);
		free(o->zbuf);
		free(o->name);
	}
	if (fflush(o->fp) != 0) rc = -1;
	if (!is_stdout && fclose(o->fp) != 0) rc = -1;
	return rc;
}

static void write_seq(out_t *fp, const uint8_t *seq, size_t len, int width)
{
	static char *line = NULL;
	static size_t cap = 0;
	if (len + 1 > cap) {
		cap = len + 1;
		line = realloc(line, cap);
		if (!line) { fprintf(stderr, "out of memory\n"); exit(EXIT_FAILURE); }
	}
	for (size_t i = 0; i < len; i++) line[i] = NT[seq[i]];
	if (width <= 0) {
		gz_put(fp, line, len);
		gz_put(fp, "\n", 1);
		return;
	}
	for (size_t i = 0; i < len; i += (size_t)width) {
		size_t m = len - i < (size_t)width ? len - i : (size_t)width;
		gz_put(fp, line + i, m);
		gz_put(fp, "\n", 1);
	}
}

static int cmp_u32_desc(const void *a, const void *b)
{
	uint32_t x = *(const uint32_t *)a, y = *(const uint32_t *)b;
	return (x < y) - (x > y);
}

/* ── verify: every FASTA k-mer is in the db, every db k-mer seen once ── */

static int verify(kset_t *s, const char *path)
{
	gzFile fp = gzopen(path, "rb");   /* reads plain files transparently */
	if (!fp) { fprintf(stderr, "verify: cannot open %s: %s\n", path, strerror(errno)); return 1; }
	gzbuffer(fp, 1u << 20);
	memset(s->vis, 0, (s->n + 7) / 8);
	const int k = s->k;
	uint64_t fwd = 0, seen = 0, foreign = 0, dup = 0, windows = 0;
	int len = 0, bol = 1, in_hdr = 0;
	enum { RB = 1 << 20 };
	char *buf = xmalloc(RB);
	int r;
	while ((r = gzread(fp, buf, RB)) > 0) {
		for (int i = 0; i < r; i++) {
			char ch = buf[i];
			if (ch == '\n') { bol = 1; in_hdr = 0; continue; }
			if (bol && ch == '>') { in_hdr = 1; len = 0; }   /* new record */
			bol = 0;
			if (in_hdr || ch == '\r') continue;
			uint8_t c = KMDB_NT2BIT[(unsigned char)ch];
			if (c > 3) { len = 0; continue; }
			fwd = ((fwd << 2) | c) & s->mask;
			if (++len < k) continue;
			windows++;
			int64_t j = lookup(s, fwd);
			if (j < 0) { foreign++; continue; }
			if (is_vis(s, (size_t)j)) dup++;
			else { set_vis(s, (size_t)j); seen++; }
		}
	}
	int gerr = 0;
	if (r < 0) fprintf(stderr, "verify: read error in %s: %s\n", path, gzerror(fp, &gerr));
	free(buf);
	gzclose(fp);
	if (r < 0) return 1;
	int ok = foreign == 0 && dup == 0 && seen == s->n;
	fprintf(stderr, "verify %s: %" PRIu64 " k-mer windows, %" PRIu64 " distinct db k-mers "
	        "seen of %zu, %" PRIu64 " not in db, %" PRIu64 " repeated -> %s\n",
	        path, windows, seen, s->n, foreign, dup, ok ? "OK" : "FAILED");
	return ok ? 0 : 1;
}

/* ── main ────────────────────────────────────────────────────────────── */

static void usage(int st)
{
	fprintf(stderr,
	"Usage: kmdb2fasta [-o out.fasta] [-z] [-l level] [-w width] [-L] [-V] in.kmdb\n"
	"  Writes the database's k-mers as unitigs (overlapping k-mers glued into\n"
	"  maximal non-branching paths). Each k-mer occurs exactly once in the\n"
	"  output, so re-counting it (kmer_scrub_count -B) gives the same k-mer\n"
	"  set. Counts are not written.\n"
	"  -o <file>  output (default <in>.fasta, <in>.fasta.gz with -z; '-' = stdout);\n"
	"             a name ending in .gz turns on gzip\n"
	"  -z         gzip the output (level 6)\n"
	"  -l <1-9>   gzip level (implies -z)\n"
	"  -w <n>     wrap sequence lines at n bases (default 0 = no wrapping)\n"
	"  -L         long headers with unitig length and k-mer count\n"
	"  -V         re-read the output (plain or .gz) and verify it spells exactly\n"
	"             the db k-mers, each once\n");
	exit(st);
}

int main(int argc, char **argv)
{
	const char *out_path = NULL;
	int width = 0, do_verify = 0, long_hdr = 0, gz = 0, level = 6, c;
	while ((c = getopt(argc, argv, "o:w:zl:LVh")) != -1) {
		switch (c) {
		case 'o': out_path = optarg; break;
		case 'w': width = atoi(optarg); break;
		case 'V': do_verify = 1; break;
		case 'L': long_hdr = 1; break;
		case 'z': gz = 1; break;
		case 'l':
			level = atoi(optarg); gz = 1;
			if (level < 1 || level > 9) { fprintf(stderr, "error: -l must be 1..9\n"); return 1; }
			break;
		case 'h': usage(0); break;
		default:  usage(1);
		}
	}
	if (optind != argc - 1) usage(1);
	const char *in = argv[optind];
	char *def_out = NULL;
	if (!out_path) {
		def_out = xmalloc(strlen(in) + 10);
		sprintf(def_out, "%s.fasta%s", in, gz ? ".gz" : "");
		out_path = def_out;
	} else {
		size_t ol = strlen(out_path);
		if (ol > 3 && strcmp(out_path + ol - 3, ".gz") == 0) gz = 1;
	}
	int to_stdout = strcmp(out_path, "-") == 0;
	if (do_verify && to_stdout) { fprintf(stderr, "error: -V needs a file (-o)\n"); return 1; }

	double t0 = now_s();
	kmdb_t *db = kmdb_open(in);
	if (!db) return 1;

	kset_t s = { 0 };
	s.k = (int)db->k;
	s.mask = (s.k == 32) ? UINT64_MAX : ((1ull << (2 * s.k)) - 1);
	s.n = db->n_kmers;
	s.keys = xmalloc(s.n * sizeof(uint64_t));
	if (!(db->flags & KMDB_FLAG_CANON_MAX)) {
		fprintf(stderr, "error: %s does not store canonical (max) k-mers\n", in);
		return 1;
	}

	kmdb_block_buf bb;
	kmdb_block_buf_init(&bb);
	size_t o = 0;
	for (uint64_t b = 0; b < db->n_blocks; b++) {
		long n = kmdb_read_block(db, b, &bb);
		if (n < 0) return 1;
		memcpy(s.keys + o, bb.kmers, (size_t)n * sizeof(uint64_t));
		o += (size_t)n;
	}
	kmdb_block_buf_free(&bb);
	kmdb_close(db);
	build_index(&s);
	s.vis = calloc((s.n + 7) / 8, 1);
	if (!s.vis) { fprintf(stderr, "out of memory\n"); return 1; }
	fprintf(stderr, "loaded %zu k-mers (k=%d) from %s (%.1fs)\n", s.n, s.k, in, now_s() - t0);

	FILE *ofp = to_stdout ? stdout : fopen(out_path, "wb");
	if (!ofp) { fprintf(stderr, "error: cannot create %s: %s\n", out_path, strerror(errno)); return 1; }
	out_t out;
	out_open(&out, ofp, gz, level, to_stdout ? "stdout" : out_path);
	out_t *fp = &out;

	buf_t left = { 0 }, right = { 0 }, seq = { 0 };
	uint32_t *lens = NULL;
	size_t n_unitigs = 0, lens_cap = 0;
	uint64_t total_bases = 0;
	const int k = s.k;

	for (size_t i = 0; i < s.n; i++) {
		if (is_vis(&s, i)) continue;
		set_vis(&s, i);
		uint64_t u = s.keys[i];
		left.n = right.n = 0;
		extend(&s, kmdb_revcomp(u, k), &left);   /* bases left of u, as rc */
		extend(&s, u, &right);

		/* seq = revcomp(left) + u + right */
		seq.n = 0;
		for (size_t j = left.n; j-- > 0; ) buf_push(&seq, (uint8_t)(3 - left.b[j]));
		for (int j = k - 1; j >= 0; j--) buf_push(&seq, (uint8_t)((u >> (2 * j)) & 3));
		for (size_t j = 0; j < right.n; j++) buf_push(&seq, right.b[j]);

		n_unitigs++;
		total_bases += seq.n;
		if (n_unitigs > lens_cap) {
			lens_cap = lens_cap ? lens_cap * 2 : 1 << 16;
			lens = realloc(lens, lens_cap * sizeof(uint32_t));
			if (!lens) { fprintf(stderr, "out of memory\n"); return 1; }
		}
		lens[n_unitigs - 1] = (uint32_t)(seq.n > UINT32_MAX ? UINT32_MAX : seq.n);

		char hdr[96];
		int hl = long_hdr
			? snprintf(hdr, sizeof hdr, ">u%zu len=%zu kmers=%zu\n", n_unitigs, seq.n, seq.n - (size_t)k + 1)
			: snprintf(hdr, sizeof hdr, ">u%zu\n", n_unitigs);
		gz_put(fp, hdr, (size_t)hl);
		write_seq(fp, seq.b, seq.n, width);
	}
	if (out_close(fp, to_stdout) != 0) {
		fprintf(stderr, "error: writing %s failed\n", out_path);
		return 1;
	}

	/* stats */
	uint64_t n50 = 0, acc = 0;
	qsort(lens, n_unitigs, sizeof(uint32_t), cmp_u32_desc);
	for (size_t j = 0; j < n_unitigs; j++) {
		acc += lens[j];
		if (acc * 2 >= total_bases) { n50 = lens[j]; break; }
	}
	fprintf(stderr, "wrote %zu unitigs, %" PRIu64 " bases (%.3f bases/k-mer), "
	        "mean %.1f bp, N50 %" PRIu64 " bp, longest %u bp -> %s%s (%.1fs)\n",
	        n_unitigs, total_bases, s.n ? (double)total_bases / s.n : 0.0,
	        n_unitigs ? (double)total_bases / n_unitigs : 0.0, n50,
	        n_unitigs ? lens[0] : 0, to_stdout ? "stdout" : out_path,
	        gz ? " (gzip)" : "", now_s() - t0);

	int rc = 0;
	if (do_verify) rc = verify(&s, out_path);

	free(left.b); free(right.b); free(seq.b); free(lens);
	free(s.keys); free(s.pidx); free(s.vis); free(def_out);
	return rc;
}