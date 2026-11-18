/* genome_compare_presence -- serial inverted k-mer presence index over a set
 * of genomes, plus the pairwise k-mer identity matrix that falls out of it.
 *
 *   genome_compare_presence -A <list|dir> [-k 31] [-t threads] [-O dir] [-n name | -o prefix]
 *
 * THREADS (-t, default 4) never change the output: the default sorted mode
 * builds the index in hash shards (one writer per shard, genomes inserted in
 * strain_id order, so every id list is identical); -U keeps one table filled
 * in the original order (readers parse ahead, one thread inserts), so even the
 * hash order is the same. Sorting and the sweep run in parallel parts that
 * are written in order; see test/speedup/run_tests_presence.sh.
 *
 * WHAT IT DOES
 * ------------
 * One serial pass over the genome list. Genome i gets strain_id i (assigned in
 * list order). Every distinct canonical k-mer of that genome is looked up in a
 * single global hash; a k-mer seen for the first time gets a new entry, one
 * already present just gets i appended to its id list. That is the inverted
 * index:
 *
 *      #kmer                              n_strains   list_strain_ids
 *      AAAA...AAA                         2           0,2
 *      ACGT...TTG                         1           1
 *
 * Because genomes are processed in id order, each id list is already sorted
 * ascending and de-duplication is a single comparison against the last element
 * -- no per-genome scratch set, no sorting of the lists at the end.
 *
 * OUTPUTS (<prefix> = <-O>/<-n>, default ./presence, or -o verbatim)
 *   <prefix>.presence.tsv[.gz|.zst]  inverted index, one row per k-mer
 *   <prefix>.manifest.tsv            strain_id -> file, with per-genome stats
 *   <prefix>.identity.tsv            square N x N Jaccard matrix
 *   <prefix>.pairs.tsv               long form: shared/union/jaccard/containment
 *
 * K-MER CONVENTION
 * ----------------
 * 2-bit packed into an unsigned __int128, so k = 1..64 (31 and 51 both fit;
 * up2bit.c is capped at 32 and uses A=0,C=1,T=2,G=3, which is why the codes
 * here are local: with A=0,C=1,G=2,T=3 the numeric order of a packed k-mer is
 * the lexicographic order of its string, so canonicalisation is one integer
 * compare). Canonical form is the lexicographically LARGER of the forward and
 * reverse-complement strings, matching orient_string()/GEN_hash_sequences() in
 * genome_compare.c so the k-mer strings printed here are the same strings the
 * scrub tools hash. -F disables canonicalisation (forward strand only).
 *
 * Any position containing a non-ACGT character (N, IUPAC codes) is skipped and
 * the rolling encoder restarts after it -- same effect as the contains_N()
 * check in the older code, but without re-scanning each window.
 *
 * IDENTITY IS OVER DISTINCT K-MERS
 * --------------------------------
 *      jaccard(i,j)     = |Ki n Kj| / |Ki u Kj|
 *      containment(i,j) = |Ki n Kj| / |Ki|
 *
 * A k-mer repeated 7x inside an rRNA operon counts once, not seven times.
 * GEN_calculate_coverage() in the old path counts k-mer POSITIONS instead, so
 * the two agree on repeat-free sequence and drift apart on repeat-rich
 * assemblies: do not carry a clustering threshold across from that tool
 * without recalibrating it on a known-ANI ladder.
 *
 * HOW THE MATRIX IS COMPUTED
 * --------------------------
 * Intersections come out of the index itself, in one sweep, split by list
 * length. A k-mer present in L strains contributes to L(L-1)/2 pairs; doing
 * that directly costs L^2/2 increments, while carrying it as one column of a
 * per-strain bitset costs N^2/128 words of popcount-AND. Rare k-mers (the long
 * tail) are cheaper the first way, near-core k-mers the second, and the two
 * are equal at L = N/8 -- which is the default crossover (-x). Both paths are
 * exact and interchangeable; -x 0 forces everything through the bitset path,
 * which is what the test script uses to cross-check them against each other.
 *
 * MEMORY
 * ------
 * ~46 B per distinct k-mer for the hash (32 B entry at 0.7 load), plus 4 B per
 * (k-mer, strain) hit beyond the second -- the first two ids live inline in the
 * entry, so k-mers private to one or two genomes never touch the allocator.
 * The bitset side is n_dense_kmers x N bits and is usually the small term.
 */

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <ctype.h>
#include <errno.h>
#include <unistd.h>
#include <dirent.h>
#include <time.h>
#include <sys/stat.h>
#include <zlib.h>
#include <pthread.h>

#ifdef HAVE_ZSTD
#include <zstd.h>
#endif

#include "kseq.h"
#include "cli_common.h"
KSEQ_INIT(gzFile, gzread)

#define DEFAULT_K              31
#define MAX_K                  64
#define DEFAULT_PREFIX         "presence"
#define INITIAL_TABLE_SLOTS    (1u << 16)   /* grows by doubling; costs ~2x the
                                             * final size in total rehash work */
#define OW_BUFFER              (1u << 20)

typedef unsigned __int128 kmer_t;

/* ---------------------------------------------------------------- helpers */

static double now_s(void)
{
	struct timespec t;
	clock_gettime(CLOCK_MONOTONIC, &t);
	return (double)t.tv_sec + (double)t.tv_nsec / 1e9;
}

static void die(const char *what, const char *detail)
{
	if (detail)
		fprintf(stderr, "genome_compare_presence: %s: %s\n", what, detail);
	else
		fprintf(stderr, "genome_compare_presence: %s\n", what);
	exit(EXIT_FAILURE);
}

static void *xmalloc(size_t n)
{
	void *p = malloc(n ? n : 1);
	if (!p) die("out of memory", NULL);
	return p;
}

static void *xcalloc(size_t n, size_t sz)
{
	void *p = calloc(n ? n : 1, sz ? sz : 1);
	if (!p) die("out of memory", NULL);
	return p;
}

/* ------------------------------------------------------- 2-bit k-mer codes */

/* A=0 C=1 G=2 T=3 so that packed-integer order == string order. */
static uint8_t BASE_CODE[256];
static const char CODE_BASE[4] = { 'A', 'C', 'G', 'T' };

static void init_codes(void)
{
	for (int i = 0; i < 256; i++) BASE_CODE[i] = 255;
	BASE_CODE[(unsigned char)'A'] = BASE_CODE[(unsigned char)'a'] = 0;
	BASE_CODE[(unsigned char)'C'] = BASE_CODE[(unsigned char)'c'] = 1;
	BASE_CODE[(unsigned char)'G'] = BASE_CODE[(unsigned char)'g'] = 2;
	BASE_CODE[(unsigned char)'T'] = BASE_CODE[(unsigned char)'t'] = 3;
}

static void kmer_to_string(kmer_t key, int k, char *out)
{
	for (int i = 0; i < k; i++) {
		unsigned b = (unsigned)((key >> (2 * (k - 1 - i))) & 3);
		out[i] = CODE_BASE[b];
	}
	out[k] = '\0';
}

/* ------------------------------------------------------------- hash table */

#define INLINE_IDS 2

typedef struct {
	kmer_t   key;
	uint32_t n;      /* number of strains holding this k-mer; 0 == empty slot */
	uint32_t cap;    /* == INLINE_IDS while the ids live inline */
	union {
		uint32_t  inl[INLINE_IDS];
		uint32_t *ptr;
	} ids;
} kentry;

typedef struct {
	kentry  *slot;
	uint64_t cap;      /* power of two */
	uint64_t mask;
	uint64_t used;
	uint64_t limit;    /* grow when used exceeds this */
} ktable;

static inline uint32_t *entry_ids(kentry *e)
{
	return e->cap == INLINE_IDS ? e->ids.inl : e->ids.ptr;
}

static inline uint64_t mix64(uint64_t x)
{
	x += 0x9E3779B97F4A7C15ULL;
	x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
	x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
	return x ^ (x >> 31);
}

static inline uint64_t key_hash(kmer_t key)
{
	return mix64((uint64_t)key ^ mix64((uint64_t)(key >> 64)));
}

static void table_init(ktable *t, uint64_t slots)
{
	uint64_t cap = 1;
	while (cap < slots) cap <<= 1;
	t->slot  = (kentry *)xcalloc(cap, sizeof *t->slot);
	t->cap   = cap;
	t->mask  = cap - 1;
	t->used  = 0;
	t->limit = (cap / 10) * 7;          /* 0.7 load factor */
}

static void table_grow(ktable *t)
{
	kentry  *old = t->slot;
	uint64_t oldcap = t->cap;

	t->cap <<= 1;
	t->mask  = t->cap - 1;
	t->limit = (t->cap / 10) * 7;
	t->slot  = (kentry *)xcalloc(t->cap, sizeof *t->slot);

	for (uint64_t i = 0; i < oldcap; i++) {
		if (old[i].n == 0) continue;
		uint64_t j = key_hash(old[i].key) & t->mask;
		while (t->slot[j].n != 0) j = (j + 1) & t->mask;
		t->slot[j] = old[i];         /* moves the ids union verbatim */
	}
	free(old);
}

/* Record that strain `id` contains `key`. Returns 1 if this was a k-mer the
 * strain had not already contributed (i.e. a new distinct k-mer for it). */
static int table_add(ktable *t, kmer_t key, uint32_t id)
{
	uint64_t i = key_hash(key) & t->mask;

	while (t->slot[i].n != 0) {
		if (t->slot[i].key == key) {
			kentry *e = &t->slot[i];
			uint32_t *ids = entry_ids(e);
			if (ids[e->n - 1] == id) return 0;   /* lists are id-ordered */
			if (e->n == e->cap) {
				uint32_t ncap = e->cap * 2;
				uint32_t *heap = (uint32_t *)xmalloc((size_t)ncap * sizeof *heap);
				memcpy(heap, ids, (size_t)e->n * sizeof *heap);
				if (e->cap != INLINE_IDS) free(e->ids.ptr);
				e->ids.ptr = heap;
				e->cap = ncap;
				ids = heap;
			}
			ids[e->n++] = id;
			return 1;
		}
		i = (i + 1) & t->mask;
	}

	t->slot[i].key = key;
	t->slot[i].cap = INLINE_IDS;
	t->slot[i].ids.inl[0] = id;
	t->slot[i].n = 1;
	if (++t->used > t->limit) table_grow(t);
	return 1;
}

static void table_free(ktable *t)
{
	for (uint64_t i = 0; i < t->cap; i++)
		if (t->slot[i].n != 0 && t->slot[i].cap != INLINE_IDS)
			free(t->slot[i].ids.ptr);
	free(t->slot);
	t->slot = NULL;
}

/* --------------------------------------------------------- output writers */

typedef enum { OW_PLAIN, OW_GZ, OW_ZSTD } ow_kind;

typedef struct {
	ow_kind  kind;
	FILE    *fp;
	gzFile   gz;
	char    *buf;
	size_t   len, cap;
	char    *path;
#ifdef HAVE_ZSTD
	ZSTD_CCtx *cctx;
	char      *zbuf;
	size_t     zcap;
#endif
} ow;

static int has_suffix(const char *s, const char *suf)
{
	size_t ls = strlen(s), lf = strlen(suf);
	return ls >= lf && strcmp(s + ls - lf, suf) == 0;
}

/* compression level for .zst / .gz outputs, set from -z */
static int g_zlevel = 3;

static ow *ow_open(const char *path)
{
	ow *w = (ow *)xcalloc(1, sizeof *w);
	w->path = strdup(path);
	w->cap  = OW_BUFFER;
	w->buf  = (char *)xmalloc(w->cap);

	if (has_suffix(path, ".zst")) {
#ifdef HAVE_ZSTD
		w->kind = OW_ZSTD;
		w->fp = fopen(path, "wb");
		if (!w->fp) die("cannot open output", path);
		w->cctx = ZSTD_createCCtx();
		if (!w->cctx) die("ZSTD_createCCtx failed", path);
		ZSTD_CCtx_setParameter(w->cctx, ZSTD_c_compressionLevel, g_zlevel);
		w->zcap = ZSTD_CStreamOutSize();
		w->zbuf = (char *)xmalloc(w->zcap);
#else
		die("built without zstd support (rebuild with libzstd, "
		    "or use a .gz / plain output name)", path);
#endif
	} else if (has_suffix(path, ".gz")) {
		w->kind = OW_GZ;
		char gzmode[8];
		snprintf(gzmode, sizeof gzmode, "wb%d",
		         g_zlevel < 1 ? 1 : g_zlevel > 9 ? 9 : g_zlevel);
		w->gz = gzopen(path, gzmode);
		if (!w->gz) die("cannot open output", path);
	} else {
		w->kind = OW_PLAIN;
		w->fp = fopen(path, "w");
		if (!w->fp) die("cannot open output", path);
	}
	return w;
}

static void ow_drain(ow *w, int final)
{
	switch (w->kind) {
	case OW_PLAIN:
		if (w->len && fwrite(w->buf, 1, w->len, w->fp) != w->len)
			die("short write", w->path);
		w->len = 0;
		break;
	case OW_GZ:
		if (w->len && gzwrite(w->gz, w->buf, (unsigned)w->len) <= 0)
			die("short write", w->path);
		w->len = 0;
		break;
	case OW_ZSTD:
#ifdef HAVE_ZSTD
	{
		ZSTD_inBuffer in = { w->buf, w->len, 0 };
		ZSTD_EndDirective mode = final ? ZSTD_e_end : ZSTD_e_continue;
		for (;;) {
			ZSTD_outBuffer out = { w->zbuf, w->zcap, 0 };
			size_t rem = ZSTD_compressStream2(w->cctx, &out, &in, mode);
			if (ZSTD_isError(rem)) die("zstd error", ZSTD_getErrorName(rem));
			if (out.pos && fwrite(w->zbuf, 1, out.pos, w->fp) != out.pos)
				die("short write", w->path);
			if (final) { if (rem == 0) break; }
			else if (in.pos == in.size) break;
		}
		w->len = 0;
	}
#endif
		break;
	}
	(void)final;
}

static inline void ow_reserve(ow *w, size_t need)
{
	if (w->len + need > w->cap) {
		ow_drain(w, 0);
		if (need > w->cap) {
			w->cap = need * 2;
			w->buf = (char *)realloc(w->buf, w->cap);
			if (!w->buf) die("out of memory", NULL);
		}
	}
}

static inline void ow_puts(ow *w, const char *s, size_t n)
{
	ow_reserve(w, n);
	memcpy(w->buf + w->len, s, n);
	w->len += n;
}

static inline void ow_u32(ow *w, uint32_t v)
{
	char tmp[12];
	int i = 12;
	do { tmp[--i] = (char)('0' + v % 10); v /= 10; } while (v);
	ow_puts(w, tmp + i, (size_t)(12 - i));
}

static inline void ow_ch(ow *w, char c)
{
	ow_reserve(w, 1);
	w->buf[w->len++] = c;
}

static void ow_close(ow *w)
{
	ow_drain(w, 1);
	switch (w->kind) {
	case OW_GZ:   gzclose(w->gz); break;
	case OW_ZSTD:
#ifdef HAVE_ZSTD
		ZSTD_freeCCtx(w->cctx);
		free(w->zbuf);
#endif
		fclose(w->fp);
		break;
	case OW_PLAIN: fclose(w->fp); break;
	}
	free(w->buf);
	free(w->path);
	free(w);
}

/* ------------------------------------------------------------ input paths */

typedef struct {
	char **path;
	size_t n, cap;
} pathvec;

static void pv_push(pathvec *v, const char *p)
{
	if (v->n == v->cap) {
		v->cap = v->cap ? v->cap * 2 : 64;
		v->path = (char **)realloc(v->path, v->cap * sizeof *v->path);
		if (!v->path) die("out of memory", NULL);
	}
	v->path[v->n++] = strdup(p);
}

static int is_fasta_name(const char *name)
{
	static const char *suf[] = {
		".fna.gz", ".fa.gz", ".fasta.gz", ".fsa.gz", ".fna", ".fa",
		".fasta", ".fsa", ".ffn", ".ffn.gz", NULL
	};
	for (int i = 0; suf[i]; i++)
		if (has_suffix(name, suf[i])) return 1;
	return 0;
}

static int cmp_str(const void *a, const void *b)
{
	char *const *x = (char *const *)a;
	char *const *y = (char *const *)b;
	return strcmp(*x, *y);
}

/* Append every FASTA file directly inside `dir`, sorted, so a run is
 * reproducible regardless of readdir order. */
static void add_directory(pathvec *v, const char *dir)
{
	DIR *d = opendir(dir);
	if (!d) { fprintf(stderr, "warning: cannot open directory %s: %s\n",
	                  dir, strerror(errno)); return; }

	pathvec local = { NULL, 0, 0 };
	struct dirent *de;
	while ((de = readdir(d)) != NULL) {
		if (de->d_name[0] == '.') continue;
		if (!is_fasta_name(de->d_name)) continue;
		char buf[4096];
		snprintf(buf, sizeof buf, "%s/%s", dir, de->d_name);
		struct stat st;
		if (stat(buf, &st) == 0 && S_ISDIR(st.st_mode)) continue;
		pv_push(&local, buf);
	}
	closedir(d);

	qsort(local.path, local.n, sizeof *local.path, cmp_str);
	for (size_t i = 0; i < local.n; i++) {
		pv_push(v, local.path[i]);
		free(local.path[i]);
	}
	free(local.path);

	if (local.n == 0)
		fprintf(stderr, "warning: no FASTA files found in %s\n", dir);
}

/* -A accepts either a directory of genomes or a text file listing genome
 * paths, where each listed path may itself be a file or a directory. */
static pathvec collect_inputs(const char *arg)
{
	pathvec v = { NULL, 0, 0 };
	struct stat st;

	if (stat(arg, &st) != 0) die("cannot stat input", arg);

	if (S_ISDIR(st.st_mode)) {
		add_directory(&v, arg);
		return v;
	}

	FILE *fp = fopen(arg, "r");
	if (!fp) die("cannot open genome list", arg);

	char *line = NULL;
	size_t cap = 0;
	ssize_t len;
	while ((len = getline(&line, &cap, fp)) != -1) {
		while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r' ||
		                   line[len - 1] == ' '  || line[len - 1] == '\t'))
			line[--len] = '\0';
		char *p = line;
		while (*p == ' ' || *p == '\t') p++;
		if (*p == '\0' || *p == '#') continue;

		struct stat es;
		if (stat(p, &es) != 0) {
			fprintf(stderr, "warning: skipping unreadable path %s\n", p);
			continue;
		}
		if (S_ISDIR(es.st_mode)) add_directory(&v, p);
		else                     pv_push(&v, p);
	}
	free(line);
	fclose(fp);
	return v;
}

static const char *base_name(const char *path)
{
	const char *s = strrchr(path, '/');
	return s ? s + 1 : path;
}

/* basename with one or two FASTA extensions removed: foo.fna.gz -> foo */
static char *sample_id_of(const char *path)
{
	char *s = strdup(base_name(path));
	for (int round = 0; round < 2; round++) {
		char *dot = strrchr(s, '.');
		if (!dot) break;
		if (strcmp(dot, ".gz") == 0 || strcmp(dot, ".zst") == 0 ||
		    strcmp(dot, ".bz2") == 0 || strcmp(dot, ".fna") == 0 ||
		    strcmp(dot, ".fa") == 0 || strcmp(dot, ".fasta") == 0 ||
		    strcmp(dot, ".fsa") == 0 || strcmp(dot, ".ffn") == 0)
			*dot = '\0';
		else break;
	}
	return s;
}

/* -------------------------------------------------------- genome scanning */

typedef struct {
	uint64_t bases;
	uint64_t positions;      /* k-mer positions actually encoded (non-N) */
	uint64_t distinct;       /* distinct canonical k-mers contributed */
	uint32_t contigs;
	uint32_t short_contigs;  /* contigs shorter than k */
} gstat;

static int scan_genome(const char *path, int k, int canonical, uint32_t id,
                       ktable *t, gstat *gs)
{
	gzFile fp = gzopen(path, "r");
	if (!fp) { fprintf(stderr, "warning: cannot open %s\n", path); return -1; }

	kseq_t *seq = kseq_init(fp);
	const kmer_t mask = (k == 64) ? (~(kmer_t)0)
	                              : ((((kmer_t)1) << (2 * k)) - 1);
	const unsigned shift = (unsigned)(2 * (k - 1));
	int l;

	while ((l = kseq_read(seq)) >= 0) {
		gs->contigs++;
		gs->bases += (uint64_t)seq->seq.l;
		if ((size_t)l < (size_t)k) { gs->short_contigs++; continue; }

		kmer_t fwd = 0, rc = 0;
		int filled = 0;
		const char *s = seq->seq.s;

		for (size_t i = 0; i < seq->seq.l; i++) {
			uint8_t code = BASE_CODE[(unsigned char)s[i]];
			if (code > 3) { filled = 0; fwd = rc = 0; continue; }

			fwd = ((fwd << 2) | (kmer_t)code) & mask;
			rc  = (rc >> 2) | (((kmer_t)(3 - code)) << shift);

			if (filled < k) filled++;
			if (filled < k) continue;

			gs->positions++;
			kmer_t key = fwd;
			if (canonical && rc > fwd) key = rc;
			gs->distinct += (uint64_t)table_add(t, key, id);
		}
	}

	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}

/* ------------------------------------------------------------------ popcnt */

static inline uint64_t and_popcount(const uint64_t *a, const uint64_t *b,
                                    uint64_t words)
{
	uint64_t n = 0;
	for (uint64_t w = 0; w < words; w++)
		n += (uint64_t)__builtin_popcountll(a[w] & b[w]);
	return n;
}


/* ------------------------------------------------------- parallel helpers */

typedef void (*par_fn)(void *ctx, int tid, int nth);
typedef struct { par_fn fn; void *ctx; int tid, nth; } par_arg;

static void *par_tramp(void *a)
{
	par_arg *p = (par_arg *)a;
	p->fn(p->ctx, p->tid, p->nth);
	return NULL;
}

/* run fn(ctx, tid, nth) on nth threads (tid 0 on the caller) and join all */
static void par_run(int nth, par_fn fn, void *ctx)
{
	if (nth <= 1) { fn(ctx, 0, 1); return; }
	pthread_t *th = (pthread_t *)xmalloc(sizeof(pthread_t) * (size_t)nth);
	par_arg   *pa = (par_arg *)xmalloc(sizeof(par_arg) * (size_t)nth);
	for (int t = 0; t < nth; t++) { pa[t].fn = fn; pa[t].ctx = ctx; pa[t].tid = t; pa[t].nth = nth; }
	for (int t = 1; t < nth; t++)
		if (pthread_create(&th[t], NULL, par_tramp, &pa[t]) != 0) die("pthread_create failed", NULL);
	fn(ctx, 0, nth);
	for (int t = 1; t < nth; t++) pthread_join(th[t], NULL);
	free(th); free(pa);
}

/* ------------------------------------------- pipelined genome k-mer reader
 * Worker threads read, 2-bit encode and hash genomes ahead of the main
 * thread (at most W genomes in flight). The main thread inserts them into the
 * single global table strictly in strain_id order and in sequence position
 * order -- the exact sequence of table_add() calls of the serial version, so
 * the table (and with it every output, including -U hash order) is identical.
 */

typedef struct {
	kmer_t   *key;
	uint64_t *hv;            /* key_hash(key), computed by the reader */
	size_t    n, cap;
	gstat     gs;
	uint32_t  id;
	int       ready;
} gbuf;

typedef struct {
	const pathvec  *files;
	int             k, canonical;
	uint32_t        N, W;
	gbuf           *ring;
	uint32_t        next, consumed;
	pthread_mutex_t m;
	pthread_cond_t  cv;
} gprod;

static void gbuf_push(gbuf *b, kmer_t key)
{
	if (b->n == b->cap) {
		b->cap = b->cap ? b->cap * 2 : (1u << 20);
		b->key = (kmer_t *)realloc(b->key, b->cap * sizeof *b->key);
		b->hv  = (uint64_t *)realloc(b->hv, b->cap * sizeof *b->hv);
		if (!b->key || !b->hv) die("out of memory", NULL);
	}
	b->key[b->n] = key;
	b->hv[b->n]  = key_hash(key);
	b->n++;
}

/* same parsing and statistics as scan_genome(); k-mers go to the buffer */
static void read_genome(const char *path, int k, int canonical, gbuf *b)
{
	memset(&b->gs, 0, sizeof b->gs);
	b->n = 0;
	gzFile fp = gzopen(path, "r");
	if (!fp) { fprintf(stderr, "warning: cannot open %s\n", path); return; }

	kseq_t *seq = kseq_init(fp);
	const kmer_t mask = (k == 64) ? (~(kmer_t)0)
	                              : ((((kmer_t)1) << (2 * k)) - 1);
	const unsigned shift = (unsigned)(2 * (k - 1));
	int l;
	gstat *gs = &b->gs;

	while ((l = kseq_read(seq)) >= 0) {
		gs->contigs++;
		gs->bases += (uint64_t)seq->seq.l;
		if ((size_t)l < (size_t)k) { gs->short_contigs++; continue; }

		kmer_t fwd = 0, rc = 0;
		int filled = 0;
		const char *s = seq->seq.s;

		for (size_t i = 0; i < seq->seq.l; i++) {
			uint8_t code = BASE_CODE[(unsigned char)s[i]];
			if (code > 3) { filled = 0; fwd = rc = 0; continue; }

			fwd = ((fwd << 2) | (kmer_t)code) & mask;
			rc  = (rc >> 2) | (((kmer_t)(3 - code)) << shift);

			if (filled < k) filled++;
			if (filled < k) continue;

			gs->positions++;
			kmer_t key = fwd;
			if (canonical && rc > fwd) key = rc;
			gbuf_push(b, key);
		}
	}
	kseq_destroy(seq);
	gzclose(fp);
}

static void reader_thread(void *ctx, int tid, int nth)
{
	gprod *g = (gprod *)ctx;
	(void)nth;
	if (tid == 0) return;                 /* tid 0 is the inserting main thread */
	for (;;) {
		pthread_mutex_lock(&g->m);
		while (g->next < g->N && g->next >= g->consumed + g->W)
			pthread_cond_wait(&g->cv, &g->m);
		if (g->next >= g->N) { pthread_mutex_unlock(&g->m); break; }
		uint32_t id = g->next++;
		pthread_mutex_unlock(&g->m);

		gbuf *b = &g->ring[id % g->W];    /* free: genome id - W was consumed */
		read_genome(g->files->path[id], g->k, g->canonical, b);

		pthread_mutex_lock(&g->m);
		b->id = id; b->ready = 1;
		pthread_cond_broadcast(&g->cv);
		pthread_mutex_unlock(&g->m);
	}
}

/* table_add() with the hash already computed */
static inline int table_add_h(ktable *t, kmer_t key, uint64_t hv, uint32_t id)
{
	uint64_t i = hv & t->mask;

	while (t->slot[i].n != 0) {
		if (t->slot[i].key == key) {
			kentry *e = &t->slot[i];
			uint32_t *ids = entry_ids(e);
			if (ids[e->n - 1] == id) return 0;   /* lists are id-ordered */
			if (e->n == e->cap) {
				uint32_t ncap = e->cap * 2;
				uint32_t *heap = (uint32_t *)xmalloc((size_t)ncap * sizeof *heap);
				memcpy(heap, ids, (size_t)e->n * sizeof *heap);
				if (e->cap != INLINE_IDS) free(e->ids.ptr);
				e->ids.ptr = heap;
				e->cap = ncap;
				ids = heap;
			}
			ids[e->n++] = id;
			return 1;
		}
		i = (i + 1) & t->mask;
	}

	t->slot[i].key = key;
	t->slot[i].cap = INLINE_IDS;
	t->slot[i].ids.inl[0] = id;
	t->slot[i].n = 1;
	if (++t->used > t->limit) table_grow(t);
	return 1;
}

typedef struct {
	gprod   *g;
	ktable  *t;
	gstat   *gs;
} build_ctx;

static void inserter(build_ctx *bc)
{
	gprod *g = bc->g;
	ktable *t = bc->t;
	const size_t D = 16;                  /* prefetch distance (k-mers) */
	for (uint32_t i = 0; i < g->N; i++) {
		gbuf *b = &g->ring[i % g->W];
		pthread_mutex_lock(&g->m);
		while (!(b->ready && b->id == i)) pthread_cond_wait(&g->cv, &g->m);
		pthread_mutex_unlock(&g->m);

		uint64_t distinct = 0;
		for (size_t j = 0; j < b->n; j++) {
			if (j + D < b->n) __builtin_prefetch(&t->slot[b->hv[j + D] & t->mask]);
			distinct += (uint64_t)table_add_h(t, b->key[j], b->hv[j], i);
		}
		bc->gs[i] = b->gs;
		bc->gs[i].distinct = distinct;

		pthread_mutex_lock(&g->m);
		b->ready = 0;
		g->consumed = i + 1;
		pthread_cond_broadcast(&g->cv);
		pthread_mutex_unlock(&g->m);

		if (g->N > 20 && (i + 1) % 50 == 0)
			fprintf(stderr, "[presence]   %u/%u genomes, %" PRIu64
			        " distinct k-mers so far\n", i + 1, g->N, t->used);
	}
}

static void build_thread(void *ctx, int tid, int nth)
{
	build_ctx *bc = (build_ctx *)ctx;
	if (tid == 0) inserter(bc);
	else reader_thread(bc->g, tid, nth);
}

/* --------------------------------------------- sharded build (sorted mode)
 * The table is split into NSHARD independent tables by the top hash bits.
 * Genomes are read in batches (one per thread); then every thread inserts the
 * whole batch, genome by genome in id order, into the shards it owns. Each
 * k-mer therefore sees exactly the same sequence of table_add() calls as in
 * the serial build (its strain list and the per-genome distinct counts are
 * identical); only the table layout differs, which is invisible once the
 * entries are sorted by k-mer. Threads never touch another thread's shard. */

#define NSHARD 64
#define SHARD_OF(hv) ((unsigned)((hv) >> 58))

typedef struct {
	const pathvec *files;
	int            k, canonical;
	ktable        *shard;
	gbuf          *buf;        /* one per thread */
	uint64_t     **dist;       /* [thread][genome] distinct k-mers */
	gstat         *gs;
	uint32_t       b0, nb;     /* current batch */
} shard_ctx;

static void shard_read(void *c, int tid, int nth)
{
	shard_ctx *s = (shard_ctx *)c;
	(void)nth;
	if ((uint32_t)tid >= s->nb) return;
	read_genome(s->files->path[s->b0 + (uint32_t)tid], s->k, s->canonical, &s->buf[tid]);
	s->gs[s->b0 + (uint32_t)tid] = s->buf[tid].gs;
}

static void shard_insert(void *c, int tid, int nth)
{
	shard_ctx *s = (shard_ctx *)c;
	const size_t D = 16;
	for (uint32_t g = 0; g < s->nb; g++) {
		const gbuf *b = &s->buf[g];
		const uint32_t id = s->b0 + g;
		uint64_t d = 0;
		for (size_t j = 0; j < b->n; j++) {
			if (j + D < b->n) {
				uint64_t h2 = b->hv[j + D];
				if ((int)(SHARD_OF(h2) % (unsigned)nth) == tid) {
					const ktable *t2 = &s->shard[SHARD_OF(h2)];
					__builtin_prefetch(&t2->slot[h2 & t2->mask]);
				}
			}
			uint64_t hv = b->hv[j];
			unsigned sh = SHARD_OF(hv);
			if ((int)(sh % (unsigned)nth) != tid) continue;
			d += (uint64_t)table_add_h(&s->shard[sh], b->key[j], hv, id);
		}
		s->dist[tid][id] += d;
	}
}

static void build_index_sharded(const pathvec *files, int k, int canonical, uint32_t N,
                                ktable *shard, gstat *gs, int nth)
{
	shard_ctx s;
	memset(&s, 0, sizeof s);
	s.files = files; s.k = k; s.canonical = canonical; s.shard = shard; s.gs = gs;
	s.buf  = (gbuf *)xcalloc((size_t)nth, sizeof(gbuf));
	s.dist = (uint64_t **)xmalloc(sizeof(uint64_t *) * (size_t)nth);
	for (int t = 0; t < nth; t++) s.dist[t] = (uint64_t *)xcalloc(N, sizeof(uint64_t));
	for (uint32_t b0 = 0; b0 < N; b0 += (uint32_t)nth) {
		s.b0 = b0;
		s.nb = (N - b0 < (uint32_t)nth) ? N - b0 : (uint32_t)nth;
		par_run(nth, shard_read, &s);
		par_run(nth, shard_insert, &s);
		if (N > 20) {
			uint64_t used = 0;
			for (int x = 0; x < NSHARD; x++) used += shard[x].used;
			for (uint32_t i = b0; i < b0 + s.nb; i++)
				if ((i + 1) % 50 == 0)
					fprintf(stderr, "[presence]   %u/%u genomes, %" PRIu64
					        " distinct k-mers so far\n", i + 1, N, used);
		}
	}
	for (uint32_t i = 0; i < N; i++) {
		uint64_t d = 0;
		for (int t = 0; t < nth; t++) d += s.dist[t][i];
		gs[i].distinct = d;
	}
	for (int t = 0; t < nth; t++) { free(s.dist[t]); free(s.buf[t].key); free(s.buf[t].hv); }
	free(s.dist); free(s.buf);
}

/* build the index with `nth` threads (1 inserter + nth-1 readers) */
static void build_index(const pathvec *files, int k, int canonical, uint32_t N,
                        ktable *t, gstat *gs, int nth)
{
	if (nth < 2) {                         /* serial: original loop */
		for (uint32_t i = 0; i < N; i++) {
			scan_genome(files->path[i], k, canonical, i, t, &gs[i]);
			if (N > 20 && (i + 1) % 50 == 0)
				fprintf(stderr, "[presence]   %u/%u genomes, %" PRIu64
				        " distinct k-mers so far\n", i + 1, N, t->used);
		}
		return;
	}
	gprod g;
	memset(&g, 0, sizeof g);
	g.files = files; g.k = k; g.canonical = canonical; g.N = N;
	g.W = (uint32_t)(nth + 1);
	g.ring = (gbuf *)xcalloc(g.W, sizeof *g.ring);
	pthread_mutex_init(&g.m, NULL);
	pthread_cond_init(&g.cv, NULL);
	build_ctx bc = { &g, t, gs };
	par_run(nth, build_thread, &bc);
	for (uint32_t w = 0; w < g.W; w++) { free(g.ring[w].key); free(g.ring[w].hv); }
	free(g.ring);
	pthread_mutex_destroy(&g.m);
	pthread_cond_destroy(&g.cv);
}

/* ------------------------------------------------------ sorted entry list
 * The sweep walks `sk`: slot index + list length per entry, in output order
 * (k-mer order, or table slot order with -U). With k <= 32 the k-mer is kept
 * inline and sorted by a parallel bucket sort; k-mers are unique, so any
 * correct sort gives the same order as the original qsort. */

typedef struct { uint64_t key; kentry *e; } skent;

static int cmp_skent(const void *a, const void *b)
{
	uint64_t x = ((const skent *)a)->key, y = ((const skent *)b)->key;
	return (x > y) - (x < y);
}

/* k > 32: full 128-bit key, same order as cmp_entry_ptr */
static int cmp_skent128(const void *a, const void *b)
{
	kmer_t x = ((const skent *)a)->e->key, y = ((const skent *)b)->e->key;
	return (x > y) - (x < y);
}

typedef struct {
	skent   *src, *dst;
	size_t   n;
	int      shift, nb;
	size_t  *hist;          /* nth x nb */
	size_t  *bstart;        /* nb + 1 */
	int      next_bucket;   /* atomic */
} sort_ctx;

static void sort_hist(void *c, int tid, int nth)
{
	sort_ctx *s = (sort_ctx *)c;
	size_t lo = s->n * (size_t)tid / (size_t)nth, hi = s->n * (size_t)(tid + 1) / (size_t)nth;
	size_t *h = s->hist + (size_t)tid * (size_t)s->nb;
	for (size_t i = lo; i < hi; i++) h[s->src[i].key >> s->shift]++;
}

static void sort_scatter(void *c, int tid, int nth)
{
	sort_ctx *s = (sort_ctx *)c;
	size_t lo = s->n * (size_t)tid / (size_t)nth, hi = s->n * (size_t)(tid + 1) / (size_t)nth;
	size_t *h = s->hist + (size_t)tid * (size_t)s->nb;     /* holds write offsets now */
	for (size_t i = lo; i < hi; i++) s->dst[h[s->src[i].key >> s->shift]++] = s->src[i];
}

static void sort_buckets(void *c, int tid, int nth)
{
	sort_ctx *s = (sort_ctx *)c;
	(void)tid; (void)nth;
	for (;;) {
		int b = __atomic_fetch_add(&s->next_bucket, 1, __ATOMIC_RELAXED);
		if (b >= s->nb) break;
		size_t lo = s->bstart[b], hi = s->bstart[b + 1];
		if (hi - lo > 1) qsort(s->dst + lo, hi - lo, sizeof(skent), cmp_skent);
	}
}

/* sorts sk[0..n) by key; returns the array holding the result */
static skent *parallel_sort(skent *sk, size_t n, int bits, int nth)
{
	skent *tmp = (skent *)xmalloc(n * sizeof *tmp);
	sort_ctx s;
	memset(&s, 0, sizeof s);
	s.src = sk; s.dst = tmp; s.n = n;
	int bb = bits < 12 ? bits : 12;
	s.nb = 1 << bb;
	s.shift = bits - bb;
	s.hist = (size_t *)xcalloc((size_t)nth * (size_t)s.nb, sizeof(size_t));
	s.bstart = (size_t *)xcalloc((size_t)s.nb + 1, sizeof(size_t));
	par_run(nth, sort_hist, &s);
	size_t off = 0;
	for (int b = 0; b < s.nb; b++) {
		s.bstart[b] = off;
		for (int t = 0; t < nth; t++) {
			size_t c = s.hist[(size_t)t * (size_t)s.nb + (size_t)b];
			s.hist[(size_t)t * (size_t)s.nb + (size_t)b] = off;
			off += c;
		}
	}
	s.bstart[s.nb] = off;
	par_run(nth, sort_scatter, &s);
	par_run(nth, sort_buckets, &s);
	free(s.hist); free(s.bstart);
	free(sk);
	return tmp;
}

/* ------------------------------------------------------------ parallel sweep
 * Rounds of entries are split into one contiguous part per thread. Each part
 * formats its presence rows into its own buffer (written to the output in
 * part order afterwards), counts its statistics and adds its direct-path
 * intersections to a private N x N matrix (summed at the end). Dense entries
 * get their global bitset column from a prefix count, and bits are set with
 * an atomic OR, since two parts can touch the same 64-bit word. */

typedef struct {
	char  *buf;
	size_t len, cap;
} pbuf;

static inline void pb_reserve(pbuf *b, size_t need)
{
	if (b->len + need > b->cap) {
		b->cap = (b->len + need) * 2;
		b->buf = (char *)realloc(b->buf, b->cap);
		if (!b->buf) die("out of memory", NULL);
	}
}

static inline void pb_u32(pbuf *b, uint32_t v)
{
	char tmp[12];
	int i = 12;
	do { tmp[--i] = (char)('0' + v % 10); v /= 10; } while (v);
	memcpy(b->buf + b->len, tmp + i, (size_t)(12 - i));
	b->len += (size_t)(12 - i);
}

typedef struct {
	skent     *sk;
	size_t     lo, hi;               /* this round */
	int        k;
	uint32_t   N;
	uint64_t   thresh, words;
	uint64_t  *bits;
	long       min_strains, max_strains;
	uint64_t **inter_t;              /* per thread, N x N */
	pbuf      *out;                  /* per thread */
	uint64_t  *dcount, *dbase;       /* per thread */
	uint64_t  *singletons, *core, *written;   /* per thread */
} sweep_ctx;

static inline void part_range(const sweep_ctx *s, int tid, int nth, size_t *lo, size_t *hi)
{
	size_t n = s->hi - s->lo;
	*lo = s->lo + n * (size_t)tid / (size_t)nth;
	*hi = s->lo + n * (size_t)(tid + 1) / (size_t)nth;
}

static void sweep_count_dense(void *c, int tid, int nth)
{
	sweep_ctx *s = (sweep_ctx *)c;
	size_t lo, hi; part_range(s, tid, nth, &lo, &hi);
	uint64_t d = 0;
	for (size_t i = lo; i < hi; i++) if (s->sk[i].e->n > s->thresh) d++;
	s->dcount[tid] = d;
}

static void sweep_part(void *c, int tid, int nth)
{
	sweep_ctx *s = (sweep_ctx *)c;
	size_t lo, hi; part_range(s, tid, nth, &lo, &hi);
	const uint32_t N = s->N;
	uint64_t *inter = s->inter_t[tid];
	uint64_t dense_col = s->dbase[tid];
	pbuf *pb = &s->out[tid];
	pb->len = 0;
	char kbuf[MAX_K + 1];
	/* thread-local counters: the per-thread result slots share cache lines,
	   so they are written once at the end of the part */
	uint64_t n_single = 0, n_core = 0, n_written = 0;

	for (size_t i = lo; i < hi; i++) {
		if (i + 8 < hi) __builtin_prefetch(s->sk[i + 8].e);
		kentry *e = s->sk[i].e;
		uint32_t *ids = entry_ids(e);

		if (e->n == 1) n_single++;
		if (e->n == N) n_core++;

		if (e->n > s->thresh) {
			uint64_t w = dense_col >> 6;
			uint64_t b = 1ULL << (dense_col & 63);
			for (uint32_t a = 0; a < e->n; a++)
				__atomic_fetch_or(&s->bits[(size_t)ids[a] * s->words + w], b, __ATOMIC_RELAXED);
			dense_col++;
		} else {
			for (uint32_t a = 0; a < e->n; a++)
				for (uint32_t b = a + 1; b < e->n; b++)
					inter[(size_t)ids[a] * N + ids[b]]++;
		}

		if ((long)e->n < s->min_strains) continue;
		if (s->max_strains > 0 && (long)e->n > s->max_strains) continue;

		pb_reserve(pb, (size_t)s->k + 2 + 11 + (size_t)e->n * 11 + 1);
		kmer_to_string(e->key, s->k, kbuf);
		memcpy(pb->buf + pb->len, kbuf, (size_t)s->k); pb->len += (size_t)s->k;
		pb->buf[pb->len++] = '\t';
		pb_u32(pb, e->n);
		pb->buf[pb->len++] = '\t';
		for (uint32_t a = 0; a < e->n; a++) {
			if (a) pb->buf[pb->len++] = ',';
			pb_u32(pb, ids[a]);
		}
		pb->buf[pb->len++] = '\n';
		n_written++;
	}
	s->singletons[tid] += n_single;
	s->core[tid]       += n_core;
	s->written[tid]    += n_written;
}

typedef struct {
	uint64_t *inter, *bits, words;
	uint32_t  N;
	int       next_row;       /* atomic */
} pop_ctx;

static void popcount_rows(void *c, int tid, int nth)
{
	pop_ctx *p = (pop_ctx *)c;
	(void)tid; (void)nth;
	for (;;) {
		int i = __atomic_fetch_add(&p->next_row, 1, __ATOMIC_RELAXED);
		if (i >= (int)p->N) break;
		for (uint32_t j = (uint32_t)i + 1; j < p->N; j++)        /* cell (i,j): this thread only */
			p->inter[(size_t)i * p->N + j] +=
				and_popcount(p->bits + (size_t)i * p->words,
				             p->bits + (size_t)j * p->words, p->words);
	}
}

typedef struct { uint64_t *dst; uint64_t **src; int ns; size_t n; } sum_ctx;

static void sum_matrices(void *c, int tid, int nth)
{
	sum_ctx *s = (sum_ctx *)c;
	size_t lo = s->n * (size_t)tid / (size_t)nth, hi = s->n * (size_t)(tid + 1) / (size_t)nth;
	for (size_t x = lo; x < hi; x++) {
		uint64_t v = 0;
		for (int t = 0; t < s->ns; t++) v += s->src[t][x];
		s->dst[x] += v;
	}
}

/* ------------------------------------------------------------------- main */

static void usage(int status)
{
	fprintf(stderr,
"Usage: genome_compare_presence -A <genome list or directory> [options]\n"
"\n"
"  -A <path>   REQUIRED. Either a directory holding .fna.gz files, or a text\n"
"              file listing one genome path per line; a listed path may itself\n"
"              be a directory (expanded, sorted). Blank lines and #comments are\n"
"              ignored. Accepts .fna/.fa/.fasta/.ffn, plain or gzipped.\n"
"  -k <int>    k-mer length, 1..%d (default %d; 51 works).\n"
"  -O <dir>    output directory (default ., created if missing).\n"
"  -n <name>   output basename (default \"%s\").\n"
"  -o <prefix> output path prefix; overrides -O/-n.\n"
"  -z <level>  compress the presence file at this level (.zst 1..22 if built\n"
"              with libzstd, otherwise .gz 1..9). Default 0 = plain TSV.\n"
"  -F          forward strand only; do not fold k-mers to their canonical\n"
"              (lexicographically larger) orientation.\n"
"  -m <int>    write only k-mers present in at least this many strains\n"
"              (default 1). Filters the presence FILE only.\n"
"  -P <int>    write only k-mers present in at most this many strains\n"
"              (default 0 = no cap). Filters the presence FILE only; the\n"
"              identity matrix always uses every k-mer.\n"
"  -U          do not sort the presence file by k-mer (faster, hash order).\n"
"  -x <int>    strain-count crossover between the direct and bitset paths for\n"
"              the intersection sweep (default auto = N/8). Tuning only; both\n"
"              paths give identical numbers. -x 0 forces all-bitset.\n"
"  -t <int>    threads (default 4). Output is identical for any value.\n"
"  -h          this message.\n"
"  (renamed: -B -> -A, -M -> -P, -z now takes a level; -s = -k)\n"
"\n"
"Outputs:\n"
"  <prefix>.presence.tsv[.gz|.zst]  #kmer, n_strains, list_strain_ids\n"
"  <prefix>.manifest.tsv            strain_id, sample_id, filename, path,\n"
"                                   n_contigs, n_bases, n_kmer_positions,\n"
"                                   n_kmers_distinct\n"
"  <prefix>.identity.tsv            square Jaccard matrix (rows = strain_id)\n"
"  <prefix>.pairs.tsv               strain_a, strain_b, n_shared, n_union,\n"
"                                   jaccard, containment_a_in_b,\n"
"                                   containment_b_in_a\n"
"\n"
"Load the presence file with polars:\n"
"  pl.scan_csv('presence.tsv', separator='\\t')\n"
"    .with_columns(pl.col('list_strain_ids').str.split(','))\n",
	MAX_K, DEFAULT_K, DEFAULT_PREFIX);
	exit(status);
}

int main(int argc, char *argv[])
{
	const char *b_arg = NULL;
	const char *out_dir = NULL, *name = DEFAULT_PREFIX, *o_arg = NULL;
	int k = DEFAULT_K;
	int canonical = 1, compress = 0, sorted = 1;
	long min_strains = 1, max_strains = 0, crossover = -1;
	int nthreads = 4;
	int c;

	while ((c = getopt(argc, argv, "A:k:s:O:n:o:m:P:x:z:t:FUhBM")) != -1) {
		switch (c) {
		case 'A': b_arg = optarg; break;
		case 'k':
		case 's': k = atoi(optarg); break;      /* -s: same letter the other tools use */
		case 'O': out_dir = optarg; break;
		case 'n': name = optarg; break;
		case 'o': o_arg = optarg; break;
		case 'm': min_strains = atol(optarg); break;
		case 'P': max_strains = atol(optarg); break;
		case 'x': crossover = atol(optarg); break;
		case 't': nthreads = atoi(optarg); break;
		case 'z': {
			long lvl = cli_parse_nonneg("genome_compare_presence", 'z', optarg);
			compress = lvl > 0;
			if (lvl > 0) g_zlevel = (int)lvl;
			break;
		}
		case 'F': canonical = 0; break;
		case 'U': sorted = 0; break;
		case 'B': cli_renamed_flag("genome_compare_presence", 'B', "-A"); break;
		case 'M': cli_renamed_flag("genome_compare_presence", 'M', "-P"); break;
		case 'h': usage(0); break;
		default:  usage(1);
		}
	}

	if (!b_arg) usage(1);
	if (max_strains < 0) die("-P must be >= 0", NULL);
#ifdef HAVE_ZSTD
	if (compress && g_zlevel > 22) die("-z must be 0..22 for .zst output", NULL);
#else
	if (compress && g_zlevel > 9)  die("-z must be 0..9 for .gz output", NULL);
#endif
	/* <O>/<n>, or -o verbatim; creates the directory */
	char *prefix = cli_resolve_output("genome_compare_presence", out_dir,
	                                  name, o_arg, "");
	if (k < 1 || k > MAX_K)
		die("-k must be between 1 and 64", NULL);
	if (min_strains < 1) min_strains = 1;
	if (nthreads < 1) nthreads = 1;

	init_codes();
	double t0 = now_s();

	pathvec files = collect_inputs(b_arg);
	if (files.n == 0) die("no genome files found for -A", b_arg);
	if (files.n > UINT32_MAX) die("too many genomes", NULL);
	const uint32_t N = (uint32_t)files.n;

	fprintf(stderr, "[presence] %u genomes, k=%d, %s orientation\n",
	        N, k, canonical ? "canonical" : "forward-only");

	/* ---- pass over the genomes: build the inverted index -------------- */
	/* sorted output (default): sharded parallel build; -U (hash order) or one
	   thread: one table filled in the original order, so its layout is the same */
	int ntab = (sorted && nthreads > 1) ? NSHARD : 1;
	ktable *tabs = (ktable *)xcalloc((size_t)ntab, sizeof *tabs);
	for (int x = 0; x < ntab; x++)
		table_init(&tabs[x], ntab == 1 ? INITIAL_TABLE_SLOTS : INITIAL_TABLE_SLOTS / 16);
	gstat *gs = (gstat *)xcalloc(N, sizeof *gs);

	if (ntab == 1) build_index(&files, k, canonical, N, &tabs[0], gs, nthreads);
	else           build_index_sharded(&files, k, canonical, N, tabs, gs, nthreads);
	ktable table;                       /* totals, for the log line */
	memset(&table, 0, sizeof table);
	for (int x = 0; x < ntab; x++) { table.used += tabs[x].used; table.cap += tabs[x].cap; }

	double t1 = now_s();
	fprintf(stderr, "[presence] index: %" PRIu64 " distinct %d-mers "
	        "(%.1f MB hash) [%.2fs]\n",
	        table.used, k,
	        (double)table.cap * sizeof(kentry) / 1e6, t1 - t0);

	/* ---- manifest ----------------------------------------------------- */
	char path[4096];
	snprintf(path, sizeof path, "%s.manifest.tsv", prefix);
	FILE *mf = fopen(path, "w");
	if (!mf) die("cannot open manifest", path);
	fprintf(mf, "strain_id\tsample_id\tfilename\tpath\tn_contigs\tn_bases\t"
	            "n_kmer_positions\tn_kmers_distinct\n");
	for (uint32_t i = 0; i < N; i++) {
		char *sid = sample_id_of(files.path[i]);
		fprintf(mf, "%u\t%s\t%s\t%s\t%u\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\n",
		        i, sid, base_name(files.path[i]), files.path[i],
		        gs[i].contigs, gs[i].bases, gs[i].positions, gs[i].distinct);
		if (gs[i].short_contigs)
			fprintf(stderr, "warning: %s has %u contig(s) shorter than k=%d\n",
			        base_name(files.path[i]), gs[i].short_contigs, k);
		if (gs[i].distinct == 0)
			fprintf(stderr, "warning: %s contributed no k-mers\n",
			        base_name(files.path[i]));
		free(sid);
	}
	fclose(mf);
	fprintf(stderr, "[presence] wrote %s.manifest.tsv\n", prefix);

	/* ---- presence file + intersections ------------------------------- */
	uint64_t *inter = NULL;
	uint64_t written = 0, singletons = 0, core = 0;
	size_t nent = 0;
	{
	/* ---- entry list in output order (k-mer order, or slot order with -U) */
	nent = (size_t)table.used;
	skent *sk = (skent *)xmalloc((nent ? nent : 1) * sizeof *sk);
	uint64_t n_dense = 0;
	uint64_t thresh;
	if (crossover >= 0)      thresh = (uint64_t)crossover;
	else if (N < 64)         thresh = N;              /* never dense for tiny N */
	else                     thresh = N / 8 < 16 ? 16 : N / 8;
	{
		size_t m = 0;
		for (int x = 0; x < ntab; x++)
			for (uint64_t i = 0; i < tabs[x].cap; i++)
				if (tabs[x].slot[i].n != 0) {
					sk[m].key = (uint64_t)tabs[x].slot[i].key;
					sk[m].e   = &tabs[x].slot[i];
					if (tabs[x].slot[i].n > thresh) n_dense++;
					m++;
				}
	}
	if (sorted) {
		if (k <= 32) {
			sk = parallel_sort(sk, nent, 2 * k, nthreads);
		} else {                       /* 128-bit keys */
			qsort(sk, nent, sizeof *sk, cmp_skent128);
		}
	}

	uint64_t words = (n_dense + 63) / 64;
	uint64_t *bits = NULL;
	if (n_dense) {
		double gb = (double)N * (double)words * 8.0 / 1e9;
		if (gb > 64.0)
			fprintf(stderr, "warning: bitset needs %.1f GB; raise -x to push "
			        "more k-mers through the direct path\n", gb);
		bits = (uint64_t *)xcalloc((size_t)N * words, sizeof *bits);
		fprintf(stderr, "[presence] %" PRIu64 " k-mers in >%" PRIu64
		        " strains -> bitset (%.1f MB), rest direct\n",
		        n_dense, thresh, (double)N * words * 8 / 1e6);
	}

	inter = (uint64_t *)xcalloc((size_t)N * N, sizeof *inter);

	/* private intersection matrices: only while they stay small */
	int sth = nthreads;
	while (sth > 1 && (double)N * N * 8.0 * sth > 2e9) sth--;
	uint64_t **inter_t = (uint64_t **)xmalloc(sizeof(uint64_t *) * (size_t)sth);
	inter_t[0] = inter;
	for (int t = 1; t < sth; t++) inter_t[t] = (uint64_t *)xcalloc((size_t)N * N, sizeof(uint64_t));

	snprintf(path, sizeof path, "%s.presence.tsv%s", prefix,
	         compress ?
#ifdef HAVE_ZSTD
	         ".zst"
#else
	         ".gz"
#endif
	         : "");
	ow *pw = ow_open(path);
	const char *phdr = "#kmer\tn_strains\tlist_strain_ids\n";
	ow_puts(pw, phdr, strlen(phdr));

	sweep_ctx sc;
	memset(&sc, 0, sizeof sc);
	sc.sk = sk; sc.k = k; sc.N = N; sc.thresh = thresh; sc.words = words;
	sc.bits = bits; sc.min_strains = min_strains; sc.max_strains = max_strains;
	sc.inter_t = inter_t;
	sc.out = (pbuf *)xcalloc((size_t)sth, sizeof(pbuf));
	sc.dcount = (uint64_t *)xcalloc((size_t)sth, sizeof(uint64_t));
	sc.dbase  = (uint64_t *)xcalloc((size_t)sth, sizeof(uint64_t));
	sc.singletons = (uint64_t *)xcalloc((size_t)sth, sizeof(uint64_t));
	sc.core       = (uint64_t *)xcalloc((size_t)sth, sizeof(uint64_t));
	sc.written    = (uint64_t *)xcalloc((size_t)sth, sizeof(uint64_t));
	const size_t ROUND = (size_t)sth << 18;
	uint64_t dense_col = 0;
	for (size_t r = 0; r < nent; r += ROUND) {
		sc.lo = r; sc.hi = r + ROUND < nent ? r + ROUND : nent;
		if (n_dense) {
			par_run(sth, sweep_count_dense, &sc);
			for (int t = 0; t < sth; t++) { sc.dbase[t] = dense_col; dense_col += sc.dcount[t]; }
		}
		par_run(sth, sweep_part, &sc);
		for (int t = 0; t < sth; t++) ow_puts(pw, sc.out[t].buf ? sc.out[t].buf : "", sc.out[t].len);
	}
	for (int t = 0; t < sth; t++) {
		written += sc.written[t]; singletons += sc.singletons[t]; core += sc.core[t];
		free(sc.out[t].buf);
	}
	ow_close(pw);
	fprintf(stderr, "[presence] wrote %s (%" PRIu64 " of %zu k-mers)\n",
	        path, written, nent);

	if (sth > 1) {
		sum_ctx su = { inter, inter_t + 1, sth - 1, (size_t)N * N };
		par_run(nthreads, sum_matrices, &su);
		for (int t = 1; t < sth; t++) free(inter_t[t]);
	}
	free(inter_t);
	free(sc.out); free(sc.dcount); free(sc.dbase); free(sc.singletons); free(sc.core); free(sc.written);

	if (n_dense) {
		pop_ctx pc = { inter, bits, words, N, 0 };
		par_run(nthreads, popcount_rows, &pc);
		free(bits);
	}
	free(sk);
	}

	double t2 = now_s();
	fprintf(stderr, "[presence] intersections done [%.2fs]\n", t2 - t1);

	/* ---- identity matrix + long-form pairs ---------------------------- */
	char **sid = (char **)xmalloc((size_t)N * sizeof *sid);
	for (uint32_t i = 0; i < N; i++) sid[i] = sample_id_of(files.path[i]);

	snprintf(path, sizeof path, "%s.identity.tsv", prefix);
	FILE *im = fopen(path, "w");
	if (!im) die("cannot open identity matrix", path);
	fprintf(im, "strain_id\tsample_id");
	for (uint32_t j = 0; j < N; j++) fprintf(im, "\t%s", sid[j]);
	fprintf(im, "\n");

	for (uint32_t i = 0; i < N; i++) {
		fprintf(im, "%u\t%s", i, sid[i]);
		for (uint32_t j = 0; j < N; j++) {
			double v;
			if (i == j) {
				v = gs[i].distinct ? 1.0 : 0.0;
			} else {
				uint32_t a = i < j ? i : j, b = i < j ? j : i;
				uint64_t sh = inter[(size_t)a * N + b];
				uint64_t un = gs[i].distinct + gs[j].distinct - sh;
				v = un ? (double)sh / (double)un : 0.0;
			}
			fprintf(im, "\t%.6f", v);
		}
		fprintf(im, "\n");
	}
	fclose(im);
	fprintf(stderr, "[presence] wrote %s.identity.tsv\n", prefix);

	snprintf(path, sizeof path, "%s.pairs.tsv", prefix);
	FILE *pf = fopen(path, "w");
	if (!pf) die("cannot open pairs file", path);
	fprintf(pf, "strain_a\tstrain_b\tsample_a\tsample_b\tn_kmers_a\tn_kmers_b\t"
	            "n_shared\tn_union\tjaccard\tcontainment_a_in_b\t"
	            "containment_b_in_a\n");
	for (uint32_t i = 0; i < N; i++)
		for (uint32_t j = i + 1; j < N; j++) {
			uint64_t sh = inter[(size_t)i * N + j];
			uint64_t un = gs[i].distinct + gs[j].distinct - sh;
			fprintf(pf, "%u\t%u\t%s\t%s\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64
			            "\t%" PRIu64 "\t%.6f\t%.6f\t%.6f\n",
			        i, j, sid[i], sid[j], gs[i].distinct, gs[j].distinct,
			        sh, un,
			        un ? (double)sh / (double)un : 0.0,
			        gs[i].distinct ? (double)sh / (double)gs[i].distinct : 0.0,
			        gs[j].distinct ? (double)sh / (double)gs[j].distinct : 0.0);
		}
	fclose(pf);
	fprintf(stderr, "[presence] wrote %s.pairs.tsv\n", prefix);

	fprintf(stderr, "[presence] %" PRIu64 " k-mers in exactly 1 strain, "
	        "%" PRIu64 " in all %u; total %.2fs\n",
	        singletons, core, N, now_s() - t0);

	/* ---- cleanup ------------------------------------------------------ */
	for (uint32_t i = 0; i < N; i++) { free(sid[i]); free(files.path[i]); }
	free(sid);
	free(files.path);
	free(inter);
	free(gs);
	for (int x = 0; x < ntab; x++) table_free(&tabs[x]);
	free(tabs);
	free(prefix);
	return 0;
}
