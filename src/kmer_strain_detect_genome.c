/*
 * kmer_strain_detect_genome.c
 *
 * Takes a kmer TSV (with a '#kmer' column) and a targets genome list,
 * counts how many times each kmer appears in each genome, and writes
 * a gzipped TSV: rows = kmers, columns = genome basenames.
 *
 * Performance design (output is identical to the original string-hash
 * implementation; see test/speedup/ for the differential + race tests):
 *   - k-mers are 2-bit packed (k <= 32) and rolled along each read; the
 *     canonical form is max(fw, rc), which equals the string code's
 *     strcmp-larger orientation. Windows touching IUPAC/U bytes are evaluated
 *     with the original string logic, so exotic input stays exact.
 *   - lookups: L2-resident bit prefilter, then an open-addressing table of
 *     {key,row} slots; survivors are batched and prefetched.
 *   - query sets that cannot be 2-bit encoded (k > 32, lowercase/N in a
 *     query, mixed lengths) use the original string path for the whole run.
 *
 * Thread safety: each genome is one job, and each job writes only its own
 * cache-line-aligned count column plus its own total_kmer_counts slot
 * (written once at the end). The k-mer index and hash are read-only during
 * scanning. Workers are joined before the pool is destroyed.
 *
 * Usage:
 *   kmer_strain_detect_genome -k <kmer.tsv[.gz]> -B <targets.txt> -o <out.kmer_hits.tsv.gz> [-t threads] [-z level]
 */

#include "zlib.h"
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <inttypes.h>
#include <pthread.h>
#include <time.h>
#include <ctype.h>
#include "kseq.h"
#include "BIO_sequence.h"
#include "BIO_hash.h"
#include "genome_compare.h"
#include "cli_common.h"

#define NOT_PAIRED_END         0
#define MAX_LINE           65536
#define MAX_KMER_LEN         512

KSEQ_INIT(gzFile, gzread)

/* ── data structures ───────────────────────────────────────────────── */

typedef struct {
    char *basename;
    char *file1;
    char *file2;   /* NULL for SE / PEI */
    int   type;
} Metagenome;

typedef struct {
    int        mg_idx;
    Metagenome *mg;
    BIO_hash   h;
    const struct kmer_index_s *ki;
    unsigned int *col;               /* this metagenome's private count column */
    int        seed;
    uint64_t  *total_kmer_counts;
} ksd_job;

typedef struct {
    ksd_job       **queue;
    int            head, tail, size, capacity;
    int            active, shutdown;
    pthread_mutex_t mutex;
    pthread_cond_t  has_work;
    pthread_cond_t  all_done;
    pthread_t      *threads;
    int             nthreads;
} ksd_pool;

/* ── forward declarations ──────────────────────────────────────────── */

static char     *strip_basename(const char *path);
static int       parse_targets(const char *file, Metagenome **out);
static int       parse_kmer_tsv(const char *file, BIO_hash h, int n_mg,
                                 int *seed_out, char ***kmer_list_out);
static void     *ksd_worker(void *arg);
static ksd_pool *ksd_pool_create(int nthreads, int capacity);
static void      ksd_pool_submit(ksd_pool *pool, ksd_job *job);
static void      ksd_pool_wait_and_destroy(ksd_pool *pool);
static void      usage(void);

/* ── file type helpers ─────────────────────────────────────────────── */

static char *strip_basename(const char *path) {
    const char *base = strrchr(path, '/');
    base = base ? base + 1 : path;
    char *b = strdup(base);
    const char *exts[] = {
        ".fastq.gz", ".fasta.gz", ".fa.gz", ".fq.gz", ".fna.gz",
        ".fastq",    ".fasta",    ".fa",    ".fq",    ".fna",
        NULL
    };
    for (int i = 0; exts[i]; i++) {
        size_t elen = strlen(exts[i]);
        size_t blen = strlen(b);
        if (blen > elen && strcmp(b + blen - elen, exts[i]) == 0) {
            b[blen - elen] = '\0';
            break;
        }
    }
    return b;
}

/* ── targets parser ────────────────────────────────────────────────── */

static int valid_genome_ext(const char *path) {
    const char *exts[] = { ".fna", ".fna.gz", ".fasta", ".fasta.gz", NULL };
    size_t plen = strlen(path);
    for (int i = 0; exts[i]; i++) {
        size_t elen = strlen(exts[i]);
        if (plen >= elen && strcmp(path + plen - elen, exts[i]) == 0) return 1;
    }
    return 0;
}

static int parse_targets(const char *file, Metagenome **out) {
    FILE *fp = fopen(file, "r");
    if (!fp) {
        fprintf(stderr, "cannot open targets file %s: %s\n", file, strerror(errno));
        exit(1);
    }

    int cap = 64, n = 0;
    Metagenome *mgs = malloc(cap * sizeof(Metagenome));

    char *line = NULL;
    size_t len = 0;
    while (getline(&line, &len, fp) != -1) {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
        char *nl = strchr(line, '\n'); if (nl) *nl = '\0';
        char *cr = strchr(line, '\r'); if (cr) *cr = '\0';
        if (line[0] == '\0') continue;
        if (!valid_genome_ext(line)) {
            fprintf(stderr, "skipping unsupported file type: %s\n", line);
            continue;
        }

        if (n == cap) { cap *= 2; mgs = realloc(mgs, cap * sizeof(Metagenome)); }
        mgs[n].type     = NOT_PAIRED_END;
        mgs[n].file1    = strdup(line);
        mgs[n].file2    = NULL;
        mgs[n].basename = strip_basename(line);
        n++;
    }
    free(line);
    fclose(fp);
    *out = mgs;
    return n;
}

/* ── kmer TSV parser ───────────────────────────────────────────────── */

static int parse_kmer_tsv(const char *file, BIO_hash h, int n_mg,
                           int *seed_out, char ***kmer_list_out) {
    gzFile fp = gzopen(file, "r");
    if (!fp) {
        fprintf(stderr, "cannot open kmer file %s: %s\n", file, strerror(errno));
        exit(1);
    }

    char line[MAX_LINE];
    int kmer_col = -1;
    int seed = -1;

    /* find the '#kmer' column index from the header */
    if (!gzgets(fp, line, MAX_LINE)) {
        fprintf(stderr, "kmer file is empty: %s\n", file);
        exit(1);
    }
    {
        char hdr[MAX_LINE];
        strncpy(hdr, line, MAX_LINE);
        hdr[MAX_LINE - 1] = '\0';
        char *tok = strtok(hdr, "\t\n\r");
        int col = 0;
        while (tok) {
            if (strcmp(tok, "#kmer") == 0) { kmer_col = col; break; }
            tok = strtok(NULL, "\t\n\r");
            col++;
        }
    }
    if (kmer_col < 0) {
        fprintf(stderr, "no '#kmer' column found in %s\n", file);
        exit(1);
    }

    int cap = 4096, n = 0;
    char **kmer_list = malloc(cap * sizeof(char *));
    char kmer_rc[MAX_KMER_LEN];

    while (gzgets(fp, line, MAX_LINE)) {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
        char *nl = strchr(line, '\n'); if (nl) *nl = '\0';
        char *cr = strchr(line, '\r'); if (cr) *cr = '\0';

        /* advance to the kmer column */
        char *tok = strtok(line, "\t");
        for (int c = 0; c < kmer_col && tok; c++)
            tok = strtok(NULL, "\t");
        if (!tok || tok[0] == '\0') continue;

        if (seed < 0) {
            seed = (int)strlen(tok);
            *seed_out = seed;
        }

        /* compute canonical (oriented) form */
        strncpy(kmer_rc, tok, MAX_KMER_LEN - 1);
        kmer_rc[MAX_KMER_LEN - 1] = '\0';
        BIO_reverseComplement(kmer_rc);
        char *oriented = orient_string(tok, kmer_rc, seed);

        /* skip duplicates (same canonical kmer) */
        if (BIO_searchHash(h, oriented) != NULL) continue;

        (void)n_mg;
        unsigned int *row = malloc(sizeof(unsigned int));   /* row index = output order */
        *row = (unsigned int)n;
        BIO_addHashData(h, oriented, row);

        if (n == cap) { cap *= 2; kmer_list = realloc(kmer_list, cap * sizeof(char *)); }
        kmer_list[n++] = strdup(tok);   /* store original for output */
    }

    gzclose(fp);
    *kmer_list_out = kmer_list;
    return n;
}

/* ── 2-bit k-mer index ─────────────────────────────────────────────── */
/*
 * Fast path (k <= 32, every query k-mer is uppercase ACGT of length k):
 * k-mers are packed 2 bits/base (A=0 C=1 G=2 T=3) so lexicographic order ==
 * numeric order and the canonical form used by the string code
 * (strcmp(fw, rc) > 0 ? fw : rc, i.e. the lexicographic max) is max(fw, rc).
 * Lookups go to an open-addressing table keyed by that uint64.
 *
 * Exactness: a read containing any byte outside ACGTN/acgtn (IUPAC codes, U,
 * other junk) is scanned with the original string loop, because the string
 * code's COMPLEMENT table maps U->A, so such windows can still hit. Reads that
 * are pure ACGTN take the 2-bit loop, where windows with N are skipped exactly
 * like contains_N() did.
 */
/* 0-3 ACGT/acgt, 4 = N/n, 5 = other printable ASCII (IUPAC, U, ...),
   6 = NUL / non-ASCII (whole read goes through the original string loop) */
static uint8_t NT2BIT[256];

static void init_nt2bit(void) {
    memset(NT2BIT, 5, 128);
    memset(NT2BIT + 128, 6, 128);
    NT2BIT[0] = 6;
    NT2BIT['A'] = NT2BIT['a'] = 0; NT2BIT['C'] = NT2BIT['c'] = 1;
    NT2BIT['G'] = NT2BIT['g'] = 2; NT2BIT['T'] = NT2BIT['t'] = 3;
    NT2BIT['N'] = NT2BIT['n'] = 4;
}

typedef struct {                     /* one 16-byte slot: key and row together */
    uint64_t key;
    uint32_t row1;                   /* row index + 1; 0 = empty slot */
} kidx_slot;

typedef struct kmer_index_s {
    kidx_slot *slots;
    uint64_t   mask;                 /* capacity - 1 (power of two) */
    /* prefilter: one bit per hash value, sized to stay cache-resident
       (>= 16 bits per k-mer, at most 2^24 bits = 2 MB). A clear bit proves
       the k-mer is absent, so most read k-mers never touch the table. */
    uint64_t  *bits;
    int        bshift;               /* 64 - log2(number of bits) */
} kmer_index;

static inline uint64_t mix64(uint64_t x) {   /* murmur3 fmix64 */
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33; return x;
}

static kmer_index *kidx_create(size_t n) {
    size_t cap = 1024;
    while (cap < 2 * n) cap <<= 1;               /* load factor <= 0.5 */
    kmer_index *ki = malloc(sizeof *ki);
    ki->slots = calloc(cap, sizeof(kidx_slot));
    ki->mask = cap - 1;
    int lb = 16;
    while (lb < 24 && ((size_t)1 << lb) < 16 * n) lb++;
    ki->bits   = calloc(((size_t)1 << lb) / 64, sizeof(uint64_t));
    ki->bshift = 64 - lb;
    return ki;
}

static void kidx_put(kmer_index *ki, uint64_t key, uint32_t row) {
    uint64_t hv = mix64(key), b = hv >> ki->bshift;
    ki->bits[b >> 6] |= 1ULL << (b & 63);
    uint64_t i = hv & ki->mask;
    while (ki->slots[i].row1 && ki->slots[i].key != key) i = (i + 1) & ki->mask;
    ki->slots[i].key = key; ki->slots[i].row1 = row + 1;
}

/* batched lookup, stage 1: prefilter; on pass, prefetch the slot */
static inline int kidx_probe(const kmer_index *ki, uint64_t key, uint64_t *hv_out) {
    uint64_t hv = mix64(key), b = hv >> ki->bshift;
    if (!(ki->bits[b >> 6] & (1ULL << (b & 63)))) return 0;
    __builtin_prefetch(&ki->slots[hv & ki->mask]);
    *hv_out = hv;
    return 1;
}

/* batched lookup, stage 2: resolve a key that passed the prefilter */
static inline uint32_t kidx_resolve(const kmer_index *ki, uint64_t key, uint64_t hv) {
    uint64_t i = hv & ki->mask;
    for (;;) {
        const kidx_slot *sl = &ki->slots[i];
        if (!sl->row1) return 0;
        if (sl->key == key) return sl->row1;
        i = (i + 1) & ki->mask;
    }
}

#define KSD_BATCH 32

static void kidx_destroy(kmer_index *ki) {
    if (!ki) return;
    free(ki->slots); free(ki->bits); free(ki);
}

/* encode an uppercase ACGT string of length k; returns 0 if not encodable */
static int encode_kmer(const char *s, int k, uint64_t *out) {
    uint64_t v = 0;
    for (int i = 0; i < k; i++) {
        uint8_t c = NT2BIT[(unsigned char)s[i]];
        if (c > 3 || s[i] > 'Z') return 0;   /* lowercase is NOT the same key in the string code */
        v = (v << 2) | c;
    }
    *out = v; return 1;
}

/* ── read scanner ──────────────────────────────────────────────────── */

/* original string-based loop for one read (exact reference semantics) */
static void scan_read_legacy(kseq_t *seq, char **rcbuf, size_t *rc_cap,
                             BIO_hash h, int seed, unsigned int *col) {
    unsigned int i;
    BIO_stringToUpper(seq->seq.s);

    if (seq->seq.l + 1 > *rc_cap) {
        *rc_cap = seq->seq.l + 1;
        *rcbuf = realloc(*rcbuf, *rc_cap);
    }
    char *seedStrRevComp = *rcbuf;
    strcpy(seedStrRevComp, seq->seq.s);
    BIO_reverseComplement(seedStrRevComp);

    char *seed_seq    = seq->seq.s;
    char *seed_seq_rc = &seedStrRevComp[seq->seq.l - seed];
    int   has_N       = contains_N(seed_seq);

    for (i = 0; i < seq->seq.l - (unsigned int)seed + 1; i++) {
        char temp_nuc     = seed_seq[seed];
        seed_seq[seed]    = '\0';
        seed_seq_rc[seed] = '\0';

        char *orientStr = (strcmp(seed_seq, seed_seq_rc) > 0)
                          ? seed_seq : seed_seq_rc;

        if (!has_N || !contains_N(orientStr)) {
            unsigned int *row = (unsigned int *)BIO_searchHash(h, orientStr);
            if (row) col[*row]++;         /* col is private to this job */
        }

        seed_seq[seed] = temp_nuc;
        seed_seq++;
        seed_seq_rc--;
    }
}

/* one window [p, p+k) containing an exotic ASCII byte, evaluated exactly like
   the string loop: uppercase, reverse complement via COMPLEMENT[], keep the
   strcmp-larger orientation, skip if it contains N, look it up by string */
static void legacy_window(const unsigned char *p, int k, char *fw, char *rc,
                          BIO_hash h, unsigned int *col) {
    for (int i = 0; i < k; i++) fw[i] = (char)toupper(p[i]);
    fw[k] = '\0';
    for (int i = 0; i < k; i++) rc[i] = COMPLEMENT[(int)fw[k - 1 - i]];
    rc[k] = '\0';
    char *o = (strcmp(fw, rc) > 0) ? fw : rc;
    if (contains_N(o)) return;
    unsigned int *row = (unsigned int *)BIO_searchHash(h, o);
    if (row) col[*row]++;
}

static void scan_file(const char *filepath, BIO_hash h, const kmer_index *ki,
                      int seed, int mg_idx, unsigned int *col,
                      uint64_t *total_kmer_counts) {
    gzFile fp = gzopen(filepath, "r");
    if (!fp) {
        fprintf(stderr, "cannot open %s: %s\n", filepath, strerror(errno));
        return;
    }

    kseq_t *seq = kseq_init(fp);
    uint64_t n_positions = 0;   /* thread-local; written to the shared array once */
    size_t rc_cap = 0;
    char *rcbuf = NULL;

    const int      k     = seed;
    const uint64_t kmask = (k == 32) ? ~0ULL : ((1ULL << (2 * k)) - 1);
    const int      rcsh  = 2 * (k - 1);

    while (kseq_read(seq) >= 0) {
        if ((int)seq->seq.l < seed) continue;
        n_positions += seq->seq.l - (unsigned int)seed + 1;

        const unsigned char *s = (const unsigned char *)seq->seq.s;
        const size_t L = seq->seq.l;

        int plain = (ki != NULL);
        if (plain) {
            for (size_t j = 0; j < L; j++)
                if (NT2BIT[s[j]] == 6) { plain = 0; break; }
        }
        if (!plain) {               /* NUL/non-ASCII bytes or no fast index */
            scan_read_legacy(seq, &rcbuf, &rc_cap, h, seed, col);
            continue;
        }

        /* 2-bit rolling scan. N resets the window (exactly the contains_N
           skip). An exotic ASCII byte poisons the next k windows, which are
           evaluated with the string logic instead (legacy_window). */
        uint64_t fw = 0, rc = 0;
        int valid = 0;                /* bases since last N */
        size_t exotic_until = 0;      /* windows ending before this index are exotic */
        uint64_t pk[KSD_BATCH], ph[KSD_BATCH];   /* keys that passed the filter */
        int np = 0;
        for (size_t j = 0; j < L; j++) {
            uint8_t c = NT2BIT[s[j]];
            if (c == 4) { valid = 0; fw = rc = 0; continue; }     /* N */
            if (c == 5) { exotic_until = j + (size_t)k; c = 0; }  /* code irrelevant */
            fw = ((fw << 2) | c) & kmask;
            rc = (rc >> 2) | ((uint64_t)(3 - c) << rcsh);
            if (++valid >= k) {
                if (j < exotic_until) {
                    if (rc_cap < (size_t)(2 * k + 2)) {
                        rc_cap = (size_t)(2 * k + 2);
                        rcbuf = realloc(rcbuf, rc_cap);
                    }
                    legacy_window(s + j + 1 - k, k, rcbuf, rcbuf + k + 1, h, col);
                    continue;
                }
                uint64_t key = fw > rc ? fw : rc;
                if (kidx_probe(ki, key, &ph[np])) {
                    pk[np++] = key;
                    if (np == KSD_BATCH) {    /* slots were prefetched while we kept scanning */
                        for (int q = 0; q < np; q++) {
                            uint32_t r1 = kidx_resolve(ki, pk[q], ph[q]);
                            if (r1) col[r1 - 1]++;   /* col is private to this job */
                        }
                        np = 0;
                    }
                }
            }
        }
        for (int q = 0; q < np; q++) {
            uint32_t r1 = kidx_resolve(ki, pk[q], ph[q]);
            if (r1) col[r1 - 1]++;
        }
    }

    total_kmer_counts[mg_idx] = n_positions;
    free(rcbuf);
    kseq_destroy(seq);
    gzclose(fp);
}

/* ── thread pool ───────────────────────────────────────────────────── */

static void *ksd_worker(void *arg) {
    ksd_pool *pool = (ksd_pool *)arg;
    for (;;) {
        pthread_mutex_lock(&pool->mutex);
        while (pool->size == 0 && !pool->shutdown)
            pthread_cond_wait(&pool->has_work, &pool->mutex);
        if (pool->shutdown && pool->size == 0) {
            pthread_mutex_unlock(&pool->mutex);
            break;
        }
        ksd_job *job = pool->queue[pool->head];
        pool->head = (pool->head + 1) % pool->capacity;
        pool->size--;
        pool->active++;
        pthread_mutex_unlock(&pool->mutex);

        scan_file(job->mg->file1, job->h, job->ki, job->seed, job->mg_idx, job->col, job->total_kmer_counts);
        free(job);

        pthread_mutex_lock(&pool->mutex);
        pool->active--;
        if (pool->size == 0 && pool->active == 0)
            pthread_cond_signal(&pool->all_done);
        pthread_mutex_unlock(&pool->mutex);
    }
    return NULL;
}

static ksd_pool *ksd_pool_create(int nthreads, int capacity) {
    ksd_pool *pool = malloc(sizeof(ksd_pool));
    pool->queue    = malloc(sizeof(ksd_job *) * capacity);
    pool->head = pool->tail = pool->size = pool->active = pool->shutdown = 0;
    pool->capacity = capacity;
    pthread_mutex_init(&pool->mutex, NULL);
    pthread_cond_init(&pool->has_work, NULL);
    pthread_cond_init(&pool->all_done, NULL);
    pool->threads  = malloc(sizeof(pthread_t) * nthreads);
    pool->nthreads = nthreads;
    for (int i = 0; i < nthreads; i++)
        pthread_create(&pool->threads[i], NULL, ksd_worker, pool);
    return pool;
}

static void ksd_pool_submit(ksd_pool *pool, ksd_job *job) {
    pthread_mutex_lock(&pool->mutex);
    pool->queue[pool->tail] = job;
    pool->tail = (pool->tail + 1) % pool->capacity;
    pool->size++;
    pthread_cond_signal(&pool->has_work);
    pthread_mutex_unlock(&pool->mutex);
}

static void ksd_pool_wait_and_destroy(ksd_pool *pool) {
    pthread_mutex_lock(&pool->mutex);
    while (pool->size > 0 || pool->active > 0)
        pthread_cond_wait(&pool->all_done, &pool->mutex);
    pool->shutdown = 1;
    pthread_cond_broadcast(&pool->has_work);
    pthread_mutex_unlock(&pool->mutex);
    /* join: every worker has left the pool before it is destroyed
       (the old 10 ms sleep could free the mutex under a running worker) */
    for (int i = 0; i < pool->nthreads; i++)
        pthread_join(pool->threads[i], NULL);
    free(pool->threads);
    pthread_mutex_destroy(&pool->mutex);
    pthread_cond_destroy(&pool->has_work);
    pthread_cond_destroy(&pool->all_done);
    free(pool->queue);
    free(pool);
}

/* ── usage ─────────────────────────────────────────────────────────── */

static void usage(void) {
    fprintf(stderr, "Usage: kmer_strain_detect_genome -k <kmer_tsv> -B <targets.txt> [-O <dir>] [-n <name> | -o <out.kmer_hits.tsv.gz>] [-G <background.txt>] [-t threads] [-z level]\n\n");
    fprintf(stderr, "  -k  kmer TSV with a '#kmer' column (plain or gzipped)\n");
    fprintf(stderr, "  -B  target genomes file (one genome path per line)\n");
    fprintf(stderr, "  -G  background metagenomes file (same format as -B); columns prefixed 'backmeta_'\n");
    fprintf(stderr, "      types: PE, SE, PEI\n");
    fprintf(stderr, "  -O  output directory (default: ., created if missing)\n");
    fprintf(stderr, "  -n  output basename (default: -k filename without extensions);\n");
    fprintf(stderr, "      writes <-O>/<-n>.kmer_hits.tsv.gz\n");
    fprintf(stderr, "  -o  output file; overrides -O/-n (e.g. sample.kmer_hits.tsv.gz)\n");
    fprintf(stderr, "  -t  threads (default: 4)   [was -j]\n");
    fprintf(stderr, "  -z  gzip level 1-9 (default: 9). Lower is much faster on large\n");
    fprintf(stderr, "      k-mer sets; the TSV content is identical, only the .gz bytes differ\n\n");
    fprintf(stderr, "Output: gzipped TSV — first column is #kmer, then target columns,\n");
    fprintf(stderr, "        then backmeta_* columns for background metagenomes.\n");
}

/* ── main ──────────────────────────────────────────────────────────── */

int main(int argc, char *argv[]) {
    char *kmer_file       = NULL;
    char *targets_file    = NULL;
    char *background_file = NULL;
    char *outfile         = NULL;
    char *out_dir         = NULL;
    char *out_name        = NULL;
    int   num_threads     = 4;
    int   zlevel          = 9;     /* default 9 = byte-identical to previous releases */
    int   c;

    while ((c = getopt(argc, argv, "k:B:G:O:n:o:t:j:z:h")) != EOF) {
        switch (c) {
            case 'k': kmer_file       = strdup(optarg); break;
            case 'B': targets_file    = strdup(optarg); break;
            case 'G': background_file = strdup(optarg); break;
            case 'O': out_dir         = strdup(optarg); break;
            case 'n': out_name        = strdup(optarg); break;
            case 'o': outfile         = strdup(optarg); break;
            case 't': num_threads     = atoi(optarg);   break;
            case 'z': zlevel          = atoi(optarg);   break;
            case 'j': cli_renamed_flag("kmer_strain_detect_genome", 'j', "-t"); break;
            case 'h': usage(); return 0;
            default:  usage(); return 1;
        }
    }

    if (!kmer_file || !targets_file) { usage(); return 1; }
    if (num_threads < 1) num_threads = 1;

    /* output: -o verbatim, else <-O>/<-n>.kmer_hits.tsv.gz with -n derived
       from the -k file name when not given */
    if (!outfile && !out_name) {
        out_name = cli_derive_basename(kmer_file);
        if (!out_name) {
            fprintf(stderr, "error: could not derive an output name from -k %s; "
                    "pass -n or -o\n", kmer_file);
            return 1;
        }
    }
    {
        char *resolved = cli_resolve_output("kmer_strain_detect_genome", out_dir, out_name,
                                            outfile, ".kmer_hits.tsv.gz");
        free(outfile);
        outfile = resolved;
    }
    fprintf(stderr, "output: %s\n", outfile);

    /* 1. parse metagenome targets */
    Metagenome *mgs;
    int n_mg = parse_targets(targets_file, &mgs);
    if (n_mg == 0) { fprintf(stderr, "no valid metagenomes in targets file\n"); return 1; }
    fprintf(stderr, "loaded %d metagenome(s)\n", n_mg);

    /* 1b. parse background metagenomes and append with backmeta_ prefix */
    int n_total = n_mg;
    if (background_file) {
        Metagenome *bg;
        int n_bg = parse_targets(background_file, &bg);
        fprintf(stderr, "loaded %d background metagenome(s)\n", n_bg);
        mgs = realloc(mgs, (n_mg + n_bg) * sizeof(Metagenome));
        for (int i = 0; i < n_bg; i++) {
            char *prefixed = malloc(strlen(bg[i].basename) + 10); /* "backmeta_" (9) + NUL */
            sprintf(prefixed, "backmeta_%s", bg[i].basename);
            free(bg[i].basename);
            bg[i].basename = prefixed;
            mgs[n_mg + i] = bg[i];
        }
        free(bg);
        n_total = n_mg + n_bg;
    }

    /* 2. load kmers into hash; count array has one slot per metagenome (target + background) */
    BIO_hash h = BIO_initHash(DEFAULT_GENOME_HASH_SIZE);
    char **kmer_list;
    int seed = 0, n_kmers;
    n_kmers = parse_kmer_tsv(kmer_file, h, n_total, &seed, &kmer_list);
    fprintf(stderr, "loaded %d kmer(s) (k=%d)\n", n_kmers, seed);

    /* 2b. 2-bit integer index (fast path) when every query k-mer is ACGT of length k<=32 */
    init_nt2bit();
    kmer_index *ki = NULL;
    if (seed >= 1 && seed <= 32) {
        ki = kidx_create((size_t)n_kmers);
        char rcb[MAX_KMER_LEN];
        for (int i = 0; i < n_kmers && ki; i++) {
            uint64_t code;
            if ((int)strlen(kmer_list[i]) != seed) { kidx_destroy(ki); ki = NULL; break; }
            strncpy(rcb, kmer_list[i], MAX_KMER_LEN - 1); rcb[MAX_KMER_LEN - 1] = '\0';
            BIO_reverseComplement(rcb);
            char *oriented = orient_string(kmer_list[i], rcb, seed);
            if (!encode_kmer(oriented, seed, &code)) { kidx_destroy(ki); ki = NULL; break; }
            kidx_put(ki, code, (uint32_t)i);   /* kmer_list[i] is row i */
        }
    }
    fprintf(stderr, "scan mode: %s\n", ki ? "2-bit integer index" : "string hash (legacy)");

    /* 3. scan all metagenomes in parallel */
    uint64_t *total_kmer_counts = calloc(n_total, sizeof(uint64_t));
    /* one count column per metagenome, each on its own cache lines:
       a worker only ever touches the column of the job it is running */
    size_t col_bytes = (((size_t)n_kmers * sizeof(unsigned int)) + 63) & ~(size_t)63;
    if (col_bytes == 0) col_bytes = 64;
    unsigned int **cols = malloc(n_total * sizeof(unsigned int *));
    for (int i = 0; i < n_total; i++) {
        cols[i] = aligned_alloc(64, col_bytes);
        memset(cols[i], 0, col_bytes);
    }
    ksd_pool *pool = ksd_pool_create(num_threads, n_total + 1);
    for (int i = 0; i < n_total; i++) {
        ksd_job *job = malloc(sizeof(ksd_job));
        job->mg_idx           = i;
        job->mg               = &mgs[i];
        job->h                = h;
        job->ki               = ki;
        job->col              = cols[i];
        job->seed             = seed;
        job->total_kmer_counts = total_kmer_counts;
        ksd_pool_submit(pool, job);
    }
    ksd_pool_wait_and_destroy(pool);
    fprintf(stderr, "done scanning\n");

    /* 4. write gzipped TSV output */
    if (zlevel < 1 || zlevel > 9) { fprintf(stderr, "-z must be 1..9\n"); return 1; }
    char zmode[4] = { 'w', 'b', (char)('0' + zlevel), '\0' };
    gzFile out = gzopen(outfile, zmode);
    if (!out) {
        fprintf(stderr, "cannot open output %s: %s\n", outfile, strerror(errno));
        return 1;
    }

    /* header row */
    gzprintf(out, "#kmer");
    for (int i = 0; i < n_total; i++) gzprintf(out, "\t%s", mgs[i].basename);
    gzprintf(out, "\n");

    /* data rows — in the same order as the input TSV.
       Each row is formatted into one buffer (same bytes gzprintf produced:
       kmer, then "\t%u" per column, then "\n") and handed to zlib in one call. */
    size_t rowcap = MAX_KMER_LEN + (size_t)n_total * 11 + 2;
    char  *rowbuf = malloc(rowcap);
    for (int i = 0; i < n_kmers; i++) {
        size_t klen = strlen(kmer_list[i]);
        if (klen + (size_t)n_total * 11 + 2 > rowcap) {
            rowcap = klen + (size_t)n_total * 11 + 2;
            rowbuf = realloc(rowbuf, rowcap);
        }
        char *w = rowbuf;
        memcpy(w, kmer_list[i], klen); w += klen;
        for (int j = 0; j < n_total; j++) {
            unsigned int v = cols[j][i];
            char tmp[10]; int n = 0;            /* UINT_MAX has 10 digits */
            do { tmp[n++] = (char)('0' + v % 10); v /= 10; } while (v);
            *w++ = '\t';
            while (n) *w++ = tmp[--n];
        }
        *w++ = '\n';
        if (gzwrite(out, rowbuf, (unsigned)(w - rowbuf)) != (int)(w - rowbuf)) {
            fprintf(stderr, "write error on %s\n", outfile);
            return 1;
        }
        free(kmer_list[i]);
    }
    free(rowbuf);

    /* summary row — total kmer positions evaluated per metagenome */
    gzprintf(out, "total_evaluated");
    for (int i = 0; i < n_total; i++)
        gzprintf(out, "\t%" PRIu64, total_kmer_counts[i]);
    gzprintf(out, "\n");

    gzclose(out);
    fprintf(stderr, "output written to %s\n", outfile);

    /* cleanup */
    kidx_destroy(ki);
    for (int i = 0; i < n_total; i++) free(cols[i]);
    free(cols);
    free(total_kmer_counts);
    free(kmer_list);
    BIO_destroyHashD(h);
    for (int i = 0; i < n_total; i++) {
        free(mgs[i].basename);
        free(mgs[i].file1);
        free(mgs[i].file2);
    }
    free(mgs);
    free(kmer_file);
    free(targets_file);
    free(background_file);
    free(outfile);
    free(out_dir);
    free(out_name);

    return 0;
}
