/*
 * kmer_strain_detect.c
 *
 * Takes a kmer file and a targets metagenome list, counts how many times each
 * kmer appears in each metagenome, and writes a gzipped TSV:
 * rows = kmers, columns = metagenome basenames.
 *
 * Accepted kmer file formats (-k, plain or gzipped):
 *   1. TSV with a header containing a '#kmer' (or 'kmer') column; all other
 *      columns are ignored (rare_kmers.mapped.tsv.gz, kmer_blocks route).
 *   2. Classic headerless list (kmer_scrub_filter.py output, *.scrubbed_kmers):
 *      kmer in the first whitespace-separated column. '#' lines (scrub stats)
 *      are skipped wherever they appear.
 *   An empty / comment-only file is not an error: it yields a header-only
 *   output with a zeroed total_evaluated row (everything was scrubbed away).
 *
 * Targets can also be per-sample count databases (.kmdb) built with
 *   kmer_metagenome_db -m sample -i R1.fq.gz -i R2.fq.gz -n <sample>
 * A .kmdb target is recognised by its extension (or type KMDB). Instead of
 * re-reading the sample, the query k-mers are packed to 2 bits, sorted, and
 * looked up block by block; blocks of the database that cannot contain a
 * query k-mer are never decompressed. Counts and total_evaluated are
 * identical to scanning the reads the database was built from (provided it
 * was built with -c 1, the sample-mode default). Read files and .kmdb
 * targets can be mixed in one targets file. Needs k <= 31 and the same k as
 * the database.
 *
 * total_evaluated row: the sum of the counts of ALL k-mers of the sample
 * (not only the query k-mers), i.e. the number of k-mer increments made.
 * That is the denominator for relative abundance. FASTQ targets: every
 * ACGT-only k-mer window of every read with length >= k, minus the windows
 * skipped by pair dedup. .kmdb targets: the database's total_evaluated
 * metadata (same definition; databases from before that key existed fall
 * back to valid_kmer_positions, their sum of counts). Windows holding any
 * non-ACGT base are not counted (until now they were).
 *
 * Pair dedup (--pair-dedup, default; --no-pair-dedup for the old counts):
 * for PE targets (two files, read in lockstep) and PEI targets
 * (consecutive records are mates) every distinct canonical k-mer is counted
 * at most once per read pair, whether it occurs in R1, R2 or both, and
 * however often it repeats inside the pair (sylph's rule). This applies to
 * every k-mer, so it also lowers total_evaluated. SE targets are unchanged.
 * Mates must be in sync (same number of records; names equal after
 * dropping a trailing /1 /2 .1 .2 _1 _2), otherwise the program stops.
 * Same rule as kmer_metagenome_db -m sample --pair-dedup, so FASTQ and
 * .kmdb targets built from the same reads give identical counts and
 * total_evaluated when both use the same setting.
 *
 * Each metagenome is processed in its own thread. Because each thread only
 * writes to its own column index in the per-kmer count array, no mutex is
 * needed for count updates. The hash is read-only during scanning.
 *
 * Usage:
 *   kmer_strain_detect -k <kmers[.gz]> -B <targets.txt> -o <out.kmer_hits.tsv.gz> [-j threads]
 *
 * Read hash (--read-hash): also writes <name>.kmer_readhash.tsv.gz, same rows
 * and columns as the hits file, holding per k-mer and sample the smallest
 * read (fragment) hash over the fragments that contain the k-mer (readhash.h;
 * 0 = not seen). Hit k-mers that share a value were (at count 1: exactly)
 * produced by the same read or read pair, so the number of distinct values
 * among a strain's hit k-mers estimates the reads behind the hits. FASTQ
 * targets: fragments are numbered while scanning (read pairs with pair dedup
 * and PE / PEI input, else single reads; only fragments with a read of
 * length >= k), exactly as kmer_metagenome_db -m sample --read-hash numbers
 * them, so both give identical hashes for the same reads. .kmdb targets: the
 * database's stored hashes (kmdb version 2); databases built without
 * --read-hash give 0 (with a warning). Last row n_fragments: fragments
 * numbered in the sample (FASTQ: counted here; .kmdb: its metadata, 0 if
 * missing).
 *
 * Build without libzstd (no -DHAVE_ZSTD): .kmdb targets are rejected with an
 * error, everything else works as before.
 */
#define _GNU_SOURCE
#include "zlib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <strings.h>
#include <errno.h>
#include <inttypes.h>
#include <pthread.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>
#include "kseq.h"
#include "BIO_sequence.h"
#include "BIO_hash.h"
#include "genome_compare.h"
#include "cli_common.h"
#include "readhash.h"
#ifdef HAVE_ZSTD
#include "kmdb.h"
#endif

#define NOT_PAIRED_END         0
#define IS_PAIRED_END          1
#define IS_PAIRED_END_INTERLEAVE 2
#define IS_KMDB                3   /* kmer_metagenome_db -m sample output */
#define UNKNOWN_FILE_TYPE     -1
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

/* query k-mer, 2-bit packed canonical key -> its count vector.
   With --read-hash the vector has 2 x n_total slots: counts[0..n_total) and
   min read hashes [n_total..2 n_total). */
typedef struct {
    uint64_t      key;
    unsigned int *counts;
} QKmer;

typedef struct {
    int        mg_idx;
    Metagenome *mg;
    BIO_hash   h;
    int        seed;
    uint64_t  *total_kmer_counts;
    const QKmer *qk;     /* sorted by key; NULL unless a .kmdb target exists */
    size_t     n_qk;
    int        pair_dedup;   /* --pair-dedup */
    int        have_fastq;   /* any FASTQ target (for the .kmdb mismatch warning) */
    int        read_hash;    /* --read-hash: min hashes at counts[n_total + mg_idx] */
    int        n_total;      /* columns (targets + background) */
    uint64_t  *n_frags;      /* per column: fragments numbered (read hash) */
} ksd_job;

typedef struct {
    ksd_job       **queue;
    int            head, tail, size, capacity;
    int            active, shutdown;
    pthread_mutex_t mutex;
    pthread_cond_t  has_work;
    pthread_cond_t  all_done;
} ksd_pool;

typedef struct {
    BIO_hash h;
    int      n_mg;
    int      seed;
    char   **list;
    int      n, cap;
    long     n_dup, n_bad, n_len;
} KmerLoad;

/* ── forward declarations ──────────────────────────────────────────── */

static int       get_file_type(const char *s);
static char     *strip_basename(const char *path);
static char     *pair_sample_name(const char *f1, const char *f2);
static int       parse_targets(const char *file, Metagenome **out);
static int       parse_kmer_tsv(const char *file, BIO_hash h, int n_mg,
                                 int *seed_out, char ***kmer_list_out);
static void      scan_target(ksd_job *job);
static void     *ksd_worker(void *arg);
static ksd_pool *ksd_pool_create(int nthreads, int capacity);
static void      ksd_pool_submit(ksd_pool *pool, ksd_job *job);
static void      ksd_pool_wait_and_destroy(ksd_pool *pool);
static void      usage(void);

/* ── file type helpers ─────────────────────────────────────────────── */

static int get_file_type(const char *s) {
    if (strcmp(s,"SE")==0 || strcmp(s,"se")==0) return NOT_PAIRED_END;
    if (strcmp(s,"PE")==0 || strcmp(s,"pe")==0) return IS_PAIRED_END;
    if (strcmp(s,"PEI")==0||strcmp(s,"pei")==0||strcmp(s,"IPE")==0||strcmp(s,"ipe")==0)
        return IS_PAIRED_END_INTERLEAVE;
    if (strcasecmp(s,"KMDB")==0) return IS_KMDB;
    return UNKNOWN_FILE_TYPE;
}

static int ends_with_kmdb(const char *p) {
    size_t n = p ? strlen(p) : 0;
    return n > 5 && strcmp(p + n - 5, ".kmdb") == 0;
}

static char *strip_basename(const char *path) {
    const char *base = strrchr(path, '/');
    base = base ? base + 1 : path;
    char *b = strdup(base);
    const char *exts[] = {
        ".fastq.gz", ".fasta.gz", ".fa.gz", ".fq.gz",
        ".fastq",    ".fasta",    ".fa",    ".fq",
        ".kmdb",
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

/* ── sample naming for PE entries ──────────────────────────────────────
 * Take the extension-stripped basenames of both mates, keep their longest
 * common prefix (i.e. cut where the two names start to disagree), then
 * clean up what is left:
 *   1. drop trailing connectors  (_ - . : space)
 *   2. drop a trailing read tag  (R, PE, read) if it sits behind a
 *      connector, so "x_R1_001"/"x_R2_001" -> "x_R" -> "x"
 *   3. drop connectors again
 * Falls back to the stripped basename of f1 if nothing sensible is left.
 *   1001099B_150804_B6_s09_PE1 / ..._PE2  -> 1001099B_150804_B6_s09
 *   DK-D32-3524-250106-H3_R1_001 / _R2_001 -> DK-D32-3524-250106-H3
 */
static int is_connector(char c) {
    return c == '_' || c == '-' || c == '.' || c == ':' || c == ' ';
}

static void trim_connectors(char *s) {
    size_t n = strlen(s);
    while (n > 0 && is_connector(s[n-1])) s[--n] = '\0';
}

static char *pair_sample_name(const char *f1, const char *f2) {
    char *a = strip_basename(f1);
    char *b = strip_basename(f2);

    size_t i = 0;
    while (a[i] && b[i] && a[i] == b[i]) i++;
    a[i] = '\0';                       /* a is now the common prefix */
    free(b);

    trim_connectors(a);

    /* strip a dangling read tag such as "_R", "_PE", "_read" */
    static const char *tags[] = { "read", "Read", "READ", "PE", "pe", "R", "r", NULL };
    size_t n = strlen(a);
    for (int t = 0; tags[t]; t++) {
        size_t tl = strlen(tags[t]);
        if (n > tl + 1 && strcmp(a + n - tl, tags[t]) == 0 && is_connector(a[n - tl - 1])) {
            a[n - tl] = '\0';
            trim_connectors(a);
            break;
        }
    }

    if (a[0] == '\0') {                /* nothing in common: keep old behaviour */
        free(a);
        return strip_basename(f1);
    }
    return a;
}

/* ── targets parser ──────────────────────────────────────────────────
 * Two layouts are accepted:
 *
 *  1. Headerless (original):   type <TAB> file1 [<TAB> file2]
 *     Sample name is derived from the file name(s) (see pair_sample_name).
 *
 *  2. With a header line (first non-blank line, optional leading '#').
 *     Columns are found by name, in any order, case-insensitive:
 *        type | sequencing_type      PE / SE / PEI
 *        file1, file2                 or  seq_file (comma-separated mates)
 *        sample_name                  optional; used verbatim as the column
 *                                     name in the output. Empty cell -> the
 *                                     automatic name is used for that row.
 *     Other columns (sample_id, isolates_to_track, ...) are ignored.
 *
 *  .kmdb targets: a file1/seq_file ending in .kmdb is a sample count database
 *  whatever the type column says (type may also be KMDB or empty). In the
 *  headerless layout a line holding only a .kmdb path is accepted too.
 *  Column name = sample_name, else the file name without .kmdb.
 */

#define MAX_FIELDS 64

static int split_tabs(char *line, char **f, int max) {
    int n = 0;
    char *p = line;
    while (n < max) {
        f[n++] = p;
        char *t = strchr(p, '\t');
        if (!t) break;
        *t = '\0';
        p = t + 1;
    }
    return n;
}

static char *trim_ws(char *s) {
    while (*s == ' ' || *s == '\t') s++;
    size_t n = strlen(s);
    while (n > 0 && (s[n-1] == ' ' || s[n-1] == '\t')) s[--n] = '\0';
    return s;
}

/* NULL if the column is absent or the cell is empty */
static char *cell(char **f, int nf, int col) {
    if (col < 0 || col >= nf) return NULL;
    char *s = trim_ws(f[col]);
    return *s ? s : NULL;
}

static int parse_targets(const char *file, Metagenome **out) {
    FILE *fp = fopen(file, "r");
    if (!fp) {
        fprintf(stderr, "cannot open targets file %s: %s\n", file, strerror(errno));
        exit(1);
    }

    int cap = 64, n = 0;
    Metagenome *mgs = malloc(cap * sizeof(Metagenome));

    /* column indices; defaults = headerless layout */
    int c_type = 0, c_f1 = 1, c_f2 = 2, c_seq = -1, c_name = -1;
    int first_line = 1;

    char *line = NULL;
    size_t len = 0;
    while (getline(&line, &len, fp) != -1) {
        char *nl = strchr(line, '\n'); if (nl) *nl = '\0';
        char *cr = strchr(line, '\r'); if (cr) *cr = '\0';
        if (line[0] == '\0') continue;

        char *f[MAX_FIELDS];

        /* header detection: first non-blank line only */
        if (first_line) {
            first_line = 0;
            char *hl = strdup(line);
            char *h = hl;
            while (*h == '#') h++;
            int nf = split_tabs(h, f, MAX_FIELDS);
            int ht = -1, h1 = -1, h2 = -1, hs = -1, hn = -1;
            for (int i = 0; i < nf; i++) {
                char *c = trim_ws(f[i]);
                if      (!strcasecmp(c, "type") || !strcasecmp(c, "sequencing_type")) ht = i;
                else if (!strcasecmp(c, "file1"))       h1 = i;
                else if (!strcasecmp(c, "file2"))       h2 = i;
                else if (!strcasecmp(c, "seq_file"))    hs = i;
                else if (!strcasecmp(c, "sample_name")) hn = i;
            }
            free(hl);
            if (ht >= 0 || h1 >= 0 || hs >= 0 || hn >= 0) {
                if (ht < 0 || (h1 < 0 && hs < 0)) {
                    fprintf(stderr, "targets header needs a type/sequencing_type column and "
                                    "file1 or seq_file column: %s\n", file);
                    exit(1);
                }
                c_type = ht; c_f1 = h1; c_f2 = h2; c_seq = hs; c_name = hn;
                continue;   /* header consumed */
            }
        }

        if (line[0] == '#') continue;

        int nf = split_tabs(line, f, MAX_FIELDS);
        char *type_s = cell(f, nf, c_type);
        char *f1 = NULL, *f2 = NULL;
        if (c_seq >= 0) {
            f1 = cell(f, nf, c_seq);
            if (f1) {
                char *comma = strchr(f1, ',');
                if (comma) { *comma = '\0'; f2 = trim_ws(comma + 1); if (!*f2) f2 = NULL; }
                f1 = trim_ws(f1);
            }
        } else {
            f1 = cell(f, nf, c_f1);
            f2 = cell(f, nf, c_f2);
        }
        /* headerless line that is just a .kmdb path */
        if (!f1 && c_seq < 0 && c_type == 0 && ends_with_kmdb(type_s)) {
            f1 = type_s;
            type_s = NULL;
        }
        if (!f1) continue;

        int type;
        if (ends_with_kmdb(f1)) {
            type = IS_KMDB;
            f2 = NULL;
        } else {
            if (!type_s) continue;
            type = get_file_type(type_s);
            if (type == IS_KMDB) {
                fprintf(stderr, "type KMDB but %s does not end in .kmdb, skipping\n", f1);
                continue;
            }
        }
        if (type == UNKNOWN_FILE_TYPE) {
            fprintf(stderr, "unknown type '%s', skipping line\n", type_s);
            continue;
        }

        if (type == IS_PAIRED_END && !f2) {
            fprintf(stderr, "PE entry missing second file for %s, skipping\n", f1);
            continue;
        }
        if (type != IS_PAIRED_END) f2 = NULL;

        if (n == cap) { cap *= 2; mgs = realloc(mgs, cap * sizeof(Metagenome)); }
        mgs[n].type  = type;
        mgs[n].file1 = strdup(f1);
        mgs[n].file2 = f2 ? strdup(f2) : NULL;

        char *given = cell(f, nf, c_name);
        if (given)      mgs[n].basename = strdup(given);                 /* explicit sample_name */
        else if (f2)    mgs[n].basename = pair_sample_name(f1, f2);      /* PE: common prefix    */
        else            mgs[n].basename = strip_basename(f1);            /* SE/PEI: strip ext    */
        n++;
    }
    free(line);
    fclose(fp);
    *out = mgs;
    return n;
}

/* ── kmer file parser ──────────────────────────────────────────────── */

static int is_acgt(const char *s) {
    if (!*s) return 0;
    for (; *s; s++) {
        switch (*s) {
            case 'A': case 'C': case 'G': case 'T':
            case 'a': case 'c': case 'g': case 't':
                break;
            default:
                return 0;
        }
    }
    return 1;
}

/* column index of '#kmer' / 'kmer' in a tab-separated header, else -1 */
static int find_header_col(const char *line) {
    char hdr[MAX_LINE];
    strncpy(hdr, line, MAX_LINE - 1);
    hdr[MAX_LINE - 1] = '\0';
    char *save = NULL;
    int col = 0;
    for (char *tok = strtok_r(hdr, "\t", &save); tok;
         tok = strtok_r(NULL, "\t", &save), col++) {
        if (strcmp(tok, "#kmer") == 0 || strcmp(tok, "kmer") == 0) return col;
    }
    return -1;
}

static void add_kmer_line(KmerLoad *L, char *line, int kmer_col, const char *delim) {
    char *save = NULL;
    char *tok = strtok_r(line, delim, &save);
    for (int c = 0; c < kmer_col && tok; c++)
        tok = strtok_r(NULL, delim, &save);
    if (!tok || !*tok) { L->n_bad++; return; }

    size_t klen = strlen(tok);
    if (klen >= MAX_KMER_LEN || !is_acgt(tok)) { L->n_bad++; return; }
    for (char *p = tok; *p; p++) *p = (char)toupper((unsigned char)*p);

    if (L->seed < 0) L->seed = (int)klen;
    else if ((int)klen != L->seed) { L->n_len++; return; }

    /* canonical (oriented) form */
    char rc[MAX_KMER_LEN];
    memcpy(rc, tok, klen + 1);
    BIO_reverseComplement(rc);
    char *oriented = orient_string(tok, rc, L->seed);

    /* skip duplicates (same canonical kmer) */
    if (BIO_searchHash(L->h, oriented) != NULL) { L->n_dup++; return; }

    unsigned int *counts = calloc(L->n_mg, sizeof(unsigned int));
    BIO_addHashData(L->h, oriented, counts);

    if (L->n == L->cap) {
        L->cap *= 2;
        L->list = realloc(L->list, L->cap * sizeof(char *));
    }
    L->list[L->n++] = strdup(tok);   /* original orientation for output */
}

static int parse_kmer_tsv(const char *file, BIO_hash h, int n_mg,
                           int *seed_out, char ***kmer_list_out) {
    gzFile fp = gzopen(file, "r");   /* zlib reads plain files transparently */
    if (!fp) {
        fprintf(stderr, "cannot open kmer file %s: %s\n", file, strerror(errno));
        exit(1);
    }

    KmerLoad L = { h, n_mg, -1, malloc(4096 * sizeof(char *)), 0, 4096, 0, 0, 0 };
    char line[MAX_LINE];
    int kmer_col = -1, pending = 0;
    const char *delim = "\t";

    /* detect format: header with '#kmer', or headerless list */
    while (gzgets(fp, line, MAX_LINE)) {
        line[strcspn(line, "\r\n")] = '\0';
        if (line[0] == '\0') continue;

        int c = find_header_col(line);
        if (c >= 0) { kmer_col = c; break; }
        if (line[0] == '#') continue;          /* stats / comment line */

        char probe[MAX_LINE];
        strncpy(probe, line, MAX_LINE - 1);
        probe[MAX_LINE - 1] = '\0';
        char *save = NULL;
        char *t = strtok_r(probe, " \t", &save);
        if (t && is_acgt(t)) {
            kmer_col = 0;
            delim    = " \t";
            pending  = 1;                       /* this line is already data */
            break;
        }
        fprintf(stderr, "%s: no '#kmer' header and first data line does not "
                "start with a kmer:\n  %.80s\n", file, line);
        exit(1);
    }

    if (kmer_col < 0) {
        fprintf(stderr, "WARNING: %s contains no kmers (empty or comments only)\n", file);
        gzclose(fp);
        *seed_out = -1;
        *kmer_list_out = L.list;
        return 0;
    }
    fprintf(stderr, "kmer file %s: %s, kmer column %d\n", file,
            pending ? "headerless list" : "header", kmer_col);

    if (pending) add_kmer_line(&L, line, kmer_col, delim);
    while (gzgets(fp, line, MAX_LINE)) {
        line[strcspn(line, "\r\n")] = '\0';
        if (line[0] == '\0' || line[0] == '#') continue;
        add_kmer_line(&L, line, kmer_col, delim);
    }
    gzclose(fp);

    if (L.n_dup || L.n_bad || L.n_len)
        fprintf(stderr, "kmer file: skipped %ld duplicate, %ld invalid, "
                "%ld wrong-length kmer line(s)\n", L.n_dup, L.n_bad, L.n_len);
    if (L.n == 0)
        fprintf(stderr, "WARNING: %s contains no valid kmers\n", file);

    *seed_out = L.seed;
    *kmer_list_out = L.list;
    return L.n;
}

/* ── read scanner ──────────────────────────────────────────────────── */

/* A-0 C-1 G-2 T-3 (either case), anything else 4: same as KMDB_NT2BIT,
   duplicated so builds without zstd/kmdb.c have it too. */
static uint8_t NT2BIT[256];
static void nt2bit_init(void) {
    memset(NT2BIT, 4, sizeof NT2BIT);
    NT2BIT['A'] = NT2BIT['a'] = 0; NT2BIT['C'] = NT2BIT['c'] = 1;
    NT2BIT['G'] = NT2BIT['g'] = 2; NT2BIT['T'] = NT2BIT['t'] = 3;
}

static inline uint64_t mix64(uint64_t x) {
    x ^= x >> 33; x *= 0xff51afd7ed558ccdull;
    x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ull;
    x ^= x >> 33;
    return x;
}

/* Per-thread set of the canonical k-mers counted in the current read pair.
   Generation-stamped: a new pair is ++cur, no clearing. Keys are the exact
   2-bit canonical k-mer for k <= 31; for longer k a 64-bit hash of the
   canonical string (a collision would merge two k-mers of one pair, at
   ~n^2/2^65 per pair, i.e. never in practice). */
typedef struct {
    uint64_t *keys;
    uint32_t *gen;
    size_t    cap;
    uint32_t  cur;
} fset_t;

static void fset_begin(fset_t *s, size_t n_windows) {
    size_t need = 1024;
    while (need < n_windows * 2) need <<= 1;
    if (need > s->cap) {
        free(s->keys); free(s->gen);
        s->cap  = need;
        s->keys = malloc(need * sizeof(uint64_t));
        s->gen  = calloc(need, sizeof(uint32_t));
        if (!s->keys || !s->gen) { fprintf(stderr, "out of memory\n"); exit(1); }
        s->cur  = 0;
    }
    if (++s->cur == 0) {
        memset(s->gen, 0, s->cap * sizeof(uint32_t));
        s->cur = 1;
    }
}

/* 1 = new in this pair, 0 = already counted */
static inline int fset_add(fset_t *s, uint64_t key) {
    size_t m = s->cap - 1, i = (size_t)mix64(key) & m;
    for (;;) {
        if (s->gen[i] != s->cur) { s->gen[i] = s->cur; s->keys[i] = key; return 1; }
        if (s->keys[i] == key) return 0;
        i = (i + 1) & m;
    }
}

typedef struct {
    BIO_hash  h;
    int       k, mg_idx;
    uint64_t  total;            /* increments made = total_evaluated */
    uint64_t  hits;             /* query-k-mer increments made */
    uint64_t  skip_total;       /* windows skipped by pair dedup */
    uint64_t  skip_hits;        /* query-k-mer increments skipped */
    uint64_t  n_pairs;
    char     *rc;  size_t rc_cap;
    fset_t    fs;
    int       rh_off;           /* --read-hash: offset of the hash slots (n_total), else -1 */
    uint64_t  frag;             /* fragments numbered so far (read hash) */
    uint32_t  cur_h;            /* hash of the current fragment */
} scan_t;

/* start a fragment: number it if it holds a read of length >= k (as
   kmer_metagenome_db does), and set its hash */
static inline void frag_begin(scan_t *S, size_t len1, size_t len2) {
    if (S->rh_off < 0) return;
    if (len1 >= (size_t)S->k || len2 >= (size_t)S->k) S->cur_h = read_hash(S->frag++);
}

/* Count the k-mers of one read. fs != NULL: pair dedup (skip k-mers already
   counted in this pair). Only ACGT-only windows count, as in
   kmer_metagenome_db. */
static void scan_read(scan_t *S, char *s, size_t len, fset_t *fs) {
    const int k = S->k;
    if (len < (size_t)k) return;
    BIO_stringToUpper(s);
    if (len + 1 > S->rc_cap) {
        S->rc_cap = len + 1;
        S->rc = realloc(S->rc, S->rc_cap);
    }
    memcpy(S->rc, s, len + 1);
    BIO_reverseComplement(S->rc);

    const int packed = k <= 31;
    const int shift = 2 * (k - 1);
    const uint64_t mask = packed ? ((1ull << (2 * k)) - 1) : 0;
    uint64_t fwd = 0, rev = 0;
    int run = 0;
    for (size_t i = 0; i < len; i++) {
        uint8_t x = NT2BIT[(unsigned char)s[i]];
        if (x > 3) { run = 0; fwd = rev = 0; continue; }
        if (packed) {
            fwd = ((fwd << 2) | x) & mask;
            rev = (rev >> 2) | ((uint64_t)(3 - x) << shift);
        }
        if (++run < k) continue;

        size_t j = i + 1 - (size_t)k;              /* window s[j, j+k) */
        char *w  = s + j;
        char *wr = S->rc + (len - j - (size_t)k);  /* its reverse complement */
        char tw = w[k], tr = wr[k];
        w[k] = '\0'; wr[k] = '\0';
        char *orient = strcmp(w, wr) > 0 ? w : wr; /* max == max(fwd,rev) */
        unsigned int *count = (unsigned int *)BIO_searchHash(S->h, orient);

        int counted = 1;
        if (fs) {
            uint64_t key;
            if (packed) key = fwd > rev ? fwd : rev;
            else {                                 /* FNV-1a 64 */
                key = 1469598103934665603ull;
                for (const char *p = orient; *p; p++)
                    key = (key ^ (unsigned char)*p) * 1099511628211ull;
            }
            counted = fset_add(fs, key);
        }
        w[k] = tw; wr[k] = tr;

        if (counted) {
            S->total++;
            if (count) {                                   /* only this thread writes this column */
                count[S->mg_idx]++; S->hits++;
                if (S->rh_off >= 0) {
                    unsigned int *mh = &count[S->rh_off + S->mg_idx];
                    *mh = read_hash_min(*mh, S->cur_h);
                }
            }
        } else {
            S->skip_total++;
            if (count) S->skip_hits++;
        }
    }
}

static size_t n_windows(size_t len, int k) {
    return len >= (size_t)k ? len - (size_t)k + 1 : 0;
}

/* mate names, ignoring a trailing /1 /2 .1 .2 _1 _2 */
static size_t mate_stem(const kstring_t *n) {
    size_t l = n->l;
    if (l >= 2 && (n->s[l-1] == '1' || n->s[l-1] == '2') &&
        (n->s[l-2] == '/' || n->s[l-2] == '.' || n->s[l-2] == '_'))
        l -= 2;
    return l;
}

static void check_mates(const kstring_t *a, const kstring_t *b, const char *what,
                        uint64_t pair_no) {
    size_t la = mate_stem(a), lb = mate_stem(b);
    if (la == lb && memcmp(a->s, b->s, la) == 0) return;
    fprintf(stderr, "error: %s: mates out of sync at pair %" PRIu64 " (%.*s vs %.*s); "
            "use --no-pair-dedup to count them as single-end reads\n",
            what, pair_no, (int)a->l, a->s, (int)b->l, b->s);
    exit(1);
}

static void scan_pair(scan_t *S, kstring_t *r1, kstring_t *r2) {
    frag_begin(S, r1->l, r2->l);
    fset_begin(&S->fs, n_windows(r1->l, S->k) + n_windows(r2->l, S->k));
    scan_read(S, r1->s, r1->l, &S->fs);
    scan_read(S, r2->s, r2->l, &S->fs);
    S->n_pairs++;
}

static gzFile open_or_die(const char *path) {
    gzFile fp = gzopen(path, "r");
    if (!fp) {
        fprintf(stderr, "error: cannot open %s: %s\n", path, strerror(errno));
        exit(1);
    }
    gzbuffer(fp, 1 << 20);
    return fp;
}

static void die_malformed(const char *path, int l) {
    fprintf(stderr, "error: %s is truncated or malformed (kseq %d)\n", path, l);
    exit(1);
}

/* single-end (also every file with --no-pair-dedup): no dedup */
static void scan_se(scan_t *S, const char *path) {
    gzFile fp = open_or_die(path);
    kseq_t *seq = kseq_init(fp);
    int l;
    while ((l = kseq_read(seq)) >= 0) {
        frag_begin(S, seq->seq.l, 0);
        scan_read(S, seq->seq.s, seq->seq.l, NULL);
    }
    if (l < -1) die_malformed(path, l);
    kseq_destroy(seq);
    gzclose(fp);
}

static void scan_pe(scan_t *S, const char *f1, const char *f2) {
    gzFile a = open_or_die(f1), b = open_or_die(f2);
    kseq_t *s1 = kseq_init(a), *s2 = kseq_init(b);
    for (;;) {
        int l1 = kseq_read(s1), l2 = kseq_read(s2);
        if (l1 < -1) die_malformed(f1, l1);
        if (l2 < -1) die_malformed(f2, l2);
        if (l1 < 0 && l2 < 0) break;
        if (l1 < 0 || l2 < 0) {
            fprintf(stderr, "error: %s and %s hold different numbers of reads "
                    "(%s ends after %" PRIu64 " pairs)\n", f1, f2,
                    l1 < 0 ? f1 : f2, S->n_pairs);
            exit(1);
        }
        check_mates(&s1->name, &s2->name, f1, S->n_pairs + 1);
        scan_pair(S, &s1->seq, &s2->seq);
    }
    kseq_destroy(s1); kseq_destroy(s2);
    gzclose(a); gzclose(b);
}

static void kstr_copy(kstring_t *d, const kstring_t *s) {
    if (d->m < s->l + 1) { d->m = s->l + 1; d->s = realloc(d->s, d->m); }
    memcpy(d->s, s->s, s->l);
    d->s[s->l] = '\0';
    d->l = s->l;
}

static void scan_pei(scan_t *S, const char *path) {
    gzFile fp = open_or_die(path);
    kseq_t *seq = kseq_init(fp);
    kstring_t name1 = {0, 0, NULL}, seq1 = {0, 0, NULL};
    for (;;) {
        int l = kseq_read(seq);
        if (l < -1) die_malformed(path, l);
        if (l < 0) break;
        kstr_copy(&name1, &seq->name);
        kstr_copy(&seq1, &seq->seq);
        l = kseq_read(seq);
        if (l < -1) die_malformed(path, l);
        if (l < 0) {
            fprintf(stderr, "error: %s (PEI) has an odd number of records; the "
                    "last read has no mate\n", path);
            exit(1);
        }
        check_mates(&name1, &seq->name, path, S->n_pairs + 1);
        scan_pair(S, &seq1, &seq->seq);
    }
    free(name1.s); free(seq1.s);
    kseq_destroy(seq);
    gzclose(fp);
}

static void scan_target(ksd_job *job) {
    const Metagenome *mg = job->mg;
    scan_t S;
    memset(&S, 0, sizeof S);
    S.h = job->h; S.k = job->seed; S.mg_idx = job->mg_idx;
    S.rh_off = job->read_hash ? job->n_total : -1;

    const char *how;
    if (!job->pair_dedup) {
        scan_se(&S, mg->file1);                         /* the old behaviour */
        if (mg->type == IS_PAIRED_END && mg->file2) scan_se(&S, mg->file2);
        how = "no pair dedup";
    } else if (mg->type == IS_PAIRED_END) {
        scan_pe(&S, mg->file1, mg->file2);
        how = "PE, pair dedup";
    } else if (mg->type == IS_PAIRED_END_INTERLEAVE) {
        scan_pei(&S, mg->file1);
        how = "PEI, pair dedup";
    } else {
        scan_se(&S, mg->file1);
        how = "SE";
    }
    job->total_kmer_counts[job->mg_idx] += S.total;
    if (job->read_hash) job->n_frags[job->mg_idx] = S.frag;

    uint64_t all = S.total + S.skip_total, qall = S.hits + S.skip_hits;
    if (job->pair_dedup && mg->type != NOT_PAIRED_END)
        fprintf(stderr, "  %s (%s): %" PRIu64 " pairs; skipped %" PRIu64 " of %" PRIu64
                " k-mer windows (%.2f%%) and %" PRIu64 " of %" PRIu64 " query k-mer hits "
                "(%.2f%%); total_evaluated %" PRIu64 "\n", mg->basename, how, S.n_pairs,
                S.skip_total, all, all ? 100.0 * S.skip_total / all : 0.0,
                S.skip_hits, qall, qall ? 100.0 * S.skip_hits / qall : 0.0, S.total);
    else
        fprintf(stderr, "  %s (%s): %" PRIu64 " query k-mer hits; total_evaluated %"
                PRIu64 "\n", mg->basename, how, S.hits, S.total);
    free(S.rc); free(S.fs.keys); free(S.fs.gen);
}

/* ── .kmdb target lookup ───────────────────────────────────────────────
 * qk[] holds the query k-mers as 2-bit canonical keys, sorted. Both qk[] and
 * the database are sorted, so a single forward pass over the blocks works:
 * for each block, find the first query >= block.first_kmer; if it is past
 * block.last_kmer the block is skipped without reading it, otherwise the
 * block is decoded and the queries inside its range are matched against it.
 * Any error aborts the program: silently missing counts would look like
 * absent strains.
 */

#ifdef HAVE_ZSTD
static int cmp_qkmer(const void *a, const void *b) {
    uint64_t x = ((const QKmer *)a)->key, y = ((const QKmer *)b)->key;
    return (x > y) - (x < y);
}

/* first index in [lo, n) with qk[i].key >= key */
static size_t qk_lower_bound(const QKmer *qk, size_t lo, size_t n, uint64_t key) {
    size_t hi = n;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (qk[mid].key < key) lo = mid + 1; else hi = mid;
    }
    return lo;
}

static void scan_kmdb(const char *path, const QKmer *qk, size_t n_qk, int seed,
                      int mg_idx, uint64_t *total_kmer_counts, int pair_dedup,
                      int have_fastq, int rh_off, uint64_t *n_frags) {
    kmdb_t *db = kmdb_open(path);
    if (!db) exit(1);

    if ((int)db->k != seed) {
        fprintf(stderr, "error: %s has k=%u but the kmer file has k=%d\n",
                path, db->k, seed);
        exit(1);
    }
    char *po = kmdb_meta_get(db, "presence_only");
    if (po && strcmp(po, "1") == 0) {
        fprintf(stderr, "error: %s is a presence-only (scrub) database; build "
                "target databases with kmer_metagenome_db -m sample\n", path);
        exit(1);
    }
    free(po);
    char *mc = kmdb_meta_get(db, "min_count");
    if (mc && strtoull(mc, NULL, 10) > 1)
        fprintf(stderr, "WARNING: %s was built with min_count=%s; k-mers seen "
                "fewer times are reported as 0\n", path, mc);
    free(mc);
    /* total_evaluated = sum of all counts of the sample. Databases written
       before that key existed were never pair-deduplicated, so their sum
       of counts is valid_kmer_positions. */
    char *te = kmdb_meta_get(db, "total_evaluated");
    if (!te) te = kmdb_meta_get(db, "valid_kmer_positions");
    if (!te) {
        fprintf(stderr, "error: %s has neither total_evaluated nor "
                "valid_kmer_positions metadata\n", path);
        exit(1);
    }
    total_kmer_counts[mg_idx] += strtoull(te, NULL, 10);
    free(te);
    char *pd = kmdb_meta_get(db, "pair_dedup");
    int db_dedup = pd && strcmp(pd, "1") == 0;
    free(pd);
    if (have_fastq && db_dedup != pair_dedup)
        fprintf(stderr, "WARNING: %s was built %s pair dedup but FASTQ targets in "
                "this run are counted %s it; their counts are not comparable\n",
                path, db_dedup ? "with" : "without", pair_dedup ? "with" : "without");

    /* --read-hash: use the database's hashes if it has them */
    const int use_rh = rh_off >= 0 && db->has_read_hash;
    if (rh_off >= 0) {
        if (!db->has_read_hash)
            fprintf(stderr, "WARNING: %s was built without --read-hash; its read "
                    "hashes are reported as 0 (rebuild with kmer_metagenome_db "
                    "-m sample --read-hash)\n", path);
        char *nf = kmdb_meta_get(db, "n_fragments");
        n_frags[mg_idx] = nf ? strtoull(nf, NULL, 10) : 0;
        free(nf);
    }

    kmdb_block_buf buf;
    kmdb_block_buf_init(&buf);
    size_t j = 0;
    uint64_t decoded = 0, hits = 0;
    for (uint64_t b = 0; b < db->n_blocks && j < n_qk; b++) {
        const kmdb_block_info *bi = &db->blocks[b];
        j = qk_lower_bound(qk, j, n_qk, bi->first_kmer);
        if (j == n_qk) break;
        if (qk[j].key > bi->last_kmer) continue;          /* nothing here */

        size_t jend = qk_lower_bound(qk, j, n_qk, bi->last_kmer);
        if (jend < n_qk && qk[jend].key == bi->last_kmer) jend++;

        long n = kmdb_read_block(db, b, &buf);
        if (n < 0) {
            fprintf(stderr, "error: cannot read %s\n", path);
            exit(1);
        }
        decoded++;
        const uint64_t *keys = buf.kmers;
        size_t m = jend - j;
        if (m * 24 < (size_t)n) {
            /* few queries: binary-search each one in the block */
            for (; j < jend; j++) {
                size_t lo = 0, hi = (size_t)n;
                while (lo < hi) {
                    size_t mid = lo + (hi - lo) / 2;
                    if (keys[mid] < qk[j].key) lo = mid + 1; else hi = mid;
                }
                if (lo < (size_t)n && keys[lo] == qk[j].key) {
                    qk[j].counts[mg_idx] += buf.counts[lo];
                    if (use_rh) qk[j].counts[rh_off + mg_idx] = buf.min_hash[lo];
                    hits++;
                }
            }
        } else {
            /* many queries: linear merge */
            size_t i = 0;
            while (j < jend && i < (size_t)n) {
                if      (keys[i] < qk[j].key) i++;
                else if (keys[i] > qk[j].key) j++;
                else {
                    qk[j].counts[mg_idx] += buf.counts[i];
                    if (use_rh) qk[j].counts[rh_off + mg_idx] = buf.min_hash[i];
                    hits++; i++; j++;
                }
            }
            j = jend;
        }
    }
    fprintf(stderr, "  %s: decoded %" PRIu64 " / %" PRIu64 " blocks, %" PRIu64
            " query k-mers present\n", path, decoded, db->n_blocks, hits);
    kmdb_block_buf_free(&buf);
    kmdb_close(db);
}
#endif

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

        if (job->mg->type == IS_KMDB) {
#ifdef HAVE_ZSTD
            scan_kmdb(job->mg->file1, job->qk, job->n_qk, job->seed,
                      job->mg_idx, job->total_kmer_counts, job->pair_dedup,
                      job->have_fastq, job->read_hash ? job->n_total : -1,
                      job->n_frags);
#endif
        } else {
            scan_target(job);
        }
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
    pthread_t *threads = malloc(sizeof(pthread_t) * nthreads);
    for (int i = 0; i < nthreads; i++)
        pthread_create(&threads[i], NULL, ksd_worker, pool);
    free(threads);
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
    struct timespec ts = {0, 10000000};
    nanosleep(&ts, NULL);
    pthread_mutex_destroy(&pool->mutex);
    pthread_cond_destroy(&pool->has_work);
    pthread_cond_destroy(&pool->all_done);
    free(pool->queue);
    free(pool);
}

/* ── usage ─────────────────────────────────────────────────────────── */

static void usage(void) {
    fprintf(stderr, "Usage: kmer_strain_detect -k <kmers> -B <targets.txt> [-O <dir>] [-n <name> | -o <out.kmer_hits.tsv.gz>] [-G <background.txt>] [-t threads]\n\n");
    fprintf(stderr, "  -k  kmers (plain or gzipped): TSV with a '#kmer' header column, or a\n");
    fprintf(stderr, "      headerless list with the kmer in the first column ('#' lines skipped)\n");
    fprintf(stderr, "  -B  target metagenomes file (tab-separated: type file1 [file2]), or with a header\n");
    fprintf(stderr, "      line naming columns: type|sequencing_type, file1/file2 or seq_file (comma-sep),\n");
    fprintf(stderr, "      and optional sample_name (used as the output column name)\n");
    fprintf(stderr, "  -G  background metagenomes file (same format as -B); columns prefixed 'b_'\n");
    fprintf(stderr, "      types: PE, SE, PEI, KMDB\n");
    fprintf(stderr, "      A file ending in .kmdb (kmer_metagenome_db -m sample) is looked up\n");
    fprintf(stderr, "      instead of scanned; a targets file may simply list .kmdb paths.\n");
    fprintf(stderr, "  -O  output directory (default: ., created if missing)\n");
    fprintf(stderr, "  -n  output basename (default: -k filename without extensions);\n");
    fprintf(stderr, "      writes <-O>/<-n>.kmer_hits.tsv.gz\n");
    fprintf(stderr, "  -o  output file; overrides -O/-n (e.g. sample.kmer_hits.tsv.gz)\n");
    fprintf(stderr, "  -t  threads (default: 4)   [was -j]\n");
    fprintf(stderr, "  --pair-dedup (default) / --no-pair-dedup\n");
    fprintf(stderr, "      PE / PEI FASTQ targets: count each k-mer at most once per read pair\n");
    fprintf(stderr, "      (in R1, R2 or both, repeats included). SE targets are unchanged.\n");
    fprintf(stderr, "      --no-pair-dedup counts every mate on its own (the old counts).\n");
    fprintf(stderr, "      .kmdb targets carry the setting they were built with.\n");
    fprintf(stderr, "  --read-hash\n");
    fprintf(stderr, "      also write <name>.kmer_readhash.tsv.gz: per k-mer and sample the\n");
    fprintf(stderr, "      smallest hash of the reads (pairs) containing it, 0 = not seen;\n");
    fprintf(stderr, "      last row n_fragments. Equal values = same read (exact at count 1).\n");
    fprintf(stderr, "      .kmdb targets need kmer_metagenome_db -m sample --read-hash.\n\n");
    fprintf(stderr, "Output: gzipped TSV — first column is #kmer, then target columns,\n");
    fprintf(stderr, "        then b_* columns for background metagenomes. Last row\n");
    fprintf(stderr, "        total_evaluated = sum of the counts of all k-mers of the sample\n");
    fprintf(stderr, "        (ACGT-only windows, after pair dedup): the relative-abundance\n");
    fprintf(stderr, "        denominator.\n");
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
    int   pair_dedup      = 1;
    int   read_hash_opt   = 0;
    int   c;

    enum { OPT_PAIR_DEDUP = 1000, OPT_NO_PAIR_DEDUP, OPT_READ_HASH };
    static const struct option longopts[] = {
        { "read-hash",     no_argument, NULL, OPT_READ_HASH },
        { "pair-dedup",    no_argument, NULL, OPT_PAIR_DEDUP },
        { "no-pair-dedup", no_argument, NULL, OPT_NO_PAIR_DEDUP },
        { "help",          no_argument, NULL, 'h' },
        { NULL, 0, NULL, 0 }
    };
    nt2bit_init();
    while ((c = getopt_long(argc, argv, "k:B:G:O:n:o:t:j:h", longopts, NULL)) != -1) {
        switch (c) {
            case OPT_PAIR_DEDUP:    pair_dedup = 1; break;
            case OPT_NO_PAIR_DEDUP: pair_dedup = 0; break;
            case OPT_READ_HASH:     read_hash_opt = 1; break;
            case 'k': kmer_file       = strdup(optarg); break;
            case 'B': targets_file    = strdup(optarg); break;
            case 'G': background_file = strdup(optarg); break;
            case 'O': out_dir         = strdup(optarg); break;
            case 'n': out_name        = strdup(optarg); break;
            case 'o': outfile         = strdup(optarg); break;
            case 't': num_threads     = atoi(optarg);   break;
            case 'j': cli_renamed_flag("kmer_strain_detect", 'j', "-t"); break;
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
        char *resolved = cli_resolve_output("kmer_strain_detect", out_dir, out_name,
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

    /* 1b. parse background metagenomes and append with b_ prefix */
    int n_total = n_mg;
    if (background_file) {
        Metagenome *bg;
        int n_bg = parse_targets(background_file, &bg);
        fprintf(stderr, "loaded %d background metagenome(s)\n", n_bg);
        mgs = realloc(mgs, (n_mg + n_bg) * sizeof(Metagenome));
        for (int i = 0; i < n_bg; i++) {
            char *prefixed = NULL;
            if (asprintf(&prefixed, "b_%s", bg[i].basename) < 0) {
                perror("asprintf"); exit(1);
            }
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
    /* with --read-hash every k-mer's vector also holds n_total min hashes */
    const int n_slots = read_hash_opt ? 2 * n_total : n_total;
    n_kmers = parse_kmer_tsv(kmer_file, h, n_slots, &seed, &kmer_list);
    fprintf(stderr, "loaded %d kmer(s) (k=%d)\n", n_kmers, seed);

    /* 2b. .kmdb targets: sorted 2-bit index of the query k-mers */
    int n_kmdb = 0;
    for (int i = 0; i < n_total; i++) n_kmdb += mgs[i].type == IS_KMDB;
    QKmer *qk = NULL;
    size_t n_qk = 0;
    if (n_kmdb > 0) {
#ifndef HAVE_ZSTD
        fprintf(stderr, "error: %d target(s) are .kmdb files, but this build has no "
                "zstd support; rebuild with libzstd\n", n_kmdb);
        return 1;
#else
        if (n_kmers > 0 && seed > KMDB_MAX_K) {
            fprintf(stderr, "error: .kmdb targets need k <= %d (kmer file has k=%d)\n",
                    KMDB_MAX_K, seed);
            return 1;
        }
        qk = malloc(((size_t)n_kmers + 1) * sizeof(QKmer));
        char rc[MAX_KMER_LEN];
        for (int i = 0; i < n_kmers; i++) {
            strncpy(rc, kmer_list[i], MAX_KMER_LEN - 1);
            rc[MAX_KMER_LEN - 1] = '\0';
            BIO_reverseComplement(rc);
            char *oriented = orient_string(kmer_list[i], rc, seed);
            uint64_t key;
            if (kmdb_encode(oriented, seed, &key) != 0) continue;   /* cannot happen: ACGT only */
            qk[n_qk].key    = key;
            qk[n_qk].counts = (unsigned int *)BIO_searchHash(h, oriented);
            n_qk++;
        }
        qsort(qk, n_qk, sizeof(QKmer), cmp_qkmer);
        fprintf(stderr, "%d .kmdb target(s): indexed %zu query k-mer(s)\n", n_kmdb, n_qk);
#endif
    }

    /* 3. scan all metagenomes in parallel (nothing to count -> skip; seed is -1) */
    int n_fastq = 0, n_paired = 0;
    for (int i = 0; i < n_total; i++) {
        if (mgs[i].type == IS_KMDB) continue;
        n_fastq++;
        n_paired += mgs[i].type == IS_PAIRED_END || mgs[i].type == IS_PAIRED_END_INTERLEAVE;
    }
    if (n_fastq > 0)
        fprintf(stderr, "pair dedup: %s for FASTQ targets (%d of %d are PE/PEI)\n",
                pair_dedup ? "ON (each k-mer at most once per read pair)"
                           : "OFF (--no-pair-dedup: every mate counted on its own)",
                n_paired, n_fastq);
    uint64_t *total_kmer_counts = calloc(n_total, sizeof(uint64_t));
    uint64_t *n_frags = calloc(n_total, sizeof(uint64_t));
    if (n_kmers > 0) {
        ksd_pool *pool = ksd_pool_create(num_threads, n_total + 1);
        for (int i = 0; i < n_total; i++) {
            ksd_job *job = malloc(sizeof(ksd_job));
            job->mg_idx           = i;
            job->mg               = &mgs[i];
            job->h                = h;
            job->seed             = seed;
            job->total_kmer_counts = total_kmer_counts;
            job->qk               = qk;
            job->n_qk             = n_qk;
            job->pair_dedup       = pair_dedup;
            job->have_fastq       = n_fastq > 0;
            job->read_hash        = read_hash_opt;
            job->n_total          = n_total;
            job->n_frags          = n_frags;
            ksd_pool_submit(pool, job);
        }
        ksd_pool_wait_and_destroy(pool);
        fprintf(stderr, "done scanning\n");
    } else {
        fprintf(stderr, "WARNING: no kmers to count, writing header-only output\n");
    }

    /* 4. write gzipped TSV output */
    gzFile out = gzopen(outfile, "wb9");
    if (!out) {
        fprintf(stderr, "cannot open output %s: %s\n", outfile, strerror(errno));
        return 1;
    }

    /* header row */
    gzprintf(out, "#kmer");
    for (int i = 0; i < n_total; i++) gzprintf(out, "\t%s", mgs[i].basename);
    gzprintf(out, "\n");

    /* --read-hash: second file, same rows and columns, min read hashes */
    gzFile rh_out = NULL;
    char *rh_path = NULL;
    if (read_hash_opt) {
        const char *suf = ".kmer_hits.tsv.gz";
        size_t ol = strlen(outfile), sl = strlen(suf);
        if (ol > sl && strcmp(outfile + ol - sl, suf) == 0) {
            if (asprintf(&rh_path, "%.*s.kmer_readhash.tsv.gz", (int)(ol - sl), outfile) < 0) {
                perror("asprintf"); return 1;
            }
        } else if (asprintf(&rh_path, "%s.readhash.tsv.gz", outfile) < 0) {
            perror("asprintf"); return 1;
        }
        rh_out = gzopen(rh_path, "wb9");
        if (!rh_out) {
            fprintf(stderr, "cannot open output %s: %s\n", rh_path, strerror(errno));
            return 1;
        }
        gzprintf(rh_out, "#kmer");
        for (int i = 0; i < n_total; i++) gzprintf(rh_out, "\t%s", mgs[i].basename);
        gzprintf(rh_out, "\n");
    }

    /* data rows — in the same order as the input file */
    char kmer_rc[MAX_KMER_LEN];
    for (int i = 0; i < n_kmers; i++) {
        strncpy(kmer_rc, kmer_list[i], MAX_KMER_LEN - 1);
        kmer_rc[MAX_KMER_LEN - 1] = '\0';
        BIO_reverseComplement(kmer_rc);
        char *oriented = orient_string(kmer_list[i], kmer_rc, seed);
        unsigned int *counts = (unsigned int *)BIO_searchHash(h, oriented);

        gzprintf(out, "%s", kmer_list[i]);
        for (int j = 0; j < n_total; j++)
            gzprintf(out, "\t%u", counts ? counts[j] : 0);
        gzprintf(out, "\n");
        if (rh_out) {
            gzprintf(rh_out, "%s", kmer_list[i]);
            for (int j = 0; j < n_total; j++)
                gzprintf(rh_out, "\t%u", counts ? counts[n_total + j] : 0);
            gzprintf(rh_out, "\n");
        }
        free(kmer_list[i]);
    }

    /* summary row: total_evaluated = sum of the counts of ALL k-mers of the
       sample (ACGT-only windows, after pair dedup), see the header comment */
    gzprintf(out, "total_evaluated");
    for (int i = 0; i < n_total; i++)
        gzprintf(out, "\t%" PRIu64, total_kmer_counts[i]);
    gzprintf(out, "\n");

    gzclose(out);
    fprintf(stderr, "output written to %s\n", outfile);
    if (rh_out) {
        gzprintf(rh_out, "n_fragments");
        for (int i = 0; i < n_total; i++) gzprintf(rh_out, "\t%" PRIu64, n_frags[i]);
        gzprintf(rh_out, "\n");
        if (gzclose(rh_out) != Z_OK) {
            fprintf(stderr, "error: writing %s failed\n", rh_path);
            return 1;
        }
        fprintf(stderr, "read hashes written to %s\n", rh_path);
    }
    free(rh_path);
    free(n_frags);

    /* cleanup */
    free(qk);
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