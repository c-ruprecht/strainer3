/*
 * kmer_block_detect.c
 *
 * Block-aware strain detection directly from FASTQ (no kmdb).
 *
 * Input k-mers come from kmer_scrub_filter.py --reference (the .kmer_loci.tsv.gz
 * table): every rare k-mer with every locus it occupies in the reference, plus
 * segment_id / block_id / family_id. A segment is a run of k-mers at consecutive
 * reference positions, so a read that truly comes from that region must hit ALL
 * segment k-mers inside its footprint, at read offsets that line up on one
 * diagonal.
 *
 * Per read (each mate separately; split at N):
 *   1. every k-mer hit proposes one candidate per locus:
 *        strand + : diag = pos - j        strand - : diag = pos + j
 *      (pos = reference start of the k-mer, j = offset of the k-mer in the read)
 *   2. hits are grouped by (segment, strand, diag). A group is ACCEPTED only if
 *        n_obs == n_exp   (complete: every segment k-mer inside the read
 *                          footprint was seen; a read may only be cut by its
 *                          own ends, never in the middle)
 *        n_obs >= min(min_hits, segment length)
 *      Anything else (mismatch, sequencing error, relative strain) is rejected.
 *   3. placement: accepted groups are summed per (contig, strand, diag). The
 *      placement with the most k-mers wins. If several placements tie (read lies
 *      completely inside a repeat), the read is AMBIGUOUS and is only counted per
 *      repeat family.
 *
 * Per fragment (R1 + R2), counts are de-duplicated, so overlapping mates are
 * never counted twice:
 *   kmer_hits      fragments with an accepted hit on the k-mer
 *   block frags    fragments placed on the block
 *   block starts   reads (mates) whose leftmost reference coordinate falls on a
 *                  k-mer position of the block (exposure = n_kmers; Poisson)
 *   block best     every placed read counted ONCE, for one block: single-copy blocks first (a read
 *                  that touches a unique block counts there, not on a repeat copy), then the block it
 *                  covers most completely (largest share of its k-mers), ties -> the block closest
 *                  to the read's 5' end. Exposure per block follows from the block layout and the
 *                  read-length histogram (block_coverage_depth.py), so the counts stay Poisson
 *                  with ~(n_kmers + L - k) start positions per isolated block instead of n_kmers
 *   family amb.    fragments that hit a repeat family but could not be placed
 *
 * Outputs (-O dir, -n basename):
 *   <n>.kmer_hits.tsv.gz       #kmer, one column per sample, last row total_evaluated
 *                              (same layout as kmer_strain_detect)
 *   <n>.block_hits.tsv.gz      block_id sample frags starts best (full grid)
 *   <n>.read_lengths.tsv.gz    sample length n   (N-free read pieces >= k; exposure of 'best')
 *   <n>.family_hits.tsv.gz     family_id sample frags_ambiguous best_ambiguous (repeat families only)
 *                              best_ambiguous: reads inside a repeat that fit several copies equally,
 *                              one count per read -> family count = sum(best of its copy blocks) + this
 *   <n>.detect_stats.tsv       per-sample read / group statistics
 *
 * Usage:
 *   kmer_block_detect -k <name.kmer_loci.tsv.gz> -s <samples.tsv> -O <outdir> -n <name>
 *                     [-t threads] [-m min_hits]
 *
 *   samples.tsv: sample_name<TAB>R1.fastq.gz[<TAB>R2.fastq.gz]   (# = comment)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <errno.h>
#include <unistd.h>
#include <pthread.h>
#include <sys/stat.h>
#include <zlib.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

#define MAX_LINE 65536

/* ── types ──────────────────────────────────────────────────────────────── */

typedef struct {
    uint32_t seg, block, family, contig, kmer;
    int64_t  pos;          /* 0-based start of the k-mer on the reference forward strand */
    uint8_t  fwd_canon;    /* 1 = ref forward k-mer is canonical, 0 = its rc is, 2 = palindrome */
} Locus;

typedef struct {
    uint32_t contig, family;
    int64_t  start, end;   /* first / last k-mer start position (inclusive) */
    uint32_t first;        /* index of the locus at position 'start' in loci[] */
    int      used;
} Segment;

typedef struct {
    uint64_t canon;
    uint32_t loc_off, loc_n;
    char    *str;          /* as written in the input file */
} Kmer;

typedef struct {
    uint32_t locus, kmer;
    int64_t  diag;
    int8_t   strand;
} Hit;

typedef struct {
    uint32_t contig, seg;
    int8_t   strand;
    int64_t  diag;
    uint32_t h0, h1;       /* hit range [h0, h1) */
    uint32_t nobs;
} Group;

typedef struct {
    char *name, *r1, *r2;
} Sample;

typedef struct {
    uint32_t *kmer_hits;     /* [n_kmers] */
    uint32_t *block_frags;   /* [max_block + 1] */
    uint32_t *block_starts;  /* [max_block + 1] */
    uint32_t *block_best;    /* [max_block + 1] */
    uint64_t *len_hist;      /* [MAX_READ_LEN + 1], N-free pieces >= k */
    uint32_t *fam_amb;       /* [max_family + 1] */
    uint32_t *fam_best;      /* [max_family + 1] ambiguous READS (pieces), one count each */
    uint64_t reads, fragments, kmer_positions, frags_with_hits, frags_placed,
             groups_accepted, groups_rejected, pieces_ambiguous, pieces_placed;
    int ok;
} SampleRes;

typedef struct {
    Hit    *hits;   size_t nh, hcap;
    Group  *grp;    size_t ng, gcap;
    uint32_t *fk;   size_t nfk, fkcap;   /* fragment k-mers */
    uint32_t *fb;   size_t nfb, fbcap;   /* fragment blocks */
    uint32_t *ff;   size_t nff, ffcap;   /* fragment ambiguous families */
    struct { uint32_t b, cnt; int64_t mn, mx; } *bs; size_t nbs, bscap;   /* per-block stats of a read */
} Work;

/* ── globals (read-only after loading) ──────────────────────────────────── */

static int       K = 0;
static uint64_t  KMASK = 0;
static unsigned char NT[256];

static Kmer     *kmers = NULL;   static uint32_t n_kmers = 0, cap_kmers = 0;
static Locus    *loci  = NULL;   static uint32_t n_loci = 0;
static uint32_t *kmer_loci = NULL;               /* CSR: locus indices per k-mer */
static Segment  *segs  = NULL;   static uint32_t max_seg = 0;
static uint32_t  max_block = 0, max_family = 0;
static uint8_t  *block_used = NULL, *family_repeat = NULL;
static uint32_t *block_nk = NULL;                /* k-mers (loci) per block */
static uint8_t  *block_multi = NULL;             /* 1 = block of multi-copy k-mers (repeat) */
#define MAX_READ_LEN 4096
static char    **contigs = NULL; static uint32_t n_contigs = 0, cap_contigs = 0;
static int       min_hits = 1;

/* k-mer hash: canonical 2-bit code -> k-mer index */
static uint64_t *hkeys = NULL;
static int32_t  *hvals = NULL;
static uint64_t  hcap = 0, hmask = 0, hcount = 0;

/* ── helpers ────────────────────────────────────────────────────────────── */

static void die(const char *msg) { fprintf(stderr, "error: %s\n", msg); exit(1); }

static void *xrealloc(void *p, size_t n) {
    void *q = realloc(p, n ? n : 1);
    if (!q) die("out of memory");
    return q;
}

#define GROW(arr, n, cap, type) do {                                   \
    if ((n) >= (cap)) { (cap) = (cap) ? (cap) * 2 : 256;               \
        (arr) = (type *)xrealloc((arr), (cap) * sizeof(type)); }       \
} while (0)

static void init_nt(void) {
    memset(NT, 4, sizeof(NT));
    NT['A'] = NT['a'] = 0; NT['C'] = NT['c'] = 1;
    NT['G'] = NT['g'] = 2; NT['T'] = NT['t'] = 3;
}

static inline uint64_t mix64(uint64_t x) {
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33; return x;
}

static inline int32_t hget(uint64_t key) {
    uint64_t i = mix64(key) & hmask;
    while (hvals[i] >= 0) {
        if (hkeys[i] == key) return hvals[i];
        i = (i + 1) & hmask;
    }
    return -1;
}

static void hrehash(uint64_t newcap) {
    uint64_t *ok = hkeys; int32_t *ov = hvals; uint64_t oc = hcap;
    hcap = newcap; hmask = newcap - 1;
    hkeys = (uint64_t *)xrealloc(NULL, hcap * sizeof(uint64_t));
    hvals = (int32_t *)xrealloc(NULL, hcap * sizeof(int32_t));
    for (uint64_t i = 0; i < hcap; i++) hvals[i] = -1;
    for (uint64_t i = 0; i < oc; i++) if (ov[i] >= 0) {
        uint64_t j = mix64(ok[i]) & hmask;
        while (hvals[j] >= 0) j = (j + 1) & hmask;
        hkeys[j] = ok[i]; hvals[j] = ov[i];
    }
    free(ok); free(ov);
}

/* returns existing index, or inserts v and returns -1 */
static int32_t hput(uint64_t key, int32_t v) {
    if ((hcount + 1) * 2 > hcap) hrehash(hcap ? hcap * 2 : (1u << 16));
    uint64_t i = mix64(key) & hmask;
    while (hvals[i] >= 0) {
        if (hkeys[i] == key) return hvals[i];
        i = (i + 1) & hmask;
    }
    hkeys[i] = key; hvals[i] = v; hcount++;
    return -1;
}

static int encode(const char *s, uint64_t *out) {
    uint64_t x = 0;
    for (int i = 0; i < K; i++) {
        int c = NT[(unsigned char)s[i]];
        if (c > 3) return 0;
        x = (x << 2) | (uint64_t)c;
    }
    *out = x;
    return 1;
}

static uint64_t rc_code(uint64_t x) {
    uint64_t r = 0;
    for (int i = 0; i < K; i++) { r = (r << 2) | (3 - (x & 3)); x >>= 2; }
    return r;
}

static uint32_t contig_index(const char *name) {
    static uint32_t last = 0;
    if (n_contigs && strcmp(contigs[last], name) == 0) return last;
    for (uint32_t i = 0; i < n_contigs; i++)
        if (strcmp(contigs[i], name) == 0) { last = i; return i; }
    GROW(contigs, n_contigs, cap_contigs, char *);
    contigs[n_contigs] = strdup(name);
    last = n_contigs;
    return n_contigs++;
}

static int parse_bool(const char *s) {
    return (s[0] == 'T' || s[0] == 't' || s[0] == '1' || s[0] == 'Y' || s[0] == 'y');
}

static void mkdir_p(const char *path) {
    char tmp[4096];
    snprintf(tmp, sizeof(tmp), "%s", path);
    for (char *p = tmp + 1; *p; p++) {
        if (*p == '/') { *p = '\0'; mkdir(tmp, 0755); *p = '/'; }
    }
    if (mkdir(tmp, 0755) != 0 && errno != EEXIST) {
        fprintf(stderr, "cannot create %s: %s\n", path, strerror(errno));
        exit(1);
    }
}

/* ── loading the k-mer loci table ───────────────────────────────────────── */

static int cmp_locus_seg_pos(const void *a, const void *b) {
    const Locus *x = (const Locus *)a, *y = (const Locus *)b;
    if (x->seg != y->seg) return x->seg < y->seg ? -1 : 1;
    if (x->pos != y->pos) return x->pos < y->pos ? -1 : 1;
    return 0;
}

static void load_loci(const char *file) {
    gzFile fp = gzopen(file, "r");
    if (!fp) { fprintf(stderr, "cannot open %s: %s\n", file, strerror(errno)); exit(1); }

    char *line = (char *)xrealloc(NULL, MAX_LINE);
    if (!gzgets(fp, line, MAX_LINE)) die("k-mer loci file is empty");

    /* header -> column indices */
    const char *want[] = { "#kmer", "contig_id", "kmer_position", "ref_rc",
                           "segment_id", "block_id", "family_id" };
    int col[7] = { -1, -1, -1, -1, -1, -1, -1 };
    {
        int c = 0;
        char *save = NULL;
        for (char *tok = strtok_r(line, "\t\r\n", &save); tok; tok = strtok_r(NULL, "\t\r\n", &save), c++)
            for (int w = 0; w < 7; w++) if (strcmp(tok, want[w]) == 0) col[w] = c;
    }
    for (int w = 0; w < 7; w++) if (col[w] < 0) {
        fprintf(stderr, "error: column '%s' missing in %s\n"
                        "       (generate it with kmer_scrub_filter.py --reference)\n", want[w], file);
        exit(1);
    }

    uint32_t cap_loci = 0;
    char *f[64];
    while (gzgets(fp, line, MAX_LINE)) {
        if (line[0] == '\n' || line[0] == '\r' || line[0] == '\0') continue;
        int nf = 0;
        char *save = NULL;
        for (char *tok = strtok_r(line, "\t\r\n", &save); tok && nf < 64; tok = strtok_r(NULL, "\t\r\n", &save))
            f[nf++] = tok;
        int need = 0;
        for (int w = 0; w < 7; w++) if (col[w] > need) need = col[w];
        if (nf <= need) continue;

        const char *ks = f[col[0]];
        if (K == 0) {
            K = (int)strlen(ks);
            if (K < 1 || K > 32) die("k must be between 1 and 32");
            KMASK = (K == 32) ? ~0ULL : ((1ULL << (2 * K)) - 1);
        }
        if ((int)strlen(ks) != K) die("k-mers of different length in loci file");

        uint64_t w_code, rc_w;
        if (!encode(ks, &w_code)) {
            fprintf(stderr, "error: non-ACGT k-mer in %s: %s\n"
                            "       (ambiguity codes in the reference; re-run kmer_scrub_filter.py, which now drops them)\n",
                    file, ks);
            exit(1);
        }
        rc_w = rc_code(w_code);
        uint64_t canon = w_code < rc_w ? w_code : rc_w;

        int32_t ki = hput(canon, (int32_t)n_kmers);
        if (ki < 0) {
            GROW(kmers, n_kmers, cap_kmers, Kmer);
            kmers[n_kmers].canon = canon;
            kmers[n_kmers].str = strdup(ks);
            kmers[n_kmers].loc_n = 0;
            ki = (int32_t)n_kmers++;
        }

        int ref_rc = parse_bool(f[col[3]]);
        uint64_t ref_fwd = ref_rc ? rc_w : w_code;
        uint64_t ref_rev = ref_rc ? w_code : rc_w;

        GROW(loci, n_loci, cap_loci, Locus);
        Locus *L = &loci[n_loci++];
        L->kmer      = (uint32_t)ki;
        L->contig    = contig_index(f[col[1]]);
        L->pos       = strtoll(f[col[2]], NULL, 10);
        L->seg       = (uint32_t)strtoul(f[col[4]], NULL, 10);
        L->block     = (uint32_t)strtoul(f[col[5]], NULL, 10);
        L->family    = (uint32_t)strtoul(f[col[6]], NULL, 10);
        L->fwd_canon = (ref_fwd == ref_rev) ? 2 : (ref_fwd < ref_rev ? 1 : 0);
        kmers[ki].loc_n++;

        if (L->seg > max_seg) max_seg = L->seg;
        if (L->block > max_block) max_block = L->block;
        if (L->family > max_family) max_family = L->family;
    }
    gzclose(fp);
    free(line);
    if (n_loci == 0) {
        /* strain without detectable blocks (kmer_scrub_filter.py wrote an empty loci file): still scan the
           reads so detect_stats / read_lengths are complete, and write header-only hit tables */
        fprintf(stderr, "WARNING: no loci in %s: writing empty hit tables (strain not detectable)\n", file);
        if (K == 0) { K = 31; KMASK = (1ULL << (2 * K)) - 1; }
        hrehash(1u << 4);
    }

    /* sort loci by segment, position and build segment geometry */
    qsort(loci, n_loci, sizeof(Locus), cmp_locus_seg_pos);
    segs = (Segment *)calloc(max_seg + 1, sizeof(Segment));
    block_used = (uint8_t *)calloc(max_block + 1, 1);
    family_repeat = (uint8_t *)calloc(max_family + 1, 1);
    uint32_t *fam_nseg = (uint32_t *)calloc(max_family + 1, sizeof(uint32_t));

    block_nk = (uint32_t *)calloc(max_block + 1, sizeof(uint32_t));
    for (uint32_t i = 0; i < n_loci; i++) {
        Locus *L = &loci[i];
        Segment *S = &segs[L->seg];
        block_used[L->block] = 1;
        block_nk[L->block]++;
        if (!S->used) {
            S->used = 1; S->contig = L->contig; S->family = L->family;
            S->start = S->end = L->pos; S->first = i;
            fam_nseg[L->family]++;
        } else {
            if (L->contig != S->contig || L->pos != S->end + 1) {
                fprintf(stderr, "error: segment %u is not contiguous (contig %s pos %" PRId64 ")\n",
                        L->seg, contigs[L->contig], L->pos);
                exit(1);
            }
            S->end = L->pos;
        }
    }
    for (uint32_t f2 = 0; f2 <= max_family; f2++) family_repeat[f2] = fam_nseg[f2] > 1;
    free(fam_nseg);

    /* CSR: k-mer -> locus indices (after sorting) */
    uint32_t off = 0;
    for (uint32_t i = 0; i < n_kmers; i++) { kmers[i].loc_off = off; off += kmers[i].loc_n; kmers[i].loc_n = 0; }
    kmer_loci = (uint32_t *)xrealloc(NULL, (size_t)n_loci * sizeof(uint32_t));
    for (uint32_t i = 0; i < n_loci; i++) {
        Kmer *km = &kmers[loci[i].kmer];
        kmer_loci[km->loc_off + km->loc_n++] = i;
    }

    block_multi = (uint8_t *)calloc(max_block + 1, 1);
    for (uint32_t i = 0; i < n_loci; i++)
        if (kmers[loci[i].kmer].loc_n > 1) block_multi[loci[i].block] = 1;

    uint32_t n_segs = 0, n_blocks = 0, n_rep = 0, multi = 0;
    for (uint32_t s = 0; s <= max_seg; s++) n_segs += segs[s].used;
    for (uint32_t b = 0; b <= max_block; b++) n_blocks += block_used[b];
    for (uint32_t f2 = 0; f2 <= max_family; f2++) n_rep += family_repeat[f2];
    for (uint32_t i = 0; i < n_kmers; i++) multi += kmers[i].loc_n > 1;
    fprintf(stderr, "loaded %u k-mers (k=%d), %u loci, %u multi-copy k-mers, %u segments, "
                    "%u blocks, %u repeat families, %u contigs\n",
            n_kmers, K, n_loci, multi, n_segs, n_blocks, n_rep, n_contigs);
}

/* ── per-read logic ─────────────────────────────────────────────────────── */

static int cmp_hit(const void *a, const void *b) {
    const Hit *x = (const Hit *)a, *y = (const Hit *)b;
    uint32_t sx = loci[x->locus].seg, sy = loci[y->locus].seg;
    if (sx != sy) return sx < sy ? -1 : 1;
    if (x->strand != y->strand) return x->strand < y->strand ? -1 : 1;
    if (x->diag != y->diag) return x->diag < y->diag ? -1 : 1;
    return 0;
}

static int cmp_group_place(const void *a, const void *b) {
    const Group *x = (const Group *)a, *y = (const Group *)b;
    if (x->contig != y->contig) return x->contig < y->contig ? -1 : 1;
    if (x->strand != y->strand) return x->strand < y->strand ? -1 : 1;
    if (x->diag != y->diag) return x->diag < y->diag ? -1 : 1;
    return 0;
}

static inline int same_place(const Group *a, const Group *b) {
    return a->contig == b->contig && a->strand == b->strand && a->diag == b->diag;
}

static inline void footprint(int8_t strand, int64_t diag, int plen, int64_t *lo, int64_t *hi) {
    int64_t span = plen - K;
    if (strand > 0) { *lo = diag;        *hi = diag + span; }
    else            { *lo = diag - span; *hi = diag; }
}

static inline void push_u32(uint32_t **arr, size_t *n, size_t *cap, uint32_t v) {
    if (*n >= *cap) { *cap = *cap ? *cap * 2 : 256; *arr = (uint32_t *)xrealloc(*arr, *cap * sizeof(uint32_t)); }
    (*arr)[(*n)++] = v;
}

/* evaluate the hits of one N-free piece of a read */
static void process_piece(Work *w, int plen, int whole_read, SampleRes *r) {
    if (w->nh == 0) return;
    qsort(w->hits, w->nh, sizeof(Hit), cmp_hit);

    /* 1. groups by (segment, strand, diag) and completeness check */
    w->ng = 0;
    size_t i = 0;
    while (i < w->nh) {
        size_t j = i + 1;
        uint32_t seg = loci[w->hits[i].locus].seg;
        while (j < w->nh && loci[w->hits[j].locus].seg == seg &&
               w->hits[j].strand == w->hits[i].strand && w->hits[j].diag == w->hits[i].diag) j++;

        const Segment *S = &segs[seg];
        int64_t lo, hi;
        footprint(w->hits[i].strand, w->hits[i].diag, plen, &lo, &hi);
        int64_t a = lo > S->start ? lo : S->start;
        int64_t b = hi < S->end ? hi : S->end;
        int64_t n_exp = b >= a ? b - a + 1 : 0;
        int64_t n_obs = (int64_t)(j - i);
        int64_t seg_n = S->end - S->start + 1;
        int64_t need  = min_hits < seg_n ? min_hits : seg_n;

        if (n_obs == n_exp && n_obs >= need) {
            if (w->ng >= w->gcap) { w->gcap = w->gcap ? w->gcap * 2 : 64;
                                    w->grp = (Group *)xrealloc(w->grp, w->gcap * sizeof(Group)); }
            Group *g = &w->grp[w->ng++];
            g->contig = S->contig; g->seg = seg;
            g->strand = w->hits[i].strand; g->diag = w->hits[i].diag;
            g->h0 = (uint32_t)i; g->h1 = (uint32_t)j; g->nobs = (uint32_t)n_obs;
            r->groups_accepted++;
        } else {
            r->groups_rejected++;
        }
        i = j;
    }
    if (w->ng == 0) return;

    /* 2. placement: total accepted k-mers per (contig, strand, diag) */
    qsort(w->grp, w->ng, sizeof(Group), cmp_group_place);
    uint32_t best = 0; int n_best = 0;
    size_t g0 = 0;
    while (g0 < w->ng) {
        size_t g1 = g0 + 1;
        uint32_t score = w->grp[g0].nobs;
        while (g1 < w->ng && same_place(&w->grp[g1], &w->grp[g0])) score += w->grp[g1++].nobs;
        if (score > best) { best = score; n_best = 1; }
        else if (score == best) n_best++;
        g0 = g1;
    }

    /* 3. record */
    int ambiguous = n_best > 1;
    if (ambiguous) r->pieces_ambiguous++; else r->pieces_placed++;
    int amb_counted = 0;

    g0 = 0;
    while (g0 < w->ng) {
        size_t g1 = g0 + 1;
        uint32_t score = w->grp[g0].nobs;
        while (g1 < w->ng && same_place(&w->grp[g1], &w->grp[g0])) score += w->grp[g1++].nobs;
        if (score == best) {
            int64_t flo = 0, fhi = 0;
            footprint(w->grp[g0].strand, w->grp[g0].diag, plen, &flo, &fhi);
            for (size_t g = g0; g < g1; g++) {
                const Group *G = &w->grp[g];
                for (uint32_t h = G->h0; h < G->h1; h++) {
                    const Hit *H = &w->hits[h];
                    push_u32(&w->fk, &w->nfk, &w->fkcap, H->kmer);
                    if (!ambiguous) push_u32(&w->fb, &w->nfb, &w->fbcap, loci[H->locus].block);
                }
                if (ambiguous) {
                    push_u32(&w->ff, &w->nff, &w->ffcap, segs[G->seg].family);
                    if (!amb_counted) { r->fam_best[segs[G->seg].family]++; amb_counted = 1; }
                } else if (whole_read) {
                    const Segment *S = &segs[G->seg];
                    if (flo >= S->start && flo <= S->end)
                        r->block_starts[loci[S->first + (uint32_t)(flo - S->start)].block]++;
                }
            }
            if (!ambiguous) {
                /* one block per read: single-copy blocks first (a read touching a unique block counts there),
                   then the most completely covered block, ties -> closest to the read's 5' end */
                w->nbs = 0;
                for (size_t g = g0; g < g1; g++) {
                    const Group *G = &w->grp[g];
                    for (uint32_t h = G->h0; h < G->h1; h++) {
                        const Locus *Lc = &loci[w->hits[h].locus];
                        size_t q = 0;
                        while (q < w->nbs && w->bs[q].b != Lc->block) q++;
                        if (q == w->nbs) {
                            if (w->nbs >= w->bscap) {
                                w->bscap = w->bscap ? w->bscap * 2 : 16;
                                w->bs = xrealloc(w->bs, w->bscap * sizeof(*w->bs));
                            }
                            w->bs[q].b = Lc->block; w->bs[q].cnt = 0;
                            w->bs[q].mn = Lc->pos; w->bs[q].mx = Lc->pos;
                            w->nbs++;
                        }
                        w->bs[q].cnt++;
                        if (Lc->pos < w->bs[q].mn) w->bs[q].mn = Lc->pos;
                        if (Lc->pos > w->bs[q].mx) w->bs[q].mx = Lc->pos;
                    }
                }
                int8_t strand = w->grp[g0].strand;
                size_t bi = 0;
                for (size_t q = 1; q < w->nbs; q++) {
                    int uq = !block_multi[w->bs[q].b], ub = !block_multi[w->bs[bi].b];
                    if (uq != ub) { if (uq) bi = q; continue; }
                    /* compare cnt_q / nk_q with cnt_bi / nk_bi without division */
                    uint64_t lhs = (uint64_t)w->bs[q].cnt * block_nk[w->bs[bi].b];
                    uint64_t rhs = (uint64_t)w->bs[bi].cnt * block_nk[w->bs[q].b];
                    if (lhs > rhs) { bi = q; continue; }
                    if (lhs == rhs) {
                        int64_t dq = strand > 0 ? w->bs[q].mn - flo : fhi - w->bs[q].mx;
                        int64_t db = strand > 0 ? w->bs[bi].mn - flo : fhi - w->bs[bi].mx;
                        if (dq < db) bi = q;
                    }
                }
                if (w->nbs) r->block_best[w->bs[bi].b]++;
            }
        }
        g0 = g1;
    }
}

static void scan_read(Work *w, const char *s, int L, SampleRes *r) {
    r->reads++;
    if (L < K) return;
    uint64_t fwd = 0, rcv = 0;
    int l = 0, ps = 0;
    const int shift = 2 * (K - 1);
    w->nh = 0;

    for (int i = 0; i < L; i++) {
        int c = NT[(unsigned char)s[i]];
        if (c > 3) {
            if (i - ps >= K) r->len_hist[i - ps > MAX_READ_LEN ? MAX_READ_LEN : i - ps]++;
            process_piece(w, i - ps, 0, r);
            w->nh = 0; l = 0; ps = i + 1;
            continue;
        }
        fwd = ((fwd << 2) | (uint64_t)c) & KMASK;
        rcv = (rcv >> 2) | ((uint64_t)(3 - c) << shift);
        if (++l < K) continue;

        r->kmer_positions++;
        uint64_t canon = fwd < rcv ? fwd : rcv;
        int32_t ki = hget(canon);
        if (ki < 0) continue;

        int read_fc = (fwd == rcv) ? 2 : (fwd < rcv ? 1 : 0);
        int64_t j = (int64_t)(i - K + 1 - ps);
        const Kmer *km = &kmers[ki];
        for (uint32_t q = 0; q < km->loc_n; q++) {
            uint32_t li = kmer_loci[km->loc_off + q];
            const Locus *Lc = &loci[li];
            int8_t strands[2]; int ns = 0;
            if (Lc->fwd_canon == 2 || read_fc == 2) { strands[ns++] = 1; strands[ns++] = -1; }
            else strands[ns++] = (Lc->fwd_canon == read_fc) ? 1 : -1;
            for (int t = 0; t < ns; t++) {
                if (w->nh >= w->hcap) { w->hcap = w->hcap ? w->hcap * 2 : 256;
                                        w->hits = (Hit *)xrealloc(w->hits, w->hcap * sizeof(Hit)); }
                Hit *H = &w->hits[w->nh++];
                H->locus = li; H->kmer = (uint32_t)ki; H->strand = strands[t];
                H->diag = strands[t] > 0 ? Lc->pos - j : Lc->pos + j;
            }
        }
    }
    if (L - ps >= K) r->len_hist[L - ps > MAX_READ_LEN ? MAX_READ_LEN : L - ps]++;
    process_piece(w, L - ps, ps == 0, r);
    w->nh = 0;
}

static int cmp_u32(const void *a, const void *b) {
    uint32_t x = *(const uint32_t *)a, y = *(const uint32_t *)b;
    return (x > y) - (x < y);
}

/* add the de-duplicated fragment counts and reset the fragment buffers */
static void finish_fragment(Work *w, SampleRes *r) {
    r->fragments++;
    if (w->nfk) r->frags_with_hits++;
    if (w->nfb) r->frags_placed++;
    if (w->nfk) {
        qsort(w->fk, w->nfk, sizeof(uint32_t), cmp_u32);
        for (size_t i = 0; i < w->nfk; i++) if (i == 0 || w->fk[i] != w->fk[i - 1]) r->kmer_hits[w->fk[i]]++;
    }
    if (w->nfb) {
        qsort(w->fb, w->nfb, sizeof(uint32_t), cmp_u32);
        for (size_t i = 0; i < w->nfb; i++) if (i == 0 || w->fb[i] != w->fb[i - 1]) r->block_frags[w->fb[i]]++;
    }
    if (w->nff) {
        qsort(w->ff, w->nff, sizeof(uint32_t), cmp_u32);
        for (size_t i = 0; i < w->nff; i++) if (i == 0 || w->ff[i] != w->ff[i - 1]) r->fam_amb[w->ff[i]]++;
    }
    w->nfk = w->nfb = w->nff = 0;
}

/* ── per-sample job ─────────────────────────────────────────────────────── */

static void run_sample(const Sample *sm, SampleRes *r) {
    r->kmer_hits    = (uint32_t *)calloc(n_kmers, sizeof(uint32_t));
    r->block_frags  = (uint32_t *)calloc(max_block + 1, sizeof(uint32_t));
    r->block_starts = (uint32_t *)calloc(max_block + 1, sizeof(uint32_t));
    r->block_best   = (uint32_t *)calloc(max_block + 1, sizeof(uint32_t));
    r->len_hist     = (uint64_t *)calloc(MAX_READ_LEN + 1, sizeof(uint64_t));
    r->fam_amb      = (uint32_t *)calloc(max_family + 1, sizeof(uint32_t));
    r->fam_best     = (uint32_t *)calloc(max_family + 1, sizeof(uint32_t));

    gzFile f1 = gzopen(sm->r1, "r");
    if (!f1) { fprintf(stderr, "[%s] cannot open %s: %s\n", sm->name, sm->r1, strerror(errno)); return; }
    gzFile f2 = NULL;
    if (sm->r2) {
        f2 = gzopen(sm->r2, "r");
        if (!f2) { fprintf(stderr, "[%s] cannot open %s: %s\n", sm->name, sm->r2, strerror(errno)); gzclose(f1); return; }
    }
    gzbuffer(f1, 1 << 17);
    if (f2) gzbuffer(f2, 1 << 17);
    kseq_t *s1 = kseq_init(f1);
    kseq_t *s2 = f2 ? kseq_init(f2) : NULL;
    Work w; memset(&w, 0, sizeof(w));

    for (;;) {
        int l1 = kseq_read(s1);
        int l2 = s2 ? kseq_read(s2) : -1;
        if (l1 < 0) {
            if (s2 && l2 >= 0) fprintf(stderr, "[%s] WARNING: R2 has more reads than R1\n", sm->name);
            break;
        }
        if (s2 && l2 < 0) { fprintf(stderr, "[%s] WARNING: R1 has more reads than R2\n", sm->name); break; }
        scan_read(&w, s1->seq.s, l1, r);
        if (s2) scan_read(&w, s2->seq.s, l2, r);
        finish_fragment(&w, r);
    }

    kseq_destroy(s1); gzclose(f1);
    if (s2) { kseq_destroy(s2); gzclose(f2); }
    free(w.hits); free(w.grp); free(w.fk); free(w.fb); free(w.ff); free(w.bs);
    r->ok = 1;
    fprintf(stderr, "[%s] %" PRIu64 " fragments, %" PRIu64 " with accepted hits, "
                    "%" PRIu64 " groups accepted, %" PRIu64 " rejected, %" PRIu64 " ambiguous pieces\n",
            sm->name, r->fragments, r->frags_with_hits, r->groups_accepted, r->groups_rejected,
            r->pieces_ambiguous);
}

/* ── thread pool (atomic job counter) ───────────────────────────────────── */

typedef struct { Sample *samples; SampleRes *res; int n; int next; } Pool;

static void *worker(void *arg) {
    Pool *p = (Pool *)arg;
    for (;;) {
        int i = __sync_fetch_and_add(&p->next, 1);
        if (i >= p->n) break;
        run_sample(&p->samples[i], &p->res[i]);
    }
    return NULL;
}

/* ── samples file ───────────────────────────────────────────────────────── */

static int load_samples(const char *file, Sample **out) {
    FILE *fp = fopen(file, "r");
    if (!fp) { fprintf(stderr, "cannot open samples file %s: %s\n", file, strerror(errno)); exit(1); }
    Sample *s = NULL; size_t n = 0, cap = 0;
    char *line = NULL; size_t len = 0;
    while (getline(&line, &len, fp) != -1) {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
        char *f[3] = { NULL, NULL, NULL }; int nf = 0;
        char *save = NULL;
        for (char *tok = strtok_r(line, "\t\r\n", &save); tok && nf < 3; tok = strtok_r(NULL, "\t\r\n", &save))
            f[nf++] = tok;
        if (nf < 2) { fprintf(stderr, "skipping malformed samples line (need name<TAB>R1[<TAB>R2])\n"); continue; }
        GROW(s, n, cap, Sample);
        s[n].name = strdup(f[0]);
        s[n].r1 = strdup(f[1]);
        s[n].r2 = (nf > 2 && f[2][0]) ? strdup(f[2]) : NULL;
        n++;
    }
    free(line);
    fclose(fp);
    *out = s;
    return (int)n;
}

/* ── output ─────────────────────────────────────────────────────────────── */

static gzFile open_out(const char *dir, const char *name, const char *suffix) {
    char path[4096];
    snprintf(path, sizeof(path), "%s/%s%s", dir, name, suffix);
    gzFile g = gzopen(path, "wb6");
    if (!g) { fprintf(stderr, "cannot open output %s: %s\n", path, strerror(errno)); exit(1); }
    return g;
}

static void write_outputs(const char *dir, const char *name, Sample *sm, SampleRes *res, int ns) {
    /* k-mer hits, legacy layout */
    gzFile g = open_out(dir, name, ".kmer_hits.tsv.gz");
    gzprintf(g, "#kmer");
    for (int s = 0; s < ns; s++) gzprintf(g, "\t%s", sm[s].name);
    gzprintf(g, "\n");
    for (uint32_t i = 0; i < n_kmers; i++) {
        gzputs(g, kmers[i].str);
        for (int s = 0; s < ns; s++) gzprintf(g, "\t%u", res[s].kmer_hits[i]);
        gzputc(g, '\n');
    }
    gzprintf(g, "total_evaluated");
    for (int s = 0; s < ns; s++) gzprintf(g, "\t%" PRIu64, res[s].kmer_positions);
    gzprintf(g, "\n");
    gzclose(g);

    /* blocks, long format, full grid */
    g = open_out(dir, name, ".block_hits.tsv.gz");
    gzprintf(g, "block_id\tsample\tfrags\tstarts\tbest\n");
    for (int s = 0; s < ns; s++)
        for (uint32_t b = 0; b <= max_block; b++) if (block_used[b])
            gzprintf(g, "%u\t%s\t%u\t%u\t%u\n", b, sm[s].name, res[s].block_frags[b], res[s].block_starts[b],
                     res[s].block_best[b]);
    gzclose(g);

    /* read-length histogram (N-free pieces >= k), the exposure model of 'best' needs it */
    g = open_out(dir, name, ".read_lengths.tsv.gz");
    gzprintf(g, "sample\tlength\tn\n");
    for (int s = 0; s < ns; s++)
        for (int l = K; l <= MAX_READ_LEN; l++) if (res[s].len_hist[l])
            gzprintf(g, "%s\t%d\t%" PRIu64 "\n", sm[s].name, l, res[s].len_hist[l]);
    gzclose(g);

    /* repeat families: fragments that could not be placed on one copy */
    g = open_out(dir, name, ".family_hits.tsv.gz");
    gzprintf(g, "family_id\tsample\tfrags_ambiguous\tbest_ambiguous\n");
    for (int s = 0; s < ns; s++)
        for (uint32_t f2 = 0; f2 <= max_family; f2++) if (family_repeat[f2])
            gzprintf(g, "%u\t%s\t%u\t%u\n", f2, sm[s].name, res[s].fam_amb[f2], res[s].fam_best[f2]);
    gzclose(g);

    /* stats */
    char path[4096];
    snprintf(path, sizeof(path), "%s/%s.detect_stats.tsv", dir, name);
    FILE *fp = fopen(path, "w");
    if (!fp) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
    fprintf(fp, "sample\tfragments\treads\tkmer_positions\tfrags_with_hits\tfrags_placed\t"
                "groups_accepted\tgroups_rejected\tpieces_placed\tpieces_ambiguous\n");
    for (int s = 0; s < ns; s++) {
        SampleRes *r = &res[s];
        fprintf(fp, "%s\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64
                    "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\n",
                sm[s].name, r->fragments, r->reads, r->kmer_positions, r->frags_with_hits, r->frags_placed,
                r->groups_accepted, r->groups_rejected, r->pieces_placed, r->pieces_ambiguous);
    }
    fclose(fp);
}

/* ── main ───────────────────────────────────────────────────────────────── */

static void usage(void) {
    fprintf(stderr,
        "Usage: kmer_block_detect -k <name.kmer_loci.tsv.gz> -s <samples.tsv> -O <outdir> -n <name>\n"
        "                         [-t threads (4)] [-m min_hits (1)]\n\n"
        "  -k  k-mer loci table from kmer_scrub_filter.py --reference\n"
        "  -s  samples: name<TAB>R1.fastq[.gz][<TAB>R2.fastq[.gz]] per line (# = comment)\n"
        "  -O  output directory\n"
        "  -n  output basename\n"
        "  -t  threads (one sample per thread)\n"
        "  -m  min k-mers a read must hit in a segment (capped at the segment size).\n"
        "      Reads cut by their own end that only clip a segment edge are\n"
        "      rejected below this. Default 1.\n\n"
        "A read counts for a segment only if it hits every segment k-mer inside its\n"
        "footprint on one diagonal (no internal gaps). See the header of the source\n"
        "for the output columns.\n");
}

int main(int argc, char *argv[]) {
    char *kfile = NULL, *sfile = NULL, *outdir = NULL, *name = NULL;
    int threads = 4, c;

    while ((c = getopt(argc, argv, "k:s:O:n:t:m:h")) != -1) {
        switch (c) {
            case 'k': kfile = optarg; break;
            case 's': sfile = optarg; break;
            case 'O': outdir = optarg; break;
            case 'n': name = optarg; break;
            case 't': threads = atoi(optarg); break;
            case 'm': min_hits = atoi(optarg); break;
            case 'h': usage(); return 0;
            default:  usage(); return 1;
        }
    }
    if (!kfile || !sfile || !outdir || !name) { usage(); return 1; }
    if (threads < 1) threads = 1;
    if (min_hits < 1) min_hits = 1;

    init_nt();
    load_loci(kfile);

    Sample *samples = NULL;
    int ns = load_samples(sfile, &samples);
    if (ns == 0) die("no samples");
    fprintf(stderr, "loaded %d sample(s); min_hits=%d; threads=%d\n", ns, min_hits, threads);

    SampleRes *res = (SampleRes *)calloc(ns, sizeof(SampleRes));
    Pool pool = { samples, res, ns, 0 };
    int nt = threads < ns ? threads : ns;
    pthread_t *tid = (pthread_t *)xrealloc(NULL, nt * sizeof(pthread_t));
    for (int i = 0; i < nt; i++) pthread_create(&tid[i], NULL, worker, &pool);
    for (int i = 0; i < nt; i++) pthread_join(tid[i], NULL);
    free(tid);

    for (int s = 0; s < ns; s++) if (!res[s].ok) {
        fprintf(stderr, "error: sample %s failed\n", samples[s].name);
        return 1;
    }

    mkdir_p(outdir);
    write_outputs(outdir, name, samples, res, ns);
    fprintf(stderr, "outputs written to %s/%s.*\n", outdir, name);

    for (int s = 0; s < ns; s++) {
        free(res[s].kmer_hits); free(res[s].block_frags); free(res[s].block_starts); free(res[s].fam_amb);
        free(res[s].block_best); free(res[s].len_hist); free(res[s].fam_best);
        free(samples[s].name); free(samples[s].r1); free(samples[s].r2);
    }
    free(res); free(samples);
    for (uint32_t i = 0; i < n_kmers; i++) free(kmers[i].str);
    free(kmers); free(loci); free(kmer_loci); free(segs);
    free(block_used); free(family_repeat); free(block_nk); free(block_multi);
    for (uint32_t i = 0; i < n_contigs; i++) free(contigs[i]);
    free(contigs); free(hkeys); free(hvals);
    return 0;
}
