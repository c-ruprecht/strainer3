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
 * Optional, off by default (the default output is unchanged):
 *   -V  fragment veto. A rejected group is a CONTRADICTION if it is incomplete (k-mers missing inside the
 *       read) and holds a single-copy k-mer: nothing but a mismatch, indel or a relative's variant explains
 *       it. With -V a fragment with a contradiction in either mate is dropped as a whole; without it only
 *       that group is dropped and the rest of the read still counts.
 *   -c  context loci (kmer_scrub_filter.py --context_bp): single-copy reference k-mers next to the rare
 *       segments. They are never counted or reported, but they join the rare k-mers into longer validation
 *       segments, so a read is checked over its whole length, not only inside one 31-61 k-mer chunk.
 *   -Q  quality-aware rejection: a gap in the run of k-mers that one wrong base can make, next to a base
 *       below the given quality, is a sequencing error: the read keeps its other k-mers, no contradiction.
 *   -P  pair mode for unmerged R1 / R2 (fastp --dedup only). Mates overlapping by >= k bases are booked as one
 *       read of the fragment's length (k-mer positions, read lengths, one best / start per fragment); a
 *       broken run in one mate is rescued if the other mate reads the same reference k-mers cleanly; a mate
 *       with context k-mers only is still checked and placed; mates on one contig must face each other
 *       within -I or the fragment is a contradiction. In pair mode compare samples on
 *       frags_contradicted / frags_hit. Recommended: -P -V -Q 20 -c <context of 300 bp>.
 *   Several sample-sheet rows with one name are pooled: give the merged reads on one row and the unmerged
 *   R1 / R2 pairs on another, so unmerged mates stay one fragment (counted once, vetoed together).
 *
 * Outputs (-O dir, -n basename):
 *   <n>.kmer_hits.tsv.gz       #kmer, one column per sample, last row total_evaluated
 *                              (same layout as kmer_strain_detect)
 *   <n>.block_hits.tsv.gz      block_id sample frags starts best (full grid)
 *   <n>.read_lengths.tsv.gz    sample length n   (N-free read pieces >= k; exposure of 'best')
 *   <n>.family_hits.tsv.gz     family_id sample frags_ambiguous best_ambiguous (repeat families only)
 *                              best_ambiguous: reads inside a repeat that fit several copies equally,
 *                              one count per read -> family count = sum(best of its copy blocks) + this
 *   <n>.kmer_hits_any.tsv.gz   same layout as kmer_hits, counting every fragment that holds the k-mer whether
 *                              its read was accepted or not (the non-aware count)
 *   block_hits 'refused'       fragments that hit the block's rare k-mers but contradict the strain. A block
 *                              with refused > 0 and frags == 0 has reads, and all of them disagree with the
 *                              reference nearby: evidence against the strain that an empty block does not give
 *   <n>.detect_stats.tsv       per-sample read / group statistics. reads_contradicted / reads_hit is the
 *                              share of reads with rare hits that contradict the reference;
 *                              reads_contra_expected / reads_hit is what sequencing errors alone would give
 *                              (from the base qualities), so the ratio of the two is comparable between
 *                              libraries. pairs_* : mates both placed, concordant = one contig, opposite
 *                              strands, outer distance <= -I
 *
 * Usage:
 *   kmer_block_detect -k <name.kmer_loci.tsv.gz> -s <samples.tsv> -O <outdir> -n <name>
 *                     [-t threads] [-m min_hits] [-V] [-c <name.context_loci.tsv.gz>] [-I max_insert]
 *
 *   samples.tsv: sample_name<TAB>R1.fastq.gz[<TAB>R2.fastq.gz]   (# = comment; same name = same sample)
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
    uint8_t  ctx;          /* 1 = context locus (-c): validates reads, is never counted or reported */
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
    uint8_t  ctx;          /* context k-mer (single copy, not rare): validation only */
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
    char *name;
    char **r1, **r2;       /* one or more inputs per sample (e.g. merged reads + unmerged pairs) */
    int   nf;
} Sample;

typedef struct {
    uint32_t *block_contra;         /* fragments that touched the block with rare k-mers but were refused */
    uint32_t *kmer_hits;     /* [n_kmers] */
    uint32_t *kmer_any;      /* [n_kmers] fragments holding the k-mer at all, accepted or not (non-aware count) */
    uint32_t *block_frags;   /* [max_block + 1] */
    uint32_t *block_refused; /* [max_block + 1] fragments that hit the block's k-mers but contradict the strain */
    uint32_t *block_starts;  /* [max_block + 1] */
    uint32_t *block_best;    /* [max_block + 1] */
    uint64_t *len_hist;      /* [MAX_READ_LEN + 1], N-free pieces >= k */
    uint32_t *fam_amb;       /* [max_family + 1] */
    uint32_t *fam_best;      /* [max_family + 1] ambiguous READS (pieces), one count each */
    uint64_t reads, fragments, kmer_positions, frags_with_hits, frags_placed,
             groups_accepted, groups_rejected, pieces_ambiguous, pieces_placed;
    /* why groups were rejected, read / fragment level counts, mate pairs */
    uint64_t groups_incomplete, groups_short, groups_contradict, reads_hit, reads_contradicted,
             frags_contradicted, frags_vetoed, pairs_placed, pairs_concordant, pairs_discordant,
             pairs_split_contig, groups_tolerated, groups_rescued, pairs_overlap, frags_hit;
    double   reads_contra_expected;   /* sum over reads_hit of P(at least one base error), from the qualities */
    int ok;
} SampleRes;

typedef struct {
    Hit    *hits;   size_t nh, hcap;
    Group  *grp;    size_t ng, gcap;
    uint32_t *fk;   size_t nfk, fkcap;   /* fragment k-mers */
    uint32_t *fb;   size_t nfb, fbcap;   /* fragment blocks */
    uint32_t *ff;   size_t nff, ffcap;   /* fragment ambiguous families */
    struct { uint32_t b, cnt; int64_t mn, mx; } *bs; size_t nbs, bscap;   /* per-block stats of a read */
    /* per-read counts are held back until the fragment is complete, so a contradicted fragment (-V)
       can be dropped as a whole */
    uint32_t *dst;  size_t ndst, dstcap; /* block starts */
    uint32_t *dbe;  size_t ndbe, dbecap; /* block best */
    uint32_t *dfb;  size_t ndfb, dfbcap; /* family best (ambiguous reads) */
    uint32_t *fa;   size_t nfa, facap;   /* rare k-mers seen anywhere in the fragment */
    uint32_t *cb;   size_t ncb, cbcap;   /* blocks whose k-mers sit in a contradicted group of this fragment */
    int      nrare;                      /* rare (non-context) hits in the current piece */
    const char *pq;                      /* qualities of the current piece (NULL = none) */
    /* -P */
    struct { uint32_t seg, o0, o1; int64_t a, b; uint8_t mate, rare, single; } *pend; size_t npend, pendcap;
    uint32_t *pobs; size_t npobs, pobscap;          /* loci seen by the pending (incomplete) groups */
    uint32_t *acc[2]; size_t nacc[2], acccap[2];    /* loci of the accepted groups, per mate */
    uint8_t *dstm; int64_t *dstlo; size_t dstmcap;  /* mate and leftmost coordinate of each 'starts' entry */
    uint8_t *dbem; size_t dbemcap;                  /* mate of each 'best' entry */
    int      mcontra[2];                            /* mate carries an unresolved contradiction */
    int      whole[2], mlen[2];                     /* mate was one N-free piece >= k, and its length */
    int      overlap;                               /* mates overlap by >= k bases (sequence overlap) */
    int      frag_rare;                             /* a mate of the fragment holds a rare k-mer */
    uint64_t *rcbuf; size_t rccap;                  /* 2-bit packed mates for the overlap search */
    uint8_t *seen; size_t seencap;       /* -Q: which expected k-mers of a group were observed */
    int      contra;                     /* fragment has a contradicted group */
    int      mplaced[2];                 /* mate placed uniquely as a whole read */
    struct { uint32_t contig; int8_t strand; int64_t lo, hi; } mp[2];
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
static int       frag_veto = 0;                  /* -V: drop fragments with a contradicted group */
static int64_t   max_insert = 1000;              /* -I: largest outer distance of a concordant pair */
static uint32_t  n_ctx_kmers = 0;
static int       min_qual = 0;                   /* -q: bases below this phred quality split the read like N */
static int       pair_mode = 0;                  /* -P: mates are one fragment: overlap bookkeeping, mate
                                                    rescue, discordant pairs contradict */
static int       tol_qual = 0;                   /* -Q: a broken run next to a base below this quality is
                                                    put down to a sequencing error, not to the strain */
static double    PERR[256];                      /* error probability of a phred+33 quality character */

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
    for (int q = 0; q < 256; q++) {
        int ph = q - 33;
        if (ph < 0) ph = 0;
        PERR[q] = 1.0;
        for (int t = 0; t < ph; t++) PERR[q] *= 0.7943282347242815;   /* 10^(-1/10) */
        if (PERR[q] > 0.75) PERR[q] = 0.75;
    }
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

static uint32_t cap_loci = 0;

/* read a loci table. context = 1: validation-only loci (#kmer contig_id kmer_position ref_rc) */
static void read_loci_table(const char *file, int context) {
    gzFile fp = gzopen(file, "r");
    if (!fp) { fprintf(stderr, "cannot open %s: %s\n", file, strerror(errno)); exit(1); }

    char *line = (char *)xrealloc(NULL, MAX_LINE);
    if (!gzgets(fp, line, MAX_LINE)) die("k-mer loci file is empty");

    /* header -> column indices */
    const char *want[] = { "#kmer", "contig_id", "kmer_position", "ref_rc",
                           "segment_id", "block_id", "family_id" };
    int col[7] = { -1, -1, -1, -1, -1, -1, -1 };
    const int n_need = context ? 4 : 7;
    {
        int c = 0;
        char *save = NULL;
        for (char *tok = strtok_r(line, "\t\r\n", &save); tok; tok = strtok_r(NULL, "\t\r\n", &save), c++)
            for (int w = 0; w < 7; w++) if (strcmp(tok, want[w]) == 0) col[w] = c;
    }
    for (int w = 0; w < n_need; w++) if (col[w] < 0) {
        fprintf(stderr, "error: column '%s' missing in %s\n"
                        "       (generate it with kmer_scrub_filter.py --reference)\n", want[w], file);
        exit(1);
    }

    char *f[64];
    while (gzgets(fp, line, MAX_LINE)) {
        if (line[0] == '\n' || line[0] == '\r' || line[0] == '\0') continue;
        int nf = 0;
        char *save = NULL;
        for (char *tok = strtok_r(line, "\t\r\n", &save); tok && nf < 64; tok = strtok_r(NULL, "\t\r\n", &save))
            f[nf++] = tok;
        int need = 0;
        for (int w = 0; w < n_need; w++) if (col[w] > need) need = col[w];
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
            kmers[n_kmers].ctx = (uint8_t)context;
            if (context) n_ctx_kmers++;
            ki = (int32_t)n_kmers++;
        } else if (context || kmers[ki].ctx) {
            fprintf(stderr, "error: context k-mer %s occurs more than once (context k-mers must be single "
                            "copy and must not be in the rare k-mer table)\n", ks);
            exit(1);
        }

        int ref_rc = parse_bool(f[col[3]]);
        uint64_t ref_fwd = ref_rc ? rc_w : w_code;
        uint64_t ref_rev = ref_rc ? w_code : rc_w;

        GROW(loci, n_loci, cap_loci, Locus);
        Locus *L = &loci[n_loci++];
        L->kmer      = (uint32_t)ki;
        L->contig    = contig_index(f[col[1]]);
        L->pos       = strtoll(f[col[2]], NULL, 10);
        L->seg       = context ? 0 : (uint32_t)strtoul(f[col[4]], NULL, 10);
        L->block     = context ? 0 : (uint32_t)strtoul(f[col[5]], NULL, 10);
        L->family    = context ? 0 : (uint32_t)strtoul(f[col[6]], NULL, 10);
        L->ctx       = (uint8_t)context;
        if (!context && L->block == 0) die("block_id 0 is reserved (context loci)");
        L->fwd_canon = (ref_fwd == ref_rev) ? 2 : (ref_fwd < ref_rev ? 1 : 0);
        kmers[ki].loc_n++;

        if (L->seg > max_seg) max_seg = L->seg;
        if (L->block > max_block) max_block = L->block;
        if (L->family > max_family) max_family = L->family;
    }
    gzclose(fp);
    free(line);
}

static int cmp_locus_contig_pos(const void *a, const void *b) {
    const Locus *x = (const Locus *)a, *y = (const Locus *)b;
    if (x->contig != y->contig) return x->contig < y->contig ? -1 : 1;
    if (x->pos != y->pos) return x->pos < y->pos ? -1 : 1;
    return 0;
}

static void load_loci(const char *file, const char *ctx_file) {
    read_loci_table(file, 0);
    uint32_t n_rare_loci = n_loci;
    if (n_loci == 0) {
        /* strain without detectable blocks (kmer_scrub_filter.py wrote an empty loci file): still scan the
           reads so detect_stats / read_lengths are complete, and write header-only hit tables */
        fprintf(stderr, "WARNING: no loci in %s: writing empty hit tables (strain not detectable)\n", file);
        if (K == 0) { K = 31; KMASK = (1ULL << (2 * K)) - 1; }
        hrehash(1u << 4);
    }

    /* repeat families = families with more than one segment of the rare table (before context is added) */
    family_repeat = (uint8_t *)calloc(max_family + 1, 1);
    uint32_t *fam_nseg = (uint32_t *)calloc(max_family + 1, sizeof(uint32_t));
    {
        uint8_t *seen = (uint8_t *)calloc(max_seg + 1, 1);
        for (uint32_t i = 0; i < n_loci; i++)
            if (!seen[loci[i].seg]) { seen[loci[i].seg] = 1; fam_nseg[loci[i].family]++; }
        free(seen);
    }
    for (uint32_t f2 = 0; f2 <= max_family; f2++) family_repeat[f2] = fam_nseg[f2] > 1;
    free(fam_nseg);

    /* context: validation segments become the runs of consecutive positions over rare + context loci, so
       a read is checked over everything it covers there, not only inside one rare chunk */
    if (ctx_file && n_rare_loci) {
        read_loci_table(ctx_file, 1);
        qsort(loci, n_loci, sizeof(Locus), cmp_locus_contig_pos);
        uint32_t sid = 0;
        for (uint32_t i = 0; i < n_loci; i++) {
            if (i && loci[i].contig == loci[i - 1].contig && loci[i].pos == loci[i - 1].pos)
                die("context and rare loci overlap (same contig and position)");
            if (i == 0 || loci[i].contig != loci[i - 1].contig || loci[i].pos != loci[i - 1].pos + 1) sid++;
            loci[i].seg = sid;
        }
        max_seg = sid;
    }

    /* sort loci by segment, position and build segment geometry */
    qsort(loci, n_loci, sizeof(Locus), cmp_locus_seg_pos);
    segs = (Segment *)calloc(max_seg + 1, sizeof(Segment));
    block_used = (uint8_t *)calloc(max_block + 1, 1);

    block_nk = (uint32_t *)calloc(max_block + 1, sizeof(uint32_t));
    for (uint32_t i = 0; i < n_loci; i++) {
        Locus *L = &loci[i];
        Segment *S = &segs[L->seg];
        if (!L->ctx) { block_used[L->block] = 1; block_nk[L->block]++; }
        if (!S->used) {
            S->used = 1; S->contig = L->contig; S->family = L->family;
            S->start = S->end = L->pos; S->first = i;
        } else {
            if (L->contig != S->contig || L->pos != S->end + 1) {
                fprintf(stderr, "error: segment %u is not contiguous (contig %s pos %" PRId64 ")\n",
                        L->seg, contigs[L->contig], L->pos);
                exit(1);
            }
            S->end = L->pos;
        }
    }

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
        if (!loci[i].ctx && kmers[loci[i].kmer].loc_n > 1) block_multi[loci[i].block] = 1;

    uint32_t n_segs = 0, n_blocks = 0, n_rep = 0, multi = 0;
    for (uint32_t s = 0; s <= max_seg; s++) n_segs += segs[s].used;
    for (uint32_t b = 0; b <= max_block; b++) n_blocks += block_used[b];
    for (uint32_t f2 = 0; f2 <= max_family; f2++) n_rep += family_repeat[f2];
    for (uint32_t i = 0; i < n_kmers; i++) multi += kmers[i].loc_n > 1;
    fprintf(stderr, "loaded %u rare k-mers (k=%d), %u context k-mers, %u loci, %u multi-copy k-mers, "
                    "%u validation segments, %u blocks, %u repeat families, %u contigs\n",
            n_kmers - n_ctx_kmers, K, n_ctx_kmers, n_loci, multi, n_segs, n_blocks, n_rep, n_contigs);
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

/* -Q: is every gap in the group's run of k-mers explained by a low-quality base? One wrong base at read
   position b removes the k-mers starting at b-K+1 .. b, so for a gap of k-mer starts m0 .. m1 the wrong base
   lies in m1 .. m0+K-1. The gap is explained if a base in that window is below tol_qual. A gap longer than
   K needs two or more differences close together and is never excused. A variant of a relative is usually
   read with good quality and stays a contradiction. */
static int gaps_explained(Work *w, size_t i, size_t j, int64_t a, int64_t b, int plen) {
    const int8_t strand = w->hits[i].strand;
    const int64_t diag = w->hits[i].diag;
    int64_t jlo = strand > 0 ? a - diag : diag - b;       /* read coordinate of the first expected k-mer */
    int64_t n = b - a + 1;
    if (n <= 0) return 0;
    if ((size_t)n > w->seencap) { w->seencap = (size_t)n * 2; w->seen = (uint8_t *)xrealloc(w->seen, w->seencap); }
    memset(w->seen, 0, (size_t)n);
    for (size_t h = i; h < j; h++) {
        int64_t p = loci[w->hits[h].locus].pos;
        int64_t jr = (strand > 0 ? p - diag : diag - p) - jlo;
        if (jr >= 0 && jr < n) w->seen[jr] = 1;
    }
    int64_t m = 0;
    while (m < n) {
        if (w->seen[m]) { m++; continue; }
        int64_t m0 = m;
        while (m < n && !w->seen[m]) m++;
        int64_t m1 = m - 1;
        if (m1 - m0 + 1 > K) return 0;       /* longer than one error can make: two or more differences */
        int64_t lo = m1 + jlo, hi = m0 + K - 1 + jlo;
        if (lo < 0) lo = 0;
        if (hi > plen - 1) hi = plen - 1;
        int low = 0;
        for (int64_t q = lo; q <= hi && !low; q++) low = (int)(unsigned char)w->pq[q] - 33 < tol_qual;
        if (!low) return 0;
    }
    return 1;
}

/* evaluate the hits of one N-free piece of a read. Counts that belong to the read (starts, best, family
   best) are buffered in the Work and only added by finish_fragment, so a contradicted fragment can be
   dropped as a whole (-V). Returns 1 if the piece carries a contradicted group. */
static int process_piece(Work *w, int plen, int whole_read, int mate, SampleRes *r) {
    /* a piece with context k-mers only says nothing about the strain on its own. In pair mode it is still
       checked and placed: if its mate carries rare k-mers, it validates the other end of the fragment */
    if (w->nh == 0 || (w->nrare == 0 && !(pair_mode && mate >= 0))) return 0;
    qsort(w->hits, w->nh, sizeof(Hit), cmp_hit);

    /* 1. groups by (segment, strand, diag) and completeness check */
    w->ng = 0;
    int contra = 0;
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

        /* rare = group holds a rare k-mer (groups of context k-mers only are not counted in the statistics);
           single = it holds a single-copy k-mer: such a k-mer has no other locus that could explain it, so an
           incomplete group is a contradiction (mismatch, indel, variant of a relative). An incomplete group
           of multi-copy k-mers only is usually the read seen at the wrong copy of a repeat. */
        int rare = 0, single = 0;
        for (size_t h = i; h < j; h++) {
            const Kmer *km = &kmers[w->hits[h].kmer];
            if (!km->ctx) rare = 1;
            if (km->loc_n == 1 && loci[w->hits[h].locus].fwd_canon != 2) single = 1;
        }

        /* -Q: an incomplete group whose gaps all sit next to a low-quality base keeps its observed k-mers and
           is no contradiction */
        int tolerated = 0;
        if (tol_qual && w->pq && n_obs < n_exp && n_obs >= need) tolerated = gaps_explained(w, i, j, a, b, plen);

        if ((n_obs == n_exp || tolerated) && n_obs >= need) {
            if (w->ng >= w->gcap) { w->gcap = w->gcap ? w->gcap * 2 : 64;
                                    w->grp = (Group *)xrealloc(w->grp, w->gcap * sizeof(Group)); }
            Group *g = &w->grp[w->ng++];
            g->contig = S->contig; g->seg = seg;
            g->strand = w->hits[i].strand; g->diag = w->hits[i].diag;
            g->h0 = (uint32_t)i; g->h1 = (uint32_t)j; g->nobs = (uint32_t)n_obs;
            if (rare) { if (tolerated) r->groups_tolerated++; else r->groups_accepted++; }
        } else if (pair_mode && mate >= 0 && n_obs < n_exp && n_obs >= need) {
            /* -P: decided in finish_fragment, where the other mate may show the missing k-mers */
            if (w->npend >= w->pendcap) { w->pendcap = w->pendcap ? w->pendcap * 2 : 16;
                                          w->pend = xrealloc(w->pend, w->pendcap * sizeof(*w->pend)); }
            w->pend[w->npend].seg = seg; w->pend[w->npend].a = a; w->pend[w->npend].b = b;
            w->pend[w->npend].mate = (uint8_t)mate; w->pend[w->npend].rare = (uint8_t)rare;
            w->pend[w->npend].single = (uint8_t)single;
            w->pend[w->npend].o0 = (uint32_t)w->npobs;
            for (size_t h = i; h < j; h++) push_u32(&w->pobs, &w->npobs, &w->pobscap, w->hits[h].locus);
            w->pend[w->npend].o1 = (uint32_t)w->npobs;
            w->npend++;
        } else {
            if (rare) {
                r->groups_rejected++;
                if (n_obs < n_exp) r->groups_incomplete++; else r->groups_short++;
            }
            if (n_obs < n_exp && single) {
                contra = 1;
                if (rare) r->groups_contradict++;
                for (size_t h = i; h < j; h++)
                    if (!loci[w->hits[h].locus].ctx) push_u32(&w->cb, &w->ncb, &w->cbcap, loci[w->hits[h].locus].block);
            }
        }
        if (pair_mode && mate >= 0 && w->ng && w->grp[w->ng - 1].h0 == (uint32_t)i && w->grp[w->ng - 1].seg == seg)
            for (size_t h = i; h < j; h++) push_u32(&w->acc[mate], &w->nacc[mate], &w->acccap[mate], w->hits[h].locus);
        i = j;
    }
    if (contra) {
        w->contra = 1;
        if (pair_mode && mate >= 0) w->mcontra[mate] = 1; else r->reads_contradicted++;
    }
    if (w->ng == 0) return contra;

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
    int counted = 0, amb_counted = 0;

    g0 = 0;
    while (g0 < w->ng) {
        size_t g1 = g0 + 1;
        uint32_t score = w->grp[g0].nobs;
        while (g1 < w->ng && same_place(&w->grp[g1], &w->grp[g0])) score += w->grp[g1++].nobs;
        if (score == best) {
            int64_t flo = 0, fhi = 0;
            int n_rare_hits = 0;
            footprint(w->grp[g0].strand, w->grp[g0].diag, plen, &flo, &fhi);
            for (size_t g = g0; g < g1; g++) {
                const Group *G = &w->grp[g];
                for (uint32_t h = G->h0; h < G->h1; h++) {
                    const Hit *H = &w->hits[h];
                    const Locus *Lc = &loci[H->locus];
                    if (Lc->ctx) continue;
                    n_rare_hits++;
                    push_u32(&w->fk, &w->nfk, &w->fkcap, H->kmer);
                    if (!ambiguous) push_u32(&w->fb, &w->nfb, &w->fbcap, Lc->block);
                    else {
                        push_u32(&w->ff, &w->nff, &w->ffcap, Lc->family);
                        if (!amb_counted) { push_u32(&w->dfb, &w->ndfb, &w->dfbcap, Lc->family); amb_counted = 1; }
                    }
                }
                if (!ambiguous && whole_read) {
                    const Segment *S = &segs[G->seg];
                    if (flo >= S->start && flo <= S->end) {
                        const Locus *Ls = &loci[S->first + (uint32_t)(flo - S->start)];
                        if (!Ls->ctx) {
                            if (w->ndst >= w->dstmcap) {
                                w->dstmcap = w->dstmcap ? w->dstmcap * 2 : 16;
                                w->dstm = (uint8_t *)xrealloc(w->dstm, w->dstmcap);
                                w->dstlo = (int64_t *)xrealloc(w->dstlo, w->dstmcap * sizeof(int64_t));
                            }
                            w->dstm[w->ndst] = (uint8_t)(mate > 0); w->dstlo[w->ndst] = flo;
                            push_u32(&w->dst, &w->ndst, &w->dstcap, Ls->block);
                        }
                    }
                }
            }
            if (n_rare_hits && !counted) {
                if (ambiguous) r->pieces_ambiguous++; else r->pieces_placed++;
                counted = 1;
            }
            if (!ambiguous && n_rare_hits) {
                /* one block per read: single-copy blocks first (a read touching a unique block counts there),
                   then the most completely covered block, ties -> closest to the read's 5' end */
                w->nbs = 0;
                for (size_t g = g0; g < g1; g++) {
                    const Group *G = &w->grp[g];
                    for (uint32_t h = G->h0; h < G->h1; h++) {
                        const Locus *Lc = &loci[w->hits[h].locus];
                        if (Lc->ctx) continue;
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
                if (w->nbs) {
                    if (w->ndbe >= w->dbemcap) { w->dbemcap = w->dbemcap ? w->dbemcap * 2 : 16;
                                                 w->dbem = (uint8_t *)xrealloc(w->dbem, w->dbemcap); }
                    w->dbem[w->ndbe] = (uint8_t)(mate > 0);
                    push_u32(&w->dbe, &w->ndbe, &w->dbecap, w->bs[bi].b);
                }
                if (whole_read && mate >= 0 && mate < 2) {
                    w->mplaced[mate] = 1;
                    w->mp[mate].contig = w->grp[g0].contig; w->mp[mate].strand = strand;
                    w->mp[mate].lo = flo; w->mp[mate].hi = fhi + K - 1;
                }
            } else if (pair_mode && !ambiguous && whole_read && mate >= 0 && mate < 2) {
                w->mplaced[mate] = 1;                 /* context-only mate: placed, never counted */
                w->mp[mate].contig = w->grp[g0].contig; w->mp[mate].strand = w->grp[g0].strand;
                w->mp[mate].lo = flo; w->mp[mate].hi = fhi + K - 1;
            }
        }
        g0 = g1;
    }
    return contra;
}

static void scan_read(Work *w, const char *s, const char *qual, int L, int mate, SampleRes *r) {
    r->reads++;
    if (mate >= 0 && mate < 2) { w->whole[mate] = 0; w->mlen[mate] = L; }
    if (L < K) return;
    uint64_t fwd = 0, rcv = 0;
    int l = 0, ps = 0, read_rare = 0;
    const int shift = 2 * (K - 1);
    w->nh = 0; w->nrare = 0;

    for (int i = 0; i < L; i++) {
        int c = NT[(unsigned char)s[i]];
        /* -q: a low-quality base is treated like N. The read is cut there and the pieces are checked on
           their own, so an unreliable base can neither support nor contradict the strain. */
        if (min_qual && qual && (int)(unsigned char)qual[i] - 33 < min_qual) c = 4;
        if (c > 3) {
            if (i - ps >= K) r->len_hist[i - ps > MAX_READ_LEN ? MAX_READ_LEN : i - ps]++;
            read_rare |= w->nrare > 0;
            w->pq = qual ? qual + ps : NULL;
            process_piece(w, i - ps, 0, mate, r);
            w->nh = 0; w->nrare = 0; l = 0; ps = i + 1;
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
        if (!km->ctx) { w->nrare++; push_u32(&w->fa, &w->nfa, &w->facap, (uint32_t)ki); }
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
    read_rare |= w->nrare > 0;
    w->pq = qual ? qual + ps : NULL;
    if (mate >= 0 && mate < 2) { w->whole[mate] = (ps == 0 && L >= K); w->mlen[mate] = L; }
    process_piece(w, L - ps, ps == 0, mate, r);
    if (read_rare) {
        w->frag_rare = 1;
        r->reads_hit++;
        if (qual) {
            /* chance that this read carries at least one sequencing error: the contradiction rate a true
               strain would show if the read were validated over its whole length (as with -c). It makes
               reads_contradicted / reads_hit comparable between libraries of different quality and length. */
            double ok = 1.0;
            for (int i = 0; i < L; i++)
                if (!min_qual || (int)(unsigned char)qual[i] - 33 >= min_qual) ok *= 1.0 - PERR[(unsigned char)qual[i]];
            r->reads_contra_expected += 1.0 - ok;
        }
    }
    w->nh = 0; w->nrare = 0;
}

static int cmp_u32(const void *a, const void *b) {
    uint32_t x = *(const uint32_t *)a, y = *(const uint32_t *)b;
    return (x > y) - (x < y);
}

/* -P: do the two mates overlap? R2 is reverse-complemented and slid along R1. offset o >= 0: R2rc starts at
   base o of R1; o < 0: it starts before R1 (insert shorter than the reads). The first offset whose overlap is
   at least K bases with at most 5 mismatches and at most 20 % mismatches wins (fastp's merge limits),
   largest overlaps first.
   Returns the overlap length (0 = none); *uni = bases covered by the two mates together. */
static inline uint64_t get32(const uint64_t *a, int pos) {
    int w0 = pos >> 5, sh = (pos & 31) << 1;
    return sh ? (a[w0] >> sh) | (a[w0 + 1] << (64 - sh)) : a[w0];
}

static int pair_overlap(Work *w, const char *s1, int L1, const char *s2, int L2, int *uni) {
    /* both reads packed 2 bits per base (callers pass N-free reads), R2 reverse-complemented; 32 bases are
       compared per step, so a wrong offset is dropped after one or two steps */
    size_t n1 = (size_t)(L1 >> 5) + 2, n2 = (size_t)(L2 >> 5) + 2;
    if (n1 + n2 > w->rccap) { w->rccap = (n1 + n2) * 2; w->rcbuf = (uint64_t *)xrealloc(w->rcbuf, w->rccap * sizeof(uint64_t)); }
    uint64_t *p1 = w->rcbuf, *p2 = w->rcbuf + n1;
    memset(p1, 0, (n1 + n2) * sizeof(uint64_t));
    for (int i = 0; i < L1; i++) p1[i >> 5] |= (uint64_t)(NT[(unsigned char)s1[i]] & 3) << ((i & 31) << 1);
    for (int i = 0; i < L2; i++)
        p2[i >> 5] |= (uint64_t)(3 - (NT[(unsigned char)s2[L2 - 1 - i]] & 3)) << ((i & 31) << 1);
    const uint64_t odd = 0x5555555555555555ULL;
    for (int pass = 0; pass < 2; pass++) {
        int o_lo = pass == 0 ? 0 : 1, o_hi = pass == 0 ? L1 - K : L2 - K;
        for (int oo = o_lo; oo <= o_hi; oo++) {
            int o = pass == 0 ? oo : -oo;
            int a1 = o > 0 ? o : 0, a2 = o > 0 ? 0 : -o;
            int ov = L1 - a1 < L2 - a2 ? L1 - a1 : L2 - a2;
            if (ov < K) break;
            int lim = ov / 5 < 5 ? ov / 5 : 5, mm = 0;
            for (int q = 0; q < ov && mm <= lim; q += 32) {
                uint64_t x = get32(p1, a1 + q) ^ get32(p2, a2 + q);
                x = (x | (x >> 1)) & odd;
                int n = ov - q;
                if (n < 32) x &= (1ULL << (n << 1)) - 1;
                mm += __builtin_popcountll(x);
            }
            if (mm <= lim) { *uni = L1 + L2 - ov; return ov; }
        }
    }
    *uni = L1 + L2;
    return 0;
}

static int has_u32(const uint32_t *a, size_t n, uint32_t v) {
    size_t lo = 0, hi = n;
    while (lo < hi) { size_t m = (lo + hi) / 2; if (a[m] < v) lo = m + 1; else hi = m; }
    return lo < n && a[lo] == v;
}

/* -P: settle the incomplete groups of a pair. A group is RESCUED if the other mate holds an accepted group
   with every k-mer this one misses: the two mates cover the same piece of the reference, one of them reads
   it without a difference, so the difference is a sequencing error. Its observed k-mers then count. A group
   that is not rescued is rejected, and is a contradiction if it holds a single-copy k-mer. */
static void settle_pending(Work *w, SampleRes *r, int paired) {
    if (w->npend == 0) return;
    if (paired) for (int m = 0; m < 2; m++) if (w->nacc[m] > 1) qsort(w->acc[m], w->nacc[m], sizeof(uint32_t), cmp_u32);
    for (size_t p = 0; p < w->npend; p++) {
        const Segment *S = &segs[w->pend[p].seg];
        const int mate = w->pend[p].mate, other = 1 - mate;
        int rescued = paired && w->nacc[other] > 0;
        for (int64_t pos = w->pend[p].a; rescued && pos <= w->pend[p].b; pos++) {
            uint32_t li = S->first + (uint32_t)(pos - S->start);
            int seen = 0;
            for (uint32_t q = w->pend[p].o0; q < w->pend[p].o1 && !seen; q++) seen = w->pobs[q] == li;
            if (!seen && !has_u32(w->acc[other], w->nacc[other], li)) rescued = 0;
        }
        if (rescued) {
            if (w->pend[p].rare) r->groups_rescued++;
            for (uint32_t q = w->pend[p].o0; q < w->pend[p].o1; q++) {
                const Locus *Lc = &loci[w->pobs[q]];
                if (Lc->ctx) continue;
                push_u32(&w->fk, &w->nfk, &w->fkcap, Lc->kmer);
                push_u32(&w->fb, &w->nfb, &w->fbcap, Lc->block);
            }
        } else {
            if (w->pend[p].rare) { r->groups_rejected++; r->groups_incomplete++; }
            if (w->pend[p].single) {
                w->contra = 1; w->mcontra[mate] = 1;
                if (w->pend[p].rare) r->groups_contradict++;
                for (uint32_t q = w->pend[p].o0; q < w->pend[p].o1; q++)
                    if (!loci[w->pobs[q]].ctx) push_u32(&w->cb, &w->ncb, &w->cbcap, loci[w->pobs[q]].block);
            }
        }
    }
}

/* add the de-duplicated fragment counts and reset the fragment buffers. With -V a fragment that carries a
   contradicted group anywhere (either mate) is dropped as a whole. paired = the fragment has two mates. */
static void finish_fragment(Work *w, SampleRes *r, int paired) {
    r->fragments++;
    if (w->frag_rare) r->frags_hit++;
    if (w->nfa) {
        qsort(w->fa, w->nfa, sizeof(uint32_t), cmp_u32);
        for (size_t i = 0; i < w->nfa; i++) if (i == 0 || w->fa[i] != w->fa[i - 1]) r->kmer_any[w->fa[i]]++;
        w->nfa = 0;
    }
    if (pair_mode) {
        if (!w->frag_rare) {                 /* context k-mers only: not a fragment of interest */
            w->npend = 0; w->contra = 0; w->mcontra[0] = w->mcontra[1] = 0; w->mplaced[0] = w->mplaced[1] = 0;
            w->ncb = 0;
        }
        settle_pending(w, r, paired);
        r->reads_contradicted += (uint64_t)(w->mcontra[0] + w->mcontra[1]);
    }
    if (w->mplaced[0] && w->mplaced[1]) {
        r->pairs_placed++;
        if (w->mp[0].contig != w->mp[1].contig) r->pairs_split_contig++;
        else {
            int64_t lo = w->mp[0].lo < w->mp[1].lo ? w->mp[0].lo : w->mp[1].lo;
            int64_t hi = w->mp[0].hi > w->mp[1].hi ? w->mp[0].hi : w->mp[1].hi;
            int ok = w->mp[0].strand != w->mp[1].strand && hi - lo + 1 <= max_insert;
            if (ok && pair_mode) {
                /* forward mate on the left, reverse mate on the right (3 bp slack for trimmed ends) */
                int f = w->mp[0].strand > 0 ? 0 : 1, v = 1 - f;
                ok = w->mp[f].lo <= w->mp[v].lo + 3 && w->mp[f].hi <= w->mp[v].hi + 3;
            }
            if (ok) r->pairs_concordant++;
            else { r->pairs_discordant++; if (pair_mode) w->contra = 1; }
        }
    }
    if (w->contra) r->frags_contradicted++;
    if (w->contra) {
        /* blocks this fragment touched but does not support: those in its contradicted groups, and with -V
           every block it would otherwise have counted for */
        if (frag_veto) for (size_t i = 0; i < w->nfb; i++) push_u32(&w->cb, &w->ncb, &w->cbcap, w->fb[i]);
        if (w->ncb) {
            qsort(w->cb, w->ncb, sizeof(uint32_t), cmp_u32);
            for (size_t i = 0; i < w->ncb; i++) if (i == 0 || w->cb[i] != w->cb[i - 1]) r->block_refused[w->cb[i]]++;
        }
    }
    w->ncb = 0;
    if (frag_veto && w->contra) {
        if (w->nfk) r->frags_vetoed++;
        w->nfk = w->nfb = w->nff = 0;
        w->ndst = w->ndbe = w->ndfb = 0;
    }
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
    if (pair_mode && w->overlap) {
        /* overlapping mates are one read of the fragment's length: one 'best' block and one start, as a
           merged read would give. best: the first mate that has one. start: the mate that reaches furthest
           left (the first mate if only one mate is placed). */
        if (w->ndbe) r->block_best[w->dbe[0]]++;
        int left = (w->mplaced[0] && w->mplaced[1]) ? (w->mp[1].lo < w->mp[0].lo) : (w->mplaced[1] && !w->mplaced[0]);
        for (size_t i = 0; i < w->ndst; i++) if (w->dstm[i] == left) { r->block_starts[w->dst[i]]++; break; }
    } else {
        for (size_t i = 0; i < w->ndst; i++) r->block_starts[w->dst[i]]++;
        for (size_t i = 0; i < w->ndbe; i++) r->block_best[w->dbe[i]]++;
    }
    for (size_t i = 0; i < w->ndfb; i++) r->fam_best[w->dfb[i]]++;
    w->nfk = w->nfb = w->nff = 0;
    w->ndst = w->ndbe = w->ndfb = 0;
    w->contra = 0; w->mplaced[0] = w->mplaced[1] = 0;
    w->npend = w->npobs = 0; w->nacc[0] = w->nacc[1] = 0;
    w->mcontra[0] = w->mcontra[1] = 0; w->overlap = 0; w->frag_rare = 0;
}

/* ── per-sample job ─────────────────────────────────────────────────────── */

static int run_files(const Sample *sm, const char *r1, const char *r2, Work *w, SampleRes *r) {
    /* testing aid: KBD_OVERLAP_DUMP=<prefix> writes read name, overlap and fragment length of every pair
       booked as overlapping to <prefix>.<sample>.tsv */
    FILE *odump = NULL;
    if (pair_mode && r2 && getenv("KBD_OVERLAP_DUMP")) {
        char path[4096];
        snprintf(path, sizeof(path), "%s.%s.tsv", getenv("KBD_OVERLAP_DUMP"), sm->name);
        odump = fopen(path, "a");
    }
    gzFile f1 = gzopen(r1, "r");
    if (!f1) { fprintf(stderr, "[%s] cannot open %s: %s\n", sm->name, r1, strerror(errno)); return 0; }
    gzFile f2 = NULL;
    if (r2) {
        f2 = gzopen(r2, "r");
        if (!f2) { fprintf(stderr, "[%s] cannot open %s: %s\n", sm->name, r2, strerror(errno)); gzclose(f1); return 0; }
    }
    gzbuffer(f1, 1 << 17);
    if (f2) gzbuffer(f2, 1 << 17);
    kseq_t *s1 = kseq_init(f1);
    kseq_t *s2 = f2 ? kseq_init(f2) : NULL;

    for (;;) {
        int l1 = kseq_read(s1);
        int l2 = s2 ? kseq_read(s2) : -1;
        if (l1 < 0) {
            if (s2 && l2 >= 0) fprintf(stderr, "[%s] WARNING: R2 has more reads than R1\n", sm->name);
            break;
        }
        if (s2 && l2 < 0) { fprintf(stderr, "[%s] WARNING: R1 has more reads than R2\n", sm->name); break; }
        scan_read(w, s1->seq.s, s1->qual.l == s1->seq.l ? s1->qual.s : NULL, l1, s2 ? 0 : -1, r);
        if (s2) scan_read(w, s2->seq.s, s2->qual.l == s2->seq.l ? s2->qual.s : NULL, l2, 1, r);
        if (pair_mode && s2 && w->whole[0] && w->whole[1]) {
            /* overlapping mates are booked as ONE read of the fragment's length, so k-mer positions, read
               lengths (the exposure model) and total bases count every sequenced base once */
            int uni = 0, ov = pair_overlap(w, s1->seq.s, l1, s2->seq.s, l2, &uni);
            if (ov >= K) {
                w->overlap = 1;
                r->pairs_overlap++;
                r->kmer_positions -= (uint64_t)(ov - K + 1);
                r->reads--;
                r->len_hist[l1 > MAX_READ_LEN ? MAX_READ_LEN : l1]--;
                r->len_hist[l2 > MAX_READ_LEN ? MAX_READ_LEN : l2]--;
                r->len_hist[uni > MAX_READ_LEN ? MAX_READ_LEN : uni]++;
                if (odump) fprintf(odump, "%s\t%d\t%d\n", s1->name.s, ov, uni);
            }
        }
        finish_fragment(w, r, s2 != NULL);
    }

    kseq_destroy(s1); gzclose(f1);
    if (s2) { kseq_destroy(s2); gzclose(f2); }
    if (odump) fclose(odump);
    return 1;
}

static void run_sample(const Sample *sm, SampleRes *r) {
    r->kmer_hits    = (uint32_t *)calloc(n_kmers, sizeof(uint32_t));
    r->kmer_any     = (uint32_t *)calloc(n_kmers, sizeof(uint32_t));
    r->block_frags  = (uint32_t *)calloc(max_block + 1, sizeof(uint32_t));
    r->block_refused = (uint32_t *)calloc(max_block + 1, sizeof(uint32_t));
    r->block_starts = (uint32_t *)calloc(max_block + 1, sizeof(uint32_t));
    r->block_best   = (uint32_t *)calloc(max_block + 1, sizeof(uint32_t));
    r->len_hist     = (uint64_t *)calloc(MAX_READ_LEN + 1, sizeof(uint64_t));
    r->fam_amb      = (uint32_t *)calloc(max_family + 1, sizeof(uint32_t));
    r->fam_best     = (uint32_t *)calloc(max_family + 1, sizeof(uint32_t));

    Work w; memset(&w, 0, sizeof(w));
    int ok = 1;
    for (int f = 0; f < sm->nf && ok; f++) ok = run_files(sm, sm->r1[f], sm->r2[f], &w, r);
    free(w.hits); free(w.grp); free(w.fk); free(w.fb); free(w.ff); free(w.bs);
    free(w.dst); free(w.dbe); free(w.dfb); free(w.seen); free(w.cb); free(w.fa);
    free(w.pend); free(w.pobs); free(w.acc[0]); free(w.acc[1]); free(w.dstm); free(w.dstlo); free(w.dbem); free(w.rcbuf);
    if (!ok) return;
    r->ok = 1;
    fprintf(stderr, "[%s] %" PRIu64 " fragments, %" PRIu64 " with accepted hits, "
                    "%" PRIu64 " groups accepted, %" PRIu64 " rejected, %" PRIu64 " ambiguous pieces, "
                    "%" PRIu64 " fragments contradicted%s\n",
            sm->name, r->fragments, r->frags_with_hits, r->groups_accepted, r->groups_rejected,
            r->pieces_ambiguous, r->frags_contradicted, frag_veto ? " (dropped)" : "");
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

/* rows with the same sample name are pooled into one sample (e.g. merged reads as one row and the
   unmerged R1 / R2 pairs as a second row) */
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
        size_t q = 0;
        while (q < n && strcmp(s[q].name, f[0]) != 0) q++;
        if (q == n) {
            GROW(s, n, cap, Sample);
            s[n].name = strdup(f[0]); s[n].r1 = NULL; s[n].r2 = NULL; s[n].nf = 0;
            n++;
        }
        s[q].r1 = (char **)xrealloc(s[q].r1, (size_t)(s[q].nf + 1) * sizeof(char *));
        s[q].r2 = (char **)xrealloc(s[q].r2, (size_t)(s[q].nf + 1) * sizeof(char *));
        s[q].r1[s[q].nf] = strdup(f[1]);
        s[q].r2[s[q].nf] = (nf > 2 && f[2][0]) ? strdup(f[2]) : NULL;
        s[q].nf++;
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
    /* k-mer hits, legacy layout (rare k-mers only: context k-mers are never reported) */
    gzFile g = open_out(dir, name, ".kmer_hits.tsv.gz");
    gzprintf(g, "#kmer");
    for (int s = 0; s < ns; s++) gzprintf(g, "\t%s", sm[s].name);
    gzprintf(g, "\n");
    for (uint32_t i = 0; i < n_kmers; i++) {
        if (kmers[i].ctx) continue;
        gzputs(g, kmers[i].str);
        for (int s = 0; s < ns; s++) gzprintf(g, "\t%u", res[s].kmer_hits[i]);
        gzputc(g, '\n');
    }
    gzprintf(g, "total_evaluated");
    for (int s = 0; s < ns; s++) gzprintf(g, "\t%" PRIu64, res[s].kmer_positions);
    gzprintf(g, "\n");
    gzclose(g);

    /* the same table for every fragment that holds the k-mer, accepted or not: what a detector that does not
       look at reads would count. Use it for sensitivity at low depth; the difference to kmer_hits is what
       the read checks removed. */
    g = open_out(dir, name, ".kmer_hits_any.tsv.gz");
    gzprintf(g, "#kmer");
    for (int s = 0; s < ns; s++) gzprintf(g, "\t%s", sm[s].name);
    gzprintf(g, "\n");
    for (uint32_t i = 0; i < n_kmers; i++) {
        if (kmers[i].ctx) continue;
        gzputs(g, kmers[i].str);
        for (int s = 0; s < ns; s++) gzprintf(g, "\t%u", res[s].kmer_any[i]);
        gzputc(g, '\n');
    }
    gzprintf(g, "total_evaluated");
    for (int s = 0; s < ns; s++) gzprintf(g, "\t%" PRIu64, res[s].kmer_positions);
    gzprintf(g, "\n");
    gzclose(g);

    /* blocks, long format, full grid */
    g = open_out(dir, name, ".block_hits.tsv.gz");
    gzprintf(g, "block_id\tsample\tfrags\tstarts\tbest\trefused\n");
    for (int s = 0; s < ns; s++)
        for (uint32_t b = 0; b <= max_block; b++) if (block_used[b])
            gzprintf(g, "%u\t%s\t%u\t%u\t%u\t%u\n", b, sm[s].name, res[s].block_frags[b], res[s].block_starts[b],
                     res[s].block_best[b], res[s].block_refused[b]);
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

    /* stats: the first ten columns are unchanged; the rest say why groups were rejected (incomplete =
       k-mers missing inside the read, short = complete but below -m), how many of the incomplete ones
       are contradictions (hold a single-copy k-mer), and what happened per read, fragment and mate pair.
       reads_contra_expected: contradicted reads expected from sequencing errors alone (base qualities of
       the reads with rare hits, whole read validated); 0 for FASTA input.
       frags_vetoed: fragments that had accepted hits and lost them to -V (with -c most contradictions sit in
       the same validation segment as the rare k-mers, so the group is rejected directly and never gets here).
       With -c the groups_* columns count the longer validation segments: compare runs on the reads_* columns.
       With -P: reads = overlapping mates counted as one read; groups_rescued, pairs_overlap; frags_hit =
       fragments with a rare k-mer, the denominator for frags_contradicted */
    char path[4096];
    snprintf(path, sizeof(path), "%s/%s.detect_stats.tsv", dir, name);
    FILE *fp = fopen(path, "w");
    if (!fp) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
    fprintf(fp, "sample\tfragments\treads\tkmer_positions\tfrags_with_hits\tfrags_placed\t"
                "groups_accepted\tgroups_rejected\tpieces_placed\tpieces_ambiguous\t"
                "groups_incomplete\tgroups_short\tgroups_contradict\treads_hit\treads_contradicted\t"
                "frags_contradicted\tfrags_vetoed\tpairs_placed\tpairs_concordant\tpairs_discordant\t"
                "pairs_split_contig\treads_contra_expected\tgroups_tolerated\tgroups_rescued\tpairs_overlap\tfrags_hit\tpair_mode\n");
    for (int s = 0; s < ns; s++) {
        SampleRes *r = &res[s];
        fprintf(fp, "%s\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64
                    "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64,
                sm[s].name, r->fragments, r->reads, r->kmer_positions, r->frags_with_hits, r->frags_placed,
                r->groups_accepted, r->groups_rejected, r->pieces_placed, r->pieces_ambiguous);
        fprintf(fp, "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64
                    "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%.2f\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%d\n",
                r->groups_incomplete, r->groups_short, r->groups_contradict, r->reads_hit,
                r->reads_contradicted, r->frags_contradicted, r->frags_vetoed, r->pairs_placed,
                r->pairs_concordant, r->pairs_discordant, r->pairs_split_contig, r->reads_contra_expected,
                r->groups_tolerated, r->groups_rescued, r->pairs_overlap, r->frags_hit, pair_mode);
    }
    fclose(fp);
}

/* ── main ───────────────────────────────────────────────────────────────── */

static void usage(void) {
    fprintf(stderr,
        "Usage: kmer_block_detect -k <name.kmer_loci.tsv.gz> -s <samples.tsv> -O <outdir> -n <name>\n"
        "                         [-t threads (4)] [-m min_hits (1)] [-V] [-c <name.context_loci.tsv.gz>]\n"
        "                         [-I max_insert (1000)] [-q min_base_quality (0 = off)]\n"
        "                         [-Q error_quality (0 = off)] [-P]\n\n"
        "  -k  k-mer loci table from kmer_scrub_filter.py --reference\n"
        "  -s  samples: name<TAB>R1.fastq[.gz][<TAB>R2.fastq[.gz]] per line (# = comment).\n"
        "      Rows with the same name are pooled into one sample, e.g. the merged reads on one\n"
        "      row and the unmerged R1 / R2 pairs on a second row.\n"
        "  -O  output directory\n"
        "  -n  output basename\n"
        "  -t  threads (one sample per thread)\n"
        "  -m  min k-mers a read must hit in a segment (capped at the segment size).\n"
        "      Reads cut by their own end that only clip a segment edge are\n"
        "      rejected below this. Default 1.\n"
        "  -V  fragment veto: a fragment (both mates) is dropped completely if any of its groups is a\n"
        "      contradiction, i.e. incomplete and holding a single-copy k-mer. Default off: only the\n"
        "      contradicted group is dropped and the rest of the read still counts.\n"
        "  -c  context loci from kmer_scrub_filter.py --context_bp: single-copy reference k-mers next to\n"
        "      the rare segments. They are never counted or reported, but a read has to match them too,\n"
        "      so it is checked over its whole length instead of only inside the rare segment.\n"
        "  -q  bases below this phred quality are treated like N: the read is cut there and each piece is\n"
        "      checked on its own. Replaces quality trimming in fastp (adapter removal and merging are still\n"
        "      needed). Reads cut this way do not count for 'starts' or the mate-pair statistics. Default 0.\n"
        "  -P  pair mode, for unmerged R1 / R2 input (fastp without --merge):\n"
        "      - mates that overlap by >= k bases are booked as one read of the fragment's length (k-mer\n"
        "        positions, read lengths, one 'best' and one 'start' per fragment);\n"
        "      - a broken run in one mate is rescued when the other mate reads the same reference k-mers\n"
        "        without a difference (groups_rescued): a sequencing error, not a variant;\n"
        "      - mates placed on one contig must face each other within -I, otherwise the fragment is a\n"
        "        contradiction. Use with -V and -c.\n"
        "  -Q  quality-aware rejection: if every gap in a read's run of k-mers lies next to a base below this\n"
        "      phred quality, the gap is put down to a sequencing error. The read keeps the k-mers it did hit\n"
        "      and is not a contradiction (groups_tolerated in detect_stats). A gap with only good bases\n"
        "      around it still rejects the group. Needs FASTQ qualities. Default 0 = every gap rejects.\n"
        "  -I  largest outer distance (bp) of two mates on one contig that still counts as a concordant\n"
        "      pair (pairs_concordant / pairs_discordant in detect_stats; diagnostic only). Default 1000.\n\n"
        "A read counts for a segment only if it hits every segment k-mer inside its\n"
        "footprint on one diagonal (no internal gaps). See the header of the source\n"
        "for the output columns.\n");
}

int main(int argc, char *argv[]) {
    char *kfile = NULL, *sfile = NULL, *outdir = NULL, *name = NULL, *cfile = NULL;
    int threads = 4, c;

    while ((c = getopt(argc, argv, "k:s:O:n:t:m:c:I:q:Q:VPh")) != -1) {
        switch (c) {
            case 'k': kfile = optarg; break;
            case 's': sfile = optarg; break;
            case 'O': outdir = optarg; break;
            case 'n': name = optarg; break;
            case 't': threads = atoi(optarg); break;
            case 'm': min_hits = atoi(optarg); break;
            case 'c': cfile = optarg; break;
            case 'I': max_insert = atoll(optarg); break;
            case 'q': min_qual = atoi(optarg); break;
            case 'Q': tol_qual = atoi(optarg); break;
            case 'V': frag_veto = 1; break;
            case 'P': pair_mode = 1; break;
            case 'h': usage(); return 0;
            default:  usage(); return 1;
        }
    }
    if (!kfile || !sfile || !outdir || !name) { usage(); return 1; }
    if (threads < 1) threads = 1;
    if (min_hits < 1) min_hits = 1;

    init_nt();
    load_loci(kfile, cfile);

    Sample *samples = NULL;
    int ns = load_samples(sfile, &samples);
    if (ns == 0) die("no samples");
    fprintf(stderr, "loaded %d sample(s); min_hits=%d; fragment veto %s; threads=%d\n", ns, min_hits,
            frag_veto ? "on" : "off", threads);

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
        free(res[s].kmer_hits); free(res[s].kmer_any); free(res[s].block_frags); free(res[s].block_refused); free(res[s].block_starts); free(res[s].fam_amb);
        free(res[s].block_best); free(res[s].len_hist); free(res[s].fam_best);
        for (int f = 0; f < samples[s].nf; f++) { free(samples[s].r1[f]); free(samples[s].r2[f]); }
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
