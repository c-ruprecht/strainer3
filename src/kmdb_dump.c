/*
	kmdb_dump — print the contents of a .kmdb k-mer database.

	Default: one line per k-mer, "<kmer>\t<count>", sorted ascending
	(canonical 2-bit order), preceded by a "#kmer\tcount" header line.

	  kmdb_dump x.kmdb                    all k-mers with counts
	  kmdb_dump -i x.kmdb                 header + metadata only
	  kmdb_dump -n 10 x.kmdb              first 10 k-mers
	  kmdb_dump -c 5 x.kmdb               k-mers with count >= 5
	  kmdb_dump -C 2 x.kmdb               k-mers with count <= 2
	  kmdb_dump -r x.kmdb                 packed uint64 instead of ACGT
	  kmdb_dump -K x.kmdb                 k-mers only (no count column)
	  kmdb_dump -q kmers.txt x.kmdb       look up k-mers (one per line, '-' =
	                                      stdin; any strand); prints
	                                      "<query>\t<count>", 0 if absent
	  kmdb_dump -s x.kmdb                 summary: n, count histogram stats

	Uses kmdb_read_block() only, so it reads exactly what every other tool
	reads. Block decoding is the only real work; output is fully buffered.
*/

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <inttypes.h>
#include <errno.h>
#include "kmdb.h"

#define OUT_BUF (8u << 20)

static void usage(int status)
{
	fprintf(stderr,
	"Usage: kmdb_dump [options] <file.kmdb>\n"
	"  -i          info only: header fields and metadata, no k-mers\n"
	"  -s          summary only: k-mer count and count statistics (reads all blocks)\n"
	"  -n <N>      print at most N k-mers (after filtering)\n"
	"  -c <MIN>    only k-mers with count >= MIN\n"
	"  -C <MAX>    only k-mers with count <= MAX\n"
	"  -r          print the packed uint64 instead of the ACGT string\n"
	"  -K          k-mers only, no count column\n"
	"  -H          no '#' header line\n"
	"  -q <file>   look up the k-mers in <file> (one per line, '-' = stdin,\n"
	"              either strand, must be length k). Prints query<TAB>count,\n"
	"              0 if absent; lines with non-ACGT bases or the wrong length\n"
	"              print NA.\n"
	"  -o <file>   write to <file> instead of stdout\n");
	exit(status);
}

static void print_info(FILE *out, const kmdb_t *db)
{
	fprintf(out, "#path\t%s\n", db->path);
	fprintf(out, "#k\t%u\n", db->k);
	fprintf(out, "#n_kmers\t%" PRIu64 "\n", db->n_kmers);
	fprintf(out, "#n_blocks\t%" PRIu64 "\n", db->n_blocks);
	fprintf(out, "#block_n\t%u\n", db->block_n);
	fprintf(out, "#flags\t%" PRIu64 "%s\n", db->flags,
	        (db->flags & KMDB_FLAG_CANON_MAX) ? " (canonical=max(fwd,rc))" : "");
	fprintf(out, "#max_comp_size\t%u\n", db->max_comp_size);
	if (db->meta) {
		const char *line = db->meta;
		while (*line) {
			const char *nl = strchr(line, '\n');
			size_t len = nl ? (size_t)(nl - line) : strlen(line);
			const char *eq = memchr(line, '=', len);
			if (eq)
				fprintf(out, "%.*s\t%.*s\n", (int)(eq - line), line,
				        (int)(len - (size_t)(eq - line) - 1), eq + 1);
			else if (len)
				fprintf(out, "%.*s\n", (int)len, line);
			if (!nl) break;
			line = nl + 1;
		}
	}
}

/* ── lookups ─────────────────────────────────────────────────────────── */

/* first block whose last_kmer >= x, or n_blocks */
static uint64_t find_block(const kmdb_t *db, uint64_t x)
{
	uint64_t lo = 0, hi = db->n_blocks;
	while (lo < hi) {
		uint64_t mid = lo + (hi - lo) / 2;
		if (db->blocks[mid].last_kmer < x) lo = mid + 1;
		else hi = mid;
	}
	return lo;
}

/* count of canonical x; 0 if absent, -1 on read error. Caches the block. */
static long long lookup(const kmdb_t *db, kmdb_block_buf *bb, int64_t *cached,
                        long *cached_n, uint64_t x)
{
	uint64_t b = find_block(db, x);
	if (b >= db->n_blocks || x < db->blocks[b].first_kmer) return 0;
	if (*cached != (int64_t)b) {
		long n = kmdb_read_block(db, b, bb);
		if (n < 0) return -1;
		*cached = (int64_t)b;
		*cached_n = n;
	}
	long lo = 0, hi = *cached_n;
	while (lo < hi) {
		long mid = lo + (hi - lo) / 2;
		if (bb->kmers[mid] < x) lo = mid + 1;
		else hi = mid;
	}
	return (lo < *cached_n && bb->kmers[lo] == x) ? (long long)bb->counts[lo] : 0;
}

static int run_queries(FILE *out, const kmdb_t *db, const char *qpath, int header)
{
	FILE *in = strcmp(qpath, "-") ? fopen(qpath, "r") : stdin;
	if (!in) {
		fprintf(stderr, "error: cannot open %s: %s\n", qpath, strerror(errno));
		return 1;
	}
	kmdb_block_buf bb;
	kmdb_block_buf_init(&bb);
	int64_t cached = -1;
	long cached_n = 0;
	int k = (int)db->k, rc = 0;
	char *line = NULL;
	size_t cap = 0;
	ssize_t len;
	if (header) fprintf(out, "#query\tcount\n");
	while ((len = getline(&line, &cap, in)) != -1) {
		while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r' ||
		                   line[len-1] == ' '  || line[len-1] == '\t'))
			line[--len] = '\0';
		if (len == 0 || line[0] == '#') continue;
		uint64_t x;
		if (len != k || kmdb_encode(line, k, &x) != 0) {
			fprintf(out, "%s\tNA\n", line);
			continue;
		}
		if (db->flags & KMDB_FLAG_CANON_MAX) x = kmdb_canonical(x, k);
		long long c = lookup(db, &bb, &cached, &cached_n, x);
		if (c < 0) { rc = 1; break; }
		fprintf(out, "%s\t%lld\n", line, c);
	}
	free(line);
	kmdb_block_buf_free(&bb);
	if (in != stdin) fclose(in);
	return rc;
}

/* ── main ────────────────────────────────────────────────────────────── */

int main(int argc, char *argv[])
{
	int info = 0, summary = 0, raw = 0, keys_only = 0, header = 1;
	long long limit = -1;
	unsigned long long cmin = 0, cmax = UINT32_MAX;
	char *qpath = NULL, *opath = NULL;
	int c;

	while ((c = getopt(argc, argv, "isn:c:C:rKHq:o:h")) != -1) {
		switch (c) {
		case 'i': info = 1; break;
		case 's': summary = 1; break;
		case 'n': limit = atoll(optarg); break;
		case 'c': cmin = strtoull(optarg, NULL, 10); break;
		case 'C': cmax = strtoull(optarg, NULL, 10); break;
		case 'r': raw = 1; break;
		case 'K': keys_only = 1; break;
		case 'H': header = 0; break;
		case 'q': qpath = optarg; break;
		case 'o': opath = optarg; break;
		case 'h': usage(0); break;
		default:  usage(1);
		}
	}
	if (optind != argc - 1) usage(1);
	if (limit < -1) limit = -1;

	kmdb_t *db = kmdb_open(argv[optind]);
	if (!db) return 1;

	FILE *out = opath ? fopen(opath, "w") : stdout;
	if (!out) {
		fprintf(stderr, "error: cannot open %s: %s\n", opath, strerror(errno));
		return 1;
	}
	setvbuf(out, NULL, _IOFBF, OUT_BUF);

	int rc = 0;
	if (info) {
		print_info(out, db);
	} else if (qpath) {
		rc = run_queries(out, db, qpath, header);
	} else {
		kmdb_block_buf bb;
		kmdb_block_buf_init(&bb);
		const int k = (int)db->k;
		char s[KMDB_MAX_K + 1];
		uint64_t printed = 0, n_seen = 0, n_pass = 0, sum = 0;
		uint32_t vmin = UINT32_MAX, vmax = 0;
		if (header && !summary)
			fprintf(out, keys_only ? "#kmer\n" : "#kmer\tcount\n");
		int done = 0;
		for (uint64_t b = 0; b < db->n_blocks && !done; b++) {
			long n = kmdb_read_block(db, b, &bb);
			if (n < 0) { rc = 1; break; }
			for (long j = 0; j < n; j++) {
				uint32_t v = bb.counts[j];
				n_seen++;
				if (v < cmin || v > cmax) continue;
				if (summary) {
					n_pass++; sum += v;
					if (v < vmin) vmin = v;
					if (v > vmax) vmax = v;
					continue;
				}
				if (limit >= 0 && printed >= (uint64_t)limit) { done = 1; break; }
				if (raw) fprintf(out, "%" PRIu64, bb.kmers[j]);
				else { kmdb_decode(bb.kmers[j], k, s); fputs(s, out); }
				if (keys_only) fputc('\n', out);
				else fprintf(out, "\t%u\n", v);
				printed++;
			}
		}
		kmdb_block_buf_free(&bb);
		if (summary && rc == 0) {
			fprintf(out, "n_kmers\t%" PRIu64 "\n", n_seen);
			fprintf(out, "n_kmers_in_count_range\t%" PRIu64 "\n", n_pass);
			fprintf(out, "count_sum\t%" PRIu64 "\n", sum);
			fprintf(out, "count_min\t%u\n", n_pass ? vmin : 0);
			fprintf(out, "count_max\t%u\n", vmax);
			fprintf(out, "count_mean\t%.4f\n", n_pass ? (double)sum / n_pass : 0.0);
		}
	}

	if (fflush(out) != 0 || (opath && fclose(out) != 0)) {
		fprintf(stderr, "error: writing output failed\n");
		rc = 1;
	}
	kmdb_close(db);
	return rc;
}
