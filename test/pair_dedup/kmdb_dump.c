/*
	kmdb_dump — test helper: print a .kmdb as "<kmer>\t<count>" lines
	(k-mer decoded to ACGT, sorted by key), or with -m a metadata key.

	  kmdb_dump x.kmdb            all k-mers with counts
	  kmdb_dump -m kmer_positions x.kmdb

	Built by run_tests.sh against ../../src/kmdb.c.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "kmdb.h"

int main(int argc, char **argv)
{
	if (argc == 4 && strcmp(argv[1], "-m") == 0) {
		kmdb_t *db = kmdb_open(argv[3]);
		if (!db) return 1;
		char *v = kmdb_meta_get(db, argv[2]);
		printf("%s\n", v ? v : "NA");
		free(v);
		kmdb_close(db);
		return 0;
	}
	if (argc != 2) {
		fprintf(stderr, "usage: kmdb_dump [-m key] <x.kmdb>\n");
		return 1;
	}
	kmdb_t *db = kmdb_open(argv[1]);
	if (!db) return 1;
	kmdb_block_buf b;
	kmdb_block_buf_init(&b);
	char s[40];
	const char *nt = "ACGT";
	for (uint64_t i = 0; i < db->n_blocks; i++) {
		long n = kmdb_read_block(db, i, &b);
		if (n < 0) return 1;
		for (long j = 0; j < n; j++) {
			uint64_t x = b.kmers[j];
			for (int p = (int)db->k - 1; p >= 0; p--) { s[p] = nt[x & 3]; x >>= 2; }
			s[db->k] = '\0';
			printf("%s\t%" PRIu32 "\n", s, b.counts[j]);
		}
	}
	kmdb_block_buf_free(&b);
	kmdb_close(db);
	return 0;
}
