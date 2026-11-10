#ifndef READHASH_H
#define READHASH_H

/*
	Read (fragment) hash shared by kmer_metagenome_db -m sample --read-hash and
	kmer_strain_detect --read-hash.

	Every fragment (a single-end read, or a read pair when pairs are counted
	together) gets a number in input order: 0, 1, 2, ... Only fragments with at
	least one read of length >= k are numbered, in both programs, so a .kmdb and
	a FASTQ scan of the same reads give the same numbers. The number is mixed
	(splitmix64 finaliser) into a uniform 32-bit hash.

	Per k-mer the programs keep min_hash = the smallest hash over the fragments
	that contain it. k-mers hit by one fragment share that fragment's hash, so
	the number of distinct min_hash values among a strain's hit k-mers estimates
	how many reads produced the hits (exactly, for k-mers with count 1).

	0 is reserved for "no hash" (k-mer not seen, or a database built without
	--read-hash), so a fragment never hashes to 0.
*/

#include <stdint.h>

#define READHASH_NONE 0u

static inline uint32_t read_hash(uint64_t frag)
{
	uint64_t x = frag + 0x9e3779b97f4a7c15ull;
	x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
	x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
	x ^= x >> 31;
	uint32_t h = (uint32_t)(x >> 32);
	return h ? h : 1u;
}

/* min with READHASH_NONE meaning "unset" */
static inline uint32_t read_hash_min(uint32_t cur, uint32_t h)
{
	return (cur == READHASH_NONE || h < cur) ? h : cur;
}

#endif /* READHASH_H */
