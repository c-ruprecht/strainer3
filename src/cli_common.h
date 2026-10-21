#ifndef CLI_COMMON_H
#define CLI_COMMON_H
/*
	Shared command-line conventions for the strainer C tools.

	Flags with the same meaning in every tool:
	  -t <n>     threads
	  -A <path>  genome list (one path per line) or a single genome
	  -B <path>  metagenome / target list
	  -O <dir>   output directory (default ".", created if missing)
	  -n <name>  output basename (default: derived from the main input)
	  -o <path>  explicit output, overrides -O/-n. Single-output tools take
	             a full file path; multi-output tools take a path prefix
	             that every output file name is appended to.
	  -P <n>     keep only k-mers present in at most n samples (0 = no cap)
	  -z <lvl>   compression level (0 = off where the output may be plain)
	  -h         usage to stderr, exit 0 (bad usage exits 1)

	Header-only (static functions) so tools that link nothing else, like
	genome_compare_presence, can use it too.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>

/* Filename of `path` with directories and known sequence/table/compression
   extensions stripped: ref/GCF_000123.4_genomic.fna.gz -> GCF_000123.4_genomic.
   Extensions are whitelisted rather than "strip after the last dot" so that
   accession-style names keep their version suffix. Returns NULL if nothing
   usable is left. Caller frees. */
static inline char *cli_derive_basename(const char *path)
{
	static const char *comp_ext[] = { ".gz", ".zst", ".bz2", ".xz", NULL };
	static const char *seq_ext[]  = { ".fa", ".fna", ".fasta", ".ffn", ".faa",
	                                  ".fas", ".fsa", ".seq", ".txt", ".list",
	                                  ".fastq", ".fq", ".kmdb", ".tsv", ".csv",
	                                  NULL };

	const char *slash = strrchr(path, '/');
	const char *name  = slash ? slash + 1 : path;
	if (*name == '\0') return NULL;

	char *base = strdup(name);
	if (!base) { perror("strdup basename"); exit(EXIT_FAILURE); }

	/* Strip at most one compression suffix, then one sequence suffix. */
	for (int round = 0; round < 2; round++) {
		const char **table = (round == 0) ? comp_ext : seq_ext;
		size_t len = strlen(base);
		for (const char **e = table; *e; e++) {
			size_t elen = strlen(*e);
			if (len > elen && strcasecmp(base + len - elen, *e) == 0) {
				base[len - elen] = '\0';
				break;
			}
		}
	}

	if (base[0] == '\0') { free(base); return NULL; }
	return base;
}

/* <dir>/<base><suffix>, with a trailing slash on dir collapsed. Caller frees. */
static inline char *cli_build_path(const char *dir, const char *base,
                                   const char *suffix)
{
	size_t dir_len = strlen(dir);
	while (dir_len > 1 && dir[dir_len - 1] == '/') dir_len--;

	size_t need = dir_len + 1 + strlen(base) + strlen(suffix) + 1;
	char *path = malloc(need);
	if (!path) { perror("malloc path"); exit(EXIT_FAILURE); }

	snprintf(path, need, "%.*s/%s%s", (int)dir_len, dir, base, suffix);
	return path;
}

/* mkdir -p. Returns 0 on success (including "already a directory"). */
static inline int cli_ensure_dir(const char *dir)
{
	struct stat st;
	if (stat(dir, &st) == 0)
		return S_ISDIR(st.st_mode) ? 0 : (errno = ENOTDIR, -1);

	char *tmp = strdup(dir);
	if (!tmp) return -1;

	size_t len = strlen(tmp);
	while (len > 1 && tmp[len - 1] == '/') tmp[--len] = '\0';

	for (char *p = tmp + 1; *p; p++) {
		if (*p != '/') continue;
		*p = '\0';
		if (mkdir(tmp, 0777) != 0 && errno != EEXIST) { free(tmp); return -1; }
		*p = '/';
	}
	if (mkdir(tmp, 0777) != 0 && errno != EEXIST) { free(tmp); return -1; }

	free(tmp);
	return 0;
}

/* Directory part of `path` ("." if none). Caller frees. */
static inline char *cli_dirname(const char *path)
{
	const char *slash = strrchr(path, '/');
	if (!slash) return strdup(".");
	if (slash == path) return strdup("/");
	return strndup(path, (size_t)(slash - path));
}

/* Resolve the -O/-n/-o convention to one output path or prefix:
     -o given        -> o_arg verbatim (its directory is created)
     otherwise       -> <out_dir or ".">/<name><suffix>
   `name` must already be set (explicit -n, or derived by the caller).
   Exits with a message on a bad -n or an uncreatable directory.
   Caller frees. */
static inline char *cli_resolve_output(const char *tool,
                                       const char *out_dir,
                                       const char *name,
                                       const char *o_arg,
                                       const char *suffix)
{
	if (o_arg) {
		char *d = cli_dirname(o_arg);
		if (cli_ensure_dir(d) != 0) {
			fprintf(stderr, "%s: error: could not create directory %s: %s\n",
			        tool, d, strerror(errno));
			exit(EXIT_FAILURE);
		}
		free(d);
		return strdup(o_arg);
	}
	if (!name || !*name) {
		fprintf(stderr, "%s: error: no output name; pass -n or -o\n", tool);
		exit(EXIT_FAILURE);
	}
	if (strchr(name, '/') != NULL) {
		fprintf(stderr, "%s: error: -n must be a bare basename, not a path "
		        "(got %s); use -O for the directory or -o for a full path\n",
		        tool, name);
		exit(EXIT_FAILURE);
	}
	const char *dir = out_dir ? out_dir : ".";
	if (cli_ensure_dir(dir) != 0) {
		fprintf(stderr, "%s: error: could not create output directory %s: %s\n",
		        tool, dir, strerror(errno));
		exit(EXIT_FAILURE);
	}
	return cli_build_path(dir, name, suffix);
}

/* Parse a non-negative integer flag value strictly (whole string must be
   digits). Exits with a message otherwise — catches e.g. "-z -F" where an
   old boolean flag now takes a value. */
static inline long cli_parse_nonneg(const char *tool, char flag, const char *s)
{
	char *end = NULL;
	errno = 0;
	long v = strtol(s, &end, 10);
	if (errno || end == s || *end != '\0' || v < 0) {
		fprintf(stderr, "%s: error: -%c expects a non-negative integer "
		        "(got \"%s\")\n", tool, flag, s);
		exit(EXIT_FAILURE);
	}
	return v;
}

/* For a flag that was renamed: explain and exit 1. */
static inline void cli_renamed_flag(const char *tool, char old_flag,
                                    const char *new_flag)
{
	fprintf(stderr, "%s: error: -%c was renamed to %s (see -h)\n",
	        tool, old_flag, new_flag);
	exit(EXIT_FAILURE);
}

#endif
