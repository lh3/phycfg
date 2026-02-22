#include <stdlib.h>
#include <stdio.h>
#include <zlib.h>
#include "phycfg.h"
#include "kommon.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

pc_tree_t *pc_tree_read(const char *fn)
{
	gzFile fp;
	kstream_t *ks;
	kstring_t buf = {0, 0, 0};
	pc_tree_t *tree;
	fp = fn && strcmp(fn, "-") ? gzopen(fn, "r") : gzdopen(0, "r");
	if (fp == NULL) {
		fprintf(stderr, "[E::pc_tree_read] failed to open '%s'\n", fn);
		return NULL;
	}
	ks = ks_init(fp);
	while (ks_getuntil2(ks, KS_SEP_LINE, &buf, NULL, 1) >= 0) {}
	ks_destroy(ks);
	gzclose(fp);
	tree = buf.l > 0 ? pc_tree_parse(buf.s, NULL) : NULL;
	free(buf.s);
	return tree;
}

pc_msa_t *pc_msa_read(const char *fn)
{
	gzFile fp;
	kseq_t *ks;
	int32_t i, j, n = 0, m = 0, n_pos;
	char **name = NULL, **seq = NULL;
	pc_msa_t *msa = NULL;

	fp = fn && strcmp(fn, "-") ? gzopen(fn, "r") : gzdopen(0, "r");
	if (fp == NULL) {
		fprintf(stderr, "[E::pc_msa_read] failed to open '%s'\n", fn);
		return NULL;
	}
	ks = kseq_init(fp);
	while (kseq_read(ks) >= 0) {
		kom_grow(char*, name, n, m);
		seq = kom_realloc(char*, seq, m);
		name[n] = kom_strdup(ks->name.s);
		seq[n] = kom_strdup(ks->seq.s);
		++n;
	}
	kseq_destroy(ks);
	gzclose(fp);

	if (n == 0) goto end_msa_read;
	n_pos = (int32_t)strlen(seq[0]);
	for (i = 1; i < n; ++i) {
		if ((int32_t)strlen(seq[i]) != n_pos) {
			fprintf(stderr, "[E::pc_msa_read] '%s' has length %d, expected %d\n",
			        name[i], (int)strlen(seq[i]), (int)n_pos);
			goto end_msa_read;
		}
	}
	msa = kom_calloc(pc_msa_t, 1);
	msa->n_seq = n;
	msa->len = n_pos;
	msa->rt = PC_RT_UNKNOWN;
	msa->m = 256;
	msa->name = name;
	name = NULL; /* ownership transferred */
	msa->msa = kom_malloc(uint8_t*, n_pos);
	for (i = 0; i < n_pos; ++i) {
		msa->msa[i] = kom_malloc(uint8_t, n);
		for (j = 0; j < n; ++j)
			msa->msa[i][j] = (uint8_t)seq[j][i];
	}

end_msa_read:
	if (name) {
		for (i = 0; i < n; ++i) free(name[i]);
		free(name);
	}
	for (i = 0; i < n; ++i) free(seq[i]);
	free(seq);
	return msa;
}

char **pc_list_read(const char *o, int *n_)
{
	int n = 0, m = 0;
	char **s = 0;
	*n_ = 0;
	if (*o != '@') {
		const char *q = o, *p;
		for (p = q;; ++p) {
			if (*p == ',' || *p == ' ' || *p == '\t' || *p == 0) {
				if (p - q > 0) {
					kom_grow(char*, s, n, m);
					s[n++] = kom_strndup(q, p - q);
				}
				if (*p == 0) break;
				q = p + 1;
			}
		}
	} else {
		gzFile fp;
		kstream_t *ks;
		kstring_t str = {0,0,0};
		int dret;

		fp = gzopen(o + 1, "r");
		if (fp == 0) return 0;
		ks = ks_init(fp);
		while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
			char *p;
			for (p = str.s; *p && !isspace(*p); ++p);
			kom_grow(char*, s, n, m);
			s[n++] = kom_strndup(str.s, p - str.s);
		}
		ks_destroy(ks);
		gzclose(fp);
	}
	if (s) s = (char**)realloc(s, n * sizeof(char*));
	*n_ = n;
	return s;
}
