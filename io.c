#include <stdlib.h>
#include <stdio.h>
#include <zlib.h>
#include "phycfg.h"
#include "kommon.h"
#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 65536)

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
