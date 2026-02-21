#include <stdlib.h>
#include <stdio.h>
#include <zlib.h>
#include "phycfg.h"
#include "knhx.h"
#include "kseq.h"

KSTREAM_INIT(gzFile, gzread, 65536)

pc_tree_t *pc_tree_read(const char *fn)
{
	gzFile fp;
	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(0, "r");
	if (fp == NULL) {
		fprintf(stderr, "[E::pc_tree_read] failed to open '%s'\n", fn);
		return NULL;
	}

	kstream_t *ks = ks_init(fp);
	kstring_t line = {0, 0, 0};
	pc_tree_t *tree = NULL;

	while (ks_getuntil(ks, KS_SEP_LINE, &line, NULL) >= 0) {
		if (line.l == 0) continue;
		int n, max, error;
		knhx1_t *nodes = kn_parse(line.s, &n, &max, &error, 0);
		if (error)
			fprintf(stderr, "[W::pc_tree_read] parse error (bits: %d)\n", error);
		tree = (pc_tree_t *)calloc(1, sizeof(pc_tree_t));
		tree->n = tree->m = n;
		tree->a = nodes;
		break;  // read the first tree only
	}

	free(line.s);
	ks_destroy(ks);
	gzclose(fp);
	return tree;
}
