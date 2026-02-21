#include <stdlib.h>
#include <stdio.h>
#include <zlib.h>
#include "phycfg.h"
#include "knhx.h"
#include "kseq.h"

KSTREAM_INIT(gzFile, gzread, 65536)

int main_view(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "Usage: phycfg view <tree.nhx.gz>\n");
		return 1;
	}

	gzFile fp = gzopen(argv[1], "r");
	if (fp == NULL) {
		fprintf(stderr, "[E::main_view] failed to open '%s'\n", argv[1]);
		return 1;
	}

	kstream_t *ks = ks_init(fp);
	kstring_t line = {0, 0, 0};

	while (ks_getuntil(ks, KS_SEP_LINE, &line, NULL) >= 0) {
		if (line.l == 0) continue;

		int n, error;
		knhx1_t *nodes = kn_parse(line.s, &n, &error);
		if (error)
			fprintf(stderr, "[W::main_view] parse error (bits: %d)\n", error);

		kstring_t out = {0, 0, 0};
		kn_format(nodes, n - 1, &out);
		puts(out.s);
		free(out.s);

		for (int i = 0; i < n; i++) {
			if (nodes[i].n > 0) free(nodes[i].child);
			free(nodes[i].name);
		}
		free(nodes);
	}

	free(line.s);
	ks_destroy(ks);
	gzclose(fp);
	return 0;
}
