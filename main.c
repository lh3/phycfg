#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "kommon.h"
#include "phycfg.h"
#include "ketopt.h"

int main_view(int argc, char *argv[]);

static int usage(FILE *fp)
{
	fprintf(fp, "Usage: phycfg <command> <arguments>\n");
	fprintf(fp, "Commands:\n");
	fprintf(fp, "  view       view phylogenetic tree\n");
	fprintf(fp, "  version    print the version number\n");
	return fp == stdout? 0 : 1;
}

int main(int argc, char *argv[])
{
	int ret = 0;
	kom_realtime();
	if (argc == 1) return usage(stdout);
	else if (strcmp(argv[1], "view") == 0) ret = main_view(argc-1, argv+1);
	else if (strcmp(argv[1], "version") == 0) {
		printf("%s\n", PC_VERSION);
		return 0;
	} else {
		fprintf(stderr, "ERROR: unknown command '%s'\n", argv[1]);
		return 1;
	}

	if (kom_verbose >= 3 && argc > 2 && ret == 0) {
		int i;
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, PC_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, kom_realtime(), kom_cputime(), kom_peakrss() / 1024.0 / 1024.0 / 1024.0);
	}
	return 0;
}

int main_view(int argc, char *argv[])
{
	const char *list_fn = NULL;
	ketopt_t o = KETOPT_INIT;
	int c;
	while ((c = ketopt(&o, argc, argv, 1, "l:", 0)) >= 0)
		if (c == 'l') list_fn = o.arg;

	if (o.ind == argc) {
		fprintf(stderr, "Usage: phycfg view [options] <tree.nhx.gz>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -l STR   list of leaf names to mark (comma/space-sep, or @file)\n");
		return 1;
	}

	pc_tree_t *tree = pc_tree_read(argv[o.ind]);
	if (tree == NULL) return 1;

	if (list_fn) {
		int nl;
		char **list = pc_list_read(list_fn, &nl);
		pc_tree_mark_leaf(tree, nl, list);
		for (int i = 0; i < nl; i++) free(list[i]);
		free(list);
	}

	pc_tree_t *out = list_fn ? pc_tree_reduce(tree) : tree;
	if (out) {
		char *s = NULL;
		int32_t max = 0;
		pc_tree_format(out, &s, &max);
		puts(s);
		free(s);
		if (out != tree) pc_tree_destroy(out);
	}
	pc_tree_destroy(tree);
	return 0;
}
