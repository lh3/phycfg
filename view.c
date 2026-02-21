#include <stdlib.h>
#include <stdio.h>
#include "phycfg.h"
#include "ketopt.h"

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
		kstring_t s = {0, 0, 0};
		pc_tree_format(&s, out);
		puts(s.s);
		free(s.s);
		if (out != tree) pc_tree_destroy(out);
	}
	pc_tree_destroy(tree);
	return 0;
}
