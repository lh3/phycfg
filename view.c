#include <stdlib.h>
#include <stdio.h>
#include "phycfg.h"
#include "knhx.h"
#include "khashl.h"
#include "ketopt.h"

KHASHL_SET_INIT(KH_LOCAL, strset_t, strset, kh_cstr_t, kh_hash_str, kh_eq_str)

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
		int i, nl;
		char **list = pc_list_read(list_fn, &nl);
		strset_t *h = strset_init();
		for (i = 0; i < nl; i++) {
			int absent;
			strset_put(h, list[i], &absent);
		}
		for (i = 0; i < tree->n; i++) {
			knhx1_t *p = &tree->a[i];
			if (p->n == 0)
				p->aux = (strset_get(h, p->name) != kh_end(h)) ? 1 : 0;
		}
		strset_destroy(h);
		for (i = 0; i < nl; i++) free(list[i]);
		free(list);
	}

	kstring_t out = {0, 0, 0};
	kn_format(tree->a, tree->n - 1, &out);
	puts(out.s);
	free(out.s);

	for (int i = 0; i < tree->n; i++) {
		if (tree->a[i].n > 0) free(tree->a[i].child);
		free(tree->a[i].name);
	}
	free(tree->a);
	free(tree);
	return 0;
}
