#include <stdlib.h>
#include <stdio.h>
#include "phycfg.h"
#include "knhx.h"

int main_view(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "Usage: phycfg view <tree.nhx.gz>\n");
		return 1;
	}

	pc_tree_t *tree = pc_tree_read(argv[1]);
	if (tree == NULL) return 1;

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
