#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "kommon.h"
#include "phycfg.h"
#include "ketopt.h"

int main_view(int argc, char *argv[]);
int main_msaflt(int argc, char *argv[]);
int main_scfg(int argc, char *argv[]);

static int usage(FILE *fp)
{
	fprintf(fp, "Usage: phycfg <command> <arguments>\n");
	fprintf(fp, "Commands:\n");
	fprintf(fp, "  view       view phylogenetic tree\n");
	fprintf(fp, "  msaflt     filter MSA columns by non-gap residue count\n");
	fprintf(fp, "  scfg       SCFG-based phylogenetic analysis\n");
	fprintf(fp, "  version    print the version number\n");
	return fp == stdout? 0 : 1;
}

int main(int argc, char *argv[])
{
	int ret = 0;
	kom_realtime();
	if (argc == 1) return usage(stdout);
	else if (strcmp(argv[1], "view") == 0) ret = main_view(argc-1, argv+1);
	else if (strcmp(argv[1], "msaflt") == 0) ret = main_msaflt(argc-1, argv+1);
	else if (strcmp(argv[1], "scfg") == 0) ret = main_scfg(argc-1, argv+1);
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
	int print_leaf = 0;
	ketopt_t o = KETOPT_INIT;
	int c;
	while ((c = ketopt(&o, argc, argv, 1, "Ll:", 0)) >= 0) {
		if      (c == 'l') list_fn = o.arg;
		else if (c == 'L') print_leaf = 1;
	}

	if (o.ind == argc) {
		fprintf(stderr, "Usage: phycfg view [options] <tree.nhx.gz>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -L       print leaf names, one per line\n");
		fprintf(stderr, "  -l STR   list of leaf names to mark (comma/space-sep, or @file)\n");
		return 1;
	}

	pc_tree_t *tree = pc_tree_read(argv[o.ind]);
	if (tree == NULL) return 1;

	uint8_t *mark = NULL;
	if (list_fn) {
		int nl;
		char **list = pc_list_read(list_fn, &nl);
		mark = kom_calloc(uint8_t, tree->n_node);
		pc_tree_mark_leaf(tree, nl, list, mark);
		for (int i = 0; i < nl; i++) free(list[i]);
		free(list);
	}

	pc_tree_t *out = list_fn ? pc_tree_reduce(tree, mark) : tree;
	if (out) {
		if (print_leaf) {
			int32_t i;
			for (i = 0; i < out->n_node; ++i)
				if (out->node[i]->n_child == 0)
					puts(out->node[i]->name);
		} else {
			char *s = NULL;
			int32_t max = 0;
			pc_tree_format(out, &s, &max);
			puts(s);
			free(s);
		}
		if (out != tree) pc_tree_destroy(out);
	}
	free(mark);
	pc_tree_destroy(tree);
	return 0;
}

int main_msaflt(int argc, char *argv[])
{
	int32_t min_cnt = 1, is_cds = 0;
	ketopt_t o = KETOPT_INIT;
	int c;
	while ((c = ketopt(&o, argc, argv, 1, "m:c", 0)) >= 0) {
		if      (c == 'm') min_cnt = atoi(o.arg);
		else if (c == 'c') is_cds = 1;
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: phycfg msaflt [options] <aln.mfa.gz>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -m INT   min non-gap/non-ambiguous residues per column [1]\n");
		fprintf(stderr, "  -c       treat as CDS; filter whole codons\n");
		return 1;
	}

	pc_msa_t *msa = pc_msa_read(argv[o.ind]);
	if (msa == NULL) return 1;

	pc_restype_t rt = pc_msa_infer_rt(msa);
	pc_msa_encode(msa, rt);
	pc_msa_filter(msa, min_cnt, is_cds);

	int32_t i, j;
	for (j = 0; j < msa->n_seq; ++j) {
		printf(">%s\n", msa->name[j]);
		for (i = 0; i < msa->len; ++i) {
			uint8_t v = msa->msa[i][j];
			char ch;
			if (rt == PC_RT_NT)
				ch = v < 4 ? "ACGT"[v] : v == 4 ? 'N' : '-';
			else if (rt == PC_RT_AA)
				ch = v <= 21 ? kom_aa_i2c[v] : v == 22 ? 'X' : '-';
			else
				ch = (char)v; /* RT_UNKNOWN: still raw ASCII */
			putchar(ch);
		}
		putchar('\n');
	}
	return 0;
}
