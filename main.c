#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "kommon.h"
#include "phycfg.h"
#include "ketopt.h"

int main_view(int argc, char *argv[]);
int main_msaflt(int argc, char *argv[]);
int main_reroot(int argc, char *argv[]);
int main_scfg(int argc, char *argv[]);

static int usage(FILE *fp)
{
	fprintf(fp, "Usage: phycfg <command> <arguments>\n");
	fprintf(fp, "Commands:\n");
	fprintf(fp, "  view       view phylogenetic tree\n");
	fprintf(fp, "  msaflt     filter MSA columns by non-gap residue count\n");
	fprintf(fp, "  reroot     reroot a tree (mid-point by default)\n");
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
	else if (strcmp(argv[1], "reroot") == 0) ret = main_reroot(argc-1, argv+1);
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
	int32_t min_cnt = 1, codon_flag = 0;
	ketopt_t o = KETOPT_INIT;
	int c;
	while ((c = ketopt(&o, argc, argv, 1, "m:123", 0)) >= 0) {
		if      (c == 'm') min_cnt = atoi(o.arg);
		else if (c == '1') codon_flag |= 1;
		else if (c == '2') codon_flag |= 2;
		else if (c == '3') codon_flag |= 4;
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: phycfg msaflt [options] <aln.mfa.gz>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -m INT   min non-gap/non-ambiguous residues per column [1]\n");
		fprintf(stderr, "  -1/-2/-3 select 1st/2nd/3rd codon position (combinable)\n");
		return 1;
	}

	pc_msa_t *msa = pc_msa_read(argv[o.ind]);
	if (msa == NULL) return 1;

	pc_restype_t rt = pc_msa_infer_rt(msa);
	if (kom_verbose >= 3)
		fprintf(stderr, "[M::%s] %s alignment\n", __func__,
		        rt == PC_RT_CODON ? "codon" : rt == PC_RT_NT ? "nucleotide" : rt == PC_RT_AA ? "amino acid" : "unknown");
	pc_msa_encode(msa, rt);
	pc_msa_filter(msa, min_cnt);
	if (codon_flag) pc_msa_select_codon(msa, codon_flag);

	int32_t i, j;
	for (j = 0; j < msa->n_seq; ++j) {
		printf(">%s\n", msa->name[j]);
		for (i = 0; i < msa->len; ++i) {
			uint8_t v = msa->msa[i][j];
			char ch;
			if (rt == PC_RT_NT || rt == PC_RT_CODON)
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

int main_reroot(int argc, char *argv[])
{
	const char *list_fn = NULL;
	ketopt_t o = KETOPT_INIT;
	int c;
	while ((c = ketopt(&o, argc, argv, 1, "l:", 0)) >= 0) {
		if (c == 'l') list_fn = o.arg;
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: phycfg reroot [options] <tree.nhx.gz>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -l STR   root at midpoint of branch to LCA of listed leaves []\n");
		return 1;
	}

	pc_tree_t *tree = pc_tree_read(argv[o.ind]);
	if (tree == NULL) return 1;

	int32_t nid;
	double dist = -1.0;  /* -1 tells pc_tree_reroot_core to use branch midpoint */

	if (list_fn) {
		/* Root at midpoint of the branch leading to the LCA of the given leaves */
		int nl;
		char **list = pc_list_read(list_fn, &nl);
		uint8_t *mark = kom_calloc(uint8_t, tree->n_node);
		pc_tree_mark_leaf(tree, nl, list, mark);
		for (int i = 0; i < nl; i++) free(list[i]);
		free(list);
		nid = pc_tree_lca(tree, mark);
		free(mark);
		if (nid < 0) {
			fprintf(stderr, "[E::%s] no leaves matched\n", __func__);
			pc_tree_destroy(tree);
			return 1;
		}
	} else {
		/* Global midpoint root: find the node on the diameter midpoint branch */
		nid = pc_tree_mid_longest(tree, &dist);
		if (nid < 0) {
			fprintf(stderr, "[E::%s] tree has fewer than two leaves\n", __func__);
			pc_tree_destroy(tree);
			return 1;
		}
	}

	pc_tree_reroot(tree, nid, dist);

	char *s = NULL;
	int32_t max = 0;
	pc_tree_format(tree, &s, &max);
	puts(s);
	free(s);
	pc_tree_destroy(tree);
	return 0;
}

int main_scfg(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int32_t i, max_iter = 100, max_iter_br = 50, nni = 0, skip_dist = 0, max_str = 0, test_mode = 0;
	pc_model_t md = PC_MD_FULL, md_test = PC_MD_UNDEF, md_EM = PC_MD_UNDEF;
	pc_scfg_buf_t *sd;
	double loglk;
	char *str = 0;

	while (ketopt(&o, argc, argv, 1, "e:b:n:m:t:Dx", 0) >= 0) {
		if (o.opt == 'n') nni = atoi(o.arg);
		else if (o.opt == 'e') max_iter = atoi(o.arg);
		else if (o.opt == 'b') max_iter_br = atoi(o.arg);
		else if (o.opt == 'm') md = pc_model_from_str(o.arg);
		else if (o.opt == 't') md_test = pc_model_from_str(o.arg);
		else if (o.opt == 'D') skip_dist = 1;
		else if (o.opt == 'x') test_mode = 1;
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: phycfg scfg [options] <tree.nhx.gz> <aln.mfa.gz>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -m STR    model: full, rev/GTR or TN93 [full]\n");
		fprintf(stderr, "  -n INT    max NNI topology search rounds [%d]\n", nni);
		fprintf(stderr, "  -e INT    number of EM iterations for the whole tree [%d]\n", max_iter);
		fprintf(stderr, "  -b INT    number of EM iterations per branch [%d]\n", max_iter_br);
		fprintf(stderr, "  -t STR    test model; if set, use this model for initial EM []\n");
		fprintf(stderr, "  -D        don't recalculate branch lengths for DNA sequences\n");
		return 1;
	}

	pc_tree_t *t = pc_tree_read(argv[o.ind]);
	if (t == NULL) return 1;
	pc_msa_t *msa = pc_msa_read(argv[o.ind + 1]);
	if (msa == NULL) { pc_tree_destroy(t); return 1; }

	pc_msa_encode(msa, pc_msa_infer_rt(msa));
	assert(msa->rt == PC_RT_NT || msa->rt == PC_RT_CODON); // only for nucleotide for now
	pc_tree_match_msa(t, msa);

	sd = pc_scfg_buf_new(t->n_node, t->m);
	pc_transmat_init(t);
	md_EM = md_test != PC_MD_UNDEF? md_test : md;
	for (i = 0; i < max_iter && !test_mode; ++i) {
		loglk = pc_scfg_em(t, msa, md_EM, sd);
		fprintf(stderr, "LK\t%d\t%.6f\n", i, loglk);
	}

	if (test_mode) {
		pc_scfg_alloc(t, msa->len);
		pc_scfg_init_par(t);
		for (i = 0; i < max_iter; ++i) {
			loglk = pc_scfg_em2(t, msa, md_EM);
			//fprintf(stderr, "LK2\t%d\t%.6f\n", i, loglk);
		}
		if (nni > 0) {
			int32_t k;
			for (k = 0; k < nni; ++k) {
				double diff = pc_scfg_nni4(t, msa, md, max_iter_br);
				if (diff == 0.0) break;
				for (i = 0; i < max_iter; ++i) {
					loglk = pc_scfg_em2(t, msa, md);
					fprintf(stderr, "XX2\t%d\t%.6f\n", i, loglk);
				}
				fprintf(stderr, "NI2\t%d\t%.6f\t%.6f\n", k + 1, loglk, diff);
			}
		} else if (md_test != PC_MD_UNDEF) {
			double *diff;
			diff = kom_calloc(double, t->n_node);
			pc_scfg_model_cmp2(t, msa, md, md_test, max_iter_br, diff);
			for (i = 0; i < t->n_node; ++i) {
				const pc_node_t *v = t->node[i];
				fprintf(stderr, "CD\t%d\t%d\t%d\t%s\t%.6f\t%.2e\t%.2f\n", i, v->n_child, v->parent ? v->parent->ftime : -1,
						v->name && v->name[0] ? v->name : ".", diff[i], pc_model_lrt(md_test, md, t->m, diff[i]),
						pc_model_BIC(md_test, md, t->m, msa->len, diff[i]));
			}
			free(diff);
		}
	} else if (nni > 0) {
		int32_t k;
		for (k = 0; k < nni; ++k) {
			double diff = pc_scfg_nni(t, msa, md, max_iter_br);
			if (diff == 0.0) break;
			for (i = 0; i < max_iter; ++i)
				loglk = pc_scfg_em(t, msa, md, sd);
			fprintf(stderr, "NI\t%d\t%.6f\t%.6f\n", k + 1, loglk, diff);
		}
	} else if (md_test != PC_MD_UNDEF) {
		double *diff;
		diff = kom_calloc(double, t->n_node);
		pc_scfg_model_cmp(t, msa, md, md_test, max_iter_br, diff);
		for (i = 0; i < t->n_node; ++i) {
			const pc_node_t *v = t->node[i];
			fprintf(stderr, "CD\t%d\t%d\t%d\t%s\t%.6f\t%.2e\t%.2f\n", i, v->n_child, v->parent ? v->parent->ftime : -1,
			        v->name && v->name[0] ? v->name : ".", diff[i], pc_model_lrt(md_test, md, t->m, diff[i]),
					pc_model_BIC(md_test, md, t->m, msa->len, diff[i]));
		}
		free(diff);
	}
	if (!skip_dist && t->m == 4)
		pc_model_dist(t, msa, PC_MD_TN93);
	pc_tree_format(t, &str, &max_str);
	puts(str);
	free(str);

	free(sd);
	pc_tree_destroy(t);
	pc_msa_destroy(msa);
	return 0;
}
