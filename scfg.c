#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "kommon.h"
#include "ketopt.h"
#include "phycfg.h"

pc_scfg_t *pc_scfg_new(int32_t n_node, int32_t m)
{
	pc_scfg_t *s;
	int32_t i, n_dbl = (sizeof(pc_scfg_t) * n_node + sizeof(double) - 1) / sizeof(double);
	double *x;
	x = kom_calloc(double, n_dbl + n_node * 3 * m);
	s = (pc_scfg_t*)x;
	x += n_dbl;
	for (i = 0; i < n_node; ++i) {
		s[i].alpha  = x, x += m;
		s[i].alpha2 = x, x += m;
		s[i].beta   = x, x += m;
	}
	return s;
}

double **pc_mat2d_new(int32_t n_row, int32_t n_col)
{
	double *x, **p;
	int32_t i, n_dbl = (sizeof(double*) * n_row + sizeof(double) - 1) / sizeof(double);
	x = kom_calloc(double, n_dbl + n_row * n_col);
	p = (double**)x;
	x += n_dbl;
	for (i = 0; i < n_row; ++i)
		p[i] = x, x += n_col;
	return p;
}

void pc_transmat_init(double **p, int32_t m, const pc_tree_t *t)
{
	int32_t k;
	assert(m > 1);
	for (k = 0; k < t->n_node; ++k) {
		int32_t i, j;
		double diag, off, *pk = p[k];
		if (k < t->n_node - 1) { // not root
			double d = t->node[k]->d;
			d = d > 0.0? d : 1e-3;
			diag = exp(-d), off = (1.0 - diag) / (m - 1);
		} else diag = off = 1.0 / m; // flat prior at the root
		for (i = 0; i < m; ++i)
			for (j = 0; j < m; ++j)
				pk[i * m + j] = i == j ? diag : off;
	}
}

static inline void pc_scfg_emit(int32_t m, int32_t c, double *alpha)
{
	int32_t i;
	if (c < m) {
		for (i = 0; i < m; ++i) alpha[i] = 0.0;
		alpha[c] = 1.0;
	} else {
		for (i = 0; i < m; ++i) alpha[i] = 1.0;
	}
}

double pc_scfg_inside(const pc_tree_t *t, double **p, const pc_msa_t *msa, int32_t pos, pc_scfg_t *sd)
{
	int32_t i, k, a, b, m = msa->m;
	double logh = 0.0, sum;
	for (i = 0; i < t->n_node; ++i) {
		const pc_node_t *v = t->node[i]; /* post-order */
		double *alpha = sd[i].alpha, *alpha2 = sd[i].alpha2, h;
		assert(v->n_child == 0 || v->n_child == 2); // only works for binary trees
		if (v->n_child == 0) { /* leaf: initialize alpha from emission */
			assert(v->seq_id >= 0 && v->seq_id < msa->n_seq);
			pc_scfg_emit(m, msa->msa[pos][v->seq_id], alpha);
		} else { /* internal: alpha(v,a) = prod_k alpha'(child_k, a) */
			for (a = 0; a < m; ++a) alpha[a] = 1.0;
			for (k = 0; k < v->n_child; ++k) {
				double *a2k = sd[v->child[k]->ftime].alpha2;
				for (a = 0; a < m; ++a) alpha[a] *= a2k[a];
			}
		}
		h = 0.0;
		for (a = 0; a < m; ++a) h += alpha[a];
		for (a = 0; a < m; ++a) alpha[a] /= h;
		/* alpha'(v,a) = sum_b p(b|a) * alpha~(v,b); stored in alpha2 for parent */
		for (a = 0; a < m; ++a) { // NB: alpha2[] is not defined at the root, though it is calculated anyway
			double s = 0.0, *pi = p[i];
			for (b = 0; b < m; ++b) s += pi[a*m + b] * alpha[b];
			alpha2[a] = s;
		}
		sd[i].h = h;
		logh += log(h);
	}
	for (a = 0, sum = 0.0; a < m; ++a)
		sum += sd[t->n_node - 1].alpha[a] * p[t->n_node - 1][a];
	return logh + log(sum); // this is equal to logh + log(h(root) * \sum_a alpha~(root,a) * beta~(root,a))
}

void pc_scfg_outside(const pc_tree_t *t, double **p, int32_t m, pc_scfg_t *sd)
{
	int32_t i, a, b, root_idx = t->n_node - 1;
	double *sib = kom_malloc(double, m);

	/* beta~(root,a) = q(a)/h_root; q(a)=p[root_idx][a] */
	for (a = 0; a < m; ++a) sd[root_idx].beta[a] = p[root_idx][a] / sd[root_idx].h;

	/* Pre-order: root (highest index) down to leaves; parent always visited first */
	for (i = root_idx - 1; i >= 0; --i) {
		pc_node_t *u = t->node[i];
		pc_node_t *v = u->parent;
		int32_t k;
		double *beta_u = sd[i].beta, *beta_v = sd[v->ftime].beta;

		/* sib[a] = prod_k alpha'~(sibling_k, a) for all siblings of u */
		for (a = 0; a < m; ++a) sib[a] = 1.0;
		for (k = 0; k < v->n_child; ++k)
			if (v->child[k] != u) {
				double *a2k = sd[v->child[k]->ftime].alpha2;
				for (a = 0; a < m; ++a) sib[a] *= a2k[a];
			}

		/* beta~(u,b) = (1/h_u) * sum_a p(b|a) * beta~(v,a) * sib[a] */
		for (b = 0; b < m; ++b) {
			double s = 0.0;
			for (a = 0; a < m; ++a)
				s += p[i][a*m + b] * beta_v[a] * sib[a];
			beta_u[b] = s / sd[i].h;
		}
	}
	free(sib);
}

double pc_scfg_em_basic(const pc_tree_t *t, double **p, const pc_msa_t *msa, pc_scfg_t *sd)
{
	int32_t i, j, k, a, b, m = msa->m;
	double loglk = 0.0;
	double *cnt = kom_calloc(double, (size_t)t->n_node * m * m);
	double *tmp = kom_calloc(double, m * m);
	double *sib = kom_malloc(double, m);

	/* E step: accumulate sufficient statistics over all alignment columns */
	for (i = 0; i < msa->len; ++i) {
		double sum, *cnt_j;
		loglk += pc_scfg_inside(t, p, msa, i, sd);
		pc_scfg_outside(t, p, m, sd);
		/* for each non-root branch, accumulate cnt[j][a*m+b] =
		 * p(b|a) * alpha~(u,b) * beta~(par,a) * prod_k alpha'~(sib_k,a) */
		for (j = 0; j < t->n_node - 1; ++j) {
			const pc_node_t *u = t->node[j], *v = u->parent;
			double *alpha_u = sd[j].alpha, *beta_v = sd[v->ftime].beta;
			cnt_j = cnt + (size_t)j * m * m;
			for (a = 0; a < m; ++a) sib[a] = 1.0;
			for (k = 0; k < v->n_child; ++k) { // prod_k alpha'~(sib_k,a)
				if (v->child[k] != u) {
					double *a2k = sd[v->child[k]->ftime].alpha2;
					for (a = 0; a < m; ++a) sib[a] *= a2k[a];
				}
			}
			for (a = 0, sum = 0.0; a < m; ++a) {
				double x = beta_v[a] * sib[a];
				for (b = 0; b < m; ++b)
					sum += (tmp[a*m + b] = p[j][a*m + b] * alpha_u[b] * x);
			}
			sum = 1.0 / sum;
			for (a = 0; a < m; ++a)
				for (b = 0; b < m; ++b)
					cnt_j[a*m + b] += tmp[a*m + b] * sum;
		}
		// special case: the root
		assert(j == t->n_node - 1);
		cnt_j = cnt + (size_t)j * m * m;
		for (a = 0, sum = 0.0; a < m; ++a)
			sum += sd[j].alpha[a] * sd[j].beta[a];
		sum = 1.0 / sum;
		for (a = 0; a < m; ++a)
			for (b = 0; b < m; ++b)
				cnt_j[a*m + b] += sd[j].alpha[b] * sd[j].beta[b] * sum;
	}

	/* M step: renormalize each row of p from accumulated counts */
	for (j = 0; j < t->n_node; ++j) {
		double *cnt_j = cnt + (size_t)j * m * m, *p_j = p[j];
		for (a = 0; a < m; ++a) {
			double s = 0.0;
			for (b = 0; b < m; ++b) s += cnt_j[a*m + b];
			if (s > 0.0)
				for (b = 0; b < m; ++b) p_j[a*m + b] = cnt_j[a*m + b] / s;
		}
	}

	free(sib); free(tmp); free(cnt);
	return loglk;
}

int main_scfg(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int32_t i;
	pc_scfg_t *sd;
	double **p;

	while (ketopt(&o, argc, argv, 1, "", 0) >= 0);
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: phycfg scfg [options] <tree.nhx.gz> <aln.mfa.gz>\n");
		return 1;
	}

	pc_tree_t *t = pc_tree_read(argv[o.ind]);
	if (t == NULL) return 1;
	pc_msa_t *msa = pc_msa_read(argv[o.ind + 1]);
	if (msa == NULL) { pc_tree_destroy(t); return 1; }

	pc_msa_encode(msa, pc_msa_infer_rt(msa));
	assert(msa->rt == PC_RT_NT); // only for nucleotide for now
	pc_tree_match_msa(t, msa);
	sd = pc_scfg_new(t->n_node, msa->m);
	p = pc_mat2d_new(t->n_node, msa->m * msa->m);
	pc_transmat_init(p, msa->m, t);
	#if 1
	for (i = 0; i < 10; ++i) {
		double loglk = pc_scfg_em_basic(t, p, msa, sd);
		printf("LK\t%.6f\n", loglk);
	}
	printf("BF");
	for (i = 0; i < msa->m; ++i) printf("\t%.4f", p[t->n_node-1][i]);
	printf("\n");
	#else
	//for (i = 0; i < msa->len; ++i) {
	{
		i = msa->len>>1;
		int32_t j, a;
		pc_scfg_inside(t, p, msa, i, sd);
		pc_scfg_outside(t, p, msa->m, sd);
		for (j = 0; j < t->n_node; ++j) {
			double x = 0.0;
			for (a = 0; a < msa->m; ++a)
				x += sd[j].alpha[a] * sd[j].beta[a];
			printf("%d\t%d\t%g\t%g\n", i, j, x * sd[j].h, sd[j].h);
		}
	}
	#endif

	pc_tree_destroy(t);
	pc_msa_destroy(msa);
	return 0;
}
