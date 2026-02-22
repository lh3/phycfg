#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "kommon.h"
#include "phycfg.h"

pc_scfg_t *pc_scfg_new(int32_t m)
{
	pc_scfg_t *s;
	s = (pc_scfg_t*)malloc(sizeof(pc_scfg_t) + 3 * m * sizeof(double));
	s->alpha  = s->x;
	s->alpha2 = s->x + m;
	s->beta   = s->x + 2 * m;
	return s;
}

pc_transmat_t *pc_transmat_new(int32_t m, double d)
{
	int32_t i, j;
	double diag, off;
	pc_transmat_t *q;
	assert(m > 1);
	diag = exp(-d);
	off  = (1.0 - diag) / (m - 1);
	q = (pc_transmat_t*)malloc(sizeof(pc_transmat_t) + m * m * sizeof(double));
	for (i = 0; i < m; ++i)
		for (j = 0; j < m; ++j)
			q->p[i*m + j] = (i == j) ? diag : off;
	return q;
}

void pc_scfg_emit(int32_t m, int32_t c, double *alpha)
{
	int32_t i;
	if (c < m) {
		for (i = 0; i < m; ++i) alpha[i] = 0.0;
		alpha[c] = 1.0;
	} else {
		for (i = 0; i < m; ++i) alpha[i] = 1.0;
	}
}

double pc_scfg_inside(const pc_tree_t *t, const pc_transmat_t *tm, const pc_msa_t *msa, int32_t pos, pc_scfg_t *sd)
{
	int32_t i, k, a, m = msa->m;
	double log_lik = 0.0;
	for (i = 0; i < t->n_node; ++i) {
		pc_node_t *v = t->node[i]; /* post-order */
		double *alpha = sd[i].alpha, *alpha2 = sd[i].alpha2, h;
		if (v->n_child == 0) { /* leaf: initialize alpha from emission */
			int32_t c = v->seq_id >= 0 ? (int32_t)msa->msa[pos][v->seq_id] : m;
			pc_scfg_emit(m, c, alpha);
		} else { /* internal: alpha(v,a) = prod_k alpha'(child_k, a) */
			for (a = 0; a < m; ++a) alpha[a] = 1.0;
			for (k = 0; k < v->n_child; ++k) {
				double *a2k = sd[v->child[k]->ftime].alpha2;
				for (a = 0; a < m; ++a) alpha[a] *= a2k[a];
			}
		}
		if (v->parent == NULL) /* root: apply flat prior q(a) = 1/m */
			for (a = 0; a < m; ++a) alpha[a] /= m;
		h = 0.0;
		for (a = 0; a < m; ++a) h += alpha[a];
		for (a = 0; a < m; ++a) alpha[a] /= h;
		/* alpha'(v,a) = sum_b p(b|a) * alpha~(v,b); stored in alpha2 for parent */
		for (a = 0; a < m; ++a) {
			int32_t b;
			double s = 0.0;
			for (b = 0; b < m; ++b) s += tm->p[a*m + b] * alpha[b];
			alpha2[a] = s;
		}
		sd[i].h = h;
		log_lik += log(h);
	}
	return log_lik;
}

double pc_scfg_outside(const pc_tree_t *t, const pc_transmat_t *tm, int32_t m, pc_scfg_t *sd)
{
	int32_t i, a, b, root_idx = t->n_node - 1;
	double log_lik = 0.0;
	double *sib = kom_malloc(double, m);

	/* beta~(root,a) = q(a)/h_root = 1/(m*h_root); uniform since q is flat */
	{
		double inv_mh = 1.0 / (m * sd[root_idx].h);
		for (a = 0; a < m; ++a) sd[root_idx].beta[a] = inv_mh;
	}

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
				s += tm->p[a*m + b] * beta_v[a] * sib[a];
			beta_u[b] = s / sd[i].h;
		}
	}
	free(sib);

	/* log P(T) = sum_v log h(v); h_root already incorporates prior 1/m */
	for (i = 0; i < t->n_node; ++i) log_lik += log(sd[i].h);
	return log_lik;
}
