#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "pcpriv.h"
#include "kommon.h"

static pc_scfg_data_t *pc_scfg_data_new(int32_t m, int32_t len)
{
	pc_scfg_data_t *s;
	int32_t n_dbl = m * m * 2 + len * m + len * m * 3;
	double *p;
	s = (pc_scfg_data_t*)calloc(1, sizeof(pc_scfg_data_t) + n_dbl * sizeof(double));
	p = s->x;
	s->p = p, p += m * m;
	s->jc = p, p += m * m;
	s->h = p, p += len;
	s->alpha = p, p += len * m;
	s->alpha2 = p, p += len * m;
	s->beta = p;
	return s;
}

void pc_scfg_alloc(pc_tree_t *t, int32_t len)
{
	int32_t i;
	for (i = 0; i < t->n_node; ++i)
		if (t->node[i]->q == 0)
			t->node[i]->q = pc_scfg_data_new(t->m, len);
}

void pc_scfg_free(pc_tree_t *t)
{
	int32_t i;
	for (i = 0; i < t->n_node; ++i) {
		free(t->node[i]->q);
		t->node[i]->q = 0;
	}
}

void pc_scfg_init_par(pc_tree_t *t) // pc_node_t::q MUST BE allocated
{
	int32_t k, m = t->m;
	for (k = 0; k < t->n_node; ++k) {
		int32_t i, j;
		double diag, off, *pk = t->node[k]->q->p;
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

double pc_scfg_inside2(pc_tree_t *t, const pc_msa_t *msa, int32_t pos)
{
	int32_t i, k, a, b, m = t->m;
	double logh = 0.0, sum;
	for (i = 0; i < t->n_node; ++i) {
		const pc_node_t *v = t->node[i]; /* post-order */
		double *alpha = &v->q->alpha[pos * m], *alpha2 = &v->q->alpha2[pos * m], h;
		const double *pi = v->q->p;
		assert(v->n_child == 0 || v->n_child == 2); // only works for binary trees
		if (v->n_child == 0) { /* leaf: initialize alpha from emission */
			assert(v->seq_id >= 0 && v->seq_id < msa->n_seq);
			pc_scfg_emit(m, msa->msa[pos][v->seq_id], alpha);
		} else { /* internal: alpha(v,a) = prod_k alpha'(child_k, a) */
			for (a = 0; a < m; ++a) alpha[a] = 1.0;
			for (k = 0; k < v->n_child; ++k) {
				const double *a2k = &v->child[k]->q->alpha2[pos * m];
				for (a = 0; a < m; ++a) alpha[a] *= a2k[a];
			}
		}
		h = 0.0;
		for (a = 0; a < m; ++a) h += alpha[a];
		for (a = 0; a < m; ++a) alpha[a] /= h;
		/* alpha'(v,a) = sum_b p(b|a) * alpha~(v,b); stored in alpha2 for parent */
		for (a = 0; a < m; ++a) { // NB: alpha2[] is not defined at the root, though it is calculated anyway
			double s = 0.0;
			for (b = 0; b < m; ++b) s += pi[a*m + b] * alpha[b];
			alpha2[a] = s;
		}
		v->q->h[pos] = h;
		logh += log(h);
	}
	for (a = 0, sum = 0.0; a < m; ++a)
		sum += t->node[t->n_node - 1]->q->alpha[pos * m + a] * t->node[t->n_node - 1]->q->p[a];
	return logh + log(sum); // this is equal to logh + log(h(root) * \sum_a alpha~(root,a) * beta~(root,a))
}

void pc_scfg_outside2(pc_tree_t *t, int32_t pos)
{
	int32_t i, a, b, m = t->m;
	double *sib = kom_malloc(double, m);

	/* beta~(root,a) = q(a)/h_root; q(a)=p[root_idx][a] */
	for (a = 0; a < m; ++a)
		t->root->q->beta[pos * m + a] = t->root->q->p[a] / t->root->q->h[pos];

	/* Pre-order: root (highest index) down to leaves; parent always visited first */
	for (i = t->n_node - 2; i >= 0; --i) {
		pc_node_t *u = t->node[i];
		pc_node_t *v = u->parent;
		int32_t k;
		double *beta_u = &u->q->beta[pos * m], *beta_v = &v->q->beta[pos * m];
		const double *pi = u->q->p;

		/* sib[a] = prod_k alpha'~(sibling_k, a) for all siblings of u */
		for (a = 0; a < m; ++a) sib[a] = 1.0;
		for (k = 0; k < v->n_child; ++k)
			if (v->child[k] != u) {
				const double *a2k = &v->child[k]->q->alpha2[pos * m];
				for (a = 0; a < m; ++a) sib[a] *= a2k[a];
			}

		/* beta~(u,b) = (1/h_u) * sum_a p(b|a) * beta~(v,a) * sib[a] */
		for (b = 0; b < m; ++b) {
			double s = 0.0;
			for (a = 0; a < m; ++a)
				s += pi[a*m + b] * beta_v[a] * sib[a];
			beta_u[b] = s / u->q->h[pos];
		}
	}
	free(sib);
}

double pc_scfg_post_cnt2(pc_tree_t *t, const pc_msa_t *msa)
{
	int32_t i, j, k, a, b, m = t->m, m2 = m * m;
	double loglk = 0.0;
	double *tmp = kom_malloc(double, m2);
	double *sib = kom_malloc(double, m);
	for (j = 0; j < t->n_node; ++j)
		memset(t->node[j]->q->jc, 0, m2 * sizeof(double));
	for (i = 0; i < msa->len; ++i) {
		double sum, *jc;
		loglk += pc_scfg_inside2(t, msa, i);
		pc_scfg_outside2(t, i);
		for (j = 0; j < t->n_node - 1; ++j) {
			const pc_node_t *u = t->node[j], *v = u->parent;
			const double *alpha_u = &u->q->alpha[i * m], *beta_v = &v->q->beta[i * m];
			const double *pj = u->q->p;
			jc = u->q->jc;
			for (a = 0; a < m; ++a) sib[a] = 1.0;
			for (k = 0; k < v->n_child; ++k) {
				if (v->child[k] != u) {
					const double *a2k = &v->child[k]->q->alpha2[i * m];
					for (a = 0; a < m; ++a) sib[a] *= a2k[a];
				}
			}
			for (a = 0, sum = 0.0; a < m; ++a) {
				double x = beta_v[a] * sib[a];
				for (b = 0; b < m; ++b)
					sum += (tmp[a*m + b] = pj[a*m + b] * alpha_u[b] * x);
			}
			sum = 1.0 / sum;
			for (a = 0; a < m; ++a)
				for (b = 0; b < m; ++b)
					jc[a*m + b] += tmp[a*m + b] * sum;
		}
		// special case: the root
		assert(j == t->n_node - 1);
		{
			const pc_node_t *r = t->node[j];
			const double *alpha_r = &r->q->alpha[i * m], *beta_r = &r->q->beta[i * m];
			jc = r->q->jc;
			for (a = 0, sum = 0.0; a < m; ++a)
				sum += alpha_r[a] * beta_r[a];
			sum = 1.0 / sum;
			for (a = 0; a < m; ++a)
				for (b = 0; b < m; ++b)
					jc[a*m + b] += alpha_r[b] * beta_r[b] * sum;
		}
	}
	free(sib); free(tmp);
	return loglk;
}

static void pc_scfg_c2p(int32_t m, const double *c, double *p)
{
	int32_t a, b;
	for (a = 0; a < m; ++a) {
		double s = 0.0;
		for (b = 0; b < m; ++b) s += c[a*m + b];
		if (s > 0.0)
			for (b = 0, s = 1.0 / s; b < m; ++b)
				p[a * m + b] = c[a * m + b] * s;
	}
}

static void pc_scfg_cnt2p(pc_tree_t *t, pc_model_t ct)
{
	int32_t u, m = t->m;
	double *tmp = kom_calloc(double, m * m);
	for (u = 0; u < t->n_node; ++u) {
		pc_model_matrix(t->node[u]->q->jc, m, ct, tmp);
		pc_scfg_c2p(m, tmp, t->node[u]->q->p);
	}
	free(tmp);
}

double pc_scfg_em2(pc_tree_t *t, const pc_msa_t *msa, pc_model_t ct)
{
	double loglk = pc_scfg_post_cnt2(t, msa);
	pc_scfg_cnt2p(t, ct);
	return loglk;
}
