#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "kommon.h"
#include "pcpriv.h"

pc_scfg_buf_t *pc_scfg_buf_new(int32_t n_node, int32_t m)
{
	pc_scfg_buf_t *s;
	int32_t i, n_dbl = (sizeof(pc_scfg_buf_t) * n_node + sizeof(double) - 1) / sizeof(double);
	double *x;
	x = kom_calloc(double, n_dbl + n_node * 3 * m);
	s = (pc_scfg_buf_t*)x;
	x += n_dbl;
	for (i = 0; i < n_node; ++i) {
		s[i].alpha  = x, x += m;
		s[i].alpha2 = x, x += m;
		s[i].beta   = x, x += m;
	}
	return s;
}

// shape of t->p: (n,m,m), where n is the number of nodes and m is the size of the alphabet
void pc_transmat_init(pc_tree_t *t)
{
	int32_t k, m = t->m;
	assert(m > 1);
	if (t->p == NULL)
		t->p = kom_calloc(double, (size_t)t->n_node * m * m);
	for (k = 0; k < t->n_node; ++k) {
		int32_t i, j;
		double diag, off, *pk = t->p + (size_t)k * m * m;
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

double pc_scfg_inside(const pc_tree_t *t, const pc_msa_t *msa, int32_t pos, pc_scfg_buf_t *sd)
{
	int32_t i, k, a, b, m = t->m;
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
			const double *pi = t->p + (size_t)i * m * m;
			double s = 0.0;
			for (b = 0; b < m; ++b) s += pi[a*m + b] * alpha[b];
			alpha2[a] = s;
		}
		sd[i].h = h;
		logh += log(h);
	}
	for (a = 0, sum = 0.0; a < m; ++a)
		sum += sd[t->n_node - 1].alpha[a] * t->p[(size_t)(t->n_node - 1) * m * m + a];
	return logh + log(sum); // this is equal to logh + log(h(root) * \sum_a alpha~(root,a) * beta~(root,a))
}

void pc_scfg_outside(const pc_tree_t *t, pc_scfg_buf_t *sd)
{
	int32_t i, a, b, m = t->m, root_idx = t->n_node - 1;
	double *sib = kom_malloc(double, m);

	/* beta~(root,a) = q(a)/h_root; q(a)=t->p[root_idx][a] */
	for (a = 0; a < m; ++a) sd[root_idx].beta[a] = t->p[(size_t)root_idx * m * m + a] / sd[root_idx].h;

	/* Pre-order: root (highest index) down to leaves; parent always visited first */
	for (i = root_idx - 1; i >= 0; --i) {
		pc_node_t *u = t->node[i];
		pc_node_t *v = u->parent;
		int32_t k;
		double *beta_u = sd[i].beta, *beta_v = sd[v->ftime].beta;
		const double *pi = t->p + (size_t)i * m * m;

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
				s += pi[a*m + b] * beta_v[a] * sib[a];
			beta_u[b] = s / sd[i].h;
		}
	}
	free(sib);
}

/* Compute eta~[n_node*m*m] from inside/outside values already stored in sd.
 * For non-root u with parent v and sibling(s) w (except the root):
 *   eta~(u,b|a) = beta~(v,a) * alpha~(u,b) * prod_k alpha'~(sib_k,a) */
void pc_scfg_eta(const pc_tree_t *t, const pc_scfg_buf_t *sd, double *eta)
{
	int32_t i, k, a, b, m = t->m, root_idx = t->n_node - 1;
	double *sib = kom_malloc(double, m);

	for (i = 0; i < root_idx; ++i) { // eta~ is not defined at the root
		const pc_node_t *u = t->node[i], *v = u->parent;
		double *eta_i = eta + (size_t)i * m * m;
		const double *alpha_u = sd[i].alpha, *beta_v = sd[v->ftime].beta;
		for (a = 0; a < m; ++a) sib[a] = 1.0;
		for (k = 0; k < v->n_child; ++k)
			if (v->child[k] != u) {
				const double *a2k = sd[v->child[k]->ftime].alpha2;
				for (a = 0; a < m; ++a) sib[a] *= a2k[a];
			}
		for (a = 0; a < m; ++a)
			for (b = 0; b < m; ++b)
				eta_i[a*m + b] = beta_v[a] * alpha_u[b] * sib[a];
	}
	free(sib);
}

typedef struct {
	int32_t rotation; // 0: ((x,y)u,w)v; 1: ((w,y)u,x)v; 2: ((x,w)u,y)v
	int32_t u;
	double loglk;
	double p[];
} pc_nni_t;

/* Compute eta~ for all three NNI rotations at each eligible node u (internal,
 * non-root). eta has shape (n_node, 3, m, m); only entries for eligible nodes
 * are written. For the original topology ((x,y)u, w)v and the scaled quantity
 * q = beta~(v,a) / h(u):
 *   rotation 0 (original):        eta0(a,b) = q * alpha'~(w,a) * alpha'~(x,b) * alpha'~(y,b)
 *   rotation 1 ((w,y)u, x)v:      eta1(a,b) = q * alpha'~(x,a) * alpha'~(w,b) * alpha'~(y,b)
 *   rotation 2 ((x,w)u, y)v:      eta2(a,b) = q * alpha'~(y,a) * alpha'~(w,b) * alpha'~(x,b)
 * Nodes that are leaves or the root are skipped (left uninitialized). */
void pc_scfg_eta3_nni(const pc_tree_t *t, const pc_scfg_buf_t *sd, double *eta3)
{ // eta shape: (n,3,m,m)
	int32_t u, m = t->m;
	for (u = 0; u < t->n_node - 1; ++u) {
		const pc_node_t *up = t->node[u], *vp = up->parent;
		int32_t v, x, y, w, a, b; // original topology: ((x,y)u,w)v
		double *eta0_u, *eta1_u, *eta2_u;
		if (vp == 0 || up->n_child == 0) continue;
		assert(up->n_child == 2);
		v = vp->ftime;
		w = vp->child[(vp->child[0] == up)]->ftime;
		x = up->child[0]->ftime;
		y = up->child[1]->ftime;
		eta0_u = &eta3[(u * 3 + 0) * m * m];
		eta1_u = &eta3[(u * 3 + 1) * m * m];
		eta2_u = &eta3[(u * 3 + 2) * m * m];
		for (a = 0; a < m; ++a) {
			double s, q = sd[v].beta[a] / sd[u].h;
			for (b = 0, s = q * sd[w].alpha2[a]; b < m; ++b) // rotation 0
				eta0_u[a * m + b] = s * sd[x].alpha2[b] * sd[y].alpha2[b];
			for (b = 0, s = q * sd[x].alpha2[a]; b < m; ++b) // rotation 1: ((w,y)u,x)v
				eta1_u[a * m + b] = s * sd[w].alpha2[b] * sd[y].alpha2[b];
			for (b = 0, s = q * sd[y].alpha2[a]; b < m; ++b) // rotation 2: ((x,w)u,y)v
				eta2_u[a * m + b] = s * sd[w].alpha2[b] * sd[x].alpha2[b];
		}
	}
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

static void pc_scfg_cnt2p(pc_tree_t *t, const double *cnt, pc_model_t ct)
{
	int32_t u, m = t->m;
	double *tmp = kom_calloc(double, m * m);
	for (u = 0; u < t->n_node; ++u) {
		const double *cnt_u = cnt + (size_t)u * m * m;
		double *p_u = t->p + (size_t)u * m * m;
		pc_model_matrix(cnt_u, m, ct, tmp);
		pc_scfg_c2p(m, tmp, p_u);
	}
	free(tmp);
}

/* Run EM on the transition matrix of branch u under a given NNI rotation,
 * using the precomputed eta~[len][n_node*3*m*m] as sufficient statistics.
 * For each MSA column l, the branch log-likelihood contribution is
 *   log( sum_{a,b} p(b|a) * eta~(u,rotation,a,b)[l] )
 * and the posterior count is the normalized product p(b|a)*eta~(u,rotation,a,b).
 * Returns an allocated pc_nni_t with the optimized p[m*m] and loglk, or NULL
 * if u is a leaf or the root. Caller must free the result. */
pc_nni_t *pc_scfg_em_branch(const pc_tree_t *t, pc_model_t ct, int32_t len, double **eta3, int32_t u, int32_t rotation, int32_t max_itr)
{
	int32_t i, l, a, b, off, m = t->m, is3 = (rotation >= 0 && rotation < 3);
	pc_nni_t *q;
	const pc_node_t *up = t->node[u], *vp = up->parent;
	double *cnt, *tmp, *tmp2;
	if (vp == 0) return 0; // ignore root
	if (is3 && up->n_child == 0) return 0;
	assert(up->n_child == 2 || up->n_child == 0); // binary tree only
	q = (pc_nni_t*)calloc(1, sizeof(pc_nni_t) + sizeof(double) * m * m);
	q->u = u, q->rotation = rotation;
	memcpy(q->p, t->p + (size_t)u * m * m, sizeof(double) * m * m);
	cnt = kom_calloc(double, 3 * m * m);
	tmp = cnt + m * m;
	tmp2 = tmp + m * m;
	off = (is3? u * 3 + rotation : u) * m * m; // whether eta3 is calculated by _eta3_nni() or _eta()
	for (i = 0; i < max_itr; ++i) {
		double loglk = 0.0;
		memset(cnt, 0, sizeof(double) * m * m);
		for (l = 0; l < len; ++l) {
			double s = 0.0;
			for (a = 0; a < m; ++a)
				for (b = 0; b < m; ++b)
					s += (tmp[a * m + b] = q->p[a * m + b] * eta3[l][off + a * m + b]);
			loglk += log(s);
			for (a = 0, s = 1.0 / s; a < m; ++a)
				for (b = 0; b < m; ++b)
					cnt[a * m + b] += tmp[a * m + b] * s;
		}
		pc_model_matrix(cnt, m, ct, tmp2);
		pc_scfg_c2p(m, tmp2, q->p);
		q->loglk = loglk;
	}
	free(cnt);
	return q;
}

/* Compute posterior counts into cnt[n_node*m*m] (zeroed on entry) and return
 * the total log likelihood summed over all alignment columns. */
double pc_scfg_post_cnt(const pc_tree_t *t, const pc_msa_t *msa, pc_scfg_buf_t *sd, double *cnt)
{
	int32_t i, j, k, a, b, m = t->m;
	double loglk = 0.0;
	double *tmp = kom_malloc(double, m * m);
	double *sib = kom_malloc(double, m);
	memset(cnt, 0, (size_t)t->n_node * m * m * sizeof(double));
	for (i = 0; i < msa->len; ++i) {
		double sum, *cnt_j;
		loglk += pc_scfg_inside(t, msa, i, sd);
		pc_scfg_outside(t, sd);
		/* for each non-root branch, accumulate cnt[j][a*m+b] =
		 * p(b|a) * alpha~(u,b) * beta~(par,a) * prod_k alpha'~(sib_k,a) */
		for (j = 0; j < t->n_node - 1; ++j) {
			const pc_node_t *u = t->node[j], *v = u->parent;
			double *alpha_u = sd[j].alpha, *beta_v = sd[v->ftime].beta;
			const double *pj = t->p + (size_t)j * m * m;
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
					sum += (tmp[a*m + b] = pj[a*m + b] * alpha_u[b] * x);
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
	free(sib); free(tmp);
	return loglk;
}

double pc_scfg_em(pc_tree_t *t, const pc_msa_t *msa, pc_model_t ct, pc_scfg_buf_t *sd)
{
	int32_t m = t->m;
	double *cnt = kom_malloc(double, (size_t)t->n_node * m * m);
	double loglk = pc_scfg_post_cnt(t, msa, sd, cnt);
	pc_scfg_cnt2p(t, cnt, ct);
	free(cnt);
	return loglk;
}

/* Apply the best single NNI move. Uses caller-supplied p (already initialized).
 * Computes eta~ for all three rotations at each MSA column, then runs
 * pc_scfg_em_branch for every eligible node. The best (u, rotation) pair is
 * the one that maximizes loglk(rotation) - loglk(rotation 0). If any
 * improvement is found, applies pc_tree_rotate and rearranges p[] to match the
 * new post-order, replacing p[u] with the optimized branch matrix.
 * Returns the improvement (>=0); 0 means no move. */
double pc_scfg_nni(pc_tree_t *t, const pc_msa_t *msa, pc_model_t ct, int32_t max_iter_br)
{
	int32_t l, u, r, m = t->m, m2 = m * m, best_u = -1, best_r = 0;
	double **eta3, best_delta = 0.0;
	pc_scfg_buf_t *sd;
	pc_nni_t **nni;

	sd = pc_scfg_buf_new(t->n_node, m);

	eta3 = kom_malloc(double*, msa->len);
	eta3[0] = kom_malloc(double, (size_t)msa->len * t->n_node * 3 * m2);
	for (l = 1; l < msa->len; ++l)
		eta3[l] = eta3[l-1] + t->n_node * 3 * m2;
	for (l = 0; l < msa->len; ++l) {
		pc_scfg_inside(t, msa, l, sd);
		pc_scfg_outside(t, sd);
		pc_scfg_eta3_nni(t, sd, eta3[l]);
	}

	nni = kom_calloc(pc_nni_t*, 3 * t->n_node);
	for (u = 0; u < t->n_node; ++u) {
		nni[u*3+0] = pc_scfg_em_branch(t, ct, msa->len, eta3, u, 0, max_iter_br);
		nni[u*3+1] = pc_scfg_em_branch(t, ct, msa->len, eta3, u, 1, max_iter_br);
		nni[u*3+2] = pc_scfg_em_branch(t, ct, msa->len, eta3, u, 2, max_iter_br);
		if (nni[u*3+0] == NULL) continue;
		for (r = 1; r <= 2; ++r) {
			double delta = nni[u*3+r]->loglk - nni[u*3+0]->loglk;
			if (delta > best_delta)
				best_delta = delta, best_u = u, best_r = r;
		}
	}

	if (best_u >= 0) {
		/* Stamp tmp = old ftime on every node; after pc_tree_sync, tmp still
		 * holds the old ftime while ftime has the new value. */
		pc_node_t *up = t->node[best_u];
		int32_t xi = up->child[best_r - 1]->ftime; // child[0] for r=1, child[1] for r=2
		double *p_tmp = kom_malloc(double, (size_t)t->n_node * m * m);
		for (u = 0; u < t->n_node; ++u) t->node[u]->tmp = u;
		memcpy(t->p + (size_t)up->ftime * m2, nni[best_u*3+best_r]->p, sizeof(double) * m2); // update the matrix at u
		memcpy(p_tmp, t->p, sizeof(double) * (size_t)t->n_node * m * m);
		pc_tree_rotate(t, xi);
		for (u = 0; u < t->n_node; ++u) // u is new ftime; tmp is old ftime
			memcpy(t->p + (size_t)u * m2, p_tmp + (size_t)t->node[u]->tmp * m2, sizeof(double) * m2);
		free(p_tmp);
	}

	for (u = 0; u < 3 * t->n_node; ++u) free(nni[u]);
	free(nni); free(eta3[0]); free(eta3); free(sd);
	return best_delta;
}

void pc_scfg_model_cmp(const pc_tree_t *t, const pc_msa_t *msa, pc_model_t md0, pc_model_t md1, int32_t max_iter_br, double *diff)
{
	int32_t l, u, m = t->m, m2 = m * m;
	double **eta;
	pc_scfg_buf_t *sd;

	sd = pc_scfg_buf_new(t->n_node, m);
	eta = kom_malloc(double*, msa->len);
	eta[0] = kom_malloc(double, (size_t)msa->len * t->n_node * m2);
	for (l = 1; l < msa->len; ++l)
		eta[l] = eta[l-1] + t->n_node * m2;
	for (l = 0; l < msa->len; ++l) {
		pc_scfg_inside(t, msa, l, sd);
		pc_scfg_outside(t, sd);
		pc_scfg_eta(t, sd, eta[l]);
	}
	for (u = 0; u < t->n_node; ++u) {
		pc_nni_t *nni0 = pc_scfg_em_branch(t, md0, msa->len, eta, u, -1, max_iter_br);
		pc_nni_t *nni1 = pc_scfg_em_branch(t, md1, msa->len, eta, u, -1, max_iter_br);
		diff[u] = nni0 == NULL? 0.0 : nni0->loglk - nni1->loglk;
		free(nni0); free(nni1);
	}
	free(eta[0]); free(eta); free(sd);
}
