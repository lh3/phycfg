#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
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

/****************************
 * Inside-outside algorithm *
 ****************************/

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

double pc_scfg_inside(pc_tree_t *t, const pc_msa_t *msa, int32_t pos)
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

void pc_scfg_outside(pc_tree_t *t, int32_t pos)
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

/*****************
 * All-branch EM *
 *****************/

static inline double pc_scfg_add_jc(int32_t m, const double *ua, const double *w2, const double *vb, const double *uq, double *uc, double *tmp)
{
	int32_t a, b;
	double s, t;
	for (a = 0, s = 0.0; a < m; ++a) {
		double t = vb[a] * w2[a];
		for (b = 0; b < m; ++b)
			s += (tmp[a * m + b] = uq[a * m + b] * t * ua[b]); // \eta_u(b|a) = ua[b] * t
	}
	for (a = 0, t = 1.0 / s; a < m; ++a)
		for (b = 0; b < m; ++b)
			uc[a * m + b] += tmp[a * m + b] * t;
	return s; // \sum_{a,b} p_u(b|a)\eta_u(b|a)
}

// fill pc_node_t::q->jc (posterior joint count)
double pc_scfg_post_cnt(pc_tree_t *t, const pc_msa_t *msa)
{
	int32_t i, j, k, a, b, m = t->m, m2 = m * m;
	double loglk = 0.0;
	double *tmp = kom_malloc(double, m2 + m), *sib = tmp + m2;
	for (j = 0; j < t->n_node; ++j)
		memset(t->node[j]->q->jc, 0, m2 * sizeof(double));
	for (i = 0; i < msa->len; ++i) {
		loglk += pc_scfg_inside(t, msa, i);
		pc_scfg_outside(t, i);
		for (j = 0; j < t->n_node - 1; ++j) {
			pc_node_t *u = t->node[j], *v = u->parent;
			for (a = 0; a < m; ++a) sib[a] = 1.0;
			for (k = 0; k < v->n_child; ++k) {
				if (v->child[k] != u) {
					const double *a2k = &v->child[k]->q->alpha2[i * m];
					for (a = 0; a < m; ++a) sib[a] *= a2k[a];
				}
			}
			pc_scfg_add_jc(m, &u->q->alpha[i * m], sib, &v->q->beta[i * m], u->q->p, u->q->jc, tmp);
		}
		// special case: the root
		assert(j == t->n_node - 1);
		{
			const pc_node_t *r = t->node[j];
			const double *alpha_r = &r->q->alpha[i * m], *beta_r = &r->q->beta[i * m];
			double sum, *jc = r->q->jc;
			for (a = 0, sum = 0.0; a < m; ++a)
				sum += alpha_r[a] * beta_r[a];
			sum = 1.0 / sum;
			for (a = 0; a < m; ++a)
				for (b = 0; b < m; ++b)
					jc[a*m + b] += alpha_r[b] * beta_r[b] * sum;
		}
	}
	free(tmp);
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

double pc_scfg_em_all(pc_tree_t *t, const pc_msa_t *msa, pc_model_t ct)
{
	double loglk = pc_scfg_post_cnt(t, msa);
	pc_scfg_cnt2p(t, ct);
	return loglk;
}

/***************
 * 1-branch EM *
 ***************/

double pc_scfg_em1(int32_t m, int32_t len, pc_model_t ct, const pc_node_t *xp, const pc_node_t *yp, const pc_node_t *up, const pc_node_t *wp, const pc_node_t *vp, int32_t max_itr, double *p)
{ // ((x,y)u,w)v
	const double eps = 1e-6;
	int32_t i, l, a;
	double *cnt, *tmp, *ua, s, loglk = 0.0, loglk0 = 0.0;
	memcpy(p, up->q->p, sizeof(double) * m * m);
	cnt = kom_calloc(double, 3 * m * m);
	tmp = cnt + m * m;
	ua = tmp + m * m;
	for (i = 0; i < max_itr; ++i) {
		memset(cnt, 0, sizeof(double) * m * m);
		for (l = 0, loglk = 0.0; l < len; ++l) {
			if (xp && yp) { // with NNI, ua is not &up->q->alpha[l * m]
				const double *a2x = &xp->q->alpha2[l * m];
				const double *a2y = &yp->q->alpha2[l * m];
				double hu1 = 1.0 / up->q->h[l];
				for (a = 0; a < m; ++a)
					ua[a] = a2x[a] * a2y[a] * hu1;
			} else {
				memcpy(ua, &up->q->alpha[l * m], sizeof(double) * m);
			}
			s = pc_scfg_add_jc(m, ua, &wp->q->alpha2[l * m], &vp->q->beta[l * m], p, cnt, tmp);
			loglk += log(s);
		}
		pc_model_matrix(cnt, m, ct, tmp);
		pc_scfg_c2p(m, tmp, p);
		if (i >= 3 && loglk > loglk0 && loglk - loglk0 < eps)
			break;
		loglk0 = loglk;
		//fprintf(stderr, "YY\t((%d,%d)%d,%d)\t%d\t%f\n", xp->ftime, yp->ftime, up->ftime, wp->ftime, i, loglk);
	}
	free(cnt);
	return loglk;
}

double pc_scfg_nni1(pc_tree_t *t, const pc_msa_t *msa, pc_model_t ct, int32_t max_iter_br)
{
	int32_t l, u, r, m = t->m, m2 = m * m, best_u = -1, best_r = 0;
	double best_delta = 0.0, *best_p, *p;
	best_p = kom_calloc(double, m2 * 4);
	p = best_p + m2;
	for (l = 0; l < msa->len; ++l) {
		pc_scfg_inside(t, msa, l);
		pc_scfg_outside(t, l);
	}
	for (u = 0; u < t->n_node - 1; ++u) { // skip the root
		const pc_node_t *up = t->node[u], *vp = up->parent, *wp, *xp, *yp;
		double loglk[3];
		if (up->n_child != 2 || vp->n_child != 2) continue;
		wp = vp->child[(vp->child[0] == up)];
		xp = up->child[0], yp = up->child[1];
		loglk[0] = pc_scfg_em1(m, msa->len, ct, xp, yp, up, wp, vp, max_iter_br, &p[0 * m2]);
		loglk[1] = pc_scfg_em1(m, msa->len, ct, wp, yp, up, xp, vp, max_iter_br, &p[1 * m2]);
		loglk[2] = pc_scfg_em1(m, msa->len, ct, xp, wp, up, yp, vp, max_iter_br, &p[2 * m2]);
		for (r = 1; r <= 2; ++r) {
			double delta = loglk[r] - loglk[0];
			if (delta > best_delta) {
				best_delta = delta, best_u = u, best_r = r;
				memcpy(best_p, &p[r * m2], m2 * sizeof(double));
			}
		}
	}
	if (best_u >= 0) {
		pc_node_t *up = t->node[best_u];
		memcpy(up->q->p, best_p, sizeof(double) * m2);
		pc_tree_rotate(t, up->child[best_r - 1]->ftime);
	}
	free(best_p);
	return best_delta;
}

void pc_scfg_model_cmp(pc_tree_t *t, const pc_msa_t *msa, pc_model_t md0, pc_model_t md1, int32_t max_iter_br, double *diff)
{
	int32_t l, u, m = t->m, m2 = m * m;
	double *p = kom_calloc(double, m2);
	for (l = 0; l < msa->len; ++l) {
		pc_scfg_inside(t, msa, l);
		pc_scfg_outside(t, l);
	}
	for (u = 0; u < t->n_node - 1; ++u) {
		const pc_node_t *up = t->node[u], *vp = up->parent, *wp = vp->child[(vp->child[0] == up)];
		double lk0 = pc_scfg_em1(m, msa->len, md0, 0, 0, up, wp, vp, max_iter_br, p);
		double lk1 = pc_scfg_em1(m, msa->len, md1, 0, 0, up, wp, vp, max_iter_br, p);
		diff[u] = lk0 - lk1;
	}
	diff[u] = 0.0; // root
	free(p);
}

/***************
 * 5-branch EM *
 ***************/

double pc_scfg_em5(int32_t m, int32_t len, pc_model_t ct, const pc_node_t *xp, const pc_node_t *yp, const pc_node_t *up, const pc_node_t *wp, const pc_node_t *vp, int32_t max_itr, double *q)
{ // topology: (((x,y)u,w)v,z)p
	const double eps = 1e-6;
	int32_t i, l, a, b, m2 = m * m;
	double *tmp, loglk = 0.0, loglk0 = 0.0;
	double *xq = q, *yq = xq + m2, *uq = yq + m2, *wq = uq + m2, *vq = wq + m2;
	double *ua, *ub, *va, *vb, *x2, *y2, *u2, *w2, *xc, *yc, *uc, *wc, *vc;
	memcpy(xq, xp->q->p, sizeof(double) * m2);
	memcpy(yq, yp->q->p, sizeof(double) * m2);
	memcpy(uq, up->q->p, sizeof(double) * m2);
	memcpy(wq, wp->q->p, sizeof(double) * m2);
	memcpy(vq, vp->q->p, sizeof(double) * m2);
	tmp = kom_calloc(double, 6 * m2 + 8 * m);
	xc = tmp + m2, yc = xc + m2, uc = yc + m2, wc = uc + m2, vc = wc + m2;
	ua = vc + m2, ub = ua + m, va = ub + m, vb = va + m;
	x2 = vb + m, y2 = x2 + m, u2 = y2 + m, w2 = u2 + m;
	for (i = 0; i < max_itr; ++i) {
		memset(xc, 0, sizeof(double) * m2 * 5); // clear all counts as they are allocated together
		for (l = 0, loglk = 0.0; l < len; ++l) {
			const double *xa = &xp->q->alpha[l * m], *ya = &yp->q->alpha[l * m], *wa = &wp->q->alpha[l * m]; // not affected
			double s, hu1 = 1.0 / up->q->h[l]; // 1/h(u), the scaling factor at u
			for (a = 0; a < m; ++a) // alpha2~(x,a)
				for (b = 0, x2[a] = 0.0; b < m; ++b)
					x2[a] += xq[a * m + b] * xa[b];
			for (a = 0; a < m; ++a) // alpha2~(y,a)
				for (b = 0, y2[a] = 0.0; b < m; ++b)
					y2[a] += yq[a * m + b] * ya[b];
			for (a = 0; a < m; ++a) // alpha2~(w,a)
				for (b = 0, w2[a] = 0.0; b < m; ++b)
					w2[a] += wq[a * m + b] * wa[b];
			for (a = 0; a < m; ++a) // alpha~(u,a)
				ua[a] = x2[a] * y2[a] * hu1;
			for (a = 0; a < m; ++a) // alpha2~(u,a)
				for (b = 0, u2[a] = 0.0; b < m; ++b)
					u2[a] += uq[a * m + b] * ua[b];
			if (vp->parent != 0) { // if vp is not the root
				const pc_node_t *pp = vp->parent, *zp = pp->child[(pp->child[0] == vp)];
				const double hv1 = 1.0 / vp->q->h[l], *pb = &pp->q->beta[l * m], *z2 = &zp->q->alpha2[l * m];
				for (a = 0; a < m; ++a) // alpha~(v,a)
					va[a] = u2[a] * w2[a] * hv1;
				for (b = 0; b < m; ++b) { // beta~(v,b)
					for (a = 0, vb[b] = 0.0; a < m; ++a)
						vb[b] += pb[a] * z2[a] * vq[a * m + b];
					vb[b] *= hv1;
				}
			} else {
				vb = &vp->q->beta[l * m];
			}
			for (b = 0, s = 0.0; b < m; ++b) { // beta~(u,b)
				for (a = 0, ub[b] = 0.0; a < m; ++a)
					ub[b] += vb[a] * w2[a] * uq[a * m + b];
				s += ua[b] * ub[b];
				ub[b] *= hu1;
			}
			loglk += log(s);
			pc_scfg_add_jc(m, xa, y2, ub, xq, xc, tmp);
			pc_scfg_add_jc(m, ya, x2, ub, yq, yc, tmp);
			pc_scfg_add_jc(m, ua, w2, vb, uq, uc, tmp);
			pc_scfg_add_jc(m, wa, u2, vb, wq, wc, tmp);
			if (vp->parent != 0) {
				const pc_node_t *pp = vp->parent, *zp = pp->child[(pp->child[0] == vp)];
				pc_scfg_add_jc(m, va, &zp->q->alpha2[l * m], &pp->q->beta[l * m], vq, vc, tmp);
			}
		}
		pc_model_matrix(xc, m, ct, tmp); pc_scfg_c2p(m, tmp, xq);
		pc_model_matrix(yc, m, ct, tmp); pc_scfg_c2p(m, tmp, yq);
		pc_model_matrix(uc, m, ct, tmp); pc_scfg_c2p(m, tmp, uq);
		pc_model_matrix(wc, m, ct, tmp); pc_scfg_c2p(m, tmp, wq);
		if (vp->parent != 0) {
			pc_model_matrix(vc, m, ct, tmp); pc_scfg_c2p(m, tmp, vq);
		}
		if (i >= 3 && loglk > loglk0 && loglk - loglk0 < eps)
			break;
		loglk0 = loglk;
		//fprintf(stderr, "YY\t((%d,%d)%d,%d)\t%d\t%f\n", xp->ftime, yp->ftime, up->ftime, wp->ftime, i, loglk);
	}
	free(tmp);
	return loglk;
}

double pc_scfg_nni5(pc_tree_t *t, const pc_msa_t *msa, pc_model_t ct, int32_t max_iter_br)
{
	int32_t l, u, r, m = t->m, m2 = m * m, best_u = -1, best_r = 0;
	double best_delta = 0.0, *best_p, *p;
	best_p = kom_calloc(double, 4 * m2 * 5);
	p = best_p + m2 * 5;
	for (l = 0; l < msa->len; ++l) {
		pc_scfg_inside(t, msa, l);
		pc_scfg_outside(t, l);
	}
	for (u = 0; u < t->n_node - 1; ++u) { // skip the root
		const pc_node_t *up = t->node[u], *vp = up->parent, *wp, *xp, *yp;
		double loglk[3];
		if (up->n_child != 2 || vp->n_child != 2) continue;
		wp = vp->child[(vp->child[0] == up)];
		xp = up->child[0], yp = up->child[1];
		loglk[0] = pc_scfg_em5(m, msa->len, ct, xp, yp, up, wp, vp, max_iter_br, &p[0 * m2 * 5]);
		loglk[1] = pc_scfg_em5(m, msa->len, ct, wp, yp, up, xp, vp, max_iter_br, &p[1 * m2 * 5]);
		loglk[2] = pc_scfg_em5(m, msa->len, ct, xp, wp, up, yp, vp, max_iter_br, &p[2 * m2 * 5]);
		for (r = 1; r <= 2; ++r) {
			double delta = loglk[r] - loglk[0];
			if (delta > best_delta) {
				best_delta = delta, best_u = u, best_r = r;
				memcpy(best_p, &p[r * m2 * 5], 5 * m2 * sizeof(double));
			}
		}
	}
	if (best_u >= 0) {
		pc_node_t *up = t->node[best_u], *vp = up->parent, *wp = vp->child[(vp->child[0] == up)];
		pc_node_t *xp = up->child[0], *yp = up->child[1];
		memcpy(up->q->p, &best_p[2 * m2], sizeof(double) * m2);
		if (vp->parent != 0)
			memcpy(vp->q->p, &best_p[4 * m2], sizeof(double) * m2);
		if (best_r == 1) {
			memcpy(wp->q->p, &best_p[0 * m2], sizeof(double) * m2);
			memcpy(yp->q->p, &best_p[1 * m2], sizeof(double) * m2);
			memcpy(xp->q->p, &best_p[3 * m2], sizeof(double) * m2);
		} else {
			memcpy(xp->q->p, &best_p[0 * m2], sizeof(double) * m2);
			memcpy(wp->q->p, &best_p[1 * m2], sizeof(double) * m2);
			memcpy(yp->q->p, &best_p[3 * m2], sizeof(double) * m2);
		}
		pc_tree_rotate(t, up->child[best_r - 1]->ftime);
	}
	free(best_p);
	return best_delta;
}
