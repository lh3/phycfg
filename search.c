#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "pcpriv.h"
#include "kavl-lite.h"
#include "kommon.h"

#define PC_NEG_INF (-1e100)

typedef struct pc_avln_s {
	double s, lk;
	int32_t no_nni;
	pc_node_t *p;
	KAVLL_HEAD(struct pc_avln_s) head;
	double p5[];
} pc_avln_t;

#define avlcmp(a, b) ((a)->s != (b)->s? ((a)->s < (b)->s) - ((b)->s < (a)->s) : ((b)->p < (a)->p) - ((a)->p < (b)->p)) // descending in s
KAVLL_INIT(pc_avl, pc_avln_t, head, avlcmp)

typedef struct {
	int32_t m, n_node, len;
	const int32_t *ucnt;
	double best_lk;
	uint64_t rng;
	pc_tree_t *best;
	pc_avln_t **avln, *root;
} pc_search_buf_t;

void pc_search_opt_init(pc_search_opt_t *opt)
{
	memset(opt, 0, sizeof(*opt));
	opt->md = PC_MD_FULL;
	opt->n_perturb_round = 5;
	opt->perturb_frac = 0.2;
	opt->eps = 0.01;
	opt->max_iter_br = 20;
	opt->max_iter_deep = 25;
}

pc_search_buf_t *pc_search_buf_alloc(pc_tree_t *t, const pc_msa_t *msa)
{
	int32_t u;
	pc_search_buf_t *sb;
	sb = kom_calloc(pc_search_buf_t, 1);
	sb->rng = 11;
	sb->m = t->m;
	sb->n_node = t->n_node;
	sb->len = msa->len_uniq;
	sb->best_lk = PC_NEG_INF;
	sb->ucnt = msa->ucnt;
	sb->avln = kom_calloc(pc_avln_t*, sb->n_node);
	for (u = 0; u < sb->n_node; ++u)
		sb->avln[u] = (pc_avln_t*)calloc(1, sizeof(pc_avln_t) + sizeof(double) * 5 * t->m * t->m);
	return sb;
}

void pc_search_buf_destroy(pc_search_buf_t *sb)
{
	int32_t u;
	for (u = 0; u < sb->n_node; ++u) free(sb->avln[u]);
	free(sb->avln); free(sb);
}

static pc_avln_t *pc_search_update_avl(pc_search_buf_t *sb, const pc_node_t *xp, pc_model_t md, double lk0, double eps, int32_t max_iter_br)
{
	pc_avln_t *xa = sb->avln[xp->ftime], *ret;
	assert(xp->n_child == 0 || xp->n_child == 2); // binary tree only
	ret = pc_avl_erase(&sb->root, xa);
	assert(ret == 0 || ret == xa);
	if (xp->parent && xp->parent->parent) {
		const pc_node_t *yp, *up, *wp, *vp;
		up = xp->parent, yp = up->child[(up->child[0] == xp)];
		vp = up->parent, wp = vp->child[(vp->child[0] == up)];
		xa->lk = pc_scfg_em5(sb->m, sb->len, sb->ucnt, md, wp, yp, up, xp, vp, max_iter_br, eps, xa->p5);
		xa->s = xa->lk - lk0;
		//fprintf(stderr, "YY\t%d\t%f\t%f -> %f\n", xp->ftime, xa->s, lk0, xa->lk);
	} else {
		xa->s = xa->lk = PC_NEG_INF;
	}
	pc_avl_insert(&sb->root, xa);
	return xa;
}

void pc_search_prepare(pc_tree_t *t, pc_search_buf_t *sb, pc_model_t md, double eps, int32_t max_iter_br)
{ // prerequisite: inside-outside applied
	int32_t u, x, l, a, m = sb->m;
	double lk0;
	for (u = 0; u < sb->n_node; ++u) {
		sb->avln[t->node[u]->ftime]->p = t->node[u];
		sb->avln[t->node[u]->ftime]->no_nni = -1;
	}
	for (l = 0, lk0 = 0.0; l < sb->len; ++l) {
		double s = 0.0;
		u = sb->n_node - 1; // the result should be the same regardless of the node
		for (a = 0; a < m; ++a)
			s += t->node[u]->q->alpha[l * m + a] * t->node[u]->q->beta[l * m + a];
		lk0 += log(s * t->node[u]->q->h[l]) * sb->ucnt[l];
	}
	sb->root = 0;
	for (x = 0, u = 0; x < sb->n_node; ++x)
		pc_search_update_avl(sb, t->node[x], md, lk0, eps, max_iter_br);
}

static void pc_search_update_tree(pc_search_buf_t *sb, const pc_avln_t *xa, pc_model_t md, double eps, int32_t max_iter_br)
{
	int32_t m2 = sb->m * sb->m;
	pc_node_t *xp = xa->p, *yp, *up, *wp, *vp;
	double lk0;

	assert(xp->parent != 0 && xp->parent->parent != 0);
	up = xp->parent, yp = up->child[(up->child[0] == xp)];
	vp = up->parent, wp = vp->child[(vp->child[0] == up)];

	memcpy(wp->q->p, &xa->p5[0 * m2], sizeof(double) * m2);
	memcpy(yp->q->p, &xa->p5[1 * m2], sizeof(double) * m2);
	memcpy(up->q->p, &xa->p5[2 * m2], sizeof(double) * m2);
	memcpy(xp->q->p, &xa->p5[3 * m2], sizeof(double) * m2);
	if (vp->parent)
		memcpy(vp->q->p, &xa->p5[4 * m2], sizeof(double) * m2);

	sb->avln[wp->ftime]->no_nni = xp->ftime;
	sb->avln[xp->ftime]->no_nni = wp->ftime;
	vp->child[vp->child[0] == wp? 0 : 1] = xp, xp->parent = vp;
	up->child[up->child[0] == xp? 0 : 1] = wp, wp->parent = up;
	pc_scfg_update5(sb->m, sb->len, sb->ucnt, up);

	lk0 = xa->lk;
	if (vp->parent && vp->parent->parent) {
		pc_search_update_avl(sb, vp, md, lk0, eps, max_iter_br);
	}
	pc_search_update_avl(sb, xp, md, lk0, eps, max_iter_br);
	pc_search_update_avl(sb, up, md, lk0, eps, max_iter_br);
	pc_search_update_avl(sb, wp, md, lk0, eps, max_iter_br);
	pc_search_update_avl(sb, yp, md, lk0, eps, max_iter_br);
	if (xp->n_child == 2) {
		pc_search_update_avl(sb, xp->child[0], md, lk0, eps, max_iter_br);
		pc_search_update_avl(sb, xp->child[1], md, lk0, eps, max_iter_br);
	}
	if (yp->n_child == 2) {
		pc_search_update_avl(sb, yp->child[0], md, lk0, eps, max_iter_br);
		pc_search_update_avl(sb, yp->child[1], md, lk0, eps, max_iter_br);
	}
	if (wp->n_child == 2) {
		pc_search_update_avl(sb, wp->child[0], md, lk0, eps, max_iter_br);
		pc_search_update_avl(sb, wp->child[1], md, lk0, eps, max_iter_br);
	}
}

int32_t pc_search_nni_greedy(pc_search_buf_t *sb, pc_model_t md, double eps, int32_t max_iter_br)
{
	int32_t n_nni = 0;
	while (1) {
		pc_avl_itr_t itr;
		const pc_avln_t *xa = 0;
		pc_avl_itr_first(sb->root, &itr);
		do {
			const pc_avln_t *xt = kavll_at(&itr);
			const pc_node_t *xp = xt->p, *up, *vp, *wp;
			if (xt->s < eps) break;
			up = xp->parent, vp = up->parent, wp = vp->child[(vp->child[0] == up)];
			if (wp->ftime != xt->no_nni) xa = xt;
			else break;
		} while (pc_avl_itr_next(&itr));
		if (xa == 0) break;
		++n_nni;
		if (kom_verbose >= 4) fprintf(stderr, "NI\t%d\t%f\n", xa->p->ftime, xa->s);
		pc_search_update_tree(sb, xa, md, eps, max_iter_br);
	}
	return n_nni;
}

int32_t pc_search_nni_perturb(pc_search_buf_t *sb, double perturb_frac, int32_t is_random)
{
	pc_avl_itr_t itr;
	pc_node_t **p;
	int32_t i, n_p = 0, n_perturb;
	assert(sb->n_node > 3);
	n_perturb = (int32_t)((sb->n_node - 3) * perturb_frac + 1.0);
	p = kom_calloc(pc_node_t*, n_perturb);
	pc_avl_itr_first(sb->root, &itr);
	do {
		const pc_avln_t *xa = kavll_at(&itr);
		pc_node_t *xp = xa->p;
		if (xp->parent == 0 || xp->parent->parent == 0) continue;
		if (!is_random || kom_u64todbl(kom_splitmix64(&sb->rng)) < perturb_frac)
			p[n_p++] = xp;
	} while (n_p < n_perturb && pc_avl_itr_next(&itr));
	for (i = 0; i < n_p; ++i) {
		pc_node_t *xp = p[i], *up, *vp, *wp;
		up = xp->parent, vp = up->parent, wp = vp->child[(vp->child[0] == up)];
		vp->child[vp->child[0] == wp? 0 : 1] = xp, xp->parent = vp;
		up->child[up->child[0] == xp? 0 : 1] = wp, wp->parent = up;
	}
	free(p);
	return n_p;
}

static double pc_search_recover(const pc_tree_t *src, const pc_msa_t *msa, pc_tree_t *t)
{
	int32_t u, m = src->m, m2 = m * m;
	pc_tree_copy(src, t);
	for (u = 0; u < t->n_node; ++u) {
		pc_node_t *up = t->node[u];
		pc_scfg_data_t *q0 = up->q;
		up->q = pc_scfg_data_new(m, msa->len_uniq);
		memcpy(up->q->p, q0->p, sizeof(double) * m2);
		free(q0);
	}
	return pc_scfg_post_cnt(t, msa);
}

static double pc_search_em_all(pc_tree_t *t, const pc_msa_t *msa, pc_model_t md, double eps, int32_t max_iter)
{
	int32_t k;
	double lk0, lk;
	for (k = 0, lk0 = PC_NEG_INF; k < max_iter; ++k) {
		lk = pc_scfg_em_all(t, msa, md);
		if (lk - lk0 < eps) break;
		lk0 = lk;
	}
	return lk;
}

void pc_search(pc_tree_t *t, const pc_msa_t *msa, const pc_search_opt_t *opt)
{
	int32_t k, l, is_random = 0;
	pc_search_buf_t *sb = pc_search_buf_alloc(t, msa);
	pc_search_em_all(t, msa, opt->md, opt->eps, opt->max_iter_deep);
	for (l = 0; l < opt->n_perturb_round + 1; ++l, is_random = !is_random) { // is_random is alternating
		double lk;
		if (l > 0) {
			pc_search_nni_perturb(sb, opt->perturb_frac, is_random);
			pc_tree_sync(t);
			lk = pc_search_em_all(t, msa, opt->md, opt->eps, opt->max_iter_br);
		}
		if (kom_verbose >= 4) fprintf(stderr, "TL\t%d\t%f\n", l, lk);
		for (k = 0; k < 5; ++k) {
			int32_t n_nni;
			if (kom_verbose >= 4) fprintf(stderr, "RD\t%d\t%d\n", l, k + 1);
			pc_search_prepare(t, sb, opt->md, opt->eps, opt->max_iter_br);
			n_nni = pc_search_nni_greedy(sb, opt->md, opt->eps, opt->max_iter_br);
			if (n_nni == 0) break;
			pc_tree_sync(t);
			lk = pc_search_em_all(t, msa, opt->md, opt->eps, opt->max_iter_br);
			if (kom_verbose >= 4) fprintf(stderr, "NL\t%d\t%d\t%f\n", l, k + 1, lk);
		}
		if (sb->best == 0 || sb->best_lk < lk) {
			if (sb->best) pc_tree_destroy(sb->best);
			sb->best = pc_tree_clone(t);
			sb->best_lk = lk;
		} else if (sb->best) {
			pc_search_recover(sb->best, msa, t);
		}
	}
	if (kom_verbose >= 4) fprintf(stderr, "BL\t%f\n", sb->best_lk);
	pc_search_buf_destroy(sb);
	if (msa->rt == PC_RT_NT || msa->rt == PC_RT_CODON)
		pc_model_dist(t, msa, PC_MD_TN93);
}
