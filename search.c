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
	pc_node_t *p;
	KAVLL_HEAD(struct pc_avln_s) head;
	double p5[];
} pc_avln_t;

#define avlcmp(a, b) ((a)->s != (b)->s? ((a)->s < (b)->s) - ((b)->s < (a)->s) : ((b)->p < (a)->p) - ((a)->p < (b)->p)) // descending in s
KAVLL_INIT(pc_avl, pc_avln_t, head, avlcmp)

struct pc_search_buf_s {
	int32_t m, n_node, len;
	pc_node_t **node;
	pc_avln_t **avln, *root;
};

void pc_search_opt_init(pc_search_opt_t *opt)
{
	memset(opt, 0, sizeof(*opt));
	opt->eps_nni_init = 0.001;
}

pc_search_buf_t *pc_search_buf_init(pc_tree_t *t, int32_t len)
{
	int32_t u;
	pc_search_buf_t *sb;
	sb = kom_calloc(pc_search_buf_t, 1);
	sb->m = t->m;
	sb->n_node = t->n_node;
	sb->len = len;
	sb->node = kom_calloc(pc_node_t*, sb->n_node);
	sb->avln = kom_calloc(pc_avln_t*, sb->n_node);
	for (u = 0; u < sb->n_node; ++u) {
		sb->node[u] = t->node[u];
		sb->avln[u] = (pc_avln_t*)calloc(1, sizeof(pc_avln_t) + sizeof(double) * 5);
		sb->avln[t->node[u]->ftime]->p = t->node[u];
	}
	return sb;
}

void pc_search_buf_destroy(pc_search_buf_t *sb)
{
	int32_t u;
	for (u = 0; u < sb->n_node; ++u) free(sb->avln[u]);
	free(sb->node); free(sb->avln); free(sb);
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
		xa->lk = pc_scfg_em5(sb->m, sb->len, md, wp, yp, up, xp, vp, max_iter_br, eps, xa->p5);
		xa->s = xa->lk - lk0;
	} else {
		xa->s = xa->lk = PC_NEG_INF;
	}
	pc_avl_insert(&sb->root, xa);
	return xa;
}

static void pc_search_update_tree(pc_search_buf_t *sb, const pc_avln_t *xa, pc_model_t md, double eps, int32_t max_iter_br)
{
	int32_t m2 = sb->m * sb->m;
	pc_node_t *xp = xa->p, *yp, *up, *wp, *vp;

	assert(xp->parent != 0 && xp->parent->parent != 0);
	up = xp->parent, yp = up->child[(up->child[0] == xp)];
	vp = up->parent, wp = vp->child[(vp->child[0] == up)];

	memcpy(wp->q->p, &xa->p5[0 * m2], sizeof(double) * m2);
	memcpy(yp->q->p, &xa->p5[1 * m2], sizeof(double) * m2);
	memcpy(up->q->p, &xa->p5[2 * m2], sizeof(double) * m2);
	memcpy(xp->q->p, &xa->p5[3 * m2], sizeof(double) * m2);
	if (vp->parent)
		memcpy(vp->q->p, &xa->p5[4 * m2], sizeof(double) * m2);

	vp->child[vp->child[0] == wp? 0 : 1] = xp, xp->parent = vp;
	up->child[up->child[0] == xp? 0 : 1] = wp, wp->parent = up;

	if (xp->n_child == 2) {
		pc_search_update_avl(sb, xp->child[0], md, xa->lk, eps, max_iter_br);
		pc_search_update_avl(sb, xp->child[1], md, xa->lk, eps, max_iter_br);
	}
	if (yp->n_child == 2) {
		pc_search_update_avl(sb, yp->child[0], md, xa->lk, eps, max_iter_br);
		pc_search_update_avl(sb, yp->child[1], md, xa->lk, eps, max_iter_br);
	}
	if (wp->n_child == 2) {
		pc_search_update_avl(sb, wp->child[0], md, xa->lk, eps, max_iter_br);
		pc_search_update_avl(sb, wp->child[1], md, xa->lk, eps, max_iter_br);
	}
}

void pc_search_prepare(pc_search_buf_t *sb, pc_model_t md, double eps, int32_t max_iter_br)
{ // prerequisite: inside-outside applied
	int32_t u, x, l, a, m = sb->m;
	double lk0;
	for (u = 0; u < sb->n_node; ++u)
		sb->node[u]->q->flag = 0;
	for (l = 0, lk0 = 0.0; l < sb->len; ++l) {
		double s = 0.0;
		u = sb->n_node - 1; // the result should be the same regardless of the node
		for (a = 0; a < m; ++a)
			s += sb->node[u]->q->alpha[l * m + a] * sb->node[u]->q->beta[l * m + a];
		lk0 += log(s * sb->node[u]->q->h[l]);
	}
	for (x = 0, u = 0; x < sb->n_node; ++x)
		pc_search_update_avl(sb, sb->node[x], md, lk0, eps, max_iter_br);
}

void pc_search_nni_greedy(pc_search_buf_t *sb, pc_model_t md, double eps, int32_t max_iter_br)
{
	for (;;) {
	}
}
