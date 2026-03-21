#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "pcpriv.h"
#include "kavl-lite.h"
#include "kommon.h"

typedef struct pc_avlnode_s {
	double s;
	pc_node_t *p;
	KAVLL_HEAD(struct pc_avlnode_s) head;
} pc_avlnode_t;

#define avlcmp(a, b) ((a)->s != (b)->s? ((a)->s < (b)->s) - ((b)->s < (a)->s) : ((b)->p < (a)->p) - ((a)->p < (b)->p)) // descending in s
KAVLL_INIT(pc_avl, pc_avlnode_t, head, avlcmp)

struct pc_search_buf_s {
	int32_t m, n_node, len;
	pc_node_t **node;
	pc_avlnode_t *avln, *root;
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
	sb->avln = kom_calloc(pc_avlnode_t, sb->n_node);
	for (u = 0; u < sb->n_node; ++u)
		sb->node[u] = t->node[u];
	return sb;
}

void pc_search_buf_destroy(pc_search_buf_t *sb)
{
	free(sb->node); free(sb->avln); free(sb);
}

void pc_search_prepare(pc_search_buf_t *sb, pc_model_t md, double eps, int32_t max_iter_br)
{ // prerequisite: inside-outside applied
	int32_t u, x, l, a, m = sb->m, m2 = m * m;
	double *p5, lk0;
	p5 = kom_calloc(double, 5 * m2);
	for (u = 0; u < sb->n_node; ++u)
		sb->node[u]->q->flag = 0;
	for (l = 0, lk0 = 0.0; l < sb->len; ++l) {
		double s = 0.0;
		u = sb->n_node - 1; // the result should be the same regardless of the node
		for (a = 0; a < m; ++a)
			s += sb->node[u]->q->alpha[l * m + a] * sb->node[u]->q->beta[l * m + a];
		lk0 += log(s * sb->node[u]->q->h[l]);
	}
	for (x = 0, u = 0; x < sb->n_node - 1; ++x) {
		pc_node_t *xp = sb->node[x], *yp, *up, *wp, *vp;
		double lk;
		pc_avlnode_t *p;
		if (xp->parent == 0 || xp->parent->parent == 0) continue;
		up = xp->parent, yp = up->child[(up->child[0] == xp)];
		vp = up->parent, wp = vp->child[(vp->child[0] == up)];
		lk = pc_scfg_em5(m, sb->len, md, wp, yp, up, xp, vp, max_iter_br, eps, p5);
		p = &sb->avln[u++];
		p->s = lk - lk0;
		p->p = xp;
		pc_avl_insert(&sb->root, p);
	}
	free(p5);
}
