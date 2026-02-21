#include <stdlib.h>
#include "kommon.h"
#include "knhx.h"
#include "phycfg.h"

pc_node_t *pc_node_new(int32_t n_child)
{
	pc_node_t *p;
	p = (pc_node_t*)calloc(1, sizeof(pc_node_t) + sizeof(void*) * n_child);
	p->n_child = n_child;
	return p;
}

pc_tree_t *pc_tree_parse(const char *str, char **en)
{
	int32_t i, k, n, max, err;
	knhx1_t *a;
	pc_tree_t *t = 0;
	a = kn_parse(str, &n, &max, &err, en);
	if (a == 0) return 0;
	if (err) goto end_parse; // TODO: check if there is a leak in knhx
	t = kom_calloc(pc_tree_t, 1);
	t->n_node = n;
	t->node = kom_calloc(pc_node_t*, n);
	for (i = 0; i < n; ++i) {
		knhx1_t *q = &a[i];
		pc_node_t *p;
		p = t->node[i] = pc_node_new(q->n);
		p->ftime = i, p->d = q->d;
		if (q->n == 0) t->n_leaf++;
	}
	for (i = n - 1; i >= 0; --i) {
		knhx1_t *q = &a[i];
		pc_node_t *p = t->node[i];
		for (k = 0; k < p->n_child; ++k) {
			p->child[k] = t->node[q->child[k]];
			p->child[k]->parent = p;
		}
	}
	t->root = t->node[t->n_node - 1];
end_parse:
	kn_destroy(n, a);
	return t;
}

typedef struct {
	pc_node_t *p;
	int32_t i;
} fstack_t;
