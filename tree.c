#include <stdlib.h>
#include <stdio.h>
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
		p->name = kom_strdup(q->name);
		p->d = q->d;
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
	const pc_node_t *p;
	int32_t i;
} fstack_t;

int pc_tree_expand(const pc_node_t *root, pc_node_t **node)
{
	int ftime = 0, max = 256, top = 0;
	fstack_t *stack, *p;
	if (root == 0) return 0;
	stack = kom_malloc(fstack_t, max);
	p = &stack[top++];
	p->i = 0, p->p = root;
	while (top) { // p == stack[top-1]
		while (p->i != p->p->n_child) {
			kom_grow(fstack_t, stack, top, max);
			p = &stack[top++];
			p->i = 0, p->p = (p-1)->p->child[(p-1)->i];
		}
		if (node) node[ftime++] = (pc_node_t*)p->p;
		p = &stack[--top];
		if (top >= 0) p->i++;
	}
	free(stack);
	return ftime + 1;
}

static void format_node_recur(kstring_t *s, pc_node_t *p, char *numbuf)
{
	if (p->n_child > 0) {
		int i;
		kom_sprintf_lite(s, "(");
		for (i = 0; i < p->n_child; ++i) {
			if (i) kom_sprintf_lite(s, ",");
			format_node_recur(s, p->child[i], numbuf);
		}
		kom_sprintf_lite(s, ")");
		if (p->name) kom_sprintf_lite(s, "%s", p->name);
		if (p->d >= 0.0) {
			snprintf(numbuf, 127, ":%g", p->d);
			kom_sprintf_lite(s, "%s", numbuf);
		}
	} else {
		kom_sprintf_lite(s, "%s", p->name);
		if (p->d >= 0.0) {
			snprintf(numbuf, 127, ":%g", p->d);
			kom_sprintf_lite(s, "%s", numbuf);
		}
	}
}

void pc_tree_format(kstring_t *s, const pc_tree_t *t) // TODO: get rid of recursion
{
	char numbuf[128];
	format_node_recur(s, t->root, numbuf);
}
