#include <stdlib.h>
#include <stdio.h>
#include "kommon.h"
#include "knhx.h"
#include "phycfg.h"
#include "khashl.h"
KHASHL_SET_INIT(KH_LOCAL, strset_t, strset, kh_cstr_t, kh_hash_str, kh_eq_str)

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
		p->ftime = i;
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
	while (top > 0) {
		p = &stack[top - 1];
		while (p->i != p->p->n_child) {
			kom_grow(fstack_t, stack, top, max);
			p = &stack[top++];
			p->i = 0, p->p = (p-1)->p->child[(p-1)->i];
		}
		if (node) node[ftime++] = (pc_node_t*)p->p;
		else ++ftime;
		if (--top > 0) {
			p = &stack[top - 1];
			p->i++;
		}
	}
	free(stack);
	return ftime;
}

void pc_tree_sync(pc_tree_t *t)
{
	int32_t i;
	free(t->node);
	t->n_node = pc_tree_expand(t->root, 0);
	t->node = kom_malloc(pc_node_t*, t->n_node);
	pc_tree_expand(t->root, t->node);
	for (i = 0; i < t->n_node; ++i)
		t->node[i]->ftime = i;
}

void pc_tree_destroy(pc_tree_t *t)
{
	int32_t i;
	if (t == NULL) return;
	for (i = 0; i < t->n_node; ++i) {
		free(t->node[i]->name);
		free(t->node[i]);
	}
	free(t->node);
	free(t);
}

void pc_tree_mark_leaf(pc_tree_t *t, int32_t n, char **leaf)
{
	int32_t i;
	strset_t *h;
	h = strset_init();
	for (i = 0; i < n; ++i) {
		int absent;
		strset_put(h, leaf[i], &absent);
	}
	for (i = 0; i < t->n_node; ++i) {
		pc_node_t *p = t->node[i];
		if (p->n_child == 0)
			p->aux = (strset_get(h, p->name) != kh_end(h)) ? 1 : 0;
	}
	strset_destroy(h);
}

pc_tree_t *pc_tree_reduce(pc_tree_t *t)
{
	int32_t i, n_marked = 0, max_child = 0, *a;
	pc_tree_t *s = 0;
	pc_node_t **sub, *p, *q;
	sub = kom_malloc(pc_node_t*, t->n_node);
	for (i = 0; i < t->n_node; ++i) {
		sub[i] = 0;
		max_child = max_child > t->node[i]->n_child? max_child : t->node[i]->n_child;
		if (t->node[i]->n_child == 0 && t->node[i]->aux > 0)
			++n_marked;
	}
	if (n_marked == 0) {
		free(sub);
		return 0;
	}
	a = kom_malloc(int32_t, max_child);
	for (i = 0; i < t->n_node; ++i) {
		p = t->node[i];
		if (p->n_child == 0) {
			if (p->aux <= 0) continue;
			q = pc_node_new(0);
			q->name = kom_strdup(p->name);
			q->d = p->d;
			sub[i] = q;
		} else {
			int32_t j, k;
			for (j = k = 0; j < p->n_child; ++j)
				if (sub[p->child[j]->ftime])
					a[k++] = p->child[j]->ftime;
			if (k == 1) { // just one node
				q = sub[i] = sub[a[0]];
				q->d += p->d;
			} else if (k >= 2) {
				q = pc_node_new(k);
				q->d = p->d;
				if (p->name) q->name = kom_strdup(p->name);
				for (j = 0; j < k; ++j) {
					sub[a[j]]->parent = q;
					q->child[j] = sub[a[j]];
				}
				sub[i] = q;
			} // if k == 0, do nothing
		}
	}
	q = sub[t->n_node - 1];
	free(a);
	free(sub);
	s = kom_calloc(pc_tree_t, 1);
	s->root = q;
	pc_tree_sync(s);
	return s;
}

/**********
 * Format *
 **********/

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
		if (p->name) kom_sprintf_lite(s, "%s", p->name);
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
