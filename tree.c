#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "kommon.h"
#include "knhx.h"
#include "pcpriv.h"
#define kh_packed
#include "khashl.h"
KHASHL_MAP_INIT(KH_LOCAL, strmap_t, strmap, kh_cstr_t, int32_t, kh_hash_str, kh_eq_str)

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
		p->seq_id = -1;
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

void pc_tree_copy(const pc_tree_t *src, pc_tree_t *dst)
{
	int32_t i, k;
	pc_node_t **map;

	/* Free existing nodes in dst */
	if (dst->node) {
		for (i = 0; i < dst->n_node; ++i) {
			free(dst->node[i]->q);
			free(dst->node[i]->name);
			free(dst->node[i]);
		}
		free(dst->node);
	}

	dst->n_node = src->n_node;
	dst->m = src->m;
	map = dst->node = kom_malloc(pc_node_t*, src->n_node);

	/* Allocate and copy each node's scalar fields */
	for (i = 0; i < src->n_node; ++i) {
		const pc_node_t *s = src->node[i];
		pc_node_t *d = map[i] = pc_node_new(s->n_child);
		d->ftime  = s->ftime;
		d->seq_id = s->seq_id;
		d->d      = s->d;
		d->name   = s->name ? kom_strdup(s->name) : NULL;
		d->q      = s->q ? pc_scfg_data_new(src->m, -1) : NULL;
	}

	/* Wire up parent and child pointers using ftime as the map key */
	for (i = 0; i < src->n_node; ++i) {
		const pc_node_t *s = src->node[i];
		pc_node_t *d = map[i];
		d->parent = s->parent ? map[s->parent->ftime] : NULL;
		for (k = 0; k < s->n_child; ++k)
			d->child[k] = map[s->child[k]->ftime];
	}

	dst->root = map[src->root->ftime];
}

pc_tree_t *pc_tree_clone(const pc_tree_t *t)
{
	pc_tree_t *c = kom_calloc(pc_tree_t, 1);
	pc_tree_copy(t, c);
	return c;
}

void pc_tree_destroy(pc_tree_t *t)
{
	int32_t i;
	if (t == NULL) return;
	pc_scfg_free(t);
	for (i = 0; i < t->n_node; ++i) {
		free(t->node[i]->q);
		free(t->node[i]->name);
		free(t->node[i]);
	}
	free(t->node);
	free(t);
}

void pc_tree_strip_iname(pc_tree_t *t)
{
	int32_t i;
	for (i = 0; i < t->n_node; ++i) {
		pc_node_t *p = t->node[i];
		if (p->n_child > 0) {
			free(p->name);
			p->name = NULL;
		}
	}
}

void pc_tree_mark_leaf(const pc_tree_t *t, int32_t n, char **leaf, uint8_t *mark)
{
	int32_t i;
	strmap_t *h;
	h = strmap_init();
	for (i = 0; i < n; ++i) {
		int absent;
		strmap_put(h, leaf[i], &absent);
	}
	for (i = 0; i < t->n_node; ++i) {
		const pc_node_t *p = t->node[i];
		mark[i] = (p->n_child == 0 && p->name && strmap_get(h, p->name) != kh_end(h));
	}
	strmap_destroy(h);
}

pc_tree_t *pc_tree_reduce(const pc_tree_t *t, const uint8_t *mark)
{
	int32_t i, n_marked = 0, max_child = 0, *a;
	pc_tree_t *s = 0;
	pc_node_t **sub, *p, *q;
	sub = kom_calloc(pc_node_t*, t->n_node);
	for (i = 0; i < t->n_node; ++i) {
		max_child = max_child > t->node[i]->n_child? max_child : t->node[i]->n_child;
		if (t->node[i]->n_child == 0 && mark[i])
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
			if (!mark[i]) continue;
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
				if (p->d >= 0.0) {
					if (q->d < 0.0) q->d = p->d;
					else q->d += p->d;
				}
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
	if (q) q->parent = 0, q->d = -1.0;
	free(a);
	free(sub);
	s = kom_calloc(pc_tree_t, 1);
	s->root = q;
	pc_tree_sync(s);
	return s;
}

int32_t pc_tree_match_msa(pc_tree_t *t, const pc_msa_t *msa)
{
	int32_t i, n_unmatch = 0;
	strmap_t *h;
	h = strmap_init();
	for (i = 0; i < msa->n_seq; ++i) {
		int absent;
		khint_t k = strmap_put(h, msa->name[i], &absent);
		kh_val(h, k) = i;
	}
	for (i = 0; i < t->n_node; ++i) {
		pc_node_t *p = t->node[i];
		if (p->n_child == 0) {
			khint_t k = strmap_get(h, p->name);
			if (k != kh_end(h)) p->seq_id = kh_val(h, k);
			else { p->seq_id = -1; ++n_unmatch; }
		}
	}
	strmap_destroy(h);
	t->m = msa->m;
	return n_unmatch;
}

int32_t pc_tree_lca(const pc_tree_t *t, const uint8_t *mark)
{
	int32_t i, k, total = 0, lca = -1;
	int32_t *cnt;
	for (i = 0; i < t->n_node; ++i) total += mark[i];
	if (total == 0) return -1;
	cnt = kom_calloc(int32_t, t->n_node);
	for (i = 0; i < t->n_node && lca < 0; ++i) {
		const pc_node_t *v = t->node[i];
		if (v->n_child == 0) {
			cnt[i] = mark[i];
		} else {
			for (k = 0; k < v->n_child; ++k)
				cnt[i] += cnt[v->child[k]->ftime];
		}
		if (cnt[i] == total) lca = i;
	}
	free(cnt);
	return lca;
}

/* Rotate: given x, change ((x,y)u, w)v to ((w,y)u, x)v.
 * x's parent u must exist, u's parent v must exist, and v must have exactly 2
 * children. Returns 0 on success, -1 if not possible. Calls pc_tree_sync on
 * success to restore post-order. */
int32_t pc_tree_rotate(pc_tree_t *t, int32_t xi)
{
	pc_node_t *x, *u, *v, *w;
	int32_t k;
	if (xi < 0 || xi >= t->n_node) return -1;
	x = t->node[xi];
	u = x->parent;
	if (u == NULL) return -1; // x is root
	v = u->parent;
	if (v == NULL) return -1; // u is root
	if (v->n_child != 2) return -1;
	w = v->child[v->child[0] == u ? 1 : 0]; // sibling of u under v
	for (k = 0; k < u->n_child; ++k) // replace x with w in u->child[]
		if (u->child[k] == x) { u->child[k] = w; break; }
	for (k = 0; k < v->n_child; ++k) // replace w with x in v->child[]
		if (v->child[k] == w) { v->child[k] = x; break; }
	w->parent = u;
	x->parent = v;
	pc_tree_sync(t);
	return 0;
}

/**********
 * Reroot *
 **********/

/* Find the node p whose incoming branch contains the midpoint of the longest
 * leaf-to-leaf path (tree diameter). Returns p's ftime, or -1 if the tree
 * has fewer than two leaves. Writes the distance from p to the midpoint
 * (along p's incoming branch) to *dist_to_mid. */
int32_t pc_tree_mid_longest(const pc_tree_t *t, double *dist_to_mid)
{
	int32_t i, k, lca_idx = -1, b1_idx = -1, b2_idx = -1;
	double diam = -1.0;
	/* max_down[i]: max distance from node i to any leaf in its subtree;
	 * best_child[i]: child ftime on that deepest path (-1 for leaves) */
	double  *max_down   = kom_calloc(double,   t->n_node);
	int32_t *best_child = kom_malloc(int32_t,  t->n_node);
	for (i = 0; i < t->n_node; ++i) best_child[i] = -1;

	/* Single post-order pass: compute max_down and track the diameter. */
	for (i = 0; i < t->n_node; ++i) {
		const pc_node_t *v = t->node[i];
		double top1 = -1.0, top2 = -1.0;   /* two deepest child depths */
		int32_t c1 = -1, c2 = -1;
		for (k = 0; k < v->n_child; ++k) {
			const pc_node_t *ck = v->child[k];
			double d   = ck->d > 0.0 ? ck->d : 0.0;
			double dep = d + max_down[ck->ftime];
			if (dep >= top1) { top2=top1; c2=c1; top1=dep; c1=ck->ftime; }
			else if (dep > top2) { top2=dep; c2=ck->ftime; }
		}
		if (c1 >= 0) { max_down[i] = top1; best_child[i] = c1; }
		/* A path through v using its two deepest subtrees is a candidate diameter. */
		if (c1 >= 0 && c2 >= 0 && top1 + top2 > diam) {
			diam=top1+top2; lca_idx=i; b1_idx=c1; b2_idx=c2;
		}
	}

	if (lca_idx < 0) {   /* fewer than two leaves */
		free(max_down); free(best_child);
		if (dist_to_mid) *dist_to_mid = 0.0;
		return -1;
	}

	/* d1: distance from LCA to the leaf1 end of the diameter path */
	double d1   = (t->node[b1_idx]->d > 0.0 ? t->node[b1_idx]->d : 0.0) + max_down[b1_idx];
	double half = diam / 2.0;
	/* Descend from the LCA toward whichever leaf is farther from the midpoint. */
	int32_t curr   = (half <= d1) ? b1_idx : b2_idx;
	double  target = (half <= d1) ? (d1 - half) : (half - d1);
	double  accum  = 0.0, edge = 0.0;
	/* Walk down best_child[] until we reach the branch containing the midpoint. */
	while (1) {
		edge = t->node[curr]->d > 0.0 ? t->node[curr]->d : 0.0;
		if (accum + edge >= target) break;
		accum += edge;
		curr = best_child[curr];
	}
	/* The midpoint is (target - accum) from curr's parent; so it is
	 * (accum + edge - target) from curr itself along the same branch. */
	if (dist_to_mid) *dist_to_mid = accum + edge - target;
	free(max_down); free(best_child);
	return curr;
}

/* Place the root on the branch leading to p0 and return the new root. The
 * distance from p0 to the new root is dist */
pc_node_t *pc_tree_reroot_core(pc_node_t *root, pc_node_t *p0, double dist)
{
	int32_t i, k;
	double d, tmp;
	pc_node_t *p, *q, *r, *s, *new_root;

	if (p0 == root || p0->parent == NULL) return root; // no need to root
	if (dist > p0->d) return 0;
	if (dist < 0.0) dist = p0->d * 0.5;

	tmp = p0->d;  /* save original branch length of p0 */

	/* Create the new root with two children: p0 and p0's former parent */
	q = new_root = pc_node_new(2);
	new_root->seq_id = -1;
	new_root->d = -1.0;

	q->child[0] = p0;
	p0->d = dist;
	p = p0->parent;
	p0->parent = q;

	/* Find p0's slot in its former parent's child array */
	for (i = 0; i < p->n_child; ++i)
		if (p->child[i] == p0) break;

	q->child[1] = p;
	d = p->d;           /* save p's branch length before overwriting */
	p->d = tmp - dist;  /* remaining half of the split branch */
	r = p->parent;
	p->parent = q;

	/* Walk up to the old root reversing each parent->child edge.
	 * At each step: make r a child of p (reversing the edge),
	 * rotate branch lengths upward, then advance p, q, r. */
	while (r != NULL) {
		s = r->parent;       /* save r's parent before modifying */
		p->child[i] = r;     /* reverse edge: r becomes child of p */
		for (i = 0; i < r->n_child; ++i)  /* find p's slot in r */
			if (r->child[i] == p) break;
		r->parent = p;
		tmp = r->d; r->d = d; d = tmp;  /* rotate branch length upward */
		q = p; p = r; r = s;
	}

	/* p is now the old root; i is the slot of q in p's child array.
	 * Suppress the old root to maintain a proper rooted binary (or
	 * reduce arity by one for multifurcating roots). */
	if (p->n_child == 2) {
		/* Binary old root: bypass it, absorbing its branch into the
		 * remaining child so the total path length is preserved */
		r = p->child[1 - i];
		for (i = 0; i < q->n_child; ++i)  /* find p's slot in q */
			if (q->child[i] == p) break;
		r->d += p->d;
		r->parent = q;
		q->child[i] = r;
		free(p->name); free(p);
	} else {
		/* Multifurcating old root: just remove the reversed child */
		for (k = i; k < p->n_child - 1; ++k)
			p->child[k] = p->child[k + 1];
		--p->n_child;
	}

	return new_root;
}

void pc_tree_reroot(pc_tree_t *t, int32_t nid, double dist)
{
	pc_node_t *new_root;
	new_root = pc_tree_reroot_core(t->root, t->node[nid], dist);
	if (new_root == NULL) return;
	t->root = new_root;
	pc_tree_sync(t);
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

int32_t pc_tree_format(const pc_tree_t *t, char **s, int32_t *max)
{
	kstring_t str;
	char numbuf[128];
	str.s = *s; str.l = 0; str.m = (size_t)*max;
	format_node_recur(&str, t->root, numbuf);
	*s = str.s; *max = (int32_t)str.m;
	return (int32_t)str.l;
}
