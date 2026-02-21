#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "knhx.h"

/*******************
 * kommon routines *
 *******************/

#define kom_malloc(type, cnt)       ((type*)malloc((cnt) * sizeof(type)))
#define kom_calloc(type, cnt)       ((type*)calloc((cnt), sizeof(type)))
#define kom_realloc(type, ptr, cnt) ((type*)realloc((ptr), (cnt) * sizeof(type)))
#define kom_grow(type, ptr, __i, __m) do { \
		if ((__i) >= (__m)) { \
			(__m) = (__i) + 1; \
			(__m) += ((__m)>>1) + 16; \
			(ptr) = kom_realloc(type, (ptr), (__m)); \
		} \
	} while (0)

static char *kom_strdup(const char *src) // strdup() doesn't conform to C99
{
	size_t len = strlen(src);
	char *dst = kom_malloc(char, len + 1);
	memcpy(dst, src, len + 1);
	return dst;
}

/*************
 * NH parser *
 *************/

typedef struct {
	int n, max, err;
	knhx1_t *a;
} kn_aux_t;

static inline char *add_node(const char *s, kn_aux_t *aux, int x)
{
	char *p, *nbeg, *nend = 0;
	knhx1_t *r;
	kom_grow(knhx1_t, aux->a, aux->n, aux->max);
	r = aux->a + (aux->n++);
	r->n = x; r->parent = -1;
	for (p = (char*)s, nbeg = p, r->d = -1.0; *p && *p != ',' && *p != ')'; ++p) {
		if (*p == '[') {
			if (nend == 0) nend = p;
			do { ++p; } while (*p && *p != ']');
			if (*p == 0) {
				aux->err |= KN_ERR_BRACKET;
				break;
			}
		} else if (*p == ':') {
			if (nend == 0) nend = p;
			r->d = strtod(p + 1, &p);
			--p;
		} else if (!isgraph(*p)) if (nend == 0) nend = p;
	}
	if (nend == 0) nend = p;
	if (nend != nbeg) {
		r->name = (char*)calloc(nend - nbeg + 1, 1);
		strncpy(r->name, nbeg, nend - nbeg);
	} else r->name = strdup("");
	return p;
}

knhx1_t *kn_parse(const char *nhx, int *_n, int *_max, int *_error, char **en)
{
	char *p = (char*)nhx;
	int *stack, top, max;
	kn_aux_t aux0, *aux = &aux0;

	#define __push_back(y) do { kom_grow(int, stack, top, max); stack[top++] = (y); } while (0)

	*_n = *_max = *_error = 0;
	if (en) *en = 0;
	memset(aux, 0, sizeof(*aux));
	while (*p && *p != '(') {} // search for the first (
	if (*p == 0) return 0; // no tree found

	stack = 0; top = max = 0;
	while (*p) {
		while (*p && !isgraph(*p)) ++p;
		if (*p == 0) break;
		if (*p == ',') ++p;
		else if (*p == '(') {
			__push_back(-1);
			++p;
		} else if (*p == ')') {
			int x = aux->n, m, i, fin = 0;
			for (i = top - 1; i >= 0; --i)
				if (stack[i] < 0) break;
			if (i < 0) {
				aux->err |= KN_ERR_NO_RGHT;
				break;
			} else if (i == 0) {
				if (en) *en = p + 1;
				fin = 1;
			}
			m = top - 1 - i;
			p = add_node(p + 1, aux, m);
			aux->a[x].child = (int*)calloc(m, sizeof(int));
			for (i = top - 1, m = m - 1; m >= 0; --m, --i) {
				aux->a[x].child[m] = stack[i];
				aux->a[stack[i]].parent = x;
			}
			top = i;
			__push_back(x);
			if (fin) break;
		} else {
			__push_back(aux->n);
			p = add_node(p, aux, 0);
		}
	}
	if (en && *p == 0) *en = p;
	*_n = aux->n, *_max = aux->max, *_error = aux->err;
	free(stack);
	return aux0.a;
}

/*************
 * Utilities *
 *************/

/****************
 * NH formatter *
 ****************/

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

static inline int kputsn(const char *p, int l, kstring_t *s)
{
	if (s->l + l + 1 >= s->m) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	memcpy(s->s + s->l, p, l);
	s->l += l; s->s[s->l] = 0;
	return l;
}

static inline int kputc(int c, kstring_t *s)
{
	if (s->l + 1 >= s->m) {
		s->m = s->l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	s->s[s->l++] = c; s->s[s->l] = 0;
	return c;
}

static void format_node_recur(const knhx1_t *node, const knhx1_t *p, kstring_t *s, char *numbuf)
{
	if (p->n) {
		int i;
		kputc('(', s);
		for (i = 0; i < p->n; ++i) {
			if (i) kputc(',', s);
			format_node_recur(node, &node[p->child[i]], s, numbuf);
		}
		kputc(')', s);
		if (p->name) kputsn(p->name, strlen(p->name), s);
		if (p->d >= 0) {
			snprintf(numbuf, 127, ":%g", p->d);
			kputsn(numbuf, strlen(numbuf), s);
		}
	} else {
	  kputsn(p->name, strlen(p->name), s);
	  if (p->d >= 0) {
	    snprintf(numbuf, 127, ":%g", p->d);
	    kputsn(numbuf, strlen(numbuf), s);
	  }
	}
}

void kn_format(const knhx1_t *node, int root, kstring_t *s) // TODO: get rid of recursion
{
	char numbuf[128];
	format_node_recur(node, &node[root], s, numbuf);
}

#ifdef KNHX_MAIN
int main(int argc, char *argv[])
{
	char *s = "((a[abc],d1)x:0.5,((b[&&NHX:S=MOUSE],h2)[&&NHX:S=HUMAN:B=99][blabla][&&NHX:K=foo],c))";
	knhx1_t *node;
	int i, j, n, max, error;
	kstring_t str;
	node = kn_parse(s, &n, &max, &error, 0);
	for (i = 0; i < n; ++i) {
		knhx1_t *p = node + i;
		printf("[%d] %s\t%d\t%d\t%g", i, p->name, p->parent, p->n, p->d);
		for (j = 0; j < p->n; ++j)
			printf("\t%d", p->child[j]);
		putchar('\n');
	}
	str.l = str.m = 0; str.s = 0;
	kn_format(node, n-1, &str);
	puts(str.s);
	free(str.s);
	return 0;
}
#endif
