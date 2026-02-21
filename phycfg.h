#ifndef PHYCFG_H
#define PHYCFG_H

#define PC_VERSION "0.0"

#include <stdint.h>
#include "kommon.h"

typedef struct pc_node_s {
	int32_t n_child, ftime;
	int64_t aux;
	double d;
	char *name;
	void *ptr;
	struct pc_node_s *parent;
	struct pc_node_s *child[];
} pc_node_t, *pc_node_p;

typedef struct {
	int32_t n_node;
	pc_node_t *root, **node;
} pc_tree_t;

#ifdef __cplusplus
extern "C" {
#endif

void pc_tree_format(kstring_t *s, const pc_tree_t *t);
pc_tree_t *pc_tree_parse(const char *str, char **en);
int pc_tree_expand(const pc_node_t *root, pc_node_t **node);
void pc_tree_sync(pc_tree_t *t);
void pc_tree_destroy(pc_tree_t *t);
void pc_tree_mark_leaf(pc_tree_t *t, int32_t n, char **leaf);
pc_tree_t *pc_tree_reduce(pc_tree_t *t);

pc_tree_t *pc_tree_read(const char *fn);
char **pc_list_read(const char *o, int *n_);

#ifdef __cplusplus
}
#endif

#endif
