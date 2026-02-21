#ifndef PHYCFG_H
#define PHYCFG_H

#define PC_VERSION "0.0"

#include <stdint.h>
#include "knhx.h"

typedef struct {
	int n;      // numbder of nodes
	int m;      // capacity of a
	int err;    // error code
	knhx1_t *a; // list of nodes, of size n
} pc_itree_t;

typedef struct pc_node_s {
	int32_t n_child, aux;
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

pc_itree_t *pc_tree_read(const char *fn);
char **pc_list_read(const char *o, int *n_);

#ifdef __cplusplus
}
#endif

#endif
