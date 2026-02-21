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
} pc_tree_t;

#ifdef __cplusplus
extern "C" {
#endif

pc_tree_t *pc_tree_read(const char *fn);

#ifdef __cplusplus
}
#endif

#endif
