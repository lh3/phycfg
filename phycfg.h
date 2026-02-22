#ifndef PHYCFG_H
#define PHYCFG_H

#define PC_VERSION "0.0"

#include <stdint.h>

typedef struct {
	double *p, *alpha, *alpha2, *beta, h;
} pc_scfg_t;

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

typedef enum { PC_RT_UNKNOWN, PC_RT_NT, PC_RT_AA } pc_restype_t;

typedef struct {
	int32_t n_pos, n_seq; // number of positions and sequences
	pc_restype_t rt; // residue type
	int32_t m; // size of the alphabet
	char **name; // sequence names
	uint8_t **msa; // n_pos rows and n_seq columns
} pc_msa_t;

#ifdef __cplusplus
extern "C" {
#endif

pc_tree_t *pc_tree_parse(const char *str, char **en);
int pc_tree_expand(const pc_node_t *root, pc_node_t **node);
void pc_tree_sync(pc_tree_t *t);
void pc_tree_destroy(pc_tree_t *t);
void pc_tree_mark_leaf(pc_tree_t *t, int32_t n, char **leaf);
pc_tree_t *pc_tree_reduce(pc_tree_t *t);

pc_tree_t *pc_tree_read(const char *fn);
int32_t pc_tree_format(const pc_tree_t *t, char **s, int32_t *max);

pc_msa_t *pc_msa_read(const char *fn);
pc_restype_t pc_msa_infer_rt(const pc_msa_t *msa);
void pc_msa_encode(pc_msa_t *msa, pc_restype_t rt);

char **pc_list_read(const char *o, int *n_);

#ifdef __cplusplus
}
#endif

#endif
