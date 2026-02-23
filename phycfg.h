#ifndef PHYCFG_H
#define PHYCFG_H

#define PC_VERSION "0.0"

#include <stdint.h>

typedef struct pc_node_s {
	int32_t n_child, ftime;
	int32_t seq_id;
	double d;
	char *name;
	struct pc_node_s *parent;
	struct pc_node_s *child[];
} pc_node_t, *pc_node_p;

typedef struct {
	int32_t n_node;
	pc_node_t *root, **node;
} pc_tree_t;

typedef enum { PC_RT_UNKNOWN, PC_RT_NT, PC_RT_AA } pc_restype_t;

#define PC_GAP_NT 5
#define PC_GAP_AA 23

typedef struct {
	int32_t len, n_seq; // alignment length and number of sequences
	pc_restype_t rt; // residue type
	int32_t m; // size of the alphabet
	char **name; // sequence names
	uint8_t **msa; // len rows and n_seq columns
} pc_msa_t;

typedef struct {
	double h, *alpha, *alpha2, *beta; // pointers point to x[]
} pc_scfg_t;

#ifdef __cplusplus
extern "C" {
#endif

pc_tree_t *pc_tree_parse(const char *str, char **en);
int pc_tree_expand(const pc_node_t *root, pc_node_t **node);
void pc_tree_sync(pc_tree_t *t);
pc_tree_t *pc_tree_clone(const pc_tree_t *t);
void pc_tree_destroy(pc_tree_t *t);
void pc_tree_mark_leaf(const pc_tree_t *t, int32_t n, char **leaf, uint8_t *mark);
pc_tree_t *pc_tree_reduce(const pc_tree_t *t, const uint8_t *mark);
int32_t pc_tree_lca(const pc_tree_t *t, const uint8_t *mark);
int32_t pc_tree_mid_longest(const pc_tree_t *t, double *dist_to_mid);
pc_node_t *pc_tree_reroot_core(pc_node_t *root, pc_node_t *p0, double dist);
void pc_tree_reroot(pc_tree_t *t, int32_t nid, double dist);
int32_t pc_tree_match_msa(pc_tree_t *t, const pc_msa_t *msa);

pc_tree_t *pc_tree_read(const char *fn);
int32_t pc_tree_format(const pc_tree_t *t, char **s, int32_t *max);

pc_msa_t *pc_msa_read(const char *fn);
void pc_msa_destroy(pc_msa_t *msa);
pc_restype_t pc_msa_infer_rt(const pc_msa_t *msa);
void pc_msa_encode(pc_msa_t *msa, pc_restype_t rt);
void pc_msa_filter(pc_msa_t *msa, int32_t min_cnt, int32_t is_cds);

char **pc_list_read(const char *o, int *n_);

double pc_scfg_inside(const pc_tree_t *t, double **tm, const pc_msa_t *msa, int32_t pos, pc_scfg_t *sd);
void pc_scfg_outside(const pc_tree_t *t, double **tm, int32_t m, pc_scfg_t *sd);

#ifdef __cplusplus
}
#endif

#endif
