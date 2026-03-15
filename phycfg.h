#ifndef PHYCFG_H
#define PHYCFG_H

#define PC_VERSION "0.0-dirty-r100"

#include <stdint.h>

typedef struct {
	double *p; // transition matrix (m,m)
	double *jc; // posterior joint count matrix (m,m)
	double *h; // scaling factor (l)
	double *alpha, *alpha2, *beta; // inside-outside functions (l,m) each
	double x[];
} pc_scfg_data_t;

typedef struct pc_node_s {
	int32_t n_child; // number of children
	int32_t ftime; // finish time
	int32_t seq_id; // sequence index in pc_msa_t; -1 for internal nodes or if missing
	int32_t tmp; // temporary field
	double d; // branch length
	pc_scfg_data_t *q;
	char *name; // node name
	struct pc_node_s *parent; // parent node; NULL for the root
	struct pc_node_s *child[]; // children; allocated along with the node
} pc_node_t, *pc_node_p;

typedef struct {
	int32_t n_node; // number of nodes
	int32_t m; // size of the alphabet
	pc_node_t *root; // the root node
	pc_node_t **node; // node[] in post order (i.e. descendants ahead of ancestors)
	double *p; // transition matrix of shape (n_node, m, m), if allocated
} pc_tree_t;

// residue type
typedef enum { PC_RT_UNKNOWN, PC_RT_NT, PC_RT_CODON, PC_RT_AA } pc_restype_t;

#define PC_GAP_NT 5
#define PC_GAP_AA 23

typedef struct {
	int32_t len; // alignment length
	int32_t n_seq; // number of sequences
	pc_restype_t rt; // residue type
	int32_t m; // size of the alphabet
	char **name; // sequence names
	uint8_t **msa; // the msa of shape (len, n_seq)
} pc_msa_t;

typedef enum { PC_MD_UNDEF = -1, PC_MD_FULL = 0, PC_MD_REV, PC_MD_TN93 } pc_model_t;

typedef struct {
	double h, *alpha, *alpha2, *beta;
} pc_scfg_buf_t;

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Parse tree in the NH format
 *
 * @param str      tree string
 * @param en       [out] end of the tree in _str_
 *
 * @return the tree, or NULL on error
 */
pc_tree_t *pc_tree_parse(const char *str, char **en);

void pc_tree_sync(pc_tree_t *t);
pc_tree_t *pc_tree_clone(const pc_tree_t *t);
void pc_tree_destroy(pc_tree_t *t);
void pc_tree_mark_leaf(const pc_tree_t *t, int32_t n, char **leaf, uint8_t *mark);
pc_tree_t *pc_tree_reduce(const pc_tree_t *t, const uint8_t *mark);
int32_t pc_tree_lca(const pc_tree_t *t, const uint8_t *mark);
int32_t pc_tree_rotate(pc_tree_t *t, int32_t xi);
int32_t pc_tree_mid_longest(const pc_tree_t *t, double *dist_to_mid);
void pc_tree_reroot(pc_tree_t *t, int32_t nid, double dist);
int32_t pc_tree_match_msa(pc_tree_t *t, const pc_msa_t *msa);

pc_tree_t *pc_tree_read(const char *fn);
int32_t pc_tree_format(const pc_tree_t *t, char **s, int32_t *max);

pc_msa_t *pc_msa_read(const char *fn);
void pc_msa_destroy(pc_msa_t *msa);
pc_restype_t pc_msa_infer_rt(const pc_msa_t *msa);
void pc_msa_encode(pc_msa_t *msa, pc_restype_t rt);
void pc_msa_filter(pc_msa_t *msa, int32_t min_cnt);
void pc_msa_select_codon(pc_msa_t *msa, int32_t codon_flag);

/**
 * Read a list of strings
 *
 * @param o        comma/space separated list, or a file if _o_ starts with '@'
 * @param n        [out] number of strings
 *
 * @return list of strings, allocated from malloc
 */
char **pc_list_read(const char *o, int *n_);

// Allocate buffer for SCFG calculation
pc_scfg_buf_t *pc_scfg_buf_new(int32_t n_node, int32_t m);

void pc_transmat_init(pc_tree_t *t);
double pc_scfg_em(pc_tree_t *t, const pc_msa_t *msa, pc_model_t ct, pc_scfg_buf_t *sd);
double pc_scfg_nni(pc_tree_t *t, const pc_msa_t *msa, pc_model_t ct, int32_t max_iter_br);

/**
 * Compare two models at each branch
 *
 * @param t           tree, ideally converged with the smaller model
 * @param msa         MSA
 * @param md0         1st model, typically the larger model
 * @param md1         2nd model, typically the smaller model
 * @param max_iter_br max EM rounds on each branch
 * @param diff        [out] log likelihood ratio log(P(md0)/P(md1)), of shape (n_node)
 */
void pc_scfg_model_cmp(const pc_tree_t *t, const pc_msa_t *msa, pc_model_t md0, pc_model_t md1, int32_t max_iter_br, double *diff);

// infer model from string; PC_MD_UNDEF if not defined
pc_model_t pc_model_from_str(const char *model_str);

// get the degree of freedom of a model
int32_t pc_model_df(pc_model_t md, int32_t m);

/**
 * Likelihood ratio test (typically on the output of pc_scfg_model_cmp)
 *
 * @param md0         first model
 * @param md1         second model of different degrees of freedom
 * @param m           sizeo fo the alphabet
 * @param lr          log likelihood ratio
 *
 * @return P-value under Wilks' theorem
 */
double pc_model_lrt(pc_model_t md0, pc_model_t md1, int32_t m, double lr);

// BIC difference
double pc_model_BIC(pc_model_t md0, pc_model_t md1, int32_t m, int32_t len, double lr);

// estimate branch lengths; only TN93 is supported for now
void pc_model_dist(pc_tree_t *t, const pc_msa_t *msa, pc_model_t md);

void pc_scfg_alloc(pc_tree_t *t, int32_t len);
void pc_scfg_free(pc_tree_t *t);
void pc_scfg_init_par(pc_tree_t *t);
double pc_scfg_em2(pc_tree_t *t, const pc_msa_t *msa, pc_model_t ct);
double pc_scfg_nni1(pc_tree_t *t, const pc_msa_t *msa, pc_model_t ct, int32_t max_iter_br);
double pc_scfg_nni4(pc_tree_t *t, const pc_msa_t *msa, pc_model_t ct, int32_t max_iter_br);

#ifdef __cplusplus
}
#endif

#endif
