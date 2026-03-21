#ifndef PCPRIV_H
#define PCPRIV_H

#include "phycfg.h"

#ifdef __cplusplus
extern "C" {
#endif

// in msa.c
void pc_msa_uniq(pc_msa_t *msa);

// in scfg.c
void pc_scfg_free(pc_tree_t *t);
double pc_scfg_post_cnt(pc_tree_t *t, const pc_msa_t *msa);
double pc_scfg_em1(int32_t m, int32_t len, pc_model_t ct, const pc_node_t *xp, const pc_node_t *yp, const pc_node_t *up, const pc_node_t *wp, const pc_node_t *vp, int32_t max_itr, double eps, double *p);
double pc_scfg_em5(int32_t m, int32_t len, pc_model_t ct, const pc_node_t *xp, const pc_node_t *yp, const pc_node_t *up, const pc_node_t *wp, const pc_node_t *vp, int32_t max_itr, double eps, double *q);
double pc_scfg_update5(int32_t m, int32_t len, pc_node_t *up);

// in sfunc.c
double kf_chi2_p(double df, double x);

// in model.c
void pc_model_matrix(const double *cnt, int32_t m, pc_model_t md, double *tmp);
double pc_model_dist_TN93(const double *cnt, double *kR, double *kY);

#ifdef __cplusplus
}
#endif

#endif
