#ifndef PCPRIV_H
#define PCPRIV_H

#include "phycfg.h"

#ifdef __cplusplus
extern "C" {
#endif

// in scfg.c
void pc_scfg_free(pc_tree_t *t);
double pc_scfg_post_cnt(pc_tree_t *t, const pc_msa_t *msa);

// in sfunc.c
double kf_chi2_p(double df, double x);

// in model.c
void pc_model_matrix(const double *cnt, int32_t m, pc_model_t md, double *tmp);
double pc_model_dist_TN93(const double *cnt, double *kR, double *kY);

#ifdef __cplusplus
}
#endif

#endif
