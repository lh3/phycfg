#ifndef PCPRIV_H
#define PCPRIV_H

#include "phycfg.h"

#ifdef __cplusplus
extern "C" {
#endif

// in scfg.c
double pc_scfg_inside(const pc_tree_t *t, const pc_msa_t *msa, int32_t pos, pc_scfg_buf_t *sd);
void pc_scfg_outside(const pc_tree_t *t, pc_scfg_buf_t *sd);
void pc_scfg_eta(const pc_tree_t *t, const pc_scfg_buf_t *sd, double *eta);
double pc_scfg_post_cnt(const pc_tree_t *t, const pc_msa_t *msa, pc_scfg_buf_t *sd, double *cnt);

// in model.c
void pc_model_matrix(const double *cnt, int32_t m, pc_model_t md, double *tmp);
double pc_model_dist_TN93(const double *cnt, double *kR, double *kY);

#ifdef __cplusplus
}
#endif

#endif
