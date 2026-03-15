#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "pcpriv.h"
#include "kommon.h"

pc_model_t pc_model_from_str(const char *model)
{
	if (strcmp(model, "full") == 0 || strcmp(model, "FULL") == 0 || strcmp(model, ".") == 0) return PC_MD_FULL;
	if (strcmp(model, "rev") == 0 || strcmp(model, "GTR") == 0 || strcmp(model, "gtr") == 0) return PC_MD_REV;
	if (strcmp(model, "TN93") == 0 || strcmp(model, "tn93") == 0) return PC_MD_TN93;
	return PC_MD_UNDEF;
}

int32_t pc_model_df(pc_model_t md, int32_t m)
{
	if (md == PC_MD_FULL) return m * (m - 1);
	if (md == PC_MD_REV) return (m - 1) * (m - 2);
	if (md == PC_MD_TN93) return 3;
	return -1; // undefined
}

double pc_model_lrt(pc_model_t md0, pc_model_t md1, int32_t m, double lr)
{
	double lambda = 2.0 * (lr < 0.0? -lr : lr);
	int32_t df_diff = pc_model_df(md1, m) - pc_model_df(md0, m);
	assert(df_diff != 0);
	if (df_diff < 0) df_diff = -df_diff;
	return kf_chi2_p(df_diff, lambda);
}

double pc_model_BIC(pc_model_t md0, pc_model_t md1, int32_t m, int32_t len, double lr)
{
	double lambda = 2.0 * (lr < 0.0? -lr : lr);
	int32_t df_diff = pc_model_df(md1, m) - pc_model_df(md0, m);
	assert(df_diff != 0);
	if (df_diff < 0) df_diff = -df_diff;
	return lambda - df_diff * log(len);
}

// make matrix symmetric
static inline void pc_model_symm(const double *cnt, int32_t m, double *out)
{
	int32_t a, b;
	for (a = 0; a < m; ++a) { // tmp is symmetric
		out[a * m + a] = cnt[a * m + a];
		for (b = 0; b < a; ++b)
			out[a * m + b] = out[b * m + a] = .5 * (cnt[a * m + b] + cnt[b * m + a]);
	}
}

// calculate residue frequency
static inline double pc_model_freq(const double *cnt, int32_t m, double *pi)
{
	int32_t a, b;
	double tot = 0.0, t;
	for (a = 0; a < m; ++a) {
		for (b = 0, pi[a] = 0.0; b < m; ++b)
			pi[a] += cnt[a * m + b];
		tot += pi[a];
	}
	for (a = 0, t = 1.0 / tot; a < m; ++a)
		pi[a] *= t;
	return tot;
}

static void pc_model_matrix_TN93(const double *cnt, double *out)
{
	double pi[4], tv;
	pc_model_symm(cnt, 4, out);
	pc_model_freq(out, 4, pi);
	tv = .25 * (out[1] / (pi[0] * pi[1]) + out[3] / (pi[0] * pi[3]) + out[6] / (pi[1] * pi[2]) + out[11] / (pi[2] * pi[3]));
	out[1]  = out[4]  = pi[0] * pi[1] * tv;
	out[3]  = out[12] = pi[0] * pi[3] * tv;
	out[6]  = out[9]  = pi[1] * pi[2] * tv;
	out[11] = out[14] = pi[2] * pi[3] * tv;
}

void pc_model_matrix(const double *cnt, int32_t m, pc_model_t md, double *tmp)
{
	if (md == PC_MD_REV) {
		pc_model_symm(cnt, m, tmp);
	} else if (md == PC_MD_TN93) {
		assert(m == 4);
		pc_model_matrix_TN93(cnt, tmp);
	} else { // PC_MD_FULL or undefined
		memcpy(tmp, cnt, sizeof(double) * m * m);
	}
}

double pc_model_dist_TN93(const double *cnt, double *kR, double *kY)
{
	int32_t i;
	double p[16], pi[4], sym[16], tot, t, e_1, e_R, e_Y, pi_R, pi_Y, beta_inv;
	*kR = *kY = -1.0;
	pc_model_symm(cnt, 4, sym);
	tot = pc_model_freq(sym, 4, pi);
	for (i = 0, t = 1.0 / tot; i < 16; ++i)
		p[i] = sym[i] * t;
	t = .25 * (p[1] / (pi[0] * pi[1]) + p[3] / (pi[0] * pi[3]) + p[6] / (pi[1] * pi[2]) + p[11] / (pi[2] * pi[3]));
	if (t > 1.0) return -1.0; // failed
	e_1 = 1.0 - t;
	pi_R = pi[0] + pi[2];
	pi_Y = pi[1] + pi[3];
	e_R = pi_R + pi_Y * e_1 - pi_R * p[2] / (pi[0] * pi[2]);
	e_Y = pi_Y + pi_R * e_1 - pi_Y * p[7] / (pi[1] * pi[3]);
	t = log(e_1);
	*kR = (log(e_R) / t - pi_Y) / pi_R;
	*kY = (log(e_Y) / t - pi_R) / pi_Y;
	beta_inv = 2.0 * (pi[0] * pi[2] * (*kR) + pi[1] * pi[3] * (*kY) + pi_R * pi_Y);
	return -t * beta_inv;
}

void pc_model_dist(pc_tree_t *t, const pc_msa_t *msa, pc_model_t md)
{
	int32_t j;
	double kR, kY;
	assert(md == PC_MD_TN93);
	pc_scfg_alloc(t, msa->len);
	pc_scfg_post_cnt2(t, msa);
	for (j = 0; j < t->n_node - 1; ++j)
		t->node[j]->d = pc_model_dist_TN93(t->node[j]->q->jc, &kR, &kY);
	t->node[j]->d = 0.0; // root is always 0
}
