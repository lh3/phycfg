#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "kommon.h"
#include "ketopt.h"
#include "phycfg.h"

pc_scfg_t *pc_scfg_new(int32_t n_node, int32_t m)
{
	pc_scfg_t *s;
	int32_t i, n_dbl = (sizeof(pc_scfg_t) * n_node + sizeof(double) - 1) / sizeof(double);
	double *x;
	x = kom_calloc(double, n_dbl + n_node * 3 * m);
	s = (pc_scfg_t*)x;
	x += n_dbl;
	for (i = 0; i < n_node; ++i) {
		s[i].alpha  = x, x += m;
		s[i].alpha2 = x, x += m;
		s[i].beta   = x, x += m;
	}
	return s;
}

// shape of p: (n,m,m), where n is the number of nodes and m is the size of the alphabet
void pc_transmat_init(double *p, int32_t m, const pc_tree_t *t)
{
	int32_t k;
	assert(m > 1);
	for (k = 0; k < t->n_node; ++k) {
		int32_t i, j;
		double diag, off, *pk = p + (size_t)k * m * m;
		if (k < t->n_node - 1) { // not root
			double d = t->node[k]->d;
			d = d > 0.0? d : 1e-3;
			diag = exp(-d), off = (1.0 - diag) / (m - 1);
		} else diag = off = 1.0 / m; // flat prior at the root
		for (i = 0; i < m; ++i)
			for (j = 0; j < m; ++j)
				pk[i * m + j] = i == j ? diag : off;
	}
}

static inline void pc_scfg_emit(int32_t m, int32_t c, double *alpha)
{
	int32_t i;
	if (c < m) {
		for (i = 0; i < m; ++i) alpha[i] = 0.0;
		alpha[c] = 1.0;
	} else {
		for (i = 0; i < m; ++i) alpha[i] = 1.0;
	}
}

double pc_scfg_inside(const pc_tree_t *t, const double *p, const pc_msa_t *msa, int32_t pos, pc_scfg_t *sd)
{
	int32_t i, k, a, b, m = msa->m;
	double logh = 0.0, sum;
	for (i = 0; i < t->n_node; ++i) {
		const pc_node_t *v = t->node[i]; /* post-order */
		double *alpha = sd[i].alpha, *alpha2 = sd[i].alpha2, h;
		assert(v->n_child == 0 || v->n_child == 2); // only works for binary trees
		if (v->n_child == 0) { /* leaf: initialize alpha from emission */
			assert(v->seq_id >= 0 && v->seq_id < msa->n_seq);
			pc_scfg_emit(m, msa->msa[pos][v->seq_id], alpha);
		} else { /* internal: alpha(v,a) = prod_k alpha'(child_k, a) */
			for (a = 0; a < m; ++a) alpha[a] = 1.0;
			for (k = 0; k < v->n_child; ++k) {
				double *a2k = sd[v->child[k]->ftime].alpha2;
				for (a = 0; a < m; ++a) alpha[a] *= a2k[a];
			}
		}
		h = 0.0;
		for (a = 0; a < m; ++a) h += alpha[a];
		for (a = 0; a < m; ++a) alpha[a] /= h;
		/* alpha'(v,a) = sum_b p(b|a) * alpha~(v,b); stored in alpha2 for parent */
		for (a = 0; a < m; ++a) { // NB: alpha2[] is not defined at the root, though it is calculated anyway
			const double *pi = p + (size_t)i * m * m;
			double s = 0.0;
			for (b = 0; b < m; ++b) s += pi[a*m + b] * alpha[b];
			alpha2[a] = s;
		}
		sd[i].h = h;
		logh += log(h);
	}
	for (a = 0, sum = 0.0; a < m; ++a)
		sum += sd[t->n_node - 1].alpha[a] * p[(size_t)(t->n_node - 1) * m * m + a];
	return logh + log(sum); // this is equal to logh + log(h(root) * \sum_a alpha~(root,a) * beta~(root,a))
}

void pc_scfg_outside(const pc_tree_t *t, const double *p, int32_t m, pc_scfg_t *sd)
{
	int32_t i, a, b, root_idx = t->n_node - 1;
	double *sib = kom_malloc(double, m);

	/* beta~(root,a) = q(a)/h_root; q(a)=p[root_idx][a] */
	for (a = 0; a < m; ++a) sd[root_idx].beta[a] = p[(size_t)root_idx * m * m + a] / sd[root_idx].h;

	/* Pre-order: root (highest index) down to leaves; parent always visited first */
	for (i = root_idx - 1; i >= 0; --i) {
		pc_node_t *u = t->node[i];
		pc_node_t *v = u->parent;
		int32_t k;
		double *beta_u = sd[i].beta, *beta_v = sd[v->ftime].beta;
		const double *pi = p + (size_t)i * m * m;

		/* sib[a] = prod_k alpha'~(sibling_k, a) for all siblings of u */
		for (a = 0; a < m; ++a) sib[a] = 1.0;
		for (k = 0; k < v->n_child; ++k)
			if (v->child[k] != u) {
				double *a2k = sd[v->child[k]->ftime].alpha2;
				for (a = 0; a < m; ++a) sib[a] *= a2k[a];
			}

		/* beta~(u,b) = (1/h_u) * sum_a p(b|a) * beta~(v,a) * sib[a] */
		for (b = 0; b < m; ++b) {
			double s = 0.0;
			for (a = 0; a < m; ++a)
				s += pi[a*m + b] * beta_v[a] * sib[a];
			beta_u[b] = s / sd[i].h;
		}
	}
	free(sib);
}

/* Compute eta~[n_node*m*m] from inside/outside values already stored in sd.
 * For non-root u with parent v and sibling(s) w (except the root):
 *   eta~(u,b|a) = beta~(v,a) * alpha~(u,b) * prod_k alpha'~(sib_k,a) */
void pc_scfg_eta(const pc_tree_t *t, int32_t m, const pc_scfg_t *sd, double *eta)
{
	int32_t i, k, a, b, root_idx = t->n_node - 1;
	double *sib = kom_malloc(double, m);

	for (i = 0; i < root_idx; ++i) { // eta~ is not defined at the root
		const pc_node_t *u = t->node[i], *v = u->parent;
		double *eta_i = eta + (size_t)i * m * m;
		const double *alpha_u = sd[i].alpha, *beta_v = sd[v->ftime].beta;
		for (a = 0; a < m; ++a) sib[a] = 1.0;
		for (k = 0; k < v->n_child; ++k)
			if (v->child[k] != u) {
				const double *a2k = sd[v->child[k]->ftime].alpha2;
				for (a = 0; a < m; ++a) sib[a] *= a2k[a];
			}
		for (a = 0; a < m; ++a)
			for (b = 0; b < m; ++b)
				eta_i[a*m + b] = beta_v[a] * alpha_u[b] * sib[a];
	}
	free(sib);
}

void pc_scfg_eta_nni(const pc_tree_t *t, int32_t m, const pc_scfg_t *sd, double *eta[3])
{
	int32_t u;
	for (u = 0; u < t->n_node - 1; ++u) {
		const pc_node_t *up = t->node[u], *vp = up->parent;
		int32_t v, x, y, w, a, b; // original topology: ((x,y)u,w)v
		double *eta0_u, *eta1_u;
		if (vp == 0 || up->n_child == 0) return;
		assert(up->n_child == 2);
		v = vp->ftime;
		w = vp->child[(vp->child[0] == up)]->ftime;
		x = up->child[0]->ftime;
		y = up->child[1]->ftime;
		eta0_u = &eta[1][u * m * m];
		eta1_u = &eta[2][u * m * m];
		for (a = 0; a < m; ++a) {
			double s, q = sd[v].beta[a] / sd[u].h;
			for (b = 0, s = q * sd[x].alpha2[a]; b < m; ++b) // alt topology1: ((w,y)u,x)v
				eta0_u[a * m + b] = s * sd[w].alpha2[b] * sd[y].alpha2[b];
			for (b = 0, s = q * sd[y].alpha2[a]; b < m; ++b) // alt topology2: ((x,w)u,y)v
				eta1_u[a * m + b] = s * sd[w].alpha2[b] * sd[x].alpha2[b];
		}
	}
}

/* Compute posterior counts into cnt[n_node*m*m] (zeroed on entry) and return
 * the total log likelihood summed over all alignment columns. */
double pc_scfg_post_cnt(const pc_tree_t *t, const double *p, const pc_msa_t *msa, pc_scfg_t *sd, double *cnt)
{
	int32_t i, j, k, a, b, m = msa->m;
	double loglk = 0.0;
	double *tmp = kom_malloc(double, m * m);
	double *sib = kom_malloc(double, m);
	memset(cnt, 0, (size_t)t->n_node * m * m * sizeof(double));
	for (i = 0; i < msa->len; ++i) {
		double sum, *cnt_j;
		loglk += pc_scfg_inside(t, p, msa, i, sd);
		pc_scfg_outside(t, p, m, sd);
		/* for each non-root branch, accumulate cnt[j][a*m+b] =
		 * p(b|a) * alpha~(u,b) * beta~(par,a) * prod_k alpha'~(sib_k,a) */
		for (j = 0; j < t->n_node - 1; ++j) {
			const pc_node_t *u = t->node[j], *v = u->parent;
			double *alpha_u = sd[j].alpha, *beta_v = sd[v->ftime].beta;
			const double *pj = p + (size_t)j * m * m;
			cnt_j = cnt + (size_t)j * m * m;
			for (a = 0; a < m; ++a) sib[a] = 1.0;
			for (k = 0; k < v->n_child; ++k) { // prod_k alpha'~(sib_k,a)
				if (v->child[k] != u) {
					double *a2k = sd[v->child[k]->ftime].alpha2;
					for (a = 0; a < m; ++a) sib[a] *= a2k[a];
				}
			}
			for (a = 0, sum = 0.0; a < m; ++a) {
				double x = beta_v[a] * sib[a];
				for (b = 0; b < m; ++b)
					sum += (tmp[a*m + b] = pj[a*m + b] * alpha_u[b] * x);
			}
			sum = 1.0 / sum;
			for (a = 0; a < m; ++a)
				for (b = 0; b < m; ++b)
					cnt_j[a*m + b] += tmp[a*m + b] * sum;
		}
		// special case: the root
		assert(j == t->n_node - 1);
		cnt_j = cnt + (size_t)j * m * m;
		for (a = 0, sum = 0.0; a < m; ++a)
			sum += sd[j].alpha[a] * sd[j].beta[a];
		sum = 1.0 / sum;
		for (a = 0; a < m; ++a)
			for (b = 0; b < m; ++b)
				cnt_j[a*m + b] += sd[j].alpha[b] * sd[j].beta[b] * sum;
	}
	free(sib); free(tmp);
	return loglk;
}

double pc_scfg_em_basic(const pc_tree_t *t, double *p, const pc_msa_t *msa, pc_scfg_t *sd)
{
	int32_t j, a, b, m = msa->m;
	double *cnt = kom_malloc(double, (size_t)t->n_node * m * m);
	double loglk = pc_scfg_post_cnt(t, p, msa, sd, cnt);

	/* M step: renormalize each row of p from accumulated counts */
	for (j = 0; j < t->n_node; ++j) {
		double *cnt_j = cnt + (size_t)j * m * m, *p_j = p + (size_t)j * m * m;
		for (a = 0; a < m; ++a) {
			double s = 0.0;
			for (b = 0; b < m; ++b) s += cnt_j[a*m + b];
			if (s > 0.0)
				for (b = 0; b < m; ++b) p_j[a*m + b] = cnt_j[a*m + b] / s;
		}
	}

	free(cnt);
	return loglk;
}

double pc_scfg_em_iter(const pc_tree_t *t, const pc_msa_t *msa, int32_t max_iter, double *p, pc_scfg_t *sd)
{
	int32_t k;
	double loglk;
	pc_transmat_init(p, msa->m, t);
	for (k = 0; k < max_iter; ++k)
		loglk = pc_scfg_em_basic(t, p, msa, sd);
	return loglk;
}

int main_scfg(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int32_t i, reroot = 0, max_iter = 20;
	pc_scfg_t *sd;
	double *p;

	while (ketopt(&o, argc, argv, 1, "rm:", 0) >= 0) {
		if (o.opt == 'r') reroot = 1;
		else if (o.opt == 'm') max_iter = atoi(o.arg);
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: phycfg scfg [options] <tree.nhx.gz> <aln.mfa.gz>\n");
		return 1;
	}

	pc_tree_t *t = pc_tree_read(argv[o.ind]);
	if (t == NULL) return 1;
	pc_msa_t *msa = pc_msa_read(argv[o.ind + 1]);
	if (msa == NULL) { pc_tree_destroy(t); return 1; }

	pc_msa_encode(msa, pc_msa_infer_rt(msa));
	assert(msa->rt == PC_RT_NT || msa->rt == PC_RT_CODON); // only for nucleotide for now
	pc_tree_match_msa(t, msa);
	sd = pc_scfg_new(t->n_node, msa->m);
	p = kom_calloc(double, (size_t)t->n_node * msa->m * msa->m);

	if (reroot) { // try all the other roots, perform EM and calculate loglk
		int32_t max_i = -1;
		double max_lk = -1e300;
		for (i = 0; i < t->n_node; ++i) {
			double loglk;
			pc_tree_t *s = pc_tree_clone(t); // don't touch t
			pc_tree_reroot(s, i, -1.0);
			loglk = pc_scfg_em_iter(s, msa, max_iter, p, sd);
			if (max_lk < loglk) max_lk = loglk, max_i = i;
			fprintf(stderr, "RT\t%d\t%.6f\n", i, loglk);
			pc_tree_destroy(s);
		}
		if (max_i >= 0) {
			char *str = NULL;
			int32_t max = 0;
			pc_tree_t *s = pc_tree_clone(t);
			pc_tree_reroot(s, max_i, -1.0);
			pc_tree_format(s, &str, &max);
			puts(str);
			free(str);
			pc_tree_destroy(s);
		}
	} else {
		pc_transmat_init(p, msa->m, t);
		for (i = 0; i < max_iter; ++i) {
			double loglk = pc_scfg_em_basic(t, p, msa, sd);
			printf("LK\t%.6f\n", loglk);
		}
		printf("BF");
		for (i = 0; i < msa->m; ++i) printf("\t%.4f", p[(size_t)(t->n_node-1) * msa->m * msa->m + i]);
		printf("\n");
	}

	free(p); free(sd);
	pc_tree_destroy(t);
	pc_msa_destroy(msa);
	return 0;
}
