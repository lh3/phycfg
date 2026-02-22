#include <ctype.h>
#include <stdint.h>
#include <stdlib.h>
#include "kommon.h"
#include "phycfg.h"

pc_restype_t pc_msa_infer_rt(const pc_msa_t *msa)
{
	int64_t n_letter = 0, n_nt = 0, n_aa = 0;
	int32_t i, j;

	for (i = 0; i < msa->len; ++i)
		for (j = 0; j < msa->n_seq; ++j) {
			uint8_t c = msa->msa[i][j];
			if (!isalpha(c)) continue;
			++n_letter;
			if (kom_nt4_table[c] < 4)  ++n_nt;
			if (kom_aa20_table[c] < 22) ++n_aa;
		}

	if (n_letter == 0)             return PC_RT_UNKNOWN;
	if (n_nt * 2 >= n_letter)      return PC_RT_NT;      /* >= 50% A/C/G/T */
	if (n_aa * 5 >= n_letter * 4)  return PC_RT_AA;      /* >= 80% amino acid */
	return PC_RT_UNKNOWN;
}

void pc_msa_encode(pc_msa_t *msa, pc_restype_t rt)
{
	int32_t i, j;
	const uint8_t *tab;
	uint8_t gap;

	msa->rt = rt, msa->m = 256;
	if (msa->rt == PC_RT_NT)      tab = kom_nt4_table,  msa->m = 4, gap = PC_GAP_NT;
	else if (msa->rt == PC_RT_AA) tab = kom_aa20_table, msa->m = 20, gap = PC_GAP_AA;
	else return;

	for (i = 0; i < msa->len; ++i)
		for (j = 0; j < msa->n_seq; ++j) {
			uint8_t c = msa->msa[i][j];
			msa->msa[i][j] = (c == '-' || c == '.') ? gap : tab[c];
		}
}

void pc_msa_destroy(pc_msa_t *msa)
{
	int32_t i;
	if (msa == NULL) return;
	for (i = 0; i < msa->n_seq; ++i) free(msa->name[i]);
	free(msa->name);
	for (i = 0; i < msa->len; ++i) free(msa->msa[i]);
	free(msa->msa);
	free(msa);
}

void pc_msa_filter(pc_msa_t *msa, int32_t min_cnt, int32_t is_cds)
{
	int32_t i, j, w = 0;
	int32_t step = is_cds ? 3 : 1;

	for (i = 0; i + step <= msa->len; i += step) {
		int32_t cnt = 0;
		for (j = 0; j < msa->n_seq; ++j) {
			int ok = msa->msa[i][j] < (uint8_t)msa->m;
			if (is_cds)
				ok = ok && msa->msa[i+1][j] < (uint8_t)msa->m
				        && msa->msa[i+2][j] < (uint8_t)msa->m;
			if (ok) ++cnt;
		}
		if (cnt >= min_cnt) {
			int32_t k;
			for (k = 0; k < step; ++k) msa->msa[w++] = msa->msa[i+k];
		} else {
			int32_t k;
			for (k = 0; k < step; ++k) free(msa->msa[i+k]);
		}
	}
	for (; i < msa->len; ++i) free(msa->msa[i]); /* trailing incomplete codon */
	msa->len = w;
}
