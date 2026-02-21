#include <ctype.h>
#include <stdint.h>
#include "kommon.h"
#include "phycfg.h"

pc_restype_t pc_msa_infer_rt(const pc_msa_t *msa)
{
	int64_t n_letter = 0, n_nt = 0, n_aa = 0;
	int32_t i, j;

	for (i = 0; i < msa->n_pos; ++i)
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

	msa->rt = rt;
	if (msa->rt == PC_RT_NT)      tab = kom_nt4_table;
	else if (msa->rt == PC_RT_AA) tab = kom_aa20_table;
	else return;

	for (i = 0; i < msa->n_pos; ++i)
		for (j = 0; j < msa->n_seq; ++j)
			msa->msa[i][j] = tab[msa->msa[i][j]];
}
