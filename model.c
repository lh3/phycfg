#include <string.h>
#include "phycfg.h"

pc_model_t pc_model_from_str(const char *model)
{
	if (strcmp(model, "null") == 0 || strcmp(model, "NULL") == 0) return PC_MD_NULL;
	if (strcmp(model, "rev") == 0 || strcmp(model, "GTR") == 0 || strcmp(model, "gtr") == 0) return PC_MD_REV;
	if (strcmp(model, "HKY") == 0 || strcmp(model, "hky") == 0) return PC_MD_HKY;
	return PC_MD_ERR;
}
