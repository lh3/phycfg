#ifndef KNHX_H_
#define KNHX_H_

#include <stddef.h>

#define KN_ERR_NO_LEFT   0x01
#define KN_ERR_NO_RGHT   0x02
#define KN_ERR_BRACKET   0x04
#define KN_ERR_COLON     0x08

typedef struct {
	int parent, n;
	unsigned long aux; // auxiliary information; not used by knhx
	int *child;
	char *name;
	double d;
} knhx1_t;

#ifdef __cplusplus
extern "C" {
#endif

knhx1_t *kn_parse(const char *nhx, int *_n, int *_max, int *_error, char **en);
void kn_destroy(int n, knhx1_t *a);

#ifdef __cplusplus
}
#endif

#endif
