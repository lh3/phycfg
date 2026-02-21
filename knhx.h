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

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif

knhx1_t *kn_parse(const char *nhx, int *_n, int *_max, int *_error, char **en);
void kn_destroy(int n, knhx1_t *a);
void kn_format(const knhx1_t *node, int root, kstring_t *s);
knhx1_t *kn_extract_marked(const knhx1_t *a0, int n0, int *n_out_);

#ifdef __cplusplus
}
#endif

#endif
