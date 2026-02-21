#include <stdio.h>
#include <string.h>
#include "kommon.h"
#include "phycfg.h"

int main_view(int argc, char *argv[]);

static int usage(FILE *fp)
{
	fprintf(fp, "Usage: phycfg <command> <arguments>\n");
	fprintf(fp, "Commands:\n");
	fprintf(fp, "  view       view phylogenetic tree\n");
	fprintf(fp, "  version    print the version number\n");
	return fp == stdout? 0 : 1;
}

int main(int argc, char *argv[])
{
	int ret = 0;
	kom_realtime();
	if (argc == 1) return usage(stdout);
	else if (strcmp(argv[1], "view") == 0) ret = main_view(argc-1, argv+1);
	else if (strcmp(argv[1], "version") == 0) {
		printf("%s\n", PC_VERSION);
		return 0;
	} else {
		fprintf(stderr, "ERROR: unknown command '%s'\n", argv[1]);
		return 1;
	}

	if (kom_verbose >= 3 && argc > 2 && ret == 0) {
		int i;
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, PC_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, kom_realtime(), kom_cputime(), kom_peakrss() / 1024.0 / 1024.0 / 1024.0);
	}
	return 0;
}
