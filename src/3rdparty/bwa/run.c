// bwa_run.c
#include "bwa.h"
#include <errno.h>
// Include other necessary headers

int run_bwa(int argc, char *argv[])
{
		// Variables
		extern char *bwa_pg;
		int i, ret;
		double t_real;
		kstring_t pg = {0, 1024, (char*)malloc(1024)};  // Initialize with buffer size
		if (!pg.s) {
			fprintf(stderr, "[run_bwa] Memory allocation failed for kstring_t\n");
			return -1;
		}

		// Add error buffer
		char error_buffer[1024];

		// Process start
		t_real = realtime();
		ksprintf(&pg, "@PG\tID:bwa\tPN:bwa\tVN:%s\tCL:%s", 0, argv[0]);
		for (i = 1; i < argc; ++i) {
			ksprintf(&pg, " %s", argv[i]);
		}
		bwa_pg = pg.s;

		if (argc < 2) {
			fprintf(stderr, "[run_bwa] Insufficient arguments\n");
			free(pg.s);
			return usage();
		}

		// Command handling
		ret = -1;  // Initialize ret to an error state
		errno = 0;  // Reset errno before the function call

		fprintf(stderr, "[run_bwa] Executing command: %s\n", argv[1]);
		fprintf(stderr, "[run_bwa] Command arguments:\n");
		for (int j = 0; j < argc; ++j) {
			fprintf(stderr, "  argv[%d]: %s\n", j, argv[j]);
		}

		if (strcmp(argv[1], "fa2pac") == 0) ret = bwa_fa2pac(argc - 1, argv + 1);
		else if (strcmp(argv[1], "pac2bwt") == 0) ret = bwa_pac2bwt(argc - 1, argv + 1);
		else if (strcmp(argv[1], "pac2bwtgen") == 0) ret = bwt_bwtgen_main(argc - 1, argv + 1);
		else if (strcmp(argv[1], "bwtupdate") == 0) ret = bwa_bwtupdate(argc - 1, argv + 1);
		else if (strcmp(argv[1], "bwt2sa") == 0) ret = bwa_bwt2sa(argc - 1, argv + 1);
		else if (strcmp(argv[1], "index") == 0) ret = bwa_index(argc - 1, argv + 1);
		else if (strcmp(argv[1], "aln") == 0) ret = bwa_aln(argc - 1, argv + 1);
		else if (strcmp(argv[1], "samse") == 0) ret = bwa_sai2sam_se(argc - 1, argv + 1);
		else if (strcmp(argv[1], "sampe") == 0) ret = bwa_sai2sam_pe(argc - 1, argv + 1);
		else if (strcmp(argv[1], "bwtsw2") == 0) ret = bwa_bwtsw2(argc - 1, argv + 1);
		else if (strcmp(argv[1], "bwasw") == 0) ret = bwa_bwtsw2(argc - 1, argv + 1);
		else if (strcmp(argv[1], "fastmap") == 0) ret = main_fastmap(argc - 1, argv + 1);
		else if (strcmp(argv[1], "mem") == 0) ret = main_mem(argc - 1, argv + 1);
		else if (strcmp(argv[1], "shm") == 0) ret = main_shm(argc - 1, argv + 1);
		else if (strcmp(argv[1], "pemerge") == 0) ret = main_pemerge(argc - 1, argv + 1);
		else if (strcmp(argv[1], "maxk") == 0) ret = main_maxk(argc - 1, argv + 1);
		else {
			fprintf(stderr, "[run_bwa] Unrecognized command '%s'\n", argv[1]);
			free(pg.s);
			return 1;
		}

		if (ret != 0) {
			// Check if errno was set
			if (errno != 0) {
				strerror_r(errno, error_buffer, sizeof(error_buffer));
				fprintf(stderr, "[run_bwa] Command failed with error: %s\n", error_buffer);
			} else {
				fprintf(stderr, "[run_bwa] Command failed with unknown error (ret = %d)\n", ret);
			}
		}

		err_fflush(stdout);

		if (ret == 0) {
			fprintf(stderr, "[run_bwa] CMD:");
			for (i = 0; i < argc; ++i) {
				fprintf(stderr, " %s", argv[i]);
			}
			fprintf(stderr, "\n[run_bwa] Real time: %.3f sec; CPU: %.3f sec\n", realtime() - t_real, cputime());
		}

		// Free allocated memory
		free(pg.s);
		return ret;
}