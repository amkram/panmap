/* The MIT License

   Copyright (c) 2018-     Dana-Farber Cancer Institute
                 2009-2018 Broad Institute, Inc.
                 2008-2009 Genome Research Ltd. (GRL)

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/
#ifndef BWA_H_
#define BWA_H_
#include "kstring.h"   // Assuming these are part of the project
#include "utils.h"

#include <stdint.h>
#include "bntseq.h"
#include "bwt.h"

#define BWA_IDX_BWT 0x1
#define BWA_IDX_BNS 0x2
#define BWA_IDX_PAC 0x4
#define BWA_IDX_ALL 0x7

#define BWA_CTL_SIZE 0x10000

#define BWTALGO_AUTO  0
#define BWTALGO_RB2   1
#define BWTALGO_BWTSW 2
#define BWTALGO_IS    3

#define BWA_DBG_QNAME 0x1

typedef struct {
	bwt_t    *bwt; // FM-index
	bntseq_t *bns; // information on the reference sequences
	uint8_t  *pac; // the actual 2-bit encoded reference sequences with 'N' converted to a random base

	int    is_shm;
	int64_t l_mem;
	uint8_t  *mem;
} bwaidx_t;

typedef struct {
	int l_seq, id;
	char *name, *comment, *seq, *qual, *sam;
} bseq1_t;

extern int bwa_verbose, bwa_dbg;
extern char bwa_rg_id[256];

extern int run_bwa(int argc, char *argv[]);

#ifdef __cplusplus
extern "C" {
#endif

	bseq1_t *bseq_read(int chunk_size, int *n_, void *ks1_, void *ks2_);
	void bseq_classify(int n, bseq1_t *seqs, int m[2], bseq1_t *sep[2]);

	void bwa_fill_scmat(int a, int b, int8_t mat[25]);
	uint32_t *bwa_gen_cigar(const int8_t mat[25], int q, int r, int w_, int64_t l_pac, const uint8_t *pac, int l_query, uint8_t *query, int64_t rb, int64_t re, int *score, int *n_cigar, int *NM);
	uint32_t *bwa_gen_cigar2(const int8_t mat[25], int o_del, int e_del, int o_ins, int e_ins, int w_, int64_t l_pac, const uint8_t *pac, int l_query, uint8_t *query, int64_t rb, int64_t re, int *score, int *n_cigar, int *NM);

	int bwa_idx_build(const char *fa, const char *prefix, int algo_type, int block_size);

	char *bwa_idx_infer_prefix(const char *hint);
	bwt_t *bwa_idx_load_bwt(const char *hint);

	bwaidx_t *bwa_idx_load_from_shm(const char *hint);
	bwaidx_t *bwa_idx_load_from_disk(const char *hint, int which);
	bwaidx_t *bwa_idx_load(const char *hint, int which);
	void bwa_idx_destroy(bwaidx_t *idx);
	int bwa_idx2mem(bwaidx_t *idx);
	int bwa_mem2idx(int64_t l_mem, uint8_t *mem, bwaidx_t *idx);

	void bwa_print_sam_hdr(const bntseq_t *bns, const char *hdr_line);
	char *bwa_set_rg(const char *s);
	char *bwa_insert_header(const char *s, char *hdr);

	int bwa_fa2pac(int argc, char *argv[]);
	int bwa_pac2bwt(int argc, char *argv[]);
	int bwa_bwtupdate(int argc, char *argv[]);
	int bwa_bwt2sa(int argc, char *argv[]);
	int bwa_index(int argc, char *argv[]);
	int bwt_bwtgen_main(int argc, char *argv[]);
	int bwa_aln(int argc, char *argv[]);
	int bwa_sai2sam_se(int argc, char *argv[]);
	int bwa_sai2sam_pe(int argc, char *argv[]);
	int bwa_bwtsw2(int argc, char *argv[]);
	int main_fastmap(int argc, char *argv[]);
	int main_mem(int argc, char *argv[]);
	int main_shm(int argc, char *argv[]);
	int main_pemerge(int argc, char *argv[]);
	int main_maxk(int argc, char *argv[]);

	static int usage()
	{
		fprintf(stderr, "\n");
		fprintf(stderr, "Program: bwa (alignment via Burrows-Wheeler transformation)\n");
		fprintf(stderr, "Contact: Heng Li <hli@ds.dfci.harvard.edu>\n\n");
		fprintf(stderr, "Usage:   bwa <command> [options]\n\n");
		fprintf(stderr, "Command: index         index sequences in the FASTA format\n");
		fprintf(stderr, "         mem           BWA-MEM algorithm\n");
		fprintf(stderr, "         fastmap       identify super-maximal exact matches\n");
		fprintf(stderr, "         pemerge       merge overlapping paired ends (EXPERIMENTAL)\n");
		fprintf(stderr, "         aln           gapped/ungapped alignment\n");
		fprintf(stderr, "         samse         generate alignment (single ended)\n");
		fprintf(stderr, "         sampe         generate alignment (paired ended)\n");
		fprintf(stderr, "         bwasw         BWA-SW for long queries (DEPRECATED)\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "         shm           manage indices in shared memory\n");
		fprintf(stderr, "         fa2pac        convert FASTA to PAC format\n");
		fprintf(stderr, "         pac2bwt       generate BWT from PAC\n");
		fprintf(stderr, "         pac2bwtgen    alternative algorithm for generating BWT\n");
		fprintf(stderr, "         bwtupdate     update .bwt to the new format\n");
		fprintf(stderr, "         bwt2sa        generate SA from BWT and Occ\n");
		fprintf(stderr, "\n");
		fprintf(stderr,
			"Note: To use BWA, you need to first index the genome with `bwa index'.\n"
			"      There are three alignment algorithms in BWA: `mem', `bwasw', and\n"
			"      `aln/samse/sampe'. If you are not sure which to use, try `bwa mem'\n"
			"      first. Please `man ./bwa.1' for the manual.\n\n");
		return 1;
	}

	
#ifdef __cplusplus
}
#endif

#endif
