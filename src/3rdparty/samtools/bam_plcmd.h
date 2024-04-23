

#include <config.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
#include <limits.h>
#include <errno.h>
#include <sys/stat.h>
#include <getopt.h>
#include <inttypes.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/kstring.h>
#include <htslib/klist.h>
#include <htslib/khash_str2int.h>
#include <htslib/cram.h>
#include "samtools.h"
#include "bedidx.h"
#include "sam_opts.h"
#include "bam_plbuf.h"
#include "sample.h"

#define MPLP_NO_COMP    (1<<2)
#define MPLP_NO_ORPHAN  (1<<3)
#define MPLP_REALN      (1<<4)
#define MPLP_NO_INDEL   (1<<5)
#define MPLP_REDO_BAQ   (1<<6)
#define MPLP_ILLUMINA13 (1<<7)
#define MPLP_IGNORE_RG  (1<<8)
#define MPLP_SMART_OVERLAPS (1<<10)

#define MPLP_PRINT_MAPQ_CHAR (1<<11)
#define MPLP_PRINT_QPOS  (1<<12)
// Start of struct active_cols elements
#define MPLP_PRINT_QNAME (1<<13)
#define MPLP_PRINT_FLAG  (1<<14)
#define MPLP_PRINT_RNAME (1<<15)
#define MPLP_PRINT_POS   (1<<16)
#define MPLP_PRINT_MAPQ  (1<<17)
#define MPLP_PRINT_CIGAR (1<<18)
#define MPLP_PRINT_RNEXT (1<<19)
#define MPLP_PRINT_PNEXT (1<<20)
#define MPLP_PRINT_TLEN  (1<<21)
#define MPLP_PRINT_SEQ   (1<<22)
#define MPLP_PRINT_QUAL  (1<<23)
#define MPLP_PRINT_RLEN  (1<<24)
// Must occur after struct active_cols element list
#define MPLP_PRINT_MODS  (1<<25)
#define MPLP_PRINT_QPOS5 (1<<26)

#define MPLP_PRINT_LAST  (1<<27) // terminator for loop

#define MPLP_MAX_DEPTH 8000
#define MPLP_MAX_INDEL_DEPTH 250

typedef struct {
    int min_mq, flag, min_baseQ, capQ_thres, max_depth, max_indel_depth, all, rev_del;
    int rflag_require, rflag_filter;
    char *reg, *pl_list, *fai_fname, *output_fname;
    faidx_t *fai;
    void *bed, *rghash, *auxlist;
    int argc;
    char **argv;
    char sep, empty, no_ins, no_ins_mods, no_del, no_ends;
    sam_global_args ga;
} mplp_conf_t;

typedef struct {
    char *ref[2];
    int ref_id[2];
    hts_pos_t ref_len[2];
} mplp_ref_t;

#define MPLP_REF_INIT {{NULL,NULL},{-1,-1},{0,0}}

typedef struct {
    samFile *fp;
    hts_itr_t *iter;
    sam_hdr_t *h;
    mplp_ref_t *ref;
    const mplp_conf_t *conf;
} mplp_aux_t;

typedef struct {
    int n;
    int *n_plp, *m_plp;
    bam_pileup1_t **plp;
} mplp_pileup_t;


int mplp_func(void *data, bam1_t *b);
int pileup_cd_create(void *data, const bam1_t *b, bam_pileup_cd *cd);
int pileup_cd_destroy(void *data, const bam1_t *b, bam_pileup_cd *cd);
int mplp_get_ref(mplp_aux_t *ma, int tid, char **ref, hts_pos_t *ref_len);
