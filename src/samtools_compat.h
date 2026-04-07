/*
 * samtools_compat.h -- Minimal struct definitions extracted from samtools
 * for panmap's pileup module. Avoids linking against libst.a.
 *
 * Original sources: samtools/bam_plcmd.h, samtools/sam_opts.h
 * Copyright (C) 2008-2024 Genome Research Ltd.
 * License: MIT/Expat (same as samtools)
 */

#ifndef SAMTOOLS_COMPAT_H
#define SAMTOOLS_COMPAT_H

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/kstring.h>
#include <htslib/klist.h>
#include <htslib/khash_str2int.h>
#include <htslib/cram.h>

/* --- From sam_opts.h --- */

typedef struct sam_global_args {
    htsFormat in;
    htsFormat out;
    char *reference;
    int nthreads;
    int write_index;
} sam_global_args;

static inline void sam_global_args_init(sam_global_args *ga) {
    if (ga) memset(ga, 0, sizeof(*ga));
}

/* --- From bam_plcmd.h --- */

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
#define MPLP_PRINT_MODS  (1<<25)
#define MPLP_PRINT_QPOS5 (1<<26)
#define MPLP_PRINT_LAST  (1<<27)

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

/* --- Inline replacements for samtools utility functions --- */

static inline void print_error(const char *subcommand, const char *fmt, ...) {
    va_list args;
    fprintf(stderr, "[%s] ", subcommand);
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    fprintf(stderr, "\n");
}

/*
 * bed_overlap is called in pileup.c but only when conf->bed != NULL.
 * panmap never sets conf->bed, so these calls are dead code.
 * Provide a stub that always returns 0 to satisfy the linker.
 */
static inline int bed_overlap(const void *bed, const char *seq,
                              hts_pos_t start, hts_pos_t end) {
    (void)bed; (void)seq; (void)start; (void)end;
    return 0;
}

#endif /* SAMTOOLS_COMPAT_H */
