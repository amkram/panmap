/*
 * pileup.h -- Pileup interface + minimal samtools struct definitions
 *
 * Original samtools sources: bam_plcmd.h, sam_opts.h
 * Copyright (C) 2008-2024 Genome Research Ltd.
 * License: MIT/Expat (same as samtools)
 */

#ifndef PILEUP_H
#define PILEUP_H

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/kstring.h>
#include <htslib/khash_str2int.h>

/* --- From sam_opts.h --- */

typedef struct sam_global_args {
    htsFormat in;
    htsFormat out;
    char* reference;
    int nthreads;
    int write_index;
} sam_global_args;

static inline void sam_global_args_init(sam_global_args* ga) {
    if (ga) memset(ga, 0, sizeof(*ga));
}

/* --- From bam_plcmd.h --- */

#define MPLP_NO_COMP (1 << 2)
#define MPLP_NO_ORPHAN (1 << 3)
#define MPLP_REALN (1 << 4)
#define MPLP_NO_INDEL (1 << 5)
#define MPLP_REDO_BAQ (1 << 6)
#define MPLP_ILLUMINA13 (1 << 7)
#define MPLP_IGNORE_RG (1 << 8)
#define MPLP_SMART_OVERLAPS (1 << 10)

#define MPLP_MAX_DEPTH 8000
#define MPLP_MAX_INDEL_DEPTH 250

typedef struct {
    int min_mq, flag, min_baseQ, capQ_thres, max_depth, max_indel_depth, all, rev_del;
    int rflag_require, rflag_filter;
    char *reg, *pl_list, *fai_fname, *output_fname;
    faidx_t* fai;
    void *bed, *rghash, *auxlist;
    int argc;
    char** argv;
    char sep, empty, no_ins, no_ins_mods, no_del, no_ends;
    sam_global_args ga;
} mplp_conf_t;

typedef struct {
    char* ref[2];
    int ref_id[2];
    hts_pos_t ref_len[2];
} mplp_ref_t;

#define MPLP_REF_INIT {{NULL, NULL}, {-1, -1}, {0, 0}}

typedef struct {
    samFile* fp;
    hts_itr_t* iter;
    sam_hdr_t* h;
    mplp_ref_t* ref;
    const mplp_conf_t* conf;
} mplp_aux_t;

typedef struct {
    int n;
    int *n_plp, *m_plp;
    bam_pileup1_t** plp;
} mplp_pileup_t;

/* --- Inline replacements for samtools utility functions --- */

static inline void print_error(const char* subcommand, const char* fmt, ...) {
    va_list args;
    fprintf(stderr, "[%s] ", subcommand);
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    fprintf(stderr, "\n");
}

static inline int bed_overlap(const void* bed, const char* seq, hts_pos_t start, hts_pos_t end) {
    (void)bed;
    (void)seq;
    (void)start;
    (void)end;
    return 0;
}

/* --- Pileup function --- */

void bam_and_ref_to_mplp(
    sam_hdr_t* header, bam1_t** bam_lines, int nbams, char* ref_string, int lref, kstring_t* mplp_string);

#endif /* PILEUP_H */
