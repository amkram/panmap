/*  bam2bcf.h -- variant calling.

    Copyright (C) 2010-2012 Broad Institute.
    Copyright (C) 2012-2022 Genome Research Ltd.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#ifndef BAM2BCF_H
#define BAM2BCF_H

#include <stdint.h>
#include "../samtools/htslib-1.20/htslib/hts.h"
#include "../samtools/htslib-1.20/htslib/vcf.h"

/**
 *  A simplified version of Mann-Whitney U-test is calculated
 *  by default (no CDF) because it is faster and seems to work
 *  better in machine learning filtering. When enabled by setting
 *  CDF_MWU_TESTS, additional annotations will appear on mpileup's
 *  output (RPB2 in addition to RPB, etc.).
 */
#ifndef CDF_MWU_TESTS
#define CDF_MWU_TESTS 0
#endif

#define B2B_INDEL_NULL 10000

#define B2B_FMT_DP      (1<<0)
#define B2B_FMT_SP      (1<<1)
#define B2B_FMT_DV      (1<<2)
#define B2B_FMT_DP4     (1<<3)
#define B2B_FMT_DPR     (1<<4)
#define B2B_INFO_DPR    (1<<5)
#define B2B_FMT_AD      (1<<6)
#define B2B_FMT_ADF     (1<<7)
#define B2B_FMT_ADR     (1<<8)
#define B2B_INFO_AD     (1<<9)
#define B2B_INFO_ADF    (1<<10)
#define B2B_INFO_ADR    (1<<11)
#define B2B_INFO_SCR    (1<<12)
#define B2B_FMT_SCR     (1<<13)
#define B2B_INFO_VDB    (1<<14)
#define B2B_FMT_QS      (1<<15)
#define B2B_FMT_NMBZ    (1<<16) // per-sample NMBZ
#define B2B_INFO_NMBZ   (1<<17)
#define B2B_INFO_BQBZ   (1<<18)
#define B2B_INFO_MQBZ   (1<<19)
#define B2B_INFO_MQSBZ  (1<<20)
#define B2B_INFO_RPBZ   (1<<21)
#define B2B_INFO_SCBZ   (1<<22)
#define B2B_INFO_SGB    (1<<23)
#define B2B_INFO_MIN_PL_SUM (1<<24)
#define B2B_INFO_NM     (1<<25)
#define B2B_INFO_MQ0F   (1<<26)
#define B2B_INFO_IDV    (1<<27)
#define B2B_INFO_IMF    (1<<28)
#define B2B_INFO_FS     (1<<29)

#define B2B_MAX_ALLELES 5
#define B2B_N_NM 32             // number of NMBZ bins, i.e. max number of mismatches


#define B2B_DROP      0
#define B2B_INC_AD    1
#define B2B_INC_AD0   2


// Pileup "client data" for each read to cache per-read information
#define PLP_CD(x) ((plp_cd_t*)((x)->p))
#define PLP_HAS_SOFT_CLIP(cd) (PLP_CD(cd)->i & 1)
#define PLP_HAS_INDEL(cd)     (PLP_CD(cd)->i & 2)
#define PLP_IS_REALN(cd)      (PLP_CD(cd)->i & 4)
#define PLP_SAMPLE_ID(cd)     (PLP_CD(cd)->i >> 3)
#define PLP_QLEN(cd)          (PLP_CD(cd)->qlen)
#define PLP_NM(cd)            (PLP_CD(cd)->nm)
#define PLP_NM_UNSET          -2

#define PLP_SET_SOFT_CLIP(cd)     (PLP_CD(cd)->i |= 1)
#define PLP_SET_INDEL(cd)         (PLP_CD(cd)->i |= 2)
#define PLP_SET_REALN(cd)         (PLP_CD(cd)->i |= 4)
#define PLP_SET_SAMPLE_ID(cd,n)   (PLP_CD(cd)->i |= (n)<<3)

typedef struct
{
    int64_t i;      // used to store sample id and flags for presence of soft-clip and indel
    uint32_t qlen;  // cached output of bam_cigar2qlen(), 0 if unset
    int nm;         // -2 PLP_NM_UNSET; -1 not available; >=0 NM value computed by get_aux_nm()
}
plp_cd_t;


typedef struct __bcf_callaux_t {
    int fmt_flag, ambig_reads;
    int capQ, min_baseQ, max_baseQ, delta_baseQ;
    int openQ, extQ, tandemQ; // for indels
    uint32_t min_support, max_support; // for collecting indel candidates
    double min_frac; // for collecting indel candidates
    float max_frac; // for collecting indel candidates
    int per_sample_flt; // indel filtering strategy
    int *ref_pos, *alt_pos, npos, *ref_mq, *alt_mq, *ref_bq, *alt_bq, *fwd_mqs, *rev_mqs, nqual; // for bias tests
    int *iref_pos, *ialt_pos, *iref_mq, *ialt_mq; // for indels
    int ref_scl[100], alt_scl[100];   // soft-clip length bias; SNP
    int iref_scl[100], ialt_scl[100]; // soft-clip length bias; INDEL
    // for internal uses
    int max_bases;
    int indel_types[4];     // indel lengths
    int indel_win_size, indels_v20, edlib;
    int seqQ_offset; // edlib mode, seqQ=MIN(seqQ, offset - MIN(20,depth)*5);
    int maxins, indelreg, poly_mqual;
    int read_len;
    char *inscns;
    uint16_t *bases;        // 5bit: unused, 6:quality, 1:is_rev, 4:2-bit base or indel allele (index to bcf_callaux_t.indel_types)
    errmod_t *e;
    void *rghash;
    float indel_bias;  // adjusts indel score threshold; lower => call more.
    float del_bias;    // (-.9 < x < .9) error profile; >0 => more del, <0 => more ins
    float vs_ref;      // 0 to 1.  0: score vs next-best. 1: score vs ref
    int32_t *ref_nm, *alt_nm;   // pointers to bcf_call_t.{ref_nm,alt_nm}
    unsigned int nnm[2];        // number of nm observations
    float nm[2];                // cumulative count of mismatches in ref and alt reads
    void *iaux;                 // auxiliary structure for --indels-2.0 calling
    char *chr;                  // current chromosome
} bcf_callaux_t;

// per-sample values
typedef struct {
    uint32_t ori_depth;     // ori_depth = anno[0..3] but before --min-BQ is applied
    unsigned int mq0;
    int32_t *ADF, *ADR, SCR, *QS;   // FMT/QS
    int32_t *ref_nm, *alt_nm;
    // The fields are:
    //      depth fwd   .. ref (0) and non-ref (2)
    //      depth rev   .. ref (1) and non-ref (3)
    //      baseQ       .. ref (4) and non-ref (6)
    //      baseQ^2     .. ref (5) and non-ref (7)
    //      mapQ        .. ref (8) and non-ref (10)
    //      mapQ^2      .. ref (9) and non-ref (11)
    //      minDist     .. ref (12) and non-ref (14)
    //      minDist^2   .. ref (13) and non-ref (15)
    // Note that this probably needs a more thorough fix: int types in
    // bcf_call_t do overflow with high-coverage data, such as exomes, and
    // BCFv2 supports only floats which may not suffice.
    double anno[16];
    float p[25];        // phred-scaled likelihood of each genotype
} bcf_callret1_t;

// values for all samples
typedef struct {
    int tid, pos;
    bcf_hdr_t *bcf_hdr;
    int a[5]; // alleles: ref, alt, alt2, alt3
    float qsum[B2B_MAX_ALLELES];  // INFO/QS tag
    int n, n_alleles, ori_ref, unseen;
    int32_t shift;  // shift is the sum of min_PL before normalization to 0 across all samples
    int n_supp; // number of supporting non-reference reads
    double anno[16];
    unsigned int depth, ori_depth, mq0;
    int32_t *PL, *DP4, *ADR, *ADF, *SCR, *QS, *ref_nm, *alt_nm;
    uint8_t *fmt_arr;
    float vdb; // variant distance bias
    float mwu_pos, mwu_mq, mwu_bq, mwu_mqs, mwu_sc, *mwu_nm, nm[2];
    float seg_bias;
    float strand_bias; // phred-scaled fisher-exact test
    kstring_t tmp;
} bcf_call_t;


#ifdef __cplusplus
extern "C" {
#endif

    bcf_callaux_t *bcf_call_init(double theta, int min_baseQ, int max_baseQ,
                                 int delta_baseQ);
    void bcf_call_destroy(bcf_callaux_t *bca);
    int bcf_call_glfgen(int _n, const bam_pileup1_t *pl, int ref_base, bcf_callaux_t *bca, bcf_callret1_t *r);
    int bcf_call_combine(int n, const bcf_callret1_t *calls, bcf_callaux_t *bca, int ref_base /*4-bit*/, bcf_call_t *call);
    int bcf_call2bcf(bcf_call_t *bc, bcf1_t *b, bcf_callret1_t *bcr, int fmt_flag,
                     const bcf_callaux_t *bca, const char *ref);
    int bcf_call_gap_prep(int n, int *n_plp, bam_pileup1_t **plp, int pos, bcf_callaux_t *bca, const char *ref);
    int bcf_iaux_gap_prep(int n, int *n_plp, bam_pileup1_t **plp, int pos, bcf_callaux_t *bca, const char *ref);
    int bcf_edlib_gap_prep(int n, int *n_plp, bam_pileup1_t **plp, int pos,
                           bcf_callaux_t *bca, const char *ref, int ref_len);

    void bcf_callaux_clean(bcf_callaux_t *bca, bcf_call_t *call);

    int bcf_cgp_l_run(const char *ref, int pos);
    int est_indelreg(int pos, const char *ref, int l, char *ins4);

/* ----------------------------------------------------------------------
 * Shared between bam2bcf_indel.c and bam2bcf_edlib.c
 */

// Take a reference position tpos and convert to a query position (returned).
// This uses the CIGAR string plus alignment c->pos to do the mapping.
//
// *_tpos is returned as tpos if query overlaps tpos, but for deletions
// it'll be either the start (is_left) or end (!is_left) ref position.
int tpos2qpos(const bam1_core_t *c, const uint32_t *cigar, int32_t tpos, int is_left, int32_t *_tpos);

// Identify spft-clip length, position in seq, and clipped seq len
void get_pos(const bcf_callaux_t *bca, bam_pileup1_t *p,
             int *sc_len_r, int *slen_r, int *epos_r, int *end);

// Compute the consensus for this sample 's', minus indels which
// get added later.
char *bcf_cgp_calc_cons(int n, int *n_plp, bam_pileup1_t **plp,
                        int pos, int *types, int n_types,
                        int max_ins, int s);

#ifdef __cplusplus
}
#endif

#endif
