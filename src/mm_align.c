#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <pthread.h>

#include "3rdparty/minimap2/kalloc.h"
#include "3rdparty/minimap2/khash.h"
#include "3rdparty/minimap2/kseq.h"
#include "3rdparty/minimap2/minimap.h"
#include "3rdparty/minimap2/mmpriv.h"

// Types from mm_align.h — redeclared here to avoid include path issues
// when this file is compiled as part of the bwa target (which doesn't have
// minimap2/ on its include path, but mm_align.h includes <minimap2/mmpriv.h>).
typedef struct {
    int32_t pos;
    int32_t rs, re;
    int32_t qs, qe;
    uint8_t mapq;
    uint8_t rev;
    uint8_t proper_frag;
    int32_t n_cigar;
    uint32_t* cigar;
} read_align_t;

typedef struct {
    read_align_t r1;
    read_align_t r2;
    int mapped;
} align_pair_result_t;

// Forward declaration
static void free_regs(mm_reg1_t* reg, int n_reg);

//
// Preset selection (from minimap2/options.c):
//   avg_len < 500   -> "sr"       short reads  (k=21, w=11, b=8, max_gap=100)
//   avg_len < 5000  -> "map-ont"  ONT/medium   (k=15, w=10, b=4, max_gap=5000)
//   else            -> "map-hifi" HiFi/long    (k=19, w=19, b=4, max_gap=10000)
//
// After preset selection, MM_F_CIGAR is always enabled.
// mm_idx_str arg order: (w, k, is_hpc=0, bucket_bits, n, seqs, names)
//
// for_scoring: when true, use sr scoring parameters (A/B/O/E) but WITHOUT
//   the MM_F_SR/MM_F_FRAG_MODE flags that require paired-end context.
//   mm_map() maps individual reads, so those flags cause all reads to fail.
//
// for_alignment: when true and short reads detected, use the full sr preset
//   WITH MM_F_SR since mm_map_frag handles paired-end correctly.
static mm_idx_t* setup_minimap2(mm_idxopt_t* iopt,
                                mm_mapopt_t* mopt,
                                const char* reference,
                                int n_reads,
                                const int* r_lens,
                                int for_scoring,
                                int paired_end) {
    // Compute average read length
    int avg_len = 150;  // default
    if (n_reads > 0) {
        int64_t total_len = 0;
        for (int i = 0; i < n_reads; i++) total_len += r_lens[i];
        avg_len = (int)(total_len / n_reads);
    }

    mm_verbose = 0;  // suppress warnings

    if (avg_len < 500) {
        if (for_scoring) {
            // Replicate minimap2 sr preset scoring/index params exactly.
            // We DON'T set MM_F_SR because the fork's older align.c has a bug in
            // the ungapped sr alignment path (NULL r->p deref in mm_align1).
            // Without MM_F_SR, minimap2 uses gapped alignment — same scoring params,
            // correct results, no crash.  We still set MM_F_FRAG_MODE + MM_F_HEAP_SORT
            // so mm_map_frag pairs reads correctly.
            mm_set_opt(0, iopt, mopt);
            // Index params (same as sr)
            iopt->k = 21;
            iopt->w = 11;
            // Scoring params (same as sr)
            mopt->a = 2;
            mopt->b = 8;
            mopt->q = 12;
            mopt->e = 2;
            mopt->q2 = 24;
            mopt->e2 = 1;
            mopt->zdrop = 100;
            mopt->zdrop_inv = 100;
            mopt->end_bonus = 10;
            mopt->max_frag_len = 800;
            mopt->max_gap = 100;
            mopt->bw = 100;
            mopt->bw_long = 100;
            mopt->pri_ratio = 0.5f;
            mopt->min_cnt = 2;
            mopt->min_chain_score = 25;
            mopt->min_dp_max = 40;
            mopt->best_n = 20;
            mopt->mid_occ = 1000;
            mopt->max_occ = 5000;
            // PE flags (no MM_F_SR to avoid fork's ungapped-alignment bug)
            mopt->flag |= MM_F_FRAG_MODE | MM_F_HEAP_SORT;
            mopt->pe_ori = 0 << 1 | 1;  // FR
            mopt->pe_bonus = 33;
        } else {
            // For assembly: full sr preset
            mm_set_opt("sr", iopt, mopt);
        }
    } else if (avg_len < 5000) {
        mm_set_opt("map-ont", iopt, mopt);
    } else {
        mm_set_opt("map-hifi", iopt, mopt);
    }
    mopt->flag |= MM_F_CIGAR;

    const char* ref_name = "ref";
    const char* ref_names[1] = {ref_name};
    mm_idx_t* mi = mm_idx_str(iopt->w, iopt->k, 0, 14, 1, &reference, ref_names);
    if (!mi) return NULL;
    mm_mapopt_update(mopt, mi);
    return mi;
}

// Helper: emit one SAM line for a single-end read, or return NULL if unmapped.
// Sets *out_rs to the 1-based reference start position (for sorting) or INT_MAX.
static char* emit_sam_se(const mm_idx_t* mi,
                         const mm_mapopt_t* mopt,
                         int read_len,
                         const char* read_seq,
                         const char* read_qual,
                         const char* read_name,
                         mm_reg1_t* reg,
                         int n_reg,
                         int* out_rs) {
    if (n_reg == 0 || reg[0].score <= 0 || reg[0].score > read_len) {
        *out_rs = INT_MAX;
        return NULL;
    }
    kstring_t sam = {0, 0, 0};
    mm_bseq1_t t;
    t.l_seq = read_len;
    t.seq = (char*)read_seq;
    t.name = (char*)read_name;
    t.qual = (char*)read_qual;

    int n_regs_arr[1] = {n_reg};
    const mm_reg1_t* regs_arr[1] = {reg};
    mm_write_sam3(&sam, mi, &t, 0, 0, 1, n_regs_arr, regs_arr, NULL, 0, 0);
    *out_rs = reg[0].rs + 1;
    return sam.s;
}

// Align reads to reference using minimap2's native seeding.
// r_lens is INPUT read lengths but OVERWRITTEN with 1-based ref start positions.
void align_reads(const char* reference,
                 int n_reads,
                 const char** reads,
                 const char** quality,
                 const char** read_names,
                 int* r_lens,
                 char** sam_alignments,
                 bool pairedEndReads) {
    mm_idxopt_t iopt;
    mm_mapopt_t mopt;

    // for_scoring=1: Use safe sr-like config (no MM_F_SR) to avoid fork's
    // ungapped-alignment NULL deref bug in align.c.  Produces identical scoring
    // with gapped alignment and valid CIGAR strings.
    mm_idx_t* mi = setup_minimap2(&iopt, &mopt, reference, n_reads, r_lens, 1, pairedEndReads ? 1 : 0);
    if (!mi) return;
    mm_tbuf_t* tbuf = mm_tbuf_init();

    if (pairedEndReads) {
        for (int k = 0; k < n_reads / 2; k++) {
            int qlens[2] = {r_lens[k * 2], r_lens[k * 2 + 1]};
            const char* seqs[2] = {reads[k * 2], reads[k * 2 + 1]};
            mm_reg1_t* reg[2] = {NULL, NULL};
            int n_reg[2] = {0, 0};

            mm_map_frag(mi, 2, qlens, seqs, n_reg, reg, tbuf, &mopt, NULL);

            if (n_reg[0] == 0 || n_reg[1] == 0 || reg[0] == NULL || reg[1] == NULL || reg[0]->score <= 0 ||
                reg[1]->score <= 0) {
                sam_alignments[k * 2] = NULL;
                sam_alignments[k * 2 + 1] = NULL;
                r_lens[k * 2] = INT_MAX;
                r_lens[k * 2 + 1] = INT_MAX;
            } else {
                // R1
                kstring_t sam1 = {0, 0, 0};
                mm_bseq1_t t1;
                t1.l_seq = qlens[0];
                t1.seq = (char*)reads[k * 2];
                t1.name = (char*)read_names[k * 2];
                t1.qual = (char*)quality[k * 2];
                const mm_reg1_t* creg[2] = {reg[0], reg[1]};
                mm_write_sam3(&sam1, mi, &t1, 0, 0, 2, n_reg, creg, NULL, 0, 0);
                r_lens[k * 2] = reg[0]->rs + 1;
                sam_alignments[k * 2] = sam1.s;
                // R2
                kstring_t sam2 = {0, 0, 0};
                mm_bseq1_t t2;
                t2.l_seq = qlens[1];
                t2.seq = (char*)reads[k * 2 + 1];
                t2.name = (char*)read_names[k * 2 + 1];
                t2.qual = (char*)quality[k * 2 + 1];
                mm_write_sam3(&sam2, mi, &t2, 1, 0, 2, n_reg, creg, NULL, 0, 0);
                r_lens[k * 2 + 1] = reg[1]->rs + 1;
                sam_alignments[k * 2 + 1] = sam2.s;
            }

            free_regs(reg[0], n_reg[0]);
            free_regs(reg[1], n_reg[1]);
        }
    } else {
        for (int k = 0; k < n_reads; k++) {
            int n_reg = 0;
            mm_reg1_t* reg = mm_map(mi, r_lens[k], reads[k], &n_reg, tbuf, &mopt, NULL);

            sam_alignments[k] =
                emit_sam_se(mi, &mopt, r_lens[k], reads[k], quality[k], read_names[k], reg, n_reg, &r_lens[k]);
            free_regs(reg, n_reg);
        }
    }
    mm_tbuf_destroy(tbuf);
    mm_idx_destroy(mi);
}

// Helper: count errors for one alignment result
// Returns edit distance (blen - mlen + ambiguous bases) for the primary alignment,
// or full read length if unmapped (no alignments returned).
static int64_t count_read_errors(mm_reg1_t* reg, int n_reg, int read_len) {
    if (n_reg > 0 && reg != NULL) {
        mm_reg1_t* r = &reg[0];
        if (r->p && r->blen > 0) {
            return (r->blen - r->mlen + r->p->n_ambi);
        }
        return read_len;
    }
    return read_len;
}

// Helper: free alignment results
static void free_regs(mm_reg1_t* reg, int n_reg) {
    if (reg) {
        for (int j = 0; j < n_reg; j++) {
            if (reg[j].p) free(reg[j].p);
        }
        free(reg);
    }
}

// Alignment scoring function for refinement
// Returns: negative total edit distance (higher = fewer errors = better)
//
// When paired_end is true, reads must be interleaved (R1_0, R2_0, R1_1, ...)
// and n_reads must be even. Uses mm_map_frag(n_segs=2) for proper paired-end
// mapping with fragment chaining, pair bonus, and insert size estimation.
//
// When paired_end is false, maps each read independently with mm_map.
int64_t score_reads_vs_reference(
    const char* reference, int n_reads, const char** reads, const int* r_lens, int kmer_size, bool paired_end) {
    (void)kmer_size;

    if (n_reads <= 0) return 0;

    mm_idxopt_t iopt;
    mm_mapopt_t mopt;

    mm_idx_t* mi = setup_minimap2(&iopt, &mopt, reference, n_reads, r_lens, 1, paired_end ? 1 : 0);
    if (!mi) return 0;
    mm_tbuf_t* tbuf = mm_tbuf_init();

    int64_t total_errors = 0;

    if (paired_end && n_reads >= 2) {
        for (int i = 0; i + 1 < n_reads; i += 2) {
            int qlens[2] = {r_lens[i], r_lens[i + 1]};
            const char* seqs[2] = {reads[i], reads[i + 1]};
            mm_reg1_t* regs[2] = {NULL, NULL};
            int n_regs[2] = {0, 0};

            mm_map_frag(mi, 2, qlens, seqs, n_regs, regs, tbuf, &mopt, NULL);

            total_errors += count_read_errors(regs[0], n_regs[0], r_lens[i]);
            total_errors += count_read_errors(regs[1], n_regs[1], r_lens[i + 1]);

            free_regs(regs[0], n_regs[0]);
            free_regs(regs[1], n_regs[1]);
        }
        // If odd number of reads, map the last one alone
        if (n_reads % 2 != 0) {
            int last = n_reads - 1;
            int n_reg = 0;
            mm_reg1_t* reg = mm_map(mi, r_lens[last], reads[last], &n_reg, tbuf, &mopt, NULL);
            total_errors += count_read_errors(reg, n_reg, r_lens[last]);
            free_regs(reg, n_reg);
        }
    } else {
        // Single-end scoring: map each read independently
        for (int i = 0; i < n_reads; i++) {
            int n_reg = 0;
            mm_reg1_t* reg = mm_map(mi, r_lens[i], reads[i], &n_reg, tbuf, &mopt, NULL);
            total_errors += count_read_errors(reg, n_reg, r_lens[i]);
            free_regs(reg, n_reg);
        }
    }

    mm_tbuf_destroy(tbuf);
    mm_idx_destroy(mi);

    // Return negative total errors (higher = fewer errors = better)
    return -total_errors;
}

// Helper: extract primary alignment fields from mm_reg1_t into read_align_t
static void extract_align_result(mm_reg1_t* reg, int n_reg, read_align_t* out) {
    memset(out, 0, sizeof(*out));
    if (n_reg > 0 && reg && reg[0].p) {
        mm_reg1_t* r = &reg[0];
        out->pos = r->rs + 1;
        out->rs = r->rs;
        out->re = r->re;
        out->qs = r->qs;
        out->qe = r->qe;
        out->mapq = r->mapq;
        out->rev = r->rev;
        out->proper_frag = r->proper_frag;
        out->n_cigar = r->p->n_cigar;
        out->cigar = (uint32_t*)malloc(r->p->n_cigar * sizeof(uint32_t));
        memcpy(out->cigar, r->p->cigar, r->p->n_cigar * sizeof(uint32_t));
    } else {
        out->pos = INT_MAX;
    }
}

typedef struct {
    const mm_idx_t* mi;
    const mm_mapopt_t* mopt;
    int start_pair;
    int end_pair;
    const char** reads;
    const char** quality;
    const char** read_names;
    const int* r_lens;
    align_pair_result_t* results;
    int paired;
} align_worker_t;

static void* align_worker_func(void* arg) {
    align_worker_t* w = (align_worker_t*)arg;
    mm_tbuf_t* tbuf = mm_tbuf_init();

    if (w->paired) {
        for (int k = w->start_pair; k < w->end_pair; k++) {
            int i = k * 2;
            int qlens[2] = {w->r_lens[i], w->r_lens[i + 1]};
            const char* seqs[2] = {w->reads[i], w->reads[i + 1]};
            mm_reg1_t* reg[2] = {NULL, NULL};
            int n_reg[2] = {0, 0};

            mm_map_frag(w->mi, 2, qlens, seqs, n_reg, reg, tbuf, w->mopt, NULL);

            align_pair_result_t* res = &w->results[k];
            memset(res, 0, sizeof(*res));

            if (n_reg[0] > 0 && n_reg[1] > 0 && reg[0] && reg[1] && reg[0]->score > 0 && reg[1]->score > 0) {
                res->mapped = 1;
                extract_align_result(reg[0], n_reg[0], &res->r1);
                extract_align_result(reg[1], n_reg[1], &res->r2);
            } else {
                res->mapped = 0;
                res->r1.pos = INT_MAX;
                res->r2.pos = INT_MAX;
            }

            free_regs(reg[0], n_reg[0]);
            free_regs(reg[1], n_reg[1]);
        }
    } else {
        for (int k = w->start_pair; k < w->end_pair; k++) {
            int n_reg = 0;
            mm_reg1_t* reg = mm_map(w->mi, w->r_lens[k], w->reads[k], &n_reg, tbuf, w->mopt, NULL);

            align_pair_result_t* res = &w->results[k];
            memset(res, 0, sizeof(*res));

            if (n_reg > 0 && reg && reg[0].score > 0 && reg[0].score <= w->r_lens[k]) {
                res->mapped = 1;
                extract_align_result(reg, n_reg, &res->r1);
            } else {
                res->r1.pos = INT_MAX;
            }
            free_regs(reg, n_reg);
        }
    }

    mm_tbuf_destroy(tbuf);
    return NULL;
}

void align_reads_direct(const char* reference,
                        int n_reads,
                        const char** reads,
                        const char** quality,
                        const char** read_names,
                        const int* r_lens,
                        align_pair_result_t* results,
                        bool pairedEndReads,
                        int n_threads) {
    mm_idxopt_t iopt;
    mm_mapopt_t mopt;
    mm_idx_t* mi = setup_minimap2(&iopt, &mopt, reference, n_reads, r_lens, 1, pairedEndReads ? 1 : 0);
    if (!mi) return;

    if (n_threads < 1) n_threads = 1;
    int n_items = pairedEndReads ? n_reads / 2 : n_reads;
    if (n_threads > n_items) n_threads = n_items;
    if (n_threads < 1) n_threads = 1;

    align_worker_t* workers = (align_worker_t*)malloc(n_threads * sizeof(align_worker_t));

    int chunk = n_items / n_threads;
    int remainder = n_items % n_threads;
    int offset = 0;

    for (int t = 0; t < n_threads; t++) {
        workers[t].mi = mi;
        workers[t].mopt = &mopt;
        workers[t].reads = reads;
        workers[t].quality = quality;
        workers[t].read_names = read_names;
        workers[t].r_lens = r_lens;
        workers[t].results = results;
        workers[t].paired = pairedEndReads ? 1 : 0;
        workers[t].start_pair = offset;
        workers[t].end_pair = offset + chunk + (t < remainder ? 1 : 0);
        offset = workers[t].end_pair;
    }

    if (n_threads == 1) {
        align_worker_func(&workers[0]);
    } else {
        pthread_t* threads = (pthread_t*)malloc(n_threads * sizeof(pthread_t));
        for (int t = 0; t < n_threads; t++) pthread_create(&threads[t], NULL, align_worker_func, &workers[t]);
        for (int t = 0; t < n_threads; t++) pthread_join(threads[t], NULL);
        free(threads);
    }

    free(workers);
    mm_idx_destroy(mi);
}
