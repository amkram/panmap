// bwa-mem backend for the single-sample alignment stage. Provides
// bwa_align_reads_direct(), a drop-in twin of align_reads_direct() (mm_align.c)
// that produces the same aligner-agnostic align_pair_result_t results, so the
// downstream BAM builder (conversion.cpp) is unchanged. Selected by --aligner bwa.

#include <limits.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "3rdparty/bwa/bntseq.h"
#include "3rdparty/bwa/bwa.h"
#include "3rdparty/bwa/bwamem.h"

// Aligner-agnostic result structs, mirrored from mm_align.h. Redeclared (not
// included) so this TU does not pull in minimap2 headers alongside bwa's.
// Layout must stay in sync with mm_align.h.
typedef struct {
    int32_t pos;  // 1-based ref position (INT_MAX if unmapped); sort key only
    int32_t rs, re;
    int32_t qs, qe;
    uint8_t mapq;
    uint8_t rev;
    uint8_t proper_frag;
    int32_t n_cigar;
    uint32_t* cigar;  // BAM-encoded; malloc'd, caller frees
} read_align_t;

typedef struct {
    read_align_t r1;
    read_align_t r2;
    int mapped;
} align_pair_result_t;

// bwa's internal CIGAR ops are MIDSH => 0,1,2,3,4; htslib/BAM is MIDNSHP=X =>
// 0..8 (S=4, H=5). Translate op codes; lengths are unchanged.
static const int BWA_OP_TO_BAM[5] = {0, 1, 2, 4, 5};

// Align one forward read against the loaded index; fill `out`. Returns 1 if mapped.
static int bwa_align_one(const bwaidx_t* idx, const mem_opt_t* opt, const char* read, int len, read_align_t* out) {
    memset(out, 0, sizeof(*out));
    out->pos = INT_MAX;
    if (len <= 0) return 0;

    uint8_t* seq = (uint8_t*)malloc(len);
    for (int i = 0; i < len; i++) {
        int c = (unsigned char)read[i];
        seq[i] = c < 4 ? c : nst_nt4_table[c];  // ASCII -> 2-bit
    }

    mem_alnreg_v ar = mem_align1(opt, idx->bwt, idx->bns, idx->pac, len, (char*)seq);

    // Highest-scoring primary region (secondary < 0).
    int best = -1, best_score = 0;
    for (size_t i = 0; i < ar.n; i++) {
        if (ar.a[i].secondary >= 0) continue;
        if (ar.a[i].score > best_score) {
            best_score = ar.a[i].score;
            best = (int)i;
        }
    }

    int mapped = 0;
    if (best >= 0) {
        mem_aln_t a = mem_reg2aln(opt, idx->bns, idx->pac, len, (char*)seq, &ar.a[best]);
        // mem_reg2aln already appends soft-clips and gives a.pos as the 0-based
        // ref position of the first matched base (BAM-ready). Keep its full cigar
        // and set qs/qe so the BAM builder adds no further clips; just remap ops.
        if (a.rid >= 0 && a.n_cigar > 0) {
            uint32_t* cig = (uint32_t*)malloc(a.n_cigar * sizeof(uint32_t));
            int64_t ref_span = 0;
            for (int j = 0; j < a.n_cigar; j++) {
                uint32_t op = a.cigar[j] & 0xf, oplen = a.cigar[j] >> 4;
                uint32_t bam_op = op < 5 ? (uint32_t)BWA_OP_TO_BAM[op] : op;
                cig[j] = (oplen << 4) | bam_op;
                if (bam_op == 0 || bam_op == 2) ref_span += oplen;  // M/D consume ref
            }
            out->pos = (int32_t)(a.pos + 1);
            out->rs = (int32_t)a.pos;
            out->re = (int32_t)(a.pos + ref_span);
            out->qs = 0;
            out->qe = len;
            out->mapq = a.mapq;
            out->rev = a.is_rev;
            out->proper_frag = 0;  // mates aligned independently
            out->n_cigar = a.n_cigar;
            out->cigar = cig;
            mapped = 1;
        }
        free(a.cigar);
    }

    free(ar.a);
    free(seq);
    return mapped;
}

typedef struct {
    const bwaidx_t* idx;
    const mem_opt_t* opt;
    const char** reads;
    const int* r_lens;
    align_pair_result_t* results;
    int paired;
    int start_pair;
    int end_pair;
} bwa_worker_t;

static void* bwa_worker_func(void* arg) {
    bwa_worker_t* w = (bwa_worker_t*)arg;
    if (w->paired) {
        for (int k = w->start_pair; k < w->end_pair; k++) {
            int i = k * 2;
            align_pair_result_t* res = &w->results[k];
            memset(res, 0, sizeof(*res));
            int m1 = bwa_align_one(w->idx, w->opt, w->reads[i], w->r_lens[i], &res->r1);
            int m2 = bwa_align_one(w->idx, w->opt, w->reads[i + 1], w->r_lens[i + 1], &res->r2);
            if (m1 && m2) {
                res->mapped = 1;
                // Both mates mapped: flag as a proper pair. bcftools mpileup
                // discards anomalous pairs by default, so without this the
                // aligned reads never reach genotyping. (Mates are aligned
                // independently, so we don't infer insert size/orientation;
                // R2 is already reverse-complemented upstream.)
                res->r1.proper_frag = 1;
                res->r2.proper_frag = 1;
            } else {
                // Mirror the minimap2 path: a pair counts only if both mates map.
                free(res->r1.cigar);
                free(res->r2.cigar);
                memset(res, 0, sizeof(*res));
                res->r1.pos = INT_MAX;
                res->r2.pos = INT_MAX;
            }
        }
    } else {
        for (int k = w->start_pair; k < w->end_pair; k++) {
            align_pair_result_t* res = &w->results[k];
            memset(res, 0, sizeof(*res));
            res->mapped = bwa_align_one(w->idx, w->opt, w->reads[k], w->r_lens[k], &res->r1);
        }
    }
    return NULL;
}

// Same signature as align_reads_direct (mm_align.c). Builds a temporary bwa
// index for the single reference, aligns reads (mates independently), and
// fills `results`. Each mapped read's cigar is malloc'd; caller frees.
void bwa_align_reads_direct(const char* reference,
                            int n_reads,
                            const char** reads,
                            const char** quality,
                            const char** read_names,
                            const int* r_lens,
                            align_pair_result_t* results,
                            bool pairedEndReads,
                            int n_threads) {
    (void)quality;
    (void)read_names;
    if (n_reads <= 0 || !reference) return;

    bwa_verbose = 0;  // silence bwa's stderr chatter

    // bwa indexes from disk; write the reference to a unique temp FASTA and
    // build a throwaway index next to it, then unlink everything once loaded.
    char tmpl[] = "/tmp/panmap_bwa_XXXXXX";
    int fd = mkstemp(tmpl);
    if (fd < 0) return;
    FILE* fa = fdopen(fd, "w");
    if (!fa) {
        close(fd);
        unlink(tmpl);
        return;
    }
    fprintf(fa, ">ref\n%s\n", reference);
    fclose(fa);

    int build_rc = bwa_idx_build(tmpl, tmpl, 0 /* BWTALGO_AUTO */, 10000000);
    bwaidx_t* idx = (build_rc == 0) ? bwa_idx_load(tmpl, BWA_IDX_ALL) : NULL;

    const char* exts[] = {"", ".amb", ".ann", ".bwt", ".pac", ".sa"};
    char path[4096];
    for (int e = 0; e < 6; e++) {
        snprintf(path, sizeof(path), "%s%s", tmpl, exts[e]);
        unlink(path);
    }
    if (!idx) return;

    mem_opt_t* opt = mem_opt_init();  // defaults tuned for short reads

    if (n_threads < 1) n_threads = 1;
    int n_items = pairedEndReads ? n_reads / 2 : n_reads;
    if (n_threads > n_items) n_threads = n_items;
    if (n_threads < 1) n_threads = 1;

    bwa_worker_t* workers = (bwa_worker_t*)malloc(n_threads * sizeof(bwa_worker_t));
    int chunk = n_items / n_threads, remainder = n_items % n_threads, offset = 0;
    for (int t = 0; t < n_threads; t++) {
        workers[t].idx = idx;
        workers[t].opt = opt;
        workers[t].reads = reads;
        workers[t].r_lens = r_lens;
        workers[t].results = results;
        workers[t].paired = pairedEndReads ? 1 : 0;
        workers[t].start_pair = offset;
        workers[t].end_pair = offset + chunk + (t < remainder ? 1 : 0);
        offset = workers[t].end_pair;
    }

    if (n_threads == 1) {
        bwa_worker_func(&workers[0]);
    } else {
        pthread_t* threads = (pthread_t*)malloc(n_threads * sizeof(pthread_t));
        char* created = (char*)calloc(n_threads, 1);
        for (int t = 0; t < n_threads; t++) {
            // On spawn failure run the chunk synchronously so no reads are dropped,
            // and don't join an uninitialized thread id.
            if (pthread_create(&threads[t], NULL, bwa_worker_func, &workers[t]) == 0)
                created[t] = 1;
            else
                bwa_worker_func(&workers[t]);
        }
        for (int t = 0; t < n_threads; t++)
            if (created[t]) pthread_join(threads[t], NULL);
        free(created);
        free(threads);
    }

    free(workers);
    free(opt);
    bwa_idx_destroy(idx);
}
