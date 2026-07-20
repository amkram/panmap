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
#include "3rdparty/bwa/bwtaln.h"
#include "3rdparty/bwa/bwase.h"
#include "3rdparty/bwa/bntseq.h"
#include "3rdparty/bwa/bwt.h"

extern int bwa_trim_read(int trim_qual, bwa_seq_t *p);

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
    char* md;
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

    // Known, benign data race: mem_align1 -> mem_mark_primary_se calls lrand48()
    // (a non-thread-safe libc PRNG on shared global state) to break ties among
    // equal-scoring primary hits. Concurrent workers racing it only perturbs which
    // of several equally-good alignments is picked as primary; it does not crash.
    // A fix would need a per-thread RNG inside vendored bwa or serializing the
    // alignment, both worse than the nondeterminism.
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
        // and set qs/qe so the BAM builder adds no further clips; remap ops.
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
            out->md = NULL;
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
                // aligned reads never reach genotyping. Mates are aligned
                // independently, so insert size/orientation is not inferred;
                // R2 is already reverse-complemented upstream.
                res->r1.proper_frag = 1;
                res->r2.proper_frag = 1;
            } else {
                // Mirror the minimap2 path: a pair counts only if both mates map.
                free(res->r1.cigar);
                free(res->r2.cigar);
                free(res->r1.md);
                free(res->r2.md);
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
                            const char* refName,
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

static const int BWA_ALN_OP_TO_BAM[4] = {0, 1, 2, 4};  // aln M/I/D/S -> BAM M/I/D/S

static void aln_seq_reverse(int len, ubyte_t* s, int is_comp) {
  for (int i = 0; i < (len >> 1); i++) {
    ubyte_t t = s[len - 1 - i];
    if (is_comp) {
      s[len - 1 - i] = s[i] < 4 ? 3 - s[i] : 4;
      s[i] = t < 4 ? 3 - t : 4;
    } else {
      s[len - 1 - i] = s[i];
      s[i] = t;
    }
  }
  if ((len & 1) && is_comp) {
    int i = len >> 1;
    s[i] = s[i] < 4 ? 3 - s[i] : 4;
  }
}

static void aln_encode_seq(bwa_seq_t* p, const char* read, const char* qual, int len, int trim_qual) {
  memset(p, 0, sizeof(*p));
  p->tid = -1;
  p->len = p->full_len = p->clip_len = len;
  p->seq = (ubyte_t*)malloc(len);
  for (int i = 0; i < len; i++) p->seq[i] = nst_nt4_table[(unsigned char)read[i]];
  if (qual) {
    p->qual = (ubyte_t*)malloc(len + 1);
    memcpy(p->qual, qual, len);
    p->qual[len] = 0;
    if (trim_qual >= 1) bwa_trim_read(trim_qual, p);  // shrinks p->len; must precede the reversal
  }
  // Mirror bwa_read_seq exactly: seq is reversed, rseq is reverse-complemented,
  // and bwa_refine_gapped reverses seq back. Deviating here silently flips strands.
  p->rseq = (ubyte_t*)malloc(p->full_len);
  memcpy(p->rseq, p->seq, p->len);
  aln_seq_reverse(p->len, p->seq, 0);
  aln_seq_reverse(p->len, p->rseq, 1);
  p->name = strdup("");
}

static int aln_fill_from_seq(read_align_t* out, const bwa_seq_t* s) {
  memset(out, 0, sizeof(*out));
  out->pos = INT_MAX;
  if (s->type == BWA_TYPE_NO_MATCH) return 0;

  int core_n = (s->cigar && s->n_cigar > 0) ? s->n_cigar : 1;
  uint32_t* cig = (uint32_t*)malloc((core_n + 1) * sizeof(uint32_t));  // +1 for a trim soft-clip
  int nc = 0;
  int64_t qlen = 0, ref_span = 0;

  if (s->cigar && s->n_cigar > 0) {
    for (int j = 0; j < s->n_cigar; j++) {
      uint32_t bam_op = (uint32_t)BWA_ALN_OP_TO_BAM[__cigar_op(s->cigar[j])];
      uint32_t oplen = (uint32_t)__cigar_len(s->cigar[j]);
      cig[nc++] = (oplen << 4) | bam_op;
      if (bam_op == 0 || bam_op == 1) qlen += oplen;
      if (bam_op == 0 || bam_op == 2) ref_span += oplen;
    }
  } else {
    // refine_gapped emits no cigar for an ungapped hit; expand to <len>M
    cig[nc++] = ((uint32_t)s->len << 4);
    qlen += s->len;
    ref_span += s->len;
  }

  int trimmed = (int)s->full_len - (int)qlen;
  if (trimmed > 0) {
    uint32_t clip = ((uint32_t)trimmed << 4) | 4u;
    if (s->strand) {  // reverse: trimmed 3' bases land at the left of the ref-forward cigar
      memmove(cig + 1, cig, nc * sizeof(uint32_t));
      cig[0] = clip;
      nc++;
    } else {
      cig[nc++] = clip;
    }
  }

  out->rs = (int32_t)s->pos;
  out->pos = (int32_t)(s->pos + 1);
  out->re = (int32_t)(s->pos + ref_span);
  out->qs = 0;
  out->qe = (int32_t)s->full_len;
  out->mapq = s->mapQ;
  out->rev = s->strand;
  out->proper_frag = 0;
  out->n_cigar = nc;
  out->cigar = cig;
  out->md = s->md ? strdup(s->md) : NULL;
  return 1;
}

typedef struct {
  const bwaidx_t* idx;
  const gap_opt_t* opt;
  const char** reads;
  const char** quals;
  const int* r_lens;
  align_pair_result_t* results;
  int paired;
  int start_pair;
  int end_pair;
} bwa_aln_worker_t;

static void* bwa_aln_worker_func(void* arg) {
  bwa_aln_worker_t* w = (bwa_aln_worker_t*)arg;
  const int n_frag = w->end_pair - w->start_pair;
  const int n_seq = w->paired ? n_frag * 2 : n_frag;
  if (n_seq <= 0) return NULL;

  bwa_seq_t* seqs = (bwa_seq_t*)calloc(n_seq, sizeof(bwa_seq_t));
  for (int j = 0; j < n_seq; j++) {
    int r = w->paired ? (w->start_pair * 2 + j) : (w->start_pair + j);
    aln_encode_seq(&seqs[j], w->reads[r], w->quals ? w->quals[r] : NULL,
                   w->r_lens[r], w->opt->trim_qual);
  }

  // One backtracking-stack allocation for the whole chunk (see notes on batching).
  bwa_cal_sa_reg_gap(0, w->idx->bwt, n_seq, seqs, w->opt);
  for (int j = 0; j < n_seq; j++) {
    // bwa_aln2seq uses drand48() to break ties among equal hits: same benign,
    // non-thread-safe PRNG race as the MEM path's lrand48; perturbs which of
    // several equally-good positions is chosen, nothing worse.
    bwa_aln2seq(seqs[j].n_aln, seqs[j].aln, &seqs[j]);
    bwa_cal_pac_pos_core(w->idx->bns, w->idx->bwt, &seqs[j], w->opt->max_diff, w->opt->fnr);
  }

  // bwa_cal_sa_reg_gap frees seq->seq and seq->rseq once the FM-index search is
  // done. The file-based driver re-reads the sequence before refine; we have no
  // file, so rebuild both buffers here. bwa_refine_gapped reverses seq->seq in
  // place and reads rseq for reverse-strand reads, so both must be restored in the
  // SAME orientation aln_encode_seq produced (seq reversed, rseq reverse-complemented).
  // We must not shrink len again: trimming already ran in the first encode, so
  // rebuild against the CURRENT seqs[j].len / full_len, without re-trimming.
  for (int j = 0; j < n_seq; j++) {
    int r = w->paired ? (w->start_pair * 2 + j) : (w->start_pair + j);
    int L = seqs[j].len;              // post-trim length refine will reverse
    int FL = seqs[j].full_len;        // original length; rseq is allocated at full_len
    seqs[j].seq = (ubyte_t*)malloc(FL);
    for (int i = 0; i < FL; i++)
      seqs[j].seq[i] = nst_nt4_table[(unsigned char)w->reads[r][i]];
    seqs[j].rseq = (ubyte_t*)malloc(FL);
    memcpy(seqs[j].rseq, seqs[j].seq, L);
    aln_seq_reverse(L, seqs[j].seq, 0);
    aln_seq_reverse(L, seqs[j].rseq, 1);
  }

  bwa_refine_gapped(w->idx->bns, n_seq, seqs, w->idx->pac);

  for (int k = w->start_pair; k < w->end_pair; k++) {
    align_pair_result_t* res = &w->results[k];
    memset(res, 0, sizeof(*res));
    if (w->paired) {
      int a = (k - w->start_pair) * 2, b = a + 1;
      int m1 = aln_fill_from_seq(&res->r1, &seqs[a]);
      int m2 = aln_fill_from_seq(&res->r2, &seqs[b]);
      if (m1 && m2) {
        res->mapped = 1;
        res->r1.proper_frag = res->r2.proper_frag = 1;
      } else {
        free(res->r1.cigar);
        free(res->r2.cigar);
        free(res->r1.md);
        free(res->r2.md);
        memset(res, 0, sizeof(*res));
        res->r1.pos = res->r2.pos = INT_MAX;
      }
    } else {
      int j = k - w->start_pair;
      res->mapped = aln_fill_from_seq(&res->r1, &seqs[j]);
    }
  }

  bwa_free_read_seq(n_seq, seqs);
  return NULL;
}

void bwa_aln_align_reads_direct(const char* reference,
                                const char* refName,
                                int n_reads,
                                const char** reads,
                                const char** quality,
                                const char** read_names,
                                const int* r_lens,
                                align_pair_result_t* results,
                                bool pairedEndReads,
                                int n_threads) {
  (void)read_names;
  if (n_reads <= 0 || !reference) return;

  bwa_verbose = 0;

  char tmpl[] = "/tmp/panmap_bwaaln_XXXXXX";
  int fd = mkstemp(tmpl);
  if (fd < 0) return;
  FILE* fa = fdopen(fd, "w");
  if (!fa) {
    close(fd);
    unlink(tmpl);
    return;
  }
  fprintf(fa, ">%s\n%s\n", refName, reference);
  fclose(fa);

  int build_rc = bwa_idx_build(tmpl, tmpl, 0, 10000000);
  bwaidx_t* idx = (build_rc == 0) ? bwa_idx_load(tmpl, BWA_IDX_ALL) : NULL;

  const char* exts[] = {"", ".amb", ".ann", ".bwt", ".pac", ".sa"};
  char path[4096];
  for (int e = 0; e < 6; e++) {
    snprintf(path, sizeof(path), "%s%s", tmpl, exts[e]);
    unlink(path);
  }
  if (!idx) return;

  bwase_initialize();  // populates g_log_n for mapQ; required because we bypass samse

  // Ancient-DNA defaults: -l 1024 disables the seed (damage clusters at read ends);
  // -n 0.01 / -o 2 per Oliva et al. 2021 (nf-core/eager); -q 0, trimming is upstream.
  gap_opt_t* opt = gap_init_opt();
  opt->fnr = 0.01f;      // -n 0.01
  opt->max_diff = -1;    // use fnr, not a fixed edit distance
  opt->max_gapo = 2;     // -o 2
  opt->seed_len = 1024;  // -l 1024
  opt->trim_qual = 0;    // -q 0
  opt->n_threads = 1;    // we chunk; bwa's own threading off

  if (n_threads < 1) n_threads = 1;
  int n_items = pairedEndReads ? n_reads / 2 : n_reads;
  if (n_threads > n_items) n_threads = n_items;
  if (n_threads < 1) n_threads = 1;

  bwa_aln_worker_t* workers = (bwa_aln_worker_t*)malloc(n_threads * sizeof(bwa_aln_worker_t));
  int chunk = n_items / n_threads, remainder = n_items % n_threads, offset = 0;
  for (int t = 0; t < n_threads; t++) {
    workers[t].idx = idx;
    workers[t].opt = opt;
    workers[t].reads = reads;
    workers[t].quals = quality;
    workers[t].r_lens = r_lens;
    workers[t].results = results;
    workers[t].paired = pairedEndReads ? 1 : 0;
    workers[t].start_pair = offset;
    workers[t].end_pair = offset + chunk + (t < remainder ? 1 : 0);
    offset = workers[t].end_pair;
  }

  if (n_threads == 1) {
    bwa_aln_worker_func(&workers[0]);
  } else {
    pthread_t* threads = (pthread_t*)malloc(n_threads * sizeof(pthread_t));
    char* created = (char*)calloc(n_threads, 1);
    for (int t = 0; t < n_threads; t++) {
      if (pthread_create(&threads[t], NULL, bwa_aln_worker_func, &workers[t]) == 0)
        created[t] = 1;
      else
        bwa_aln_worker_func(&workers[t]);
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