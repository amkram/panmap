#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include "3rdparty/minimap2/kalloc.h"
#include "3rdparty/minimap2/khash.h"
#include "3rdparty/minimap2/kseq.h"
#include "3rdparty/minimap2/minimap.h"
#include "3rdparty/minimap2/mmpriv.h"

mm128_t *collect_seed_hits_heap(void *km, const mm_mapopt_t *opt, int max_occ,
                                const mm_idx_t *mi, const char *qname,
                                const mm128_v *mv, int qlen, int64_t *n_a,
                                int *rep_len, int *n_mini_pos,
                                uint64_t **mini_pos);
mm128_t *collect_seed_hits(void *km, const mm_mapopt_t *opt, int max_occ,
                           const mm_idx_t *mi, const char *qname,
                           const mm128_v *mv, int qlen, int64_t *n_a,
                           int *rep_len, int *n_mini_pos, uint64_t **mini_pos);
void chain_post(const mm_mapopt_t *opt, int max_chain_gap_ref,
                const mm_idx_t *mi, void *km, int qlen, int n_segs,
                const int *qlens, int *n_regs, mm_reg1_t *regs, mm128_t *a);
mm_reg1_t *align_regs(const mm_mapopt_t *opt, const mm_idx_t *mi, void *km,
                      int qlen, const char *seq, int *n_regs, mm_reg1_t *regs,
                      mm128_t *a);

//#include "mm_align.h"

// mi: holds reference info and alignment flags and options
// num_reads
// read_lengths
// read_seqs: sequence of reads, as a C string of A C G T or U
// n_regs: output of function, number of alignment hits (per read)
// regs: output of function, array of alignment hits, must free (*regs[i])->p
// using free(), must free (*regs) using free() b: thread buffer, just
// initialize using mm_tbuf_init() before all alignments, then afterward destroy
// using mm_tbuf_destroy() opt: alignment options seed_l: length of seeds
// n_seeds_0: number of seed hits in first read
// n_seeds_1: number of seed hits in second read (or 0 if there only one read)
// r_poss: position of seeds(from both reads) in reference coords (the ends of
// the seed, the last index included in the seed) q_poss: position of seeds in
// query(read) coords (the ends of the seed, the last index included in the
// seed) reversed: whether seeds are reversed
void align_read_given_seeds(const mm_idx_t *mi, int num_reads,
                            const int *read_lengths, const char **read_seqs,
                            int *n_regs, mm_reg1_t **regs, mm_tbuf_t *b,
                            const mm_mapopt_t *opt, int seed_l, int n_seeds_0,
                            int n_seeds_1, int *r_poss, int *q_poss,
                            uint8_t *reversed) {
  int n_seeds = n_seeds_0 + n_seeds_1;
  int i, j, rep_len = 0, qlen_sum = 0, n_regs0, n_mini_pos = 0;
  int max_chain_gap_qry, max_chain_gap_ref,
      is_splice = !!(opt->flag & MM_F_SPLICE), is_sr = !!(opt->flag & MM_F_SR);
  uint32_t hash;
  int64_t n_a;
  uint64_t *u, *mini_pos;
  mm128_t *a;
  mm128_v mv = {0, 0, 0};
  mm_reg1_t *regs0;
  km_stat_t kmst;
  float chn_pen_gap, chn_pen_skip;

  for (i = 0; i < num_reads; i++) {
    n_regs[i] = 0;
    regs[i] = 0;
    qlen_sum += read_lengths[i];
  }

  if (opt->max_qlen > 0 && qlen_sum > opt->max_qlen)
    return;

  hash = __ac_Wang_hash(qlen_sum) + __ac_Wang_hash(opt->seed);
  hash = __ac_Wang_hash(hash);

  a = (mm128_t *)kmalloc(b->km, n_seeds * sizeof(mm128_t));

  for (i = 0; i < n_seeds; ++i) {
    uint32_t qpos = (uint32_t)
        q_poss[i]; // query position, position on read coordinates, position of
                   // the END of the seed, including it. so if this is 6 and
                   // read[6] is A then A is the last character of the seed
    uint32_t rpos =
        (uint32_t)r_poss[i]; // reference position of the END of the seed
    uint32_t span = (uint32_t)seed_l; // seed length

    uint64_t x, y;

    x = ((uint64_t)reversed[i] << 63) | rpos;
    y = (uint64_t)span << 32 | qpos;

    a[i].x = x;
    a[i].y = y;
  }
  n_a = n_seeds;

  if (num_reads ==
      2) { // TODO check if this needs to be done non-paired end as well
    for (i = 0; i < n_seeds_0; i++) {
      if (reversed[i]) {
        a[i].y += read_lengths[1];
      }
    }
  }

  for (i = n_seeds_0; i < n_seeds; i++) {
    a[i].y |= 0x1000000000000;
    if (!reversed[i]) {
      a[i].y += read_lengths[0];
    }
  }

  if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
    fprintf(stderr, "RS\t%d\n", rep_len);
    for (i = 0; i < n_a; ++i)
      fprintf(stderr, "SD\t%s\t%d\t%c\t%d\t%d\t%d\t%lx\t%lx\n",
              mi->seq[a[i].x << 1 >> 33].name, (int32_t)a[i].x,
              "+-"[a[i].x >> 63], (int32_t)a[i].y,
              (int32_t)(a[i].y >> 32 & 0xff),
              i == 0 ? 0
                     : ((int32_t)a[i].y - (int32_t)a[i - 1].y) -
                           ((int32_t)a[i].x - (int32_t)a[i - 1].x),
              a[i].x, a[i].y);
  }

  // set max chaining gap on the query and the reference sequence
  if (is_sr) {
    max_chain_gap_qry = qlen_sum > opt->max_gap ? qlen_sum : opt->max_gap;
  } else
    max_chain_gap_qry = opt->max_gap;
  if (opt->max_gap_ref > 0) {
    max_chain_gap_ref =
        opt->max_gap_ref; // always honor mm_mapopt_t::max_gap_ref if set
  } else if (opt->max_frag_len > 0) {
    max_chain_gap_ref = opt->max_frag_len - qlen_sum;
    if (max_chain_gap_ref < opt->max_gap)
      max_chain_gap_ref = opt->max_gap;
  } else
    max_chain_gap_ref = opt->max_gap;

  chn_pen_gap = opt->chain_gap_scale * 0.01 * mi->k;
  chn_pen_skip = opt->chain_skip_scale * 0.01 * mi->k;
  if (opt->flag & MM_F_RMQ) {
    a = mg_lchain_rmq(opt->max_gap, opt->rmq_inner_dist, opt->bw,
                      opt->max_chain_skip, opt->rmq_size_cap, opt->min_cnt,
                      opt->min_chain_score, chn_pen_gap, chn_pen_skip, n_seeds,
                      a, &n_regs0, &u, b->km);
  } else {
    a = mg_lchain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw,
                     opt->max_chain_skip, opt->max_chain_iter, opt->min_cnt,
                     opt->min_chain_score, chn_pen_gap, chn_pen_skip, is_splice,
                     num_reads, n_seeds, a, &n_regs0, &u, b->km);
  }

  if (opt->bw_long > opt->bw &&
      (opt->flag & (MM_F_SPLICE | MM_F_SR | MM_F_NO_LJOIN)) == 0 &&
      num_reads == 1 && n_regs0 > 1) { // re-chain/long-join for long sequences
    int32_t st = (int32_t)a[0].y, en = (int32_t)a[(int32_t)u[0] - 1].y;
    if (qlen_sum - (en - st) > opt->rmq_rescue_size ||
        en - st > qlen_sum * opt->rmq_rescue_ratio) {
      int32_t i;
      for (i = 0, n_seeds = 0; i < n_regs0; ++i)
        n_seeds += (int32_t)u[i];
      kfree(b->km, u);
      radix_sort_128x(a, a + n_seeds);
      a = mg_lchain_rmq(opt->max_gap, opt->rmq_inner_dist, opt->bw_long,
                        opt->max_chain_skip, opt->rmq_size_cap, opt->min_cnt,
                        opt->min_chain_score, chn_pen_gap, chn_pen_skip,
                        n_seeds, a, &n_regs0, &u, b->km);
    }
  } else if (opt->max_occ > opt->mid_occ && rep_len > 0 &&
             !(opt->flag & MM_F_RMQ)) { // re-chain, mostly for short reads
    int rechain = 0;
    if (n_regs0 > 0) { // test if the best chain has all the segments
      int n_chained_segs = 1, max = 0, max_i = -1, max_off = -1, off = 0;
      for (i = 0; i < n_regs0; ++i) { // find the best chain
        if (max < (int)(u[i] >> 32))
          max = u[i] >> 32, max_i = i, max_off = off;
        off += (uint32_t)u[i];
      }
      for (i = 1; i < (int32_t)u[max_i];
           ++i) { // count the number of segments in the best chain
        if ((a[max_off + i].y & MM_SEED_SEG_MASK) !=
            (a[max_off + i - 1].y & MM_SEED_SEG_MASK))
          ++n_chained_segs;
      }
      if (n_chained_segs < num_reads)
        rechain = 1;
    } else
      rechain = 1;
    if (rechain) { // redo chaining with a higher max_occ threshold
      kfree(b->km, a);
      kfree(b->km, u);
      kfree(b->km, mini_pos);
      if (opt->flag & MM_F_HEAP_SORT)
        a = collect_seed_hits_heap(b->km, opt, opt->max_occ, mi, NULL, &mv,
                                   qlen_sum, &n_seeds, &rep_len, &n_mini_pos,
                                   &mini_pos);
      else
        a = collect_seed_hits(b->km, opt, opt->max_occ, mi, NULL, &mv, qlen_sum,
                              &n_seeds, &rep_len, &n_mini_pos, &mini_pos);
      a = mg_lchain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw,
                       opt->max_chain_skip, opt->max_chain_iter, opt->min_cnt,
                       opt->min_chain_score, chn_pen_gap, chn_pen_skip,
                       is_splice, num_reads, n_seeds, a, &n_regs0, &u, b->km);
    }
  }
  b->frag_gap = max_chain_gap_ref;
  b->rep_len = rep_len;

  regs0 = mm_gen_regs(b->km, hash, qlen_sum, n_regs0, u, a,
                      !!(opt->flag & MM_F_QSTRAND));
  if (mi->n_alt) {
    mm_mark_alt(mi, n_regs0, regs0);
    mm_hit_sort(b->km, &n_regs0, regs0,
                opt->alt_drop); // this step can be merged into mm_gen_regs();
                                // will do if this shows up in profile
  }

  chain_post(opt, max_chain_gap_ref, mi, b->km, qlen_sum, num_reads,
             read_lengths, &n_regs0, regs0, a);
  if (!is_sr && !(opt->flag & MM_F_QSTRAND)) {
    mm_est_err(mi, qlen_sum, n_regs0, regs0, a, n_mini_pos, mini_pos);
    n_regs0 = mm_filter_strand_retained(n_regs0, regs0);
  }

  if (num_reads == 1) { // uni-segment
    regs0 = align_regs(opt, mi, b->km, read_lengths[0], read_seqs[0], &n_regs0,
                       regs0, a);
    regs0 = (mm_reg1_t *)realloc(regs0, sizeof(*regs0) * n_regs0);
    mm_set_mapq(b->km, n_regs0, regs0, opt->min_chain_score, opt->a, rep_len,
                is_sr);
    n_regs[0] = n_regs0, regs[0] = regs0;
  } else { // multi-segment
    mm_seg_t *seg;
    seg =
        mm_seg_gen(b->km, hash, num_reads, read_lengths, n_regs0, regs0, n_regs,
                   regs, a); // split fragment chain to separate segment chains
    free(regs0);

    for (i = 0; i < num_reads; ++i) {
      mm_set_parent(b->km, opt->mask_level, opt->mask_len, n_regs[i], regs[i],
                    opt->a * 2 + opt->b, opt->flag & MM_F_HARD_MLEVEL,
                    opt->alt_drop); // update mm_reg1_t::parent
      regs[i] = align_regs(opt, mi, b->km, read_lengths[i], read_seqs[i],
                           &n_regs[i], regs[i], seg[i].a);
      mm_set_mapq(b->km, n_regs[i], regs[i], opt->min_chain_score, opt->a,
                  rep_len, is_sr);
    }

    mm_seg_free(b->km, num_reads, seg);
    if (num_reads == 2 && opt->pe_ori >= 0 && (opt->flag & MM_F_CIGAR)) {
      mm_pair(b->km, max_chain_gap_ref, opt->pe_bonus, opt->a * 2 + opt->b,
              opt->a, read_lengths, n_regs, regs); // pairing
    }
  }

  // kfree(b->km, mv.a);
  kfree(b->km, a);
  kfree(b->km, u);
  // kfree(b->km, mini_pos);

  if (b->km) {
    km_stat(b->km, &kmst);
    assert(kmst.n_blocks == kmst.n_cores); // otherwise, there is a memory leak
    if (kmst.largest > 1U << 28 ||
        (opt->cap_kalloc > 0 && kmst.capacity > opt->cap_kalloc)) {
      if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
        fprintf(stderr, "[W::%s] reset thread-local memory after read %s\n",
                __func__, (char *)NULL);
      km_destroy(b->km);
      b->km = km_init();
    }
  }
}

// Stores sam strings in sam_alignments, and stores reference positions of
// alignments in r_lens
void align_reads(const char *reference, int n_reads, const char **reads,
                 const char **quality, const char **read_names, int *r_lens,
                 int *seed_counts, uint8_t **reversed, int **ref_positions,
                 int **qry_positions, char **sam_alignments, int syncmer_k,
                 bool pairedEndReads) {

  mm_idxopt_t iopt;
  mm_mapopt_t mopt;
  int n_threads = 1;

  mm_verbose = 2; // disable message output to stderr
  mm_set_opt(0, &iopt, &mopt);
  mopt.flag |= MM_F_CIGAR; // perform alignment
  // mopt.flag |= MM_F_SR;
  // mopt.max_gap = 100;
  // mopt.max_frag_len = 800;
  // mopt.max_gap = 100;

  iopt.k = syncmer_k; // setting seed length

  mm_idx_t *mi = mm_idx_str(10, iopt.k, 0, 14, 1, &reference,
                            NULL); // Read reference into structure

  mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence;
  mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading,
                                    // allocate one tbuf for each thread

  if (pairedEndReads) {
    for (int k = 0; k < n_reads / 2; k++) {

      mm_reg1_t *reg[2] = {NULL, NULL};
      int n_reg[2] = {0, 0};
      align_read_given_seeds(mi, 2, &(r_lens[k * 2]), &(reads[k * 2]), n_reg,
                             reg, tbuf, &mopt, iopt.k, seed_counts[k * 2],
                             seed_counts[k * 2 + 1], ref_positions[k],
                             qry_positions[k], reversed[k]);
      mi->seq[0].name = "ref";

      if (n_reg[0] == 0 || n_reg[1] == 0 || reg[0]->score <= 0 ||
          reg[1]->score <= 0) {
        sam_alignments[k * 2] = NULL;
        sam_alignments[k * 2 + 1] = NULL;
        r_lens[k * 2] = INT_MAX;
        r_lens[k * 2 + 1] = INT_MAX;
      } else {

        kstring_t sam = {0, 0, 0};
        mm_bseq1_t t;
        t.l_seq = r_lens[k * 2];
        t.seq = reads[k * 2];
        t.name = read_names[k * 2];
        t.qual = quality[k * 2];
        mm_write_sam3(&sam, mi, &t, 0, 0, 2, n_reg, reg, NULL, 0, 0);
        r_lens[k * 2] = reg[0]->rs + 1; // Length of sam line string
        sam_alignments[k * 2] = sam.s;

        kstring_t sam2 = {0, 0, 0};
        t.l_seq = r_lens[k * 2 + 1];
        t.seq = reads[k * 2 + 1];
        t.name = read_names[k * 2 + 1];
        t.qual = quality[k * 2 + 1];
        mm_write_sam3(&sam2, mi, &t, 1, 0, 2, n_reg, reg, NULL, 0, 0);
        r_lens[k * 2 + 1] = reg[1]->rs + 1; // Length of sam line string
        sam_alignments[k * 2 + 1] = sam2.s;
      }

      for (int j = 0; j < n_reg[0]; ++j) {
        mm_reg1_t *r = &reg[0][j];
        assert(r->p); // with MM_F_CIGAR, this should not be NULL
        free(r->p);
      }
      for (int j = 0; j < n_reg[1]; ++j) {
        mm_reg1_t *r = &reg[1][j];
        assert(r->p); // with MM_F_CIGAR, this should not be NULL
        free(r->p);
      }
      free(reg[0]);
      free(reg[1]);
    }
  } else {

    for (int k = 0; k < n_reads; k++) {

      mm_reg1_t *reg[1] = {NULL}; // Alignments
      int n_reg[1] = {0};

      align_read_given_seeds(mi, 1, &(r_lens[k]), &(reads[k]), n_reg, reg, tbuf,
                             &mopt, iopt.k, seed_counts[k], 0, ref_positions[k],
                             qry_positions[k], reversed[k]);

      mi->seq[0].name = "ref"; // TODO use actual ref name

      // reg[0]->score is the number of quality aligned bp

      if (n_reg[0] == 0 || reg[0]->score <= 0 || reg[0]->score > r_lens[k]) {
        sam_alignments[k] = NULL;
        r_lens[k] = INT_MAX;
      } else {

        kstring_t sam = {0, 0, 0};
        mm_bseq1_t t;
        t.l_seq = r_lens[k];
        t.seq = reads[k];
        t.name = read_names[k];
        t.qual = quality[k];

        mm_write_sam3(&sam, mi, &t, 0, 0, 1, n_reg, reg, NULL, 0, 0);
        r_lens[k] = reg[0]->rs + 1; // Length of sam line string
        sam_alignments[k] = sam.s;
      }

      for (int j = 0; j < n_reg[0]; ++j) {
        mm_reg1_t *r = &reg[0][j];
        assert(r->p); // with MM_F_CIGAR, this should not be NULL

        free(r->p);
      }
      free(reg[0]);
    }
  }
  mm_tbuf_destroy(tbuf);
  mm_idx_destroy(mi);
}