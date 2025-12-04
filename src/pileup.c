#include "pileup.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#define dummy_free(p)
KLIST_INIT(auxlist, char *, dummy_free)

char *global_ref_string;
hts_pos_t global_ref_length;
int bam_index;
int bams;
bam1_t **bam_records;

int our_printw(int c, kstring_t *k) {
  char buf[16];
  int l, x;
  if (c == 0) {
    // fputc('0', fp);
    return kputc('0', k);
  }
  for (l = 0, x = c < 0 ? -c : c; x > 0; x /= 10)
    buf[l++] = x % 10 + '0';
  if (c < 0)
    buf[l++] = '-';
  buf[l] = 0;
  for (x = 0; x < l / 2; ++x) {
    int y = buf[x];
    buf[x] = buf[l - 1 - x];
    buf[l - 1 - x] = y;
  }
  kputs(buf, k);
  // fputs(buf, fp);
  return 0;
}

int our_mplp_get_ref(mplp_aux_t *ma, int tid, char **ref, hts_pos_t *ref_len) {
  mplp_ref_t *r = ma->ref;

  // printf("get ref %d {%d/%p, %d/%p}\n", tid, r->ref_id[0], r->ref[0],
  // r->ref_id[1], r->ref[1]);

  if (!r) {
    *ref = NULL;
    return 0;
  }

  // Do we need to reference count this so multiple mplp_aux_t can
  // track which references are in use?
  // For now we just cache the last two. Sufficient?
  if (tid == r->ref_id[0]) {
    *ref = r->ref[0];
    *ref_len = r->ref_len[0];
    return 1;
  }

  if (tid == r->ref_id[1]) {

    // Last, swap over
    int tmp_id;
    hts_pos_t tmp_len;
    tmp_id = r->ref_id[0];
    r->ref_id[0] = r->ref_id[1];
    r->ref_id[1] = tmp_id;
    tmp_len = r->ref_len[0];
    r->ref_len[0] = r->ref_len[1];
    r->ref_len[1] = tmp_len;

    char *tc;
    tc = r->ref[0];
    r->ref[0] = r->ref[1];
    r->ref[1] = tc;
    *ref = r->ref[0];
    *ref_len = r->ref_len[0];
    return 1;
  }

  // New, so migrate to old and load new
  free(r->ref[1]);
  r->ref[1] = r->ref[0];
  r->ref_id[1] = r->ref_id[0];
  r->ref_len[1] = r->ref_len[0];

  r->ref_id[0] = tid;

  r->ref[0] = global_ref_string; // Ugly use of global variables
  r->ref_len[0] = global_ref_length;

  if (!r->ref[0]) {
    r->ref[0] = NULL;
    r->ref_id[0] = -1;
    r->ref_len[0] = 0;
    *ref = NULL;
    return 0;
  }

  *ref = r->ref[0];
  *ref_len = r->ref_len[0];
  return 1;
}

int our_pileup_seq(kstring_t *mplp_string, const bam_pileup1_t *p,
                   hts_pos_t pos, hts_pos_t ref_len, const char *ref,
                   kstring_t *ks, int rev_del, int no_ins, int no_ins_mods,
                   int no_del, int no_ends) {
  no_ins_mods |= no_ins;
  int j;
  hts_base_mod_state *m = p->cd.p;
  if (!no_ends && p->is_head) {
    kputc('^', mplp_string);
    // putc('^', fp);
    kputc(p->b->core.qual > 93 ? 126 : p->b->core.qual + 33, mplp_string);
    // putc(p->b->core.qual > 93? 126 : p->b->core.qual + 33, fp);
  }
  if (!p->is_del) {
    int c = p->qpos < p->b->core.l_qseq
                ? seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos)]
                : 'N';
    if (ref) {
      int rb = pos < ref_len ? ref[pos] : 'N';
      if (c == '=' || seq_nt16_table[c] == seq_nt16_table[rb])
        c = bam_is_rev(p->b) ? ',' : '.';
      else
        c = bam_is_rev(p->b) ? tolower(c) : toupper(c);
    } else {
      if (c == '=')
        c = bam_is_rev(p->b) ? ',' : '.';
      else
        c = bam_is_rev(p->b) ? tolower(c) : toupper(c);
    }
    kputc(c, mplp_string);
    // putc(c, fp);
    if (m) {
      int nm;
      hts_base_mod mod[256];
      if ((nm = bam_mods_at_qpos(p->b, p->qpos, m, mod, 256)) > 0) {
        kputc('[', mplp_string);
        // putc('[', fp);
        int j;
        for (j = 0; j < nm && j < 256; j++) {
          char qual[20];
          if (mod[j].qual >= 0)
            sprintf(qual, "%d", mod[j].qual);
          else
            *qual = 0;
          if (mod[j].modified_base < 0) {
            // ChEBI
            ksprintf(mplp_string, "%c(%d)%s", "+-"[mod[j].strand],
                     -mod[j].modified_base, qual);
            // fprintf(fp, "%c(%d)%s", "+-"[mod[j].strand],
            //-mod[j].modified_base, qual);
          } else {
            ksprintf(mplp_string, "%c%c%s", "+-"[mod[j].strand],
                     mod[j].modified_base, qual);
            // fprintf(fp, "%c%c%s", "+-"[mod[j].strand],
            // mod[j].modified_base, qual);
          }
        }
        kputc(']', mplp_string);
        // putc(']', fp);
      }
    }
  } else {
    kputc(p->is_refskip ? (bam_is_rev(p->b) ? '<' : '>')
                        : ((bam_is_rev(p->b) && rev_del) ? '#' : '*'),
          mplp_string);
    // putc(p->is_refskip? (bam_is_rev(p->b)? '<' : '>') : ((bam_is_rev(p->b) &&
    // rev_del) ? '#' : '*'), fp);
  }
  int del_len = -p->indel;
  if (p->indel > 0) {
    int len =
        bam_plp_insertion_mod(p, m && !no_ins_mods ? m : NULL, ks, &del_len);
    if (len < 0) {
      print_error("our_mpileup", "bam_plp_insertion() failed");
      return -1;
    }
    if (no_ins < 2) {
      kputc('+', mplp_string);
      // putc('+', fp);
      our_printw(len, mplp_string);
    }
    if (!no_ins) {
      if (bam_is_rev(p->b)) {
        char pad = rev_del ? '#' : '*';
        int in_mod = 0;
        for (j = 0; j < ks->l; j++) {
          if (ks->s[j] == '[')
            in_mod = 1;
          else if (ks->s[j] == ']')
            in_mod = 0;
          kputc(ks->s[j] != '*' ? (in_mod ? ks->s[j] : tolower(ks->s[j])) : pad,
                mplp_string);
          // putc(ks->s[j] != '*'
          //? (in_mod ? ks->s[j] : tolower(ks->s[j]))
          //: pad, fp);
        }
      } else {
        int in_mod = 0;
        for (j = 0; j < ks->l; j++) {
          if (ks->s[j] == '[')
            in_mod = 1;
          if (ks->s[j] == ']')
            in_mod = 0;
          kputc(in_mod ? ks->s[j] : toupper(ks->s[j]), mplp_string);
          // putc(in_mod ? ks->s[j] : toupper(ks->s[j]), fp);
        }
      }
    }
  }
  if (del_len > 0) {
    if (no_del < 2)
      our_printw(-del_len, mplp_string);
    if (!no_del) {
      for (j = 1; j <= del_len; ++j) {
        int c = (ref && (int)pos + j < ref_len) ? ref[pos + j] : 'N';
        kputc(bam_is_rev(p->b) ? tolower(c) : toupper(c), mplp_string);
        // putc(bam_is_rev(p->b)? tolower(c) : toupper(c), fp);
      }
    }
  }
  if (!no_ends && p->is_tail) {
    kputc('$', mplp_string);
    // putc('$', fp);
  }
  return 0;
}

int our_pileup_cd_create(void *data, const bam1_t *b, bam_pileup_cd *cd) {
  int ret;
  hts_base_mod_state *m = hts_base_mod_state_alloc();
  ret = bam_parse_basemod(b, m);
  cd->p = m;
  return ret;
}

int our_pileup_cd_destroy(void *data, const bam1_t *b, bam_pileup_cd *cd) {
  hts_base_mod_state_free(cd->p);
  return 0;
}

void our_print_empty_pileup(kstring_t *k, const mplp_conf_t *conf,
                            const char *tname, hts_pos_t pos, int n,
                            const char *ref, hts_pos_t ref_len) {
  int i;
  ksprintf(k, "%s\t%" PRIhts_pos "\t%c", tname, pos + 1,
           (ref && pos < ref_len) ? ref[pos] : 'N');
  // fprintf(fp, "%s\t%"PRIhts_pos"\t%c", tname, pos+1, (ref && pos < ref_len)?
  // ref[pos] : 'N');
  for (i = 0; i < n; ++i) {
    kputs("\t0\t*\t*", k);
    // fputs("\t0\t*\t*", fp);
    int flag_value = MPLP_PRINT_MAPQ_CHAR;
    while (flag_value < MPLP_PRINT_LAST) {
      if (flag_value != MPLP_PRINT_MODS && (conf->flag & flag_value)) {
        kputs("\t*", k);
        // fputs("\t*", fp);
      }
      flag_value <<= 1;
    }
    if (conf->auxlist) {
      int t = 0;
      while (t++ < ((klist_t(auxlist) *)conf->auxlist)->size) {
        kputs("\t*", k);
        // fputs("\t*", fp);
      }
    }
  }
  kputs("\n", k);
  // putc('\n', fp);
}
// I think this function populates the bam1_t structure.

int our_mplp_func(void *data, bam1_t *b) {
  char *ref;
  mplp_aux_t *ma = (mplp_aux_t *)data;
  int ret, skip = 0;
  hts_pos_t ref_len;

  do {
    int has_ref;

    if (bam_index >= bams) {
      ret = -1;
    } else {
      ret = 1;
      (void)bam_copy1(b, bam_records[bam_index]);
    }
    bam_index += 1;

    if (ret < 0)
      break;

    // The 'B' cigar operation is not part of the specification, considering as
    // obsolete.
    //  bam_remove_B(b);
    if (b->core.tid < 0 ||
        (b->core.flag & BAM_FUNMAP)) { // exclude unmapped reads
      skip = 1;
      continue;
    }

    if (ma->conf->rflag_require && !(ma->conf->rflag_require & b->core.flag)) {
      skip = 1;
      continue;
    }

    if (ma->conf->rflag_filter && ma->conf->rflag_filter & b->core.flag) {
      skip = 1;
      continue;
    }

    if (ma->conf->bed && ma->conf->all == 0) { // test overlap
      skip = !bed_overlap(ma->conf->bed, sam_hdr_tid2name(ma->h, b->core.tid),
                          b->core.pos, bam_endpos(b));
      if (skip)
        continue;
    }

    if (ma->conf->rghash) { // exclude read groups
      uint8_t *rg = bam_aux_get(b, "RG");
      skip = (rg && khash_str2int_get(ma->conf->rghash, (const char *)(rg + 1),
                                      NULL) == 0);
      if (skip)
        continue;
    }

    if (ma->conf->flag & MPLP_ILLUMINA13) {
      int i;
      uint8_t *qual = bam_get_qual(b);
      for (i = 0; i < b->core.l_qseq; ++i)
        qual[i] = qual[i] > 31 ? qual[i] - 31 : 0;
    }

    if (b->core.tid >= 0) { // changed this line
      has_ref = our_mplp_get_ref(ma, b->core.tid, &ref, &ref_len);
      if (has_ref &&
          ref_len <=
              b->core.pos) { // exclude reads outside of the reference sequence
        fprintf(stderr,
                "[%s] Skipping because %" PRIhts_pos
                " is outside of %" PRIhts_pos " [ref:%d]\n",
                __func__, (int64_t)b->core.pos, ref_len, b->core.tid);
        skip = 1;
        continue;
      }
    } else {
      has_ref = 0;
    }

    skip = 0;
    if (has_ref && (ma->conf->flag & MPLP_REALN)) {

      sam_prob_realn(b, ref, ref_len, (ma->conf->flag & MPLP_REDO_BAQ) ? 7 : 3);
    }
    if (has_ref && ma->conf->capQ_thres > 10) {
      int q = sam_cap_mapq(b, ref, ref_len, ma->conf->capQ_thres);
      if (q < 0)
        skip = 1;
      else if (b->core.qual > q)
        b->core.qual = q;
    }
    if (b->core.qual < ma->conf->min_mq)
      skip = 1;
    else if ((ma->conf->flag & MPLP_NO_ORPHAN) &&
             (b->core.flag & BAM_FPAIRED) && !(b->core.flag & BAM_FPROPER_PAIR))
      skip = 1;
  } while (skip);
  return ret;
}

/*
 * Performs pileup
 * @param conf configuration for this pileup
 * @param n number of files specified in fn
 * @param fn filenames
 * @param fn_idx index filenames
 */
int our_mpileup(mplp_conf_t *conf, sam_hdr_t *sam_header,
                kstring_t *mplp_string) {

  mplp_aux_t **data;
  int i, tid, *n_plp, tid0 = 0, max_depth, n = 1;
  hts_pos_t pos, beg0 = 0, end0 = HTS_POS_MAX, ref_len;
  const bam_pileup1_t **plp;
  mplp_ref_t mp_ref = MPLP_REF_INIT;

  bam_mplp_t iter;
  sam_hdr_t *h = NULL; /* header of first file in input list */
  char *ref;
  // FILE *pileup_fp = NULL;

  kstring_t buf;
  mplp_pileup_t gplp;

  memset(&gplp, 0, sizeof(mplp_pileup_t));
  memset(&buf, 0, sizeof(kstring_t));
  data = calloc(n, sizeof(mplp_aux_t *));
  plp = calloc(n, sizeof(bam_pileup1_t *));
  n_plp = calloc(n, sizeof(int));

  if (n == 0) {
    fprintf(stderr, "[%s] no input file/data given\n", __func__);
    exit(EXIT_FAILURE);
  }

  // read the header of each file in the list and initialize  data
  refs_t *refs = NULL;
  for (i = 0; i < n; ++i) {
    sam_hdr_t *h_tmp = sam_header; // Providing sam header
    data[i] = calloc(1, sizeof(mplp_aux_t));

    data[i]->conf = conf;
    data[i]->ref = &mp_ref;

    // bam_smpl_add(sm, fn[i], (conf->flag&MPLP_IGNORE_RG)? 0 :
    // sam_hdr_str(h_tmp));

    data[i]->iter = NULL;

    if (i == 0)
      h = data[i]->h = h_tmp; // save the header of the first file
    else {
      sam_hdr_destroy(h_tmp);

      // we store only the first file's header; it's (alleged to be)
      // compatible with the i-th file's target_name lookup needs
      data[i]->h = h;
    }
  }

  iter = bam_mplp_init(n, our_mplp_func, (void **)data);

  if (conf->flag & MPLP_PRINT_MODS) {
    bam_mplp_constructor(iter, our_pileup_cd_create);
    bam_mplp_destructor(iter, our_pileup_cd_destroy);
  }
  if (conf->flag & MPLP_SMART_OVERLAPS)
    bam_mplp_init_overlaps(iter);
  if (!conf->max_depth) {
    max_depth = INT_MAX;
    fprintf(stderr, "[%s] Max depth set to maximum value (%d)\n", __func__,
            INT_MAX);
  } else {
    max_depth = conf->max_depth;
    if (max_depth * n > 1 << 20)
      fprintf(stderr,
              "[%s] Combined max depth is above 1M. Potential memory hog!\n",
              __func__);
  }

  bam_mplp_set_maxcnt(iter, max_depth);
  int ret;
  int last_tid = -1;
  hts_pos_t last_pos = -1;
  int one_seq = 0;

  // begin pileup
  while ((ret = bam_mplp64_auto(iter, &tid, &pos, n_plp, plp)) > 0) {

    one_seq = 1; // at least 1 output
    if (conf->reg && (pos < beg0 || pos >= end0))
      continue; // out of the region requested

    our_mplp_get_ref(data[0], tid, &ref, &ref_len);

    if (conf->all) {
      // Deal with missing portions of previous tids
      while (tid > last_tid) {
        if (last_tid >= 0 && !conf->reg) {
          while (++last_pos < sam_hdr_tid2len(h, last_tid)) {
            if (conf->bed &&
                bed_overlap(conf->bed, sam_hdr_tid2name(h, last_tid), last_pos,
                            last_pos + 1) == 0)
              continue;
            // TODO handle empty case
            our_print_empty_pileup(mplp_string, conf,
                                   sam_hdr_tid2name(h, last_tid), last_pos, n,
                                   ref, ref_len);
          }
        }
        last_tid++;
        last_pos = -1;
        if (conf->all < 2)
          break;
      }
    }
    if (conf->all) {
      // Deal with missing portion of current tid
      while (++last_pos < pos) {
        if (conf->reg && last_pos < beg0)
          continue; // out of range; skip
        if (conf->bed && bed_overlap(conf->bed, sam_hdr_tid2name(h, tid),
                                     last_pos, last_pos + 1) == 0)
          continue;
        // TODO handle empty case
        our_print_empty_pileup(mplp_string, conf, sam_hdr_tid2name(h, tid),
                               last_pos, n, ref, ref_len);
      }
      last_tid = tid;
      last_pos = pos;
    }
    if (conf->bed && tid >= 0 &&
        !bed_overlap(conf->bed, sam_hdr_tid2name(h, tid), pos, pos + 1))
      continue;
    ksprintf(mplp_string, "%s\t%" PRIhts_pos "\t%c", sam_hdr_tid2name(h, tid),
             pos + 1, (ref && pos < ref_len) ? ref[pos] : 'N');
    // fprintf(pileup_fp, "%s\t%"PRIhts_pos"\t%c", sam_hdr_tid2name(h, tid), pos
    // + 1, (ref && pos < ref_len)? ref[pos] : 'N');
    for (i = 0; i < n; ++i) {
      int j, cnt;
      for (j = cnt = 0; j < n_plp[i]; ++j) {
        const bam_pileup1_t *p = plp[i] + j;
        int c = p->qpos < p->b->core.l_qseq ? bam_get_qual(p->b)[p->qpos] : 0;
        if (c >= conf->min_baseQ)
          ++cnt;
      }
      ksprintf(mplp_string, "\t%d\t", cnt);
      // fprintf(pileup_fp, "\t%d\t", cnt);
      if (n_plp[i] == 0) {
        kputs("*\t*", mplp_string);
        // fputs("*\t*", pileup_fp);
        int flag_value = MPLP_PRINT_MAPQ_CHAR;
        while (flag_value < MPLP_PRINT_LAST) {
          if (flag_value != MPLP_PRINT_MODS && (conf->flag & flag_value)) {
            kputs("\t*", mplp_string);
            // fputs("\t*", pileup_fp);
          }
          flag_value <<= 1;
        }
        if (conf->auxlist) {
          int t = 0;
          while (t++ < ((klist_t(auxlist) *)conf->auxlist)->size) {
            kputs("\t*", mplp_string);
            // fputs("\t*", pileup_fp);
          }
        }
      } else {
        int n = 0;
        kstring_t ks = KS_INITIALIZE;
        for (j = 0; j < n_plp[i]; ++j) {
          const bam_pileup1_t *p = plp[i] + j;
          int c = p->qpos < p->b->core.l_qseq ? bam_get_qual(p->b)[p->qpos] : 0;
          if (c >= conf->min_baseQ) {
            n++;
            // TODO what do i do
            if (our_pileup_seq(mplp_string, plp[i] + j, pos, ref_len, ref, &ks,
                               conf->rev_del, conf->no_ins, conf->no_ins_mods,
                               conf->no_del, conf->no_ends) < 0) {
              ret = 1;
              goto fail;
            }
          }
        }
        if (!n) {
          kputc('*', mplp_string);
          // putc('*', pileup_fp);
        }

        /* Print base qualities */
        n = 0;
        ks_free(&ks);
        kputc('\t', mplp_string);
        // putc('\t', pileup_fp);

        for (j = 0; j < n_plp[i]; ++j) {
          const bam_pileup1_t *p = plp[i] + j;
          int c = p->qpos < p->b->core.l_qseq ? bam_get_qual(p->b)[p->qpos] : 0;
          if (c >= conf->min_baseQ) {
            c = c + 33 < 126 ? c + 33 : 126;
            kputc(c, mplp_string);
            // putc(c, pileup_fp);
            n++;
          }
        }
        if (!n) {
          kputc('*', mplp_string);
          // putc('*', pileup_fp);
        }

        /* Print selected columns */
        int flag_value = MPLP_PRINT_MAPQ_CHAR;
        while (flag_value < MPLP_PRINT_LAST) {
          if (flag_value != MPLP_PRINT_MODS && (conf->flag & flag_value)) {
            n = 0;
            kputc('\t', mplp_string);
            // putc('\t', pileup_fp);
            for (j = 0; j < n_plp[i]; ++j) {
              const bam_pileup1_t *p = &plp[i][j];
              int c =
                  p->qpos < p->b->core.l_qseq ? bam_get_qual(p->b)[p->qpos] : 0;
              if (c < conf->min_baseQ)
                continue;
              if (n > 0 && flag_value != MPLP_PRINT_MAPQ_CHAR) {
                kputc(',', mplp_string);
                // putc(',', pileup_fp);
              }
              n++;

              switch (flag_value) {
              case MPLP_PRINT_MAPQ_CHAR:
                c = p->b->core.qual + 33;
                if (c > 126)
                  c = 126;
                kputc(c, mplp_string);
                // putc(c, pileup_fp);
                break;
              case MPLP_PRINT_QPOS:
                // query position in current orientation
                ksprintf(mplp_string, "%d", p->qpos + 1);
                // fprintf(pileup_fp, "%d", p->qpos + 1);
                break;
              case MPLP_PRINT_QPOS5: {
                // query position in 5' to 3' orientation
                int pos5 = bam_is_rev(p->b)
                               ? p->b->core.l_qseq - p->qpos + p->is_del
                               : p->qpos + 1;
                ksprintf(mplp_string, "%d", pos5);
                // fprintf(pileup_fp, "%d", pos5);
                break;
              }
              case MPLP_PRINT_QNAME:
                kputs(bam_get_qname(p->b), mplp_string);
                // fputs(bam_get_qname(p->b), pileup_fp);
                break;
              case MPLP_PRINT_FLAG:
                ksprintf(mplp_string, "%d", p->b->core.flag);
                // fprintf(pileup_fp, "%d", p->b->core.flag);
                break;
              case MPLP_PRINT_RNAME:
                if (p->b->core.tid >= 0) {
                  kputs(sam_hdr_tid2name(h, p->b->core.tid), mplp_string);
                  // fputs(sam_hdr_tid2name(h, p->b->core.tid), pileup_fp);
                } else {
                  kputc('*', mplp_string);
                  // putc('*', pileup_fp);
                }
                break;
              case MPLP_PRINT_POS:
                ksprintf(mplp_string, "%" PRId64, (int64_t)p->b->core.pos + 1);
                // fprintf(pileup_fp, "%"PRId64, (int64_t) p->b->core.pos + 1);
                break;
              case MPLP_PRINT_MAPQ:
                ksprintf(mplp_string, "%d", p->b->core.qual);
                // fprintf(pileup_fp, "%d", p->b->core.qual);
                break;
              case MPLP_PRINT_RNEXT:
                if (p->b->core.mtid >= 0) {
                  kputs(sam_hdr_tid2name(h, p->b->core.mtid), mplp_string);
                  // fputs(sam_hdr_tid2name(h, p->b->core.mtid), pileup_fp);
                } else {
                  kputc('*', mplp_string);
                  // putc('*', pileup_fp);
                }
                break;
              case MPLP_PRINT_PNEXT:
                ksprintf(mplp_string, "%" PRId64, (int64_t)p->b->core.mpos + 1);
                // fprintf(pileup_fp, "%"PRId64, (int64_t) p->b->core.mpos + 1);
                break;
              case MPLP_PRINT_RLEN:
                ksprintf(mplp_string, "%d", p->b->core.l_qseq);
                // fprintf(pileup_fp, "%d", p->b->core.l_qseq);
                break;
              }
            }
            if (!n) {
              kputc('*', mplp_string);
              // putc('*', pileup_fp);
            }
          }
          flag_value <<= 1;
        }

        /* Print selected tags */
        klist_t(auxlist) *auxlist_p = ((klist_t(auxlist) *)conf->auxlist);
        if (auxlist_p && auxlist_p->size) {
          kliter_t(auxlist) * aux;
          for (aux = kl_begin(auxlist_p); aux != kl_end(auxlist_p);
               aux = kl_next(aux)) {
            n = 0;
            kputc('\t', mplp_string);
            // putc('\t', pileup_fp);
            for (j = 0; j < n_plp[i]; ++j) {
              const bam_pileup1_t *p = &plp[i][j];
              int c =
                  p->qpos < p->b->core.l_qseq ? bam_get_qual(p->b)[p->qpos] : 0;
              if (c < conf->min_baseQ)
                continue;

              if (n > 0) {
                kputc(conf->sep, mplp_string);
                // putc(conf->sep, pileup_fp);
              }
              n++;
              uint8_t *tag_u = bam_aux_get(p->b, kl_val(aux));
              if (!tag_u) {
                kputc(conf->empty, mplp_string);
                // putc(conf->empty , pileup_fp);
                continue;
              }

              int tag_supported = 0;

              /* Tag value is string */
              if (*tag_u == 'Z' || *tag_u == 'H') {
                char *tag_s = bam_aux2Z(tag_u);
                if (!tag_s)
                  continue;
                kputs(tag_s, mplp_string);
                // fputs(tag_s, pileup_fp);
                tag_supported = 1;
              }

              /* Tag value is integer */
              if (*tag_u == 'I' || *tag_u == 'i' || *tag_u == 'C' ||
                  *tag_u == 'c' || *tag_u == 'S' || *tag_u == 's') {
                int64_t tag_i = bam_aux2i(tag_u);
                ksprintf(mplp_string, "%" PRId64 "", tag_i);
                // fprintf(pileup_fp, "%" PRId64 "", tag_i);
                tag_supported = 1;
              }

              /* Tag value is float */
              if (*tag_u == 'd' || *tag_u == 'f') {
                double tag_f = bam_aux2f(tag_u);
                ksprintf(mplp_string, "%lf", tag_f);
                // fprintf(pileup_fp, "%lf", tag_f);
                tag_supported = 1;
              }

              /* Tag value is character */
              if (*tag_u == 'A') {
                char tag_c = bam_aux2A(tag_u);
                kputc(tag_c, mplp_string);
                // putc(tag_c, pileup_fp);
                tag_supported = 1;
              }

              if (!tag_supported) {
                kputc('*', mplp_string);
                // putc('*', pileup_fp);
              }
            }
            if (!n) {
              kputc('*', mplp_string);
              // putc('*', pileup_fp);
            }
          }
        }
      }
    }
    kputc('\n', mplp_string);
    // putc('\n', pileup_fp);
  }

  if (ret < 0) {
    print_error("our_mpileup", "error reading from input file");
    ret = EXIT_FAILURE;
    goto fail;
  }

  if (conf->all) {
    // Handle terminating region
    if (last_tid < 0 && conf->reg && conf->all > 1) {
      last_tid = tid0;
      last_pos = beg0 - 1;
      our_mplp_get_ref(data[0], tid0, &ref, &ref_len);
    } else if (last_tid < 0 && !one_seq && conf->all > 1) {
      last_tid = 0; // --aa on a blank file
    }
    while (last_tid >= 0 && last_tid < sam_hdr_nref(h)) {
      our_mplp_get_ref(data[0], last_tid, &ref, &ref_len);
      while (++last_pos < sam_hdr_tid2len(h, last_tid)) {
        if (last_pos >= end0)
          break;
        if (conf->bed && bed_overlap(conf->bed, sam_hdr_tid2name(h, last_tid),
                                     last_pos, last_pos + 1) == 0)
          continue;
        // TODO handle empty case
        our_print_empty_pileup(mplp_string, conf, sam_hdr_tid2name(h, last_tid),
                               last_pos, n, ref, ref_len);
      }
      last_tid++;
      last_pos = -1;
      if (conf->all < 2 || conf->reg)
        break;
    }
  }

fail:
  // clean up
  // if (pileup_fp && conf->output_fname) fclose(pileup_fp);
  free(buf.s);
  for (i = 0; i < gplp.n; ++i)
    free(gplp.plp[i]);
  free(gplp.plp);
  free(gplp.n_plp);
  free(gplp.m_plp);
  bam_mplp_destroy(iter);
  sam_hdr_destroy(h);
  for (i = 0; i < n; ++i) {
    sam_close(data[i]->fp);
    if (data[i]->iter)
      hts_itr_destroy(data[i]->iter);
    free(data[i]);
  }
  free(data);
  free(plp);
  free(n_plp);
  // free(mp_ref.ref[0]);
  // free(mp_ref.ref[1]);
  return ret;
}

// takes as input:
//   a pointer to a header type (sam_hdr_t)
//   an array of (bam1_t *)
//   an int for length of array
//   a pointer to a reference string
//   an int for ref string length
//
//   a pointer to an empty kstring_t for output (this will be populated with the
//   mpileup string)
//
// destroys the header and nothing else
void bam_and_ref_to_mplp(sam_hdr_t *header, bam1_t **bam_lines, int nbams,
                         char *ref_string, int lref, kstring_t *mplp_string) {

  mplp_conf_t mplp;
  memset(&mplp, 0, sizeof(mplp_conf_t));
  mplp.min_baseQ = 13;
  mplp.capQ_thres = 0;
  mplp.max_depth = MPLP_MAX_DEPTH;
  mplp.flag = MPLP_NO_ORPHAN | MPLP_REALN | MPLP_SMART_OVERLAPS;
  mplp.argc = 4;
  mplp.argv = NULL;
  mplp.rflag_filter = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
  mplp.output_fname = NULL;
  mplp.all = 0;
  mplp.rev_del = 0;
  mplp.sep = ',';
  mplp.empty = '*';
  sam_global_args_init(&mplp.ga);

  bam_index = 0;
  bams = nbams;
  bam_records = bam_lines;

  global_ref_string = ref_string;
  global_ref_length = lref;

  our_mpileup(&mplp, header, mplp_string);
}
