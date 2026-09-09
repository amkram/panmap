#pragma once

typedef struct {
  int mm2_k;
  int mm2_w;
  int mm2_match;
  int mm2_mismatch;
  int mm2_gap_open;
  int mm2_gap_extend;

  int mm2_k_set;
  int mm2_w_set;
  int mm2_match_set;
  int mm2_mismatch_set;
  int mm2_gap_open_set;
  int mm2_gap_extend_set;

  double bwaMaxDiff;
  int bwaMaxGapOpen;
  int bwaSeedLength;
  int bwaTrimQual;
} aligner_options_t;