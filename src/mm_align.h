#include <minimap2/mmpriv.h>
#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void align_reads(const char *reference, int n_reads, const char **reads,
                 const char **quality, const char **read_names, int *r_lens,
                 int *seed_counts, uint8_t **reversed, int **ref_positions,
                 int **qry_positions, char **sam_alignments, int syncmer_k,
                 bool pairedEndReads);

// Simple alignment scoring for refinement (no pre-computed seeds needed)
// Returns negative total edit distance (higher = fewer errors = better)
// When paired_end is true, reads are interleaved (R1_0, R2_0, R1_1, R2_1, ...)
// and n_reads must be even. Pairs are mapped together via mm_map_frag.
int64_t score_reads_vs_reference(const char *reference, int n_reads, 
                                  const char **reads, const int *r_lens,
                                  int kmer_size, bool paired_end);

#ifdef __cplusplus
}
#endif