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
// Returns sum of primary alignment scores across all reads
int64_t score_reads_vs_reference(const char *reference, int n_reads, 
                                  const char **reads, const int *r_lens,
                                  int kmer_size);

#ifdef __cplusplus
}
#endif