#include <minimap2/mmpriv.h>
#include <stdbool.h>

void align_reads(const char *reference, int n_reads, const char **reads,
                 const char **quality, const char **read_names, int *r_lens,
                 int *seed_counts, uint8_t **reversed, int **ref_positions,
                 int **qry_positions, char **sam_alignments, int syncmer_k,
                 bool pairedEndReads);