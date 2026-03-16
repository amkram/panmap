#include <minimap2/mmpriv.h>
#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// Align reads to reference using minimap2's native seeding.
// r_lens is INPUT read lengths but OVERWRITTEN with 1-based ref start positions
// (or INT_MAX for unmapped reads) on output.
void align_reads(const char *reference, int n_reads, const char **reads,
                 const char **quality, const char **read_names, int *r_lens,
                 char **sam_alignments, bool pairedEndReads);

// Simple alignment scoring for refinement (no pre-computed seeds needed)
// Returns negative total edit distance (higher = fewer errors = better)
// When paired_end is true, reads are interleaved (R1_0, R2_0, R1_1, R2_1, ...)
// and n_reads must be even. Pairs are mapped together via mm_map_frag.
int64_t score_reads_vs_reference(const char *reference, int n_reads,
                                  const char **reads, const int *r_lens,
                                  int kmer_size, bool paired_end);

// Structured alignment result for one read
typedef struct {
    int32_t pos;         // 1-based ref position (INT_MAX if unmapped)
    int32_t rs, re;      // 0-based ref start/end
    int32_t qs, qe;      // query start/end
    uint8_t mapq;
    uint8_t rev;         // reverse strand
    uint8_t proper_frag; // proper pair flag (from minimap2)
    int32_t n_cigar;
    uint32_t *cigar;     // BAM-encoded CIGAR; malloc'd, caller frees
} read_align_t;

// Result for a read pair (or single read)
typedef struct {
    read_align_t r1;
    read_align_t r2;     // only valid when paired
    int mapped;          // 1 if primary hits exist
} align_pair_result_t;

// Parallel alignment returning structured results (no SAM text).
// For paired-end: reads are interleaved (R1_0, R2_0, R1_1, R2_1, ...),
// n_reads must be even, results has n_reads/2 elements.
// For single-end: results has n_reads elements.
// Each result's cigar arrays are malloc'd; caller must free them.
void align_reads_direct(const char *reference, int n_reads,
                        const char **reads, const char **quality,
                        const char **read_names, const int *r_lens,
                        align_pair_result_t *results, bool pairedEndReads,
                        int n_threads);

#ifdef __cplusplus
}
#endif
