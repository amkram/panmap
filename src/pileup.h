#include "samtools_compat.h"

void bam_and_ref_to_mplp(sam_hdr_t *header, bam1_t **bam_lines, int nbams,
                         char *ref_string, int lref, kstring_t *mplp_string);