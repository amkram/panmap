#include <htslib/cram.h>
#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/khash_str2int.h>
#include <htslib/klist.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>
#include <samtools/bam_plcmd.h>
#include <samtools/samtools.h>

void bam_and_ref_to_mplp(sam_hdr_t *header, bam1_t **bam_lines, int nbams,
                         char *ref_string, int lref, kstring_t *mplp_string);