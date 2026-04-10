#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

extern "C" {
#include <htslib/sam.h>
}

namespace genotyping {
struct mutationMatrices;
}

void createMplp(std::string& bestMatchSequence,
                sam_hdr_t* header,
                bam1_t** bamRecords,
                int numBams,
                std::string& mpileupFileName,
                char*& mplpString);

void createMplpBcf(const std::string& prefix,
                   const std::string& refFileName,
                   const std::string& bestMatchSequence,
                   const std::string& bamFileName,
                   std::string& mpileupFileName,
                   bool baq = false);

void createVcf(char* mplpString, const genotyping::mutationMatrices& mutMat, std::string& vcfFileName, bool keep_alts);

void createVcfWithMutationMatrices(std::string& prefix,
                                   std::string& mpileupFileName,
                                   std::string& vcfFileName,
                                   const std::vector<std::vector<double>>& substMatrixPhred);

int createConsensus(const std::string& vcfFileName,
                    const std::string& refFileName,
                    const std::string& consensusFileName);

// Direct alignment-to-BAM pipeline: parallel minimap2 alignment with direct
// bam1_t construction (no SAM text intermediate). Writes sorted BAM file.
void alignAndWriteBam(std::vector<std::string>& readSequences,
                      std::vector<std::string>& readQuals,
                      std::vector<std::string>& readNames,
                      std::string& reference,
                      const std::string& bamFileName,
                      bool pairedEndReads,
                      int n_threads);
