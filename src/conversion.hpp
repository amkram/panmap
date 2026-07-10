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

void createMplpBcf(const std::string& prefix,
                   const std::string& refFileName,
                   const std::string& bestMatchSequence,
                   const std::string& bamFileName,
                   std::string& mpileupFileName,
                   bool baq = false,
                   const std::string& refName = "ref");

void createVcfWithMutationMatrices(std::string& prefix,
                                   std::string& mpileupFileName,
                                   std::string& vcfFileName,
                                   const std::vector<std::vector<double>>& substMatrixPhred);

int createConsensus(const std::string& vcfFileName,
                    const std::string& refFileName,
                    const std::string& consensusFileName,
                    const std::string& consensusHeader = "");

// Parallel minimap2 alignment with direct bam1_t construction (no SAM text
// intermediate). Writes sorted BAM file.
void alignAndWriteBam(std::vector<std::string>& readSequences,
                      std::vector<std::string>& readQuals,
                      std::vector<std::string>& readNames,
                      std::string& reference,
                      const std::string& bamFileName,
                      bool pairedEndReads,
                      int n_threads,
                      bool useBwa = false,
                      const std::string& refName = "ref");
