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

// Returns 0 on success, non-zero if bcftools mpileup fails.
int createMplpBcf(const std::string& prefix,
                  const std::string& refFileName,
                  const std::string& bestMatchSequence,
                  const std::string& bamFileName,
                  std::string& mpileupFileName,
                  bool baq = false,
                  const std::string& refName = "ref");

// Returns 0 on success, non-zero if bcftools call fails.
int createVcfWithMutationMatrices(std::string& prefix,
                                  std::string& mpileupFileName,
                                  std::string& vcfFileName,
                                  const std::vector<std::vector<double>>& substMatrixPhred,
                                  int minDepth,
                                  double minQual);

int createConsensus(const std::string& vcfFileName,
                    const std::string& refFileName,
                    const std::string& consensusFileName,
                    const std::string& consensusHeader = "");

// Parallel minimap2/bwa alignment with direct bam1_t construction (no SAM text
// intermediate). Writes a sorted BAM file. Returns 0 on success, non-zero on a
// BAM write/close error.
int alignAndWriteBam(std::vector<std::string>& readSequences,
                      std::vector<std::string>& readQuals,
                      std::vector<std::string>& readNames,
                      std::string& reference,
                      const std::string& bamFileName,
                      bool pairedEndReads,
                      int n_threads,
                      std::string aligner = "minimap2",
                      const std::string& refName = "ref");
