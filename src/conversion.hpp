#include "seeding.hpp"
#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

extern "C" {
#include <htslib/sam.h>
}

// Forward declarations
namespace genotyping {
    struct mutationMatrices;
}

// Main functions
void createSam(
    std::vector<std::vector<seeding::seed_t>> &readSeeds,
    std::vector<std::string> &readSequences,
    std::vector<std::string> &readQuals,
    std::vector<std::string> &readNames,
    std::string &bestMatchSequence,
    std::unordered_map<size_t, std::pair<std::vector<uint32_t>, std::vector<uint32_t>>> &seedToRefPositions,
    std::string &samFileName,
    int k,
    bool shortenSyncmers,
    bool pairedEndReads,
    std::vector<char *> &samAlignments,
    std::string &samHeader);

void createBam(std::vector<char *> &samAlignments,
               std::string &samHeader,
               std::string &bamFileName,
               sam_hdr_t *&header,
               bam1_t **&bamRecords);

void createMplp(std::string &bestMatchSequence,
                sam_hdr_t *header,
                bam1_t **bamRecords,
                int numBams,
                std::string &mpileupFileName,
                char *&mplpString);

void createMplpBcf(std::string &prefix,
                   std::string &refFileName,
                   std::string &bestMatchSequence,
                   std::string &bamFileName,
                   std::string &mpileupFileName);

void createVcf(char *mplpString,
               const genotyping::mutationMatrices &mutMat,
               std::string &vcfFileName,
               bool keep_alts);

void createVcfWithMutationMatrices(
    std::string &prefix,
    std::string &mpileupFileName,
    const genotyping::mutationMatrices &mutMat,
    std::string &vcfFileName,
    double mutationRate);
