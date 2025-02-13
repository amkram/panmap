#include "util.hpp"
#include "seed_annotated_tree.hpp"

extern "C" {
    #include <htslib/sam.h>
}

void getAnchors(std::vector<std::tuple<int64_t, int32_t, int>> &anchors, const std::vector<std::vector<seeding::seed>> &readSeeds,
                const std::vector<std::string> &readSequences, const std::unordered_map<size_t, std::pair<std::vector<uint32_t>, std::vector<uint32_t>>> &seedToRefPositions, int k);

void createSam(
    std::vector<std::vector<seeding::seed>> &readSeeds,
    std::vector<std::string> &readSequences,
    std::vector<std::string> &readQuals,
    std::vector<std::string> &readNames,
    std::string &bestMatchSequence,
    std::unordered_map<size_t, std::pair<std::vector<uint32_t>, std::vector<uint32_t>>> &seedToRefPositions,
    std::string &samFileName,
    int k,
    bool pairedEndReads,
    
    std::vector<char *> &samAlignments,
    std::string &samHeader
);



void createBam(
    std::vector<char *> &samAlignments,
    std::string &samHeader,
    std::string &bamFileName,

    sam_hdr_t * &header,
    bam1_t ** &bamRecords
);


void createMplp(
    std::string &bestMatchSequence,
    sam_hdr_t *header,
    bam1_t **bamRecords,
    int numBams,
    std::string &mpileupFileName,

    char * &mplpString
);

void createMplpBcf(
    std::string &prefix,
    std::string &refFileName,
    std::string &bestMatchSequence,
    std::string &bamFileName,
    std::string &mpileupFileName
);

void createVcf(
    char *mplpString,
    const seed_annotated_tree::mutationMatrices& mutMat,
    std::string &vcfFileName,
    bool keep_alts
);

void createVcfWithMutationMatrices(
  std::string &prefix,
  std::string &mpileupFileName,
  const seed_annotated_tree::mutationMatrices& mutMat,
  std::string &vcfFileName,
  double mutationRate
);
