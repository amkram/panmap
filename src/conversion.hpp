#include "util.hpp"
#include "tree.hpp"


void createSam(
    std::vector<seeding::seed> &refSeeds,
    std::vector<std::vector<seeding::seed>> &readSeeds,
    std::vector<std::string> &readSequences,
    std::vector<std::string> &readQuals,
    std::vector<std::string> &readNames,
    std::string &bestMatchSequence,
    std::unordered_map<std::string, std::vector<int32_t>> &seedToRefPositions,
    std::string &samFileName,
    int k,
    
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

void createVcf(
    char *mplpString,
    const tree::mutationMatrices& mutMat,
    std::string &vcfFileName
);