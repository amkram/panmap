#include "util.hpp"
#include "tree.hpp"
#include <htslib/sam.h>


void createSam(
    std::vector<std::vector<seeding::seed>> &readSeeds,
    std::vector<std::string> &readSequences,
    std::vector<std::string> &readQuals,
    std::vector<std::string> &readNames,
    std::string &bestMatchSequence,
    std::unordered_map<std::string, std::vector<int32_t>> &seedToRefPositions,
    const bool& makeSam,
    const std::string& samFileName,
    int k,
    bool pairedEndReads,
    
    std::vector<char *> &samAlignments,
    std::string &samHeader
);



void createBam(
    std::vector<char *> &samAlignments,
    std::string &samHeader,
    const bool& makeBam,
    const std::string& bamFileName,

    sam_hdr_t * &header,
    bam1_t ** &bamRecords
);


void createMplp(
    std::string &bestMatchSequence,
    sam_hdr_t *header,
    bam1_t **bamRecords,
    int numBams,
    const bool& makeMPileup,
    const std::string& mpileupFileName,

    char * &mplpString
);

void createVcf(
    char *mplpString,
    const tree::mutationMatrices& mutMat,
    const bool& makeVCF,
    const std::string& vcfFileName
);