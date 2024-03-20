#include "util.hpp"


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