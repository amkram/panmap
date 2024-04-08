#ifndef ALIGN_HPP
#define ALIGN_HPP

#include "seeding.hpp"
#include <vector>
#include "tree.hpp"

using namespace seeding;

namespace align {

    void mapToTarget(PangenomeMAT::Tree *T, std::unordered_map<std::string, std::vector<int32_t>> &seedToRefPositions, std::vector<std::vector<seed>> &readSeeds, std::vector<std::string> &readSequences, std::vector<std::string> &readQuals, std::vector<std::string> &readNames, std::string &bestMatchSequence, int k, std::string &samFileName, std::string &bamFileName, std::string &mpileupFileName, std::string &vcfFileName);

    void align_reads(const char *reference, int n_reads, char **reads, char **quality, char **read_names, int *r_lens, int *seed_counts, uint8_t **reversed, int **ref_positions, int **qry_positions, char** sam_alignments, int syncmer_k);
}

#endif