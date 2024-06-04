#ifndef __MGSR_HPP
#define __MGSR_HPP

#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <tbb/concurrent_vector.h>
#include "PangenomeMAT.hpp"
#include "tree.hpp"
#include "pmi.hpp"
#include <eigen3/Eigen/Dense>




typedef size_t hash_t;
typedef std::pair<std::vector<std::tuple<size_t, int32_t, int32_t, bool, int32_t>>, std::unordered_set<size_t>> readSeedmers_t;
typedef std::tuple<int32_t, int32_t, int32_t, int32_t, bool, int32_t> match_t;
namespace mgsr {
    struct seedmers {
        //       beg                hash    rev
        // std::map<int32_t, std::pair<size_t, bool>> positionMap;
        //       beg                 end      hash    rev
        std::map<int32_t, std::tuple<int32_t, size_t, bool>> positionMap;
        //                 hash                       begs
        std::unordered_map<size_t, std::unordered_set<int32_t>> seedmersMap;
    };

    // pmi index file with k = k, s = s, l = 1
    // k
    // l
    // tree
    void buildSeedmer(pmi::seedIndex &Index, PangenomeMAT::Tree *T, const size_t l, const size_t k, const size_t s, std::stringstream& seedmersOutStream);
    void buildSeedmerHelper(tree::mutableTreeData &data, seedMap_t seedMap, pmi::seedIndex &index, mgsr::seedmers seedmersIndex, std::map<int32_t, std::pair<int32_t, size_t>> curSeeds, PangenomeMAT::Tree *T, const PangenomeMAT::Node *node, const int32_t l, const int32_t k, const int32_t s, const tree::globalCoords_t &globalCoords, std::stringstream& seedmersOutStream);
    void scorePseudo(std::ifstream &indexFile, const std::string &reads1Path, const std::string &reads2Path, std::unordered_map<std::string, tbb::concurrent_vector<std::pair<int32_t, double>>>& allScores, std::vector<std::pair<int32_t, std::vector<size_t>>>& numReadDuplicates, std::unordered_map<std::string, std::string>& leastRecentIdenticalAncestor, std::unordered_map<std::string, std::unordered_set<std::string>>& identicalSets, int32_t& numReads, Tree *T, std::vector<std::string>& readSequences, std::vector<std::string>& readQuals, std::vector<std::string>& readNames, const int& maximumGap, const int& minimumCount, const int& minimumScore, const double& errorRate);
    void squaremHelper(PangenomeMAT::Tree *T, const std::unordered_map<std::string, tbb::concurrent_vector<std::pair<int32_t, double>>>& allScores, const std::vector<std::pair<int32_t, std::vector<size_t>>>& numReadDuplicates, const int32_t& numReads, const std::unordered_map<std::string, std::string>& leastRecentIdenticalAncestor, const std::unordered_map<std::string, std::unordered_set<std::string>>& identicalSets, Eigen::MatrixXd& probs, std::vector<std::string>& nodes, Eigen::VectorXd& props, double& llh, const int32_t& roundsRemove, const double& removeThreshold, std::string exclude);
    void accio(const std::unordered_map<std::string, tbb::concurrent_vector<std::pair<int32_t, double>>>& allScores, const std::vector<std::string>& nodes, const Eigen::MatrixXd& probs, const Eigen::VectorXd& props, const std::vector<std::pair<int32_t, std::vector<size_t>>>& numReadDuplicates, std::unordered_map<std::string, std::unordered_set<size_t>>& assignedReads);
    
    // for testing
    std::unordered_map<std::string, std::pair<double, double>> getReadAssignmentAccuracy( const std::unordered_map<std::string, std::unordered_set<size_t>>& assignedReads, const std::vector<std::string>& nodes, const std::vector<std::string>& readNames, const std::unordered_map<std::string, std::string>& leastRecentIdenticalAncestor);
}




#endif
