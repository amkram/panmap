#ifndef __PMI_HPP
#define __PMI_HPP

//#include <__config>
#pragma once
#include "seeding.hpp"
#include "seed_annotated_tree.hpp"
#include "index.capnp.h"
#include <unordered_map>
#include <vector>
#include <stack>
#include <tuple>
#include <functional> 


using namespace panmanUtils;
using namespace seeding;
using namespace seed_annotated_tree;

class SeedTracker {
    struct SeedState { int depth = -1; };

    struct SeedKey {
        size_t hash, pos;
        bool operator==(const SeedKey& other) const { return hash == other.hash && pos == other.pos; }
    };

    struct SeedKeyHash {
        size_t operator()(const SeedKey& key) const {
            return std::hash<size_t>()(key.hash) ^ (std::hash<size_t>()(key.pos) << 1);
        }
    };

    std::unordered_map<SeedKey, std::stack<SeedState>, SeedKeyHash> activeSeeds; // Tracks active seeds
    int currentDepth = 0;  // Current DFS depth
    int maxDepth = 0;      // Tracks maximum depth in the tree

public:
    std::unordered_map<size_t, long long> hashCounts;  // Total occurrences of each seed hash

    // Process changes when entering a node
    void enterNode(int depth, const std::vector<std::tuple<size_t, size_t, bool>>& changes) {
        currentDepth = depth;
        maxDepth = std::max(maxDepth, depth); // Update maximum depth

        for (const auto& [hash, pos, isAdded] : changes) {
            SeedKey key = {hash, pos};
            if (isAdded) {
                if (activeSeeds[key].empty() || activeSeeds[key].top().depth == -1) {
                    activeSeeds[key].push({depth});
                }
            } else {
                if (!activeSeeds[key].empty() && activeSeeds[key].top().depth != -1) {
                    int contribution = depth - activeSeeds[key].top().depth;
                    
                    // Prevent underflow when removing
                    hashCounts[key.hash] = std::max(0LL, hashCounts[key.hash] - contribution);

                    activeSeeds[key].pop();
                }
            }
        }
    }

    // Roll back changes when exiting a node
    void exitNode(int depth, const std::vector<std::tuple<size_t, size_t, bool>>& changes) {
        currentDepth = depth;

        for (const auto& [hash, pos, isAdded] : changes) {
            SeedKey key = {hash, pos};
            if (isAdded) {
                if (!activeSeeds[key].empty() && activeSeeds[key].top().depth == depth) {
                    activeSeeds[key].pop();
                }
            } else {
                if (activeSeeds[key].empty() || activeSeeds[key].top().depth == -1) {
                    activeSeeds[key].push({depth});
                }
            }
        }
    }

    // Finalize counts after DFS traversal
    void finalize() {
        for (const auto& [key, stack] : activeSeeds) {
            if (!stack.empty() && stack.top().depth != -1) {
                int contribution = maxDepth - stack.top().depth + 1;
                hashCounts[key.hash] += std::max(contribution, 0);  // Ensure non-negative count
            }
        }
    }

    // Get total occurrences of a specific seed hash
    long long count(size_t seed_hash) const {
        auto it = hashCounts.find(seed_hash);
        return (it != hashCounts.end()) ? it->second : 0;
    }
};


// Forward declarations
struct NodeMutationData;
struct PlacementObjects;

enum posWidth {pos16, pos32, pos64};

enum Metric {
  jaccard,
  weighted_jaccard,
  weighted_count
};


// Scoring weights for alignment metrics
const double MAPQ_WEIGHT = 0.3;      // Weight for mapping quality
const double COVERAGE_WEIGHT = 0.3;   // Weight for genome coverage
const double IDENTITY_WEIGHT = 0.2;    // Weight for sequence identity
const double ALIGN_RATE_WEIGHT = 0.2;  // Weight for read alignment rate
const double PENALTY_WEIGHT = 0.1;     // Weight for penalties

// Weighted Jaccard calculations
struct WeightedJaccardMetrics {
    double raw_score;           // Raw weighted Jaccard score
    double relative_score;      // Normalized relative to read set size
    double asymmetric_score;    // Asymmetric score focusing on read contribution
    double z_score;             // Z-score normalized across all nodes

    void computeZScore(double mean, double stddev) {
        z_score = stddev > 0 ? (raw_score - mean) / stddev : 0.0;
    }
};

// Custom comparator for sorting pairs by first element (string)
struct CompareByFirst {
    bool operator()(const std::pair<Node*, double>& a, const std::pair<Node*, double>& b) const {
        return a.first->identifier < b.first->identifier;  // Compare based on string (node identifier)
    }
};

struct AlignmentMetrics {
    double avgMapQ;
    double coverage;
    double identity;
    double alignmentRate;
    double penalty;
    double finalScore;

    void computeScore() {
        finalScore = MAPQ_WEIGHT * avgMapQ +
               COVERAGE_WEIGHT * coverage +
               IDENTITY_WEIGHT * identity +
               ALIGN_RATE_WEIGHT * alignmentRate -
               PENALTY_WEIGHT * penalty;
    }
};

struct PlacementResult {
  std::string simulation_id;
  std::string true_node_id;
  std::string best_node_id_jaccard;
  std::string best_node_id_cosine;
  double best_node_score_jaccard;
  double best_node_score_cosine;
  std::vector<double> all_jaccard_scores;
  std::vector<double> all_cosine_scores;
  std::vector<std::pair<std::string, double>> top_nodes_jaccard;
  std::vector<std::pair<std::string, double>> top_nodes_cosine;
  std::vector<AlignmentMetrics> top_nodes_metrics_jaccard;  // Store metrics for each top node
  std::vector<AlignmentMetrics> top_nodes_metrics_cosine;  // Store metrics for each top node
  std::vector<WeightedJaccardMetrics> jaccard_metrics; // Store Jaccard metrics for each node
  std::vector<WeightedJaccardMetrics> cosine_metrics; // Store Jaccard metrics for each node
  int32_t k;
  int32_t s;
  Metric metric;
  int32_t read_count;
  int32_t mutation_count;
  double elapsed_time_s;
};

namespace pmi { // functions and types for seed indexing

    void build(Tree *T, Index::Builder &index);
    float align(AlignmentMetrics &metrics, Node *node, Tree *T, int32_t k, int32_t s, int32_t t, bool open, int32_t l, 
                ::capnp::List<SeedMutations>::Reader &perNodeSeedMutations_Reader, 
                ::capnp::List<GapMutations>::Reader &perNodeGapMutations_Reader, 
                const std::string &reads1Path,
                const std::string &reads2Path,
                std::vector<std::vector<seed>> &readSeeds,
                std::vector<std::string> &readSequences,
                std::vector<std::string> &readQuals,
                std::vector<std::string> &readNames, 
                std::string &samFileName, std::string &bamFileName, bool pairedEndReads, std::string &refFileName, std::string aligner="minimap2");

    void genotype(std::string prefix, std::string refFileName, std::string bestMatchSequence, std::string bamFileName, std::string mpileupFileName, std::string vcfFileName, seed_annotated_tree::mutationMatrices &mutMat);
    void place(PlacementResult &result,
      Tree *T, Index::Reader &index, const std::string &reads1Path, const std::string &reads2Path,
      seed_annotated_tree::mutationMatrices &mutMat, std::string prefix, std::string refFileName, std::string samFileName,
      std::string bamFileName, std::string mpileupFileName, std::string vcfFileName, std::string aligner,
      const std::string& refNode, const bool& save_jaccard, const bool& show_time, std::ofstream& logFile,
      const float& score_proportion = 0.01, const int& max_tied_nodes = 16, const std::string& true_node_id="", const std::string& species="", const int& read_count=0, const int& mutation_count=0);

    void parallel_tester(Tree *T, Index::Reader &index, const std::string &reads1Path, const std::string &reads2Path, const std::string &prefix);

    void place_per_read(
      Tree *T, Index::Reader &index, const std::string &reads1Path, const std::string &reads2Path,
      const int& maximumGap, const int& minimumCount, const int& minimumScore, const double& errorRate,
      const int& redoReadThreshold, const bool& recalculateScore, const bool& rescueDuplicates,
      const double& rescueDuplicatesThreshold, const double& excludeDuplicatesThreshold,
      const std::string& preEMFilterMethod, const int& preEMFilterNOrder, const int& preEMFilterMBCNum, const int& emFilterRound, const int& checkFrequency,
      const int& removeIteration, const double& insigProb, const int& roundsRemove, const double& removeThreshold,
      const bool& leafNodesOnly, const bool& callSubconsensus, const std::string& prefix, const bool& save_kminmer_binary_coverage);

    void evaluate(Tree *T, std::string input_tsv, Index::Reader &index, 
              seed_annotated_tree::mutationMatrices &mutMat, std::string aligner, std::string species,
              std::chrono::high_resolution_clock::time_point start_time, std::string default_index_path,
              std::string default_mutmat_path,
              const float& score_proportion = 0.01, const int& max_tied_nodes = 16);

} // namespace pmi

/* Expose some pmi.cpp helpers for unit testing also tree.cpp uses
 * applyMutations for now */
using namespace pmi;

// void buildHelper(mutableTreeData &data, seedMap_t &seedMap, int32_t &k, int32_t &s, ::capnp::List<Mutations>::Builder &mutations,
//                  Tree *T, const Node *node, const globalCoords_t &globalCoords,
//                  CoordNavigator &navigator, int64_t &dfsIndex, posWidth &width, std::map<int64_t, int64_t> &gapRuns);
// void applyMutations(mutableTreeData &data, seedMap_t &seedMap,
//                     blockMutationInfo_t &blockMutData,
//                     std::vector<tupleRange> &recompRanges,
//                     mutationInfo_t &nucMutData, Tree *T, const Node *node,
//                     globalCoords_t &globalCoords, const ::capnp::List<Mutations>::Builder &mutations);
// void undoMutations(mutableTreeData &data, ::capnp::List<Mutations>::Builder &mutations, Tree *T,
//                     Node *node, const blockMutationInfo_t &blockMutData,
//                    const mutationInfo_t &nucMutData);
void flipCoords(int32_t blockId, globalCoords_t &globalCoords);
void updateGapMapStep(std::map<int64_t, int64_t>& gapMap, const std::pair<bool, std::pair<int64_t, int64_t>>& update, std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& backtrack, std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& gapMapUpdates, bool recordGapMapUpdates=true);
void updateGapMap(std::map<int64_t, int64_t>& gapMap, const std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& updates, std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& backtrack, std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& gapMapUpdates);
std::vector<std::pair<int64_t, int64_t>> invertRanges(const std::vector<std::pair<int64_t, int64_t>>& nucRanges, const std::pair<int64_t, int64_t>& invertRange);
void invertGapMap(std::map<int64_t, int64_t>& gapMap, const std::pair<int64_t, int64_t>& invertRange, std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& backtrack, std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& gapMapUpdates);
void makeCoordIndex(std::map<int64_t, int64_t>& coordIndex, const std::map<int64_t, int64_t>& gapMap, const std::vector<std::pair<int64_t, int64_t>>& blockRanges);

// // Go upstream until neededNongap nucleotides are seen and return the new coord.
// tupleCoord_t expandLeft(const CoordNavigator &navigator, tupleCoord_t coord,
//                         int neededNongap, blockExists_t &blockExists);

// // Go downstream until neededNongap nucleotides are seen and return the new coord.
// tupleCoord_t expandRight(const CoordNavigator &navigator, tupleCoord_t coord,
//                          int neededNongap, blockExists_t &blockExists, blockStrand_t &blockStrand);

// // Merges each range with overlapping ranges after expanding left and right
// // by `neededNongap` non-gap nucleotides.
// std::vector<tupleRange> expandAndMergeRanges(const CoordNavigator &navigator, std::vector<tupleRange> &ranges, int neededNongap, blockExists_t &blockExists);
int64_t tupleToScalarCoord(const tupleCoord_t &coord, const globalCoords_t &globalCoords);

#endif