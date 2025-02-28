#ifndef PLACEMENT_HPP
#define PLACEMENT_HPP

#include "coordinates.hpp"
#include "gap_map.hpp"
#include "seed_annotated_tree.hpp"
#include "seeding.hpp"
#include <set>
#include <tbb/concurrent_vector.h>

namespace placement {

using namespace seed_annotated_tree;

// Forward declarations
struct Jaccard;
struct Cosine;
struct ScoreGlobal;
struct CompareByFirst;
// struct LinkedNode;
struct PlacementResult;

// Helper functions for placement
void placementTraversal(
    coordinates::CoordinateTraverser &traverser,
    std::vector<SeedChange> &seedChanges,
    std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> &gapRunBacktracks,
    int64_t hitsInThisGenome, int64_t &maxHitsInAnyGenome, Node *&maxHitsNode,
    int64_t jaccardNumeratorThisGenome, int64_t jaccardDenominatorThisGenome,
    double &bestJaccardScore, Node *&bestJaccardNode,
    double cosineNumeratorThisGenome, double cosineSumOfSquaresThisGenome,
    double &bestCosineScore, Node *&bestCosineNode, Node *parent, Node *current,
    PlacementGlobalState &state, PlacementNodeState &nodeState,
    PlacementResult &result, ::capnp::List<SeedMutations>::Reader &seedIndex,
    ::capnp::List<GapMutations>::Reader &gapIndex, int seedK, int seedS,
    int seedT, bool open, int seedL, Tree *T,
    std::unordered_map<size_t, std::pair<size_t, size_t>> &readSeedCounts,
    std::unordered_map<size_t, int64_t> &currentGenomeSeedCounts,
    const size_t &totalReadSeedCount, const size_t &numGenomes,
    blockExists_t &oldBlockExists, coordinates::blockStrand_t &oldBlockStrand,
    const std::string &true_node_id, const std::string &species,
    const int &read_count, const int &mutation_count, bool firstCall,
    bool &stopReached, Node *startNode, Node *stopNode,
    const std::unordered_set<Node *> &groupNodes, int64_t &searchCount);

void backtrackNode(
    PlacementGlobalState &state, coordinates::CoordinateTraverser &traverser,
    Tree *T, Node *current, PlacementNodeState &nodeState,
    std::unordered_map<size_t, int64_t> &currentGenomeSeedCounts,
    std::vector<std::tuple<int64_t, size_t, bool, int64_t, int64_t>>
        &seedChanges,
    std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>
        &gapRunBacktracks);

void performRecursiveDFS(
    coordinates::CoordinateTraverser &traverser, PlacementGlobalState &state,
    int64_t hitsInThisGenome, int64_t &maxHitsInAnyGenome, Node *&maxHitsNode,
    int64_t jaccardNumeratorThisGenome, int64_t jaccardDenominatorThisGenome,
    double &bestJaccardScore, Node *&bestJaccardNode,
    double cosineNumeratorThisGenome, double cosineSumOfSquaresThisGenome,
    double &bestCosineScore, Node *&bestCosineNode, bool firstCall,
    bool &stopReached, Node *startNode, Node *current, Node *stopNode,
    const std::unordered_set<Node *> &groupNodes, PlacementResult &result,
    ::capnp::List<SeedMutations>::Reader &perNodeSeedMutations_Index,
    ::capnp::List<GapMutations>::Reader &perNodeGapMutations_Index, int seedK,
    int seedS, int seedT, bool open, int seedL, Tree *T,
    std::unordered_map<size_t, std::pair<size_t, size_t>> &readSeedCounts,
    std::unordered_map<size_t, int64_t> &currentGenomeSeedCounts,
    const size_t &totalReadSeedCount, const size_t &numGenomes,
    int64_t &searchCount, coordinates::blockExists_t &rootBlockExists,
    coordinates::blockStrand_t &rootBlockStrand,
    const std::string &true_node_id, const std::string &species,
    const int &read_count, const int &mutation_count);

/**
 * @brief Main placement function that places reads in the phylogenetic tree
 * @param maxHitsInAnyGenome Maximum hits in any genome
 * @param maxHitsNode Node with maximum hits
 * @param bestJaccardScore Best Jaccard score
 * @param bestJaccardNode Node with best Jaccard score
 * @param bestCosineScore Best Cosine score
 * @param bestCosineNode Node with best Cosine score
 * @param result Result structure to store placement scores
 * @param T Tree to place reads in
 * @param index Index reader containing seed and gap mutations
 * @param reads1Path Path to first read file
 * @param reads2Path Path to second read file
 * @param mutMat Mutation matrices
 * @param prefix Output file prefix
 * @param refFileName Reference file name
 * @param samFileName SAM file name
 * @param bamFileName BAM file name
 * @param mpileupFileName MPileup file name
 * @param vcfFileName VCF file name
 * @param aligner Aligner to use
 * @param refNode Reference node ID
 * @param save_jaccard Whether to save Jaccard scores
 * @param show_time Whether to show timing information
 * @param score_proportion Score proportion for filtering
 * @param max_tied_nodes Maximum number of tied nodes
 * @param true_node_id True node ID for evaluation
 * @param species Species information
 * @param read_count Read count
 * @param mutation_count Mutation count
 */
void place(int64_t &maxHitsInAnyGenome, Node *&maxHitsNode,
           double &bestJaccardScore, Node *&bestJaccardNode,
           double &bestCosineScore, Node *&bestCosineNode,
           PlacementResult &result, Tree *T, Index::Reader &index,
           const std::string &reads1Path, const std::string &reads2Path,
           seed_annotated_tree::mutationMatrices &mutMat, std::string prefix,
           std::string refFileName, std::string samFileName,
           std::string bamFileName, std::string mpileupFileName,
           std::string vcfFileName, std::string aligner,
           const std::string &refNode, const bool &save_jaccard,
           const bool &show_time, const float &score_proportion,
           const int &max_tied_nodes, const std::string &true_node_id,
           const std::string &species, const int &read_count,
           const int &mutation_count, const std::string &placementFileName);

/**
 * @brief Batch processing function that places multiple sets of reads from a
 * batch file
 * @param T Tree to place reads in
 * @param index Index reader containing seed and gap mutations
 * @param batchFilePath Path to batch file containing read file paths
 * @param mutMat Mutation matrices
 * @param prefixBase Base prefix for output files
 * @param refFileNameBase Base reference filename (empty to disable)
 * @param samFileNameBase Base SAM filename (empty to disable)
 * @param bamFileNameBase Base BAM filename (empty to disable)
 * @param mpileupFileNameBase Base MPileup filename (empty to disable)
 * @param vcfFileNameBase Base VCF filename (empty to disable)
 * @param aligner Aligner to use
 * @param refNode Reference node ID
 * @param save_jaccard Whether to save Jaccard scores
 * @param show_time Whether to show timing information
 * @param score_proportion Score proportion for filtering
 * @param max_tied_nodes Maximum number of tied nodes
 */
void placeBatch(Tree *T, Index::Reader &index, const std::string &batchFilePath,
                seed_annotated_tree::mutationMatrices &mutMat,
                std::string prefixBase, std::string refFileNameBase,
                std::string samFileNameBase, std::string bamFileNameBase,
                std::string mpileupFileNameBase, std::string vcfFileNameBase,
                std::string aligner, const std::string &refNode,
                const bool &save_jaccard, const bool &show_time,
                const float &score_proportion, const int &max_tied_nodes);

// Scoring structures
struct Jaccard {
  int64_t numerator;
  int64_t denominator;
  Node *node;

  Jaccard(int64_t n, int64_t d, Node *ln)
      : numerator(n), denominator(d), node(ln) {}
};

struct Cosine {
  double numerator;
  double sumOfSquaresGenome;
  Node *node;

  Cosine(double n, double s, Node *ln)
      : numerator(n), sumOfSquaresGenome(s), node(ln) {}
};

struct ScoreGlobal {
  double score;
  Node *node;

  ScoreGlobal(double s, Node *n) : score(s), node(n) {}
};

struct CompareByFirst {
  bool operator()(const std::pair<Node *, double> &a,
                  const std::pair<Node *, double> &b) const {
    return a.second > b.second;
  }
};

// struct LinkedNode {
//     Node* node;
//     LinkedNode* parent;
//     std::vector<LinkedNode*> children;
//     LinkedNode* next;
//     LinkedNode* prev;
//     int64_t hits;
//     std::vector<std::optional<seeding::onSeedsHash>> seeds;

//     LinkedNode(Node* n, int64_t h, const
//     std::vector<std::optional<seeding::onSeedsHash>>& s, LinkedNode* p =
//     nullptr)
//         : node(n), parent(p), next(nullptr), prev(nullptr), hits(h), seeds(s)
//         {}

//     LinkedNode(Node* n, LinkedNode* p = nullptr)
//         : node(n), parent(p), next(nullptr), prev(nullptr), hits(0) {}
// };

struct PlacementResult {
  std::vector<std::pair<Node *, float>> placementScoresJaccard;
  std::vector<std::pair<Node *, float>> placementScoresCosine;
};

} // namespace placement

#endif // PLACEMENT_HPP