#pragma once

#include "capnp/list.h"
#include "index.capnp.h"
#include "panman.hpp"
#include "progress_state.hpp"
#include "state.hpp"
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <tbb/concurrent_vector.h>
#include <unordered_map>
#include <utility>
#include <vector>
#include <absl/container/flat_hash_map.h>

// Forward declaration for indexing 
namespace indexing {
void processNodesByLevel(
    const std::vector<panmanUtils::Node*>& nodes,
    std::function<void(panmanUtils::Node*)> processFunction);
}

namespace placement {

// Forward declarations
struct PlacementResult;

// External declarations
extern std::shared_ptr<PlacementProgressState> progress_state;

// Parameters for traversal - consolidated
struct TraversalParams {
  int k = 32;              // k-mer size
  int s = 8;              // syncmer parameter s
  int t = 0;               // t-syncmer parameter
  bool open = false;        // Whether to use open syncmers
  double scoreScale = 1.0; // Scaling factor for scores
  std::string debug_node_id;
};

// Track global state during placement
struct PlacementGlobalState {
  // Reference to seed mutations from the index
  ::capnp::List<::SeedMutations>::Reader perNodeSeedMutations;

  // Reference to gap mutations from the index
  ::capnp::List<::GapMutations>::Reader perNodeGapMutations;

  // Reference to node path information from the index
  ::capnp::List<::NodePathInfo>::Reader nodePathInfo;

  // Reference to block information from the index
  ::capnp::List<::BlockInfo>::Reader blockInfo;

  // Reference to ancestor matrix from the index
  ::capnp::List<::capnp::List<bool>>::Reader ancestorMatrix;

  // K-mer dictionary and seed frequencies
  std::unordered_map<uint32_t, std::string> kmerDictionary;  // Dictionary ID -> k-mer sequence
  std::unordered_map<size_t, uint32_t> kmerHashToId;         // Hash -> Dictionary ID
  std::unordered_map<size_t, bool> hashOrientation;          // Hash -> is reverse orientation
  std::unordered_map<size_t, bool> canonicalHashes;          // Hash -> is canonical form
  
  absl::flat_hash_map<size_t, int64_t> seedFreqInReads;      // Hash -> read frequency count
  absl::flat_hash_map<size_t, std::string> hashToKmer;      // Hash -> k-mer sequence
  
  float totalReadSeedCount = 0.0f;                          // Total seeds in reads (for normalization)
  double jaccardDenominator = 0.0;                          // Denominator for Jaccard calculation
  size_t readUniqueSeedCount = 0;                           // Count of unique seed hashes in reads

  // Flag to stop the traversal
  std::atomic<bool> stopTraversal{false};
  
  // Store k-mer size for hash calculations
  int kmerSize = 32;  // Default value, will be set from params in place()

  // Helper methods for node relationships
  bool isAncestor(const std::string &potentialAncestorId, const std::string &nodeId) const;
  uint32_t getNodeLevel(const std::string &nodeId) const;
  std::vector<int32_t> getNodeActiveBlocks(const std::string &nodeId) const;
};

// Helper class to store score info for a single node during placement
class PlacementNodeScore {
public:
    int64_t hitsInThisGenome;
    int64_t currentJaccardNumerator;
    // currentJaccardDenominator is calculated on the fly based on uniqueSeedHashes.size() and state.jaccardDenominator
    double currentCosineNumerator;
    double currentCosineDenominator;
    absl::flat_hash_map<int64_t, seeding::seed_t> kmerSeedMap;
    absl::flat_hash_map<size_t, int64_t> currentGenomeSeedCounts; // For weighted Jaccard and Cosine

    // New fields for additional scoring metrics
    int64_t rawSeedMatchScore;                   // For raw seed matches based on read frequency
    int64_t jaccardPresenceNumerator;            // Intersection size for presence/absence Jaccard
    absl::flat_hash_set<size_t> currentGenomeUniqueSeedHashes; // Unique seeds in current genome (for presence/absence Jaccard denominator)
};

// Consolidated score tracking and results
struct PlacementResult {
  // Hits-based results (this seems to be the current "weighted Jaccard" numerator)
  int64_t maxHitsInAnyGenome = 0; // This might be better named to reflect its actual calculation if it's not raw hits
  panmanUtils::Node *maxHitsNode = nullptr;
  std::vector<panmanUtils::Node *> tiedMaxHitsNodes;

  // Raw Seed Match Score (New)
  int64_t bestRawSeedMatchScore = 0;
  panmanUtils::Node *bestRawSeedMatchNode = nullptr;
  std::vector<panmanUtils::Node *> tiedRawSeedMatchNodes;

  // Jaccard-based results (current "weighted" Jaccard)
  double bestJaccardScore = 0.0; // This is the weighted Jaccard
  panmanUtils::Node *bestJaccardNode = nullptr;
  std::vector<panmanUtils::Node *> tiedJaccardNodes;

  // Jaccard Index (Presence/Absence) (New)
  double bestJaccardPresenceScore = 0.0;
  panmanUtils::Node *bestJaccardPresenceNode = nullptr;
  std::vector<panmanUtils::Node *> tiedJaccardPresenceNodes;

  // Cosine-based results
  double bestCosineScore = 0.0;
  panmanUtils::Node *bestCosineNode = nullptr;
  std::vector<panmanUtils::Node *> tiedCosineNodes;

  // Weighted result
  double bestWeightedScore = 0.0;
  panmanUtils::Node *bestWeightedNode = nullptr;
  std::vector<panmanUtils::Node *> tiedWeightedNodes;

  // Current tracking values
  int64_t hitsInThisGenome = 0;
  int64_t currentJaccardNumerator = 0;
  int64_t currentJaccardDenominator = 0;
  double currentCosineNumerator = 0.0;
  double currentCosineDenominator = 0.0;
  absl::flat_hash_map<size_t, int64_t> currentGenomeSeedCounts;
  
  // Added for placement - map of node seeds without dependency on StateManager
  // This allows us to track seeds directly in the placement result
  // Format: nodeId -> position -> seed
  std::unordered_map<std::string, std::unordered_map<int64_t, seeding::seed_t>> nodeSeedMap;

  // Performance metrics
  int64_t totalReadsProcessed = 0;
  double totalTimeSeconds = 0.0;
  
  // Helper methods
  void updateHitsScore(panmanUtils::Node* node, int64_t hits); // This is for the existing "maxHitsInAnyGenome"
  void updateRawSeedMatchScore(panmanUtils::Node* node, int64_t score); // New
  void updateJaccardScore(panmanUtils::Node* node, double score); // This is for weighted Jaccard
  void updateJaccardPresenceScore(panmanUtils::Node* node, double score); // New
  void updateCosineScore(panmanUtils::Node* node, double score);
  void updateWeightedScore(panmanUtils::Node* node, double score, double scale);
};

// Core functions for placement
void processNodeMutations(panmanUtils::Node *node,
                          state::StateManager &stateManager,
                          PlacementGlobalState &state,
                          PlacementResult &result,
                          const TraversalParams &params);

void placementTraversal(state::StateManager &stateManager,
                        PlacementResult &result,
                        panmanUtils::Tree *T, 
                        PlacementGlobalState &state,
                        const TraversalParams &params,
                        const std::string &debug_node_id_param);

/**
 * @brief Load seed information from index into the StateManager for placement
 * 
 * This function decodes the quaternary-encoded seed changes from the index
 * and applies them to the StateManager's hierarchical seed stores.
 * 
 * @param stateManager StateManager to populate with seed data
 * @param index The index to read seed data from
 */
void loadSeedsFromIndex(state::StateManager& stateManager, const ::Index::Reader& index);

void place(PlacementResult &result, panmanUtils::Tree *T,
           ::Index::Reader &index, const std::string &reads1,
           const std::string &reads2,
           std::string &outputPath,
           const std::string &indexPath,
           const std::string &debug_node_id_param);

void placeBatch(panmanUtils::Tree *T, ::Index::Reader &index,
                const std::string &batchFilePath,
                std::string prefixBase, 
                std::string refFileNameBase,
                std::string samFileNameBase, 
                std::string bamFileNameBase,
                std::string mpileupFileNameBase, 
                std::string vcfFileNameBase,
                std::string aligner, 
                const std::string &refNode,
                const bool &save_jaccard, 
                const bool &show_time,
                const float &score_proportion, 
                const int &max_tied_nodes,
                const std::string &indexPath,
                const std::string &debug_node_id_param);

// Helper function
std::pair<double, double> getCosineDelta(
    bool isRemoval, bool isAddition, size_t seedHash,
    const absl::flat_hash_map<size_t, int64_t>& readSeedCounts,
    const absl::flat_hash_map<size_t, int64_t>& genomeSeedCounts);

// NEW: Declaration for placement summary dump function
void dumpPlacementSummary(const PlacementResult& result, const std::string& outputFilename);

// NEW: Declaration for kmer debug dump function 
void dumpKmerDebugData(
    const std::vector<std::string>& readSequences,
    int k,
    const std::string& outputFilename);

// Consolidated seed processing functions
void processSeedOperation(
    panmanUtils::Node* node,
    state::StateManager& stateManager,
    PlacementGlobalState& state,
    PlacementResult& result,
    absl::flat_hash_set<size_t>& uniqueSeedHashes,
    int64_t pos,
    bool isRemoval,
    bool isAddition,
    const seeding::seed_t* newSeed = nullptr);

seeding::seed_t createAndProcessSeed(
    panmanUtils::Node* node,
    state::StateManager& stateManager,
    PlacementGlobalState& state,
    PlacementResult& result,
    absl::flat_hash_set<size_t>& uniqueSeedHashes,
    int64_t pos,
    const std::string& kmerStr,
    int64_t endPos,
    int params_k,
    size_t& seedAdditions,
    size_t& dictionaryLookups);

std::tuple<double, double, size_t> calculateAndUpdateScores(
    panmanUtils::Node* node,
    PlacementGlobalState& state,
    PlacementResult& result,
    const TraversalParams& params,
    absl::flat_hash_set<size_t>& uniqueSeedHashes);

// Dumps detailed seed information for debugging placement issues
void dumpPlacementDebugData(
    state::StateManager& stateManager,
    PlacementGlobalState& state,
    const std::vector<std::string>& nodesToDump,
    size_t maxNodes,
    const TraversalParams& params,
    const std::string& outputFilename,
    const std::vector<std::string>& readSequences);

} // namespace placement