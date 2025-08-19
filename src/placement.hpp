#pragma once

#include "capnp/list.h"
#include "mgsr_index.capnp.h"
#include "panman.hpp"
#include "progress_state.hpp"
#include "seeding.hpp"
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

// Type alias for seed structures
using seed_t = seeding::seed_t;

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
  int k = 0;              // k-mer size
  int s = 0;              // syncmer parameter s
  int t = 0;               // t-syncmer parameter
  bool open = false;        // Whether to use open syncmers
  double scoreScale = 1.0; // Scaling factor for scores
  std::string debug_node_id;
};

// Track global state during placement (MGSR-only simplified)
struct PlacementGlobalState {
    // Seed frequencies in reads
    absl::flat_hash_map<size_t, int64_t> seedFreqInReads;      // Hash -> read frequency count
    absl::flat_hash_map<size_t, std::string> hashToKmer;       // Optional: hash -> k-mer sequence (if available)
    float totalReadSeedCount = 0.0f;
    double jaccardDenominator = 0.0;
    size_t readUniqueSeedCount = 0;
    std::atomic<bool> stopTraversal{false};
    int kmerSize = 32;  // set from params

    // Active seed set during traversal (map of seedIndex -> count)
    std::unordered_map<uint64_t, uint32_t> activeSeedIndices;
    
    // MGSR index data
    ::capnp::List<SeedInfo>::Reader seedInfo;
    ::capnp::List<NodeChanges>::Reader perNodeChanges;
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
// Note: Legacy processNodeMutations and placementTraversal removed 
// MGSR placement now uses inline traversal in place() function

void place(PlacementResult &result, panmanUtils::Tree *T,
           ::MGSRIndex::Reader &mgsrIndex, const std::string &reads1,
           const std::string &reads2,
           std::vector<std::vector<seeding::seed_t>>& readSeeds,
           std::vector<std::string>& readSequences,
           std::vector<std::string>& readNames,
           std::vector<std::string>& readQuals,
           std::string &outputPath,
           const std::string &indexPath,
           const std::string &debug_node_id_param);
// Batch mode now accepts path to MGSR index file instead of legacy Index::Reader
void placeBatch(panmanUtils::Tree *T,
                ::MGSRIndex::Reader &mgsrIndex,
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

// Consolidated seed processing functions (MGSR compatible)
void processSeedOperation(
    panmanUtils::Node* node,
    state::StateManager& stateManager,
    PlacementGlobalState& state,
    PlacementResult& result,
    absl::flat_hash_set<size_t>& uniqueSeedHashes,
    uint32_t seedIndex);

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