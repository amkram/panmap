#pragma once

#include "capnp/list.h"
#include "index.capnp.h"
#include "panman.hpp"
#include "panmap_utils.hpp"
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

// Forward declarations
namespace placement {
    struct PlacementGlobalState;
}

// Efficient genome state representation for traversal
// This represents the accumulated genome state at a node in phylogenetic tree
// The root starts empty (ancestral state), and each node adds phylogenetic mutations
struct GenomeState {
  // Position → seed index mapping (for handling deletions and substitutions)
  absl::flat_hash_map<uint32_t, uint32_t> positionToSeedIndex;
  
  // Core seed count data (this is the main memory consumer)
  absl::flat_hash_map<size_t, int64_t> seedCounts;
  absl::flat_hash_set<size_t> uniqueSeeds;
  
  // Cached similarity metrics (computed on-demand to avoid recalculation)
  mutable bool metricsValid = false;
  mutable int64_t weightedJaccardNumerator = 0;
  mutable int64_t jaccardNumerator = 0;
  mutable int64_t jaccardDenominator = 0;
  mutable double cosineNumerator = 0.0;
  mutable double cosineDenominator = 0.0;
  
  // Efficient copy constructor that preserves cached metrics when possible
  GenomeState() = default;
  GenomeState(const GenomeState& parent) : 
      positionToSeedIndex(parent.positionToSeedIndex),
      seedCounts(parent.seedCounts), 
      uniqueSeeds(parent.uniqueSeeds),
      metricsValid(false) {} // Always recalculate metrics for child nodes
      
  // Apply a node's phylogenetic changes to create child state
  void applyNodeChanges(const std::vector<std::pair<size_t, int64_t>>& seedDeltas);
  
  // Forward declare the method - implementation will be in .cpp with proper namespace
  void computeMetrics(const placement::PlacementGlobalState& state) const;
  
  // Estimate memory usage for monitoring
  size_t estimateMemoryUsage() const {
      return positionToSeedIndex.size() * (sizeof(uint32_t) + sizeof(uint32_t)) +
             seedCounts.size() * (sizeof(size_t) + sizeof(int64_t)) + 
             uniqueSeeds.size() * sizeof(size_t) + 
             sizeof(GenomeState);
  }
};

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
  int l = 0;               // k-minimizer window size (0 = use raw syncmers)
  bool open = false;        // Whether to use open syncmers
  bool useRawSeeds = false; // Whether to use raw syncmers instead of k-min-mers
  double scoreScale = 1.0; // Scaling factor for scores
  std::string debug_node_id;
  bool verify_scores = false; // Whether to recompute scores from scratch for verification
};

// Track global state during placement
struct PlacementGlobalState {
    // Seed frequencies in reads
    absl::flat_hash_map<size_t, int64_t> seedFreqInReads;      // Hash -> read frequency count
    absl::flat_hash_map<size_t, std::string> hashToKmer;       // Optional: hash -> k-mer sequence (if available)
    float totalReadSeedCount = 0.0f;
    int64_t totalReadSeedFrequency = 0;  // Sum of all read seed frequencies (constant denominator for Jaccard)
    double jaccardDenominator = 0.0;
    size_t readUniqueSeedCount = 0;
    double readMagnitude = 0.0;  // Precomputed read magnitude for cosine similarity
    std::atomic<bool> stopTraversal{false};
    int kmerSize = 32;  // set from params

    // Active seed set during traversal (map of seedIndex -> count)
    std::unordered_map<uint64_t, uint32_t> activeSeedIndices;
    
    // MGSR index data
    ::capnp::List<SeedInfo>::Reader seedInfo;
    ::capnp::List<NodeChanges>::Reader perNodeChanges;
    ::capnp::List<LiteNode>::Reader liteNodes;  // For node ID lookups
    
    // Root node pointer for traversal
    panmapUtils::LiteNode* root = nullptr;
    
    // Optional: Full tree for verification mode (nullptr if not using verification)
    panmanUtils::Tree* fullTree = nullptr;
};

// MGSR-style global genome state (single instance, modified in-place with backtracking)
struct MgsrGenomeState {
    // Position-based seed tracking (like MGSR)
    std::map<uint64_t, uint64_t> positionMap;  // position -> seedIndex
    absl::flat_hash_map<size_t, std::vector<std::map<uint64_t, uint64_t>::iterator>> hashToPositionMap;  // hash -> list of position iterators
    
    // Cached similarity metrics for current node (updated incrementally)
    int64_t currentJaccardNumerator = 0;
    int64_t currentJaccardDenominator = 0;  // Cached Jaccard denominator (union of read and genome)
    int64_t currentWeightedJaccardNumerator = 0;
    int64_t currentWeightedJaccardDenominator = 0;  // Cached weighted Jaccard denominator
    double currentCosineNumerator = 0.0;
    double currentGenomeMagnitude = 0.0;  // Cached genome magnitude for cosine (sqrt of magnitudeSquared)
    double currentGenomeMagnitudeSquared = 0.0;  // Cached squared magnitude (avoids rounding errors)
    size_t currentPresenceIntersectionCount = 0;  // Cached presence/absence intersection count
    size_t currentPresenceUnionCount = 0;  // Cached presence/absence union count
    absl::flat_hash_map<size_t, int64_t> currentSeedCounts;  // For frequency-based metrics
    absl::flat_hash_set<size_t> currentUniqueSeeds;          // For presence-based metrics
    
    void reset() {
        positionMap.clear();
        hashToPositionMap.clear();
        currentJaccardNumerator = 0;
        currentJaccardDenominator = 0;
        currentWeightedJaccardNumerator = 0;
        currentWeightedJaccardDenominator = 0;
        currentCosineNumerator = 0.0;
        currentGenomeMagnitude = 0.0;
        currentGenomeMagnitudeSquared = 0.0;
        currentPresenceIntersectionCount = 0;
        currentPresenceUnionCount = 0;
        currentSeedCounts.clear();
        currentUniqueSeeds.clear();
    }
    
    size_t estimateMemoryUsage() const {
        return positionMap.size() * (sizeof(uint64_t) + sizeof(uint64_t)) +
               hashToPositionMap.size() * sizeof(size_t) * 8 +  // Rough estimate for vector overhead
               currentSeedCounts.size() * (sizeof(size_t) + sizeof(int64_t)) +
               currentUniqueSeeds.size() * sizeof(size_t);
    }
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
  std::string maxHitsNodeId;  // Changed from Node* to string ID
  std::vector<std::string> tiedMaxHitsNodeIds;  // Changed from vector<Node*> to vector<string>

  // Raw Seed Match Score (New)
  int64_t bestRawSeedMatchScore = 0;
  std::string bestRawSeedMatchNodeId;  // Changed from Node* to string ID
  std::vector<std::string> tiedRawSeedMatchNodeIds;  // Changed from vector<Node*> to vector<string>

  // Jaccard-based results (current "weighted" Jaccard)
  double bestJaccardScore = 0.0; // This is the weighted Jaccard
  std::string bestJaccardNodeId;  // Changed from Node* to string ID
  std::vector<std::string> tiedJaccardNodeIds;  // Changed from vector<Node*> to vector<string>

  // Weighted Jaccard results (separate tracking)
  double bestWeightedJaccardScore = 0.0;
  std::string bestWeightedJaccardNodeId;
  std::vector<std::string> tiedWeightedJaccardNodeIds;

  // Jaccard Index (Presence/Absence) (New)
  double bestJaccardPresenceScore = 0.0;
  std::string bestJaccardPresenceNodeId;  // Changed from Node* to string ID
  std::vector<std::string> tiedJaccardPresenceNodeIds;  // Changed from vector<Node*> to vector<string>

  // Cosine-based results
  double bestCosineScore = 0.0;
  std::string bestCosineNodeId;  // Changed from Node* to string ID
  std::vector<std::string> tiedCosineNodeIds;  // Changed from vector<Node*> to vector<string>

  // Weighted result
  double bestWeightedScore = 0.0;
  std::string bestWeightedNodeId;  // Changed from Node* to string ID
  std::vector<std::string> tiedWeightedNodeIds;  // Changed from vector<Node*> to vector<string>

  // Current tracking values
  int64_t hitsInThisGenome = 0;
  int64_t currentWeightedJaccardNumerator = 0;    // For weighted Jaccard (min frequencies)
  int64_t currentJaccardNumerator = 0;            // Legacy field, now used for regular Jaccard  
  int64_t currentJaccardDenominator = 0;
  double currentCosineNumerator = 0.0;
  double currentCosineDenominator = 0.0;
  double currentGenomeMagnitudeSquared = 0.0;     // Sum of (genomeCount^2) for ALL seeds
  absl::flat_hash_map<size_t, int64_t> currentAllGenomeSeedCounts;  // ALL genome seeds (authoritative)
  
  // Added for placement - map of node seeds without dependency on StateManager
  // This allows us to track seeds directly in the placement result
  // Format: nodeId -> position -> seed
  std::unordered_map<std::string, std::unordered_map<int64_t, seeding::seed_t>> nodeSeedMap;

  // Performance metrics
  int64_t totalReadsProcessed = 0;
  double totalTimeSeconds = 0.0;
  
  // Helper methods - updated to use node IDs instead of pointers
  void resetCurrentNodeState() {
    // Reset per-node tracking values (should be called before evaluating each node)
    hitsInThisGenome = 0;
    currentWeightedJaccardNumerator = 0;
    currentJaccardNumerator = 0;
    currentJaccardDenominator = 0;
    currentCosineNumerator = 0.0;
    currentCosineDenominator = 0.0;
    currentGenomeMagnitudeSquared = 0.0;
    currentAllGenomeSeedCounts.clear();
    // Note: nodeSeedMap is preserved across nodes as it stores the genome being built
  }
  
  void updateHitsScore(const std::string& nodeId, int64_t hits); // This is for the existing "maxHitsInAnyGenome"
  void updateRawSeedMatchScore(const std::string& nodeId, int64_t score); // New
  void updateJaccardScore(const std::string& nodeId, double score); // This is for regular Jaccard
  void updateWeightedJaccardScore(const std::string& nodeId, double score); // For weighted Jaccard
  void updateJaccardPresenceScore(const std::string& nodeId, double score); // New
  void updateCosineScore(const std::string& nodeId, double score);
  void updateWeightedScore(const std::string& nodeId, double score, double scale);
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

// NEW: LiteTree-based placement (avoids loading full panman until after placement)
void placeLite(PlacementResult &result, 
               panmapUtils::LiteTree *liteTree,
               ::MGSRIndex::Reader &mgsrIndex, 
               const std::string &reads1,
               const std::string &reads2,
               std::vector<std::vector<seeding::seed_t>>& readSeeds,
               std::vector<std::string>& readSequences,
               std::vector<std::string>& readNames,
               std::vector<std::string>& readQuals,
               std::string &outputPath,
               const std::string &indexPath,
               const std::string &debug_node_id_param,
               bool verify_scores = false,
               panmanUtils::Tree *fullTree = nullptr);  // Optional: only needed for verification


// NEW: Declaration for placement summary dump function
void dumpPlacementSummary(const PlacementResult& result, const std::string& outputFilename);

// NEW: Declaration for kmer debug dump function 
void dumpKmerDebugData(
    const std::vector<std::string>& readSequences,
    int k,
    const std::string& outputFilename);

// Consolidated seed processing functions (MGSR compatible with LiteTree)
void processSeedOperation(
    const std::string& nodeId,  // Changed from Node* to string ID
    state::StateManager& stateManager,
    PlacementGlobalState& state,
    PlacementResult& result,
    absl::flat_hash_set<size_t>& uniqueSeedHashes,
    int64_t pos,
    bool isRemoval,
    bool isAddition,
    const seeding::seed_t* newSeed = nullptr);

seeding::seed_t createAndProcessSeed(
    const std::string& nodeId,  // Changed from Node* to string ID
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
    const std::string& nodeId,  // Changed from Node* to string ID
    PlacementGlobalState& state,
    PlacementResult& result,
    const TraversalParams& params,
    absl::flat_hash_set<size_t>& uniqueSeedHashes);

// LiteTree-specific helper function
std::tuple<double, double, size_t> calculateAndUpdateScoresLite(
    const std::string& nodeId,
    PlacementGlobalState& state,
    PlacementResult& result,
    const TraversalParams& params,
    absl::flat_hash_set<size_t>& uniqueSeedHashes,
    uint32_t nodeChangeIndex);

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