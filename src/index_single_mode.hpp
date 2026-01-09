#pragma once

#include "panmanUtils.hpp"
#include "capnp/message.h"
#include "capnp/serialize-packed.h"
#include "index_lite.capnp.h"
#include "panmap_utils.hpp"
#include "seeding.hpp"
#include <tbb/task_arena.h>
#include <tbb/task_group.h>
#include <tbb/concurrent_hash_map.h>
#include <tbb/concurrent_vector.h>
#include <tbb/spin_mutex.h>
#include <atomic>
#include <absl/container/btree_set.h>
#include <absl/container/flat_hash_map.h>

namespace index_single_mode {

// Type alias for ordered set - btree_set has better cache locality than std::set
using SyncmerSet = absl::btree_set<uint64_t>;

// Sparse map types for memory-efficient storage
// Instead of vector<optional<T>> sized to genome length, store only actual values
using SyncmerMap = absl::flat_hash_map<uint64_t, seeding::rsyncmer_t>;
using KminmerMap = absl::flat_hash_map<uint64_t, uint64_t>;

// ============================================================================
// NodeSeedDelta: Stores just the seed changes for a node (not full counts)
// ============================================================================
struct NodeSeedDelta {
    std::vector<uint64_t> addedHashes;      // Seeds added at this node
    std::vector<uint64_t> deletedHashes;    // Seeds deleted at this node  
    std::vector<std::pair<uint64_t, uint64_t>> substitutedHashes; // <oldHash, newHash>
    
    void clear() {
        addedHashes.clear();
        deletedHashes.clear();
        substitutedHashes.clear();
    }
    
    bool empty() const {
        return addedHashes.empty() && deletedHashes.empty() && substitutedHashes.empty();
    }
};

// ============================================================================
// BuildState: Encapsulates all mutable state for parallel DFS traversal
// Memory-optimized: uses sparse maps instead of dense vectors
// ============================================================================
struct BuildState {
    // Block sequence state (mutated during traversal)
    panmapUtils::BlockSequences blockSequences;
    std::vector<char> blockExistsDelayed;
    std::vector<char> blockStrandDelayed;
    
    // Gap and inversion tracking
    std::map<uint64_t, uint64_t> gapMap;
    std::unordered_set<uint64_t> invertedBlocks;
    
    // Genome extent tracking: first and last non-gap scalar positions
    // Seeds in flank regions (before firstNonGapScalar or after lastNonGapScalar)
    // should NOT be deleted when they become gaps - they're missing data, not true gaps
    uint64_t firstNonGapScalar = UINT64_MAX;
    uint64_t lastNonGapScalar = 0;
    
    // Syncmer state - SPARSE: only stores positions with actual syncmers
    // Memory: O(num_syncmers) instead of O(genome_length)
    SyncmerMap refOnSyncmers;     // position -> rsyncmer_t (sparse)
    SyncmerSet refOnSyncmersMap;  // btree_set of positions for ordered iteration
    std::unordered_map<uint32_t, std::unordered_set<uint64_t>> blockOnSyncmers;
    
    // K-minmer state - SPARSE: only stores positions with actual kminmers
    KminmerMap refOnKminmers;     // position -> kminmer hash (sparse)
    
    // Running seed hash counts for computing node changes on-the-fly
    // Each chunk maintains its own copy; modified during DFS traversal
    std::unordered_map<uint64_t, int64_t> runningCounts;
    
    // Shared storage for final node changes - indexed by dfsIndex (thread-safe writes)
    // Each node writes to its own slot, so no synchronization needed
    std::shared_ptr<std::vector<std::vector<std::tuple<uint64_t, int64_t, int64_t>>>> nodeChanges;
    
    // Shared storage for genome metrics - computed during main DFS (indexed by dfsIndex)
    std::shared_ptr<std::vector<double>> genomeMagnitudeSquared;
    std::shared_ptr<std::vector<uint64_t>> genomeUniqueSeedCount;
    std::shared_ptr<std::vector<int64_t>> genomeTotalSeedFrequency;
    
    // IDF tracking: For each seed hash, count how many leaf genomes contain it
    // Thread-local map, merged at end. Key = seed hash, Value = count of leaf genomes
    std::shared_ptr<std::unordered_map<uint64_t, uint32_t>> leafSeedGenomeCounts;
    
    // Track which DFS indices are leaf nodes (for IDF computation)
    std::shared_ptr<std::vector<bool>> isLeafNode;
    
    // Default constructor
    BuildState() : nodeChanges(std::make_shared<std::vector<std::vector<std::tuple<uint64_t, int64_t, int64_t>>>>()) {}
    
    // Move constructor and assignment
    BuildState(BuildState&&) = default;
    BuildState& operator=(BuildState&&) = default;
    
    // Copy constructor - shares nodeChanges storage (doesn't deep copy it)
    // But copies runningCounts since each chunk needs its own copy
    BuildState(const BuildState& other) = default;
    BuildState& operator=(const BuildState&) = default;
    
    // Clone method for explicit copying
    BuildState clone() const { return *this; }
    
    // Helper methods for sparse access
    bool hasSyncmer(uint64_t pos) const { return refOnSyncmers.contains(pos); }
    bool hasKminmer(uint64_t pos) const { return refOnKminmers.contains(pos); }
    
    // Check if a position is in the genome's valid extent (not in flank regions)
    bool isInGenomeExtent(uint64_t scalarPos) const {
        return scalarPos >= firstNonGapScalar && scalarPos <= lastNonGapScalar;
    }
};

struct BacktrackInfo {
  std::vector<std::tuple<uint32_t, bool, bool, bool, bool>> blockMutationRecord;
  std::vector<std::tuple<panmapUtils::Coordinate, char, char>> nucMutationRecord;
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>> gapRunBacktracks;
  std::vector<std::pair<uint64_t, bool>> invertedBlocksBacktracks;
  std::vector<std::tuple<uint64_t, panmapUtils::seedChangeType, seeding::rsyncmer_t>> refOnSyncmersChangeRecord;
  std::vector<std::tuple<uint64_t, uint64_t, panmapUtils::seedChangeType>> blockOnSyncmersChangeRecord;
  std::vector<std::tuple<uint64_t, panmapUtils::seedChangeType, uint64_t>> refOnKminmersChangeRecord;
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>> gapRunBlockInversionBacktracks;
  // Running count changes: (hash, delta) where delta is +1 for added, -1 for deleted
  std::vector<std::pair<uint64_t, int64_t>> runningCountChanges;
  
  // Genome extent before this node's mutations (for backtracking)
  uint64_t prevFirstNonGapScalar = UINT64_MAX;
  uint64_t prevLastNonGapScalar = 0;

  void clear() {
    blockMutationRecord.clear();
    nucMutationRecord.clear();
    gapRunBacktracks.clear();
    invertedBlocksBacktracks.clear();
    refOnSyncmersChangeRecord.clear();
    blockOnSyncmersChangeRecord.clear();
    refOnKminmersChangeRecord.clear();
    gapRunBlockInversionBacktracks.clear();
    runningCountChanges.clear();
    prevFirstNonGapScalar = UINT64_MAX;
    prevLastNonGapScalar = 0;
  }
};

// Compute genome extent from gapMap
// Returns (firstNonGapScalar, lastNonGapScalar)
// If genome is all gaps, returns (UINT64_MAX, 0)
inline std::pair<uint64_t, uint64_t> computeExtentFromGapMap(
    const std::map<uint64_t, uint64_t>& gapMap,
    uint64_t lastScalarCoord) {
  uint64_t firstNonGap = 0;
  uint64_t lastNonGap = lastScalarCoord;
  
  // Check for leading gap run (starts at position 0)
  auto it = gapMap.find(0);
  if (it != gapMap.end()) {
    // There's a gap run starting at 0
    if (it->second >= lastScalarCoord) {
      // Entire genome is gaps
      return {UINT64_MAX, 0};
    }
    firstNonGap = it->second + 1;
  }
  
  // Check for trailing gap run (ends at lastScalarCoord)
  if (!gapMap.empty()) {
    auto lastIt = std::prev(gapMap.end());
    if (lastIt->second == lastScalarCoord) {
      lastNonGap = lastIt->first - 1;
    }
  }
  
  return {firstNonGap, lastNonGap};
}

void updateGapMapStep(
    std::map<uint64_t, uint64_t>& gapMap,
    uint64_t startPos,
    uint64_t endPos,
    bool toGap,
    std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& backtrack,
    std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& gapMapUpdates,
    bool recordGapMapUpdates
);

void updateGapMap(
  panmanUtils::Node *node,
  size_t dfsIndex,
  std::map<uint64_t, uint64_t>& gapMap,
  const std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& updates,
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& backtrack,
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& gapMapUpdates
);

void invertGapMap(
  std::map<uint64_t, uint64_t>& gapMap,
  const std::pair<uint64_t, uint64_t>& invertRange,
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& backtrack,
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& gapMapUpdates
);

void revertGapMapInversions(
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& gapMapBlocksBacktracks,
  std::map<uint64_t, uint64_t>& gapMap
);

void makeCoordIndex(
  std::map<uint64_t, uint64_t>& degapCoordIndex,
  std::map<uint64_t, uint64_t>& regapCoordIndex,
  const std::map<uint64_t, uint64_t>& gapMap,
  uint64_t lastScalarCoord
);

uint64_t degapGlobal(const uint64_t& globalCoord, const std::map<uint64_t, uint64_t>& degapCoordsIndex);

uint64_t regapGlobal(const uint64_t& localCoord, const std::map<uint64_t, uint64_t>& regapCoordsIndex);

int open_file(const std::string& path);


class IndexBuilder {
  public:
    // capnp object
    ::capnp::MallocMessageBuilder outMessage;
    LiteIndex::Builder indexBuilder;
    

    // tree pointer
    panmanUtils::Tree *T;
  
    // syncmer and k-min-mer objects (legacy - used by sequential path)
    std::unordered_map<uint32_t, std::unordered_set<uint64_t>> blockOnSyncmers;
    SyncmerMap refOnSyncmers;     // SPARSE: position -> rsyncmer_t
    SyncmerSet refOnSyncmersMap;  // btree_set for ordered iteration

    // Use concurrent_vector for thread-safe growth without invalidating references
    tbb::concurrent_vector<seeding::uniqueKminmer_t> uniqueKminmers;
    KminmerMap refOnKminmers;     // SPARSE: position -> kminmer hash
    std::unordered_map<seeding::uniqueKminmer_t, uint64_t> kminmerToUniqueIndex;

    std::unordered_map<std::string, uint32_t> nodeToDfsIndex;
    
    // For tracking seed counts per node (needed for LiteIndex format)
    std::vector<std::unordered_map<uint64_t, int64_t>> nodeSeedCounts;
    
    // Thread-safe containers for parallel building
    tbb::spin_mutex uniqueKminmersMutex_;
    tbb::spin_mutex nodeToDfsIndexMutex_;
    std::atomic<uint64_t> processedNodes_{0};
    
    // Memory management for parallel cloning
    std::atomic<int> activeClones_{0};  // Track number of active state clones
    int maxConcurrentClones_{0};        // Maximum allowed concurrent clones (set based on memory)
    
    // Instrumentation stats for parallel builds
    std::atomic<uint64_t> totalClonesCreated_{0};    // Cumulative clones created
    std::atomic<uint64_t> parallelForkPoints_{0};    // Times we forked into parallel tasks
    std::atomic<uint64_t> sequentialFallbacks_{0};   // Times we fell back to sequential (clone limit)
    std::atomic<int> peakActiveClones_{0};           // Peak concurrent clones
    mutable std::chrono::steady_clock::time_point lastProgressLog_{};
    std::chrono::steady_clock::time_point buildStartTime_{};  // For nodes/sec calculation
    
    // Pre-computed subtree sizes for DFS index assignment
    std::unordered_map<std::string, uint64_t> subtreeSizes_;
    
    // Thread-safe empty nodes tracking
    tbb::spin_mutex emptyNodesMutex_;

    IndexBuilder(panmanUtils::Tree *T, int k, int s, int t, int l, bool openSyncmer, int flankMaskBp = 250) 
      : outMessage(), indexBuilder(outMessage.initRoot<LiteIndex>()), T(T), k_(k), s_(s), l_(l), flankMaskBp_(flankMaskBp)
    {
      indexBuilder.setK(k);
      indexBuilder.setS(s);
      indexBuilder.setT(t);
      indexBuilder.setL(l);
      indexBuilder.setOpen(openSyncmer);
      indexBuilder.setVersion(3);
      nodeToDfsIndex.reserve(T->allNodes.size());
      nodeSeedCounts.resize(T->allNodes.size());
    }

    void buildIndex();
    void buildIndexParallel(int numThreads = 0);  // 0 = auto-detect

    void writeIndex(const std::string& path, int numThreads = 0);  // 0 = auto-detect

    // Getters for parameters (for testing)
    int getK() const { return k_; }
    int getS() const { return s_; }
    int getL() const { return l_; }

  private:
    int k_, s_, l_;
    int flankMaskBp_;  // Hard mask first/last N bp at genome ends
    
  public:
    int getFlankMaskBp() const { return flankMaskBp_; }
    
    // Compute subtree size for a node (recursive)
    uint64_t computeSubtreeSize(panmanUtils::Node* node);
    
    // Sequential DFS helper (original)
    void buildIndexHelper(
      panmanUtils::Node *node,
      std::unordered_set<std::string_view>& emptyNodes,
      panmapUtils::BlockSequences &blockSequences,
      std::vector<char> &blockExistsDelayed,
      std::vector<char> &blockStrandDelayed,
      panmapUtils::GlobalCoords &globalCoords,
      std::map<uint64_t, uint64_t> &gapMap,
      std::unordered_set<uint64_t> &invertedBlocks,
      uint64_t &dfsIndex
    );
    
    // Process a single node without recursion or backtracking (for parallel path processing)
    void processSingleNodeNoBacktrack(
      panmanUtils::Node *node,
      std::unordered_set<std::string_view>& emptyNodes,
      panmapUtils::BlockSequences &blockSequences,
      std::vector<char> &blockExistsDelayed,
      std::vector<char> &blockStrandDelayed,
      panmapUtils::GlobalCoords &globalCoords,
      std::map<uint64_t, uint64_t> &gapMap,
      std::unordered_set<uint64_t> &invertedBlocks,
      uint64_t dfsIndex
    );
    
    // Sequential subtree processing (no cloning, uses passed-by-reference state)
    void processSubtreeSequential(
      panmanUtils::Node *node,
      BuildState& state,
      panmapUtils::GlobalCoords &globalCoords,
      std::unordered_set<std::string_view>& localEmptyNodes,
      uint64_t dfsIndex,
      BacktrackInfo* backtrackInfo = nullptr
    );
    
    // Parallel subtree processing - spawns parallel tasks at fork points
    void processSubtreeParallel(
      panmanUtils::Node *node,
      BuildState& state,
      panmapUtils::GlobalCoords &globalCoords,
      std::unordered_set<std::string_view>& localEmptyNodes,
      uint64_t dfsIndex,
      tbb::spin_mutex& emptyNodesMutex,
      size_t parallelThreshold,
      size_t depth = 0  // Recursion depth for instrumentation
    );
    
    // Parallel DFS helper (clones state only at top levels)
    void buildIndexHelperParallel(
      panmanUtils::Node *node,
      BuildState state,
      panmapUtils::GlobalCoords &globalCoords,
      std::unordered_set<std::string_view>& emptyNodes,
      uint64_t dfsIndex,
      tbb::task_group& taskGroup,
      int depth
    );
    
    // Process a single node and compute its seed counts
    void processNode(
      panmanUtils::Node *node,
      BuildState& state,
      panmapUtils::GlobalCoords &globalCoords,
      std::unordered_set<std::string_view>& localEmptyNodes,
      uint64_t dfsIndex,
      BacktrackInfo* backtrackInfo = nullptr,
      bool skipNodeChanges = false  // If true, don't compute/record nodeChanges (for path walking in parallel)
    );
    
    // Backtrack changes made by processNode using recorded changes
    void backtrackNode(
      BuildState& state,
      const BacktrackInfo& backtrackInfo
    );

    std::vector<panmapUtils::NewSyncmerRange> computeNewSyncmerRangesJump(
      panmanUtils::Node* node,
      size_t dfsIndex,
      const panmapUtils::BlockSequences& blockSequences,
      const std::vector<char>& blockExistsDelayed,
      const std::vector<char>& blockStrandDelayed,
      const panmapUtils::GlobalCoords& globalCoords,
      const std::map<uint64_t, uint64_t>& gapMap,
      std::vector<std::pair<panmapUtils::Coordinate, panmapUtils::Coordinate>>& localMutationRanges,
      std::vector<std::tuple<uint64_t, uint64_t, panmapUtils::seedChangeType>>& blockOnSyncmersBacktracks,
      const SyncmerMap& refOnSyncmers,
      std::unordered_map<uint32_t, std::unordered_set<uint64_t>>& blockOnSyncmers,
      uint64_t firstNonGapScalar = 0,
      uint64_t lastNonGapScalar = UINT64_MAX
    );

    std::vector<panmapUtils::NewSyncmerRange> computeNewSyncmerRangesWalk(
      panmanUtils::Node* node,
      size_t dfsIndex,
      const panmapUtils::BlockSequences& blockSequences,
      const std::vector<char>& blockExistsDelayed,
      const std::vector<char>& blockStrandDelayed,
      const panmapUtils::GlobalCoords& globalCoords,
      const std::map<uint64_t, uint64_t>& gapMap,
      std::vector<std::pair<panmapUtils::Coordinate, panmapUtils::Coordinate>>& localMutationRanges,
      std::vector<std::tuple<uint64_t, uint64_t, panmapUtils::seedChangeType>>& blockOnSyncmersBacktracks,
      const SyncmerMap& refOnSyncmers,
      std::unordered_map<uint32_t, std::unordered_set<uint64_t>>& blockOnSyncmers
    );

    std::vector<std::pair<SyncmerSet::iterator, SyncmerSet::iterator>> computeNewKminmerRanges(
      std::vector<std::tuple<uint64_t, panmapUtils::seedChangeType, seeding::rsyncmer_t>>& refOnSyncmersChangeRecord,
      BuildState& state,
      const uint64_t dfsIndex
    );
    
    // Overload for sequential path (uses member state)
    std::vector<std::pair<SyncmerSet::iterator, SyncmerSet::iterator>> computeNewKminmerRanges(
      std::vector<std::tuple<uint64_t, panmapUtils::seedChangeType, seeding::rsyncmer_t>>& refOnSyncmersChangeRecord,
      const uint64_t dfsIndex
    );
};

// end of namespace index_single_mode
}
