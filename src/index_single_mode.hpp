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

namespace index_single_mode {

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
// ============================================================================
struct BuildState {
    // Block sequence state (mutated during traversal)
    panmapUtils::BlockSequences blockSequences;
    std::vector<char> blockExistsDelayed;
    std::vector<char> blockStrandDelayed;
    
    // Gap and inversion tracking
    std::map<uint64_t, uint64_t> gapMap;
    std::unordered_set<uint64_t> invertedBlocks;
    
    // Syncmer state
    std::vector<std::optional<seeding::rsyncmer_t>> refOnSyncmers;
    std::set<uint64_t> refOnSyncmersMap;
    std::unordered_map<uint32_t, std::unordered_set<uint64_t>> blockOnSyncmers;
    
    // K-minmer state
    std::vector<std::optional<uint64_t>> refOnKminmers;
    
    // Per-node seed deltas - shared across clones (NOT copied during clone)
    // Stores only the changes at each node, not full counts
    std::shared_ptr<std::vector<NodeSeedDelta>> nodeSeedDeltas;
    
    // Default constructor
    BuildState() : nodeSeedDeltas(std::make_shared<std::vector<NodeSeedDelta>>()) {}
    
    // Move constructor and assignment
    BuildState(BuildState&&) = default;
    BuildState& operator=(BuildState&&) = default;
    
    // Copy constructor - shares nodeSeedDeltas (doesn't deep copy it)
    BuildState(const BuildState& other) = default;
    BuildState& operator=(const BuildState&) = default;
    
    // Clone method for explicit copying
    BuildState clone() const { return *this; }
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

  void clear() {
    blockMutationRecord.clear();
    nucMutationRecord.clear();
    gapRunBacktracks.clear();
    invertedBlocksBacktracks.clear();
    refOnSyncmersChangeRecord.clear();
    blockOnSyncmersChangeRecord.clear();
    refOnKminmersChangeRecord.clear();
    gapRunBlockInversionBacktracks.clear();
  }
};

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
    std::vector<std::optional<seeding::rsyncmer_t>> refOnSyncmers;
    std::set<uint64_t> refOnSyncmersMap;

    // Use concurrent_vector for thread-safe growth without invalidating references
    tbb::concurrent_vector<seeding::uniqueKminmer_t> uniqueKminmers;
    std::vector<std::optional<uint64_t>> refOnKminmers;
    std::unordered_map<seeding::uniqueKminmer_t, uint64_t> kminmerToUniqueIndex;

    std::unordered_map<std::string, uint32_t> nodeToDfsIndex;
    
    // For tracking seed counts per node (needed for LiteIndex format)
    std::vector<std::unordered_map<uint64_t, int64_t>> nodeSeedCounts;
    
    // Thread-safe containers for parallel building
    tbb::spin_mutex uniqueKminmersMutex_;
    tbb::spin_mutex nodeToDfsIndexMutex_;
    std::atomic<uint64_t> processedNodes_{0};
    
    // Pre-computed subtree sizes for DFS index assignment
    std::unordered_map<std::string, uint64_t> subtreeSizes_;
    
    // Thread-safe empty nodes tracking
    tbb::spin_mutex emptyNodesMutex_;

    IndexBuilder(panmanUtils::Tree *T, int k, int s, int t, int l, bool openSyncmer) 
      : outMessage(), indexBuilder(outMessage.initRoot<LiteIndex>()), T(T), k_(k), s_(s), l_(l)
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

    void writeIndex(const std::string& path);

    // Getters for parameters (for testing)
    int getK() const { return k_; }
    int getS() const { return s_; }
    int getL() const { return l_; }

  private:
    int k_, s_, l_;
    
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
      size_t parallelThreshold
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
      BacktrackInfo* backtrackInfo = nullptr
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
      const std::vector<std::optional<seeding::rsyncmer_t>>& refOnSyncmers,
      std::unordered_map<uint32_t, std::unordered_set<uint64_t>>& blockOnSyncmers
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
      const std::vector<std::optional<seeding::rsyncmer_t>>& refOnSyncmers,
      std::unordered_map<uint32_t, std::unordered_set<uint64_t>>& blockOnSyncmers
    );

    std::vector<std::pair<std::set<uint64_t>::iterator, std::set<uint64_t>::iterator>> computeNewKminmerRanges(
      std::vector<std::tuple<uint64_t, panmapUtils::seedChangeType, seeding::rsyncmer_t>>& refOnSyncmersChangeRecord,
      BuildState& state,
      const uint64_t dfsIndex
    );
    
    // Overload for sequential path (uses member state)
    std::vector<std::pair<std::set<uint64_t>::iterator, std::set<uint64_t>::iterator>> computeNewKminmerRanges(
      std::vector<std::tuple<uint64_t, panmapUtils::seedChangeType, seeding::rsyncmer_t>>& refOnSyncmersChangeRecord,
      const uint64_t dfsIndex
    );
};

// end of namespace index_single_mode
}
