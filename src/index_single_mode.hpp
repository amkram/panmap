#pragma once

#include "panmanUtils.hpp"
#include "logging.hpp"
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
#include <array>
#include <cstdint>
#include <string>

namespace index_single_mode {

// A small uncompressed header prepended to the .idx so cache-validation can read the
// seeding params without decompressing the (up to GB-scale) payload -- previously the
// validation step fully decompressed the index just to read six fields, then the
// placement step decompressed it a second time.
constexpr uint32_t kIndexMagic = 0x31494D50u;  // "PMI1" little-endian
constexpr uint32_t kIndexHeaderVersion = 1;
constexpr size_t kIndexHeaderSize = 32;
struct IndexParamsHeader {
    int32_t k = 0, s = 0, t = 0, l = 0;
    bool hpc = false, open = false;
};
std::array<uint8_t, kIndexHeaderSize> encodeIndexHeader(const IndexParamsHeader& p);
// Reads the header from the start of `path`. Returns false if absent / not this format.
bool readIndexHeader(const std::string& path, IndexParamsHeader& out);

// btree_set has better cache locality than std::set
using SyncmerSet = absl::btree_set<uint64_t>;

// Sparse: store only actual values, not a vector<optional<T>> sized to genome length
using SyncmerMap = absl::flat_hash_map<uint64_t, seeding::rsyncmer_t>;
using KminmerMap = absl::flat_hash_map<uint64_t, uint64_t>;

struct NodeSeedDelta {
    std::vector<uint64_t> addedHashes;
    std::vector<uint64_t> deletedHashes;
    std::vector<std::pair<uint64_t, uint64_t>> substitutedHashes;  // <oldHash, newHash>

    void clear() {
        addedHashes.clear();
        deletedHashes.clear();
        substitutedHashes.clear();
    }

    bool empty() const { return addedHashes.empty() && deletedHashes.empty() && substitutedHashes.empty(); }
};

struct BuildState {
    panmapUtils::BlockSequences blockSequences;
    std::vector<char> blockExistsDelayed;
    std::vector<char> blockStrandDelayed;

    std::map<uint64_t, uint64_t> gapMap;
    std::unordered_set<uint64_t> invertedBlocks;

    // first/last non-gap scalar positions. Seeds in flank regions (before
    // firstNonGapScalar or after lastNonGapScalar) are missing data, not true gaps,
    // so must not be deleted when they become gaps. Permissive by default; only
    // restricted when extentGuard is on.
    uint64_t firstNonGapScalar = 0;
    uint64_t lastNonGapScalar = UINT64_MAX;

    // Sparse: O(num_syncmers) memory, not O(genome_length)
    SyncmerMap refOnSyncmers;     // position -> rsyncmer_t
    SyncmerSet refOnSyncmersMap;  // btree_set of positions for ordered iteration
    std::unordered_map<uint32_t, std::unordered_set<uint64_t>> blockOnSyncmers;

    KminmerMap refOnKminmers;  // position -> kminmer hash

    // Running seed hash counts for node changes. Each chunk keeps its own copy.
    std::unordered_map<uint64_t, int64_t> runningCounts;

    // Shared, indexed by dfsIndex. Each node writes its own slot, so no locking needed.
    std::shared_ptr<std::vector<std::vector<std::tuple<uint64_t, int64_t, int64_t>>>> nodeChanges;

    BuildState() : nodeChanges(std::make_shared<std::vector<std::vector<std::tuple<uint64_t, int64_t, int64_t>>>>()) {}

    BuildState(BuildState&&) = default;
    BuildState& operator=(BuildState&&) = default;

    // Copy gives each chunk its own runningCounts
    BuildState(const BuildState& other) = default;
    BuildState& operator=(const BuildState&) = default;
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
    uint64_t prevFirstNonGapScalar = 0;
    uint64_t prevLastNonGapScalar = UINT64_MAX;

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
        prevFirstNonGapScalar = 0;
        prevLastNonGapScalar = UINT64_MAX;
    }
};

void makeCoordIndex(std::map<uint64_t, uint64_t>& degapCoordIndex,
                    std::map<uint64_t, uint64_t>& regapCoordIndex,
                    const std::map<uint64_t, uint64_t>& gapMap,
                    uint64_t lastScalarCoord);

uint64_t degapGlobal(const uint64_t& globalCoord, const std::map<uint64_t, uint64_t>& degapCoordsIndex);

uint64_t regapGlobal(const uint64_t& localCoord, const std::map<uint64_t, uint64_t>& regapCoordsIndex);

int open_file(const std::string& path);

class IndexBuilder {
   public:
    ::capnp::MallocMessageBuilder outMessage;
    LiteIndex::Builder indexBuilder;

    panmanUtils::Tree* T;

    // Used by the sequential path
    std::unordered_map<uint32_t, std::unordered_set<uint64_t>> blockOnSyncmers;
    SyncmerMap refOnSyncmers;     // position -> rsyncmer_t
    SyncmerSet refOnSyncmersMap;  // btree_set for ordered iteration

    // concurrent_vector: thread-safe growth without invalidating references
    tbb::concurrent_vector<seeding::uniqueKminmer_t> uniqueKminmers;
    KminmerMap refOnKminmers;  // position -> kminmer hash
    std::unordered_map<seeding::uniqueKminmer_t, uint64_t> kminmerToUniqueIndex;

    std::unordered_map<std::string, uint32_t> nodeToDfsIndex;

    // Seed counts per node (needed for LiteIndex format)
    std::vector<std::unordered_map<uint64_t, int64_t>> nodeSeedCounts;

    tbb::spin_mutex uniqueKminmersMutex_;
    tbb::spin_mutex nodeToDfsIndexMutex_;
    std::atomic<uint64_t> processedNodes_{0};

    // Instrumentation stat for parallel builds
    std::atomic<uint64_t> totalClonesCreated_{0};
    mutable std::chrono::steady_clock::time_point lastProgressLog_{};
    std::chrono::steady_clock::time_point buildStartTime_{};  // For nodes/sec calculation
    std::unique_ptr<output::ProgressBar> buildProgressBar_;

    // Pre-computed subtree sizes for DFS index assignment
    std::unordered_map<std::string, uint64_t> subtreeSizes_;

    tbb::spin_mutex emptyNodesMutex_;

    IndexBuilder(panmanUtils::Tree* T,
                 int k,
                 int s,
                 int t,
                 int l,
                 bool openSyncmer,
                 int flankMaskBp = 250,
                 bool hpc = false,
                 bool imputeAmb = false,
                 bool extentGuard = false)
        : outMessage(),
          indexBuilder(outMessage.initRoot<LiteIndex>()),
          T(T),
          k_(k),
          s_(s),
          l_(l),
          flankMaskBp_(flankMaskBp),
          hpc_(hpc),
          imputeAmb_(imputeAmb),
          extentGuard_(extentGuard) {
        indexBuilder.setK(k);
        indexBuilder.setS(s);
        indexBuilder.setT(t);
        indexBuilder.setL(l);
        indexBuilder.setOpen(openSyncmer);
        indexBuilder.setHpc(hpc);
        indexBuilder.setFormatVersion(panmapUtils::INDEX_FORMAT_VERSION);
        nodeToDfsIndex.reserve(T->allNodes.size());
        nodeSeedCounts.resize(T->allNodes.size());
    }

    void buildIndex();
    void buildIndexParallel(int numThreads = 0);  // 0 = auto-detect

    // Compute 4x4 substitution spectrum from tree mutations and store in index
    void computeSubstitutionSpectrum();

    void writeIndex(const std::string& path, int numThreads = 0, int zstdLevel = 7);

    // Parameter getters (for testing)
    int getK() const { return k_; }

    int getS() const { return s_; }

    int getL() const { return l_; }

   private:
    int k_, s_, l_;
    int flankMaskBp_;   // Hard mask first/last N bp at genome ends
    bool hpc_;          // Homopolymer-compressed seeds
    bool imputeAmb_;    // Skip N mutations (impute from parent)
    bool extentGuard_;  // Guard seed deletions at genome extent boundaries

   public:
    int getFlankMaskBp() const { return flankMaskBp_; }

    bool getHpc() const { return hpc_; }

    bool getImputeAmb() const { return imputeAmb_; }

    bool getExtentGuard() const { return extentGuard_; }

    uint64_t computeSubtreeSize(panmanUtils::Node* node);

    // Sequential DFS helper
    void buildIndexHelper(panmanUtils::Node* node,
                          std::unordered_set<std::string_view>& emptyNodes,
                          panmapUtils::BlockSequences& blockSequences,
                          std::vector<char>& blockExistsDelayed,
                          std::vector<char>& blockStrandDelayed,
                          panmapUtils::GlobalCoords& globalCoords,
                          std::map<uint64_t, uint64_t>& gapMap,
                          std::unordered_set<uint64_t>& invertedBlocks,
                          uint64_t& dfsIndex);

    // Process a single node and compute its seed counts
    void processNode(
        panmanUtils::Node* node,
        BuildState& state,
        panmapUtils::GlobalCoords& globalCoords,
        std::unordered_set<std::string_view>& localEmptyNodes,
        uint64_t dfsIndex,
        BacktrackInfo* backtrackInfo = nullptr,
        bool skipNodeChanges = false  // If true, don't compute/record nodeChanges (for path walking in parallel)
    );

    // Backtrack changes made by processNode using recorded changes
    void backtrackNode(BuildState& state, const BacktrackInfo& backtrackInfo);

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
        uint64_t lastNonGapScalar = UINT64_MAX);

    std::vector<std::pair<SyncmerSet::iterator, SyncmerSet::iterator>> computeNewKminmerRanges(
        std::vector<std::tuple<uint64_t, panmapUtils::seedChangeType, seeding::rsyncmer_t>>& refOnSyncmersChangeRecord,
        BuildState& state,
        const uint64_t dfsIndex);

    // Overload for sequential path (uses member state)
    std::vector<std::pair<SyncmerSet::iterator, SyncmerSet::iterator>> computeNewKminmerRanges(
        std::vector<std::tuple<uint64_t, panmapUtils::seedChangeType, seeding::rsyncmer_t>>& refOnSyncmersChangeRecord,
        const uint64_t dfsIndex);
};

}  // namespace index_single_mode
