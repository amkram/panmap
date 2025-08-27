#pragma once

#include "panmanUtils.hpp"
#include "capnp/message.h"
#include "capnp/serialize-packed.h"
#include "mgsr_index.capnp.h"
#include "panmap_utils.hpp"
#include "seeding.hpp"
#include <eigen3/Eigen/Dense>

namespace mgsr {


typedef uint64_t minichain_t;

void updateGapMapStep(
    std::map<uint64_t, uint64_t>& gapMap,
    const std::pair<bool, std::pair<uint64_t, uint64_t>>& update,
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


struct IteratorComparator {
  bool operator()(const std::map<uint64_t, uint64_t>::iterator& lhs,
                  const std::map<uint64_t, uint64_t>::iterator& rhs) const {
      return lhs->first < rhs->first;
  }
};

bool inline compareMinichainByBeg(const minichain_t& a, const minichain_t& b) {
  return ((a >> 1) & 0x7FFFFFFF) < ((b >> 1) & 0x7FFFFFFF);
}
bool inline compareMinichainByEnd(const minichain_t& a, const minichain_t& b) {
  return ((a >> 32) & 0x7FFFFFFF) < ((b >> 32) & 0x7FFFFFFF);
}

struct VectorHash {
  size_t operator()(const std::vector<uint32_t>& v) const {
    size_t hash = 0;
    for (const auto& val : v) {
      hash ^= std::hash<uint32_t>{}(val) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
    }
    return hash;
  }
};

struct PairHash {
  std::size_t operator()(const std::pair<uint64_t, uint64_t>& p) const {
    return std::hash<uint64_t>{}(p.first) ^ (std::hash<uint64_t>{}(p.second) << 1);
  }
};

struct readSeedmer {
  const size_t hash;
  const uint32_t begPos;
  const uint32_t endPos;
  const bool rev;
  const uint32_t iorder;
};

struct SeedmerState {
  bool match;
  bool rev;
  bool inChain;
};

struct hashCoordInfoCache {
  int32_t rGlobalBeg = -1;
  int32_t rGlobalEnd = -1;
  int32_t rLocalBeg = -1;
  int32_t rLocalEnd = -1;
};

enum RefSeedmerExistStatus : uint8_t {
  EXIST_UNIQUE,
  EXIST_DUPLICATE,
  NOT_EXIST
};

enum RefSeedmerChangeType : uint8_t {
  EXIST_UNIQUE_TO_EXIST_UNIQUE = 0,
  EXIST_UNIQUE_TO_EXIST_DUPLICATE = 1,
  EXIST_UNIQUE_TO_NOT_EXIST = 2,
  EXIST_DUPLICATE_TO_EXIST_UNIQUE = 3,
  EXIST_DUPLICATE_TO_EXIST_DUPLICATE = 4,
  EXIST_DUPLICATE_TO_NOT_EXIST = 5,
  NOT_EXIST_TO_EXIST_UNIQUE = 6,
  NOT_EXIST_TO_EXIST_DUPLICATE = 7,
  NOT_EXIST_TO_NOT_EXIST = 8
};

struct RefSeedmerChangeCountStats {
  size_t EXIST_UNIQUE_TO_EXIST_UNIQUE = 0;
  size_t EXIST_UNIQUE_TO_EXIST_DUPLICATE = 0;
  size_t EXIST_UNIQUE_TO_NOT_EXIST = 0;
  size_t EXIST_DUPLICATE_TO_EXIST_UNIQUE = 0;
  size_t EXIST_DUPLICATE_TO_EXIST_DUPLICATE = 0;
  size_t EXIST_DUPLICATE_TO_NOT_EXIST = 0;
  size_t NOT_EXIST_TO_EXIST_UNIQUE = 0;
  size_t NOT_EXIST_TO_EXIST_DUPLICATE = 0;
  size_t NOT_EXIST_TO_NOT_EXIST = 0;
  size_t TOTAL_SEEDMERS = 0;

  bool operator==(const RefSeedmerChangeCountStats& other) const noexcept {
    return
      EXIST_UNIQUE_TO_EXIST_UNIQUE == other.EXIST_UNIQUE_TO_EXIST_UNIQUE &&
      EXIST_UNIQUE_TO_EXIST_DUPLICATE == other.EXIST_UNIQUE_TO_EXIST_DUPLICATE &&
      EXIST_UNIQUE_TO_NOT_EXIST == other.EXIST_UNIQUE_TO_NOT_EXIST &&
      EXIST_DUPLICATE_TO_EXIST_UNIQUE == other.EXIST_DUPLICATE_TO_EXIST_UNIQUE &&
      EXIST_DUPLICATE_TO_EXIST_DUPLICATE == other.EXIST_DUPLICATE_TO_EXIST_DUPLICATE &&
      EXIST_DUPLICATE_TO_NOT_EXIST == other.EXIST_DUPLICATE_TO_NOT_EXIST &&
      NOT_EXIST_TO_EXIST_UNIQUE == other.NOT_EXIST_TO_EXIST_UNIQUE &&
      NOT_EXIST_TO_EXIST_DUPLICATE == other.NOT_EXIST_TO_EXIST_DUPLICATE &&
      NOT_EXIST_TO_NOT_EXIST == other.NOT_EXIST_TO_NOT_EXIST &&
      TOTAL_SEEDMERS == other.TOTAL_SEEDMERS;
  }
};

struct RefSeedmerChangeCountStatsHash {
  std::size_t operator()(const RefSeedmerChangeCountStats& stats) const noexcept {
    return static_cast<std::size_t>(stats.EXIST_UNIQUE_TO_EXIST_UNIQUE) * 100000000000000000ULL +
           static_cast<std::size_t>(stats.EXIST_UNIQUE_TO_EXIST_DUPLICATE) * 1000000000000000ULL +
           static_cast<std::size_t>(stats.EXIST_UNIQUE_TO_NOT_EXIST) * 10000000000000ULL +
           static_cast<std::size_t>(stats.EXIST_DUPLICATE_TO_EXIST_UNIQUE) * 100000000000ULL +
           static_cast<std::size_t>(stats.EXIST_DUPLICATE_TO_EXIST_DUPLICATE) * 1000000000ULL +
           static_cast<std::size_t>(stats.EXIST_DUPLICATE_TO_NOT_EXIST) * 10000000ULL +
           static_cast<std::size_t>(stats.NOT_EXIST_TO_EXIST_UNIQUE) * 100000ULL +
           static_cast<std::size_t>(stats.NOT_EXIST_TO_EXIST_DUPLICATE) * 1000ULL +
           static_cast<std::size_t>(stats.NOT_EXIST_TO_NOT_EXIST) * 10ULL +
           static_cast<std::size_t>(stats.TOTAL_SEEDMERS);
  }
};

enum readType : uint8_t {
  PASS,
  HIGH_DUPLICATES,
  IDENTICAL_SCORE_ACROSS_NODES
};

class Read {
  public:
    std::vector<readSeedmer> seedmersList;
    std::vector<SeedmerState> seedmerStates;
    std::unordered_map<size_t, std::vector<uint64_t>> uniqueSeedmers;
    std::vector<minichain_t> minichains;
    std::unordered_set<int32_t> duplicates;
};

void seedmersFromFastq(
  const std::string& readPath1, const std::string& readPath2,
  std::vector<Read>& reads,
  std::unordered_map<size_t, std::vector<std::pair<uint32_t, std::vector<uint64_t>>>>& seedmerToReads,
  std::vector<std::vector<size_t>>& readSeedmersDuplicatesIndex,
  int k, int s, int t, int l, bool open, bool fast_mode
);


class mgsrIndexBuilder {
  public:
    // capnp object
    ::capnp::MallocMessageBuilder outMessage;
    MGSRIndex::Builder indexBuilder;
    capnp::List<NodeChanges>::Builder perNodeChanges;
    

    // tree pointer
    panmanUtils::Tree *T;
  
    // syncmer and k-min-mer objects
    std::unordered_map<uint32_t, std::unordered_set<uint64_t>> blockOnSyncmers;
    std::vector<std::optional<seeding::rsyncmer_t>> refOnSyncmers;
    std::set<uint64_t> refOnSyncmersMap;

    std::vector<seeding::uniqueKminmer_t> uniqueKminmers;
    std::vector<std::optional<uint64_t>> refOnKminmers;
    std::unordered_map<seeding::uniqueKminmer_t, uint64_t> kminmerToUniqueIndex;

    mgsrIndexBuilder(panmanUtils::Tree *T, int k, int s, int t, int l, bool open) 
      : outMessage(), indexBuilder(outMessage.initRoot<MGSRIndex>()), T(T)
    {
      indexBuilder.setK(k);
      indexBuilder.setS(s);
      indexBuilder.setT(t);
      indexBuilder.setL(l);
      indexBuilder.setOpen(open);
      perNodeChanges = indexBuilder.initPerNodeChanges(T->allNodes.size());
    }

    void buildIndex();

    void writeIndex(const std::string& path);
  private:
    void buildIndexHelper(
      panmanUtils::Node *node,
      panmapUtils::BlockSequences &blockSequences,
      std::vector<char> &blockExistsDelayed,
      std::vector<char> &blockStrandDelayed,
      panmapUtils::GlobalCoords &globalCoords,
      std::map<uint64_t, uint64_t> &gapMap,
      std::unordered_set<uint64_t> &invertedBlocks,
      uint64_t &dfsIndex
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
      std::vector<std::tuple<uint64_t, uint64_t, panmapUtils::seedChangeType>>& blockOnSyncmersBacktracks
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
      std::vector<std::tuple<uint64_t, uint64_t, panmapUtils::seedChangeType>>& blockOnSyncmersBacktracks
    );

    std::vector<std::pair<std::set<uint64_t>::iterator, std::set<uint64_t>::iterator>> computeNewKminmerRanges(
      std::vector<std::tuple<uint64_t, panmapUtils::seedChangeType, seeding::rsyncmer_t>>& refOnSyncmersChangeRecord,
      const uint64_t dfsIndex
    );
};

struct readScoreDelta {
  size_t readIndex;
  uint32_t scoreDelta;
};

class mgsrPlacer {
  public:
    // mutation structures
    std::vector<seeding::uniqueKminmer_t> seedInfos;
    std::vector<std::vector<uint32_t>> seedInsubIndices;
    std::vector<std::vector<uint32_t>> seedDeletions;
    std::vector<std::vector<std::pair<uint32_t, std::optional<uint32_t>>>> coordDeltas;
    std::vector<std::vector<uint32_t>> invertedBlocks;

    // tree pointer
    panmanUtils::Tree *T;
    std::unordered_map<std::string, int64_t> nodeToDfsIndex;

    // parameters from index
    uint k = 0;
    uint s = 0;
    uint t = 0;
    uint l = 0;
    bool openSyncmer = false;

    // parameters from user input... preset for now
    double excludeDuplicatesThreshold = 0.5;
    double errorRate = 0.005;
    int64_t maximumGap = 100;
    
    // dynamic reference kminmer structures
    std::map<uint64_t, uint64_t> gapMap;
    std::map<uint64_t, uint64_t> positionMap;
    std::unordered_map<size_t, std::vector<std::map<uint64_t, uint64_t>::iterator>> hashToPositionMap;
    std::unordered_map<size_t, RefSeedmerExistStatus> delayedRefSeedmerStatus;

    // preallocated structures to prevent memory allocation and deletion in tight loops
    std::vector<std::pair<minichain_t, bool>> minichainsToUpdate;
    std::vector<std::pair<std::vector<uint64_t>, uint64_t>> readToAffectedSeedmerStorageHelper;  // might be memory intensive....

    // current query kminmer structures
    std::vector<mgsr::Read> reads;
    std::vector<mgsr::readType> readTypes;
    std::unordered_map<size_t, std::vector<std::pair<uint32_t, std::vector<uint64_t>>>> seedmerToReads;
    std::vector<std::vector<size_t>> readSeedmersDuplicatesIndex;

    // current query score index structures
    std::vector<int32_t> readScores;
    std::vector<std::pair<int64_t, int64_t>> kminmerMatches;
    std::vector<std::vector<readScoreDelta>> perNodeScoreDeltasIndex;
    std::vector<std::vector<std::tuple<size_t, int64_t, int64_t>>> perNodeKminmerMatchesDeltasIndex;
    std::vector<std::vector<minichain_t>> maxMinichains;
    std::vector<int32_t> maxScores;
    int64_t totalScore = 0;
    int64_t totalDirectionalKminmerMatches = 0;
    
    // counters for calculating overlap coefficient
    std::unordered_map<std::string, double> kminmerOverlapCoefficients;
    size_t binaryOverlapKminmerCount = 0;

    // for identical pairs of nodes
    std::unordered_map<std::string, uint64_t> totalScores;
    std::unordered_map<std::string, std::vector<std::string>> identicalGroups;
    std::unordered_map<std::string, std::string> identicalNodeToGroup;

    // misc
    uint64_t curDfsIndex = 0;
    uint64_t readMinichainsInitialized = 0;

    uint64_t readMinichainsAdded = 0;
    uint64_t readMinichainsAddedToEmpty = 0;
    uint64_t readMinichainsAddedToSingleton = 0;
    uint64_t readMinichainsAddedToMultiple = 0;

    uint64_t readMinichainsRemoved = 0;
    uint64_t readMinichainsRemovedInplace = 0;
    uint64_t readMinichainsRemovedFromMultiple = 0;
    
    mgsrPlacer(panmanUtils::Tree* tree, const std::string& path) : T(tree) {
      ::capnp::ReaderOptions readerOptions {
        .traversalLimitInWords = std::numeric_limits<uint64_t>::max(),
        .nestingLimit = 1024
      };
      int fd = open_file(path);
      ::capnp::PackedFdMessageReader reader(fd, readerOptions);
      MGSRIndex::Reader indexReader = reader.getRoot<MGSRIndex>();
      capnp::List<SeedInfo>::Reader seedInfosReader = indexReader.getSeedInfo();
      capnp::List<NodeChanges>::Reader perNodeChangesReader = indexReader.getPerNodeChanges();
      
      k = indexReader.getK();
      s = indexReader.getS();
      t = indexReader.getT();
      l = indexReader.getL();
      openSyncmer = indexReader.getOpen();

      seedInfos.resize(seedInfosReader.size());
      for (size_t i = 0; i < seedInfos.size(); i++) {
        const auto& seedReader = seedInfosReader[i];
        auto& seed = seedInfos[i];
        seed.hash = seedReader.getHash();
        seed.startPos = seedReader.getStartPos();
        seed.endPos = seedReader.getEndPos();
        seed.isReverse = seedReader.getIsReverse();
      }

      seedInsubIndices.resize(perNodeChangesReader.size());
      seedDeletions.resize(perNodeChangesReader.size());
      coordDeltas.resize(perNodeChangesReader.size());
      invertedBlocks.resize(perNodeChangesReader.size());
      for (size_t i = 0; i < perNodeChangesReader.size(); i++) {
        const auto& currentPerNodeChangeReader = perNodeChangesReader[i];
        auto& currentSeedInsubIndices = seedInsubIndices[i];
        auto& currentSeedDeletions = seedDeletions[i];
        auto& currentCoordDeltas = coordDeltas[i];
        auto& currentInvertedBlocks = invertedBlocks[i];

        const auto& currentSeedInsubIndicesReader = currentPerNodeChangeReader.getSeedInsubIndices();
        const auto& currentSeedDeletionsReader = currentPerNodeChangeReader.getSeedDeletions();
        const auto& currentCoordDeltasReader = currentPerNodeChangeReader.getCoordDeltas();
        const auto& currentInvertedBlocksReader = currentPerNodeChangeReader.getInvertedBlocks();

        currentSeedInsubIndices.resize(currentSeedInsubIndicesReader.size());
        currentSeedDeletions.resize(currentSeedDeletionsReader.size());
        currentCoordDeltas.resize(currentCoordDeltasReader.size());
        currentInvertedBlocks.resize(currentInvertedBlocksReader.size());

        for (size_t j = 0; j < currentSeedInsubIndicesReader.size(); j++) {
          currentSeedInsubIndices[j] = currentSeedInsubIndicesReader[j];
        }

        for (size_t j = 0; j < currentSeedDeletionsReader.size(); j++) {
          currentSeedDeletions[j] = currentSeedDeletionsReader[j];
        }

        for (size_t j = 0; j < currentCoordDeltasReader.size(); j++) {
          const auto& currentCoordDeltaReader = currentCoordDeltasReader[j];
          auto& currentCoordDelta = currentCoordDeltas[j];
          currentCoordDelta.first = currentCoordDeltaReader.getPos();
          currentCoordDelta.second = currentCoordDeltaReader.getEndPos().which() == CoordDelta::EndPos::VALUE ? std::optional<uint32_t>(currentCoordDeltaReader.getEndPos().getValue()) : std::nullopt;
        }

        for (size_t j = 0; j < currentInvertedBlocksReader.size(); j++) {
          currentInvertedBlocks[j] = currentInvertedBlocksReader[j];
        }
      }
    }

    void initializeQueryData(const std::string& readPath1, const std::string& readPath2, bool fast_mode = false);

    void placeReadsHelper(panmanUtils::Node* node, const panmapUtils::GlobalCoords& globalCoords);
    void placeReads();

    // for updating reference seeds and gapMap
    void updateSeeds(std::vector<std::pair<uint64_t, panmapUtils::seedChangeType>>& seedBacktracks, std::unordered_set<uint64_t>& affectedSeedmers);
    void updateGapMap(const panmapUtils::GlobalCoords& globalCoords, std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& gapMapBacktracks, std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& gapMapBlocksBacktracks);
    void addSeedAtPosition(uint64_t uniqueKminmerIndex, std::vector<std::pair<uint64_t, panmapUtils::seedChangeType>>& seedBacktracks, std::unordered_set<uint64_t>& affectedSeedmers);
    void addSeedAtPosition(uint64_t uniqueKminmerIndex);
    void delSeedAtPosition(uint64_t pos, std::vector<std::pair<uint64_t, panmapUtils::seedChangeType>>& seedBacktracks, std::unordered_set<uint64_t>& affectedSeedmers);
    void delSeedAtPosition(uint64_t pos);

    // for updating read scores and kminmer matches
    void setReadScore(size_t readIndex, const int32_t score, const size_t numDuplicates);
    void initializeReadMinichains(size_t readIndex);
    void initializeReadMinichains(mgsr::Read& curRead);
    void updateMinichains(size_t readIndex, const std::vector<uint64_t>& affectedSeedmerIndexCodes, const uint64_t affectedSeedmerIndexCodesSize, bool allUniqueToNonUnique, bool allNonUniqueToUnique);
    int64_t getReadPseudoScore(mgsr::Read& curRead, const std::map<uint64_t, uint64_t>& degapCoordIndex, const std::map<uint64_t, uint64_t>& regapCoordIndex, std::unordered_map<size_t, mgsr::hashCoordInfoCache>& hashCoordInfoCacheTable);
    inline uint64_t decodeBegFromMinichain(uint64_t minichain);
    inline uint64_t decodeEndFromMinichain(uint64_t minichain);
    inline bool decodeRevFromMinichain(uint64_t minichain);
    RefSeedmerExistStatus getCurrentRefSeedmerExistStatus(uint64_t hash);
    RefSeedmerExistStatus getDelayedRefSeedmerExistStatus(uint64_t hash);

    // For preparing for squareEM
    void updateIdenticalGroups();
    void mergeIdenticalNodes(const std::unordered_set<std::string>& identicalGroup);
    bool identicalReadScores(const std::string& node1, const std::string& node2, bool fast_mode = false) const;
    std::vector<uint32_t> getScoresAtNode(const std::string& nodeId) const;
    void getScoresAtNode(const std::string& nodeId, std::vector<uint32_t>& curNodeScores) const;
    




  private:
    int open_file(const std::string& path);

    void updateRefSeedmerStatus(size_t hash, mgsr::RefSeedmerChangeType& seedmerChangeType, mgsr::RefSeedmerExistStatus refSeedmerOldStatus, mgsr::RefSeedmerExistStatus refSeedmerNewStatus);
    void updateSeedmerChangesTypeFlag(mgsr::RefSeedmerChangeType seedmerChangeType, std::pair<bool, bool>& flags);
    void fillReadToAffectedSeedmerIndex(
      std::unordered_map<size_t, std::pair<uint64_t, std::pair<bool, bool>>>& readToAffectedSeedmerIndex,
      const std::unordered_set<uint64_t>& affectedSeedmers
    );
    uint64_t extendMinichain(std::map<uint64_t, uint64_t>::const_iterator refPositionIt, const mgsr::Read& curRead, uint64_t& curEnd, bool rev, uint64_t qidx, uint64_t c);
    void extendChainRemoval(uint64_t& c, uint64_t& curEnd, const std::vector<uint64_t>& affectedSeedmerIndexCode, const uint64_t affectedSeedmerIndexCodesSize, uint32_t lastSeedmerIndex);
    void extendChainAddition(uint64_t& c, uint64_t& curEnd, const std::vector<uint64_t>& affectedSeedmerIndexCode, const uint64_t affectedSeedmerIndexCodesSize, bool chainRev, std::map<uint64_t, uint64_t>::const_iterator refPositionIt, uint64_t readIndex);
    bool colinearAdjacent(std::map<uint64_t, uint64_t>::const_iterator from, std::map<uint64_t, uint64_t>::const_iterator to, bool fromRev, bool toRev);
    bool colinearAdjacent(std::map<uint64_t, uint64_t>::const_iterator from, std::map<uint64_t, uint64_t>::const_iterator to, bool fromRev);
    void addToMinichains(const std::vector<readSeedmer>& curSeedmerList, std::vector<minichain_t>& curMinichains, uint64_t minichain);
    void removeFromMinichains(std::vector<minichain_t>& curMinichains, uint64_t minichain);
    uint64_t getRefSeedmerBegFromHash(const size_t hash) const;
    uint64_t getRefSeedmerEndFromHash(const size_t hash) const;
    bool isColinearFromMinichains(
      mgsr::Read& curRead, const bool rev,
      const uint64_t beg1, const uint64_t end1,
      const uint64_t beg2, const uint64_t end2,
      const std::map<uint64_t, uint64_t>& degapCoordIndex,
      const std::map<uint64_t, uint64_t>& regapCoordIndex,
      std::unordered_map<size_t, mgsr::hashCoordInfoCache>& hashCoordInfoCacheTable) const;

    int64_t getReadBruteForceScore(size_t readIndex, const std::map<uint64_t, uint64_t>& degapCoordIndex, const std::map<uint64_t, uint64_t>& regapCoordIndex, std::unordered_map<size_t, mgsr::hashCoordInfoCache>& hashCoordInfoCacheTable);
};

class squareEM {
  public:
    // main prob and prop matrices
    Eigen::MatrixXd probs;
    Eigen::VectorXd props;

    // preallocate intermediate variables
    Eigen::VectorXd props0;
    Eigen::VectorXd props1;
    Eigen::VectorXd props2;
    Eigen::VectorXd propsSq;
    Eigen::VectorXd denoms;
    Eigen::VectorXd inverseDenoms;
    Eigen::VectorXd r;
    Eigen::VectorXd v;
    Eigen::VectorXd readDuplicates;
    double llh;
    double invTotalWeight;

    std::unordered_map<std::string, std::vector<std::string>> identicalGroups;
    std::unordered_map<std::string, std::string> identicalNodeToGroup;
    std::vector<std::string> nodes;
    size_t numNodes;
    size_t curIteration = 0;


    // Input parameters... hardcoded for now
    double eta = 0.00001;
    double maxChangeThreshold = 0.0001;
    double propThresholdToRemove = 0.005;

    squareEM(mgsrPlacer& placer, uint32_t overlapCoefficientCutoff);


    void runSquareEM(uint64_t maximumIterations);
    bool removeLowPropNodes();

  private:
    void updateProps1();
    void updateProps2();
    void updateProps(const Eigen::VectorXd& original, Eigen::VectorXd& result);
    void normalizeProps(Eigen::VectorXd& props);
    double getExp(const Eigen::VectorXd& props);
    void resetIteration();
};




// end of namespace mgsr
}