#pragma once

#include "panmanUtils.hpp"
#include "capnp/message.h"
#include "capnp/serialize-packed.h"
#include "mgsr_index.capnp.h"
#include "panmap_utils.hpp"
#include "seeding.hpp"


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
      std::vector<std::tuple<uint64_t, panmapUtils::seedChangeType, seeding::rsyncmer_t>>& refOnSyncmersChangeRecord
    );
};

struct readScoreDelta {
  size_t readIndex;
  uint32_t scoreDelta;
  double probDelta;
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
    std::vector<std::pair<int32_t, double>> readScores;
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
    std::unordered_map<std::string, std::vector<std::string>> identicalGroups;
    std::unordered_map<std::string, std::string> identicalNodeToGroup;

    // misc
    uint64_t curDfsIndex = 0;
    
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
    void setReadScore(size_t readIndex, const int32_t score, const double prob, const size_t numDuplicates);
    void initializeReadMinichains(size_t readIndex);
    void initializeReadMinichains(mgsr::Read& curRead);
    void updateMinichains(size_t readIndex, const std::vector<uint64_t>& affectedSeedmerIndexCodes, const uint64_t affectedSeedmerIndexCodesSize, bool allUniqueToNonUnique, bool allNonUniqueToUnique);
    int64_t getReadPseudoScore(mgsr::Read& curRead, const std::map<uint64_t, uint64_t>& degapCoordIndex, const std::map<uint64_t, uint64_t>& regapCoordIndex, std::unordered_map<size_t, mgsr::hashCoordInfoCache>& hashCoordInfoCacheTable);
    inline uint64_t decodeBegFromMinichain(uint64_t minichain);
    inline uint64_t decodeEndFromMinichain(uint64_t minichain);
    inline bool decodeRevFromMinichain(uint64_t minichain);

    RefSeedmerExistStatus getCurrentRefSeedmerExistStatus(uint64_t hash);
    RefSeedmerExistStatus getDelayedRefSeedmerExistStatus(uint64_t hash);




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




// end of namespace mgsr
}