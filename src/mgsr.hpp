#pragma once

#include "panmanUtils.hpp"
#include "capnp/message.h"
#include "capnp/serialize-packed.h"
#include "mgsr_index.capnp.h"
#include "panmap_utils.hpp"
#include "seeding.hpp"


namespace mgsr {


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

class Read {
  public:
    std::vector<readSeedmer> seedmersList;
    std::vector<SeedmerState> seedmerStates;
    std::unordered_map<size_t, std::vector<uint64_t>> uniqueSeedmers;
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
      std::vector<bool> &blockExistsDelayed,
      std::vector<bool> &blockStrandDelayed,
      panmapUtils::GlobalCoords &globalCoords,
      std::map<uint64_t, uint64_t> &gapMap,
      std::unordered_set<uint64_t> &invertedBlocks,
      uint64_t &dfsIndex
    );

    std::vector<panmapUtils::NewSyncmerRange> computeNewSyncmerRanges(
      panmanUtils::Node* node,
      size_t dfsIndex,
      const panmapUtils::BlockSequences& blockSequences,
      const panmapUtils::GlobalCoords& globalCoords,
      std::vector<std::pair<panmapUtils::Coordinate, panmapUtils::Coordinate>>& localMutationRanges,
      std::vector<std::tuple<uint64_t, uint64_t, panmapUtils::seedChangeType>>& blockOnSyncmersBacktracks
    );

    std::vector<std::pair<std::set<uint64_t>::iterator, std::set<uint64_t>::iterator>> computeNewKminmerRanges(
      std::vector<std::tuple<uint64_t, panmapUtils::seedChangeType, seeding::rsyncmer_t>>& refOnSyncmersChangeRecord
    );
};

class mgsrIndexReader {
  public:
    ::capnp::ReaderOptions readerOptions {
      .traversalLimitInWords = std::numeric_limits<uint64_t>::max(),
      .nestingLimit = 1024
    };
    ::capnp::PackedFdMessageReader reader;
    MGSRIndex::Reader indexReader;
    capnp::List<SeedInfo>::Reader seedInfos;
    capnp::List<NodeChanges>::Reader perNodeChanges;

    panmanUtils::Tree *T;

    mgsrIndexReader(panmanUtils::Tree* tree, const std::string& path)
        : mgsrIndexReader(tree, open_file(path)) {}

  private:

    mgsrIndexReader(panmanUtils::Tree* tree, int fd)
    : reader(fd, readerOptions),
      indexReader(reader.getRoot<MGSRIndex>()),
      seedInfos(indexReader.getSeedInfo()),
      perNodeChanges(indexReader.getPerNodeChanges()),
      T(tree)
    {}

    int open_file(const std::string& path);
};



class mgsrPlacer {
  public:
    mgsrIndexReader indexReader;
    
    panmanUtils::Tree *T;
    
    // dynamic objects
    std::map<uint64_t, uint64_t> gapMap;
    std::map<uint64_t, uint64_t> positionMap;
    std::unordered_map<size_t, std::set<std::map<uint64_t, uint64_t>::iterator, IteratorComparator>> hashToPositionMap;

    mgsrPlacer(panmanUtils::Tree* tree, const std::string& path)
      : indexReader(tree, path),
        T(tree) {
    }

    void initialize();

    void placeReadsHelper(panmanUtils::Node* node, uint64_t& dfsIndex, const panmapUtils::GlobalCoords& globalCoords);
    void placeReads();

    void updateSeeds(uint64_t currentDfsIndex, std::vector<std::pair<uint64_t, panmapUtils::seedChangeType>>& seedBacktracks);
    void updateGapMap(uint64_t currentDfsIndex, const panmapUtils::GlobalCoords& globalCoords, std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& gapMapBacktracks, std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& gapMapBlocksBacktracks);
    void addSeedAtPosition(uint64_t uniqueKminmerIndex, std::vector<std::pair<uint64_t, panmapUtils::seedChangeType>>& seedBacktracks);
    void addSeedAtPosition(uint64_t uniqueKminmerIndex);
    void delSeedAtPosition(uint64_t pos, std::vector<std::pair<uint64_t, panmapUtils::seedChangeType>>& seedBacktracks);
    void delSeedAtPosition(uint64_t pos);




  private:
};



// end of namespace mgsr
}