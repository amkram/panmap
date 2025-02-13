@0xf6f6f6f6f6f6f6f6;

struct MapDelta {
    pos @0 :Int32;
    maybeValue :union {
        value @1 :Int64;
        none @2 :Void;
    }
}

struct GapMutations {
    deltas @0 :List(MapDelta);
}


struct SeedMutations {
    basePositions @0 :List(Int64);
    perPosMasks @1 :List(UInt64);
}


struct DFSPositionKmer {
  dfsIndex @0 :UInt32;
  position @1 :UInt64;
  kmer @2 :UInt64;
  endPosition @3 :UInt64;
  isReverse @4 :Bool;
}
struct SeedWithEndPosition {
  seed @0 :UInt64;
  endPosition @1 :UInt64;
  isReverse @2 :Bool;
}
struct HotSeedIndexSerial {
  immutableSeeds @0 :List(SeedWithEndPosition); # Seeds that never change
  hotSeeds @1 :List(DFSPositionKmer); # Frequently mutating positions
}

struct Index {
    k @0 :Int32;
    s @1 :Int32;
    t @2 :Int32;
    l @3 :Int32;
    open @4 :Bool;
    perNodeSeedMutations @5 :List(SeedMutations);
    perNodeGapMutations @6 :List(GapMutations);
    hotSeedIndex @7 :HotSeedIndexSerial;
}
