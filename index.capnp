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
struct SeedPosition {
  pos @0 :Int64;
  dfsIndex @1 :Int64;
}

struct HotSeedEntry {
  kmer @0 :UInt64;
  endPos @1 :Int64;
  isReverse @2 :Bool;
  accessCount @3 :Int64;
  positions @4 :List(SeedPosition);
}

struct HotSeedIndexSerial {
  entries @0 :List(HotSeedEntry);
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
