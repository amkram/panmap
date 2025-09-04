@0xaeb74a09d04c9118;

struct CoordDelta {
  pos @0 :UInt32;
  endPos :union {
    value @1 :UInt32;
    none @2 :Void;
  }
}

struct SeedInfo {
  hash @0 :UInt64;
  startPos @1 :UInt32;
  endPos @2 :UInt32;
  isReverse @3 :Bool;
}


struct NodeChanges {
  nodeIndex @0 :UInt32;
  seedInsubIndices @1 :List(UInt32);
  seedDeletions @2 :List(UInt32);
  coordDeltas @3 :List(CoordDelta);
  invertedBlocks @4 :List(UInt32);
}

struct LiteNode {
  id @0: Text;
  parentIndex @1: UInt32;
}

struct BlockRange {
  rangeBeg @0: UInt32;
  rangeEnd @1: UInt32;
}

struct LiteTree {
  liteNodes @0: List(LiteNode);
  blockRanges @1: List(BlockRange);
}

struct MGSRIndex {
  k @0 :UInt16;
  s @1 :UInt16;
  t @2 :UInt16;
  l @3 :UInt16;
  open @4 :Bool;
  useRawSeeds @7 :Bool;

  liteTree @5 :LiteTree;
  seedInfo @6 :List(SeedInfo);
  perNodeChanges @7 :List(NodeChanges);
  
}