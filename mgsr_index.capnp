@0xaeb74a09d04c9118;

struct GapRunDelta {
  startPos @0 :UInt32;
  endPos @1 :UInt32;
  toGap @2 :Bool;
}

struct SeedInfo {
  hash @0 :UInt64;
  startPos @1 :UInt32;
  endPos @2 :UInt32;
  isReverse @3 :Bool;
}

struct SeedDelta {
  seedIndex @0 :UInt32;
  isDeleted @1 :Bool;
}


struct NodeChanges {
  nodeIndex @0 :UInt32;
  seedDeltas @1 :List(SeedDelta);
  gapRunDeltas @2 :List(GapRunDelta);
  invertedBlocks @3 :List(UInt32);
}

struct LiteNode {
  id @0: Text;
  parentIndex @1: UInt32;
  identicalToParent @2: Bool;
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

  liteTree @5 :LiteTree;
  seedInfo @6 :List(SeedInfo);
  perNodeChanges @7 :List(NodeChanges);
  
}