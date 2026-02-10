@0xbeb74a09d04c9119;

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

struct LiteIndex {
  k @0 :UInt16;
  s @1 :UInt16;
  t @2 :UInt16;
  l @3 :UInt16;
  open @4 :Bool;
  liteTree @5 :LiteTree;

  # Struct-of-arrays format for seed changes
  seedChangeHashes @6 :List(UInt64);
  nodeChangeOffsets @7 :List(UInt32);  # Size = numNodes + 1
  seedChangeParentCounts @8 :List(Int64);
  seedChangeChildCounts @9 :List(Int64);

  # Homopolymer compression mode: if true, seeds were extracted from HPC sequences
  hpc @10 :Bool = false;
}
