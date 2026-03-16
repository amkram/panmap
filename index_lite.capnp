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

struct GapRunDelta {
  # deprecate soon
  startPos @0 :UInt32;
  endPos @1 :UInt32;
  toGap @2 :Bool;
}

struct SeedInfo {
  # deprecate soon
  hash @0 :UInt64;
  startPos @1 :UInt32;
  endPos @2 :UInt32;
  isReverse @3 :Bool;
}

struct SeedDelta {
  # deprecate soon
  seedIndex @0 :UInt32;
  isDeleted @1 :Bool;
}


struct NodeChanges {
  # deprecate soon
  nodeIndex @0 :UInt32;
  seedDeltas @1 :List(SeedDelta);
  gapRunDeltas @2 :List(GapRunDelta);
  invertedBlocks @3 :List(UInt32);
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
  seedChangeParentCounts @8 :List(Int16);
  seedChangeChildCounts @9 :List(Int16);

  # Homopolymer compression mode: if true, seeds were extracted from HPC sequences
  hpc @10 :Bool = false;

  # 4x4 nucleotide substitution rate matrix (A=0,C=1,G=2,T=3), row-major
  # Each entry is P(col | row) per branch, computed from the pangenome tree.
  # 16 elements: [A→A, A→C, A→G, A→T, C→A, C→C, ...]
  substitutionMatrix @16 :List(Float64);

    # mgsr (deprecate soon)
  seedInfo @11 :List(SeedInfo);
  perNodeChanges @12 :List(NodeChanges);

  # Overflow arrays for large indices (>500M seed changes)
  # When present, the full array = concat(primary, overflow)
  seedChangeHashes2 @13 :List(UInt64);
  seedChangeParentCounts2 @14 :List(Int16);
  seedChangeChildCounts2 @15 :List(Int16);
}
