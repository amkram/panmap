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


struct NodeChanges {
  # deprecate soon
  nodeIndex @0 :UInt32;
  seedDeltaIndices @1 :List(UInt32);
  seedDeltaIsDeleted @2 :List(Bool);
  gapRunDeltas @3 :List(GapRunDelta);
  invertedBlocks @4 :List(UInt32);
}

struct LiteIndex {
  k @0 :UInt16;
  s @1 :UInt16;
  t @2 :UInt16;
  l @3 :UInt16;
  open @4 :Bool;
  liteTree @5 :LiteTree;

  # Struct-of-arrays format for seed changes
  seedChangeHashes @6 :List(List(UInt64));
  seedChangeParentCounts @7 :List(List(Int16));
  seedChangeChildCounts @8 :List(List(Int16));
  nodeChangeOffsets @9 :List(UInt32);  # Size = numNodes + 1


  # Homopolymer compression mode: if true, seeds were extracted from HPC sequences
  hpc @10 :Bool = false;

  # mgsr (deprecate soon)
  seedHashes @11 :List(UInt64);
  seedStartPos @12 :List(UInt32);
  seedEndPos @13 :List(UInt32);
  seedIsReverse @14 :List(Bool);

  perNodeChanges @15 :List(NodeChanges);



  # 4x4 nucleotide substitution rate matrix (A=0,C=1,G=2,T=3), row-major
  substitutionMatrix @16 :List(Float64);
}
