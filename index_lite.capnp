@0xbeb74a09d04c9119;

struct SeedInfo {
  hash @0 :UInt64;  # Only hash stored; positions not needed for lite placement
}

struct SeedChange {
  seedHash @0 :UInt64;      # Hash of the seed that changed
  parentCount @1 :Int64;    # Absolute count in parent node
  delta @2 :Int64;          # Change from parent to child (childCount = parentCount + delta)
}

struct NodeSeedChanges {
  nodeIndex @0 :UInt32;
  seedChanges @1 :List(SeedChange);  # Only seeds that differ from parent
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

struct LiteIndex {
  k @0 :UInt16;
  s @1 :UInt16;
  t @2 :UInt16;
  l @3 :UInt16;
  open @4 :Bool;

  liteTree @5 :LiteTree;
  seedInfo @6 :List(SeedInfo);
  deprecatedPerNodeChanges @7 :List(NodeSeedChanges); # DEPRECATED: Was perNodeChanges
  useRawSeeds @8 :Bool;
  version @9 :UInt16 = 3;  # Version 3 with struct-of-arrays for parallel access
  
  largestNodeChangeCount @10 :UInt32;  # Hint for pre-allocation
  totalSeedChanges @11 :UInt64;        # Total across all nodes
  
  # Pre-computed genome-only metrics (parallel arrays indexed by DFS index)
  genomeMagnitudeSquared @12 :List(Float64);  # Sum of (count^2) for all seeds in each genome
  genomeUniqueSeedCount @13 :List(UInt32);    # Total number of unique seeds in each genome
  genomeTotalSeedFrequency @14 :List(Int64);  # Sum of all seed counts in each genome

  # VERSION 3: STRUCT-OF-ARRAYS format for optimal parallel deserialization
  # All three arrays have the same length (totalSeedChanges).
  # This layout allows threads to read contiguous primitive arrays without pointer chasing.
  # Cap'n Proto stores primitive lists as contiguous memory (encoding C=5 for 64-bit values).
  seedChangeHashes @15 :List(UInt64);        # Seed hashes for all changes (flat array)
  
  # Offsets into seed change arrays for each node.
  # nodeChangeOffsets[i] is the start index for node i.
  # The number of changes for node i is nodeChangeOffsets[i+1] - nodeChangeOffsets[i].
  # Size should be numNodes + 1.
  nodeChangeOffsets @16 :List(UInt32);
  
  seedChangeParentCounts @17 :List(Int64);   # Parent counts for all changes (flat array)
  seedChangeChildCounts @18 :List(Int64);    # Child counts (parentCount + delta) for all changes (flat array)
  
  deprecatedAllSeedChanges @19 :List(SeedChange); # DEPRECATED: Was allSeedChanges
}