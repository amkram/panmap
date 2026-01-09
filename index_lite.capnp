@0xbeb74a09d04c9119;

# Panmap Lite Index Schema

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
  useRawSeeds @6 :Bool;
  version @7 :UInt16 = 3;
  largestNodeChangeCount @8 :UInt32;
  totalSeedChanges @9 :UInt64;
  
  # Pre-computed genome metrics (parallel arrays indexed by DFS index)
  genomeMagnitudeSquared @10 :List(Float64);
  genomeUniqueSeedCount @11 :List(UInt32);
  genomeTotalSeedFrequency @12 :List(Int64);

  # Struct-of-arrays format for seed changes
  # All three arrays have length totalSeedChanges.
  seedChangeHashes @13 :List(UInt64);
  nodeChangeOffsets @14 :List(UInt32);  # Size = numNodes + 1
  seedChangeParentCounts @15 :List(Int64);
  seedChangeChildCounts @16 :List(Int64);
  
  # IDF (Inverse Document Frequency) data for seed weighting
  # Seeds that appear in fewer genomes are more informative
  totalLeafGenomes @17 :UInt32;           # Total number of leaf genomes (N for IDF)
  idfSeedHashes @18 :List(UInt64);        # Unique seed hashes (sorted)
  idfGenomeCounts @19 :List(UInt32);      # Number of genomes containing each seed (parallel array)
}
