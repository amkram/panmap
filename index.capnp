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

struct KmerDictionaryEntry {
    sequence @0 :Text;    # Actual k-mer sequence
    canonicalForm @1 :Bool; # Whether this is the canonical (lexicographically smaller) form
}

struct SeedMutations {
    basePositions @0 :List(Int64);
    perPosMasks @1 :List(UInt64);
    kmerDictionaryIds @2 :List(UInt32);  # IDs into the global dictionary
    kmerPositions @3 :List(Int64);       # Positions these dictionary entries correspond to
}

struct NodePathInfo {
    nodeId @0 :Text;
    dfsIndex @1 :Int64;
    level @2 :UInt32;
    parentId @3 :Text;
    activeBlocks @4 :List(Int32);
}

struct BlockInfo {
    blockId @0 :Int32;
    rangeStart @1 :Int64;
    rangeEnd @2 :Int64;
    seedCount @3 :UInt32;
}

struct Index {
    k @0 :Int32;
    s @1 :Int32;
    t @2 :Int32;
    l @3 :Int32;
    open @4 :Bool;
    perNodeSeedMutations @5 :List(SeedMutations);
    perNodeGapMutations @6 :List(GapMutations);
    nodePathInfo @7 :List(NodePathInfo);
    blockInfo @8 :List(BlockInfo);
    ancestorMatrix @9 :List(List(Bool));
    kmerDictionary @10 :List(KmerDictionaryEntry);  # Global dictionary of unique k-mers
}
