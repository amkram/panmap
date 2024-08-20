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

struct TernaryMasks {
    masks @0 :List(UInt32);
}

struct SeedMutations {
    basePositions @0 :List(Int64);
    perPosMasks @1 :List(TernaryMasks);
}

struct Index {
    k @0 :Int32;
    s @1 :Int32;
    perNodeSeedMutations @2 :List(SeedMutations);
    perNodeGapMutations @3 :List(GapMutations);
}