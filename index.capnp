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

struct Index {
    k @0 :Int32;
    s @1 :Int32;
    t @2 :Int32;
    l @3 :Int32;
    open @4 :Bool;
    perNodeSeedMutations @5 :List(SeedMutations);
    perNodeGapMutations @6 :List(GapMutations);
}