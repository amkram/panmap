@0xaeb74a09d04c9118;

struct CoordDelta {
  pos @0 :UInt32;
  maybeValue :union {
    value @1 :UInt32;
    none @2 :Void;
  }
}

struct SeedInsub {
  startPos @0 :UInt32;
  endPos @1 :UInt32;
  seedHashIndex @2 :UInt32;
  isReversed @3 :Bool;
}


struct SeedDeletion {
  startPos @0 :UInt32;
}

struct NodeChanges {
  nodeIndex @0 :UInt32;
  seedInsubs @1 :List(SeedInsub);
  seedDeletions @2 :List(SeedDeletion);
  coordDeltas @3 :List(CoordDelta);
}

struct MGSRIndex {
  k @0 :UInt16;
  s @1 :UInt16;
  t @2 :UInt16;
  l @3 :UInt16;
  open @4 :Bool;

  seedHashSet @5 :List(UInt64);
  perNodeChanges @6 :List(NodeChanges);
  perNodeGapMutations @7 :List(List(CoordDelta));
}