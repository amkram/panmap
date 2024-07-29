@0xf6f6f6f6f6f6f6f6;

struct Insertion {
    pos :union {
        pos32 @0 :Int32;
        pos64 @1 :Int64;
    }
}
struct Deletion {
    pos :union {
        pos32 @0 :Int32;
        pos64 @1 :Int64;
    }
}

struct InsertionWithOffset {
    bitset @0 :Int64;
    pos :union {
        pos32 @1 :Int32;
        pos64 @2 :Int64;
    }
}

struct DeletionWithOffset {
    bitset @0 :Int64;
    pos :union {
        pos32 @1 :Int32;
        pos64 @2 :Int64;
    }
}

struct Mutations {
    insertions @0 :List(Insertion);
    deletions @1 :List(Deletion);
    insertionsWithOffset @2 :List(InsertionWithOffset);
    deletionsWithOffset @3 :List(DeletionWithOffset);
}

struct Index {
    k @0 :Int32;
    s @1 :Int32;
    width @2 :Int16;
    perNodeMutations @3 :List(Mutations);
}