#!/bin/sh
# Patch panman source for compatibility with oneTBB 2021+ (macOS)
# On Linux with TBB 2019, the second argument should be "false" to skip TBB compat patches.
SRCDIR="$1"
TBB_COMPAT="${2:-true}"

if [ "$TBB_COMPAT" = "true" ]; then
    # Remove tbb/task_scheduler_init.h includes and usages (removed in oneTBB 2021)
    find "$SRCDIR/src" \( -name "*.hpp" -o -name "*.cpp" -o -name "*.h" \) \
        -exec sed -i.bak '/task_scheduler_init/d' {} +

    # Add std::hash specialization for std::pair (needed by tbb::concurrent_unordered_set in oneTBB 2021+)
    cat > "$SRCDIR/src/pair_hash_compat.hpp" << 'HEREDOC'
#pragma once
#include <functional>
#include <utility>
namespace std {
template<typename A, typename B>
struct hash<pair<A,B>> {
    size_t operator()(const pair<A,B>& p) const {
        return hash<A>{}(p.first) ^ (hash<B>{}(p.second) << 1);
    }
};
}
HEREDOC

    # Include the hash compat header at the top of common.hpp
    sed -i.bak '1i\
#include "pair_hash_compat.hpp"
' "$SRCDIR/src/common.hpp"
fi
