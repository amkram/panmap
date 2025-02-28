#ifndef INDEXING_HPP
#define INDEXING_HPP

#include "coordinates.hpp"
#include "gap_map.hpp"
#include "index.capnp.h"
#include "seed_annotated_tree.hpp"
#include <capnp/common.h>
#include <capnp/message.h>
#include <capnp/serialize.h>

namespace indexing {

// Helper functions for processing node seeds
void processSeedChanges(
    seed_annotated_tree::TraversalGlobalState &state,
    seed_annotated_tree::NodeLocalState &nodeState,
    coordinates::CoordinateTraverser &traverser, const PanmapParams &params,
    panmanUtils::Tree *T, panmanUtils::Node *node, int32_t lastOnBlock,
    const coordinates::tupleCoord_t *lastOnCoord,
    std::set<seed_annotated_tree::SeedChange,
             seed_annotated_tree::SeedChangeComparator> &seedChanges);

void undoSeedChanges(
    seed_annotated_tree::TraversalGlobalState &state,
    seed_annotated_tree::NodeLocalState &nodeState,
    coordinates::CoordinateTraverser &traverser,
    const std::set<seed_annotated_tree::SeedChange,
                   seed_annotated_tree::SeedChangeComparator> &seedChanges);

void storeSeedChanges(
    ::capnp::List<SeedMutations>::Builder &perNodeSeedMutations,
    ::capnp::List<GapMutations>::Builder &perNodeGapMutations, int64_t dfsIndex,
    const seed_annotated_tree::NodeLocalState &nodeState,
    const std::set<seed_annotated_tree::SeedChange,
                   seed_annotated_tree::SeedChangeComparator> &seedChanges);

void indexingTraversal(
    seed_annotated_tree::TraversalGlobalState &state,
    ::capnp::List<SeedMutations>::Builder &perNodeSeedMutations,
    ::capnp::List<GapMutations>::Builder &perNodeGapMutations,
    const PanmapParams &params, panmanUtils::Tree *T, panmanUtils::Node *node,
    coordinates::CoordinateTraverser &traverser);

void indexingTraversalInit(const coordinates::CoordinateManager &coordManager,
                           seed_annotated_tree::TraversalGlobalState &state);

/**
 * @brief Build panmap index for a PanMAN pangenome tree
 * @param T Tree to be indexed
 * @param index Builder for the index
 * @param k Length of k-mer seeds
 * @param s Syncmer s parameter
 */
void index(panmanUtils::Tree *T, ::Index::Builder &index, int k, int s);

} // namespace indexing

#endif // INDEXING_HPP