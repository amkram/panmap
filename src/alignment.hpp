#include "placement.hpp"
#include "seed_annotated_tree.hpp"
#include <capnp/message.h>
#include <capnp/serialize.h>

namespace alignment {
void getAnchors(
    std::vector<std::tuple<int64_t, int32_t, int>> &anchors,
    const std::vector<std::vector<seeding::seed>> &readSeeds,
    const std::vector<std::string> &readSequences,
    const std::unordered_map<
        size_t, std::pair<std::vector<uint32_t>, std::vector<uint32_t>>>
        &seedToRefPositions,
    int k);
float align(
    std::vector<std::optional<seeding::onSeedsHash>> &bestNodeSeeds,
    std::unordered_map<size_t,
                       std::pair<std::vector<uint32_t>, std::vector<uint32_t>>>
        &seedToRefPositions,
    std::string &nodeSequence, panmanUtils::Node *node, panmanUtils::Tree *T,
    int32_t k, int32_t s, int32_t t, bool open, int32_t l,
    ::capnp::List<SeedMutations>::Reader &perNodeSeedMutations_Reader,
    ::capnp::List<GapMutations>::Reader &perNodeGapMutations_Reader,
    const std::string &reads1Path, const std::string &reads2Path,
    std::vector<std::vector<seed>> &readSeeds,
    std::vector<std::string> &readSequences,
    std::vector<std::string> &readQuals, std::vector<std::string> &readNames,
    std::string &samFileName, std::string &bamFileName, bool pairedEndReads,
    std::string &refFileName, std::string aligner);

} // namespace alignment