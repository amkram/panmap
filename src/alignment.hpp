#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <tuple>
#include "capnp/list.h"
#include "seeding.hpp"
#include <optional>
#include "index_lite.capnp.h"
#include "panman.hpp"
namespace alignment {
void getAnchors(
    std::vector<std::tuple<int64_t, int32_t, int>> &anchors,
    const std::vector<std::vector<seeding::seed_t>> &readSeeds,
    const std::vector<std::string> &readSequences,
    const std::unordered_map<
        size_t, std::pair<std::vector<uint32_t>, std::vector<uint32_t>>>
        &seedToRefPositions,
    int k);
float align(
    std::vector<std::optional<seeding::seed_t>> &bestNodeSeeds,
    std::unordered_map<size_t,
                       std::pair<std::vector<uint32_t>, std::vector<uint32_t>>>
        &seedToRefPositions,
    std::string &nodeSequence, panmanUtils::Node *node, panmanUtils::Tree *T,
    int32_t k, int32_t s, int32_t t, bool open, int32_t l,
    const std::string &reads1Path, const std::string &reads2Path,
    std::vector<std::vector<seeding::seed_t>> &readSeeds,
    std::vector<std::string> &readSequences,
    std::vector<std::string> &readQuals, std::vector<std::string> &readNames,
    std::string &samFileName, std::string &bamFileName, bool pairedEndReads,
    std::string &refFileName, std::string aligner);

} // namespace alignment