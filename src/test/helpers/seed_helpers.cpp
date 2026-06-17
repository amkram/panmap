#include "seed_helpers.hpp"

#include "seeding.hpp"

#include <algorithm>
#include <fstream>
#include <random>
#include <stdexcept>

namespace ts {

indexUtils::SeedCountMap extractSeeds(std::string_view seq, int k, int s, bool countDuplicates) {
    indexUtils::SeedCountMap out;
    auto syncmers = seeding::rollingSyncmers(seq, k, s, /*open=*/false, /*t=*/0, /*returnAll=*/false);
    for (const auto& [hash, isReverse, isSyncmer, pos] : syncmers) {
        if (!isSyncmer) continue;
        if (countDuplicates) {
            out[hash]++;
        } else {
            out[hash] = 1;
        }
    }
    return out;
}

std::string readFasta(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open FASTA file: " + path);
    }
    std::string sequence, line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '>') continue;
        line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
        sequence += line;
    }
    return sequence;
}

std::vector<std::string> generateReads(const std::string& genome, int readLen, int numReads, uint64_t seed) {
    std::vector<std::string> reads;
    if (genome.size() < static_cast<size_t>(readLen)) return reads;

    std::mt19937_64 gen(seed);
    std::uniform_int_distribution<size_t> posDist(0, genome.size() - readLen);
    for (int i = 0; i < numReads; i++) {
        size_t pos = posDist(gen);
        std::string read = genome.substr(pos, readLen);
        size_t nCount = std::count(read.begin(), read.end(), 'N');
        if (nCount < read.size() / 4) {
            reads.push_back(std::move(read));
        }
    }
    return reads;
}

}  // namespace ts
