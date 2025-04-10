#pragma once

#include <cstddef>
#include <functional>
#include <utility>

namespace util {

/**
 * @brief Utility hash function for std::pair to use in unordered containers
 * 
 * @tparam T First type in the pair
 * @tparam U Second type in the pair
 */
template <typename T, typename U>
struct PairHash {
    std::size_t operator()(const std::pair<T, U>& p) const {
        return std::hash<T>()(p.first) ^ 
               (std::hash<U>()(p.second) << 1);
    }
};

} // namespace util
