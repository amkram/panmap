#pragma once

#include <string>

// CMake provides TEST_DATA_DIR as a compile definition: absolute path to
// src/test/data. Fall back to a relative path for ad-hoc builds.
#ifndef TEST_DATA_DIR
#define TEST_DATA_DIR "src/test/data"
#endif

namespace ts {

inline std::string dataPath(const std::string& relative) {
    return std::string(TEST_DATA_DIR) + "/" + relative;
}

}  // namespace ts
