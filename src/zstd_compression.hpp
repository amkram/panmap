#pragma once

#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>

namespace panmap_zstd {

// numThreads 0 = auto-detect; compressionLevel 1-22.
[[nodiscard]] bool compressToFile(const void* inputData,
                                  size_t inputSize,
                                  const std::string& outputPath,
                                  int compressionLevel = 3,
                                  int numThreads = 0,
                                  size_t frameSize = 4 * 1024 * 1024,
                                  const void* headerData = nullptr,    // uncompressed prefix, written verbatim before the frames
                                  size_t headerSize = 0
);

// numThreads 0 = auto-detect; dataOffset skips a leading uncompressed header.
[[nodiscard]] bool
decompressFromFile(const std::string& inputPath, std::vector<uint8_t>& outputData, int numThreads = 0,
                   size_t dataOffset = 0);

}  // namespace panmap_zstd
