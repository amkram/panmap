#include "zstd_compression.hpp"
#include "logging.hpp"
#include <zstd.h>
#include <algorithm>
#include <atomic>
#include <fstream>
#include <thread>
#include <utility>
#include <vector>
#include <tbb/parallel_for.h>
#include <boost/iostreams/device/mapped_file.hpp>

namespace panmap_zstd {

bool compressToFile(const void* inputData,
                    size_t inputSize,
                    const std::string& outputPath,
                    int compressionLevel,
                    int numThreads,
                    size_t frameSize,
                    const void* headerData,
                    size_t headerSize) {
    if (numThreads == 0) {
        numThreads = std::thread::hardware_concurrency();
    }
    if (frameSize == 0) {
        frameSize = 64ull * 1024 * 1024;
    }

    // Independent per-chunk frames (not one nbWorkers frame) so decompressFromFile
    // can inflate them in parallel -- single-frame inflate dominates large-index load.
    const size_t numChunks = std::max<size_t>(1, (inputSize + frameSize - 1) / frameSize);

    logging::info("Compressing {} bytes to {} with ZSTD level {} using {} threads ({} frames x {} MB)",
                  inputSize,
                  outputPath,
                  compressionLevel,
                  numThreads,
                  numChunks,
                  frameSize / (1024 * 1024));

    auto startTime = std::chrono::high_resolution_clock::now();

    const uint8_t* in = static_cast<const uint8_t*>(inputData);
    std::vector<std::vector<uint8_t>> frames(numChunks);
    std::atomic<bool> anyError{false};

    tbb::parallel_for(size_t(0), numChunks, [&](size_t i) {
        if (anyError.load(std::memory_order_relaxed)) return;
        const size_t inOff = i * frameSize;
        const size_t inLen = std::min(frameSize, inputSize - inOff);
        std::vector<uint8_t> out(ZSTD_compressBound(inLen));
        ZSTD_CCtx* cctx = ZSTD_createCCtx();
        if (!cctx) {
            anyError.store(true, std::memory_order_relaxed);
            return;
        }
        ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, compressionLevel);
        ZSTD_CCtx_setParameter(cctx, ZSTD_c_checksumFlag, 1);  // one-shot => content-size header written
        const size_t r = ZSTD_compress2(cctx, out.data(), out.size(), in + inOff, inLen);
        ZSTD_freeCCtx(cctx);
        if (ZSTD_isError(r)) {
            logging::err("ZSTD compression failed on frame {}: {}", i, ZSTD_getErrorName(r));
            anyError.store(true, std::memory_order_relaxed);
            return;
        }
        out.resize(r);
        out.shrink_to_fit();  // free the compressBound slack so only ~compressed size is held
        frames[i] = std::move(out);
    });

    if (anyError.load()) {
        return false;
    }

    std::ofstream outFile(outputPath, std::ios::binary);
    if (!outFile) {
        logging::err("Failed to open output file: {}", outputPath);
        return false;
    }
    if (headerData != nullptr && headerSize > 0) {  // uncompressed prefix (e.g. index params)
        outFile.write(reinterpret_cast<const char*>(headerData), static_cast<std::streamsize>(headerSize));
    }
    size_t compressedSize = 0;
    for (const auto& f : frames) {  // concatenate frames in order
        outFile.write(reinterpret_cast<const char*>(f.data()), static_cast<std::streamsize>(f.size()));
        if (!outFile) {
            logging::err("Failed to write compressed data to file: {}", outputPath);
            return false;
        }
        compressedSize += f.size();
    }
    outFile.close();

    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);

    double ratio = 100.0 * compressedSize / inputSize;
    double throughput = inputSize / (1024.0 * 1024.0 * duration.count() / 1000.0);

    logging::info("Compression complete: {} -> {} bytes ({:.1f}% ratio) in {}ms ({:.1f} MB/s)",
                  inputSize,
                  compressedSize,
                  ratio,
                  duration.count(),
                  throughput);

    return true;
}

bool decompressFromFile(const std::string& inputPath, std::vector<uint8_t>& outputData, int numThreads,
                        size_t dataOffset) {
    if (numThreads == 0) {
        numThreads = std::thread::hardware_concurrency();
    }

    logging::debug("Decompressing {} using {} threads", inputPath, numThreads);

    auto startTime = std::chrono::high_resolution_clock::now();

    boost::iostreams::mapped_file_source mappedFile;
    try {
        mappedFile.open(inputPath);
    } catch (const std::exception& e) {
        logging::err("Failed to memory map input file: {} ({})", inputPath, e.what());
        return false;
    }

    if (!mappedFile.is_open()) {
        logging::err("Failed to open input file: {}", inputPath);
        return false;
    }

    const char* compressedData = mappedFile.data();
    size_t compressedSize = mappedFile.size();

    // Total size = sum over frames; getFrameContentSize on the whole buffer reports
    // only the first frame (indexes are written multi-frame).
    size_t numFrames = 0;
    unsigned long long decompressedSize = 0;
    size_t pos = dataOffset;  // skip an uncompressed header prefix, if any
    while (pos < compressedSize) {
        unsigned long long frameCompSize = ZSTD_findFrameCompressedSize(compressedData + pos, compressedSize - pos);
        if (ZSTD_isError(frameCompSize)) {
            break;
        }
        unsigned long long frameContentSize = ZSTD_getFrameContentSize(compressedData + pos, frameCompSize);
        if (frameContentSize == ZSTD_CONTENTSIZE_ERROR) {
            logging::err("Not a valid ZSTD frame in: {}", inputPath);
            return false;
        }
        if (frameContentSize == ZSTD_CONTENTSIZE_UNKNOWN) {
            logging::err("Cannot determine uncompressed size for a frame in: {}", inputPath);
            return false;
        }
        decompressedSize += frameContentSize;
        numFrames++;
        pos += frameCompSize;
    }
    if (numFrames == 0) {
        logging::err("No valid ZSTD frames in: {}", inputPath);
        return false;
    }

    outputData.resize(decompressedSize);

    if (numFrames > 1 && numThreads > 1) {
        logging::debug("Detected {} ZSTD frames, decompressing in parallel", numFrames);

        std::vector<std::pair<size_t, size_t>> frameRanges;
        std::vector<std::pair<size_t, size_t>> outputRanges;

        pos = dataOffset;  // frames begin after the uncompressed header
        size_t outPos = 0;
        while (pos < compressedSize) {
            unsigned long long frameCompSize = ZSTD_findFrameCompressedSize(compressedData + pos, compressedSize - pos);
            if (ZSTD_isError(frameCompSize)) break;

            unsigned long long frameDecompSize = ZSTD_getFrameContentSize(compressedData + pos, frameCompSize);
            if (ZSTD_isError(frameDecompSize) || frameDecompSize == ZSTD_CONTENTSIZE_UNKNOWN) break;

            frameRanges.push_back({pos, frameCompSize});
            outputRanges.push_back({outPos, frameDecompSize});

            pos += frameCompSize;
            outPos += frameDecompSize;
        }

        std::atomic<bool> error{false};
        tbb::parallel_for(size_t(0), frameRanges.size(), [&](size_t i) {
            if (error.load()) return;

            ZSTD_DCtx* dctx = ZSTD_createDCtx();
            if (!dctx) {
                error.store(true);
                return;
            }

            size_t result = ZSTD_decompressDCtx(dctx,
                                                outputData.data() + outputRanges[i].first,
                                                outputRanges[i].second,
                                                compressedData + frameRanges[i].first,
                                                frameRanges[i].second);

            ZSTD_freeDCtx(dctx);

            if (ZSTD_isError(result) || result != outputRanges[i].second) {
                error.store(true);
            }
        });

        if (error.load()) {
            logging::err("Parallel ZSTD decompression failed");
            mappedFile.close();
            return false;
        }

    } else {
        logging::debug("Decompressing {} frame(s) single-threaded", numFrames);

        ZSTD_DCtx* dctx = ZSTD_createDCtx();
        if (!dctx) {
            logging::err("Failed to create ZSTD decompression context");
            return false;
        }

        // ZSTD_decompressDCtx handles one frame; walk them for multi-frame input.
        size_t inPos = dataOffset;  // frames begin after the uncompressed header
        size_t outPos = 0;
        while (inPos < compressedSize) {
            unsigned long long frameCompSize =
                ZSTD_findFrameCompressedSize(compressedData + inPos, compressedSize - inPos);
            if (ZSTD_isError(frameCompSize)) break;
            unsigned long long frameContentSize = ZSTD_getFrameContentSize(compressedData + inPos, frameCompSize);
            size_t const r = ZSTD_decompressDCtx(dctx,
                                                 outputData.data() + outPos,
                                                 static_cast<size_t>(frameContentSize),
                                                 compressedData + inPos,
                                                 static_cast<size_t>(frameCompSize));
            if (ZSTD_isError(r) || r != frameContentSize) {
                logging::err("ZSTD decompression failed: {}",
                             ZSTD_isError(r) ? ZSTD_getErrorName(r) : "frame size mismatch");
                ZSTD_freeDCtx(dctx);
                mappedFile.close();
                return false;
            }
            inPos += frameCompSize;
            outPos += frameContentSize;
        }

        ZSTD_freeDCtx(dctx);

        if (outPos != decompressedSize) {
            logging::err("Decompressed size mismatch: expected {}, got {}", decompressedSize, outPos);
            mappedFile.close();
            return false;
        }
    }

    mappedFile.close();

    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);

    double throughput = decompressedSize / (1024.0 * 1024.0 * duration.count() / 1000.0);

    logging::debug("Decompression complete: {} -> {} bytes in {}ms ({:.1f} MB/s)",
                   compressedSize,
                   decompressedSize,
                   duration.count(),
                   throughput);

    return true;
}

}  // namespace panmap_zstd
