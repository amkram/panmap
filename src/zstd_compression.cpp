#include "zstd_compression.hpp"
#include "logging.hpp"
#include <zstd.h>
#include <fstream>
#include <cstring>
#include <thread>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <boost/iostreams/device/mapped_file.hpp>
#include <sys/mman.h>

namespace panmap_zstd {

// ZSTD magic number: 0xFD2FB528 (little-endian)
static const uint32_t ZSTD_MAGIC = 0x28B52FFD;

bool isSegmentedZstd(const std::string& inputPath) {
    std::ifstream file(inputPath, std::ios::binary);
    if (!file) return false;
    
    uint32_t magic;
    file.read(reinterpret_cast<char*>(&magic), sizeof(magic));
    return file.good() && magic == SEGMENTED_MAGIC;
}

size_t getUncompressedSize(const std::string& inputPath) {
    std::ifstream file(inputPath, std::ios::binary);
    if (!file) return 0;
    
    file.seekg(0, std::ios::end);
    size_t compressedSize = file.tellg();
    file.seekg(0, std::ios::beg);
    
    std::vector<uint8_t> compressedData(compressedSize);
    file.read(reinterpret_cast<char*>(compressedData.data()), compressedSize);
    
    unsigned long long const decompressedSize = ZSTD_getFrameContentSize(
        compressedData.data(), compressedData.size());
    
    if (decompressedSize == ZSTD_CONTENTSIZE_ERROR) {
        logging::err("Not a valid ZSTD frame: {}", inputPath);
        return 0;
    }
    if (decompressedSize == ZSTD_CONTENTSIZE_UNKNOWN) {
        logging::warn("Uncompressed size unknown for: {}", inputPath);
        return 0;
    }
    
    return static_cast<size_t>(decompressedSize);
}

bool compressSegmentedToFile(
    const std::vector<std::pair<const void*, size_t>>& segments,
    const std::string& outputPath,
    int compressionLevel,
    int numThreads
) {
    if (numThreads == 0) {
        numThreads = std::thread::hardware_concurrency();
    }
    
    size_t totalUncompressed = 0;
    for (const auto& [data, size] : segments) {
        totalUncompressed += size;
    }
    
    logging::info("Compressing {} segments ({} bytes total) to {} with ZSTD level {} using {} threads",
                  segments.size(), totalUncompressed, outputPath, compressionLevel, numThreads);
    
    // Debug: Log compression level received
    logging::info("DEBUG [compressSegmentedToFile]: compressionLevel parameter = {}", compressionLevel);
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // Compress all segments in parallel
    std::vector<std::vector<uint8_t>> compressedSegments(segments.size());
    std::vector<SegmentInfo> segmentTable(segments.size());
    
    // Use thread-local contexts to avoid reallocation overhead
    tbb::enumerable_thread_specific<ZSTD_CCtx*> threadContexts([compressionLevel]() {
        ZSTD_CCtx* cctx = ZSTD_createCCtx();
        ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, compressionLevel);
        ZSTD_CCtx_setParameter(cctx, ZSTD_c_checksumFlag, 1);
        return cctx;
    });

    tbb::parallel_for(0, (int)segments.size(), [&](int i) {
        const auto& [data, size] = segments[i];
        
        ZSTD_CCtx* cctx = threadContexts.local();
        
        // Compress segment
        size_t maxCompressedSize = ZSTD_compressBound(size);
        compressedSegments[i].resize(maxCompressedSize);
        
        size_t compressedSize = ZSTD_compress2(
            cctx,
            compressedSegments[i].data(), maxCompressedSize,
            data, size
        );
        
        if (ZSTD_isError(compressedSize)) {
            logging::err("ZSTD compression failed for segment {}: {}", i, ZSTD_getErrorName(compressedSize));
            compressedSegments[i].clear();
        } else {
            compressedSegments[i].resize(compressedSize);
            segmentTable[i].uncompressedSize = size;
            segmentTable[i].compressedSize = compressedSize;
        }
    });
    
    // Cleanup contexts
    for (auto cctx : threadContexts) {
        ZSTD_freeCCtx(cctx);
    }
    
    // Check for compression errors
    for (size_t i = 0; i < compressedSegments.size(); i++) {
        if (compressedSegments[i].empty()) {
            return false;
        }
    }
    
    // Write file: [Magic][NumSegments][SegmentTable][Segment0][Segment1]...
    std::ofstream outFile(outputPath, std::ios::binary);
    if (!outFile) {
        logging::err("Failed to open output file: {}", outputPath);
        return false;
    }
    
    // Write header
    uint32_t magic = SEGMENTED_MAGIC;
    uint32_t numSegments = segments.size();
    outFile.write(reinterpret_cast<const char*>(&magic), sizeof(magic));
    outFile.write(reinterpret_cast<const char*>(&numSegments), sizeof(numSegments));
    
    // Calculate offsets and write segment table
    size_t headerSize = sizeof(magic) + sizeof(numSegments) + 
                       numSegments * sizeof(SegmentInfo);
    size_t currentOffset = headerSize;
    
    for (size_t i = 0; i < numSegments; i++) {
        segmentTable[i].offset = currentOffset;
        currentOffset += segmentTable[i].compressedSize;
        outFile.write(reinterpret_cast<const char*>(&segmentTable[i]), sizeof(SegmentInfo));
    }
    
    // Write compressed segments
    size_t totalCompressed = 0;
    for (const auto& compressed : compressedSegments) {
        outFile.write(reinterpret_cast<const char*>(compressed.data()), compressed.size());
        totalCompressed += compressed.size();
    }
    
    outFile.close();
    
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    
    double ratio = 100.0 * totalCompressed / totalUncompressed;
    double throughput = totalUncompressed / (1024.0 * 1024.0 * duration.count() / 1000.0);
    
    logging::info("Segmented compression complete: {} -> {} bytes ({:.1f}% ratio) in {}ms ({:.1f} MB/s)",
                  totalUncompressed, totalCompressed, ratio, duration.count(), throughput);
    
    return true;
}

bool decompressSegmentedFromFile(
    const std::string& inputPath,
    std::vector<uint8_t>& outputData,
    int numThreads
) {
    if (numThreads == 0) {
        numThreads = std::thread::hardware_concurrency();
    }
    
    logging::info("Decompressing segmented file {} using {} threads", inputPath, numThreads);
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // Use mmap instead of reading entire file into vector
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
    
    const char* fileData = mappedFile.data();
    size_t fileSize = mappedFile.size();
    
    if (fileSize < sizeof(uint32_t) * 2) {
        logging::err("File too small: {}", inputPath);
        return false;
    }
    
    uint32_t magic = *reinterpret_cast<const uint32_t*>(fileData);
    uint32_t numSegments = *reinterpret_cast<const uint32_t*>(fileData + sizeof(uint32_t));
    
    if (magic != SEGMENTED_MAGIC) {
        logging::err("Invalid segmented ZSTD magic number: {:#x}", magic);
        return false;
    }
    
    // Read segment table
    // We copy it to a vector to be safe and easy to use
    std::vector<SegmentInfo> segmentTable(numSegments);
    const char* tablePtr = fileData + sizeof(uint32_t) * 2;
    // Check bounds
    if (fileSize < sizeof(uint32_t) * 2 + numSegments * sizeof(SegmentInfo)) {
        logging::err("File too small for segment table: {}", inputPath);
        return false;
    }
    
    std::memcpy(segmentTable.data(), tablePtr, numSegments * sizeof(SegmentInfo));
    
    // Calculate total uncompressed size
    size_t totalUncompressed = 0;
    for (const auto& seg : segmentTable) {
        totalUncompressed += seg.uncompressedSize;
    }
    
    // Allocate output buffer
    outputData.resize(totalUncompressed);
    
    // Decompress all segments in parallel
    std::atomic<bool> success{true};
    size_t outputOffset = 0;
    std::vector<size_t> outputOffsets(numSegments);
    for (uint32_t i = 0; i < numSegments; i++) {
        outputOffsets[i] = outputOffset;
        outputOffset += segmentTable[i].uncompressedSize;
    }
    
    // Thread-local decompression contexts
    tbb::enumerable_thread_specific<ZSTD_DCtx*> threadContexts([]() {
        return ZSTD_createDCtx();
    });
    
    tbb::parallel_for(0, (int)numSegments, [&](int i) {
        if (!success.load()) return;
        
        ZSTD_DCtx* dctx = threadContexts.local();
        
        // Check bounds for segment
        if (segmentTable[i].offset + segmentTable[i].compressedSize > fileSize) {
             logging::err("Segment {} out of file bounds", i);
             success.store(false);
             return;
        }

        size_t actualSize = ZSTD_decompressDCtx(
            dctx,
            outputData.data() + outputOffsets[i], segmentTable[i].uncompressedSize,
            fileData + segmentTable[i].offset, segmentTable[i].compressedSize
        );
        
        if (ZSTD_isError(actualSize)) {
            logging::err("ZSTD decompression failed for segment {}: {}", i, ZSTD_getErrorName(actualSize));
            success.store(false);
            return;
        }
        
        if (actualSize != segmentTable[i].uncompressedSize) {
            logging::err("Decompressed size mismatch for segment {}: expected {}, got {}",
                        i, segmentTable[i].uncompressedSize, actualSize);
            success.store(false);
            return;
        }
    });
    
    // Cleanup
    for (auto dctx : threadContexts) {
        ZSTD_freeDCtx(dctx);
    }
    
    mappedFile.close();
    
    if (!success.load()) {
        return false;
    }
    
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    
    double throughput = totalUncompressed / (1024.0 * 1024.0 * duration.count() / 1000.0);
    
    logging::info("Segmented decompression complete: {} segments, {} bytes in {}ms ({:.1f} MB/s)",
                  numSegments, totalUncompressed, duration.count(), throughput);
    
    return true;
}

bool decompressSegmentedFromFile(
    const std::string& inputPath,
    std::unique_ptr<uint8_t[]>& outputData,
    size_t& outputSize,
    int numThreads
) {
    if (numThreads == 0) {
        numThreads = std::thread::hardware_concurrency();
    }
    
    logging::info("Decompressing segmented file {} using {} threads (raw buffer)", inputPath, numThreads);
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // Use mmap instead of reading entire file into vector
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
    
    const char* fileData = mappedFile.data();
    size_t fileSize = mappedFile.size();
    
    // Hint to OS that we will access this sequentially (or randomly depending on access pattern)
    // Since we read segments which are sequential, but in parallel, it's mixed.
    // MADV_WILLNEED might be good if we have enough RAM.
    #ifdef MADV_WILLNEED
    madvise((void*)fileData, fileSize, MADV_WILLNEED);
    #endif
    
    if (fileSize < sizeof(uint32_t) * 2) {
        logging::err("File too small: {}", inputPath);
        return false;
    }
    
    uint32_t magic = *reinterpret_cast<const uint32_t*>(fileData);
    uint32_t numSegments = *reinterpret_cast<const uint32_t*>(fileData + sizeof(uint32_t));
    
    if (magic != SEGMENTED_MAGIC) {
        logging::err("Invalid segmented ZSTD magic number: {:#x}", magic);
        return false;
    }
    
    // Read segment table
    std::vector<SegmentInfo> segmentTable(numSegments);
    const char* tablePtr = fileData + sizeof(uint32_t) * 2;
    // Check bounds
    if (fileSize < sizeof(uint32_t) * 2 + numSegments * sizeof(SegmentInfo)) {
        logging::err("File too small for segment table: {}", inputPath);
        return false;
    }
    
    std::memcpy(segmentTable.data(), tablePtr, numSegments * sizeof(SegmentInfo));
    
    // Calculate total uncompressed size
    size_t totalUncompressed = 0;
    for (const auto& seg : segmentTable) {
        totalUncompressed += seg.uncompressedSize;
    }
    
    // Allocate output buffer using new[] to avoid zero-initialization
    try {
        outputData.reset(new uint8_t[totalUncompressed]);
        outputSize = totalUncompressed;
    } catch (const std::bad_alloc& e) {
        logging::err("Failed to allocate {} bytes for decompression: {}", totalUncompressed, e.what());
        return false;
    }
    
    // Decompress all segments in parallel
    std::atomic<bool> success{true};
    size_t outputOffset = 0;
    std::vector<size_t> outputOffsets(numSegments);
    for (uint32_t i = 0; i < numSegments; i++) {
        outputOffsets[i] = outputOffset;
        outputOffset += segmentTable[i].uncompressedSize;
    }
    
    // Thread-local decompression contexts
    tbb::enumerable_thread_specific<ZSTD_DCtx*> threadContexts([]() {
        return ZSTD_createDCtx();
    });
    
    tbb::parallel_for(0, (int)numSegments, [&](int i) {
        if (!success.load()) return;
        
        ZSTD_DCtx* dctx = threadContexts.local();
        
        // Check bounds for segment
        if (segmentTable[i].offset + segmentTable[i].compressedSize > fileSize) {
             logging::err("Segment {} out of file bounds", i);
             success.store(false);
             return;
        }

        size_t actualSize = ZSTD_decompressDCtx(
            dctx,
            outputData.get() + outputOffsets[i], segmentTable[i].uncompressedSize,
            fileData + segmentTable[i].offset, segmentTable[i].compressedSize
        );
        
        if (ZSTD_isError(actualSize)) {
            logging::err("ZSTD decompression failed for segment {}: {}", i, ZSTD_getErrorName(actualSize));
            success.store(false);
            return;
        }
        
        if (actualSize != segmentTable[i].uncompressedSize) {
            logging::err("Decompressed size mismatch for segment {}: expected {}, got {}",
                        i, segmentTable[i].uncompressedSize, actualSize);
            success.store(false);
            return;
        }
    });
    
    // Cleanup
    for (auto dctx : threadContexts) {
        ZSTD_freeDCtx(dctx);
    }
    
    mappedFile.close();
    
    if (!success.load()) {
        outputData.reset();
        outputSize = 0;
        return false;
    }
    
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    
    double throughput = totalUncompressed / (1024.0 * 1024.0 * duration.count() / 1000.0);
    
    logging::info("Segmented decompression complete: {} segments, {} bytes in {}ms ({:.1f} MB/s)",
                  numSegments, totalUncompressed, duration.count(), throughput);
    
    return true;
}

bool compressToFile(
    const void* inputData,
    size_t inputSize,
    const std::string& outputPath,
    int compressionLevel,
    int numThreads,
    size_t frameSize
) {
    if (numThreads == 0) {
        numThreads = std::thread::hardware_concurrency();
    }
    
    logging::info("Compressing {} bytes to {} with ZSTD level {} using {} threads (frame size: {} MB)",
                  inputSize, outputPath, compressionLevel, numThreads, frameSize / (1024 * 1024));
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // Open output file
    std::ofstream outFile(outputPath, std::ios::binary);
    if (!outFile) {
        logging::err("Failed to open output file: {}", outputPath);
        return false;
    }
    
    // Create compression context with multi-threading
    ZSTD_CCtx* cctx = ZSTD_createCCtx();
    if (!cctx) {
        logging::err("Failed to create ZSTD compression context");
        return false;
    }
    
    // Set compression parameters for parallel compression
    ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, compressionLevel);
    ZSTD_CCtx_setParameter(cctx, ZSTD_c_nbWorkers, numThreads);
    ZSTD_CCtx_setParameter(cctx, ZSTD_c_jobSize, frameSize);  // Set frame size for parallelism
    
    // Enable checksums for data integrity
    ZSTD_CCtx_setParameter(cctx, ZSTD_c_checksumFlag, 1);
    
    // Allocate output buffer (worst case: slightly larger than input)
    size_t const maxCompressedSize = ZSTD_compressBound(inputSize);
    std::vector<uint8_t> compressedData(maxCompressedSize);
    
    // Compress with parallel frames
    size_t const compressedSize = ZSTD_compress2(
        cctx,
        compressedData.data(), maxCompressedSize,
        inputData, inputSize
    );
    
    ZSTD_freeCCtx(cctx);
    
    if (ZSTD_isError(compressedSize)) {
        logging::err("ZSTD compression failed: {}", ZSTD_getErrorName(compressedSize));
        return false;
    }
    
    // Write to file
    outFile.write(reinterpret_cast<const char*>(compressedData.data()), compressedSize);
    if (!outFile) {
        logging::err("Failed to write compressed data to file: {}", outputPath);
        return false;
    }
    outFile.close();
    
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    
    double ratio = 100.0 * compressedSize / inputSize;
    double throughput = inputSize / (1024.0 * 1024.0 * duration.count() / 1000.0);
    
    logging::info("Compression complete: {} -> {} bytes ({:.1f}% ratio) in {}ms ({:.1f} MB/s)",
                  inputSize, compressedSize, ratio, duration.count(), throughput);
    
    return true;
}

bool decompressFromFile(
    const std::string& inputPath,
    std::vector<uint8_t>& outputData,
    int numThreads
) {
    if (numThreads == 0) {
        numThreads = std::thread::hardware_concurrency();
    }
    
    logging::info("Decompressing {} using {} threads", inputPath, numThreads);
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // Use mmap for input
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
    
    // Get uncompressed size
    unsigned long long const decompressedSize = ZSTD_getFrameContentSize(
        compressedData, compressedSize);
    
    if (decompressedSize == ZSTD_CONTENTSIZE_ERROR) {
        logging::err("Not a valid ZSTD frame: {}", inputPath);
        return false;
    }
    if (decompressedSize == ZSTD_CONTENTSIZE_UNKNOWN) {
        logging::err("Cannot determine uncompressed size: {}", inputPath);
        return false;
    }
    
    // Allocate output buffer
    outputData.resize(decompressedSize);
    
    // Check if the compressed data has multiple frames (parallel compression was used)
    // If so, we can decompress in parallel
    size_t numFrames = 0;
    size_t pos = 0;
    while (pos < compressedSize) {
        unsigned long long frameSize = ZSTD_findFrameCompressedSize(compressedData + pos, compressedSize - pos);
        if (ZSTD_isError(frameSize)) {
            break;
        }
        numFrames++;
        pos += frameSize;
    }
    
    if (numFrames > 1 && numThreads > 1) {
        // Parallel decompression of multiple frames
        logging::debug("Detected {} ZSTD frames, decompressing in parallel", numFrames);
        
        std::vector<std::pair<size_t, size_t>> frameRanges;  // (compressed offset, compressed size)
        std::vector<std::pair<size_t, size_t>> outputRanges; // (output offset, output size)
        
        pos = 0;
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
        
        // Decompress frames in parallel
        std::atomic<bool> error{false};
        tbb::parallel_for(size_t(0), frameRanges.size(), [&](size_t i) {
            if (error.load()) return;
            
            ZSTD_DCtx* dctx = ZSTD_createDCtx();
            if (!dctx) {
                error.store(true);
                return;
            }
            
            size_t result = ZSTD_decompressDCtx(
                dctx,
                outputData.data() + outputRanges[i].first, outputRanges[i].second,
                compressedData + frameRanges[i].first, frameRanges[i].second
            );
            
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
        // Single-threaded decompression (single frame or single thread requested)
        logging::debug("Using single-threaded decompression");
        
        ZSTD_DCtx* dctx = ZSTD_createDCtx();
        if (!dctx) {
            logging::err("Failed to create ZSTD decompression context");
            return false;
        }
        
        size_t const actualSize = ZSTD_decompressDCtx(
            dctx,
            outputData.data(), decompressedSize,
            compressedData, compressedSize
        );
        
        ZSTD_freeDCtx(dctx);
        
        if (ZSTD_isError(actualSize)) {
            logging::err("ZSTD decompression failed: {}", ZSTD_getErrorName(actualSize));
            mappedFile.close();
            return false;
        }
        
        if (actualSize != decompressedSize) {
            logging::err("Decompressed size mismatch: expected {}, got {}", 
                         decompressedSize, actualSize);
            mappedFile.close();
            return false;
        }
    }
    
    mappedFile.close();
    
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    
    double throughput = decompressedSize / (1024.0 * 1024.0 * duration.count() / 1000.0);
    
    logging::info("Decompression complete: {} -> {} bytes in {}ms ({:.1f} MB/s)",
                  compressedSize, decompressedSize, duration.count(), throughput);
    
    return true;
}

} // namespace panmap_zstd
