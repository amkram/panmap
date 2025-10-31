#include <boost/iostreams/categories.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <lzma.h>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <csignal>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <unistd.h>
#include <vector>
#include <stack>
#include <sys/stat.h>
// Additional headers for robust signal handling
#include <signal.h>
#include <atomic>
// Headers for backtrace functionality
#include <execinfo.h>
#include <cxxabi.h>
// Header for gprof profiling cleanup
extern "C" void _mcleanup(void);

#include "capnp/message.h"
#include "capnp/serialize-packed.h"
#include "capnp/serialize.h"
#include "kj/io.h"
#include "kj/array.h"
#include "panman.capnp.h"

#include "genotyping.hpp"
#include "index.capnp.h"
#include "mgsr_index.capnp.h"
#include "indexing.hpp"
#include "logging.hpp"
#include "panman.hpp"
#include "placement.hpp"
#include "timing.hpp"
#include "seeding.hpp"
#include "mgsr.hpp"
#include "panmap_utils.hpp"

using namespace logging;
namespace po = boost::program_options;

// Namespace aliases
namespace fs = boost::filesystem;

// Global constants
static constexpr int DEFAULT_K = 32;
static constexpr int DEFAULT_S = 8;

// Thread safety
std::mutex outputMutex, indexMutex, placementMutex;
std::atomic<bool> shouldStop{false};

// Forward declarations
std::unique_ptr<::capnp::MessageReader> readCapnp(const std::string &path);
void writeCapnp(::capnp::MallocMessageBuilder &message, const std::string &path);

// Create a custom wrapper to hold ownership of fd and reader
struct IndexReaderWrapper : public ::capnp::MessageReader {
    int fd_;
    std::unique_ptr<::capnp::PackedFdMessageReader> reader;
    
    IndexReaderWrapper(int fd, const ::capnp::ReaderOptions& opts)
        : ::capnp::MessageReader(opts), fd_(fd), reader(std::make_unique<::capnp::PackedFdMessageReader>(fd, opts)) {}
    
    ~IndexReaderWrapper() {
        // First destroy the reader (which might use the fd)
        reader.reset();
        // Then close the fd
        if (fd_ >= 0) {
            close(fd_);
            fd_ = -1;
        }
    }
    
    // Implement required MessageReader methods by delegating to the inner reader
    ::kj::ArrayPtr<const ::capnp::word> getSegment(uint id) override {
        return reader->getSegment(id);
    }
    
    // Return root object from the reader
    template <typename T>
    typename T::Reader getRoot() {
        return reader->getRoot<T>();
    }
};

void writeCapnp(::capnp::MallocMessageBuilder &message, const std::string &path) {
  // Normalize path to ensure consistent resolution
  boost::filesystem::path normalizedPath = boost::filesystem::absolute(path);
  std::string resolvedPath = normalizedPath.string();
  boost::filesystem::path parentDir = normalizedPath.parent_path();
  
  // Create a temporary file path in the same directory
  std::string tempPath = resolvedPath + ".tmp";
  logging::debug("Writing Cap'n Proto message to temporary file: {}", tempPath);

  // Create directory if it doesn't exist
  if (!parentDir.empty() && !boost::filesystem::exists(parentDir)) {
    logging::debug("Creating parent directory: {}", parentDir.string());
    if (!boost::filesystem::create_directories(parentDir)) {
      std::string error = "Failed to create directory: " + parentDir.string();
      logging::err("{}", error);
      throw std::runtime_error(error);
    }
  }

  // ---- CRITICAL: Validate message integrity before writing ----
  // Check key fields of the index to ensure they're not corrupted
  auto indexRoot = message.getRoot<Index>();
  uint32_t k = indexRoot.getK();
  uint32_t s = indexRoot.getS();
  size_t seedMutations = indexRoot.getPerNodeSeedMutations().size();
  size_t dictSize = 0;
  if (indexRoot.hasKmerDictionary()) {
    dictSize = indexRoot.getKmerDictionary().size();
  }

  // Keep the original k and s values for later reference
  const uint32_t original_k = k;
  const uint32_t original_s = s;
  
  // Check for potential corruption
  bool needsRepair = false;
  if (k == 0 || s == 0) {
    logging::critical("CRITICAL ERROR: Corrupt values detected before writing: k={}, s={}", k, s);
    needsRepair = true;
    
    // Try to guess what k and s should be
    k = k == 0 ? 32 : k; // Default to 32 if k is zero
    s = s == 0 ? 8 : s;  // Default to 8 if s is zero
    
    // Attempt to repair the message
    logging::warn("ATTEMPTING REPAIR: Setting k={}, s={}", k, s);
    indexRoot.setK(k);
    indexRoot.setS(s);
    
    // Re-read the values to verify repair
    auto indexAfterRepair = message.getRoot<Index>();
    uint32_t k_after_repair = indexAfterRepair.getK();
    uint32_t s_after_repair = indexAfterRepair.getS();
    
    logging::info("After repair: k={}, s={}", k_after_repair, s_after_repair);
    if (k_after_repair == 0 || s_after_repair == 0) {
      logging::critical("REPAIR FAILED: Unable to set non-zero values in Cap'n Proto message");
      // Continue anyway as a last resort
    }
  }
  
  // CRITICAL: If seedMutations or dictionary count is corrupted, we should abort and not write the file
  if (seedMutations == 0 && dictSize == 0) {
    logging::critical("CRITICAL ERROR: Message data is corrupted (empty seedMutations and dictionary)");
    throw std::runtime_error("Cannot write corrupted index with empty data");
  }
  
  logging::debug("Integrity check: k={}, s={}, seedMutations={}, dictionary={}", 
                indexRoot.getK(), indexRoot.getS(), seedMutations, dictSize);
  // ---- End integrity check ----

  // Write to temporary file first
  int fd = open(tempPath.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
  if (fd == -1) {
    std::string error = "Failed to open file for writing: " + tempPath + 
                       " (errno=" + std::to_string(errno) + ": " + 
                       std::strerror(errno) + ")";
    logging::err("{}", error);
    throw std::runtime_error(error);
  }

  try {
    // Write the message directly - no copying
    auto rootBeforeWrite = message.getRoot<Index>();
    uint32_t k_before_write = rootBeforeWrite.getK();
    uint32_t s_before_write = rootBeforeWrite.getS();
    size_t seedMutations_before_write = rootBeforeWrite.getPerNodeSeedMutations().size();
    size_t dictionary_before_write = rootBeforeWrite.getKmerDictionary().size();
    
    // Double check critical values one more time
    if (seedMutations_before_write == 0 && seedMutations > 0) {
      logging::critical("CRITICAL ERROR: seedMutations count was corrupted from {} to 0", seedMutations);
      throw std::runtime_error("Cannot write index with corrupted seedMutations count");
    }
    
    // Log key parameters to verify correct data is being written
    logging::info("Writing message: k={}, s={}, seedMutations={}, dictionary={}",
                 k_before_write, s_before_write, 
                 seedMutations_before_write, dictionary_before_write);
    
    // One last check - if values are still zero despite repair attempts, log and proceed
    if (k_before_write == 0 || s_before_write == 0) {
      logging::critical("CRITICAL ERROR: Zero values detected at write time, even after repair attempts. Will try to write anyway.");
    }
    
    // Use packed writing for efficiency
    ::capnp::writePackedMessageToFd(fd, message);
    
    // Explicitly sync to ensure data is written to disk
    if (fsync(fd) != 0) {
      std::string error = "Failed to fsync file: " + tempPath + 
                         " (errno=" + std::to_string(errno) + ": " + 
                         std::strerror(errno) + ")";
      logging::err("{}", error);
    close(fd);
      throw std::runtime_error(error);
    }
    
    // Close the file descriptor
    if (close(fd) != 0) {
      std::string error = "Failed to close file: " + tempPath + 
                         " (errno=" + std::to_string(errno) + ": " + 
                         std::strerror(errno) + ")";
      logging::err("{}", error);
      throw std::runtime_error(error);
    }
    
    // Verify written file by reading it back
    if (boost::filesystem::exists(tempPath) && 
        boost::filesystem::file_size(tempPath) > 0) {
        
      // Skip full validation in production to avoid performance impact
      // But we could add a flag to enable it
      
      logging::debug("Wrote {} bytes to temporary file", boost::filesystem::file_size(tempPath));
    } else {
      logging::warn("Temporary file is empty or doesn't exist after writing!");
    }
    
    // Sync the directory to ensure the file entry is updated
    int dirFd = open(parentDir.c_str(), O_RDONLY);
    if (dirFd != -1) {
      fsync(dirFd);
      close(dirFd);
    }
    
    // Atomically rename the temporary file to the target file
    if (std::rename(tempPath.c_str(), resolvedPath.c_str()) != 0) {
      std::string error = "Failed to rename file: " + tempPath + " -> " + 
                         resolvedPath + " (errno=" + std::to_string(errno) + 
                         ": " + std::strerror(errno) + ")";
      logging::err("{}", error);
      throw std::runtime_error(error);
    }
    
    // Final directory sync after rename
    dirFd = open(parentDir.c_str(), O_RDONLY);
    if (dirFd != -1) {
      fsync(dirFd);
      close(dirFd);
    }
    
    logging::info("Successfully wrote Cap'n Proto message to: {}", resolvedPath);
  } catch (const std::exception& e) {
    // Clean up temporary file on error
    if (fd != -1) {
    close(fd);
    }
    if (boost::filesystem::exists(tempPath)) {
      boost::filesystem::remove(tempPath);
    }
    std::string error = "Error writing Cap'n Proto message: " + std::string(e.what());
    logging::err("{}", error);
    throw std::runtime_error(error);
  }
}

std::unique_ptr<::capnp::MessageReader> readCapnp(const std::string &path) {
  // Normalize path to ensure consistent resolution
  boost::filesystem::path normalizedPath = boost::filesystem::absolute(path);
  std::string resolvedPath = normalizedPath.string();
  logging::debug("Reading Cap'n Proto message from: {}", resolvedPath);
  
  try {
    // Check if file exists and is not empty
    if (!fs::exists(resolvedPath)) {
      throw std::runtime_error("Index file not found: " + resolvedPath);
    }
    
    uintmax_t fileSize = fs::file_size(resolvedPath);
    if (fileSize == 0) {
      throw std::runtime_error("Index file is empty: " + resolvedPath);
    }
    logging::debug("Index file exists with size {} bytes", fileSize);
    
    // Check file permissions
    bool isReadable = access(resolvedPath.c_str(), R_OK) == 0;
    if (!isReadable) {
      int errnum = errno;
      throw std::runtime_error("Index file is not readable: " + resolvedPath + 
                             " (errno: " + std::to_string(errnum) + ", " + 
                             strerror(errnum) + ")");
    }
    
    // Open the file
    int fd = open(resolvedPath.c_str(), O_RDONLY);
    if (fd < 0) {
      int errnum = errno;
      throw std::runtime_error("Failed to open index file: " + resolvedPath + 
                            " (errno: " + std::to_string(errnum) + ", " + 
                            strerror(errnum) + ")");
    }
    
    try {
    // Configure reader options for large messages
    ::capnp::ReaderOptions opts;
    opts.traversalLimitInWords = std::numeric_limits<uint64_t>::max();
    opts.nestingLimit = 1024;
    
      // FIXED: Create wrapper that properly manages reader and fd lifetimes
      auto wrapper = std::make_unique<IndexReaderWrapper>(fd, opts);
      
      // Validate the root structure as MGSRIndex (new format)
      auto root = wrapper->getRoot<MGSRIndex>();
      
      // Validate essential fields
      uint32_t k = root.getK();
      uint32_t s = root.getS();
      auto liteTree = root.getLiteTree();
      size_t nodeCount = liteTree.getLiteNodes().size();
      size_t seedCount = root.getSeedInfo().size();
        
      if (k == 0 || s == 0 || nodeCount == 0) {
        // Don't close fd here - will be closed by wrapper destructor
        throw std::runtime_error("Invalid MGSRIndex data: k=" + std::to_string(k) + 
                              ", s=" + std::to_string(s) + 
              ", nodes=" + std::to_string(nodeCount));
        }
        
      // Log validation checks for MGSRIndex format
      logging::debug("Validated MGSRIndex structure: k={}, s={}, nodes={}, seeds={}, nodeChanges={}",
                    k, s, nodeCount, seedCount, root.getPerNodeChanges().size());
      
      logging::info("Successfully read and validated MGSRIndex from {}", resolvedPath);
      
      // Return the wrapper as a MessageReader (it will be upcast)
      return std::unique_ptr<::capnp::MessageReader>(wrapper.release());
    } catch (const ::kj::Exception& e) {
      // Clean up fd on exception
      close(fd);
      throw std::runtime_error("Cap'n Proto parsing error: " + 
                             std::string(e.getDescription().cStr()));
    } catch (const std::exception& e) {
      // Clean up fd on exception
      close(fd);
      throw; // Re-throw other exceptions
      }
  } catch (const std::exception &e) {
    logging::err("Failed to read Cap'n Proto message: {}", e.what());
    throw;
  }
}

/**
 * Demangles a C++ symbol name to make it human-readable
 * 
 * @param symbol The mangled symbol name
 * @return The demangled symbol name or the original if demangling fails
 */
std::string demangle(const char* symbol) {
    int status = 0;
    char* demangled = abi::__cxa_demangle(symbol, nullptr, nullptr, &status);
    
    if (status == 0 && demangled) {
        std::string result(demangled);
        free(demangled);
        return result;
    }
    
    return std::string(symbol);
}

/**
 * Prints the current stack trace
 * 
 * @param skip Number of frames to skip from the top of the stack
 */
void printStackTrace(int skip = 1) {
    void* callstack[128];
    int frames = backtrace(callstack, 128);
    char** symbols = backtrace_symbols(callstack, frames);
    
    std::cerr << "\nStack trace:" << std::endl;
    
    for (int i = skip; i < frames; i++) {
        std::string frame = symbols[i];
        
        // Try to extract the function name for demangling
        size_t nameStart = frame.find('(');
        size_t nameEnd = frame.find('+', nameStart);
        
        if (nameStart != std::string::npos && nameEnd != std::string::npos) {
            std::string mangledName = frame.substr(nameStart + 1, nameEnd - nameStart - 1);
            std::string prettyName = demangle(mangledName.c_str());
            
            // Replace the mangled name with the demangled one
            frame = frame.substr(0, nameStart + 1) + prettyName + frame.substr(nameEnd);
        }
        
        std::cerr << "  #" << (i - skip) << ": " << frame << std::endl;
    }
    
    free(symbols);
}

/**
 * Signal handler for handling interrupts (SIGINT, SIGTERM, etc.)
 * Simple and robust handler to ensure gmon.out gets saved
 * 
 * @param signum The signal number
 */
void signalHandler(int signum) {
    // Prevent recursive signal handling
    static volatile sig_atomic_t handling_signal = 0;
    if (handling_signal) {
        _exit(1);
    }
    handling_signal = 1;

    const char* sigName = "UNKNOWN";
    switch (signum) {
        case SIGINT:  sigName = "SIGINT (Interrupt/Ctrl+C)"; break;
        case SIGTERM: sigName = "SIGTERM (Termination request)"; break;
        case SIGSEGV: sigName = "SIGSEGV (Segmentation fault)"; break;
        case SIGABRT: sigName = "SIGABRT (Abort)"; break;
    }

    if (signum == SIGINT || signum == SIGTERM) {
        // Use write() instead of printf/fprintf for signal safety
        const char* msg = "\n*** Program interrupted by signal ***\n";
        write(STDERR_FILENO, msg, strlen(msg));
        
        // Set the flag to indicate we should stop
        shouldStop = true;
        
        const char* cleanup_msg = "Exiting gracefully to save profiling data (gmon.out)...\n";
        write(STDERR_FILENO, cleanup_msg, strlen(cleanup_msg));
        
        // Minimal cleanup - just try to shutdown spdlog quietly
        try {
            spdlog::set_level(spdlog::level::off);
            spdlog::shutdown();
        } catch (...) {
            // Ignore any exceptions during shutdown
        }
        
        // Use exit() (not _exit()) to trigger atexit handlers including gprof's profiling finalization
        // The key is that exit() will call gprof's internal cleanup which writes gmon.out
        exit(128 + signum);  // Standard Unix exit code for signal termination

    } else {
        // For fatal signals, exit immediately without cleanup
        const char* fatal_msg = "\n*** Fatal signal received, exiting immediately ***\n";
        write(STDERR_FILENO, fatal_msg, strlen(fatal_msg));
        _exit(1);
    }
}

bool isFileReadable(const std::string &path) {
  std::ifstream file(path);
  return file.good();
}

bool isFileWritable(const std::string &path) {
  if (fs::exists(path)) {
    std::ofstream file(path, std::ios::app);
    return file.good();
  }

  auto dir = fs::path(path).parent_path();
  if (!fs::exists(dir)) {
    try {
      fs::create_directories(dir);
    } catch (const fs::filesystem_error &) {
      return false;
    }
  }

  std::ofstream file(path);
  if (!file.good())
    return false;
  fs::remove(path);
  return true;
}

void validateInputFile(const std::string &path,
                       const std::string &description) {
  if (path.empty())
    throw std::runtime_error(description + " path is empty");
  if (!fs::exists(path))
    throw std::runtime_error(description + " not found: " + path);
  if (!isFileReadable(path))
    throw std::runtime_error("Cannot read " + description + ": " + path);
}

void validateOutputFile(const std::string &path,
                        const std::string &description) {
  if (path.empty())
    throw std::runtime_error(description + " path is empty");
  if (!isFileWritable(path))
    throw std::runtime_error("Cannot write to " + description + ": " + path);
}

void cleanupPanMAN(panmanUtils::TreeGroup *TG) { delete TG; }

// Function to get a random node from the TreeGroup
panmanUtils::Node *getRandomNode(panmanUtils::TreeGroup *TG,
                                 std::mt19937 &rng) {
  if (!TG || TG->trees.empty()) {
    return nullptr;
  }

  // First select a random tree
  std::uniform_int_distribution<size_t> treeDist(0, TG->trees.size() - 1);
  size_t treeIndex = treeDist(rng);
  panmanUtils::Tree &randomTree = TG->trees[treeIndex];

  // Then select a random node from that tree
  if (randomTree.allNodes.empty()) {
    return nullptr;
  }

  std::vector<panmanUtils::Node *> nodes;
  for (auto &pair : randomTree.allNodes) {
    nodes.push_back(pair.second);
  }

  std::uniform_int_distribution<size_t> nodeDist(0, nodes.size() - 1);
  size_t nodeIndex = nodeDist(rng);

  return nodes[nodeIndex];
}

// Function to save a node's sequence to a FASTA file
bool saveNodeSequence(panmanUtils::Tree *tree, panmanUtils::Node *node,
                      const std::string &outputFileName) {
  if (!tree || !node) {
    return false;
  }

  std::string sequence =
      tree->getStringFromReference(node->identifier, false, true);
  std::ofstream outFile(outputFileName);

  if (outFile.is_open()) {
    outFile << ">" << node->identifier << "\n";
    outFile << sequence << "\n";
    outFile.close();
    return true;
  }

  return false;
}

// wrapper function for timing functions
template <typename Func, typename... Args>
auto timeFunction(const std::string& name, Func&& func, Args&&... args) {
  auto start = std::chrono::high_resolution_clock::now();

    // Handle both void and non-void return types
  if constexpr (std::is_void_v<std::invoke_result_t<Func, Args...>>) {
    std::invoke(std::forward<Func>(func), std::forward<Args>(args)...);
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << name << " took: " << duration.count() << " microseconds\n";
  } else {
    auto result = std::invoke(std::forward<Func>(func), std::forward<Args>(args)...);
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << name << " took: " << duration.count() << " microseconds\n";
    return result;
  }
}

/**
 * @brief Efficient KJ InputStream wrapper for in-memory data
 * 
 * This provides a zero-copy KJ InputStream implementation that reads from
 * a memory buffer, avoiding the iostream overhead that was causing 71.6%
 * of the load time bottleneck.
 */
class MemoryInputStream: public kj::InputStream {
public:
  MemoryInputStream(const char* data, size_t size)
    : data_(data), size_(size), position_(0) {}

  size_t tryRead(void* buffer, size_t minBytes, size_t maxBytes) override {
    size_t available = size_ - position_;
    size_t amount = std::min(maxBytes, available);
    
    if (amount > 0) {
      std::memcpy(buffer, data_ + position_, amount);
      position_ += amount;
    }
    
    return amount;
  }

private:
  const char* data_;
  size_t size_;
  size_t position_;
};

/**
 * @brief Load a PanMAN tree from a file with optimized I/O
 *
 * This implementation eliminates the iostream bottleneck by:
 * 1. Using large buffers (4MB) for file reading
 * 2. Decompressing entire file to memory in large chunks (2MB)
 * 3. Creating a zero-copy KJ InputStream from the buffer
 * 4. Letting Cap'n Proto parse directly from memory
 *
 * @param filename Path to the PanMAN file
 * @return panmanUtils::TreeGroup* Loaded TreeGroup or nullptr if loading failed
 */
panmanUtils::TreeGroup *loadPanMAN(const std::string &filename) {
  const char* method_env = std::getenv("PANMAN_LOAD_METHOD");
  std::string method = method_env ? method_env : "parallel"; // Default to parallel

  if (method == "original") {
    // ========== ORIGINAL PANMAN METHOD ==========
    logging::info("Using ORIGINAL panman loading method");
    auto total_start = std::chrono::high_resolution_clock::now();
    
    std::ifstream inputFile(filename, std::ios::binary);
    if (!inputFile.is_open()) {
      logging::err("Failed to open PanMAN file: {}", filename);
      return nullptr;
    }
    
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inBuffer;
    inBuffer.push(boost::iostreams::lzma_decompressor());
    inBuffer.push(inputFile);
    std::istream inputStream(&inBuffer);
    
    panmanUtils::TreeGroup *TG = new panmanUtils::TreeGroup(inputStream);
    
    auto total_end = std::chrono::high_resolution_clock::now();
    auto total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        total_end - total_start).count();
    
    logging::info("Original method loaded PanMAN in {} ms", total_ms);
    return TG;

  } else if (method == "mmap") {
    // ========== MMAP-OPTIMIZED METHOD ==========
    logging::info("Using MMAP-OPTIMIZED panman loading method");
    try {
        auto total_start = std::chrono::high_resolution_clock::now();

        auto mmap_start = std::chrono::high_resolution_clock::now();
        boost::iostreams::mapped_file_source mappedFile(filename);
        if (!mappedFile.is_open()) {
            logging::err("Failed to memory-map PanMAN file: {}", filename);
            return nullptr;
        }
        auto mmap_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - mmap_start).count();
        logging::debug("Memory-mapped file in {} ms", mmap_ms);

  auto decompress_start = std::chrono::high_resolution_clock::now();
  boost::iostreams::filtering_streambuf<boost::iostreams::input> inBuffer;
  inBuffer.push(boost::iostreams::lzma_decompressor());
  // Wrap the mmap'd data with an array_source to avoid Boost optimal_buffer_size template issues
  boost::iostreams::array_source arr(mappedFile.data(), mappedFile.size());
  inBuffer.push(arr);
  std::istream inputStream(&inBuffer);

        auto parse_start = std::chrono::high_resolution_clock::now();
        panmanUtils::TreeGroup *TG = new panmanUtils::TreeGroup(inputStream);
        auto parse_end = std::chrono::high_resolution_clock::now();

        mappedFile.close();

        auto total_end = std::chrono::high_resolution_clock::now();
        auto total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            total_end - total_start).count();
        auto decompress_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            parse_start - decompress_start).count();
        auto parse_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            parse_end - parse_start).count();

        logging::info("MMAP-OPTIMIZED method loaded PanMAN: {} (total: {} ms, decompress: {} ms, parse: {} ms)", 
                      filename, total_ms, decompress_ms, parse_ms);
        
        return TG;
    } catch (const std::exception &e) {
        logging::err("Exception while loading PanMAN file with mmap: {}", e.what());
        return nullptr;
    }

  } else if (method == "lzma_parallel") {
    // ========== LZMA_PARALLEL METHOD ==========
    logging::info("Using LZMA_PARALLEL panman loading method");
    try {
        auto total_start = std::chrono::high_resolution_clock::now();

        // Open file
        std::ifstream inputFile(filename, std::ios::binary | std::ios::ate);
        if (!inputFile.is_open()) {
            logging::err("Failed to open PanMAN file: {}", filename);
            return nullptr;
        }
        std::streamsize fileSize = inputFile.tellg();
        inputFile.seekg(0, std::ios::beg);

        // Read the entire file into a buffer for parallel processing
        std::vector<uint8_t> file_buf(fileSize);
        inputFile.read(reinterpret_cast<char*>(file_buf.data()), fileSize);
        inputFile.close();

        // Find and decode the .xz index from the end of the file.
        // The .xz format stores the index before a 12-byte Stream Footer.
        // The Stream Footer contains a "backward size" field that indicates the size of the index.
        
        // Read the Stream Footer (last 12 bytes from the buffer)
        if (fileSize < 12) {
            logging::err("File is too small to be a valid XZ file: {}", filename);
            return nullptr;
        }
        uint8_t footer[12];
        std::memcpy(footer, file_buf.data() + fileSize - 12, 12);
        
        // Parse the backward size from the footer (bytes 4-7, little-endian)
        // The backward size is stored as (index_size_in_bytes / 4) - 1
        uint32_t backward_size_field = 0;
        for (int i = 0; i < 4; i++) {
            backward_size_field |= static_cast<uint32_t>(footer[4 + i]) << (i * 8);
        }
        
        // Calculate the actual size in bytes from the encoded field
        uint64_t index_size_bytes = (static_cast<uint64_t>(backward_size_field) + 1) * 4;
        
        // Point to the start of the index in the buffer
        if (static_cast<uint64_t>(fileSize) < index_size_bytes + 12) {
            logging::err("Invalid XZ index size in file: {}", filename);
            return nullptr;
        }
        uint8_t* index_ptr = file_buf.data() + fileSize - 12 - index_size_bytes;
        
        // Decode the index
        lzma_index *index = nullptr;
        uint64_t memlimit = UINT64_MAX;
        size_t in_pos = 0;
        
        lzma_ret ret = lzma_index_buffer_decode(&index, &memlimit, NULL, 
                                                 index_ptr, &in_pos, index_size_bytes);
        
        if (ret != LZMA_OK || index == nullptr) {
            logging::err("Failed to decode .xz index from file: {} (error code: {})", filename, ret);
            return nullptr;
        }
        
        uint64_t uncompressed_size = lzma_index_uncompressed_size(index);
        uint64_t num_blocks = lzma_index_block_count(index);
        
        logging::info("XZ file has {} blocks, uncompressed size: {} bytes", num_blocks, uncompressed_size);
        
        // Collect block information using lzma_index_iter
        struct BlockInfo {
            uint64_t compressed_offset;
            uint64_t uncompressed_offset;
            uint64_t compressed_size;
            uint64_t uncompressed_size;
            lzma_check check_type;
        };
        
        std::vector<BlockInfo> blocks;
        blocks.reserve(num_blocks);
        
        lzma_index_iter iter;
        lzma_index_iter_init(&iter, index);
        
        // First, iterate to the stream to get check type
        lzma_check check_type = LZMA_CHECK_NONE;
        if (!lzma_index_iter_next(&iter, LZMA_INDEX_ITER_STREAM)) {
            if (iter.stream.flags != nullptr) {
                check_type = iter.stream.flags->check;
            }
        }
        
        // Reset iterator to get blocks
        lzma_index_iter_init(&iter, index);
        
        while (!lzma_index_iter_next(&iter, LZMA_INDEX_ITER_BLOCK)) {
            BlockInfo info;
            info.compressed_offset = iter.block.compressed_file_offset;
            info.uncompressed_offset = iter.block.uncompressed_file_offset;
            info.compressed_size = iter.block.total_size;
            info.uncompressed_size = iter.block.uncompressed_size;
            info.check_type = check_type;
            blocks.push_back(info);
        }
        
        lzma_index_end(index, NULL);
        
        logging::debug("Collected {} blocks for parallel decompression", blocks.size());
        
        // Log first few blocks for debugging
        for (size_t i = 0; i < std::min(size_t(3), blocks.size()); ++i) {
            logging::debug("Block {}: compressed_offset={}, uncompressed_offset={}, compressed_size={}, uncompressed_size={}", 
                          i, blocks[i].compressed_offset, blocks[i].uncompressed_offset, 
                          blocks[i].compressed_size, blocks[i].uncompressed_size);
        }
        
        std::vector<char> decompressed_data(uncompressed_size);
        
        // Parallel block decompression using lzma_block_decoder()
        auto decompress_start = std::chrono::high_resolution_clock::now();
        
        std::atomic<bool> error_occurred{false};
        std::vector<std::atomic<double>> block_times(blocks.size());
        for (auto& t : block_times) {
            t.store(0.0);
        }
        
        tbb::parallel_for(size_t(0), blocks.size(), [&](size_t i) {
            if (error_occurred.load()) return;
            
            auto block_start = std::chrono::high_resolution_clock::now();
            
            const BlockInfo& block = blocks[i];
            
            // Bounds check
            if (block.compressed_offset + block.compressed_size > static_cast<uint64_t>(fileSize)) {
                logging::err("Block {} offset out of bounds: offset={}, size={}, fileSize={}", 
                           i, block.compressed_offset, block.compressed_size, fileSize);
                error_occurred = true;
                return;
            }
            
            const uint8_t* block_data = file_buf.data() + block.compressed_offset;
            
            // Decode block header size from first byte
            uint8_t header_size = lzma_block_header_size_decode(block_data[0]);
            if (header_size == 0 || header_size > LZMA_BLOCK_HEADER_SIZE_MAX) {
                logging::err("Invalid block header size for block {}", i);
                error_occurred = true;
                return;
            }
            
            // Prepare lzma_block structure with filter array
            lzma_filter filters[LZMA_FILTERS_MAX + 1];
            lzma_block lzma_block_struct;
            std::memset(&lzma_block_struct, 0, sizeof(lzma_block_struct));
            lzma_block_struct.version = 0;
            lzma_block_struct.check = block.check_type;
            lzma_block_struct.filters = filters;
            lzma_block_struct.header_size = header_size;
            
            // Decode block header to populate filters
            lzma_ret ret = lzma_block_header_decode(&lzma_block_struct, NULL, block_data);
            if (ret != LZMA_OK) {
                logging::err("Failed to decode block header for block {}: error {}", i, ret);
                error_occurred = true;
                return;
            }
            
            // Create a block decoder for this specific block
            lzma_stream strm = LZMA_STREAM_INIT;
            ret = lzma_block_decoder(&strm, &lzma_block_struct);
            if (ret != LZMA_OK) {
                logging::err("Failed to create block decoder for block {}: error {}", i, ret);
                error_occurred = true;
                return;
            }
            
            // Decompress this block
            strm.next_in = block_data + header_size;
            strm.avail_in = block.compressed_size - header_size;
            strm.next_out = reinterpret_cast<uint8_t*>(decompressed_data.data()) + block.uncompressed_offset;
            strm.avail_out = block.uncompressed_size;
            
            ret = lzma_code(&strm, LZMA_FINISH);
            lzma_end(&strm);
            
            if (ret != LZMA_STREAM_END) {
                logging::err("Failed to decompress block {}: error {}", i, ret);
                error_occurred = true;
                return;
            }
            
            auto block_end = std::chrono::high_resolution_clock::now();
            double block_time = std::chrono::duration<double>(block_end - block_start).count();
            block_times[i].store(block_time);
        });
        
        if (error_occurred) {
            logging::err("Parallel decompression failed");
            return nullptr;
        }
        
        auto decompress_end = std::chrono::high_resolution_clock::now();
        double total_decompress_time = std::chrono::duration<double>(decompress_end - decompress_start).count();
        
        // Calculate block timing statistics
        double min_time = std::numeric_limits<double>::max();
        double max_time = 0.0;
        double total_time = 0.0;
        for (size_t i = 0; i < block_times.size(); ++i) {
            double t = block_times[i].load();
            min_time = std::min(min_time, t);
            max_time = std::max(max_time, t);
            total_time += t;
        }
        double avg_time = total_time / block_times.size();
        
        logging::info("Block decompression stats: min={:.3f}s, max={:.3f}s, avg={:.3f}s, total_cpu={:.3f}s, wall={:.3f}s, speedup={:.2f}x",
                     min_time, max_time, avg_time, total_time, total_decompress_time, total_time / total_decompress_time);
        
        // Log slowest blocks if there are many
        if (blocks.size() > 10) {
            std::vector<std::pair<size_t, double>> sorted_blocks;
            for (size_t i = 0; i < blocks.size(); ++i) {
                sorted_blocks.push_back({i, block_times[i].load()});
            }
            std::sort(sorted_blocks.begin(), sorted_blocks.end(), 
                     [](const auto& a, const auto& b) { return a.second > b.second; });
            
            logging::debug("Top 5 slowest blocks:");
            for (size_t i = 0; i < std::min(size_t(5), sorted_blocks.size()); ++i) {
                size_t block_idx = sorted_blocks[i].first;
                double time = sorted_blocks[i].second;
                logging::debug("  Block {}: {:.3f}s (compressed: {} bytes, uncompressed: {} bytes)", 
                             block_idx, time, blocks[block_idx].compressed_size, blocks[block_idx].uncompressed_size);
            }
        }
        
        logging::debug("Parallel decompression completed");

        auto decompress_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - decompress_start).count();

        // Parse the final buffer
        auto parse_start = std::chrono::high_resolution_clock::now();
        kj::ArrayPtr<const capnp::word> words(
            reinterpret_cast<const capnp::word*>(decompressed_data.data()),
            decompressed_data.size() / sizeof(capnp::word));
        
        capnp::ReaderOptions readerOptions;
        readerOptions.traversalLimitInWords = std::numeric_limits<uint64_t>::max();
        readerOptions.nestingLimit = 1024;
        
        capnp::FlatArrayMessageReader messageReader(words, readerOptions);
        panman::TreeGroup::Reader TGReader = messageReader.getRoot<panman::TreeGroup>();
        
        auto treesList = TGReader.getTrees();
        size_t numTrees = treesList.size();
        std::vector<panmanUtils::Tree*> treePtrs(numTrees);
        
        // This part remains serial as it was faster
        for(size_t i = 0; i < numTrees; ++i) {
            treePtrs[i] = new panmanUtils::Tree(treesList[i]);
        }
        
        panmanUtils::TreeGroup *TG = new panmanUtils::TreeGroup(treePtrs);
        
        for (auto* tree : treePtrs) {
            delete tree;
        }
        
        auto parse_end = std::chrono::high_resolution_clock::now();
        auto parse_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            parse_end - parse_start).count();
        
        auto total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            parse_end - total_start).count();
        
        logging::info("LZMA_PARALLEL method loaded PanMAN: {} (total: {} ms, decompress: {} ms, parse: {} ms)", 
                      filename, total_ms, decompress_ms, parse_ms);
        
        return TG;

    } catch (const std::exception &e) {
        logging::err("Exception while loading PanMAN file with lzma_parallel: {}", e.what());
        return nullptr;
    }
  } else { // "parallel" or default
    // ========== PARALLEL-CONSTRUCTION-OPTIMIZED METHOD ==========
    logging::info("Using PARALLEL-CONSTRUCTION-OPTIMIZED panman loading method");
    try {
        auto total_start = std::chrono::high_resolution_clock::now();
        
        auto decompress_start = std::chrono::high_resolution_clock::now();
        std::ifstream inputFile(filename, std::ios::binary);
        if (!inputFile.is_open()) {
          logging::err("Failed to open PanMAN file: {}", filename);
          return nullptr;
        }

        inputFile.seekg(0, std::ios::end);
        size_t compressedSize = inputFile.tellg();
        inputFile.seekg(0, std::ios::beg);

        boost::iostreams::filtering_streambuf<boost::iostreams::input> inBuffer;
        inBuffer.push(boost::iostreams::lzma_decompressor());
        inBuffer.push(inputFile);
        std::istream inputStream(&inBuffer);

        std::string decompressedData;
        size_t estimatedSize = static_cast<size_t>(compressedSize * 6.5);
        decompressedData.reserve(estimatedSize);
        
        std::vector<char> buffer(4 * 1024 * 1024);
        while (inputStream.read(buffer.data(), buffer.size()) || inputStream.gcount() > 0) {
          decompressedData.append(buffer.data(), inputStream.gcount());
        }
        inputFile.close();
        auto decompress_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - decompress_start).count();

        auto parse_start = std::chrono::high_resolution_clock::now();
        kj::ArrayPtr<const capnp::word> words(
            reinterpret_cast<const capnp::word*>(decompressedData.data()),
            decompressedData.size() / sizeof(capnp::word));
        
        capnp::ReaderOptions readerOptions;
        readerOptions.traversalLimitInWords = std::numeric_limits<uint64_t>::max();
        readerOptions.nestingLimit = 1024;
        
        capnp::FlatArrayMessageReader messageReader(words, readerOptions);
        panman::TreeGroup::Reader TGReader = messageReader.getRoot<panman::TreeGroup>();
        
        auto treesList = TGReader.getTrees();
        size_t numTrees = treesList.size();
        std::vector<panmanUtils::Tree*> treePtrs(numTrees);
        
        tbb::parallel_for(tbb::blocked_range<size_t>(0, numTrees),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i != range.end(); ++i) {
                    treePtrs[i] = new panmanUtils::Tree(treesList[i]);
                }
            });
        
        panmanUtils::TreeGroup *TG = new panmanUtils::TreeGroup(treePtrs);
        
        for (auto* tree : treePtrs) {
            delete tree;
        }
        
        auto parse_end = std::chrono::high_resolution_clock::now();
        auto parse_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            parse_end - parse_start).count();
        
        auto total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            parse_end - total_start).count();
        
        logging::info("PARALLEL-CONSTRUCTION-OPTIMIZED method loaded PanMAN: {} (total: {} ms, decompress: {} ms, parse: {} ms)", 
                      filename, total_ms, decompress_ms, parse_ms);
        
        return TG;
    } catch (const std::exception &e) {
        logging::err("Exception while loading PanMAN file with parallel construction: {}", e.what());
        return nullptr;
    }
  }
}

// Function to save a node's constructed sequence to a FASTA file
bool saveNodeConstructedSequence(state::StateManager &stateManager, panmanUtils::Tree *tree, 
                               panmanUtils::Node *node, const std::string &outputFileName) {
  if (!tree || !node) {
    return false;
  }

  std::string nodeId = node->identifier;
  logging::info("Extracting constructed sequence for node {}", nodeId);
  
  try {
    // First ensure the node's blocks are properly activated
    auto activeBlocks = stateManager.getActiveBlocks(nodeId);
    if (activeBlocks.empty()) {
      logging::warn("Node {} has no active blocks. Sequence will be empty!", nodeId);
    } else {
      logging::info("Node {} has {} active blocks", nodeId, activeBlocks.size());
    }
    
    // Get the sequence using StateManager (this shows what's actually in memory)
    int64_t startPos = std::numeric_limits<int64_t>::max();
    int64_t endPos = std::numeric_limits<int64_t>::min();
    
    // Find the span of all active blocks
    for (int32_t blockId : activeBlocks) {
        auto blockRange = stateManager.getBlockRange(blockId);
        startPos = std::min(startPos, blockRange.start);
        endPos = std::max(endPos, blockRange.end);
    }
    
    // If no blocks are active, use a reasonable default
    if (startPos > endPos) {
        startPos = 0;
        endPos = 0;
    }
    
    coordinates::CoordRange fullRange = {startPos, endPos};
    auto [sequence, _, __, ___] = stateManager.extractSequence(nodeId, fullRange);
    
    logging::info("Extracted sequence of length {} for node {}", sequence.length(), nodeId);
    std::ofstream outFile(outputFileName);

    if (outFile.is_open()) {
      outFile << ">" << nodeId << "_constructed_sequence\n";
      outFile << sequence << "\n";
      outFile.close();
      return true;
    }
  } catch (const std::exception &e) {
    logging::err("Error extracting sequence for node {}: {}", nodeId, e.what());
  }

  return false;
}


// Main program entry point
int main(int argc, char *argv[]) {
    // Register signal handlers for various signals
    signal(SIGINT, signalHandler);   // Ctrl+C
    signal(SIGTERM, signalHandler);  // Termination request
    signal(SIGSEGV, signalHandler);  // Segmentation fault
    signal(SIGABRT, signalHandler);  // Abort
    
    // Register an atexit handler to ensure gprof data is written
    std::atexit([]() {
        // This ensures that any profiling data (gmon.out) gets written
        // even if we exit from a signal handler
        std::cerr << "Program cleanup completed - profiling data should be saved." << std::endl;
    });
    
    auto main_start = std::chrono::high_resolution_clock::now();
    
    // --- Boost.Program_options Setup ---
    po::options_description generic_opts("Generic options");
    generic_opts.add_options()
        ("help,h", "Show help message")
        ("version,V", "Show version")
        ("quiet,q", "Run in quiet mode")
        ("verbose,v", "Run in verbose mode")
        ("log-level", po::value<std::string>()->default_value("info"), 
         "Set logging level (trace, debug, info, warn, error, critical, off)")
        ("time", "Show time taken at each step")
        ("cpus,c", po::value<int>()->default_value(1), "Number of CPUs to use")
        ("stop-after,x", po::value<std::string>()->default_value(""), 
         "Stop after stage (indexing/i, placement/p, mapping/m, genotyping/g, assembly/a)")
        ("seed,Q", po::value<int>()->default_value(42), "Seed for random number generation")
        ("test", "Run coordinate system tests")
        ("test-nucleotide-mutations", "Run targeted test for nucleotide mutations in node_2")
    ;

    po::options_description input_opts("Input/output options");
    input_opts.add_options()
        ("prefix,p", po::value<std::string>()->default_value("panmap"), "Prefix for output files")
        ("outputs,o", po::value<std::string>()->default_value("bam,vcf,assembly"), 
         "Outputs (placement/p, assembly/a, reference/r, spectrum/c, sam/s, bam/b, mpileup/m, vcf/v, all/A)")
        ("index,i", po::value<std::string>()->default_value(""), "Path to precomputed index")
    ;
    
    po::options_description seeding_opts("Seeding/alignment options");
    seeding_opts.add_options()
        ("k,k", po::value<int>()->default_value(DEFAULT_K), "Length of k-mer seeds") // Use DEFAULT_K defined later
        ("s,s", po::value<int>()->default_value(DEFAULT_S), "Length of s-mers for syncmers") // Use DEFAULT_S defined later
        ("aligner,a", po::value<std::string>()->default_value("minimap2"), "Aligner (minimap2 or bwa-aln)")
        ("ref,r", po::value<std::string>()->default_value(""), "Reference node ID to align to (skip placement)")
        ("prior,P", "Use mutation spectrum prior for genotyping")
        ("reindex,f", "Force index rebuild")
        ("with-mm", "Build mutation matrix during indexing (off by default)")
        
    ;

    po::options_description mgsr_opts("MGSR options");
    mgsr_opts.add_options()
        ("index-mgsr", po::value<std::string>(), "Path to build/rebuild MGSR index")
        ("mgsr-index,m", po::value<std::string>(), "Path to precomputed MGSR index")
        ("l,l", po::value<int>()->default_value(1), "Length of k-min-mers (i.e. l seeds per kminmer)")
        ("no-progress", "Disable progress bars")
        ("seed-scores", "Use seed scores instead of read scores to score nodes")
        ("mask-reads", po::value<uint32_t>()->default_value(0), "mask reads containing k-min-mers with total occurrence <= threshold")
        ("mask-seeds", po::value<uint32_t>()->default_value(0), "mask k-min-mer seeds in query with total occurrence <= threshold")

        ("low-memory", "Use low memory mode")
    ;

    po::options_description dev_opts("Developer options");
    dev_opts.add_options()
        ("dump,D", "Dump all seeds to file")
        ("dump-real,X", "Dump true seeds to file")
        ("genotype-from-sam", "Generate VCF from SAM file")
        ("sam-file", po::value<std::string>(), "Path to SAM file for VCF generation")
        ("ref-file", po::value<std::string>(), "Path to reference FASTA for VCF generation")
        ("save-jaccard", "Save Jaccard index to <prefix>.jaccard.txt")
        ("save-kminmer-binary-coverage", "Save kminmer binary coverage")
        ("parallel-tester", "Run parallel tester")
        ("eval", po::value<std::string>(), "Evaluate placement accuracy (path to TSV)")
        ("dump-sequence", po::value<std::string>(), "Dump sequence for a specific node ID")
        ("dump-random-node", "Dump sequence for a random node")
        ("dump-node-cluster", po::value<std::vector<std::string>>()->multitoken(), "Dump N nearest relatives for specified node IDs (includes internal nodes)")
        ("dump-node-cluster-leaves", po::value<std::vector<std::string>>()->multitoken(), "Dump N nearest leaf-node relatives for specified node IDs (does not include internal nodes)")
        ("dump-random-node-cluster", po::value<uint32_t>(), "Dump sequences for a random node and N of its closest relatives (includes internal nodes)")
        ("dump-random-node-cluster-leaves", po::value<uint32_t>(), "Dump sequences for a random node and N of its closest leaf-node relatives (does not include internal nodes)")
        ("dump-random-node-clusters", po::value<std::vector<uint32_t>>()->multitoken(), "Dump sequences for multiple random node clusters, each with N closest relatives (includes internal nodes)")
        ("dump-random-node-clusters-leaves", po::value<std::vector<uint32_t>>()->multitoken(), "Dump sequences for multiple random node clusters, each with N closest leaf-node relatives (does not include internal nodes)")
        ("debug-node-id", po::value<std::string>()->default_value(""), "Log detailed placement debug info for this specific node ID")
        ("verify-scores", "Recompute all similarity scores from scratch at each node for verification (slow)")
        ("candidate-threshold", po::value<float>()->default_value(0.01f), "Placement candidate threshold proportion") 
        ("max-candidates", po::value<int>()->default_value(16), "Maximum placement candidates")
        ("overlap-coefficients", "Output overlap coefficients then exit")
        ("true-abundance", po::value<std::string>(), "Path to true abundance TSV for comparison")
        ("use-lite-placement", "Use LiteTree-based placement (avoids loading full panman until after placement)")
        ("use-full-placement", "Use traditional full Tree-based placement (default for backward compatibility)")
    ;

    // Hidden options for positional arguments
    po::options_description hidden_opts("Hidden options");
    hidden_opts.add_options()
        ("guide-panman", po::value<std::string>(), "Pangenome reference file")
        ("reads1", po::value<std::string>(), "Reads file 1")
        ("reads2", po::value<std::string>(), "Reads file 2")
    ;

    po::options_description cmdline_options;
    cmdline_options.add(generic_opts).add(input_opts).add(mgsr_opts).add(seeding_opts).add(dev_opts).add(hidden_opts);

    po::options_description visible_options("panmap -- v0.0 \u27d7 \u27d7\nPangenome phylogenetic placement, alignment, genotyping, and assembly of reads\n\nUsage: panmap [options] <guide.panman> [<reads1.fastq>] [<reads2.fastq>]\n\nAllowed options");
    visible_options.add(generic_opts).add(input_opts).add(mgsr_opts).add(seeding_opts).add(dev_opts);

    po::positional_options_description p;
    p.add("guide-panman", 1).add("reads1", 1).add("reads2", 1);

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
        po::notify(vm);
    } catch (const po::error& e) {
        std::cerr << "Error parsing command line: " << e.what() << std::endl;
        std::cerr << visible_options << std::endl;
        return 1;
    }

    // Handle --help
    if (vm.count("help")) {
        std::cout << visible_options << std::endl;
        return 0;
    }

    // Handle --version
    if (vm.count("version")) {
        std::cout << "panmap 0.0" << std::endl;
        return 0;
    }

    // --- End Boost.Program_options Setup ---

    // Set up parallel processing
    tbb::global_control c(tbb::global_control::max_allowed_parallelism,
                          vm["cpus"].as<int>());

    // Set logging verbosity level
    bool quiet_mode = vm.count("quiet") > 0;
    bool verbose_mode = vm.count("verbose") > 0;
    std::string levelStr = vm["log-level"].as<std::string>(); // Default handled by boost

    // Determine effective log level
    if (quiet_mode) {
        logging::loggingLevel = logging::LogLevel::CRITICAL; 
        levelStr = "critical"; // Update levelStr for logging message
    } else if (verbose_mode) {
        logging::loggingLevel = logging::LogLevel::TRACE;    
        levelStr = "trace"; // Update levelStr for logging message
    } else {
        // Map string from --log-level to enum
        std::string lowerLevelStr = levelStr;
        std::transform(lowerLevelStr.begin(), lowerLevelStr.end(), lowerLevelStr.begin(), ::tolower);
        
        if (lowerLevelStr == "trace") logging::loggingLevel = logging::LogLevel::TRACE;
        else if (lowerLevelStr == "debug") logging::loggingLevel = logging::LogLevel::DEBUG;
        else if (lowerLevelStr == "info") logging::loggingLevel = logging::LogLevel::INFO;
        else if (lowerLevelStr == "warn") logging::loggingLevel = logging::LogLevel::WARN;
        else if (lowerLevelStr == "error") logging::loggingLevel = logging::LogLevel::ERROR;
        else if (lowerLevelStr == "critical") logging::loggingLevel = logging::LogLevel::CRITICAL;
        else if (lowerLevelStr == "off") logging::loggingLevel = logging::LogLevel::OFF;
        else { 
            logging::loggingLevel = logging::LogLevel::WARN; // Default to WARN (quieter)
            levelStr = "warn"; // Ensure levelStr is consistent
            std::cerr << "Warning: Invalid log level '" << vm["log-level"].as<std::string>() << "' provided. Using default 'warn'." << std::endl;
        }
    }

    // Initialize spdlog (uses the static logging::loggingLevel variable)
    logging::initSpdlog(); 
    
    // Initialize node tracking only if verbose
    if (logging::loggingLevel <= logging::LogLevel::DEBUG) { 
        logging::initNodeTracking();
    }
    // Log the final effective level
    // Use levelStr which reflects quiet/verbose overrides
    logging::critical("Effective log level set to: {}", levelStr); 

    logging::critical("Starting panmap with verbosity level: {}",
                      quiet_mode ? "quiet"
                                 : (verbose_mode ? "verbose" : "normal"));

    // Get basic parameters
    // Check for mandatory positional argument
    if (!vm.count("guide-panman")) {
      err("Missing required pangenome file argument <guide.panman>");
      std::cerr << visible_options << std::endl; // Show help on error
      return 1;
    }
    std::string guide = vm["guide-panman"].as<std::string>();
    
    // Optional positional arguments
    std::string reads1 = vm.count("reads1") ? vm["reads1"].as<std::string>() : "";
    std::string reads2 = vm.count("reads2") ? vm["reads2"].as<std::string>() : "";
    
    // Get other parameters from vm (using defaults where applicable)
    std::string prefix = vm["prefix"].as<std::string>();
    std::string outputs = vm["outputs"].as<std::string>();
    std::string aligner = vm["aligner"].as<std::string>(); // Default handled by boost
    std::string refNode = vm["ref"].as<std::string>();     // Default handled by boost
    std::string eval = vm.count("eval") ? vm["eval"].as<std::string>() : ""; // Optional

  // Placement/diagnostic flags from CLI
  bool verify_scores_flag = vm.count("verify-scores") > 0;
  std::string debug_node_id_flag = vm.count("debug-node-id") ? vm["debug-node-id"].as<std::string>() : std::string();

    // Validate input file exists
    if (!fs::exists(guide)) {
      err("Pangenome file not found: {}", guide);
      return 1;
    }

    // Parse output options
    std::vector<std::string> outputs_seperated;
    std::stringstream ss(outputs);
    std::string token;
    while (std::getline(ss, token, ',')) {
      outputs_seperated.push_back(token);
    }

    // Initialize output filenames
    std::string refFileName, samFileName, bamFileName, mpileupFileName,
        vcfFileName;

    for (const auto &output : outputs_seperated) {
      if (output.size() == 1) {
        switch (output[0]) {
        case 'r':
          refFileName = prefix + ".reference.fa";
          break;
        case 's':
          samFileName = prefix + ".sam";
          break;
        case 'b':
          bamFileName = prefix + ".bam";
          break;
        case 'm':
          mpileupFileName = prefix + ".mpileup";
          break;
        case 'v':
          vcfFileName = prefix + ".vcf";
          break;
        case 'A':
          refFileName = prefix + ".reference.fa";
          samFileName = prefix + ".sam";
          bamFileName = prefix + ".bam";
          mpileupFileName = prefix + ".mpileup";
          vcfFileName = prefix + ".vcf";
          break;
        }
      } else {
        if (output == "reference")
          refFileName = prefix + ".reference.fa";
        else if (output == "sam")
          samFileName = prefix + ".sam";
        else if (output == "bam") {
          samFileName = prefix + ".sam";
          bamFileName = prefix + ".bam";
        } else if (output == "mpileup")
          mpileupFileName = prefix + ".mpileup";
        else if (output == "vcf")
          vcfFileName = prefix + ".vcf";
        else if (output == "all") {
          refFileName = prefix + ".reference.fa";
          samFileName = prefix + ".sam";
          bamFileName = prefix + ".bam";
          mpileupFileName = prefix + ".mpileup";
          vcfFileName = prefix + ".vcf";
        }
      }
    }

    std::mt19937 rng;
    if (vm.count("random-seed")) {
      std::string seed_str = vm["random-seed"].as<std::string>();
      std::hash<std::string> hasher;
      rng = std::mt19937(hasher(seed_str));
    } else {
      std::random_device rd;
      rng = std::mt19937(rd());
    }
    
    if (vm.count("dump-node-cluster") || vm.count("dump-node-cluster-leaves")
     || vm.count("dump-random-node-cluster") || vm.count("dump-random-node-cluster-leaves")
     || vm.count("dump-random-node-clusters") || vm.count("dump-random-node-clusters-leaves")
    ) {
    // Get timing flag early for use in all code paths
    bool show_time = vm.count("time") > 0;

    if (vm.count("mgsr-index")) {
      std::string mgsr_index_path = vm["mgsr-index"].as<std::string>();
      int fd = mgsr::open_file(mgsr_index_path);
      ::capnp::ReaderOptions readerOptions {.traversalLimitInWords = std::numeric_limits<uint64_t>::max(), .nestingLimit = 1024};
      ::capnp::PackedFdMessageReader reader(fd, readerOptions);
      MGSRIndex::Reader indexReader = reader.getRoot<MGSRIndex>();
      LiteTree::Reader liteTreeReader = indexReader.getLiteTree();
      size_t numThreads = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
      bool lowMemory = vm.count("low-memory") > 0;
      
      mgsr::MgsrLiteTree T;
      T.initialize(indexReader, numThreads, lowMemory, false);

      if (vm.count("dump-node-cluster") || vm.count("dump-node-cluster-leaves")) {
        if (!vm.count("mgsr-index")) {
          err("--mgsr-index must be provided to use --dump-node-cluster");
          return 1;
        }
        auto inputParameters = vm.count("dump-node-cluster") 
                            ? vm["dump-node-cluster"].as<std::vector<std::string>>()
                            : vm["dump-node-cluster-leaves"].as<std::vector<std::string>>();
        if (inputParameters.size() != 2) {
          err("Expected 2 parameters for --dump-node-cluster: <nodeID> <numNodes>");
          return 1;
        }

        uint32_t numNodes;
        try {
          numNodes = std::stoi(inputParameters[1]);
          // Proceed with using numNodes
        } catch (const std::invalid_argument& e) {
          err("The second parameter is not convertible to an integer");
          return 1;
        }

        std::string nodeId = inputParameters[0];
        if (T.allLiteNodes.find(nodeId) == T.allLiteNodes.end()) {
          err("Node ID {} not found in the tree", nodeId);
          return 1;
        }

        std::vector<mgsr::MgsrLiteNode*> nearestNodes = vm.count("dump-node-cluster")
                                                    ? mgsr::getNearestNodes(T.allLiteNodes.find(nodeId)->second, numNodes, false)
                                                    : mgsr::getNearestNodes(T.allLiteNodes.find(nodeId)->second, numNodes, true);
        std::ofstream outFile(prefix + ".clusterIDs.tsv");
        outFile << "Strain\tClusterID" << std::endl;
        for (const auto& node : nearestNodes) {
          outFile << node->identifier << "\t0" << std::endl;
        }
        outFile.close();
        msg("Cluster node IDs written to {}", prefix + ".clusterIDs.tsv");
        exit(0);
      } else if (vm.count("dump-random-node-cluster") || vm.count("dump-random-node-cluster-leaves")) {
        if (!vm.count("mgsr-index")) {
          err("--mgsr-index must be provided to use --dump-node-cluster");
          return 1;
        }
        uint32_t numNodes = vm.count("dump-random-node-cluster")
                            ? vm["dump-random-node-cluster"].as<uint32_t>()
                            : vm["dump-random-node-cluster-leaves"].as<uint32_t>();
        std::vector<std::string_view> allNodeIDs;
        allNodeIDs.reserve(T.allLiteNodes.size());
        for (const auto& [nodeID, node] : T.allLiteNodes) {
          if (node->children.empty()) {
            allNodeIDs.push_back(nodeID);
          }
        }
        allNodeIDs.shrink_to_fit();

        std::shuffle(allNodeIDs.begin(), allNodeIDs.end(), rng);

        std::vector<mgsr::MgsrLiteNode*> nearestNodes = vm.count("dump-random-node-cluster")
                                                    ? mgsr::getNearestNodes(T.allLiteNodes.find(std::string(allNodeIDs[0]))->second, numNodes, false)
                                                    : mgsr::getNearestNodes(T.allLiteNodes.find(std::string(allNodeIDs[0]))->second, numNodes, true);

        std::ofstream outFile(prefix + ".randomClusterIDs.tsv");
        outFile << "Strain\tClusterID" << std::endl;
        for (const auto& node : nearestNodes) {
          outFile << node->identifier << "\t0" << std::endl;
        }
        outFile.close();
        msg("Random node IDs written to {}", prefix + ".randomClusterIDs.tsv");
        exit(0);
      } else if (vm.count("dump-random-node-clusters") || vm.count("dump-random-node-clusters-leaves")) {
        if (!vm.count("mgsr-index")) {
          err("--mgsr-index must be provided to use --dump-node-cluster");
          return 1;
        }
        std::vector<uint32_t> nodeClusterNodes = vm.count("dump-random-node-clusters")
                                              ? vm["dump-random-node-clusters"].as<std::vector<uint32_t>>()
                                              : vm["dump-random-node-clusters-leaves"].as<std::vector<uint32_t>>();

        std::sort(nodeClusterNodes.begin(), nodeClusterNodes.end(), std::greater<uint32_t>());

        std::vector<std::string_view> allNodeIDs;
        allNodeIDs.reserve(T.allLiteNodes.size());
        for (const auto& [nodeID, node] : T.allLiteNodes) {
          if (node->children.empty()) {
            allNodeIDs.push_back(nodeID);
          }
        }
        allNodeIDs.shrink_to_fit();

        std::shuffle(allNodeIDs.begin(), allNodeIDs.end(), rng);

        std::vector<std::vector<std::string_view>> nodeClusters(nodeClusterNodes.size());
        std::unordered_set<std::string_view> selectedNodes;
        std::unordered_set<std::string_view> paddedNodes;
        size_t allNodeIDsIndex = 0;
        for (size_t i = 0; i < nodeClusterNodes.size(); i++) {
          uint32_t curClusterSize = nodeClusterNodes[i];
          uint32_t nextClusterSize = (i + 1 < nodeClusterNodes.size()) ? nodeClusterNodes[i + 1] : 1;

          std::string_view curNode;
          while (true) {
            curNode = allNodeIDs[allNodeIDsIndex];
            ++allNodeIDsIndex;
            if (allNodeIDsIndex >= allNodeIDs.size()) {
              err("Not enough nodes to fulfill the requested cluster sizes");
              exit(1);
            }
            if (paddedNodes.find(curNode) == paddedNodes.end()) {
              break;
            }
          }

          std::vector<mgsr::MgsrLiteNode*> nearestNodes = vm.count("dump-random-node-clusters")
              ? mgsr::getNearestNodes(T.allLiteNodes.find(std::string(curNode))->second, selectedNodes, curClusterSize + nextClusterSize - 1, false)
              : mgsr::getNearestNodes(T.allLiteNodes.find(std::string(curNode))->second, selectedNodes, curClusterSize + nextClusterSize - 1, true);
          for (size_t j = 0; j < nearestNodes.size(); j++) {
            if (j < curClusterSize) {
              selectedNodes.insert(nearestNodes[j]->identifier);
              nodeClusters[i].push_back(nearestNodes[j]->identifier);
            }
            paddedNodes.insert(nearestNodes[j]->identifier);
          }
        }

        std::ofstream outFile(prefix + ".randomClusterIDs.tsv");
        outFile << "Strain\tClusterID" << std::endl;
        for (size_t i = 0; i < nodeClusters.size(); i++) {
          for (const auto& node : nodeClusters[i]) {
            outFile << node << "\t" << (i + 1) << std::endl;
          }
        }
        outFile.close();
        msg("Random node IDs written to {}", prefix + ".randomClusterIDs.tsv");
        exit(0);
      }
    }



    if (vm.count("mgsr-index") && !reads1.empty()) {
      std::string mgsr_index_path = vm["mgsr-index"].as<std::string>();
      int fd = mgsr::open_file(mgsr_index_path);
      ::capnp::ReaderOptions readerOptions {.traversalLimitInWords = std::numeric_limits<uint64_t>::max(), .nestingLimit = 1024};
      ::capnp::PackedFdMessageReader reader(fd, readerOptions);
      MGSRIndex::Reader indexReader = reader.getRoot<MGSRIndex>();
      LiteTree::Reader liteTreeReader = indexReader.getLiteTree();
      size_t numThreads = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
      bool lowMemory = vm.count("low-memory") > 0;
      
      mgsr::MgsrLiteTree liteTree;
      if (vm.count("true-abundance")) {
        std::string trueAbundancePath = vm["true-abundance"].as<std::string>();
        liteTree.loadTrueAbundances(trueAbundancePath);
      }
      liteTree.initialize(indexReader, numThreads, lowMemory, true);
      liteTree.debugNodeID = vm["debug-node-id"].as<std::string>();
      std::vector<std::string> readSequences;
      mgsr::extractReadSequences(reads1, reads2, readSequences);

      uint32_t maskReads = vm.count("mask-reads") ? vm["mask-reads"].as<uint32_t>() : 0;
      uint32_t maskSeedThreshold = vm["mask-seeds"].as<uint32_t>();
      bool progressBar = true;
      if (vm.count("no-progress")) progressBar = false;
      mgsr::ThreadsManager threadsManager(&liteTree, numThreads, maskReads, progressBar, lowMemory);
      threadsManager.initializeMGSRIndex(indexReader);
      close(fd);
      threadsManager.initializeQueryData(readSequences, maskSeedThreshold);

      if (vm.count("overlap-coefficients")) {
        auto start_time_computeOverlapCoefficients = std::chrono::high_resolution_clock::now();
        mgsr::mgsrPlacer placer(&liteTree, threadsManager, lowMemory, 0);
        auto overlapCoefficients = placer.computeOverlapCoefficients(threadsManager.allSeedmerHashesSet);
        std::ofstream overlapCoefficientsFile(prefix + ".overlapCoefficients.txt");
        std::sort(overlapCoefficients.begin(), overlapCoefficients.end(), [](const auto& a, const auto& b) {
          return a.second > b.second;
        });
        uint32_t rank = 0;
        double currentOverlapCoefficient = overlapCoefficients[0].second;
        for (const auto& [nodeId, overlapCoefficient] : overlapCoefficients) {
          if (overlapCoefficient != currentOverlapCoefficient) {
            currentOverlapCoefficient = overlapCoefficient;
            ++rank;
          }
          overlapCoefficientsFile << nodeId << "\t" << std::fixed << std::setprecision(6) << overlapCoefficient << "\t" << rank <<  std::endl;
        }
        overlapCoefficientsFile.close();
        auto end_time_computeOverlapCoefficients = std::chrono::high_resolution_clock::now();
        auto duration_computeOverlapCoefficients = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_computeOverlapCoefficients - start_time_computeOverlapCoefficients);
        std::cout << "Computed overlap coefficients in " << static_cast<double>(duration_computeOverlapCoefficients.count()) / 1000.0 << "s\n" << std::endl;
        exit(0);
      }

      liteTree.collapseIdenticalScoringNodes(threadsManager.allSeedmerHashesSet);
      // liteTree.collapseEmptyNodes(true);

      if (vm.count("seed-scores")) {
        threadsManager.countSeedNodesFrequency();
        threadsManager.computeNodeSeedScores();
        std::ofstream nodeSeedScoresOut(prefix + ".nodeSeedScores.tsv");
        nodeSeedScoresOut << "NodeId\tDistance\tnodeSeedScores\tnodeSeedScoresCorrected\tselected\tselectedNeighbor\tinSample\tcollapsedNodes" << std::endl;
        for (const auto& [node, nodeSeedScore] : threadsManager.nodeSeedScores) {
          nodeSeedScoresOut << node->identifier
                            << "\t" << node->seedDistance
                            << "\t" << nodeSeedScore
                            << "\t" << threadsManager.nodeSeedScoresCorrected.find(node)->second
                            << "\t" << node->selected
                            << "\t" << node->selectedNeighbor
                            << "\t" << node->inSample << "\t";
          if (node->identicalNodeIdentifiers.empty()) {
            nodeSeedScoresOut << "." << std::endl;
          } else {
            for (size_t i = 0; i < node->identicalNodeIdentifiers.size(); ++i) {
              nodeSeedScoresOut << node->identicalNodeIdentifiers[i];
              if (i != node->identicalNodeIdentifiers.size() - 1) {
                nodeSeedScoresOut << ",";
              }
            }
            nodeSeedScoresOut << std::endl;
          }

        }
        nodeSeedScoresOut.close();
        exit(0);
      }
            
      
      panmapUtils::LiteTree liteTree;
      liteTree.initialize(liteTreeReader);

      std::vector<std::string> readSequences;
      mgsr::extractReadSequences(reads1, reads2, readSequences);

      bool skipSingleton = vm.count("skip-singleton") > 0;
      bool lowMemory = vm.count("low-memory") > 0;
      mgsr::ThreadsManager threadsManager(&liteTree, numThreads, skipSingleton, lowMemory);
      threadsManager.initializeMGSRIndex(indexReader);
      close(fd);
      threadsManager.initializeQueryData(readSequences);
      std::cout << "Total unique kminmers: " << threadsManager.allSeedmerHashesSet.size() << std::endl;
      
      std::vector<uint64_t> totalNodesPerThread(numThreads, 0);
      for (size_t i = 0; i < numThreads; ++i) {
        totalNodesPerThread[i] = liteTree.allLiteNodes.size();
      }
      ProgressTracker progressTracker(numThreads, totalNodesPerThread);

      std::cout << "Using " << numThreads << " threads" << std::endl;

      // mgsr::mgsrPlacer placer(&liteTree, threadsManager, lowMemory);
      // placer.setProgressTracker(&progressTracker, 0);
      // auto start_time_traverseTree = std::chrono::high_resolution_clock::now();
      // placer.traverseTree();
      // auto end_time_traverseTree = std::chrono::high_resolution_clock::now();
      // auto duration_traverseTree = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_traverseTree - start_time_traverseTree);
      // std::cout << "Traversed tree in " << std::fixed << std::setprecision(3) << static_cast<double>(duration_traverseTree.count()) / 1000.0 << "s\n" << std::endl;
      // exit(0);


      // auto start_time_computeOverlapCoefficients = std::chrono::high_resolution_clock::now();
      // mgsr::mgsrPlacer placer(&liteTree, threadsManager);
      // placer.computeOverlapCoefficients(threadsManager.allSeedmerHashesSet);
      // auto end_time_computeOverlapCoefficients = std::chrono::high_resolution_clock::now();
      // auto duration_computeOverlapCoefficients = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_computeOverlapCoefficients - start_time_computeOverlapCoefficients);
      // std::cout << "Computed overlap coefficients in " << static_cast<double>(duration_computeOverlapCoefficients.count()) / 1000.0 << "s\n" << std::endl;
      // exit(0);

      auto start_time_place = std::chrono::high_resolution_clock::now();
      std::atomic<size_t> numGroupsUpdate = 0;
      std::atomic<size_t> numReadsUpdate = 0;
      tbb::parallel_for(tbb::blocked_range<size_t>(0, threadsManager.threadRanges.size()), [&](const tbb::blocked_range<size_t>& rangeIndex){
        for (size_t i = rangeIndex.begin(); i != rangeIndex.end(); ++i) {
          auto [start, end] = threadsManager.threadRanges[i];
        
          std::span<mgsr::Read> curThreadReads(threadsManager.reads.data() + start, end - start);
          mgsr::mgsrPlacer curThreadPlacer(&liteTree, threadsManager, lowMemory);
          curThreadPlacer.initializeQueryData(curThreadReads);
          curThreadPlacer.setAllSeedmerHashesSet(threadsManager.allSeedmerHashesSet);
          curThreadPlacer.setShowTime(show_time);

          curThreadPlacer.setProgressTracker(&progressTracker, i);

          curThreadPlacer.placeReads();

          // Move score deltas from placer to thread manager
          threadsManager.perNodeScoreDeltasIndexByThreadId[i] = std::move(curThreadPlacer.perNodeScoreDeltasIndex);
          threadsManager.readMinichainsInitialized[i] = curThreadPlacer.readMinichainsInitialized;
          threadsManager.readMinichainsAdded[i] = curThreadPlacer.readMinichainsAdded;
          threadsManager.readMinichainsRemoved[i] = curThreadPlacer.readMinichainsRemoved;
          threadsManager.readMinichainsUpdated[i] = curThreadPlacer.readMinichainsUpdated;
          if (i == 0) {
            threadsManager.identicalGroups = std::move(curThreadPlacer.identicalGroups);
            threadsManager.identicalNodeToGroup = std::move(curThreadPlacer.identicalNodeToGroup);
          }
          // Note: numGroupsUpdate and numReadsUpdate removed from mgsrPlacer
        }
      });
      auto end_time_place = std::chrono::high_resolution_clock::now();
      auto duration_place = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_place - start_time_place);
      std::cerr << "\n\nPlaced reads in " << static_cast<double>(duration_place.count()) / 1000.0 << "s\n" << std::endl;

      // threadsManager.scoreNodes();
      threadsManager.scoreNodesMultithreaded();
      
      mgsr::mgsrPlacer placerOC(&liteTree, threadsManager, lowMemory, 0);
      auto overlapCoefficients = placerOC.computeOverlapCoefficients(threadsManager.allSeedmerHashesSet);
      std::unordered_map<std::string, double>().swap(threadsManager.kminmerOverlapCoefficients);
      for (const auto& [nodeId, overlapCoefficient] : overlapCoefficients) {
        threadsManager.kminmerOverlapCoefficients[nodeId] = overlapCoefficient;
      }

      std::string collapsedNewick = liteTree.toNewick(true);
      std::ofstream of(prefix + ".collapsed.newick");
      of << collapsedNewick;
      of.close();
      std::ofstream scoresOut(prefix + ".nodeScores.tsv");
      scoresOut << "NodeId\tdistance\tWEPPScore\tWEPPScoreCorrected\tWEPPScoreCorrectedSelected\tSelectedNeighbor\tinSample\tcollapsedNodes" << std::endl;
      for (const auto& [nodeId, node] : liteTree.allLiteNodes) {
        if (liteTree.detachedNodes.find(node) != liteTree.detachedNodes.end()) {
          continue;
        }
        scoresOut << nodeId
                  << "\t" << node->seedDistance
                  << "\t" << node->sumWEPPScore.sum
                  << "\t" << node->sumWEPPScoreCorrected.sum
                  << "\t" << node->selected
                  << "\t" << node->selectedNeighbor
                  << "\t" << node->inSample << "\t";
        if (node->identicalNodeIdentifiers.empty()) {
          scoresOut << "." << std::endl;
        } else {
          for (size_t i = 0; i < node->identicalNodeIdentifiers.size(); ++i) {
            scoresOut << node->identicalNodeIdentifiers[i];
            if (i != node->identicalNodeIdentifiers.size() - 1) {
              scoresOut << ",";
            }
          }
          scoresOut << std::endl;
        }
      }
      scoresOut.close();
      exit(0);

      mgsr::squareEM squareEM(threadsManager, liteTree, prefix, 1000);
      liteTree.cleanup(); // no longer needed. clear memory to prep for EM.
      
      auto start_time_squareEM = std::chrono::high_resolution_clock::now();
      for (size_t i = 0; i < 5; ++i) {
        squareEM.runSquareEM(1000);
        std::cout << "\nRound " << i << " of squareEM completed... nodes size changed from " << squareEM.nodes.size() << " to ";
        bool removed = squareEM.removeLowPropNodes();
        std::cout << squareEM.nodes.size() << std::endl;
        if (!removed) {
          break;
        }
      }
      
      std::vector<uint64_t> indices(squareEM.nodes.size());
      std::iota(indices.begin(), indices.end(), 0);
      std::sort(indices.begin(), indices.end(), [&squareEM](uint64_t i, uint64_t j) {
        return squareEM.props[i] > squareEM.props[j];
      });
      std::cout << std::endl;

      std::ofstream abundanceOutput(prefix + ".mgsr.abundance.out");
      abundanceOutput << std::setprecision(5) << std::fixed;
      for (size_t i = 0; i < indices.size(); ++i) {
        size_t index = indices[i];
        abundanceOutput << squareEM.nodes[index];
        if (squareEM.identicalGroups.find(squareEM.nodes[index]) != squareEM.identicalGroups.end()) {
          for (const auto& member : squareEM.identicalGroups[squareEM.nodes[index]]) {
            abundanceOutput << "," << member;
          }
        }
        abundanceOutput << "\t" << squareEM.props[index] << std::endl;
      }
      abundanceOutput.close();
      auto end_time_squareEM = std::chrono::high_resolution_clock::now();
      auto duration_squareEM = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_squareEM - start_time_squareEM);
      std::cout << "SquareEM completed in " << static_cast<double>(duration_squareEM.count()) / 1000.0 << "s\n" << std::endl;

      exit(0);
    }

    // Get remaining parameters
    bool reindex = vm.count("reindex") > 0;
    bool prior = vm.count("prior") > 0;
    bool genotype_from_sam = vm.count("genotype-from-sam") > 0;
    bool save_jaccard = vm.count("save-jaccard") > 0;
    // show_time already declared earlier for use in MGSR section
    bool dump_random_node = vm.count("dump-random-node") > 0;

    int k = vm["k"].as<int>(); // Default handled by boost
    int s = vm["s"].as<int>(); // Default handled by boost
    std::string index_path = vm["index"].as<std::string>(); // Default handled by boost


    // Load pangenome
    msg("Loading reference pangenome from: {}", guide);
    std::random_device rd;
    std::mt19937 rng(rd());

    // OPTIMIZATION: Defer loading the full panman Tree until actually needed.
    // LiteTree-based placement doesn't require the full Tree, so we only load it for:
    // - Index building (--index-mgsr or when index doesn't exist)
    // - Diagnostic operations (--dump-sequence, --dump-random-node)
    // - Post-placement operations (alignment, genotyping)
    // This significantly reduces memory usage and startup time for placement-only workflows.
    
    // Helper lambda to load Tree when needed
    std::unique_ptr<panmanUtils::TreeGroup> TG;
    panmanUtils::Tree* T_ptr = nullptr;
    
    auto ensureTreeLoaded = [&]() -> panmanUtils::Tree& {
      if (!TG) {
        auto time_pangenome_start = std::chrono::high_resolution_clock::now();
        msg("Loading reference pangenome from: {}", guide);
        
        try {
          TG = std::unique_ptr<panmanUtils::TreeGroup>(loadPanMAN(guide));
          
          if (!TG || TG->trees.empty()) {
            throw std::runtime_error("No valid trees found in the loaded PanMAN file.");
          }
          
          // Debug log for tree structure
          for (size_t i = 0; i < TG->trees.size(); i++) {
            panmanUtils::Tree &T = TG->trees[i];
            logging::debug(
                "[DEBUG] Tree {} details: blocks.size={}, gaps.size={}, "
                "blockGaps.blockPosition.size={}, allNodes.size={}, root={}",
                i, T.blocks.size(), T.gaps.size(), T.blockGaps.blockPosition.size(),
                T.allNodes.size(), T.root ? T.root->identifier : "null");
          }
          
          msg("Successfully loaded pangenome with {} trees:", TG->trees.size());
          T_ptr = &TG->trees[0];
          
          auto time_pangenome_end = std::chrono::high_resolution_clock::now();
          auto duration_pangenome = std::chrono::duration_cast<std::chrono::milliseconds>(
              time_pangenome_end - time_pangenome_start);
          if (show_time) {
            msg("[TIME] Pangenome loading: {}ms", duration_pangenome.count());
          }
        } catch (const std::exception &e) {
          err("Failed to load reference PanMAN: {}", e.what());
          throw;
        }
      }
      return *T_ptr;
    };

    if (vm.count("index-mgsr")) {
      panmanUtils::Tree &T = ensureTreeLoaded();
      std::string mgsr_index_path = vm["index-mgsr"].as<std::string>();
      int mgsr_t = 0;
      int mgsr_l = 3;
      bool open = false;
      mgsr::mgsrIndexBuilder mgsrIndexBuilder(&T, 19, 8, mgsr_t, mgsr_l, open);
      mgsrIndexBuilder.buildIndex();
      mgsrIndexBuilder.writeIndex(mgsr_index_path);
      msg("MGSR index written to: {}", mgsr_index_path);
      return 0;
    }

    // Handle --dump-random-node parameter if provided
    if (dump_random_node) {
      ensureTreeLoaded(); // Ensure TG is loaded
      panmanUtils::Node *randomNode = getRandomNode(TG.get(), rng);
      if (!randomNode) {
        err("Failed to select a random node");
        return 1;
      }

      // Find which tree the node belongs to
      panmanUtils::Tree *nodeTree = nullptr;
      for (auto &tree : TG->trees) {
        if (tree.allNodes.find(randomNode->identifier) != tree.allNodes.end()) {
          nodeTree = &tree;
          break;
        }
      }

      if (!nodeTree) {
        err("Could not determine which tree the random node belongs to");
        return 1;
      }

      // Sanitize node identifier for use in filename (replace / with _)
      std::string sanitizedNodeId = randomNode->identifier;
      std::replace(sanitizedNodeId.begin(), sanitizedNodeId.end(), '/', '_');
      
      std::string outputFileName =
          guide + ".random." + sanitizedNodeId + ".fa";
      if (saveNodeSequence(nodeTree, randomNode, outputFileName)) {
        msg("Random node {} sequence written to {}", randomNode->identifier,
            outputFileName);
        return 0; // Exit after dumping sequence
      } else {
        err("Failed to save random node sequence to {}", outputFileName);
        return 1;
      }
    }

    // Handle --dump-sequence parameter if provided
    if (vm.count("dump-sequence")) {
      panmanUtils::Tree &T = ensureTreeLoaded();
      std::string nodeID = vm["dump-sequence"].as<std::string>();
      if (T.allNodes.find(nodeID) == T.allNodes.end()) {
        err("Node ID {} not found in the tree", nodeID);
        return 1;
      }

      std::string sequence = T.getStringFromReference(nodeID, false, true);
      std::string outputFileName = guide + "." + nodeID + ".fa";
      std::ofstream outFile(outputFileName);

      if (outFile.is_open()) {
        outFile << ">" << nodeID << "\n";
        outFile << sequence << "\n";
        outFile.close();
        msg("Sequence for node {} written to {}", nodeID, outputFileName);
        return 0; // Exit after dumping sequence
      } else {
        err("Failed to open file {} for writing", outputFileName);
        return 1;
      }
    }

    // Log settings
    msg("--- Settings ---");
    msg("Reads: {}",
        (reads1.empty() ? "<none>"
                        : reads1 + (reads2.empty() ? "" : " + " + reads2)));
    msg("Reference PanMAN: {}", guide);
    msg("Using {} threads", vm["cpus"].as<int>());

    // Handle indexing
    bool build = true;
    ::capnp::MallocMessageBuilder outMessage;
    std::unique_ptr<::capnp::MessageReader> inMessage;
    
    // Ensure paths are consistent by normalizing to absolute paths
    boost::filesystem::path guidePath = boost::filesystem::absolute(guide);
    std::string default_index_path = guidePath.string() + ".pmi";
    
    // If an explicit index path was provided, also normalize it
    std::string normalized_index_path = "";
    if (!index_path.empty()) {
      boost::filesystem::path explicitIndexPath = boost::filesystem::absolute(index_path);
      normalized_index_path = explicitIndexPath.string();
    }
    
    // Use the normalized paths for all operations
    std::string effective_index_path = normalized_index_path.empty() ? default_index_path : normalized_index_path;
    
    logging::debug("DEBUG-INDEX: Default index path: {}", default_index_path);
    logging::debug("DEBUG-INDEX: Explicit index path: {}", normalized_index_path);
    logging::debug("DEBUG-INDEX: Effective index path: {}", effective_index_path);
    logging::debug("DEBUG-INDEX: Reindex flag: {}", reindex ? "true" : "false");

    // CRITICAL FIX: Add checkpoint to check for and remove corrupted index
    if (fs::exists(effective_index_path) && !reindex) {
      msg("Checking if existing index {} is valid...", effective_index_path);
      // Try to read the index file to validate it
      try {
        logging::debug("DEBUG-INDEX: Attempting to load existing index from effective path");
        inMessage = readCapnp(effective_index_path);
        
        // If we get here, the index loaded successfully
        build = false;
        msg("Successfully loaded existing index from: {}", effective_index_path);
        logging::debug("DEBUG-INDEX: Successfully loaded index from effective path");
      } catch (const std::exception &e) {
        err("Existing index appears to be corrupted: {}", e.what());
        logging::debug("DEBUG-INDEX: Failed to load existing index: {}", e.what());
        msg("Will rebuild the index. Use -f/--reindex to skip this check.");
        logging::debug("DEBUG-INDEX: Removing corrupted index file: {}", effective_index_path);
        fs::remove(effective_index_path);
        build = true;
      }
    } 
    else {
      msg(reindex ? "Will re-build index (-f)"
                  : "No index found, will build new index");
      logging::debug("DEBUG-INDEX: Will build new index. Reason: {}", 
                  reindex ? "Reindex flag set" : "No existing index found");
    }

    int mgsr_t = 0;
    int mgsr_l = 0;
    bool open = false;
    bool use_raw_seeds = true; // true for panmap, false for panmama
    
    // Only load tree if we need to build the index
    if (build) {
      panmanUtils::Tree &T = ensureTreeLoaded();
      auto time_index_build_start = std::chrono::high_resolution_clock::now();
      mgsr::mgsrIndexBuilder mgsrIndexBuilder(&T, k, s, mgsr_t, mgsr_l, open, use_raw_seeds);
      mgsrIndexBuilder.buildIndex();
      auto time_index_build_end = std::chrono::high_resolution_clock::now();
      auto duration_index_build = std::chrono::duration_cast<std::chrono::milliseconds>(
          time_index_build_end - time_index_build_start);
      
      if (show_time) {
        msg("[TIME] Index building: {}ms", duration_index_build.count());
      }
      
      // Log unique seed/k-mer count
      if (use_raw_seeds) {
        msg("Built raw seeds index with {} unique seeds", mgsrIndexBuilder.uniqueSyncmers.size());
      } else {
        msg("Built k-minmers index with {} unique k-minmers", mgsrIndexBuilder.uniqueKminmers.size());
      }
      
      // Use the effective output path for MGSR index
      auto time_index_write_start = std::chrono::high_resolution_clock::now();
      mgsrIndexBuilder.writeIndex(effective_index_path);
      auto time_index_write_end = std::chrono::high_resolution_clock::now();
      auto duration_index_write = std::chrono::duration_cast<std::chrono::milliseconds>(
          time_index_write_end - time_index_write_start);
      
      if (show_time) {
               msg("[TIME] Index writing: {}ms", duration_index_write.count());
      }
    }

    // Generate mutation matrices during indexing if requested
    std::string mm_path = guide + ".mm";
    bool build_mm = vm.count("with-mm") > 0;
    
    if (!build_mm) {
      msg("Skipping mutation matrix generation (use --with-mm to enable)");
    } else if (reindex || !fs::exists(mm_path)) {
      msg("=== Generating Mutation Matrices ===");
      panmanUtils::Tree &T = ensureTreeLoaded();
      auto time_mutmat_start = std::chrono::high_resolution_clock::now();
      genotyping::mutationMatrices mutMat;
      msg("Building mutation matrices from tree");
      genotyping::fillMutationMatricesFromTree_test(mutMat, &T, mm_path);
      auto time_mutmat_end = std::chrono::high_resolution_clock::now();
      auto duration_mutmat = std::chrono::duration_cast<std::chrono::milliseconds>(
          time_mutmat_end - time_mutmat_start);
      
      msg("Mutation matrices written to: {}", mm_path);
      if (show_time) {
        msg("[TIME] Mutation matrix generation: {}ms", duration_mutmat.count());
      }
    }

    // Exit if no reads were provided - indexing is complete
    if (reads1.empty()) {
      msg("Index building complete. No reads provided, skipping placement.");
      
      auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
          std::chrono::high_resolution_clock::now() - main_start);
      if (show_time) {
        msg("=== TIMING SUMMARY ===");
        msg("Total runtime: {}ms", total_duration.count());
      }
      
      return 0;
    }

    // Prepare MGSRIndex reader for placement
    int fd_mgsr = ::open(effective_index_path.c_str(), O_RDONLY);
    if (fd_mgsr < 0) {
      throw std::runtime_error("Failed to open MGSR index file: " + effective_index_path);
    }
    ::capnp::ReaderOptions opts; opts.traversalLimitInWords = std::numeric_limits<uint64_t>::max(); opts.nestingLimit = 1024;
    ::capnp::PackedFdMessageReader mgsrMsg(fd_mgsr, opts);
    ::MGSRIndex::Reader mgsrIndexRoot = mgsrMsg.getRoot<MGSRIndex>();

    // Log index parameters
    msg("Index parameters: k={}, s={}, t={}, l={}, open={}", 
        mgsrIndexRoot.getK(), mgsrIndexRoot.getS(), mgsrIndexRoot.getT(), 
        mgsrIndexRoot.getL(), mgsrIndexRoot.getOpen());

    // Run placement using MGSR index format
    // Control placement method via command line options
    bool use_lite_placement = vm.count("use-lite-placement") > 0;
    bool use_full_placement = vm.count("use-full-placement") > 0;
    
    // Default to LiteTree placement if neither option is specified (more efficient)
    if (!use_lite_placement && !use_full_placement) {
        use_lite_placement = true;
    }
    
    // If both options are specified, prefer full placement for backward compatibility
    if (use_lite_placement && use_full_placement) {
        logging::warn("Both --use-lite-placement and --use-full-placement specified. Using full placement.");
        use_lite_placement = false;
    }
    
    placement::PlacementResult result;
    std::vector<std::vector<seeding::seed_t>> readSeeds;
    std::vector<std::string> readSequences;
    std::vector<std::string> readNames;
    std::vector<std::string> readQuals;
    std::string placementFileName = guide + ".placement.tsv";
    
    if (use_lite_placement) {
        msg("=== Using LiteTree-based placement (efficient, no full tree loading) ===");
        
        // Initialize LiteTree from MGSR index
        panmapUtils::LiteTree liteTree;
        auto liteTreeReader = mgsrIndexRoot.getLiteTree();
        liteTree.initialize(liteTreeReader);
        
        msg("Initialized LiteTree with {} nodes and {} block ranges", 
            liteTree.allLiteNodes.size(), liteTree.blockScalarRanges.size());
        
        // Load full tree ONLY if verification mode is enabled
        panmanUtils::Tree* treeForVerification = nullptr;
        if (verify_scores_flag) {
            msg("Verification mode enabled - loading full Tree for getStringFromReference()");
            panmanUtils::Tree &T = ensureTreeLoaded();
            treeForVerification = &T;
        }
        
        // Use LiteTree-based placement
        placement::placeLite(result, &liteTree, mgsrIndexRoot, reads1, reads2,
                readSeeds, readSequences, readNames, readQuals,
                placementFileName, effective_index_path, debug_node_id_flag, 
                verify_scores_flag, treeForVerification);
        
        msg("LiteTree-based placement completed successfully");
        
    } else {
        msg("=== Using LiteTree-based placement (fallback) ===");
        
        // Initialize LiteTree from MGSR index since original place function was removed
        panmapUtils::LiteTree liteTree;
        auto liteTreeReader = mgsrIndexRoot.getLiteTree();
        liteTree.initialize(liteTreeReader);
        
        // Load full tree ONLY if verification mode is enabled
        panmanUtils::Tree* treeForVerification = nullptr;
        if (verify_scores_flag) {
            msg("Verification mode enabled - loading full Tree for getStringFromReference()");
            panmanUtils::Tree &T = ensureTreeLoaded();
            treeForVerification = &T;
        }
        
        // Use LiteTree-based placement as fallback
        placement::placeLite(result, &liteTree, mgsrIndexRoot, reads1, reads2,
                readSeeds, readSequences, readNames, readQuals,
                placementFileName, effective_index_path, debug_node_id_flag, 
                verify_scores_flag, treeForVerification);
    }
    
    return 0;

    // Build index if needed
    if (build) {
      logging::debug("DEBUG-INDEX: Starting index build process");
      try {
        panmanUtils::Tree &T = ensureTreeLoaded();
        Index::Builder index = outMessage.initRoot<Index>();
        logging::debug("DEBUG-INDEX: Initialized root for new index");
        
        if (k <= 0 || s <= 0 || s >= k) {
          logging::err("DEBUG-INDEX: Invalid k-mer or s-mer parameters: k={}, s={}", k, s);
          throw std::runtime_error("Invalid k-mer or s-mer parameters");
        }
        
        logging::debug("DEBUG-INDEX: Building index with k={}, s={}", k, s);
        auto start_indexing = std::chrono::high_resolution_clock::now(); // Renamed variable
        // This function now only builds the index in outMessage, it does not write it.
        indexing::index(&T, index, k, s, outMessage, effective_index_path); 
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start_indexing); // Use renamed variable
        
        logging::debug("DEBUG-INDEX: Index built in memory in {}ms", duration.count());
        if (show_time) {
          msg("[TIME] Index building: {}ms", duration.count());
        }

        // The index function now writes the index file directly when given a path
        // No need to call writeCapnp here anymore
        msg("Index written to: {}", effective_index_path);
        
        // Try to immediately read back the index to confirm it's valid
        logging::debug("DEBUG-INDEX: Verifying newly written index");
        try {
          // Add a small delay to ensure file system operations complete
          std::this_thread::sleep_for(std::chrono::milliseconds(500));
          inMessage = readCapnp(effective_index_path);
          if (inMessage) {
            auto root = inMessage->getRoot<Index>();
            bool indexValid = (root.getK() == k && root.getS() == s);
            logging::debug("DEBUG-INDEX: Verification result: index is {}valid", indexValid ? "" : "in");
          }
        } catch (const std::exception &e) {
          err("WARNING: Failed to verify newly written index: {}", e.what());
          logging::debug("DEBUG-INDEX: Verification failed: {}", e.what());
        }
       
      } catch (const std::exception &e) {
        err("ERROR during indexing or writing: {}", e.what());
        logging::debug("DEBUG-INDEX: Exception during index build/write: {}", e.what());
        return 1;
      }
    }

    // Load mutation matrices
    msg("=== Loading Mutation Matrices ===");
    genotyping::mutationMatrices mutMat;
    // Assuming --mutmat isn't needed with boost (wasn't in original docopt?)
    // If needed, add: ("mutmat", po::value<std::string>(), "Path to mutation matrix") to input_opts
    std::string mutmat_path = ""; // Replace with vm access if --mutmat option is added
    std::string default_mutmat_path = guide + ".mm";

    if (!mutmat_path.empty()) {
      msg("Loading mutation matrices from: {}", mutmat_path);
      std::ifstream mutmat_file(mutmat_path);
      genotyping::fillMutationMatricesFromFile(mutMat, mutmat_file);
    } else if (fs::exists(default_mutmat_path)) {
      msg("Loading default mutation matrices from: {}", default_mutmat_path);
      std::ifstream mutmat_file(default_mutmat_path);
      genotyping::fillMutationMatricesFromFile(mutMat, mutmat_file);
    } else {
      panmanUtils::Tree &T = ensureTreeLoaded();
      msg("Building new mutation matrices");
      genotyping::fillMutationMatricesFromTree_test(mutMat, &T, default_mutmat_path);
    }

    if (!eval.empty()) {
      msg("[Developer mode] --- Evaluate placement accuracy ---");
      std::exit(0);
    }

    // Placement phase
    msg("=== Starting Read Placement ===");
    auto start = std::chrono::high_resolution_clock::now();

    try {
      // Ensure index is loaded
      msg("Reading index...");
      auto time_index_load_start = std::chrono::high_resolution_clock::now();
      logging::debug("DEBUG-PLACE: inMessage is {}", inMessage ? "valid" : "nullptr");
      
      if (!inMessage) {
        try {
          logging::warn("Index was not successfully loaded, attempting to use existing index...");
          logging::debug("DEBUG-PLACE: Looking for index at {}", effective_index_path);
          // Only try to load from the default index path
          if (fs::exists(effective_index_path)) {
            logging::info("Trying to load index from: {}", effective_index_path);
            logging::debug("DEBUG-PLACE: Default index exists with size: {} bytes", fs::file_size(effective_index_path));
            
            // Check file permissions
            bool isReadable = access(effective_index_path.c_str(), R_OK) == 0;
            uintmax_t fileSize = fs::file_size(effective_index_path);
            if (!isReadable) {
              logging::err("DEBUG-PLACE: Index file exists but cannot be read (permission denied): {}", effective_index_path);
              throw std::runtime_error("Index file permission denied: " + effective_index_path);
            }
            
            if (fileSize == 0) {
              logging::err("DEBUG-PLACE: Index file exists but is empty: {}", effective_index_path);
              throw std::runtime_error("Empty index file: " + effective_index_path);
            }
            
            // Print detailed file info
            logging::debug("DEBUG-PLACE: Index file size: {} bytes", fileSize);
            
            // Try to load the index with detailed error reporting
            try {
              logging::debug("DEBUG-PLACE: Calling readCapnp on default index path");
              inMessage = readCapnp(effective_index_path);
              logging::debug("DEBUG-PLACE: Successfully loaded index from default path");
              logging::info("Successfully loaded index from: {}", effective_index_path);
            } catch (const std::exception& e) {
              logging::err("DEBUG-PLACE: Error loading index: {}", e.what());
              throw;
            }
          }
          // If index still not loaded and reindex isn't set, error out
          else if (!vm.count("reindex")) {
            logging::debug("DEBUG-PLACE: Default index path doesn't exist and reindex not set");
            throw std::runtime_error("Failed to load index from " + effective_index_path + 
                                   " and reindex option not specified");
          } else {
            logging::debug("DEBUG-PLACE: Default index doesn't exist but reindex flag set - will rebuild");
          }
        } catch (const std::exception& e) {
          err("ERROR loading index: {}. Run with -f/--reindex to rebuild the index or specify a working index with -i.", e.what());
          logging::debug("DEBUG-PLACE: Index loading failed: {}", e.what());
          return 1;
        }
      }
      
      // CRITICAL FIX: We must ensure inMessage is valid throughout the index_input usage
      // Get a reference to the inMessage as it's needed for index_input to stay valid
      if (!inMessage) {
        logging::err("DEBUG-PLACE: inMessage is still null after loading attempts");
        throw std::runtime_error("Critical error: Index message is null after loading");
      }
      
      auto time_index_load_end = std::chrono::high_resolution_clock::now();
      auto duration_index_load = std::chrono::duration_cast<std::chrono::milliseconds>(
          time_index_load_end - time_index_load_start);
      if (show_time) {
        msg("[TIME] Index loading: {}ms", duration_index_load.count());
      }
      
      logging::debug("DEBUG-PLACE: inMessage is valid, proceeding to extract root");
      
      // Scope to ensure inMessage stays alive as long as index_input is used
      {
        // Extract the root from the message - this reader depends on inMessage staying alive!
        Index::Reader index_input = inMessage->getRoot<Index>();
        logging::debug("DEBUG-PLACE: Successfully extracted root from inMessage");
        
        // Let's check for critical fields
        uint32_t k = index_input.getK();
        uint32_t s = index_input.getS();
        bool open = index_input.getOpen();
        uint32_t t = index_input.getT();
        size_t nodeCount = index_input.getPerNodeSeedMutations().size();
        logging::debug("DEBUG-PLACE: Index contains k={}, nodes={}", k, nodeCount);
        
        if (genotype_from_sam) {
          msg("Genotyping from SAM file");
          // Get required args for this mode
          if (!vm.count("sam-file") || !vm.count("ref-file")) {
              throw std::runtime_error(
                "--sam-file and --ref-file are required for --genotype-from-sam");
          }
          samFileName = vm["sam-file"].as<std::string>();
          refFileName = vm["ref-file"].as<std::string>();
          // Note: bamFileName, mpileupFileName, vcfFileName might need prefix logic here
          // depending on whether --prefix applies to this mode.
          // Assuming the existing output logic handles prefixing correctly.
          genotyping::genotype(prefix, refFileName, "", bamFileName,
                              mpileupFileName, vcfFileName, mutMat);
        } else {
          if (!refNode.empty()) {
            panmanUtils::Tree &T = ensureTreeLoaded();
            msg("Using reference node: {}", refNode);
            if (T.allNodes.find(refNode) == T.allNodes.end()) {
              throw std::runtime_error("Reference node not found in pangenome");
            }
          }
          std::string debug_specific_node_id = "";
          try {
            debug_specific_node_id = vm["debug-node-id"].as<std::string>();
          } catch (const std::exception &e) {
            debug_specific_node_id = "";
          }

          // Initialize placement variables and read information needed for alignment
          placement::PlacementResult result;
          std::vector<std::vector<seeding::seed_t>> readSeeds;
          std::vector<std::string> readSequences;
          std::vector<std::string> readNames;
          std::vector<std::string> readQuals;

          std::string placementFileName = prefix + "/placements.tsv";
          
          // Initialize LiteTree from MGSR index since original place function was removed
          auto time_placement_start = std::chrono::high_resolution_clock::now();
          panmapUtils::LiteTree liteTree;
          auto liteTreeReader = mgsrIndexRoot.getLiteTree();
          liteTree.initialize(liteTreeReader);
          
          // Perform placement using LiteTree-based approach
          placement::placeLite(result, &liteTree, mgsrIndexRoot, reads1, reads2,
                             readSeeds, readSequences, readNames, readQuals,
                             placementFileName, effective_index_path, debug_specific_node_id, verify_scores_flag);

          auto time_placement_end = std::chrono::high_resolution_clock::now();
          auto duration_placement = std::chrono::duration_cast<std::chrono::milliseconds>(
              time_placement_end - time_placement_start);
          if (show_time) {
            msg("[TIME] Placement computation: {}ms", duration_placement.count());
          }


          /* @Alan this is the end of placement
            -> next step is pass seeds of target node to Nico's alignment code
            -> and genotyping

            placementResult should have nodeSeedMap[targetId] with the target seed set
            
            TODO: I'm not sure k-mer end positions are correct yet
          */
          
          // Now load the tree for post-placement operations (alignment, genotyping)
          panmanUtils::Tree &T = ensureTreeLoaded();
          std::string bestMatchSequence = panmapUtils::getStringFromReference(&T, result.bestJaccardPresenceNodeId, false);
          logging::info("Best match sequence built for {} with length {}", result.bestJaccardPresenceNodeId, bestMatchSequence.size());


          
          std::vector<std::tuple<size_t, bool, bool, int64_t>> refSyncmers;
          std::unordered_map<size_t, std::pair<std::vector<uint32_t>, std::vector<uint32_t>>> seedToRefPositions;
          bool shortenSyncmers = false;
          
          if (k > 28 && !shortenSyncmers) {
            logging::warn("k > 28, setting k = 19, s = 10, t = 0, open = {} for minimap alignment", open);
            int k_minimap = 28;
            int s_minimap = 15;
            bool open_minimap = open;
            int t_minimap = 0;
            seeding::recalculateReadSeeds(k_minimap, s_minimap, open_minimap, t_minimap, readSequences, readSeeds);
            refSyncmers = seeding::rollingSyncmers(bestMatchSequence, k_minimap, s_minimap, open_minimap, t_minimap, false);
          } else {
            refSyncmers = seeding::rollingSyncmers(bestMatchSequence, k, s, open, t, false);
          }

          // going to build ref seed from scratch for now until Alex corrects k-mer end positions.
          for (const auto &[kmerHash, isReverse, isSyncmer, startPos] : refSyncmers) {
            if (!isSyncmer) {
              continue;
            }
            if (seedToRefPositions.find(kmerHash) == seedToRefPositions.end()) {
              seedToRefPositions[kmerHash] = std::make_pair(std::vector<uint32_t>(), std::vector<uint32_t>());
            }
            if (isReverse) {
              seedToRefPositions[kmerHash].second.push_back(startPos);
            } else {
              seedToRefPositions[kmerHash].first.push_back(startPos);
            }
          }

          bool pairedEndReads = reads1.size() > 0 && reads2.size() > 0;
          std::vector<char *> samAlignments;
          std::string samHeader;
          auto time_alignment_start = std::chrono::high_resolution_clock::now();
          if (k > 28) {
            if (shortenSyncmers) {
              createSam(readSeeds, readSequences, readQuals, readNames, bestMatchSequence, seedToRefPositions, samFileName, 28, shortenSyncmers, pairedEndReads, samAlignments, samHeader);
            } else {
              createSam(readSeeds, readSequences, readQuals, readNames, bestMatchSequence, seedToRefPositions, samFileName, 19, shortenSyncmers, pairedEndReads, samAlignments, samHeader);
            }
          } else {
            createSam(readSeeds, readSequences, readQuals, readNames, bestMatchSequence, seedToRefPositions, samFileName, k, shortenSyncmers, pairedEndReads, samAlignments, samHeader);
          }

          sam_hdr_t *header;
          bam1_t **bamRecords;
          createBam(samAlignments, samHeader, bamFileName, header, bamRecords);

          auto time_alignment_end = std::chrono::high_resolution_clock::now();
          auto duration_alignment = std::chrono::duration_cast<std::chrono::milliseconds>(
              time_alignment_end - time_alignment_start);
          if (show_time) {
            msg("[TIME] SAM/BAM creation: {}ms", duration_alignment.count());
          }

          auto time_genotyping_start = std::chrono::high_resolution_clock::now();
          createMplpBcf(prefix, refFileName, bestMatchSequence, bamFileName, mpileupFileName);

          createVcfWithMutationMatrices(prefix, mpileupFileName, mutMat, vcfFileName, 0.0011);

          auto time_genotyping_end = std::chrono::high_resolution_clock::now();
          auto duration_genotyping = std::chrono::duration_cast<std::chrono::milliseconds>(
              time_genotyping_end - time_genotyping_start);
          if (show_time) {
            msg("[TIME] Genotyping/VCF creation: {}ms", duration_genotyping.count());
          }

          std::string placementSummaryFileName = prefix + ".placement.summary.md";
          placement::dumpPlacementSummary(result, placementSummaryFileName);

          // Report the top 3 metrics as requested
          msg("=== TOP PLACEMENT RESULTS ===");
          
          // 1. Raw number of matches (set, no frequency) - using raw count
          msg("Top by raw seed matches (unique count): node {} with {} matches (Jaccard score {:.4f})",
              result.bestJaccardPresenceNodeId.empty() ? "none" : result.bestJaccardPresenceNodeId, 
              result.bestJaccardPresenceScore,
              result.bestJaccardPresenceScore);
          
          // 2. Number of matches scaled by read frequency 
          msg("Top by weighted seed matches (frequency-scaled): node {} with score {}",
              result.bestRawSeedMatchNodeId.empty() ? "none" : result.bestRawSeedMatchNodeId, 
              result.bestRawSeedMatchScore);
          
          // 3. Cosine similarity
          msg("Top by cosine similarity: node {} with score {:.4f}",
              result.bestCosineNodeId.empty() ? "none" : result.bestCosineNodeId, 
              result.bestCosineScore);
          
          msg("=== ADDITIONAL METRICS ===");
          msg("Best Weighted Jaccard: node {} with score {:.4f}",
              result.bestJaccardNodeId.empty() ? "none" : result.bestJaccardNodeId, result.bestJaccardScore); // bestJaccardScore is Weighted Jaccard
          msg("Best Cosine Similarity: node {} with score {:.4f}",
              result.bestCosineNodeId.empty() ? "none" : result.bestCosineNodeId, result.bestCosineScore);
          msg("Overall Best Weighted Score (Jaccard*scale + Cosine*(1-scale)): node {} with score {:.4f}",
              result.bestWeightedNodeId.empty() ? "none" : result.bestWeightedNodeId, result.bestWeightedScore);
        }
        
        // Validate the loaded index
        try {
          if (!inMessage) {
            throw std::runtime_error("Index could not be loaded");
          }
          
          // Try to access the index root to verify it's valid
          auto indexRoot = inMessage->getRoot<Index>();
          uint32_t k = indexRoot.getK();
          size_t nodeCount = indexRoot.getPerNodeSeedMutations().size();
          
          if (k == 0 || nodeCount == 0) {
            throw std::runtime_error("Invalid index: k=" + std::to_string(k) + 
                ", nodes=" + std::to_string(nodeCount) + 
                ". Index appears to be corrupted or incomplete.");
          }
          
          logging::info("Successfully validated index with {} nodes, k={}, s={}",
                      nodeCount, k, indexRoot.getS());
        } catch (const std::exception& e) {
          err("ERROR validating index: {}", e.what());
          if (vm.count("reindex")) {
            msg("Attempting to rebuild index due to --reindex flag being set...");
          } else {
            err("Index validation failed. Run with -f/--reindex to rebuild the index.");
            return 1;
          }
        }
      } // End of scope - index_input is no longer used after this point
      
    } catch (const std::exception &e) {
      err("ERROR during placement: {}", e.what());
      return 1;
    }

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start);
    msg("Placement completed in {}ms", duration.count());

    auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - main_start);
    
    if (show_time) {
      msg("=== TIMING SUMMARY ===");
    }
    msg("=== panmap run completed ===");
    msg("Total runtime: {}ms", total_duration.count());

    
    return 0;
}
