#include <boost/iostreams/categories.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
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

#include "capnp/message.h"
#include "capnp/serialize-packed.h"

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
      
      // Validate the root structure
      auto root = wrapper->getRoot<Index>();
      
      // Validate essential fields
      uint32_t k = root.getK();
      uint32_t s = root.getS();
      size_t nodeCount = root.getPerNodeSeedMutations().size();
        
      if (k == 0 || s == 0 || nodeCount == 0) {
        // Don't close fd here - will be closed by wrapper destructor
        throw std::runtime_error("Invalid index data: k=" + std::to_string(k) + 
                              ", s=" + std::to_string(s) + 
              ", nodes=" + std::to_string(nodeCount));
        }
        
      // Log additional validation checks
      logging::debug("Validated index structure: k={}, s={}, nodes={}, gapMutations={}, dictionary={}",
                    k, s, nodeCount, root.getPerNodeGapMutations().size(), root.getKmerDictionary().size());
      
      // Check optional fields if they're expected to be set
      if (root.hasNodePathInfo() && root.getNodePathInfo().size() == 0) {
        logging::warn("Index has empty nodePathInfo list");
      }
      
      if (root.hasBlockInfo() && root.getBlockInfo().size() == 0) {
        logging::warn("Index has empty blockInfo list");
      }
      
      logging::info("Successfully read and validated index from {}", resolvedPath);
      
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
        
        // Minimal cleanup - just try to shutdown spdlog quietly
        try {
            spdlog::set_level(spdlog::level::off);
            spdlog::shutdown();
        } catch (...) {
            // Ignore any exceptions during shutdown
        }
        
        const char* cleanup_msg = "Attempting graceful exit to preserve profiling data...\n";
        write(STDERR_FILENO, cleanup_msg, strlen(cleanup_msg));
        
        // Use exit() (not _exit()) to ensure atexit handlers run, including gprof's
        exit(0);  // Exit with 0 for gprof compatibility

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
 * @brief Load a PanMAN tree from a file
 *
 * @param filename Path to the PanMAN file
 * @return panmanUtils::TreeGroup* Loaded TreeGroup or nullptr if loading failed
 */
panmanUtils::TreeGroup *loadPanMAN(const std::string &filename) {
  try {
    // Open the input file
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
      logging::err("Failed to open PanMAN file: {}", filename);
      return nullptr;
    }

    // Set up filtering stream for decompression
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inBuffer;
    inBuffer.push(boost::iostreams::lzma_decompressor());
    inBuffer.push(inputFile);

    // Create the input stream that will be passed to the TreeGroup constructor
    std::istream inputStream(&inBuffer);

    // Create the TreeGroup (default isOld=false for newer format)
    panmanUtils::TreeGroup *TG = new panmanUtils::TreeGroup(inputStream);

    inputFile.close();
    logging::info("Successfully loaded PanMAN file: {}", filename);
    return TG;
  } catch (const std::exception &e) {
    logging::err("Exception while loading PanMAN file: {}", e.what());
    return nullptr;
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
        ("batch,b", po::value<std::string>(), "Path to batch TSV file")
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
        
    ;

    po::options_description mgsr_opts("MGSR options");
    mgsr_opts.add_options()
        ("index-mgsr", po::value<std::string>(), "Path to build/rebuild MGSR index")
        ("mgsr-index,m", po::value<std::string>(), "Path to precomputed MGSR index")
        ("l,l", po::value<int>()->default_value(1), "Length of k-min-mers (i.e. l seeds per kminmer)")
        ("skip-singleton", "Skip singleton reads")
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
        ("random-seed", po::value<std::string>(), "Seed for rng (read in as string then hashed). If not provided, a random seed will be used.")
        ("dump-sequence", po::value<std::string>(), "Dump sequence for a specific node ID")
        ("dump-sequences",po::value<std::vector<std::string>>()->multitoken(), "Dump sequences for a list of node IDs")
        ("simulate-snps",po::value<std::vector<uint32_t>>()->multitoken(), "Simulate number of SNPs for node IDs, parameter position is relative to dump-sequences")
        ("dump-random-nodeIDs", po::value<uint32_t>(), "Dump specified number of random node IDs from the tree")
        ("dump-random-node", "Dump sequence for a random node")
        ("debug-node-id", po::value<std::string>()->default_value(""), "Log detailed placement debug info for this specific node ID")
        ("candidate-threshold", po::value<float>()->default_value(0.01f), "Placement candidate threshold proportion") 
        ("max-candidates", po::value<int>()->default_value(16), "Maximum placement candidates")
        ("overlap-coefficients", "Output overlap coefficients then exit")
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

    if (vm.count("mgsr-index")) {
      std::string mgsr_index_path = vm["mgsr-index"].as<std::string>();
      int fd = mgsr::open_file(mgsr_index_path);
      ::capnp::ReaderOptions readerOptions {.traversalLimitInWords = std::numeric_limits<uint64_t>::max(), .nestingLimit = 1024};
      ::capnp::PackedFdMessageReader reader(fd, readerOptions);
      MGSRIndex::Reader indexReader = reader.getRoot<MGSRIndex>();
      LiteTree::Reader liteTreeReader = indexReader.getLiteTree();
      size_t numThreads = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
      bool lowMemory = vm.count("low-memory") > 0;
      
      mgsr::MgsrLiteTree liteTree;
      liteTree.initialize(indexReader, numThreads, lowMemory);
      
      std::vector<std::string> readSequences;
      mgsr::extractReadSequences(reads1, reads2, readSequences);

      bool skipSingleton = vm.count("skip-singleton") > 0;
      mgsr::ThreadsManager threadsManager(&liteTree, numThreads, skipSingleton, lowMemory);
      threadsManager.initializeMGSRIndex(indexReader);
      close(fd);
      threadsManager.initializeQueryData(readSequences);

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
            
      std::vector<uint64_t> totalNodesPerThread(numThreads, 0);
      for (size_t i = 0; i < numThreads; ++i) {
        totalNodesPerThread[i] = liteTree.getNumActiveNodes();
      }
      ProgressTracker progressTracker(numThreads, totalNodesPerThread);
      std::cout << "Using " << numThreads << " threads" << std::endl;

      auto start_time_place = std::chrono::high_resolution_clock::now();
      std::atomic<size_t> numGroupsUpdate = 0;
      std::atomic<size_t> numReadsUpdate = 0;
      std::vector<size_t> numUniqueKminmersPerThread(numThreads, 0);
      std::vector<size_t> numNodesPostCollapsePerThread(numThreads, 0);
      tbb::parallel_for(tbb::blocked_range<size_t>(0, threadsManager.threadRanges.size()), [&](const tbb::blocked_range<size_t>& rangeIndex){
        for (size_t i = rangeIndex.begin(); i != rangeIndex.end(); ++i) {
          auto [start, end] = threadsManager.threadRanges[i];
        
          std::span<mgsr::Read> curThreadReads(threadsManager.reads.data() + start, end - start);
          mgsr::mgsrPlacer curThreadPlacer(&liteTree, threadsManager, lowMemory, i);
          curThreadPlacer.initializeQueryData(curThreadReads);

          curThreadPlacer.setAllSeedmerHashesSet(threadsManager.allSeedmerHashesSet);

          curThreadPlacer.setProgressTracker(&progressTracker, i);
          // curThreadPlacer.placeReads();

          curThreadPlacer.scoreReads();

          threadsManager.readMinichainsInitialized[i] = curThreadPlacer.readMinichainsInitialized;
          threadsManager.readMinichainsAdded[i] = curThreadPlacer.readMinichainsAdded;
          threadsManager.readMinichainsRemoved[i] = curThreadPlacer.readMinichainsRemoved;
          threadsManager.readMinichainsUpdated[i] = curThreadPlacer.readMinichainsUpdated;
          if (i == 0) {
            threadsManager.identicalGroups = std::move(curThreadPlacer.identicalGroups);
            threadsManager.identicalNodeToGroup = std::move(curThreadPlacer.identicalNodeToGroup);
            // threadsManager.kminmerOverlapCoefficients = std::move(curThreadPlacer.kminmerOverlapCoefficients);
          }
          numGroupsUpdate += curThreadPlacer.numGroupsUpdate;
          numReadsUpdate += curThreadPlacer.numReadsUpdate;
        }
      });
      auto end_time_place = std::chrono::high_resolution_clock::now();
      auto duration_place = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_place - start_time_place);
      std::cerr << "\n\nPlaced reads in " << static_cast<double>(duration_place.count()) / 1000.0 << "s\n" << std::endl;

      mgsr::mgsrPlacer placerOC(&liteTree, threadsManager, lowMemory, 0);
      auto overlapCoefficients = placerOC.computeOverlapCoefficients(threadsManager.allSeedmerHashesSet);
      std::unordered_map<std::string, double>().swap(threadsManager.kminmerOverlapCoefficients);
      for (const auto& [nodeId, overlapCoefficient] : overlapCoefficients) {
        threadsManager.kminmerOverlapCoefficients[nodeId] = overlapCoefficient;
      }

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
    bool show_time = vm.count("time") > 0;
    bool dump_random_node = vm.count("dump-random-node") > 0;

    int k = vm["k"].as<int>(); // Default handled by boost
    int s = vm["s"].as<int>(); // Default handled by boost
    std::string index_path = vm["index"].as<std::string>(); // Default handled by boost

    std::mt19937 rng;
    if (vm.count("random-seed")) {
      std::string seed_str = vm["random-seed"].as<std::string>();
      std::hash<std::string> hasher;
      rng = std::mt19937(hasher(seed_str));
    } else {
      std::random_device rd;
      rng = std::mt19937(rd());
    }


    // Load pangenome
    msg("Loading reference pangenome from: {}", guide);

    std::unique_ptr<panmanUtils::TreeGroup> TG;

    try {
      TG = std::unique_ptr<panmanUtils::TreeGroup>(loadPanMAN(guide));

      if (!TG || TG->trees.empty()) {
        throw std::runtime_error(
            "No valid trees found in the loaded PanMAN file.");
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
    } catch (const std::exception &e) {
      err("Failed to load reference PanMAN: {}", e.what());
      return 1;
    }

    msg("Using first tree as reference.");
    panmanUtils::Tree &T = TG->trees[0];

    if (vm.count("index-mgsr")) {
      std::string mgsr_index_path = vm["index-mgsr"].as<std::string>();
      int mgsr_t = 0;
      int mgsr_l = 2;
      bool open = false;
      mgsr::mgsrIndexBuilder mgsrIndexBuilder(&T, 19, 8, mgsr_t, mgsr_l, open);
      mgsrIndexBuilder.buildIndex();
      mgsrIndexBuilder.writeIndex(mgsr_index_path);
      msg("MGSR index written to: {}", mgsr_index_path);
      return 0;
    }



    // Handle --dump-random-node parameter if provided
    if (dump_random_node) {
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

      std::string outputFileName =
          guide + ".random." + randomNode->identifier + ".fa";
      if (saveNodeSequence(nodeTree, randomNode, outputFileName)) {
        msg("Random node {} sequence written to {}", randomNode->identifier,
            outputFileName);
        return 0; // Exit after dumping sequence
      } else {
        err("Failed to save random node sequence to {}", outputFileName);
        return 1;
      }
    }


    if (vm.count("dump-random-nodeIDs")) {
      uint32_t num_nodes = vm["dump-random-nodeIDs"].as<uint32_t>();
      std::vector<std::string_view> allNodeIDs;
      allNodeIDs.reserve(T.allNodes.size());
      for (const auto& [nodeID, node] : T.allNodes) {
        if (node->children.empty()) {
          allNodeIDs.push_back(nodeID);
        }
      }
      allNodeIDs.shrink_to_fit();

      std::shuffle(allNodeIDs.begin(), allNodeIDs.end(), rng);

      std::ofstream outFile(prefix + ".randomNodeIDs.txt");
      for (size_t i = 0; i < num_nodes; i++) {
        outFile << allNodeIDs[i] << std::endl;
      }
      outFile.close();
      msg("Random node IDs written to {}", prefix + ".randomNodeIDs.txt");
      exit(0);
    }

    // Handle --dump-sequence parameter if provided
    if (vm.count("dump-sequence")) {
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

    if (vm.count("dump-sequences")) {
      auto nodeIDs = vm["dump-sequences"].as<std::vector<std::string>>();
      std::vector<uint32_t> numsnps;
      if (vm.count("simulate-snps")) {
        numsnps = vm["simulate-snps"].as<std::vector<uint32_t>>();
        if (numsnps.size() != nodeIDs.size()) {
          err("Number of SNP parameters does not match number of node IDs");
          return 1;
        }
      }

      for (size_t i = 0; i < nodeIDs.size(); i++) {
        const auto& nodeID = nodeIDs[i];
        uint32_t numsnp = (numsnps.empty() ? 0 : numsnps[i]);
        if (T.allNodes.find(nodeID) == T.allNodes.end()) {
          err("Node ID {} not found in the tree", nodeID);
          return 1;
        }

        std::string sequence = panmapUtils::getStringFromReference(&T, nodeID, false);
        std::vector<std::tuple<char, char, uint32_t>> snpRecords;
        panmapUtils::simulateSNPsOnSequence(sequence, snpRecords, numsnp, rng);
        std::string nodeIDClean = nodeID;
        std::replace(nodeIDClean.begin(), nodeIDClean.end(), '/', '_');
        std::replace(nodeIDClean.begin(), nodeIDClean.end(), '|', '_');
        std::string outputFileName = prefix + "." + nodeIDClean + "." + std::to_string(numsnp) + "snps.fa";
        std::ofstream outFile(outputFileName);

        if (outFile.is_open()) {
          outFile << ">" << nodeID << " ";
          for (const auto& [ref, alt, pos] : snpRecords) {
            outFile << ref << pos << alt << " ";
          }
          outFile << "\n";
          for (size_t i = 0; i < sequence.size(); i += 80) {
            outFile << sequence.substr(i, 80) << "\n";
          }
          outFile.close();
          msg("Sequence for node {} with {} SNPs written to {}", nodeID, numsnp, outputFileName);
        } else {
          err("Failed to open file {} for writing", outputFileName);
          return 1;
        }
      }
      exit(0);
    }

    // Log settings
    msg("--- Settings ---");
    msg("Reads: {}",
        (reads1.empty() ? "<none>"
                        : reads1 + (reads2.empty() ? "" : " + " + reads2)));
    msg("Reference PanMAN: {} ({} nodes)", guide, T.allNodes.size());
    msg("Using {} threads", vm["cpus"].as<int>());
    msg("k={}, s={}", k, s);

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


    // Build index if needed
    if (build) {
      logging::debug("DEBUG-INDEX: Starting index build process");
      try {
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
      msg("Building new mutation matrices");
      genotyping::fillMutationMatricesFromTree_test(mutMat, &T,
                                                    default_mutmat_path);
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

          // Handle batch file if provided
          if (vm.count("batch")) {
            std::string batchFilePath = vm["batch"].as<std::string>();
            try {
              msg("=== Starting Batch Placement ===");
              // Process batch file - keep inMessage alive during the whole process
              placement::placeBatch(&T, index_input, batchFilePath, prefix,
                                    refFileName, samFileName, bamFileName,
                                    mpileupFileName, vcfFileName, aligner, refNode,
                                    vm.count("save-jaccard") > 0, vm.count("time") > 0,
                                    vm["candidate-threshold"].as<float>(), 
                                    vm["max-candidates"].as<int>(),      
                                    effective_index_path, 
                                    debug_specific_node_id);
              msg("Batch placement completed.");
            } catch (const std::exception &e) {
              err("ERROR during batch processing: {}", e.what());
            }
            return 0;
          }

          // Initialize placement variables and read information needed for alignment
          placement::PlacementResult result;
          std::vector<std::vector<seeding::seed_t>> readSeeds;
          std::vector<std::string> readSequences;
          std::vector<std::string> readNames;
          std::vector<std::string> readQuals;

          std::string placementFileName = prefix + ".placement.tsv";
          // Perform placement - keep inMessage alive during the whole process
          placement::place(result, &T, index_input, reads1, reads2,
                          readSeeds, readSequences, readNames, readQuals,
                         placementFileName, effective_index_path, debug_specific_node_id);


          /* @Alan this is the end of placement
            -> next step is pass seeds of target node to Nico's alignment code
            -> and genotyping

            placementResult should have nodeSeedMap[targetId] with the target seed set
            
            TODO: I'm not sure k-mer end positions are correct yet
          */
          

          std::string bestMatchSequence = panmapUtils::getStringFromReference(&T, result.bestJaccardPresenceNode->identifier, false);
          logging::info("Best match sequence built for {} with length {}", result.bestJaccardPresenceNode->identifier, bestMatchSequence.size());


          
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

          createMplpBcf(prefix, refFileName, bestMatchSequence, bamFileName, mpileupFileName);

          createVcfWithMutationMatrices(prefix, mpileupFileName, mutMat, vcfFileName, 0.0011);

          std::string placementSummaryFileName = prefix + ".placement.summary.md";
          placement::dumpPlacementSummary(result, placementSummaryFileName);

          // Report the top 3 metrics as requested
          msg("=== TOP PLACEMENT RESULTS ===");
          
          // 1. Raw number of matches (set, no frequency) - using raw count
          msg("Top by raw seed matches (unique count): node {} with {} matches (Jaccard score {:.4f})",
              result.bestJaccardPresenceNode ? result.bestJaccardPresenceNode->identifier : "none", 
              result.bestJaccardPresenceCount,
              result.bestJaccardPresenceScore);
          
          // 2. Number of matches scaled by read frequency 
          msg("Top by weighted seed matches (frequency-scaled): node {} with score {}",
              result.bestRawSeedMatchNode ? result.bestRawSeedMatchNode->identifier : "none", 
              result.bestRawSeedMatchScore);
          
          // 3. Cosine similarity
          msg("Top by cosine similarity: node {} with score {:.4f}",
              result.bestCosineNode ? result.bestCosineNode->identifier : "none", 
              result.bestCosineScore);
          
          msg("=== ADDITIONAL METRICS ===");
          msg("Best Weighted Jaccard: node {} with score {:.4f}",
              result.bestJaccardNode ? result.bestJaccardNode->identifier : "none", result.bestJaccardScore);
          msg("Overall Best Weighted Score (Jaccard*scale + Cosine*(1-scale)): node {} with score {:.4f}",
              result.bestWeightedNode ? result.bestWeightedNode->identifier : "none", result.bestWeightedScore);
          
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
    msg("=== panmap run completed ===");
    msg("Total runtime: {}ms", total_duration.count());

    
    return 0;
}
