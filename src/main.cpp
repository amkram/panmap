#include <boost/iostreams/categories.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <cctype>
#include <chrono>
#include <cmath>
#include <csignal>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
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
#include <unistd.h>
#include <vector>
#include <stack>

#include "capnp/message.h"
#include "capnp/serialize-packed.h"

#include "docopt.h"

#include "docopt_value.h"
#include "genotyping.hpp"
#include "index.capnp.h"
#include "indexing.hpp"
#include "logging.hpp"
#include "panman.hpp"
#include "placement.hpp"
#include "timing.hpp"
#include <atomic>
#include <boost/filesystem.hpp>
#include <boost/filesystem/exception.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/program_options.hpp>
#include <htslib/kseq.h>
#include <panmanUtils.hpp>
#include <cxxabi.h>
#include <execinfo.h>

static const char USAGE[] =
    R"(panmap -- v0.0 ⟗ ⟗

Pangenome phylogenetic placement, alignment, genotyping, and assembly of reads

Usage:
  panmap [options] <guide.panman> [<reads1.fastq>] [<reads2.fastq>]

<guide.panman>  Path to pangenome reference (PanMAN file).
<reads1.fastq>  [Optional] Path to first FASTQ file.
<reads2.fastq>  [Optional] Path to second FASTQ file (for paired-end reads).

Input/output options:
  -b --batch <tsv>          Path to tsv file with reads to process (one per line or two per line for paired-end reads).
  -p --prefix <prefix>      Prefix for output files. [default: panmap]
  -o <outputs>              List of outputs to generate.
                              Accepted values:
                                  placement / p:  Save phylogenetic placement results to <prefix>.placement.tsv (or <prefix>.placements.tsv in batch mode)
                                  assembly / a:   Save consensus assembly to <prefix>.assembly.fa
                                  reference / r:  Save panmap's selected reference to <prefix>.reference.fa
                                  spectrum / c:   Save mutation spectrum matrix to <prefix>.mm
                                  sam / s:        Save aligned reads to <prefix>.sam
                                  bam / b:        Save aligned reads to <prefix>.bam
                                  mpileup / m:    Save read pileup to <prefix>.mpileup
                                  vcf / v:        Save variant calls and likelihoods to <prefix>.vcf
                                  all / A:        Save all possible outputs, each to <prefix>.<ext>
                              [default: bam,vcf,assembly]

  -i --index <path>         Provide a precomputed panmap index. If not specified, index
                                is loaded from <guide.pmat>.pmi, if it exists, otherwise
                                it is built with parameters <k> and <s> (see below). [default: ]

  -q --quiet                Run in quiet mode (show only critical messages, errors, and warnings).
  -v --verbose              Run in verbose mode (show detailed logs for debugging).

Seeding/alignment options:
  -k <k>                             Length of k-mer seeds. [default: 32]
  -s <s>                             Length of s-mers for seed (syncmer) selection. [default: 8]
  -a --aligner <method>              The alignment algorithm to use ('minimap2' or 'bwa-aln').
                                     For very short or damaged reads, use bwa-aln. [default: minimap2]
  -r --ref <node_id>                 Provide a reference node to align reads to (skip placement) [default: ].
  -P --prior                         Compute and use a mutation spectrum prior for genotyping.
  -f --reindex                       Don't load index from disk, build it from scratch.

Other options:
  -c --cpus <num>            Number of CPUs to use. [default: 1]
  -x --stop-after <stage>    Stop after the specified stage. Accepted values:
                                    indexing / i:   Stop after seed indexing
                                    placement / p:  Stop after placement
                                    mapping / m:    Stop after read mapping
                                    genotyping / g: Stop after computing genotype likelihoods
                                    assembly / a:   Stop after consensus assembly
                                 [default: ]
  -Q <seed>                 Integer seed for random number generation. [default: 42]
  -V --version              Show version.
  --time                    Show time taken at each step
  -h --help                 Show this screen.
  -D --dump                 Dump all seeds to file.
  -X --dump-real            Dump true seeds to file.
  --test                    Run coordinate system tests.

Developer options:
  --genotype-from-sam                   Generate VCF from SAM file using mutation spectrum as prior.
  --sam-file <path>                     Path to SAM file to generate VCF from.
  --ref-file <path>                     Path to reference FASTA file to generate VCF from.
  --save-jaccard                        Save jaccard index between reads and haplotypes to <prefix>.jaccard.txt
  --save-kminmer-binary-coverage        Save kminmer binary coverage to <prefix>.kminmer_binary_coverage.txt
  --parallel-tester                     Run parallel tester.
  --eval <path>                         Evaluate accuracy of placement. <path> is a tsv file.
  --dump-sequence <nodeID>              Dump sequence for a specified node ID to <panmanfile>.<nodeID>.fa
  --dump-random-node                    Dump sequence for a random node to <panmanfile>.random.<nodeID>.fa

Placement options:
  --candidate-threshold <proportion>    Initial seed matching retains the top <proportion> of all nodes as candidates, whose scores are refined by seed extension. [default: 0.01]
  --max-candidates <count>              Maximum number of candidate nodes to consider [default: 16]
)";

using namespace logging;

// Namespace aliases
namespace fs = boost::filesystem;

// using namespace coordinates;
// using namespace seed_annotated_tree;
// using namespace placement;
// using namespace genotyping;

// Global constants
static constexpr int DEFAULT_K = 32;
static constexpr int DEFAULT_S = 8;

// Thread safety
std::mutex outputMutex, indexMutex, placementMutex;
std::atomic<bool> shouldStop{false};

void writeCapnp(::capnp::MallocMessageBuilder &message,
                const std::string &path) {
  std::lock_guard<std::mutex> lock(indexMutex);
  
  // Use the indexing namespace version with forceFlush=true
  indexing::periodicallyFlushCapnp(message, path, true, 0, 1000);
}

std::unique_ptr<::capnp::PackedFdMessageReader>
readCapnp(const std::string &path) {
  std::lock_guard<std::mutex> lock(indexMutex);
  int fd = open(path.c_str(), O_RDONLY);
  if (fd < 0) {
    throw std::runtime_error("Failed to open file: " + path);
  }

  ::capnp::ReaderOptions opts;
  opts.traversalLimitInWords = std::numeric_limits<uint64_t>::max();
  opts.nestingLimit = 1024;

  // Create the reader that manages the file descriptor
  // NOTE: Don't close the file descriptor here - PackedFdMessageReader owns it
  // and will close it when it's destroyed
  return std::unique_ptr<::capnp::PackedFdMessageReader>(
      new ::capnp::PackedFdMessageReader(fd, opts));
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
 * Captures the stack trace and sets the shouldStop flag
 * 
 * @param signum The signal number
 */
void signalHandler(int signum) {
    // Get the signal name
    const char* sigName = "UNKNOWN";
    switch (signum) {
        case SIGINT:  sigName = "SIGINT (Interrupt/Ctrl+C)"; break;
        case SIGTERM: sigName = "SIGTERM (Termination request)"; break;
        case SIGSEGV: sigName = "SIGSEGV (Segmentation fault)"; break;
        case SIGABRT: sigName = "SIGABRT (Abort)"; break;
    }
    
    std::cerr << "\n\n*** Program interrupted by signal: " << sigName << " ***" << std::endl;
    
    // Print the current stack trace
    printStackTrace(2);  // Skip signalHandler and signal handler dispatcher frames
    
    std::cerr << "\nGracefully shutting down..." << std::endl;
    shouldStop = true;
    
    // Attempt to trigger atexit handlers (like gprof's) for non-fatal signals
    if (signum == SIGINT || signum == SIGTERM) {
        std::cerr << "Attempting graceful exit to finalize profiling data (gmon.out)..." << std::endl;
        // Exit with a non-zero status to indicate termination by signal
        // This should allow atexit handlers to run.
        exit(128 + signum); 
    } else if (signum == SIGSEGV || signum == SIGABRT) {
        // For fatal errors, exit immediately after printing trace.
        // gmon.out is unlikely to be useful or complete here anyway.
        std::cerr << "Fatal error, exiting immediately." << std::endl;
        exit(128 + signum);
    }
    // If it's another signal, we just set shouldStop and let the program potentially handle it,
    // though this handler isn't registered for other signals currently.
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

// Main program entry point
int main(int argc, char *argv[]) {
    // Register signal handlers for various signals
    signal(SIGINT, signalHandler);   // Ctrl+C
    signal(SIGTERM, signalHandler);  // Termination request
    signal(SIGSEGV, signalHandler);  // Segmentation fault
    signal(SIGABRT, signalHandler);  // Abort
    
    auto main_start = std::chrono::high_resolution_clock::now();
    
    std::map<std::string, docopt::value> args =
        docopt::docopt(USAGE, {argv + 1, argv + argc}, true, "panmap 0.0");

    // Set up parallel processing
    tbb::global_control c(tbb::global_control::max_allowed_parallelism,
                          std::stoi(args["--cpus"].asString()));

    // Set logging verbosity level
    bool quiet_mode = false;
    bool verbose_mode = false;

    // Check if the quiet and verbose flags exist in the args map
    quiet_mode = args["--quiet"] && args["--quiet"].isBool()
                ? args["--quiet"].asBool()
                : false;

    verbose_mode = args["--verbose"] && args["--verbose"].isBool()
                  ? args["--verbose"].asBool()
                  : false;

    // Set logging level based on flags
    if (quiet_mode) {
      logging::setLoggingLevel(logging::LogLevel::QUIET);
      // Only critical messages will be shown
    } else if (verbose_mode) {
      logging::setLoggingLevel(logging::LogLevel::VERBOSE);
      // Initialize node tracking for detailed logging
      logging::initNodeTracking();
    } else {
      // Default is NORMAL level
      logging::setLoggingLevel(logging::LogLevel::INFO);
    }

    // Display startup message based on verbosity
    logging::critical("Starting panmap with verbosity level: {}",
                      quiet_mode ? "quiet"
                                 : (verbose_mode ? "verbose" : "normal"));

    // Get basic parameters
    if (!args["<guide.panman>"] || !args["<guide.panman>"].isString()) {
      err("Missing required pangenome file argument");
      return 1;
    }
    std::string guide = args["<guide.panman>"].asString();
    std::string reads1 =
        args["<reads1.fastq>"] ? args["<reads1.fastq>"].asString() : "";
    std::string reads2 =
        args["<reads2.fastq>"] ? args["<reads2.fastq>"].asString() : "";
    std::string prefix =
        args["--prefix"] ? args["--prefix"].asString() : "panmap";
    std::string outputs = args["-o"] && args["-o"].isString()
                              ? args["-o"].asString()
                              : "bam,vcf,assembly";
    std::string aligner =
        args["--aligner"].asString() == "bwa-aln" ? "bwa-aln" : "minimap2";
    std::string refNode = args["--ref"] ? args["--ref"].asString() : "";
    std::string eval = args["--eval"] ? args["--eval"].asString() : "";

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

    // Get remaining parameters
    bool reindex = args["--reindex"] && args["--reindex"].isBool()
                       ? args["--reindex"].asBool()
                       : false;
    bool prior = args["--prior"] && args["--prior"].isBool()
                     ? args["--prior"].asBool()
                     : false;

    bool genotype_from_sam =
        args["--genotype-from-sam"] && args["--genotype-from-sam"].isBool()
            ? args["--genotype-from-sam"].asBool()
            : false;
    bool save_jaccard = args["--save-jaccard"] && args["--save-jaccard"].isBool()
                            ? args["--save-jaccard"].asBool()
                            : false;
    bool show_time = args["--time"] && args["--time"].isBool()
                         ? args["--time"].asBool()
                         : false;
    bool dump_random_node =
        args["--dump-random-node"] && args["--dump-random-node"].isBool()
            ? args["--dump-random-node"].asBool()
            : false;

    // Use DEFAULT_K and DEFAULT_S if not specified in command line args
    int k = args["-k"] ? std::stoi(args["-k"].asString()) : DEFAULT_K;
    int s = args["-s"] ? std::stoi(args["-s"].asString()) : DEFAULT_S;
    std::string index_path = args["--index"] ? args["--index"].asString() : "";

    std::random_device rd;
    std::mt19937 rng(rd());

    // Load pangenome
    msg("Loading reference pangenome from: {}", guide);

    std::unique_ptr<panmanUtils::TreeGroup> TG;

    logging::setLoggingLevel(logging::LogLevel::VERBOSE);
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

    // Handle --dump-sequence parameter if provided
    if (args["--dump-sequence"] && args["--dump-sequence"].isString()) {
      std::string nodeID = args["--dump-sequence"].asString();
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
    msg("Reference PanMAN: {} ({} nodes)", guide, T.allNodes.size());
    msg("Using {} threads", args["--cpus"].asString());
    msg("k={}, s={}", k, s);

    // Handle indexing
    bool build = true;
    ::capnp::MallocMessageBuilder outMessage;
    std::unique_ptr<::capnp::PackedFdMessageReader> inMessage;
    std::string default_index_path = guide + ".pmi";

    if (!index_path.empty() && !reindex) {
      msg("Using provided index: {}", index_path);
    } else if (!reindex && fs::exists(default_index_path)) {
      msg("Loading existing index from: {}", default_index_path);
      inMessage = readCapnp(default_index_path);
      build = false;
    } else {
      msg(reindex ? "Will re-build index (-f)"
                  : "No index found, will build new index");
    }

    // Build index if needed
    if (build) {

      try {
        Index::Builder index = outMessage.initRoot<Index>();
        if (k <= 0 || s <= 0 || s >= k) {
          throw std::runtime_error("Invalid k-mer or s-mer parameters");
        }

        auto start = std::chrono::high_resolution_clock::now();
        indexing::index(&T, index, k, s, outMessage, default_index_path);
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start);

        writeCapnp(outMessage, default_index_path);
        msg("Index written to: {}", default_index_path);
      } catch (const std::exception &e) {
        err("ERROR during indexing: {}", e.what());
        return 1;
      }
    }

    // Load mutation matrices
    msg("=== Loading Mutation Matrices ===");
    genotyping::mutationMatrices mutMat;
    std::string mutmat_path = args["--mutmat"] ? args["--mutmat"].asString() : "";
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
      if (!inMessage) {
        inMessage = readCapnp(default_index_path);
      }
      Index::Reader index_input = inMessage->getRoot<Index>();
      
      if (genotype_from_sam) {
        msg("Genotyping from SAM file");
        if (samFileName.empty() || refFileName.empty()) {
          throw std::runtime_error(
              "SAM file and reference file are required for genotyping from SAM");
        }
        genotyping::genotype(prefix, refFileName, "", bamFileName,
                             mpileupFileName, vcfFileName, mutMat);
      } else {
        if (!refNode.empty()) {
          msg("Using reference node: {}", refNode);
          if (T.allNodes.find(refNode) == T.allNodes.end()) {
            throw std::runtime_error("Reference node not found in pangenome");
          }
        }

        // Check if batch mode is enabled
        std::string batchFilePath =
            args["--batch"] ? args["--batch"].asString() : "";
        if (!batchFilePath.empty()) {
          msg("Running in batch mode with file: {}", batchFilePath);

          // Validate batch file exists
          if (!fs::exists(batchFilePath)) {
            throw std::runtime_error("Batch file not found: " + batchFilePath);
          }

          // Process batch file - keep inMessage alive during the whole process
          placement::placeBatch(&T, index_input, batchFilePath, prefix,
                                refFileName, samFileName, bamFileName,
                                mpileupFileName, vcfFileName, aligner, refNode,
                                save_jaccard, show_time, 0.01f, 16);

          msg("Batch processing completed");
        } else {
          // Process single sample mode
          // Initialize placement variables
          placement::PlacementResult result;

          std::string placementFileName = prefix + ".placement.tsv";
          // Perform placement - keep inMessage alive during the whole process
          placement::place(result, &T, index_input, reads1, reads2,
                           placementFileName);

          // Get placement results from result struct
          panmanUtils::Node *maxHitsNode = result.maxHitsNode;
          int64_t maxHitsInAnyGenome = result.maxHitsInAnyGenome;

          msg("Placed sample at {} ({} seed matches)",
              maxHitsNode ? maxHitsNode->identifier : "none", maxHitsInAnyGenome);
        }
      }
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

    if (show_time) {
      timing::Timer::report();
    }

    // Final cleanup 
    
    // Flush any remaining log entries for node_1 and node_2
    logging::flushNodeTracking();
    
    return 0;
}
