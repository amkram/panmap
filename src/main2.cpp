#include "indexing.hpp"
#include "placement.hpp"
#include "genotyping.hpp"
#include "logging.hpp"
#include <panmanUtils.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/program_options.hpp>
#include <minimap2/kseq.h>
#include <tbb/global_control.h>
#include "docopt.h"
#include "index.capnp.h"
#include <fcntl.h>


static const char USAGE[] = 
R"(panmap -- v0.0 ⟗ ⟗

Pangenome phylogenetic placement, alignment, genotyping, and assembly of reads

Usage:
  panmap <guide.panman> [<reads1.fastq>] [<reads2.fastq>] [options]

<guide.panman>  Path to pangenome reference (PanMAN file).
<reads1.fastq>  [Optional] Path to first FASTQ file.
<reads2.fastq>  [Optional] Path to second FASTQ file (for paired-end reads).

Input/output options:
  -p --prefix <prefix>      Prefix for output files. [default: panmap]
  -o <outputs>              List of outputs to generate.
                              Accepted values:
                                  placement / p:  Save phylogenetic placement results to <prefix>.placement
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
  
  -m --mutmat <path>        Provide a precomputed mutation spectrum matrix instead of computing
                                one from the tree. if not specified, one is computed from the tree.
                                Overrides --prior. [default: ]

Seeding/alignment options:
  -k <k>                             Length of k-mer seeds. [default: 19]
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
  --hot-threshold <int>     Mutation frequency threshold for hot positions to index. [default: 1000]

Developer options:
  --genotype-from-sam                   Generate VCF from SAM file using mutation spectrum as prior.
  --sam-file <path>                     Path to SAM file to generate VCF from.
  --ref-file <path>                     Path to reference FASTA file to generate VCF from.
  --save-jaccard                        Save jaccard index between reads and haplotypes to <prefix>.jaccard.txt
  --save-kminmer-binary-coverage        Save kminmer binary coverage to <prefix>.kminmer_binary_coverage.txt
  --parallel-tester                     Run parallel tester.
  --eval <path>                         Evaluate accuracy of placement. <path> is a tsv file.

Placement options:
  --candidate-threshold <proportion>    Initial seed matching retains the top <proportion> of all nodes as candidates, whose scores are refined by seed extension. [default: 0.01]
  --max-candidates <count>              Maximum number of candidate nodes to consider [default: 16]
)";


// Namespace aliases
namespace fs = boost::filesystem;
namespace po = boost::program_options;
namespace sat = seed_annotated_tree;
using namespace logging;

using namespace coordinates;  // For coordinate types and traversal
using namespace seed_annotated_tree;
using namespace placement;
using namespace genotyping;

// Global constants
static constexpr int DEFAULT_K = 31;
static constexpr int DEFAULT_S = 15;
static constexpr int64_t DEFAULT_HOT_THRESHOLD = 1000;
static constexpr bool DEFAULT_REINDEX = false;

// Thread safety
namespace {
    std::mutex outputMutex, indexMutex, placementMutex;
    std::atomic<bool> shouldStop{false};

    void log(const std::string& prefix, const std::string& message) {
        std::lock_guard<std::mutex> lock(outputMutex);
        std::cout << "[" << prefix << "] " << message << std::endl;
    }

    void writeCapnp(::capnp::MallocMessageBuilder& message, const std::string& path) {
        std::lock_guard<std::mutex> lock(indexMutex);
        int fd = open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if (fd < 0) throw std::runtime_error("Failed to open file for writing: " + path);
        writePackedMessageToFd(fd, message);
        close(fd);
    }

    std::unique_ptr<::capnp::PackedFdMessageReader> readCapnp(const std::string& path) {
        std::lock_guard<std::mutex> lock(indexMutex);
        int fd = open(path.c_str(), O_RDONLY);
        if (fd < 0) throw std::runtime_error("Failed to open file for reading: " + path);
        auto reader = std::make_unique<::capnp::PackedFdMessageReader>(fd);
        close(fd);
        return reader;
    }

    void signalHandler(int signum) {
        shouldStop = true;
    }

    bool isFileReadable(const std::string& path) {
        std::ifstream file(path);
        return file.good();
    }

    bool isFileWritable(const std::string& path) {
        if (fs::exists(path)) {
            std::ofstream file(path, std::ios::app);
            return file.good();
        }
        
        auto dir = fs::path(path).parent_path();
        if (!fs::exists(dir)) {
            try {
                fs::create_directories(dir);
            } catch (const fs::filesystem_error&) {
                return false;
            }
        }
        
        std::ofstream file(path);
        if (!file.good()) return false;
        fs::remove(path);
        return true;
    }

    void validateInputFile(const std::string& path, const std::string& description) {
        if (path.empty()) throw std::runtime_error(description + " path is empty");
        if (!fs::exists(path)) throw std::runtime_error(description + " not found: " + path);
        if (!isFileReadable(path)) throw std::runtime_error("Cannot read " + description + ": " + path);
    }

    void validateOutputFile(const std::string& path, const std::string& description) {
        if (path.empty()) throw std::runtime_error(description + " path is empty");
        if (!isFileWritable(path)) throw std::runtime_error("Cannot write to " + description + ": " + path);
    }
}

// Forward declarations
static void writeCapnp(::capnp::MallocMessageBuilder& message, const std::string& path);
static std::unique_ptr<::capnp::PackedFdMessageReader> readCapnp(const std::string& path);

using namespace std;

void writeCapnp(::capnp::MallocMessageBuilder &message, std::string &filename) {
    int fd = open(filename.c_str(), O_WRONLY | O_CREAT, 0644);
    if (fd < 0) {
        err("Failed to open proto file for writing");
        return;
    }
    try {
        capnp::writePackedMessageToFd(fd, message);
    } catch (const std::exception &e) {
        err("Failed to write Cap'n Proto message: {}", e.what());
    }
    close(fd);
}

std::unique_ptr<::capnp::PackedFdMessageReader> readCapnp(std::string &filename) {
    int fd = open(filename.c_str(), O_RDONLY);
    ::capnp::ReaderOptions options = {(uint64_t) -1, 64};
    return std::make_unique<::capnp::PackedFdMessageReader>(fd, options);
}

panmanUtils::Tree* loadPanmanOrPanmat(const std::string &pmatFile) {
    panmanUtils::Tree *T = nullptr;
    panmanUtils::TreeGroup *TG = nullptr;

    try {
        // Try loading as PanMAN
        std::ifstream inputFile(pmatFile);
        boost::iostreams::filtering_streambuf<boost::iostreams::input> inBuffer;
        inBuffer.push(boost::iostreams::lzma_decompressor());
        inBuffer.push(inputFile);
        std::istream inputStream(&inBuffer);
        TG = new panmanUtils::TreeGroup(inputStream);
        inputFile.close();
    } catch (const std::exception &e) {
        msg("Attempting to load as PanMAT...");
        try {
            std::ifstream inputFile(pmatFile);
            boost::iostreams::filtering_streambuf<boost::iostreams::input> inBuffer;
            inBuffer.push(boost::iostreams::lzma_decompressor());
            inBuffer.push(inputFile);
            std::istream inputStream(&inBuffer);
            T = new panmanUtils::Tree(inputStream);
        } catch (const std::exception &e) {
            return nullptr;
        }
    }
    return TG ? &(TG->trees[0]) : T;
}

// Memory management helper
class ScopedLinkedNodeCleaner {
    std::vector<LinkedNode*> nodes;

    void cleanupLinkedNodes(LinkedNode* node) {
        if (!node) return;
        for (auto* child : node->children) {
            cleanupLinkedNodes(child);
        }
        delete node;
    }

public:
    void add(LinkedNode* node) {
        if (node) nodes.push_back(node);
    }

    ~ScopedLinkedNodeCleaner() {
        for (auto* node : nodes) {
            cleanupLinkedNodes(node);
        }
    }
};

int main(int argc, const char** argv) {
    auto main_start = std::chrono::high_resolution_clock::now();
    
    // Parse command line arguments
    std::map<std::string, docopt::value> args = docopt::docopt(USAGE, { argv + 1, argv + argc }, true, "panmap 0.0");
    
    // Set up parallel processing
    tbb::global_control c(tbb::global_control::max_allowed_parallelism, std::stoi(args["--cpus"].asString()));

    msg("=== Starting panmap execution ===");
    msg("Command line arguments parsed successfully");
    msg("Parallel processing initialized with {} threads", args["--cpus"].asString());

    // Get basic parameters
    if (!args["<guide.panman>"] || !args["<guide.panman>"].isString()) {
        err("Missing required pangenome file argument");
        return 1;
    }
    std::string guide = args["<guide.panman>"].asString();
    std::string reads1 = args["<reads1.fastq>"] ? args["<reads1.fastq>"].asString() : "";
    std::string reads2 = args["<reads2.fastq>"] ? args["<reads2.fastq>"].asString() : "";
    std::string prefix = args["--prefix"] ? args["--prefix"].asString() : "panmap";
    std::string outputs = args["-o"] && args["-o"].isString() ? args["-o"].asString() : "bam,vcf,assembly";
    std::string aligner = args["--aligner"].asString() == "bwa-aln" ? "bwa-aln" : "minimap2";
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
    std::string refFileName, samFileName, bamFileName, mpileupFileName, vcfFileName;
    
    for (const auto& output : outputs_seperated) {
        if (output.size() == 1) {
            switch(output[0]) {
                case 'r': refFileName = prefix + ".reference.fa"; break;
                case 's': samFileName = prefix + ".sam"; break;
                case 'b': bamFileName = prefix + ".bam"; break;
                case 'm': mpileupFileName = prefix + ".mpileup"; break;
                case 'v': vcfFileName = prefix + ".vcf"; break;
                case 'A': 
                    refFileName = prefix + ".reference.fa";
                    samFileName = prefix + ".sam";
                    bamFileName = prefix + ".bam";
                    mpileupFileName = prefix + ".mpileup";
                    vcfFileName = prefix + ".vcf";
                    break;
            }
        } else {
            if (output == "reference") refFileName = prefix + ".reference.fa";
            else if (output == "sam") samFileName = prefix + ".sam";
            else if (output == "bam") {
                samFileName = prefix + ".sam";
                bamFileName = prefix + ".bam";
            }
            else if (output == "mpileup") mpileupFileName = prefix + ".mpileup";
            else if (output == "vcf") vcfFileName = prefix + ".vcf";
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
    bool reindex = args["--reindex"] && args["--reindex"].isBool() ? args["--reindex"].asBool() : false;
    bool prior = args["--prior"] && args["--prior"].isBool() ? args["--prior"].asBool() : false;
    bool placement_per_read = args["--place-per-read"] && args["--place-per-read"].isBool() ? args["--place-per-read"].asBool() : false;
    bool genotype_from_sam = args["--genotype-from-sam"] && args["--genotype-from-sam"].isBool() ? args["--genotype-from-sam"].asBool() : false;
    bool save_jaccard = args["--save-jaccard"] && args["--save-jaccard"].isBool() ? args["--save-jaccard"].asBool() : false;
    bool show_time = args["--time"] && args["--time"].isBool() ? args["--time"].asBool() : false;

    int k = std::stoi(args["-k"].asString());
    int s = std::stoi(args["-s"].asString());
    std::string index_path = args["--index"] ? args["--index"].asString() : "";

    // Load pangenome
    msg("Loading pangenome from: {}", guide);
    panmanUtils::Tree *T = loadPanmanOrPanmat(guide);
    if (!T) {
        err("Failed to load guide panMAN/panMAT");
        return 1;
    }
    msg("Successfully loaded pangenome with {} nodes", T->allNodes.size());

    // Log settings
    msg("--- Settings ---");
    msg("Pangenome: {} ({} nodes)", guide, T->allNodes.size());
    msg("Reads: {}", (reads1.empty() ? "<none>" : reads1 + (reads2.empty() ? "" : " + " + reads2)));
    msg("Output prefix: {}", prefix);
    msg("Outputs: {}", outputs);
    msg("Reindex: {}", reindex);
    msg("k-mer length: {}", k);
    msg("s-mer length: {}", s);

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
        msg(reindex ? "Reindexing requested" : "No existing index found, will build new index");
    }

    if (args["--stop-after"]) {
        msg("Will stop after stage: {}", args["--stop-after"].asString());
    }

    int64_t hot_threshold = args["--hot-threshold"] ? std::stoi(args["--hot-threshold"].asString()) : 1000;
    msg("Hot threshold set to: {}", hot_threshold);

    // Build index if needed
    if (build) {
        msg("=== Building Index ===");
        msg("Initializing index parameters:");
        msg("  k-mer length: {}", k);
        msg("  s-mer length: {}", s);

        try {
            Index::Builder index = outMessage.initRoot<Index>();
            if (k <= 0 || s <= 0 || s >= k) {
                throw std::runtime_error("Invalid k-mer or s-mer parameters");
            }
            
            auto start = std::chrono::high_resolution_clock::now();
            indexing::index(T, index, hot_threshold, k, s);
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - start
            );
            msg("Index build completed in {}ms", duration.count());
            
            writeCapnp(outMessage, default_index_path);
            msg("Index written to: {}", default_index_path);
            
            // If no reads were provided, we're done after building the index
            if (reads1.empty() && reads2.empty()) {
                msg("No reads provided - index building complete, exiting");
                return 0;
            }
        } catch (const std::exception& e) {
            err("ERROR during indexing: {}", e.what());
            return 1;
        }
    }

    // Load mutation matrices
    msg("=== Loading Mutation Matrices ===");
    sat::mutationMatrices mutMat;
    std::string mutmat_path = args["--mutmat"] ? args["--mutmat"].asString() : "";
    std::string default_mutmat_path = guide + ".mm";
    
    if (!mutmat_path.empty()) {
        msg("Loading mutation matrices from: {}", mutmat_path);
        std::ifstream mutmat_file(mutmat_path);
        sat::fillMutationMatricesFromFile(mutMat, mutmat_file);
    } else if (fs::exists(default_mutmat_path)) {
        msg("Loading default mutation matrices from: {}", default_mutmat_path);
        std::ifstream mutmat_file(default_mutmat_path);
        sat::fillMutationMatricesFromFile(mutMat, mutmat_file);
    } else {
        msg("Building new mutation matrices");
        sat::fillMutationMatricesFromTree_test(mutMat, T, default_mutmat_path);
    }

    msg("Reading index...");
    inMessage = readCapnp(default_index_path);
    Index::Reader index_input = inMessage->getRoot<Index>();

    if (!eval.empty()) {
        msg("[Developer mode] --- Evaluate placement accuracy ---");
        std::exit(0);
    }

    // Placement phase
    msg("=== Starting Read Placement ===");
    auto start = std::chrono::high_resolution_clock::now();
    
    try {
        if (genotype_from_sam) {
            msg("Genotyping from SAM file");
            if (samFileName.empty() || refFileName.empty()) {
                throw std::runtime_error("SAM file and reference file are required for genotyping from SAM");
            }
            genotype(prefix, refFileName, "", bamFileName, mpileupFileName, vcfFileName, mutMat);
        } else if (placement_per_read) {
            throw std::runtime_error("PLACEMENT PER READ REMOVED IN THIS BRANCH");
        } else {
            if (!refNode.empty()) {
                msg("Using reference node: {}", refNode);
                if (T->allNodes.find(refNode) == T->allNodes.end()) {
                    throw std::runtime_error("Reference node not found in pangenome");
                }
            }
            
            // Initialize placement variables
            placement::PlacementResult result;
            placement::LinkedNode *maxHitsNode = nullptr, *bestJaccardNode = nullptr,
                                *bestCosineNode = nullptr;
            int64_t maxHitsInAnyGenome = 0;
            double bestJaccardScore = 0.0, bestCosineScore = 0.0;
            
            // Perform placement
            place(
                maxHitsInAnyGenome, maxHitsNode,
                bestJaccardScore, bestJaccardNode,
                bestCosineScore, bestCosineNode,
                result, T, index_input,
                reads1, reads2,
                mutMat, prefix, refFileName, samFileName,
                bamFileName, mpileupFileName, vcfFileName,
                aligner, refNode,
                save_jaccard, show_time,
                0.01f, 16,
                "", "",
                0, 0,
                hot_threshold
            );

            msg("Placement completed successfully");
            if (bestJaccardNode && bestJaccardNode->node) {
                msg("Best placement node (Jaccard): {}", bestJaccardNode->node->identifier);
                msg("Jaccard score: {}", std::to_string(bestJaccardScore));
            }
        }
    } catch (const std::exception& e) {
        err("ERROR during placement: {}", e.what());
        return 1;
    }
    
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start
    );
    msg("Placement completed in {}ms", duration.count());

    auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - main_start
    );
    msg("=== panmap execution completed ===");
    msg("Total runtime: {}ms", total_duration.count());

    if (show_time) {
        timing::Timer::report();
    }

    return 0;
}
