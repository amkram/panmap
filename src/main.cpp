/**
 * @file main.cpp
 * @brief Panmap - Pangenome-based sequence placement, alignment, and genotyping
 *
 * Supports single-sample (isolate) and metagenomic workflows.
 */

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <boost/algorithm/string.hpp>
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#include <tbb/tbb.h>
// Compatibility shim: oneTBB 2021+ uses tbb::filter_mode enum class.
// TBB 2019 uses tbb::filter::serial_in_order nested enum.
#if defined(TBB_VERSION_MAJOR) && TBB_VERSION_MAJOR >= 2021
#define TBB_FILTER_SERIAL_IN_ORDER tbb::filter_mode::serial_in_order
#define TBB_FILTER_PARALLEL tbb::filter_mode::parallel
#define TBB_FILTER_SERIAL_OUT_OF_ORDER tbb::filter_mode::serial_out_of_order
#else
#define TBB_FILTER_SERIAL_IN_ORDER tbb::filter::serial_in_order
#define TBB_FILTER_PARALLEL tbb::filter::parallel
#define TBB_FILTER_SERIAL_OUT_OF_ORDER tbb::filter::serial_out_of_order
#endif
#include <fstream>
#include <sstream>
#include <iostream>
#include <memory>
#include <random>
#include <chrono>
#include <optional>
#include <algorithm>
#include <atomic>
#include <mutex>
#include <unordered_map>
#include <unordered_set>

#include <absl/container/flat_hash_map.h>

#include "capnp/message.h"
#include "capnp/serialize.h"
#include "kj/io.h"
#include "index_lite.capnp.h"
#include "logging.hpp"
#include "panman.hpp"
#include "placement.hpp"
#include "panmap_utils.hpp"
#include "index_single_mode.hpp"
#include "zstd_compression.hpp"
#include "genotyping.hpp"
#include "conversion.hpp"
#include "seeding.hpp"
#include "mgsr.hpp"

extern "C" {
#include <htslib/sam.h>
#include <htslib/faidx.h>
}

#include "version.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

constexpr const char* VERSION = PANMAP_VERSION;
constexpr const char* PROGRAM_NAME = "panmap";

namespace color {
inline const char* reset() {
    return output::style::reset();
}

inline const char* bold() {
    return output::style::bold();
}

inline const char* dim() {
    return output::style::dim();
}

inline const char* red() {
    return output::style::red();
}

inline const char* green() {
    return output::style::green();
}

inline const char* yellow() {
    return output::style::yellow();
}

inline const char* blue() {
    return output::style::blue();
}

inline const char* cyan() {
    return output::style::cyan();
}
}  // namespace color

enum class PipelineStage {
    Index,      // Build index only
    Place,      // Placement only
    Align,      // Placement + Alignment
    Genotype,   // Placement + Alignment + Genotyping
    Consensus,  // + Consensus FASTA generation
    Full        // Full pipeline
};

struct Config {
    // Input files
    std::string panman;  // Guide pangenome (.panman)
    std::string reads1;  // First read file (FASTQ/FASTA)
    std::string reads2;  // Second read file for paired-end

    // Output
    std::string output;  // Output prefix
    std::string index;   // Index file path

    // Pipeline control
    PipelineStage stopAfter = PipelineStage::Place;
    bool forceReindex = false;

    // Mode
    bool metagenomic = false;  // Metagenomic mode (multi-sample)
    int topN = 1;              // Report top N placements
    bool dedupReads = false;   // Deduplicate reads before placement (for amplicon data)

    std::string aligner = "minimap2";

    // Index parameters
    int k = 19;                   // syncmer k
    int s = 8;                    // syncmer s
    int l = 3;                    // l-mer size
    int t = 0;                    // syncmer offset
    bool openSyncmer = false;     // Open syncmer
    int flankMaskBp = 250;        // Hard mask first/last N bp at genome ends
    double seedMaskFraction = 0;  // Mask top 0.1% most frequent seeds
    int minSeedQuality = 0;       // Min avg Phred quality for seed region (0=disabled)
    int trimStart = 0;            // Trim N bases from start of each read (primer removal)
    int trimEnd = 0;              // Trim N bases from end of each read (primer removal)
    int minReadSupport = 1;       // Min reads for a seed to be counted (2 = filter singletons)
    bool hpc = false;             // Homopolymer-compressed seeds
    bool extentGuard = false;     // Guard seed deletions at genome extent boundaries

    // Resources
    int threads = 1;
    int zstdLevel = 7;

    // Metagenomic options
    bool indexFull = false;
    bool indexPacked = false;
    bool readPacked = false;
    bool noProgress = false;
    size_t topOc = 1000;
    uint32_t maskReads = 0;
    uint32_t maskSeeds = 0;
    std::string indexMgsr;
    std::string ampliconDepth;
    double maskReadsRelativeFrequency = 0.0;
    double maskSeedsRelativeFrequency = 0.0;
    double emConvergenceThreshold = 0.00001;
    double emDeltaThreshold = 0.0;
    uint32_t emMaximumIterations = 1000;
    uint32_t emMaximumRounds = 5;
    bool emLeavesOnly = false;
    bool filterAndAssign = false;
    double dust = 100.0;
    double discard = 0.0;
    uint32_t maskReadEnds = 0;
    std::string taxonomicMetadata;
    std::string taxonomicRank;
    size_t maximumTaxonNumber = 1;
    double ambiguousScoreThresholdRatio = 0.0;
    int ambiguousScoreThreshold = 0;
    bool breadthRatio = false;
    bool pseudochain = false;
    std::string batchFilesPath;
    size_t batchSize = 1000000;

    // Utility modes
    std::string dumpNodeId;
    bool writeMetaReadScoresFiltered = false;
    bool writeMetaReadScoresUnfiltered = false;
    bool writeOCRanks = false;
    int seed = 42;

    // Batch mode
    std::string batchFile;  // Path to batch file listing samples (one per line: reads1 [reads2])

    // Leave-one-out validation mode

    // Diagnostic options
    std::string dumpAllScores;  // Dump all node scores to this file

    // Output control
    bool quiet = false;    // Minimal output (errors only)
    bool verbose = false;  // Extra debug output
    bool plain = false;    // Plain text output (no colors/unicode)

    // Leaf-only placement
    bool forceLeaf = false;  // Restrict placement to leaf nodes only

    // Consensus options
    bool impute = false;              // Impute N's from parent sequence (ignore _->N mutations)
    bool noMutationSpectrum = false;  // Skip mutation spectrum filtering in VCF
    bool baq = false;                 // Enable BAQ (Base Alignment Quality) in mpileup

    // Pre-computed substitution spectrum from index (4x4 phred-scaled matrix)
    // Empty if index has no spectrum or --no-mutation-spectrum is set
    std::vector<std::vector<double>> substMatrixPhred;

    // Alignment-based refinement options
    bool refine = false;           // Enable alignment-based refinement
    double refineTopPct = 0.01;    // Top X% of nodes to refine (default 1%)
    int refineMaxTopN = 150;       // Max nodes to align against
    int refineNeighborRadius = 2;  // Expand to neighbors within N branches
    int refineMaxNeighborN = 150;  // Max additional nodes from neighbor expansion
};

// Cap'n Proto reader with ZSTD decompression
class IndexReader : public ::capnp::MessageReader {
   public:
    std::vector<uint8_t> data;
    std::unique_ptr<::capnp::FlatArrayMessageReader> reader;

    explicit IndexReader(const std::string& path, int numThreads = 0) : ::capnp::MessageReader(makeOptions()) {
        if (!panmap_zstd::decompressFromFile(path, data, numThreads)) {
            throw std::runtime_error("Failed to decompress index: " + path);
        }

        reader = std::make_unique<::capnp::FlatArrayMessageReader>(
            kj::ArrayPtr<const capnp::word>(reinterpret_cast<const capnp::word*>(data.data()),
                                            data.size() / sizeof(capnp::word)),
            makeOptions());
    }

    kj::ArrayPtr<const capnp::word> getSegment(uint id) override { return reader->getSegment(id); }

   private:
    static ::capnp::ReaderOptions makeOptions() {
        ::capnp::ReaderOptions opts;
        opts.traversalLimitInWords = kj::maxValue;
        opts.nestingLimit = 1024;
        return opts;
    }
};

class PackedFdReader {
   public:
    int fd = -1;
    std::unique_ptr<::capnp::PackedFdMessageReader> reader;

    explicit PackedFdReader(const std::string& path) {
        ::capnp::ReaderOptions opts;
        opts.traversalLimitInWords = kj::maxValue;
        opts.nestingLimit = 1024;

        fd = mgsr::open_file(path);
        reader = std::make_unique<::capnp::PackedFdMessageReader>(fd, opts);
    }

    ~PackedFdReader() {
        reader.reset();
        if (fd >= 0) {
            close(fd);
        }
    }

    template <typename T>
    typename T::Reader getRoot() {
        return reader->getRoot<T>();
    }
};

class FdReader {
   public:
    int fd = -1;
    std::unique_ptr<::capnp::StreamFdMessageReader> reader;

    explicit FdReader(const std::string& path) {
        ::capnp::ReaderOptions opts;
        opts.traversalLimitInWords = 16ULL * 1024 * 1024 * 1024;
        opts.nestingLimit = 1024;
        fd = mgsr::open_file(path);
        reader = std::make_unique<::capnp::StreamFdMessageReader>(fd, opts);
    }

    ~FdReader() {
        reader.reset();
        if (fd >= 0) close(fd);
    }

    template <typename T>
    typename T::Reader getRoot() {
        return reader->getRoot<T>();
    }
};

// Load 4x4 substitution matrix from index and convert to phred scale.
// Returns empty matrix if index has no spectrum data.
std::vector<std::vector<double>> loadSubstMatrixFromIndex(LiteIndex::Reader& idx) {
    auto matReader = idx.getSubstitutionMatrix();
    if (matReader.size() != 16) return {};

    // Check if matrix is populated (non-identity)
    bool allZeroOffDiag = true;
    for (int i = 0; i < 4 && allZeroOffDiag; i++)
        for (int j = 0; j < 4; j++)
            if (i != j && matReader[i * 4 + j] > 0) {
                allZeroOffDiag = false;
                break;
            }
    if (allZeroOffDiag) return {};

    // Convert probability matrix to phred: -10 * log10(p)
    std::vector<std::vector<double>> phred(4, std::vector<double>(4, 0.0));
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {
            double p = matReader[i * 4 + j];
            phred[i][j] = (p > 0) ? -10.0 * log10(p) : 100.0;
        }
    return phred;
}

std::unique_ptr<panmanUtils::TreeGroup> loadPanMAN(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Cannot open PanMAN file: " + path);
    }

    auto buffer = std::make_unique<boost::iostreams::filtering_streambuf<boost::iostreams::input>>();
    buffer->push(boost::iostreams::lzma_decompressor());
    buffer->push(file);

    std::istream stream(buffer.get());
    return std::make_unique<panmanUtils::TreeGroup>(stream);
}

/**
 * Get ungapped genome length for a node (count non-gap characters).
 */

std::string sanitizeFilename(const std::string& s) {
    std::string result = s;
    for (char& c : result) {
        if (c == '/' || c == '|' || c == ' ' || c == ':' || c == '\\') c = '_';
    }
    return result;
}

void saveNodeSequence(panmanUtils::Tree* T, const std::string& nodeId, const std::string& path) {
    std::string seq = T->getStringFromReference(nodeId, false, true);

    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write: " + path);

    out << ">" << nodeId << "\n";
    for (size_t i = 0; i < seq.size(); i += 80) {
        out << seq.substr(i, 80) << "\n";
    }
    logging::msg("Saved {} ({} bp) to {}", nodeId, seq.size(), path);
}

/* INDEXING */

bool buildMgsrIndex(const Config& cfg) {
    auto tg = loadPanMAN(cfg.panman);
    if (!tg || tg->trees.empty()) {
        logging::err("Failed to load pangenome");
        return false;
    }

    panmanUtils::Tree* T = &tg->trees[0];
    mgsr::mgsrIndexBuilder mgsrIndexBuilder(T, cfg.k, cfg.s, cfg.t, cfg.l, cfg.openSyncmer, cfg.impute, cfg.indexFull);
    mgsrIndexBuilder.buildIndex();
    mgsrIndexBuilder.writeIndex(cfg.indexMgsr, cfg.indexPacked);
    std::cout << (cfg.indexFull ? "Full" : "Lite") << " MGSR index written to " << cfg.indexMgsr << std::endl;
    return true;
}

bool buildIndex(const Config& cfg) {
    if (fs::exists(cfg.index) && !cfg.forceReindex) {
        logging::msg("Using existing index: {}", cfg.index);
        return true;
    }

    logging::msg("Building index: {}", cfg.index);
    auto tg = loadPanMAN(cfg.panman);
    if (!tg || tg->trees.empty()) {
        logging::err("Failed to load pangenome");
        return false;
    }

    index_single_mode::IndexBuilder builder(
        &tg->trees[0], cfg.k, cfg.s, 0, cfg.l, false, cfg.flankMaskBp, cfg.hpc, cfg.impute, cfg.extentGuard);
    builder.buildIndexParallel(cfg.threads);
    builder.computeSubstitutionSpectrum();
    builder.writeIndex(cfg.index, cfg.threads, cfg.zstdLevel);

    logging::msg("Index built with k={}, s={}, l={}, flankMask={}bp{}{}{}",
                 cfg.k,
                 cfg.s,
                 cfg.l,
                 cfg.flankMaskBp,
                 cfg.hpc ? ", hpc=on" : "",
                 cfg.impute ? ", impute=on" : "",
                 cfg.extentGuard ? ", extentGuard=on" : "");
    return true;
}

// Compute tree distance between two nodes (number of edges)

void writeOCRanks(const std::string& outputFile,
                  const std::vector<std::pair<std::string, double>>& overlapCoefficients) {
    std::ofstream outFile(outputFile);
    uint32_t rank = 0;
    double currentOverlapCoefficient = overlapCoefficients[0].second;
    for (const auto& [nodeId, overlapCoefficient] : overlapCoefficients) {
        if (overlapCoefficient != currentOverlapCoefficient) {
            currentOverlapCoefficient = overlapCoefficient;
            ++rank;
        }
        outFile << nodeId << "\t" << std::fixed << std::setprecision(6) << overlapCoefficient << "\t" << rank
                << std::endl;
    }
    outFile.close();
}

void writeMetaReadScores(const std::string& outputFile,
                         const mgsr::ThreadsManager& threadsManager,
                         bool includeOverMaxTaxonNum) {
    std::ofstream outFile(outputFile);
    outFile << "ReadIndex\tNumDuplicates\tTotalScore\tMaxScore\tNumMaxScoreNodes\t";
    if (includeOverMaxTaxonNum) outFile << "OvermaximumTaxonNumber\t";
    outFile << "RawReadsIndices" << std::endl;
    for (size_t i = 0; i < threadsManager.reads.size(); ++i) {
        const auto& curRead = threadsManager.reads[i];
        if (curRead.maxScore == 0) continue;
        outFile << i << "\t" << threadsManager.readSeedmersDuplicatesIndex[i].size() << "\t"
                << curRead.seedmersList.size() << "\t" << curRead.maxScore << "\t" << curRead.epp << "\t";
        if (includeOverMaxTaxonNum) outFile << curRead.overMaximumTaxonNumber << "\t";
        for (size_t j = 0; j < threadsManager.readSeedmersDuplicatesIndex[i].size(); ++j) {
            if (j > 0) outFile << ",";
            outFile << threadsManager.readSeedmersDuplicatesIndex[i][j];
        }
        outFile << std::endl;
    }
    outFile.close();
}

void scoreReadsMultiThreaded(mgsr::MgsrLiteTree& T, mgsr::ThreadsManager& threadsManager, const Config& cfg) {
    std::vector<uint64_t> totalNodesPerThread(threadsManager.numThreads, 0);
    for (size_t i = 0; i < threadsManager.numThreads; ++i) {
        totalNodesPerThread[i] = T.getNumActiveNodes();
    }
    ProgressTracker progressTracker(threadsManager.numThreads, totalNodesPerThread);
    std::cout << "Placing reads with " << threadsManager.numThreads << " threads..." << std::endl;

    bool lowMemory = false;
    auto start_time_place = std::chrono::high_resolution_clock::now();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, threadsManager.threadRanges.size()),
                      [&](const tbb::blocked_range<size_t>& rangeIndex) {
                          for (size_t i = rangeIndex.begin(); i != rangeIndex.end(); ++i) {
                              auto [start, end] = threadsManager.threadRanges[i];

                              std::span<mgsr::Read> curThreadReads(threadsManager.reads.data() + start, end - start);
                              mgsr::mgsrPlacer curThreadPlacer(&T, threadsManager, lowMemory, i);
                              curThreadPlacer.initializeQueryData(curThreadReads);

                              curThreadPlacer.setAllSeedmerHashesSet(threadsManager.allSeedmerHashesSet);

                              curThreadPlacer.setProgressTracker(&progressTracker, i);
                              if (cfg.pseudochain) {
                                  curThreadPlacer.placeReads();
                              } else {
                                  curThreadPlacer.scoreReads();
                              }

                              if (i == 0) {
                                  threadsManager.identicalGroups = std::move(curThreadPlacer.identicalGroups);
                                  threadsManager.identicalNodeToGroup = std::move(curThreadPlacer.identicalNodeToGroup);
                              }
                          }
                      });
    auto end_time_place = std::chrono::high_resolution_clock::now();
    auto duration_place = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_place - start_time_place);
    std::cout << "\n\nPlaced reads in " << static_cast<double>(duration_place.count()) / 1000.0 << "s\n" << std::endl;
}

bool isFileReadable(const std::string& path) {
    std::ifstream file(path);
    return file.good();
}

void validateInputFile(const std::string& path, const std::string& description) {
    if (path.empty()) throw std::runtime_error(description + " path is empty");
    if (!fs::exists(path)) throw std::runtime_error(description + " not found: " + path);
    if (!isFileReadable(path)) throw std::runtime_error("Cannot read " + description + ": " + path);
}

void filterAndAssignBatch(mgsr::ThreadsManager& threadsManager,
                          mgsr::MgsrLiteTree& T,
                          const std::string& readPath1,
                          const std::string& readPath2,
                          double dustThreshold,
                          double discardThreshold,
                          uint32_t maskReadsEnds,
                          int ambiguousScoreThreshold,
                          double ambiguousScoreThresholdRatio,
                          size_t maximumTaxonNumber,
                          double breadthRatio,
                          size_t batchSize,
                          const std::string& prefix) {
    size_t maxLiveTokens = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism) + 2;
    mgsr::FastqFile fq1(readPath1);
    kseq_t* seq1 = kseq_init(fileno(fq1.fp));

    std::unique_ptr<mgsr::FastqFile> fq2;
    kseq_t* seq2 = nullptr;
    bool pairedEnd = !readPath2.empty();
    if (pairedEnd) {
        fq2 = std::make_unique<mgsr::FastqFile>(readPath2);
        seq2 = kseq_init(fileno(fq2->fp));
    }

    std::ofstream assignedReadsFastq(prefix + ".mgsr.assignedReads.fastq");
    if (!assignedReadsFastq.is_open()) {
        logging::err("Failed to open assigned reads fastq file: {}", prefix + ".mgsr.assignedReads.fastq");
        return;
    }

    std::unordered_map<mgsr::MgsrLiteNode*, std::vector<size_t>> allAssignedReadsByNode;
    std::unordered_map<mgsr::MgsrLiteNode*, std::vector<size_t>> allAssignedReadsByLCANode;

    struct Batch {
        size_t batchIndex;
        std::vector<std::string> sequences;
        std::vector<std::string> names;
        std::vector<std::string> quals;
        std::vector<mgsr::Read> reads;
        std::vector<std::vector<size_t>> readDuplicatesIndex;

        std::unordered_map<mgsr::MgsrLiteNode*, std::vector<size_t>> assignedReadsByNode;
        std::unordered_map<mgsr::MgsrLiteNode*, std::vector<size_t>> assignedReadsByLCANode;

        size_t numReadsAssigned;

        size_t numReadsUnmapped;
        size_t numReadsDiscarded;
        size_t numReadsOvermaximumTaxonNumber;

        double readProcessingTime;
        double scoringTime;
        double assigningTime;
        double postProcessingTime;
    };

    size_t batchesCompleted = 0;
    size_t readsCompleted = 0;
    size_t totalReadsUnmapped = 0;
    size_t totalReadsDiscarded = 0;
    size_t totalFastqReadsWritten = 0;

    size_t batchIndex = 0;
    tbb::parallel_pipeline(
        maxLiveTokens,
        tbb::make_filter<void, Batch*>(TBB_FILTER_SERIAL_IN_ORDER,
                                       [&](tbb::flow_control& fc) -> Batch* {
                                           Batch* batch = new Batch();
                                           batch->batchIndex = batchIndex++;
                                           batch->sequences.reserve(batchSize + 1);
                                           batch->names.reserve(batchSize + 1);
                                           batch->quals.reserve(batchSize + 1);
                                           int l1;
                                           if (pairedEnd) {
                                               int l2;
                                               while ((l1 = kseq_read(seq1)) >= 0 && (l2 = kseq_read(seq2)) >= 0) {
                                                   batch->sequences.emplace_back(seq1->seq.s);
                                                   batch->names.emplace_back(seq1->name.s);
                                                   batch->quals.emplace_back(seq1->qual.s);
                                                   batch->sequences.emplace_back(seq2->seq.s);
                                                   batch->names.emplace_back(seq2->name.s);
                                                   batch->quals.emplace_back(seq2->qual.s);
                                                   if (batch->sequences.size() >= batchSize) break;
                                               }
                                           } else {
                                               while ((l1 = kseq_read(seq1)) >= 0) {
                                                   batch->sequences.emplace_back(seq1->seq.s);
                                                   batch->names.emplace_back(seq1->name.s);
                                                   batch->quals.emplace_back(seq1->qual.s);
                                                   if (batch->sequences.size() >= batchSize) break;
                                               }
                                           }
                                           if (batch->sequences.empty()) {
                                               delete batch;
                                               fc.stop();
                                               return nullptr;
                                           }
                                           batch->sequences.shrink_to_fit();
                                           batch->names.shrink_to_fit();
                                           batch->quals.shrink_to_fit();
                                           return batch;
                                       }) &

            tbb::make_filter<Batch*, Batch*>(
                TBB_FILTER_PARALLEL,
                [&](Batch* batch) -> Batch* {
                    mgsr::mgsrPlacer placer(&T, threadsManager, false, 0);
                    placer.readScoreDeltasBatch.resize(T.getNumActiveNodes());

                    auto startReadProcessing = std::chrono::high_resolution_clock::now();
                    placer.initializeQueryDataBatch(batch->sequences,
                                                    batch->reads,
                                                    batch->readDuplicatesIndex,
                                                    dustThreshold,
                                                    discardThreshold,
                                                    maskReadsEnds);
                    auto endReadProcessing = std::chrono::high_resolution_clock::now();

                    auto startScoring = std::chrono::high_resolution_clock::now();
                    placer.scoreReadsBatch(discardThreshold);
                    auto endScoring = std::chrono::high_resolution_clock::now();

                    size_t& numDiscarded = batch->numReadsDiscarded;
                    size_t& numUnmapped = batch->numReadsUnmapped;
                    for (size_t i = 0; i < batch->reads.size(); ++i) {
                        auto& read = batch->reads[i];
                        if (read.maxScore == 0) {
                            numUnmapped += batch->readDuplicatesIndex[i].size();
                        } else if (read.maxScore < read.discardThreshold) {
                            read.maxScore = 0;
                            numDiscarded += batch->readDuplicatesIndex[i].size();
                        }
                    }
                    auto startAssigning = std::chrono::high_resolution_clock::now();
                    placer.assignReadsBatch(batch->assignedReadsByNode,
                                            batch->assignedReadsByLCANode,
                                            maximumTaxonNumber,
                                            ambiguousScoreThreshold,
                                            ambiguousScoreThresholdRatio);

                    auto endAssigning = std::chrono::high_resolution_clock::now();

                    batch->readProcessingTime =
                        std::chrono::duration<double>(endReadProcessing - startReadProcessing).count();
                    batch->scoringTime = std::chrono::duration<double>(endScoring - startScoring).count();
                    batch->assigningTime = std::chrono::duration<double>(endAssigning - startAssigning).count();
                    return batch;
                }) &

            tbb::make_filter<Batch*, void>(TBB_FILTER_SERIAL_OUT_OF_ORDER, [&](Batch* batch) {
                auto startPostProcessing = std::chrono::high_resolution_clock::now();
                readsCompleted += batch->sequences.size();
                totalReadsUnmapped += batch->numReadsUnmapped;
                totalReadsDiscarded += batch->numReadsDiscarded;
                ++batchesCompleted;

                std::unordered_map<size_t, size_t> batchReadIndexToFastqIndex;
                for (auto& [node, uniqueReads] : batch->assignedReadsByNode) {
                    for (size_t uniqueReadIndex : uniqueReads) {
                        for (size_t readIndex : batch->readDuplicatesIndex[uniqueReadIndex]) {
                            auto [it, inserted] = batchReadIndexToFastqIndex.emplace(readIndex, totalFastqReadsWritten);
                            if (inserted) {
                                assignedReadsFastq << "@" << batch->names[readIndex] << "\n";
                                assignedReadsFastq << batch->sequences[readIndex] << "\n";
                                assignedReadsFastq << "+\n";
                                assignedReadsFastq << batch->quals[readIndex] << "\n";
                                totalFastqReadsWritten++;
                            }
                            allAssignedReadsByNode[node].push_back(it->second);
                        }
                    }
                }

                for (auto& [node, uniqueReads] : batch->assignedReadsByLCANode) {
                    for (size_t uniqueReadIndex : uniqueReads) {
                        for (size_t readIndex : batch->readDuplicatesIndex[uniqueReadIndex]) {
                            auto [it, inserted] = batchReadIndexToFastqIndex.emplace(readIndex, totalFastqReadsWritten);
                            if (inserted) {
                                std::cerr << "Read index " << readIndex << " not found in batchReadIndexToFastqIndex"
                                          << std::endl;
                                exit(1);
                            }
                            allAssignedReadsByLCANode[node].push_back(it->second);
                        }
                    }
                }
                batch->numReadsAssigned = batchReadIndexToFastqIndex.size();
                auto endPostProcessing = std::chrono::high_resolution_clock::now();
                batch->postProcessingTime =
                    std::chrono::duration<double>(endPostProcessing - startPostProcessing).count();

                std::cout << readsCompleted << " total reads processed | " << batch->sequences.size()
                          << " reads processed | " << batch->numReadsAssigned
                          << " reads assigned and written to fastq file | " << batch->numReadsUnmapped
                          << " reads unmapped | " << batch->numReadsDiscarded << " reads discarded | "
                          << batch->numReadsOvermaximumTaxonNumber << " reads over maximum taxon number | "
                          << batch->readProcessingTime << " seconds (read processing) | " << batch->scoringTime
                          << " seconds (scoring) | " << batch->assigningTime << " seconds (assigning) | "
                          << batch->postProcessingTime << " seconds (post-processing)\n";
                // collect results ...
                delete batch;
            }));

    kseq_destroy(seq1);
    if (pairedEnd) {
        kseq_destroy(seq2);
    }
    assignedReadsFastq.close();

    std::ofstream assignedReadsOut(prefix + ".mgsr.assignedReads.out");
    if (!assignedReadsOut.is_open()) {
        logging::err("Failed to open assigned reads output file: {}", prefix + ".mgsr.assignedReads.out");
        return;
    }
    for (auto& [node, readIndices] : allAssignedReadsByNode) {
        assignedReadsOut << node->identifier;
        for (const auto& identicalNodeId : node->identicalNodeIdentifiers) {
            assignedReadsOut << "," << identicalNodeId;
        }

        std::sort(readIndices.begin(), readIndices.end());
        assignedReadsOut << "\t" << readIndices.size() << "\t";

        for (size_t i = 0; i < readIndices.size(); ++i) {
            assignedReadsOut << readIndices[i];
            if (i != readIndices.size() - 1) {
                assignedReadsOut << ",";
            }
        }
        assignedReadsOut << "\n";
    }
    assignedReadsOut.close();
    std::cout << totalFastqReadsWritten << " reads written to fastq file" << std::endl;

    std::ofstream assignedReadsLCANodeOut(prefix + ".mgsr.assignedReadsLCANode.out");
    if (!assignedReadsLCANodeOut.is_open()) {
        logging::err("Failed to open assigned reads LCA node output file: {}",
                     prefix + ".mgsr.assignedReadsLCANode.out");
        return;
    }
    for (auto& [node, readIndices] : allAssignedReadsByLCANode) {
        assignedReadsLCANodeOut << node->identifier;
        for (const auto& identicalNodeId : node->identicalNodeIdentifiers) {
            assignedReadsLCANodeOut << "," << identicalNodeId;
        }

        std::sort(readIndices.begin(), readIndices.end());
        assignedReadsLCANodeOut << "\t" << readIndices.size() << "\t";

        for (size_t i = 0; i < readIndices.size(); ++i) {
            assignedReadsLCANodeOut << readIndices[i];
            if (i != readIndices.size() - 1) {
                assignedReadsLCANodeOut << ",";
            }
        }
        assignedReadsLCANodeOut << "\n";
    }

    assignedReadsLCANodeOut.close();

    // Calculate breadth from assigned reads
    if (breadthRatio) {
        std::cout << "Calculating breadth ratio from assigned reads..." << std::endl;
        FILE* fpAssigned = fopen(std::string(prefix + ".mgsr.assignedReads.fastq").c_str(), "r");
        kseq_t* seqAssigned = kseq_init(fileno(fpAssigned));
        std::vector<std::string> assignedReads;
        assignedReads.reserve(totalFastqReadsWritten);
        int l;
        while ((l = kseq_read(seqAssigned)) >= 0) {
            assignedReads.emplace_back(seqAssigned->seq.s);
        }
        kseq_destroy(seqAssigned);
        fclose(fpAssigned);
        std::cout << "Read in " << assignedReads.size() << " assigned reads" << std::endl;
        if (!assignedReads.empty()) {
            mgsr::mgsrPlacer placer(&T, threadsManager, false, 0);
            placer.readScoreDeltasBatch.resize(T.getNumActiveNodes());

            std::vector<mgsr::Read> reads;
            std::vector<std::vector<size_t>> readDuplicatesIndex;
            placer.initializeQueryDataBatch(
                assignedReads, reads, readDuplicatesIndex, dustThreshold, discardThreshold, maskReadsEnds);

            placer.scoreReadsBatch(discardThreshold);

            std::unordered_map<mgsr::MgsrLiteNode*, std::vector<size_t>> assignedReadsByNode;
            std::unordered_map<mgsr::MgsrLiteNode*, std::vector<size_t>> assignedReadsByLCANode;
            placer.assignReadsBatch(assignedReadsByNode,
                                    assignedReadsByLCANode,
                                    maximumTaxonNumber,
                                    ambiguousScoreThreshold,
                                    ambiguousScoreThresholdRatio);
            std::unordered_map<size_t, int64_t> refSeedsCount;
            std::vector<std::pair<mgsr::MgsrLiteNode*, mgsr::breadthInfo>> breadths;
            placer.calculateBreadthRatio(T.root, refSeedsCount, breadths, assignedReadsByNode);

            std::ofstream breadthsOut(prefix + ".mgsr.breadths.out");
            breadthsOut << "NodeId\tTotalRefSeeds\tObservedBreadthCount\tObservedBreadthRatio\tTotalDepth\tMeanDepth\tE"
                           "xpectedBreadthRatio\tObservedToExpectedBreadthRatio"
                        << std::endl;
            for (const auto& [node, breadthInfo] : breadths) {
                breadthsOut << node->identifier;
                for (const auto& identicalNodeId : node->identicalNodeIdentifiers) {
                    breadthsOut << "," << identicalNodeId;
                }
                breadthsOut << "\t" << breadthInfo.totalRefSeeds << "\t" << breadthInfo.observedBreadthCount << "\t"
                            << breadthInfo.observedBreadthRatio << "\t" << breadthInfo.totalDepth << "\t"
                            << breadthInfo.meanDepth << "\t" << breadthInfo.expectedBreadthRatio << "\t"
                            << breadthInfo.observedToExpectedBreadthRatio << std::endl;
            }
            breadthsOut.close();
        } else {
            std::cout << "No assigned reads found" << std::endl;
            std::ofstream breadthsOut(prefix + ".mgsr.breadths.out");
            breadthsOut << "NodeId\tTotalRefSeeds\tObservedBreadthCount\tObservedBreadthRatio\tTotalDepth\tMeanDepth\tE"
                           "xpectedBreadthRatio\tObservedToExpectedBreadthRatio"
                        << std::endl;
            breadthsOut.close();
        }
    }
}

struct BatchEntry {
  std::string reads1;
  std::string reads2;
  std::string prefix;
};

bool readBatchFiles(const std::string& batchFilesPath, std::vector<BatchEntry>& entries) {
  validateInputFile(batchFilesPath, "batch files TSV");
  std::ifstream batchIn(batchFilesPath);
  if (!batchIn.is_open()) {
    logging::err("Cannot open batch files TSV: {}", batchFilesPath);
    exit(1);
  }
  std::string line;
  size_t lineNum = 0;
  while (std::getline(batchIn, line)) {
    ++lineNum;
    boost::trim(line);
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::vector<std::string> cols;
    boost::split(cols, line, boost::is_any_of("\t"), boost::token_compress_on);
    for (auto& c : cols) {
      boost::trim(c);
    }
    while (!cols.empty() && cols.back().empty()) {
      cols.pop_back();
    }
    if (cols.size() == 2) {
      validateInputFile(cols[0], "reads1 (batch TSV)");
      entries.push_back({cols[0], "", cols[1]});
    } else if (cols.size() == 3) {
      validateInputFile(cols[0], "reads1 (batch TSV)");
      validateInputFile(cols[1], "reads2 (batch TSV)");
      entries.push_back({cols[0], cols[1], cols[2]});
    } else {
      logging::err("Batch TSV {} line {}: expected 2 or 3 tab-separated columns, got {}",
                   batchFilesPath,
                   lineNum,
                   cols.size());
      return false;
    }
  }
  return true;
}

bool runFilterAndAssign(mgsr::MgsrLiteTree& T, mgsr::ThreadsManager& threadsManager, const Config& cfg) {
  auto start_time_filterAndAssign = std::chrono::high_resolution_clock::now();

  auto start_time_buildEulerTour = std::chrono::high_resolution_clock::now();
  T.buildEulerTour();
  auto end_time_buildEulerTour = std::chrono::high_resolution_clock::now();
  auto duration_buildEulerTour =
      std::chrono::duration_cast<std::chrono::milliseconds>(end_time_buildEulerTour - start_time_buildEulerTour);
  std::cout << "Built Euler tour in " << static_cast<double>(duration_buildEulerTour.count()) / 1000.0 << "s\n"
            << std::endl;

  if (cfg.breadthRatio) {
    auto start_time_calculateRefSeedCounts = std::chrono::high_resolution_clock::now();
    T.calculateRefSeedCounts();
    auto end_time_calculateRefSeedCounts = std::chrono::high_resolution_clock::now();
    auto duration_calculateRefSeedCounts = std::chrono::duration_cast<std::chrono::milliseconds>(
        end_time_calculateRefSeedCounts - start_time_calculateRefSeedCounts);
    std::cout << "Calculated reference seed counts in "
              << static_cast<double>(duration_calculateRefSeedCounts.count()) / 1000.0 << "s\n"
              << std::endl;
  }

  auto start_time_toRefSeedDeltas = std::chrono::high_resolution_clock::now();
  T.toRefSeedDeltas();
  auto end_time_toRefSeedDeltas = std::chrono::high_resolution_clock::now();
  auto duration_toRefSeedDeltas =
      std::chrono::duration_cast<std::chrono::milliseconds>(end_time_toRefSeedDeltas - start_time_toRefSeedDeltas);
  std::cout << "Converted to reference seed deltas in "
            << static_cast<double>(duration_toRefSeedDeltas.count()) / 1000.0 << "s\n"
            << std::endl;

  if (!cfg.reads1.empty()) {
    filterAndAssignBatch(threadsManager,
                         T,
                         cfg.reads1,
                         cfg.reads2,
                         cfg.dust,
                         cfg.discard,
                         cfg.maskReadEnds,
                         cfg.ambiguousScoreThreshold,
                         cfg.ambiguousScoreThresholdRatio,
                         cfg.maximumTaxonNumber,
                         cfg.breadthRatio,
                         cfg.batchSize,
                         cfg.output);
  } else if (!cfg.batchFilesPath.empty()) {
    std::vector<BatchEntry> batchEntries;
    if (!readBatchFiles(cfg.batchFilesPath, batchEntries)) {
      return false;
    }
    for (size_t i = 0; i < batchEntries.size(); ++i) {
      const auto& entry = batchEntries[i];
      if (entry.reads2.empty()) {
        logging::info("filter-and-assign batch entry {}: reads1={}, prefix={}",
                      i + 1, entry.reads1, entry.prefix);
      } else {
        logging::info("filter-and-assign batch entry {}: reads1={}, reads2={}, prefix={}",
                      i + 1, entry.reads1, entry.reads2, entry.prefix);
      }
      filterAndAssignBatch(threadsManager,
                           T,
                           entry.reads1,
                           entry.reads2,
                           cfg.dust,
                           cfg.discard,
                           cfg.maskReadEnds,
                           cfg.ambiguousScoreThreshold,
                           cfg.ambiguousScoreThresholdRatio,
                           cfg.maximumTaxonNumber,
                           cfg.breadthRatio,
                           cfg.batchSize,
                           entry.prefix);
    }
  }
  auto end_time_filterAndAssign = std::chrono::high_resolution_clock::now();
  auto duration_filterAndAssign =
      std::chrono::duration_cast<std::chrono::milliseconds>(end_time_filterAndAssign - start_time_filterAndAssign);
  std::cout << "Filter and assign took " << static_cast<double>(duration_filterAndAssign.count()) / 1000.0 << "s\n"
            << std::endl;
  return true;
}

bool runDeconvolution(mgsr::MgsrLiteTree& T, mgsr::ThreadsManager& threadsManager, const Config& cfg) {
    auto start_time_deconvolution = std::chrono::high_resolution_clock::now();

    auto start_time_initializeQueryData = std::chrono::high_resolution_clock::now();
    threadsManager.initializeQueryData(cfg.reads1, cfg.reads2, cfg.ampliconDepth, cfg.dust, cfg.maskReadEnds);
    auto end_time_initializeQueryData = std::chrono::high_resolution_clock::now();
    auto duration_initializeQueryData = std::chrono::duration_cast<std::chrono::milliseconds>(
        end_time_initializeQueryData - start_time_initializeQueryData);
    std::cout << "Initialize query data took " << static_cast<double>(duration_initializeQueryData.count()) / 1000.0
              << "s\n"
              << std::endl;

    // compute overlap coefficients
    {
        bool lowMemory = false;
        mgsr::mgsrPlacer placer(&T, threadsManager, lowMemory, 0);
        auto overlapCoefficients = placer.computeOverlapCoefficients(threadsManager.allSeedmerHashesSet);

        T.fillOCRanks(overlapCoefficients);

        if (cfg.writeOCRanks) {
            writeOCRanks(cfg.output + ".overlapCoefficients.tsv", overlapCoefficients);
        }
    }

    auto start_time_collapseIdenticalScoringNodes = std::chrono::high_resolution_clock::now();
    T.collapseIdenticalScoringNodes(threadsManager.allSeedmerHashesSet);
    auto end_time_collapseIdenticalScoringNodes = std::chrono::high_resolution_clock::now();
    auto duration_collapseIdenticalScoringNodes = std::chrono::duration_cast<std::chrono::milliseconds>(
        end_time_collapseIdenticalScoringNodes - start_time_collapseIdenticalScoringNodes);
    std::cout << "Collapsed identical scoring nodes in "
              << static_cast<double>(duration_collapseIdenticalScoringNodes.count()) / 1000.0 << "s\n"
              << std::endl;

    scoreReadsMultiThreaded(T, threadsManager, cfg);

    if (cfg.writeMetaReadScoresUnfiltered) {
        writeMetaReadScores(cfg.output + ".read_scores_info.unfiltered.tsv", threadsManager, false);
    }

    double discard_threshold = cfg.discard;
    size_t num_discarded = 0;
    size_t num_unmapped = 0;
    for (auto& read : threadsManager.reads) {
        if (read.maxScore == 0) {
            ++num_unmapped;
        } else if (read.maxScore < static_cast<int>(read.seedmersList.size() * discard_threshold)) {
            read.maxScore = 0;
            ++num_discarded;
        }
    }
    std::cout << num_unmapped << " reads unmapped... " << std::endl;
    std::cout << num_discarded << " reads discarded due to low parsimony score... " << std::endl;
    std::cout << threadsManager.reads.size() - num_unmapped - num_discarded << " reads mapped to nodes..." << std::endl;
    if (threadsManager.reads.size() - num_unmapped - num_discarded == 0) {
        std::cerr << "No reads remain for node scoring and EM after discarding low-score reads... Exiting... "
                  << std::endl;
        return true;
    }

    T.seedInfos.clear();  // no longer needed. clear memory to prep for EM.
    mgsr::squareEM squareEM(threadsManager,
                            T,
                            cfg.output,
                            cfg.topOc,
                            cfg.emConvergenceThreshold,
                            cfg.emDeltaThreshold,
                            cfg.emMaximumIterations,
                            false,
                            cfg.emLeavesOnly);

    auto start_time_squareEM = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < cfg.emMaximumRounds; ++i) {
        squareEM.runSquareEM();
        std::cout << "\nRound " << i << " of squareEM completed... nodes size changed from " << squareEM.nodes.size()
                  << " to ";
        bool removed = squareEM.removeLowPropNodes();
        std::cout << squareEM.nodes.size() << std::endl;
        if (!removed) {
            break;
        }
    }
    std::cout << std::endl;

    auto end_time_squareEM = std::chrono::high_resolution_clock::now();
    auto duration_squareEM =
        std::chrono::duration_cast<std::chrono::milliseconds>(end_time_squareEM - start_time_squareEM);
    std::cout << "SquareEM completed in " << static_cast<double>(duration_squareEM.count()) / 1000.0 << "s\n"
              << std::endl;

    std::vector<uint64_t> indices(squareEM.nodes.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&squareEM](uint64_t i, uint64_t j) {
        return squareEM.props[i] > squareEM.props[j];
    });

    std::cout << "writing abundance file: " << cfg.output + ".mgsr.abundance.out" << std::endl;
    std::ofstream abundanceOutput(cfg.output + ".mgsr.abundance.out");
    abundanceOutput << std::setprecision(5) << std::fixed;
    for (size_t i = 0; i < indices.size(); ++i) {
        size_t index = indices[i];
        abundanceOutput << squareEM.nodes[index];
        for (const auto& member : T.allLiteNodes.at(squareEM.nodes[index])->identicalNodeIdentifiers) {
            abundanceOutput << "," << member;
        }
        if (squareEM.identicalGroups.find(squareEM.nodes[index]) != squareEM.identicalGroups.end()) {
            for (const auto& member : squareEM.identicalGroups[squareEM.nodes[index]]) {
                abundanceOutput << "," << member;
                for (const auto& identicalMember : T.allLiteNodes.at(member)->identicalNodeIdentifiers) {
                    abundanceOutput << "," << identicalMember;
                }
            }
        }
        abundanceOutput << "\t" << squareEM.props[index] << std::endl;
    }
    abundanceOutput.close();

    if (cfg.writeMetaReadScoresFiltered) {
        writeMetaReadScores(cfg.output + ".read_scores_info.filtered.tsv", threadsManager, true);
    }

    return true;
}

bool runMetagenomic(const Config& cfg) {
    std::cout << "Running metagenomic mode with index: " << cfg.index << " and threads: " << cfg.threads << std::endl;

    // Checking IO
    if (cfg.index.empty() || !fs::exists(cfg.index)) {
        std::cerr << "Error: Index file " << cfg.index << " does not exist" << std::endl;
        return false;
    }

    if (!cfg.taxonomicMetadata.empty() && !fs::exists(cfg.taxonomicMetadata)) {
        std::cerr << "Error: Taxonomic metadata file " << cfg.taxonomicMetadata << " does not exist" << std::endl;
        return false;
    }


    if (cfg.reads1.empty()) {
      if (cfg.batchFilesPath.empty()) {
        std::cerr << "Error: Reads1 file is required when not using batch mode" << std::endl;
        return false;
      } else {
        // check if batch files path exists
        std::vector<BatchEntry> dummyBatchEntries;
        if (!readBatchFiles(cfg.batchFilesPath, dummyBatchEntries)) {
          std::cerr << "Error: Failed to read batch files" << std::endl;
          return false;
        }
      }
    } else {
      if (!fs::exists(cfg.reads1)) {
        std::cerr << "Error: Reads1 file " << cfg.reads1 << " does not exist" << std::endl;
        return false;
      }
    }

    if (!cfg.reads2.empty() && !fs::exists(cfg.reads2)) {
        std::cerr << "Error: Reads2 file " << cfg.reads2 << " does not exist" << std::endl;
        return false;
    }

    if (cfg.dust > 100.0) {
        std::cerr << "Error: --dust must be <= 100" << std::endl;
        return false;
    }

    if (cfg.discard < 0.0 || cfg.discard > 1.0) {
        std::cerr << "Error: --discard must be between 0 and 1" << std::endl;
        return false;
    }

    std::unique_ptr<IndexReader> compressedReader;
    std::unique_ptr<PackedFdReader> packedFdReader;
    std::unique_ptr<FdReader> fdReader;
    ::capnp::MessageReader* baseReader = nullptr;

    try {
        compressedReader = std::make_unique<IndexReader>(cfg.index, cfg.threads);
        baseReader = compressedReader.get();
    } catch (...) {
        try {
            std::cerr << "Trying to open index file as packed file..." << std::endl;
            if (cfg.readPacked) {
                packedFdReader = std::make_unique<PackedFdReader>(cfg.index);
                baseReader = packedFdReader->reader.get();
            } else {
                fdReader = std::make_unique<FdReader>(cfg.index);
                baseReader = fdReader->reader.get();
            }
        } catch (...) {
            std::cerr << "Failed to open index file" << std::endl;
            return false;
        }
    }

    LiteIndex::Reader indexReader = baseReader->getRoot<LiteIndex>();

    bool lowMemory = false;
    // initialize tree
    mgsr::MgsrLiteTree T;
    T.initialize(indexReader,
                 cfg.taxonomicMetadata,
                 cfg.taxonomicRank,
                 cfg.taxonomicMetadata.empty() ? 0 : cfg.maximumTaxonNumber,
                 cfg.threads,
                 lowMemory,
                 true);

    // initialize threads manager
    mgsr::ThreadsManager threadsManager(&T,
                                        cfg.output,
                                        cfg.threads,
                                        cfg.maskSeeds,
                                        cfg.maskReads,
                                        cfg.maskSeedsRelativeFrequency,
                                        cfg.maskReadsRelativeFrequency,
                                        !cfg.noProgress,
                                        lowMemory);
    threadsManager.initializeMGSRIndex(T.k, T.s, T.t, T.l, T.openSyncmer);

    // filterAndAssign
    if (cfg.filterAndAssign) {
        if (!runFilterAndAssign(T, threadsManager, cfg)) {
            return false;
        }
        return true;  // great success!
    } else {
        if (!runDeconvolution(T, threadsManager, cfg)) {
            return false;
        }
        return true;  // great success!
    }

    return true;
}

int runAlignment(const Config& cfg,
                 const placement::PlacementResult& placement,
                 panmanUtils::Tree* preloadedTree = nullptr);
int runGenotyping(const Config& cfg);
int runConsensus(const Config& cfg);

int runBatchPlacement(const Config& cfg) {
    // Parse batch file
    std::ifstream batchIn(cfg.batchFile);
    if (!batchIn.is_open()) {
        output::error("Cannot open batch file: {}", cfg.batchFile);
        return 1;
    }

    struct BatchSample {
        std::string reads1, reads2, outputPrefix;
    };

    std::vector<BatchSample> samples;
    std::string line;
    int lineNum = 0;
    while (std::getline(batchIn, line)) {
        lineNum++;
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        BatchSample s;
        iss >> s.reads1;
        if (s.reads1.empty()) continue;
        // Optional second field (reads2 or output prefix)
        std::string field2, field3;
        iss >> field2 >> field3;
        if (!field2.empty()) {
            // If field3 exists: field2=reads2, field3=output
            // If only field2: check if it looks like a fastq (reads2) or output prefix
            if (!field3.empty()) {
                s.reads2 = field2;
                s.outputPrefix = field3;
            } else {
                // Heuristic: if it ends with .fastq/.fq/.gz, it's reads2
                std::string lower = field2;
                std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
                if (lower.find(".fastq") != std::string::npos || lower.find(".fq") != std::string::npos) {
                    s.reads2 = field2;
                } else {
                    s.outputPrefix = field2;
                }
            }
        }
        // Auto-derive output prefix if not specified
        if (s.outputPrefix.empty()) {
            fs::path p(s.reads1);
            std::string stem = p.stem().string();
            for (const auto& suffix : {"_R1", "_R2", "_1", "_2", ".R1", ".R2"}) {
                if (stem.size() > strlen(suffix) && stem.substr(stem.size() - strlen(suffix)) == suffix) {
                    stem = stem.substr(0, stem.size() - strlen(suffix));
                    break;
                }
            }
            for (const auto& ext : {".fastq", ".fq"}) {
                if (stem.size() > strlen(ext) && stem.substr(stem.size() - strlen(ext)) == ext) {
                    stem = stem.substr(0, stem.size() - strlen(ext));
                    break;
                }
            }
            s.outputPrefix = (p.parent_path() / stem).string();
        }
        if (!fs::exists(s.reads1)) {
            output::error("Batch line {}: reads file not found: {}", lineNum, s.reads1);
            return 1;
        }
        if (!s.reads2.empty() && !fs::exists(s.reads2)) {
            output::error("Batch line {}: reads file not found: {}", lineNum, s.reads2);
            return 1;
        }
        samples.push_back(std::move(s));
    }
    batchIn.close();

    if (samples.empty()) {
        output::error("No samples found in batch file");
        return 1;
    }

    output::info("Batch mode: {} samples", samples.size());

    // Load index ONCE
    logging::msg("Loading index...");
    IndexReader reader(cfg.index, cfg.threads);
    auto idx = reader.getRoot<LiteIndex>();
    logging::debug("Index parameters: k={}, s={}, l={}", idx.getK(), idx.getS(), idx.getL());

    panmapUtils::LiteTree tree;
    tree.initialize(idx.getLiteTree());
    logging::msg("Tree loaded: {} nodes", tree.allLiteNodes.size());

    // Load full tree if refinement is enabled or if running alignment/genotyping
    std::unique_ptr<panmanUtils::TreeGroup> tg;
    panmanUtils::Tree* fullTreePtr = nullptr;
    bool refineEnabled = cfg.refine;
    bool needFullTree = refineEnabled || cfg.stopAfter >= PipelineStage::Align;
    if (needFullTree) {
        tg = loadPanMAN(cfg.panman);
        if (tg && !tg->trees.empty()) {
            fullTreePtr = &tg->trees[0];
        } else {
            logging::warn("Failed to load full tree");
            if (refineEnabled) refineEnabled = false;
        }
    }

    int successCount = 0, failCount = 0;
    auto batchStart = std::chrono::high_resolution_clock::now();

    placement::TraversalParams batchParams;
    batchParams.seedMaskFraction = cfg.seedMaskFraction;
    batchParams.minSeedQuality = cfg.minSeedQuality;
    batchParams.dedupReads = cfg.dedupReads;
    batchParams.trimStart = cfg.trimStart;
    batchParams.trimEnd = cfg.trimEnd;
    batchParams.minReadSupport = cfg.minReadSupport;
    batchParams.refineEnabled = refineEnabled;
    batchParams.refineTopPct = cfg.refineTopPct;
    batchParams.refineMaxTopN = cfg.refineMaxTopN;
    batchParams.refineNeighborRadius = cfg.refineNeighborRadius;
    batchParams.refineMaxNeighborN = cfg.refineMaxNeighborN;
    batchParams.forceLeaf = cfg.forceLeaf;

    // Pre-load seed changes by running placement on the first sample.
    {
        const auto& s = samples[0];
        std::string outPath = s.outputPrefix + ".placement.tsv";
        fs::path outDir = fs::path(s.outputPrefix).parent_path();
        if (!outDir.empty()) fs::create_directories(outDir);

        output::config().quiet = true;
        placement::PlacementResult result;
        auto start = std::chrono::high_resolution_clock::now();
        placement::placeLite(result, &tree, reader, s.reads1, s.reads2, outPath, batchParams, fullTreePtr);
        output::config().quiet = false;

        if (result.bestLogContainmentNodeId.empty()) {
            auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - start);
            fmt::print(stderr, "[1/{}] {} -> NO PLACEMENT ({}ms)\n", samples.size(), s.outputPrefix, elapsed.count());
            failCount++;
        } else {
            if (cfg.stopAfter >= PipelineStage::Align) {
                Config sampleCfg = cfg;
                sampleCfg.output = s.outputPrefix;
                sampleCfg.reads1 = s.reads1;
                sampleCfg.reads2 = s.reads2;
                if (runAlignment(sampleCfg, result, fullTreePtr) != 0) {
                    fmt::print(stderr, "[1/{}] {} -> alignment failed\n", samples.size(), s.outputPrefix);
                    failCount++;
                } else if (cfg.stopAfter >= PipelineStage::Genotype) {
                    if (runGenotyping(sampleCfg) != 0) {
                        fmt::print(stderr, "[1/{}] {} -> genotyping failed\n", samples.size(), s.outputPrefix);
                        failCount++;
                    } else {
                        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
                            std::chrono::high_resolution_clock::now() - start);
                        fmt::print(stderr,
                                   "[1/{}] {} -> {} ({}ms)\n",
                                   samples.size(),
                                   s.outputPrefix,
                                   result.bestLogContainmentNodeId,
                                   elapsed.count());
                        successCount++;
                    }
                } else {
                    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::high_resolution_clock::now() - start);
                    fmt::print(stderr,
                               "[1/{}] {} -> {} ({}ms)\n",
                               samples.size(),
                               s.outputPrefix,
                               result.bestLogContainmentNodeId,
                               elapsed.count());
                    successCount++;
                }
            } else {
                auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::high_resolution_clock::now() - start);
                fmt::print(stderr,
                           "[1/{}] {} -> {} ({}ms)\n",
                           samples.size(),
                           s.outputPrefix,
                           result.bestLogContainmentNodeId,
                           elapsed.count());
                successCount++;
            }
        }
    }

    // Process remaining samples in parallel (seed changes already loaded, tree is read-only)
    if (samples.size() > 1) {
        std::atomic<int> atomicSuccess(0), atomicFail(0);
        std::atomic<int> completedCount(1);  // first sample already done
        std::mutex progressMutex;
        output::config().quiet = true;

        tbb::parallel_for(tbb::blocked_range<size_t>(1, samples.size()), [&](const tbb::blocked_range<size_t>& range) {
            for (size_t i = range.begin(); i < range.end(); i++) {
                if (signals::check_interrupted()) continue;
                const auto& s = samples[i];
                std::string outPath = s.outputPrefix + ".placement.tsv";

                fs::path outDir = fs::path(s.outputPrefix).parent_path();
                if (!outDir.empty()) fs::create_directories(outDir);

                placement::PlacementResult result;
                auto start = std::chrono::high_resolution_clock::now();
                placement::placeLite(result, &tree, reader, s.reads1, s.reads2, outPath, batchParams, fullTreePtr);

                if (result.bestLogContainmentNodeId.empty()) {
                    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::high_resolution_clock::now() - start);
                    int n = completedCount.fetch_add(1, std::memory_order_relaxed) + 1;
                    std::lock_guard<std::mutex> lock(progressMutex);
                    fmt::print(stderr,
                               "[{}/{}] {} -> NO PLACEMENT ({}ms)\n",
                               n,
                               samples.size(),
                               s.outputPrefix,
                               elapsed.count());
                    atomicFail.fetch_add(1, std::memory_order_relaxed);
                    continue;
                }

                bool sampleOk = true;
                if (cfg.stopAfter >= PipelineStage::Align) {
                    Config sampleCfg = cfg;
                    sampleCfg.output = s.outputPrefix;
                    sampleCfg.reads1 = s.reads1;
                    sampleCfg.reads2 = s.reads2;

                    if (runAlignment(sampleCfg, result, fullTreePtr) != 0) {
                        int n = completedCount.fetch_add(1, std::memory_order_relaxed) + 1;
                        std::lock_guard<std::mutex> lock(progressMutex);
                        fmt::print(stderr, "[{}/{}] {} -> alignment failed\n", n, samples.size(), s.outputPrefix);
                        atomicFail.fetch_add(1, std::memory_order_relaxed);
                        continue;
                    }

                    if (cfg.stopAfter >= PipelineStage::Genotype) {
                        if (runGenotyping(sampleCfg) != 0) {
                            int n = completedCount.fetch_add(1, std::memory_order_relaxed) + 1;
                            std::lock_guard<std::mutex> lock(progressMutex);
                            fmt::print(stderr, "[{}/{}] {} -> genotyping failed\n", n, samples.size(), s.outputPrefix);
                            atomicFail.fetch_add(1, std::memory_order_relaxed);
                            continue;
                        }
                    }
                }

                auto totalElapsedSample = std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::high_resolution_clock::now() - start);
                int n = completedCount.fetch_add(1, std::memory_order_relaxed) + 1;
                std::lock_guard<std::mutex> lock(progressMutex);
                fmt::print(stderr,
                           "[{}/{}] {} -> {} ({}ms)\n",
                           n,
                           samples.size(),
                           s.outputPrefix,
                           result.bestLogContainmentNodeId,
                           totalElapsedSample.count());
                atomicSuccess.fetch_add(1, std::memory_order_relaxed);
            }
        });

        output::config().quiet = false;
        successCount += atomicSuccess.load();
        failCount += atomicFail.load();
    }

    auto totalElapsed =
        std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - batchStart);
    output::done(
        fmt::format("Batch complete: {}/{} placed in {}s", successCount, samples.size(), totalElapsed.count()));
    if (failCount > 0) logging::warn("{} samples failed placement", failCount);
    return 0;
}

std::optional<placement::PlacementResult> runPlacement(const Config& cfg) {
    logging::msg("Loading index...");
    IndexReader reader(cfg.index, cfg.threads);

    auto idx = reader.getRoot<LiteIndex>();
    logging::debug("Index parameters: k={}, s={}, l={}", idx.getK(), idx.getS(), idx.getL());

    panmapUtils::LiteTree tree;
    tree.initialize(idx.getLiteTree());
    logging::msg("Tree loaded: {} nodes", tree.allLiteNodes.size());

    placement::PlacementResult result;
    std::string outPath = cfg.output + ".placement.tsv";

    bool storeDiagnostics = false;

    // Load full tree if refinement is enabled (needed for genome sequences)
    std::unique_ptr<panmanUtils::TreeGroup> tg;
    panmanUtils::Tree* fullTreePtr = nullptr;
    bool refineEnabled = cfg.refine;
    if (refineEnabled) {
        tg = loadPanMAN(cfg.panman);
        if (tg && !tg->trees.empty()) {
            fullTreePtr = &tg->trees[0];
            logging::msg("Loaded full tree for refinement ({} nodes)", fullTreePtr->allNodes.size());
        } else {
            logging::warn("Failed to load full tree for refinement - disabling refinement");
            refineEnabled = false;
        }
    }

    auto start = std::chrono::high_resolution_clock::now();
    placement::TraversalParams params;
    params.store_diagnostics = storeDiagnostics;
    params.seedMaskFraction = cfg.seedMaskFraction;
    params.minSeedQuality = cfg.minSeedQuality;
    params.dedupReads = cfg.dedupReads;
    params.trimStart = cfg.trimStart;
    params.trimEnd = cfg.trimEnd;
    params.minReadSupport = cfg.minReadSupport;
    params.refineEnabled = refineEnabled;
    params.refineTopPct = cfg.refineTopPct;
    params.refineMaxTopN = cfg.refineMaxTopN;
    params.refineNeighborRadius = cfg.refineNeighborRadius;
    params.refineMaxNeighborN = cfg.refineMaxNeighborN;
    params.forceLeaf = cfg.forceLeaf;
    placement::placeLite(result, &tree, reader, cfg.reads1, cfg.reads2, outPath, params, fullTreePtr);
    auto elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);

    logging::msg("Placement complete in {}ms", elapsed.count());

    // Check if any metric found a placement
    if (result.bestLogContainmentNodeId.empty()) {
        logging::warn("No placement found");
        return std::nullopt;
    }

    // Report results
    logging::msg("{}Best placement:{} {} (LogContainment: {:.6f}, LogRaw: {:.6f}, LogCosine: {:.6f}, Containment: {:.6f})",
                 color::green(),
                 color::reset(),
                 result.bestLogContainmentNodeId,
                 result.bestLogContainmentScore,
                 result.bestLogRawScore,
                 result.bestLogCosineScore,
                 result.bestContainmentScore);

    // Dump all scores to file if requested
    if (!cfg.dumpAllScores.empty()) {
        std::ofstream outFile(cfg.dumpAllScores);
        if (outFile) {
            outFile << "node\tlogRaw\tlogCosine\tcontainment\n";
            std::vector<std::pair<double, std::string>> allScores;
            for (auto& [id, node] : tree.allLiteNodes) {
                if (node->logRawScore > 0 || node->logCosineScore > 0 || node->containmentScore > 0) {
                    allScores.push_back({node->logRawScore, id});
                }
            }
            // Sort by logRaw descending
            std::sort(allScores.begin(), allScores.end(), std::greater<>());
            for (auto& [score, id] : allScores) {
                auto* node = tree.allLiteNodes[id];
                outFile << id << "\t" << node->logRawScore << "\t" << node->logCosineScore << "\t"
                        << node->containmentScore << "\n";
            }
            logging::msg("Dumped {} node scores to {}", allScores.size(), cfg.dumpAllScores);
        } else {
            logging::warn("Could not open {} for writing", cfg.dumpAllScores);
        }
    }

    if (cfg.metagenomic && cfg.topN > 1) {
        logging::msg("Top {} placements written to {}", cfg.topN, outPath);
    }

    return result;
}

int runAlignment(const Config& cfg, const placement::PlacementResult& placement, panmanUtils::Tree* preloadedTree) {
    std::unique_ptr<panmanUtils::TreeGroup> tg;
    panmanUtils::Tree* T = preloadedTree;
    if (!T) {
        logging::msg("Loading tree for alignment...");
        tg = loadPanMAN(cfg.panman);
        if (!tg || tg->trees.empty()) {
            logging::err("Failed to load tree");
            return 1;
        }
        T = &tg->trees[0];
    }
    // Use LogContainment as primary placement metric
    std::string nodeId = placement.bestLogContainmentNodeId;

    if (nodeId.empty()) {
        logging::err("No best node ID from placement - cannot align");
        return 1;
    }

    std::string bestMatchSequence = panmapUtils::getStringFromReference(T, nodeId, false);

    if (bestMatchSequence.empty()) {
        logging::err("Empty sequence for node '{}' - cannot align", nodeId);
        return 1;
    }

    // Output file paths
    std::string refFileName = cfg.output + ".ref.fa";
    std::string samFileName = cfg.output + ".sam";
    std::string bamFileName = cfg.output + ".bam";

    // Write reference fasta
    {
        std::ofstream outFile(refFileName);
        if (!outFile) {
            logging::err("Cannot write reference file: {}", refFileName);
            return 1;
        }
        outFile << ">ref\n" << bestMatchSequence << "\n";
        outFile.close();
    }

    // Create FASTA index (.fai) for htslib/samtools
    if (fai_build(refFileName.c_str()) != 0) {
        logging::err("Failed to create FASTA index for {}", refFileName);
        return 1;
    }

    logging::msg("Reference: {} bp from node {} -> {}", bestMatchSequence.size(), nodeId, refFileName);

    auto start = std::chrono::high_resolution_clock::now();

    // Opt 1: Lightweight FASTQ reader (no seed computation)
    std::vector<std::string> readSequences, readQuals, readNames;
    seeding::readFastqPaired(readSequences, readQuals, readNames, cfg.reads1, cfg.reads2);

    logging::msg("Loaded {} reads", readSequences.size());

    // Opts 2+3: Parallel alignment with direct BAM construction
    bool pairedEndReads = !cfg.reads2.empty();
    alignAndWriteBam(readSequences, readQuals, readNames, bestMatchSequence, bamFileName, pairedEndReads, cfg.threads);

    auto elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);

    logging::msg("Alignment complete in {}ms -> {}", elapsed.count(), bamFileName);
    return 0;
}

int runGenotyping(const Config& cfg) {
    std::string prefix = cfg.output;
    std::string refFileName = cfg.output + ".ref.fa";
    std::string bamFileName = cfg.output + ".bam";
    std::string mpileupFileName = cfg.output + ".mpileup";
    std::string vcfFileName = cfg.output + ".vcf";

    // Read reference sequence
    std::string bestMatchSequence;
    std::ifstream ref(refFileName);
    std::string line;
    while (std::getline(ref, line)) {
        if (line[0] != '>') bestMatchSequence += line;
    }

    auto start = std::chrono::high_resolution_clock::now();

    // Create mpileup and VCF
    createMplpBcf(prefix, refFileName, bestMatchSequence, bamFileName, mpileupFileName, cfg.baq);
    createVcfWithMutationMatrices(prefix, mpileupFileName, vcfFileName, cfg.substMatrixPhred);

    auto elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);

    logging::msg("Genotyping complete in {}ms -> {}", elapsed.count(), vcfFileName);
    return 0;
}

int runConsensus(const Config& cfg) {
    std::string refFileName = cfg.output + ".ref.fa";
    std::string vcfFileName = cfg.output + ".vcf";
    std::string consensusFileName = cfg.output + ".consensus.fa";

    auto start = std::chrono::high_resolution_clock::now();

    if (createConsensus(vcfFileName, refFileName, consensusFileName) != 0) {
        logging::err("bcftools consensus failed");
        return 1;
    }

    auto elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);

    logging::msg("Consensus complete in {}ms -> {}", elapsed.count(), consensusFileName);
    return 0;
}

void printUsage() {
    std::cout << color::bold() << "panmap" << color::reset() << " v" << VERSION << "\n";
    std::cout << "Pangenome-based sequence placement, alignment, and genotyping\n\n";
    std::cout << color::bold() << "Usage:" << color::reset() << "  panmap [options] <panman> [reads.fq] [reads2.fq]\n";
    std::cout << color::bold() << "Output:" << color::reset() << " <prefix>.vcf, <prefix>.bam, <prefix>.consensus.fa, ...\n"
              << "        (prefix defaults to reads filename, or use -o)\n\n";
}

int main(int argc, char** argv) {
    Config cfg;

    // Define options - organized into visible (common) and advanced groups

    // === Common options (shown in --help) ===
    po::options_description visible("Options");
    visible.add_options()
        ("help,h", "Show help (--help-all for more)")
        ("help-all", "Show all options")
        ("version,V", "Show version")
        ("output,o", po::value<std::string>(&cfg.output), "Output prefix")
        ("threads,t", po::value<int>(&cfg.threads)->default_value(1), "Threads")
        ("stop", po::value<std::string>()->default_value("consensus"), "Stop after: index|place|align|genotype|consensus")
        ("meta", po::bool_switch(&cfg.metagenomic), "Metagenomic mode (for more options, see --help-all)")
        ("aligner,a", po::value<std::string>(&cfg.aligner)->default_value("minimap2"), "Aligner: minimap2|bwa")
        ("verbose,v", po::bool_switch(&cfg.verbose), "Verbose output")
        ("quiet,q", po::bool_switch(&cfg.quiet), "Quiet output")
        ("no-color", po::bool_switch(&cfg.plain), "No colors");
    
    po::options_description advanced("Advanced");
    advanced.add_options()
        ("index,i", po::value<std::string>(&cfg.index), "Index file path")
        ("reindex,f", po::bool_switch(&cfg.forceReindex), "Force rebuild index")
        ("dedup", po::bool_switch(&cfg.dedupReads), "Deduplicate reads")
        ("impute", po::bool_switch(&cfg.impute), "Impute N's from parent (skip _->N mutations in indexing and output)")
        ("no-mutation-spectrum", po::bool_switch(&cfg.noMutationSpectrum), "Disable mutation spectrum filtering in VCF genotyping")
        ("baq", po::bool_switch(&cfg.baq), "Enable BAQ (Base Alignment Quality) in mpileup (default: off)")
        ("kmer,k", po::value<int>(&cfg.k)->default_value(19), "Syncmer k")
        ("syncmer,s", po::value<int>(&cfg.s)->default_value(8), "Syncmer s")
        ("offset", po::value<int>(&cfg.t)->default_value(0), "Syncmer offset")
        ("lmer,l", po::value<int>(&cfg.l)->default_value(3), "Syncmers per seed")
        ("open-syncmer", po::bool_switch(&cfg.openSyncmer), "Open syncmer")
        ("flank-mask", po::value<int>(&cfg.flankMaskBp)->default_value(250), "Mask bp at ends")
        ("seed-mask-fraction", po::value<double>(&cfg.seedMaskFraction)->default_value(0), "Mask top seed fraction")
        ("min-seed-quality", po::value<int>(&cfg.minSeedQuality)->default_value(0), "Min seed quality")
        ("trim-start", po::value<int>(&cfg.trimStart)->default_value(0), "Trim read start")
        ("trim-end", po::value<int>(&cfg.trimEnd)->default_value(0), "Trim read end")
        ("min-read-support", po::value<int>(&cfg.minReadSupport)->default_value(1), "Min reads for a seed (2=filter singletons)")
        ("hpc", po::bool_switch(&cfg.hpc), "Homopolymer-compressed seeds")
        ("extent-guard", po::bool_switch(&cfg.extentGuard), "Guard seed deletions at genome extent boundaries")
        ("force-leaf", po::bool_switch(&cfg.forceLeaf), "Restrict placement to leaf nodes only (default when --stop genotype)")
        ("refine", po::bool_switch(&cfg.refine), "Enable alignment-based refinement")
        ("refine-top-pct", po::value<double>(&cfg.refineTopPct)->default_value(0.01), "Top % of nodes to refine (default 1%)")
        ("refine-max-top-n", po::value<int>(&cfg.refineMaxTopN)->default_value(150), "Max nodes to align against")
        ("refine-neighbor-radius", po::value<int>(&cfg.refineNeighborRadius)->default_value(2), "Expand to neighbors within N branches")
        ("refine-max-neighbor-n", po::value<int>(&cfg.refineMaxNeighborN)->default_value(150), "Max additional nodes from neighbor expansion")
        ("zstd-level", po::value<int>(&cfg.zstdLevel)->default_value(7), "ZSTD compression level for index (1-22)")
        ("batch", po::value<std::string>(&cfg.batchFile), "Batch file listing samples (one per line: reads1 [reads2] [output_prefix])");
    
    po::options_description metagenomic("Metagenomic");
    metagenomic.add_options()
        ("index-mgsr", po::value<std::string>(&cfg.indexMgsr), "Path to build/rebuild MGSR index")
        ("index-full", po::bool_switch(&cfg.indexFull), "Build full index (default index-mgsr builds lite index)")
        ("index-packed", po::bool_switch(&cfg.indexPacked), "Build packed capnp message (default false)")
        ("read-packed", po::bool_switch(&cfg.readPacked), "Read packed capnp message (default false)")
        ("no-progress", po::bool_switch(&cfg.noProgress), "Disable progress bars");
    
    po::options_description em("Metagenomic: EM");
    em.add_options()
        ("top-oc", po::value<size_t>(&cfg.topOc)->default_value(1000), "Select top <int> nodes by overlap coefficients to send to EM")
        ("mask-reads", po::value<uint32_t>(&cfg.maskReads)->default_value(0), "mask reads containing k-min-mers with total occurrence <= threshold")
        ("mask-seeds", po::value<uint32_t>(&cfg.maskSeeds)->default_value(0), "mask k-min-mer seeds in query with total occurrence <= threshold")
        ("amplicon-depth", po::value<std::string>(&cfg.ampliconDepth), "Path to amplicon depth TSV file (if specified, will be used to mask-reads/seeds basedd)")
        ("mask-reads-relative-frequency", po::value<double>(&cfg.maskReadsRelativeFrequency)->default_value(0.0), "mask reads containing k-min-mers with relative frequency < threadshold * amplicon_depth")
        ("mask-seeds-relative-frequency", po::value<double>(&cfg.maskSeedsRelativeFrequency)->default_value(0.0), "mask k-min-mer seeds in query with with relative frequency < threadshold * amplicon_depth")
        ("em-convergence-threshold", po::value<double>(&cfg.emConvergenceThreshold)->default_value(0.00001), "EM converges when likelihood difference is less than <float> (choose em-convergence-threshold or em-delta-threshold, default is em-convergence-threshold)")
        ("em-delta-threshold", po::value<double>(&cfg.emDeltaThreshold)->default_value(0.0), "EM converges when maximum proportion change is less than <float> (choose em-convergence-threshold or em-delta-threshold, default is em-delta-threshold)")
        ("em-maximum-rounds", po::value<uint32_t>(&cfg.emMaximumRounds)->default_value(5), "EM maximum rounds")
        ("em-maximum-iterations", po::value<uint32_t>(&cfg.emMaximumIterations)->default_value(1000), "EM maximum iterations")
        ("em-leaves-only", po::bool_switch(&cfg.emLeavesOnly), "Only run EM on leaf (sample) nodes");
    
    po::options_description filterAndAssign("Metagenomic: Filter and Assign");
    filterAndAssign.add_options()
        ("filter-and-assign", po::bool_switch(&cfg.filterAndAssign), "Filter and assign reads to nodes without running EM")
        ("dust", po::value<double>(&cfg.dust)->default_value(100.0), "Discard reads with Prinseq scale dust score > <FLOAT> (default 100, i.e. no dust filtering)")
        ("discard", po::value<double>(&cfg.discard)->default_value(0.0), "Discard reads with maximum parsimony score < FLOAT * read_total_seed (default 0, i.e. no discard)")
        ("mask-read-ends", po::value<uint32_t>(&cfg.maskReadEnds)->default_value(0), "mask <int> bases from the beginning and end of reads (for ancient eDNA damage)")
        ("taxonomic-metadata", po::value<std::string>(&cfg.taxonomicMetadata), "Path to taxonomic metadata TSV file")
        ("taxonomic-rank", po::value<std::string>(&cfg.taxonomicRank)->default_value("Family"), "Taxonomic rank to use for filtering and assigning reads (should match the column name in the taxonomic metadata TSV file), only applicable if taxonomic-metadata is provided")
        ("maximum-taxon-number", po::value<size_t>(&cfg.maximumTaxonNumber)->default_value(1), "Discard reads assigned to nodes spanning more than <int> distinct taxonos at the specified taxonomic rank, only applicable if taxonomic-metadata is provided")
        ("ambiguous-score-threshold-ratio", po::value<double>(&cfg.ambiguousScoreThresholdRatio)->default_value(0.0), "Discard reads scoring max score - <double> * max score outside of the max scoring families")
        ("ambiguous-score-threshold", po::value<int>(&cfg.ambiguousScoreThreshold)->default_value(0), "Discard reads scoring max score - <int> outside of the max scoring families")
        ("breadth-ratio", po::bool_switch(&cfg.breadthRatio), "Calculate observed / expected breadth ratio")
        ("pseudochain", po::bool_switch(&cfg.pseudochain), "Use pseudo-chains for scoring reads (default: off)")
        ("batch-files-path", po::value<std::string>(&cfg.batchFilesPath), "Path to tsv file containg batch file paths")
        ("batch-size", po::value<size_t>(&cfg.batchSize)->default_value(1000000), "Batch size for filtering and assigning reads");
    
    po::options_description developer("Developer");
    developer.add_options()
        ("dump-sequence", po::value<std::string>(&cfg.dumpNodeId), "Dump node FASTA")
        ("dump-all-scores", po::value<std::string>(&cfg.dumpAllScores), "Dump all node scores to TSV file")
        ("write-meta-read-scores-filtered", po::bool_switch(&cfg.writeMetaReadScoresFiltered), "Write filtered meta read scores to TSV file")
        ("write-meta-read-scores-unfiltered", po::bool_switch(&cfg.writeMetaReadScoresUnfiltered), "Write unfiltered meta read scores to TSV file")
        ("write-ocranks", po::bool_switch(&cfg.writeOCRanks), "Write overlap coefficients info to TSV file")
        ("seed", po::value<int>(&cfg.seed)->default_value(42), "Random seed");

    // Positional arguments (always hidden)
    po::options_description positional;
    positional.add_options()("panman", po::value<std::string>(&cfg.panman), "")(
        "reads1", po::value<std::string>(&cfg.reads1), "")("reads2", po::value<std::string>(&cfg.reads2), "");

    po::positional_options_description pos;
    pos.add("panman", 1).add("reads1", 1).add("reads2", 1);

    // Combine option groups
    po::options_description all;  // For parsing
    all.add(visible).add(advanced).add(metagenomic).add(em).add(filterAndAssign).add(developer).add(positional);

    po::options_description visible_all;  // For --help-all
    visible_all.add(visible).add(advanced).add(metagenomic).add(em).add(filterAndAssign).add(developer);

    // Parse
    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).options(all).positional(pos).run(), vm);
        po::notify(vm);
    } catch (const po::error& e) {
        // Check for NO_COLOR env var
        const char* noColorEnv = std::getenv("NO_COLOR");
        bool plainMode = vm.count("no-color") > 0 || (noColorEnv && noColorEnv[0] != '\0');
        output::init(false, false, plainMode);
        output::error("{}", e.what());
        std::cerr << "\n";
        printUsage();
        std::cout << visible << "\n";
        return 1;
    }

    // Check for NO_COLOR environment variable (no-color.org standard)
    const char* noColorEnv = std::getenv("NO_COLOR");
    if (noColorEnv && noColorEnv[0] != '\0') cfg.plain = true;

    // Initialize output early for help/version formatting
    output::init(cfg.quiet, cfg.verbose, cfg.plain);

    // Handle help/version
    if (vm.count("help") || argc == 1) {
        printUsage();
        std::cout << visible << "\n";
        return 0;
    }

    if (vm.count("help-all")) {
        printUsage();
        std::cout << visible_all << "\n";
        return 0;
    }

    if (vm.count("version")) {
        std::cout << PROGRAM_NAME << " " << VERSION << "\n";
        return 0;
    }

    // Validate required args
    if (cfg.panman.empty()) {
        output::error("PanMAN file required");
        return 1;
    }

    // Auto-detect: if panman arg looks like an index file (.pmi/.idx), use it as the index
    {
        auto ext = fs::path(cfg.panman).extension().string();
        if (ext == ".pmi" || ext == ".idx") {
            if (cfg.index.empty()) {
                cfg.index = cfg.panman;
            }
            output::error("'{}' looks like an index file, not a PanMAN.", cfg.panman);
            output::error("Usage: panmap <panman> [reads.fq] -i {}", cfg.panman);
            return 1;
        }
    }

    // Parse stop stage
    std::string stopStr = vm["stop"].as<std::string>();
    if (stopStr == "index")
        cfg.stopAfter = PipelineStage::Index;
    else if (stopStr == "place")
        cfg.stopAfter = PipelineStage::Place;
    else if (stopStr == "align")
        cfg.stopAfter = PipelineStage::Align;
    else if (stopStr == "genotype")
        cfg.stopAfter = PipelineStage::Genotype;
    else if (stopStr == "consensus")
        cfg.stopAfter = PipelineStage::Consensus;
    else {
        output::error("Invalid stage '{}'", stopStr);
        return 1;
    }

    // Default --force-leaf on when running genotype (unless user explicitly set it)
    if (!vm.count("force-leaf") && cfg.stopAfter >= PipelineStage::Genotype) {
        cfg.forceLeaf = true;
    }

    // Set defaults
    // If index not explicitly set, derive from output prefix if set, otherwise from panman
    if (cfg.index.empty()) {
        if (!cfg.output.empty()) {
            cfg.index = cfg.output + ".idx";
        } else {
            cfg.index = cfg.panman + ".idx";
        }
    }
    // Only set default output if we're not in dump mode (dump modes handle their own output)
    if (cfg.output.empty() && cfg.dumpNodeId.empty()) {
        // Derive output prefix from reads filename (without path and common extensions)
        if (!cfg.reads1.empty()) {
            fs::path readsPath(cfg.reads1);
            std::string stem = readsPath.stem().string();
            // Remove common paired-end suffixes like _R1, _1, .R1, etc.
            for (const auto& suffix : {"_R1", "_R2", "_1", "_2", ".R1", ".R2", ".1", ".2"}) {
                if (stem.size() > strlen(suffix) && stem.substr(stem.size() - strlen(suffix)) == suffix) {
                    stem = stem.substr(0, stem.size() - strlen(suffix));
                    break;
                }
            }
            // Also remove .fastq, .fq extensions if present (for .fastq.gz -> .fastq stem)
            for (const auto& ext : {".fastq", ".fq"}) {
                if (stem.size() > strlen(ext) && stem.substr(stem.size() - strlen(ext)) == ext) {
                    stem = stem.substr(0, stem.size() - strlen(ext));
                    break;
                }
            }
            cfg.output = stem;
        } else {
            cfg.output = cfg.panman;
        }
    }

    // Install signal handlers for graceful interruption
    signals::install_handlers();

    // Initialize threading
    tbb::global_control tbb_ctl(tbb::global_control::max_allowed_parallelism, cfg.threads);

    // ========================================================================
    // Print Configuration Summary
    // ========================================================================

    auto printConfigSummary = [&]() {
        if (cfg.quiet) return;                // Skip in quiet mode
        if (!cfg.dumpNodeId.empty()) return;  // Skip for utility modes

        // Build stage string using arrow from box chars
        std::string arrow = output::box::arrow();
        std::string stageStr;
        switch (cfg.stopAfter) {
            case PipelineStage::Index: stageStr = "index"; break;
            case PipelineStage::Place: stageStr = fmt::format("index {} place", arrow); break;
            case PipelineStage::Align: stageStr = fmt::format("index {} place {} align", arrow, arrow); break;
            case PipelineStage::Genotype:
                stageStr = fmt::format("index {} place {} align {} genotype", arrow, arrow, arrow);
                break;
            default: stageStr = "full"; break;
        }

        // Build input string
        std::string inputStr = cfg.panman;
        if (!cfg.reads1.empty()) {
            inputStr += "  + " + cfg.reads1;
            if (!cfg.reads2.empty()) inputStr += ", " + cfg.reads2;
        }

        // Build config string
        std::string configStr = fmt::format("threads={}  k={} s={} l={}", cfg.threads, cfg.k, cfg.s, cfg.l);
        if (cfg.hpc) {
            configStr += "  hpc";
        }
        if (cfg.metagenomic) {
            configStr += fmt::format("  meta(top={})", cfg.topN);
        }
        if (cfg.forceReindex) {
            configStr += "  reindex";
        }
        if (cfg.aligner != "minimap2") {
            configStr += "  aligner=" + cfg.aligner;
        }

        output::print_header("panmap", VERSION);
        output::print_row("Input ", inputStr);
        output::print_row("Output", cfg.output + ".*");
        output::print_row("Stages", stageStr);
        output::print_row("Config", configStr);
        output::print_footer();
    };

    printConfigSummary();

    // ========================================================================
    // Run Pipeline
    // ========================================================================

    try {
        // Utility: dump specific node sequence
        if (!cfg.dumpNodeId.empty()) {
            auto tg = loadPanMAN(cfg.panman);
            // Use -o output path if specified, otherwise default to panman path
            std::string outPath =
                cfg.output.empty() ? cfg.panman + "." + sanitizeFilename(cfg.dumpNodeId) + ".fa" : cfg.output;
            saveNodeSequence(&tg->trees[0], cfg.dumpNodeId, outPath);
            std::cout << cfg.dumpNodeId << "\n";
            return 0;
        }

        // metagenomics mode related
        if (!cfg.indexMgsr.empty()) {
            if (!buildMgsrIndex(cfg)) return 1;
            return 0;
        }

        if (cfg.metagenomic) {
            if (!runMetagenomic(cfg)) return 1;
            std::cout << "Metagenomic mode run completed" << std::endl;
            return 0;
        }

        // Stage 1: Index
        if (!buildIndex(cfg)) return 1;
        if (signals::check_interrupted()) return 130;  // Standard exit code for SIGINT

        // Load substitution spectrum from index for genotyping
        if (!cfg.noMutationSpectrum && cfg.stopAfter >= PipelineStage::Genotype) {
            IndexReader specReader(cfg.index, cfg.threads);
            auto specIdx = specReader.getRoot<LiteIndex>();
            cfg.substMatrixPhred = loadSubstMatrixFromIndex(specIdx);
            if (!cfg.substMatrixPhred.empty()) {
                logging::msg("Loaded substitution spectrum from index");
            }
        }

        if (cfg.stopAfter == PipelineStage::Index) {
            output::done("Index ready: " + cfg.index);
            output::info("");
            output::info("{}Next steps:{}", color::dim(), color::reset());
            output::info("  {} panmap {} reads.fq", color::dim(), cfg.panman);
            return 0;
        }

        // Batch mode: load index once, place all samples
        if (!cfg.batchFile.empty()) {
            return runBatchPlacement(cfg);
        }

        // Check for reads
        if (cfg.reads1.empty()) {
            output::info("No reads provided. Index is ready.");
            output::info("");
            output::info("{}Next steps:{}", color::dim(), color::reset());
            output::info("  {} panmap {} reads.fq", color::dim(), cfg.panman);
            return 0;
        }

        // Stage 2: Placement
        auto placement = runPlacement(cfg);
        if (signals::check_interrupted()) return 130;
        if (!placement) {
            output::error("Placement failed");
            return 1;
        }
        if (cfg.stopAfter == PipelineStage::Place) {
            output::done("Placement: " + cfg.output + ".placement.tsv");
            output::info("");
            output::info("{}Next steps:{}", color::dim(), color::reset());
            output::info("  {} Continue to alignment: panmap {} {} --stop align", color::dim(), cfg.panman, cfg.reads1);
            return 0;
        }

        // Stage 3: Alignment
        if (runAlignment(cfg, *placement) != 0) {
            output::error("Alignment failed");
            return 1;
        }
        if (signals::check_interrupted()) return 130;
        if (cfg.stopAfter == PipelineStage::Align) {
            output::done("Alignment: " + cfg.output + ".bam");
            output::info("");
            output::info("{}Next steps:{}", color::dim(), color::reset());
            output::info("  {} View: samtools view {}.bam | head", color::dim(), cfg.output);
            output::info("  {} Stats: samtools flagstat {}.bam", color::dim(), cfg.output);
            return 0;
        }

        // Stage 4: Genotyping
        if (runGenotyping(cfg) != 0) {
            output::error("Genotyping failed");
            return 1;
        }
        if (signals::check_interrupted()) return 130;
        if (cfg.stopAfter == PipelineStage::Genotype) {
            output::done("Variants: " + cfg.output + ".vcf");
            output::info("");
            output::info("{}Next steps:{}", color::dim(), color::reset());
            output::info("  {} View variants: bcftools view {}.vcf | head", color::dim(), cfg.output);
            output::info("  {} Consensus: panmap {} {} --stop consensus", color::dim(), cfg.panman, cfg.reads1);
            return 0;
        }

        // Stage 5: Consensus
        if (runConsensus(cfg) != 0) {
            output::error("Consensus generation failed");
            return 1;
        }
        if (signals::check_interrupted()) return 130;

        output::info("");
        output::info("{}Pipeline complete.{}", color::green(), color::reset());
        output::info("  Placement:  {}.placement.tsv", cfg.output);
        output::info("  Reference:  {}.ref.fa", cfg.output);
        output::info("  Alignment:  {}.bam", cfg.output);
        output::info("  Variants:   {}.vcf", cfg.output);
        output::info("  Consensus:  {}.consensus.fa", cfg.output);
        output::info("");
        output::info("{}Next steps:{}", color::dim(), color::reset());
        output::info("  {} View variants: bcftools view {}.vcf | head", color::dim(), cfg.output);
        output::info("  {} View alignment: samtools tview {}.bam {}.ref.fa", color::dim(), cfg.output, cfg.output);

        return 0;

    } catch (const std::exception& e) {
        output::error("Fatal error: {}", e.what());
        return 1;
    }
}
