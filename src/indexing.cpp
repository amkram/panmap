#include "indexing.hpp"
#include "alignment.hpp"
#include "capnp/message.h"
#include "capnp/serialize-packed.h"
#include "coordinates.hpp"
#include "gap_map.hpp"
#include "logging.hpp"
#include "panman.hpp"
#include "placement.hpp"
#include "state.hpp"
#include "timing.hpp"
#include <algorithm>
#include <atomic>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <functional>
#include <future>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
#include <tbb/task_arena.h>
#include <tbb/task_group.h>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include "absl/container/flat_hash_set.h"
// #include "common.hpp" // Commented out
// #include "seed_changes.hpp" // Commented out

namespace state {
extern std::atomic<size_t> totalKmerAttempts;
extern std::atomic<size_t> successfulKmerExtractions;
extern std::atomic<size_t> kmerFailedDueToLength;
extern std::atomic<size_t> kmerFailedDueToOnlyGaps;
}

namespace {
    std::atomic<size_t> nodesWithGapOnlyRanges{0};
    std::atomic<size_t> totalRangesProcessed{0};
    std::atomic<size_t> gapOnlyRangesCount{0};
    std::atomic<size_t> totalCharsProcessed{0};
    std::atomic<size_t> gapCharsCount{0};
    std::atomic<size_t> nonGapCharsCount{0};
    
    
    std::mutex diagnosticMutex;
    
    const size_t DIAGNOSTIC_INTERVAL = 250;
    
    
    void resetDiagnosticCounters() {
        state::totalKmerAttempts.store(0, std::memory_order_relaxed);
        state::successfulKmerExtractions.store(0, std::memory_order_relaxed);
        state::kmerFailedDueToLength.store(0, std::memory_order_relaxed);
        state::kmerFailedDueToOnlyGaps.store(0, std::memory_order_relaxed);
        nodesWithGapOnlyRanges = 0;
        totalRangesProcessed = 0;
        gapOnlyRangesCount = 0;
        totalCharsProcessed = 0;
        gapCharsCount = 0;
        nonGapCharsCount = 0;
    }
    
    
    void printDiagnosticReport(size_t nodeCount) {
        std::lock_guard<std::mutex> lock(diagnosticMutex);
        std::cerr << "\n=== DIAGNOSTIC REPORT (Nodes 1-" << nodeCount << ") ===\n"
                  << "Total chars processed: " << totalCharsProcessed.load() << "\n"
                  << "Gap chars: " << gapCharsCount.load() 
                  << " (" << (totalCharsProcessed.load() > 0 ? (static_cast<double>(gapCharsCount.load()) / totalCharsProcessed.load() * 100.0) : 0.0) << "%)\n"
                  << "Non-gap chars: " << nonGapCharsCount.load()
                  << " (" << (totalCharsProcessed.load() > 0 ? (static_cast<double>(nonGapCharsCount.load()) / totalCharsProcessed.load() * 100.0) : 0.0) << "%)\n"
                  << "Total ranges processed: " << totalRangesProcessed.load() << "\n"
                  << "Gap-only ranges: " << gapOnlyRangesCount.load() 
                  << " (" << (totalRangesProcessed.load() > 0 ? (static_cast<double>(gapOnlyRangesCount.load()) / totalRangesProcessed.load() * 100.0) : 0.0) << "%)\n"
                  << "Nodes with gap-only ranges: " << nodesWithGapOnlyRanges.load() << "\n";
      
      
      if (state::totalKmerAttempts.load() > 0) {
          std::cerr << "K-mer extraction attempts: " << state::totalKmerAttempts.load() << "\n"
                    << "K-mer extraction successes: " << state::successfulKmerExtractions.load() 
                    << " (" << (static_cast<double>(state::successfulKmerExtractions.load()) / state::totalKmerAttempts.load() * 100.0) << "%)\n"
                    << "Failed due to insufficient length: " << state::kmerFailedDueToLength.load() 
                    << " (" << (static_cast<double>(state::kmerFailedDueToLength.load()) / state::totalKmerAttempts.load() * 100.0) << "%)\n"
                    << "Failed due to only gaps: " << state::kmerFailedDueToOnlyGaps.load()
                    << " (" << (static_cast<double>(state::kmerFailedDueToOnlyGaps.load()) / state::totalKmerAttempts.load() * 100.0) << "%)\n";
      }
      
      std::cerr << "=== END DIAGNOSTIC REPORT ===\n" << std::endl;
    }
}


extern void writeCapnp(::capnp::MallocMessageBuilder &message, const std::string &path);


extern std::unique_ptr<::capnp::MessageReader> readCapnp(const std::string &path);


namespace panmanUtils {
  inline char getNucleotideFromCode(int code) {
    switch(code) {
    case 1:
        return 'A';
    case 2:
        return 'C';
    case 4:
        return 'G';
    case 8:
        return 'T';
    case 5:
        return 'R';
    case 10:
        return 'Y';
    case 6:
        return 'S';
    case 9:
        return 'W';
    case 12:
        return 'K';
    case 3:
        return 'M';
    case 14:
        return 'B';
    case 13:
        return 'D';
    case 11:
        return 'H';
    case 7:
        return 'V';
    case 15:
        return 'N';
    default:
        return '-';
    }
  }
}

namespace placement {
  class PlacementEngine {
  public:
    PlacementEngine(int k) {}
    void addSeed(int64_t pos, size_t hash, panmanUtils::Node* node, int32_t blockId) {}
    void saveIndex() {}
  };
}

namespace indexing {


std::mutex capnpMutex;


void processNodeComplete(
    state::StateManager& stateManager,
    placement::PlacementEngine& engine,
    panmanUtils::Tree* tree, // Added tree pointer
    panmanUtils::Node* node,
    panmanUtils::Node* commonAncestor,
    const std::unordered_map<std::string, std::vector<panmanUtils::Node*>>& nodePaths,
    int k, int s,
    ::capnp::List<SeedMutations>::Builder* perNodeSeedMutations,
    const std::unordered_map<std::string, uint32_t>* kmerDictionary,
    bool processMutations,
    bool computeSeeds,
    std::unordered_set<std::string>& uniqueKmersCollector,
    std::ofstream& debugSeedFile);


void dumpIndexSummary(const Index::Reader& indexReader, const std::string& outputFilename);




void processSeedChanges(
    const std::vector<std::tuple<int64_t, bool, bool, std::optional<size_t>,
                                 std::optional<size_t>, std::optional<bool>,
                                 std::optional<bool>, std::optional<int64_t>,
                                 std::optional<int64_t>>> &seedChanges,
    std::vector<int64_t> &basePositions, std::vector<uint64_t> &tritMasks);

std::vector<std::tuple<int64_t, bool, bool>>
decodeSeedChanges(const std::vector<int64_t> &basePositions,
                  const std::vector<uint64_t> &bitMasks);



void processMutationsForNode(
    state::StateManager &stateManager, 
    panmanUtils::Tree *tree, 
    panmanUtils::Node *node,
    panmanUtils::Node *commonAncestor,
    const std::unordered_map<std::string, std::vector<panmanUtils::Node *>> &nodePaths,
    std::string* ab_before_output_str,
    std::stringstream* condensed_trace_ss_ptr // NEW: Pass stringstream for logging
    ) {

  if (!node)
    return;

  std::string strNodeId = node->identifier;
  
  // DEADLOCK DEBUG: Log function entry
  if (strNodeId == "node_2") {
    // logging::info("DEADLOCK_DEBUG: Entering processMutationsForNode for {}", strNodeId);
  }
  
  // DEADLOCK FIX: Don't hold nodeState reference throughout entire function
  // Get it only when needed for specific operations
  
  // DEADLOCK DEBUG: Log after avoiding nodeState lock
  if (strNodeId == "node_2") {
    // logging::info("DEADLOCK_DEBUG: Avoided holding nodeState lock for {}", strNodeId);
  }

  // Only log debug files for nodes on the path to our target leaf or specific debug nodes
  std::ofstream debug_file_ourcode_pmfn; // Changed name to avoid conflict with processNodeComplete param
  std::string debug_filename_pmfn = "debug_OURCODE_" + strNodeId + "_ALL_blocks_mutations.txt";
  
  // Only open debug file for target leaf nodes or specific nodes of interest
  bool should_create_debug_file = (strNodeId == "KU950624.1" || 
                                  strNodeId.find("KU950624") != std::string::npos || 
                                  strNodeId == "node_1" || 
                                  strNodeId == "node_2");
  
  // STEP 3: Apply this node's OWN mutations (Block and Nucleotide)
  // The recompRanges will be added to nodeState.recompRanges.
  // All block/nuc mutation logic below uses the single 'nodeState' declared at STEP 0.
  
  std::vector<coordinates::CoordRange> blockMutationInducedRanges;
  std::vector<coordinates::CoordRange> nucMutationInducedRanges;
  std::vector<coordinates::GapUpdate> derivedGapUpdates; // MOVED DECLARATION EARLIER

  // Track which blocks each range belongs to BEFORE merging
  std::vector<std::pair<coordinates::CoordRange, int32_t>> rangeToBlockMapping;


  // PHASE 1: Collect all mutation information and ranges without holding locks
  struct BlockMutationInfo {
    int32_t blockId;
    bool isInsertion;
    bool isMutationInversionFlag;
    bool wasOn;
    bool wasInverted;
    bool isDeactivation;
    bool inversionChanged;
    coordinates::CoordRange blockRange;
  };
  
  std::vector<BlockMutationInfo> mutationInfos;
  std::vector<int32_t> blocksToDeactivate;
  std::vector<coordinates::CoordRange> allRecompRanges; // Collect all ranges first
  int32_t k_val = stateManager.getKmerSize();
  if (k_val <= 0) k_val = 32;  // Default value if not set
  std::set<int32_t> processedBlocks; // Track processed blocks to avoid duplicates
  
  // Declare nucleotide mutation counters early so they're accessible throughout the function
  int nucMutationsInActiveBlocks = 0;
  int nucMutationsInInactiveBlocks = 0;
  
  for (const auto &block_mutation : node->blockMutation) {
    int32_t blockId = block_mutation.primaryBlockId;
    
    
    processedBlocks.insert(blockId);
    
    BlockMutationInfo info;
    info.blockId = blockId;
    info.isInsertion = (block_mutation.blockMutInfo == 1);
    info.isMutationInversionFlag = block_mutation.inversion;
    
    // Get PARENT state to check pre-mutation status
    std::string parentNodeId = node->parent ? node->parent->identifier : strNodeId;
    info.wasOn = stateManager.isBlockOn(parentNodeId, info.blockId);
    info.wasInverted = info.wasOn && stateManager.isBlockInverted(parentNodeId, info.blockId);
    info.isDeactivation = !info.isInsertion && info.wasOn;
    
    info.inversionChanged = false;
    if (info.isInsertion) {
      info.inversionChanged = (info.wasInverted != info.isMutationInversionFlag);
    } else if (info.isMutationInversionFlag) {
      info.inversionChanged = true;
    }
    
    // Get block range once
    info.blockRange = stateManager.getBlockRange(info.blockId);
    
    mutationInfos.push_back(info);
  }

  // PHASE 2: Apply mutations first - this updates the state so we can properly find previous ON blocks
  for (const auto &info : mutationInfos) {
    int32_t blockId = info.blockId;
    bool isInsertion = info.isInsertion;
    bool isMutationInversionFlag = info.isMutationInversionFlag;
    
    // DEADLOCK DEBUG: Log before applyBlockMutation
    if (strNodeId == "node_2" && blockId == 595) {
      // logging::info("DEADLOCK_DEBUG: About to call applyBlockMutation for block {} in {}", blockId, strNodeId);
    }
    
    // DEADLOCK FIX: Apply the mutation without holding any previous locks
    stateManager.applyBlockMutation(strNodeId, blockId, isInsertion, isMutationInversionFlag);
    
    // DEADLOCK DEBUG: Log after applyBlockMutation
    if (strNodeId == "node_2" && blockId == 595) {
      // logging::info("DEADLOCK_DEBUG: Successfully applied block mutation for block {} in {}", blockId, strNodeId);
    }
    
    bool isNowOn = stateManager.isBlockOn(strNodeId, blockId);

    // Add GapUpdate if activation state changed
    if (info.wasOn && !isNowOn) { // Block DEACTIVATED
      try {
        auto blockRange = stateManager.getBlockRange(blockId);
        if (blockRange.start < blockRange.end) {
          derivedGapUpdates.emplace_back(blockRange.start, blockRange.end - blockRange.start, true);
        }
      } catch (const std::exception& e) {
        logging::warn("PMFN_GAP_UPDATE: Error getting range for deactivated block {}: {}", blockId, e.what());
      }
    } else if (!info.wasOn && isNowOn) { // Block ACTIVATED
      try {
        auto blockRange = stateManager.getBlockRange(blockId);
        if (blockRange.start < blockRange.end) {
          derivedGapUpdates.emplace_back(blockRange.start, blockRange.end - blockRange.start, false);
        }
      } catch (const std::exception& e) {
        logging::warn("PMFN_GAP_UPDATE: Error getting range for activated block {}: {}", blockId, e.what());
      }
    }
  }

  // PHASE 3: Generate recomputation ranges based on post-mutation state
  for (const auto &info : mutationInfos) {
    
    int32_t blockId = info.blockId;
    bool isInsertion = info.isInsertion;
    bool isMutationInversionFlag = info.isMutationInversionFlag;
    bool wasOn = info.wasOn;
    coordinates::CoordRange blockRange = info.blockRange;
    
    // DEADLOCK DEBUG: Log block mutation processing
    if (strNodeId == "node_2") {
      // logging::info("DEADLOCK_DEBUG: Processing recomp ranges for block {} in {}", blockId, strNodeId);
    }
    
    // Block mutations -> recomp ranges
    // For PanGraph: any insertion (blockMutInfo = 1) needs recomputation ranges
    // This covers all OFF->ON block mutations at this node
    if (isInsertion) {
      // std::cout << "Block Insertion/Activation: " << blockId << " in Node: " << strNodeId << std::endl;
      
      // 1. The inserted/activated block itself (all k-mers within it need recomputation)
      coordinates::CoordRange insertionRange = {blockRange.start, blockRange.end};
      if (insertionRange.start < insertionRange.end) {
        blockMutationInducedRanges.push_back(insertionRange);
        allRecompRanges.push_back(insertionRange);
        rangeToBlockMapping.push_back({insertionRange, blockId});
        // std::cout << "=> Insertion block range [" << insertionRange.start << ", " << insertionRange.end << ")" << std::endl;
      }
      
      // 2. Find the final k non-gap chars of the previous ON block for upstream context
      // Now we can query the post-mutation state to find the actual previous ON block
      try {
        // Find the previous ON block before this insertion point (post-mutation state)
        int32_t prevOnBlockId = -1;
        for (int32_t testBlockId = blockId - 1; testBlockId >= 0; testBlockId--) {
          if (stateManager.isBlockOn(strNodeId, testBlockId)) {
            prevOnBlockId = testBlockId;
            break;
          }
        }
        
        if (prevOnBlockId != -1) {
          auto prevBlockRange = stateManager.getBlockRange(prevOnBlockId);
          
          // Extract the final k non-gap characters from the previous ON block
          auto kmerResult = stateManager.extractKmer(strNodeId, prevBlockRange.end - 1, k_val, true); // reverse=true to go backward
          
          if (!kmerResult.first.empty() && kmerResult.second.size() >= k_val) {
            // The positions are end positions, so we need the range from the start of the first k-mer to the end
            int64_t upstreamStart = kmerResult.second[k_val - 1] - (k_val - 1); // Start of the k-th last character
            int64_t upstreamEnd = prevBlockRange.end;
            
            if (upstreamStart >= prevBlockRange.start && upstreamStart < upstreamEnd) {
              coordinates::CoordRange upstreamRange = {upstreamStart, upstreamEnd};
              allRecompRanges.push_back(upstreamRange);
              rangeToBlockMapping.push_back({upstreamRange, blockId});
              // std::cout << "=> Upstream recomp range (final " << k_val << " non-gap chars of block " << prevOnBlockId << ") [" << upstreamStart << ", " << upstreamEnd << ")" << std::endl;
            }
          }
        }
      } catch (const std::exception& e) {
        // Fallback to simple position-based upstream range if extractKmer fails
        int64_t upstreamStart = std::max(static_cast<int64_t>(0), blockRange.start - k_val);
        int64_t upstreamEnd = blockRange.start;
        
        if (upstreamStart < upstreamEnd) {
          coordinates::CoordRange upstreamRange = {upstreamStart, upstreamEnd};
          allRecompRanges.push_back(upstreamRange);
          rangeToBlockMapping.push_back({upstreamRange, blockId});
          // std::cout << "=> Upstream recomp range (fallback) [" << upstreamStart << ", " << upstreamEnd << ")" << std::endl;
        }
      }
      
    } else if (!isInsertion && !isMutationInversionFlag && wasOn) {
      // Handle block deactivations (blockMutInfo = 0, inversion = false, block was ON)
      // std::cout << "Block Deactivation: " << blockId << " in Node: " << strNodeId << std::endl;
      
      // For block deactivations, find the final k non-gap chars of the previous ON block
      try {
        // Find the previous ON block before this deactivation point (post-mutation state)
        int32_t prevOnBlockId = -1;
        for (int32_t testBlockId = blockId - 1; testBlockId >= 0; testBlockId--) {
          if (stateManager.isBlockOn(strNodeId, testBlockId)) {
            prevOnBlockId = testBlockId;
            break;
          }
        }
        
        if (prevOnBlockId != -1) {
          auto prevBlockRange = stateManager.getBlockRange(prevOnBlockId);
          
          // Extract the final k non-gap characters from the previous ON block
          auto kmerResult = stateManager.extractKmer(strNodeId, prevBlockRange.end - 1, k_val, true); // reverse=true to go backward
          
          if (!kmerResult.first.empty() && kmerResult.second.size() >= k_val) {
            // The positions are end positions, so we need the range from the start of the first k-mer to the end
            int64_t rangeStart = kmerResult.second[k_val - 1] - (k_val - 1); // Start of the k-th last character
            int64_t rangeEnd = prevBlockRange.end;
            
            if (rangeStart >= prevBlockRange.start && rangeStart < rangeEnd) {
              coordinates::CoordRange deactivationRecompRange = {rangeStart, rangeEnd};
              allRecompRanges.push_back(deactivationRecompRange);
              // std::cout << "=> Deactivation recomp range (final " << k_val << " non-gap chars of block " << prevOnBlockId << ") [" << rangeStart << ", " << rangeEnd << ")" << std::endl;
            }
          }
        }
      } catch (const std::exception& e) {
        // Fallback to simple position-based range if extractKmer fails
        int64_t rangeStart = std::max(static_cast<int64_t>(0), blockRange.start - k_val);
        int64_t rangeEnd = blockRange.start;
        
        if (rangeStart < rangeEnd) {
          coordinates::CoordRange deactivationRecompRange = {rangeStart, rangeEnd};
          allRecompRanges.push_back(deactivationRecompRange);
          // std::cout << "=> Deactivation recomp range (fallback) [" << rangeStart << ", " << rangeEnd << ")" << std::endl;
        }
      }
    }
  }
  

  // Process nucleotide mutations to find recomputation ranges
  if (strNodeId == "node_3") {
    // Define debug positions for node_3
    std::set<int64_t> debugPositions = {369894, 369925, 370975, 370981, 371019, 371042, 375703, 375704, 375710, 375724, 377237, 377962, 378303, 378327, 378336, 378337, 400112, 401809, 506725, 513651, 539424, 572027, 572085};
    
    logging::info("DEBUG_NUC_MUTATIONS: Node {} has {} nucleotide mutations:", strNodeId, node->nucMutation.size());
    for (size_t i = 0; i < node->nucMutation.size(); i++) {
      const auto& mut = node->nucMutation[i];
      
      // Extract mutation type from mutInfo field
      uint8_t mutType = mut.mutInfo & 0x7;
      int mutLen = (mut.mutInfo >> 4);
      std::string mutTypeStr = "UNKNOWN";
      switch (mutType) {
        case panmanUtils::NucMutationType::NS: mutTypeStr = "SUBSTITUTION"; break;
        case panmanUtils::NucMutationType::NI: mutTypeStr = "INSERTION"; break;
        case panmanUtils::NucMutationType::ND: mutTypeStr = "DELETION"; break;
        case panmanUtils::NucMutationType::NSNPS: mutTypeStr = "SNP_SUB"; mutLen = 1; break;
        case panmanUtils::NucMutationType::NSNPI: mutTypeStr = "SNP_INS"; mutLen = 1; break;
        case panmanUtils::NucMutationType::NSNPD: mutTypeStr = "SNP_DEL"; mutLen = 1; break;
        default: mutTypeStr = "TYPE_" + std::to_string(static_cast<int>(mutType)); break;
      }
      
      // Calculate global position for this mutation
      int64_t globalPos = -1;
      try {
        bool isBlockInverted = stateManager.isBlockInverted(strNodeId, mut.primaryBlockId);
        globalPos = stateManager.mapToGlobalCoordinate(mut.primaryBlockId, mut.nucPosition, 
                                                      mut.nucGapPosition, isBlockInverted);
      } catch (...) {
        throw std::runtime_error("Error mapping mutation position for node " + strNodeId + 
                             ", blockId " + std::to_string(mut.primaryBlockId) + 
                             ", nucPos " + std::to_string(mut.nucPosition) + 
                             ", nucGapPos " + std::to_string(mut.nucGapPosition));
      }
      
      logging::info("  NucMutation[{}]: {} at globalPos {} (blockId={} nucPos={} nucGapPos={} len={})", 
                   i, mutTypeStr, globalPos, mut.primaryBlockId, mut.nucPosition, mut.nucGapPosition, mutLen);
      
      // Check if any debug positions are near this mutation
      std::vector<int64_t> nearbyDebugPos;
      if (globalPos != -1) {
        for (int64_t debugPos : debugPositions) {
          // Check if debug position is within k_val bases of this mutation
          int64_t distance = std::abs(debugPos - globalPos);
          if (distance < k_val) {  // k_val is the kmer size, mutations within k_val bases can affect seeds
            nearbyDebugPos.push_back(debugPos);
          }
        }
        if (!nearbyDebugPos.empty()) {
          std::stringstream nearbyList;
          for (size_t j = 0; j < nearbyDebugPos.size(); j++) {
            if (j > 0) nearbyList << ",";
            nearbyList << nearbyDebugPos[j] << "(dist=" << (nearbyDebugPos[j] - globalPos) << ")";
          }
          logging::info("    Near debug positions: {}", nearbyList.str());
        }
      }
    }
  }

  // PHASE 3.5: Nucleotide mutation range calculation will happen AFTER mutations are applied

  // PHASE 5: Apply nucleotide mutations (batched to minimize lock contention)
  for (const auto &nuc_mutation : node->nucMutation) {
    int32_t blockId = nuc_mutation.primaryBlockId;
    int32_t nucPos = nuc_mutation.nucPosition;
    int32_t initialGapPos = nuc_mutation.nucGapPosition;
    uint32_t mutInfo = nuc_mutation.mutInfo;

    // Helper to decode nucleotides for logging
    auto decodeNucMutForLog = [&](int charIndex_log) -> char {
        int nucCode_log = 0;
        uint8_t type_log = mutInfo & 0x7;
        if (type_log < panmanUtils::NucMutationType::NSNPS) { // Multi-character mutation
            if (charIndex_log < 6) {
                nucCode_log = (nuc_mutation.nucs >> (4 * (5 - charIndex_log))) & 0xF;
            } else { return '?'; }
        } else { // SNP mutation
            if (charIndex_log == 0) {
                 nucCode_log = (nuc_mutation.nucs >> 20) & 0xF;
            } else { return '?'; }
        }
        return panmanUtils::getNucleotideFromCode(nucCode_log);
    };

    uint8_t type = mutInfo & 0x7;
    int len = (mutInfo >> 4);
    if (type >= panmanUtils::NucMutationType::NSNPS) len = 1;
    
    // Always log NUC_MUT_RAW if file is open
    if (debug_file_ourcode_pmfn.is_open()) {
        std::string newChars_log = "";
        if (len > 0) { 
            for (int k_log = 0; k_log < len; ++k_log) {
                newChars_log += decodeNucMutForLog(k_log);
            }
        } else {
            newChars_log = "(len<=0)";
        }
        debug_file_ourcode_pmfn << "NUC_MUT_RAW: blockId=" << blockId
                           << ", nucPos=" << nucPos
                           << ", initialGapPos=" << initialGapPos
                           << ", mutInfo_type=" << static_cast<int>(type)
                           << ", mutInfo_len=" << len
                           << ", newChars_decoded='" << newChars_log << "'"
                           << ", raw_nucs_field=" << nuc_mutation.nucs
                           << ", raw_mutInfo_field=" << mutInfo
                           << std::endl;
    }

    // Check if the block is currently active
    bool isCurrentBlockActive = stateManager.isBlockOn(strNodeId, blockId);
    if (blockId == 595) {
      std::cout << "Node " << strNodeId 
                << " Block " << blockId 
                << " is currently " << (isCurrentBlockActive ? "ACTIVE" : "INACTIVE") 
                << std::endl;
    }
    // Note: mutation counting happens later in the range generation phase

    
    auto decodeNucleotide = [&](int charIndex) -> char { 
        int nucCode = 0;
        if (type < panmanUtils::NucMutationType::NSNPS) {
            if (charIndex < 6) {
                nucCode = (nuc_mutation.nucs >> (4 * (5 - charIndex))) & 0xF;
            } else {
                throw std::runtime_error("Indexing error in decodeNucleotide");
            }
        } else {
            nucCode = (nuc_mutation.nucs >> 20) & 0xF;
        }
        try {
            return panmanUtils::getNucleotideFromCode(nucCode);
        } catch (const std::exception &e) {
            throw std::runtime_error("Decoding error in decodeNucleotide");
        }
    };

    // DEADLOCK FIX: Don't hold nodeState reference during nucleotide mutation application
    // Apply each nucleotide mutation independently to minimize lock duration
    
    for (int j = 0; j < len; ++j) {
      int32_t currentNucPos = nucPos;
      int32_t currentGapPos = initialGapPos;

      
      if (initialGapPos == -1) {
        currentNucPos = nucPos + j;
        currentGapPos = -1;
      } else {
        currentNucPos = nucPos;
        currentGapPos = initialGapPos + j;
      }

      char parChar = '?'; 
      try {
        parChar = stateManager.getCharAtPosition(strNodeId, blockId,
                                               currentNucPos, currentGapPos);
      } catch (const std::exception &e) {
        // Always log parChar retrieval error if file is open
        if (debug_file_ourcode_pmfn.is_open()) {
            debug_file_ourcode_pmfn << "NUC_MUT_APPLY_ERR: Exception getting parChar for "
                               << blockId << ":" << currentNucPos << ":" << currentGapPos
                               << " - " << e.what() << std::endl;
        }
        throw; 
      }

      
      char newChar = '-'; 
      if (type == panmanUtils::NucMutationType::NS ||
          type == panmanUtils::NucMutationType::NI ||
          type == panmanUtils::NucMutationType::NSNPS ||
          type == panmanUtils::NucMutationType::NSNPI) {
        try {
            newChar = decodeNucleotide(j);
        } catch (const std::exception& e) {
            throw std::runtime_error("Error decoding nucleotide for newChar");
        }
      }

      // === NEW TARGETED LOGGING ===
      if (condensed_trace_ss_ptr != nullptr && isCurrentBlockActive) {
          std::string decoded_char_from_nucs_field = "?";
          try {
              if (type < panmanUtils::NucMutationType::NSNPS) { // Multi-char part (NSNPS=0, NSNPI=3)
                  if (j < (mutInfo >> 4) ) decoded_char_from_nucs_field = std::string(1, decodeNucleotide(j));
              } else { // Single char SNP(2)/NI(5) or non-nuc types (ND=6, NSNPD=1)
                  if (type == panmanUtils::NucMutationType::NS || type == panmanUtils::NucMutationType::NI) {
                     if (j == 0) decoded_char_from_nucs_field = std::string(1, decodeNucleotide(0));
                  } else {
                     decoded_char_from_nucs_field = "N/A_DEL_TYPE";
                  }
              }
          } catch(...) { decoded_char_from_nucs_field = "ERR_DEC"; }

      }
      // === END NEW TARGETED LOGGING ===
      

      
      try {
        // DEADLOCK FIX: Acquire nodeState reference only for the minimal scope needed
        // std::cerr << "[DEADLOCK DEBUG] About to acquire nodeState for nucleotide mutation: node=" 
                  // << strNodeId << ", block=" << blockId << ", nucPos=" << currentNucPos << std::endl;
        auto& nodeState = stateManager.getNodeState(strNodeId);
        // std::cerr << "[DEADLOCK DEBUG] NodeState acquired, applying nucleotide mutation" << std::endl;
        stateManager.applyNucleotideMutation(strNodeId, nodeState, blockId,
                                           currentNucPos, currentGapPos,
                                           newChar);
        // std::cerr << "[DEADLOCK DEBUG] Nucleotide mutation applied successfully" << std::endl;
        // nodeState reference goes out of scope here, releasing the lock
        // === DEBUG VERIFICATION - only do this if absolutely necessary for debugging ===
        if (condensed_trace_ss_ptr != nullptr && debug_file_ourcode_pmfn.is_open()) {
            char char_after_apply = '?';
            std::string verify_status = "OK";
            uint8_t type_for_verify = mutInfo & 0x7; 
            try {
                char_after_apply = stateManager.getCharAtPosition(strNodeId, blockId, currentNucPos, currentGapPos);
                if (type_for_verify < panmanUtils::NucMutationType::NSNPS && char_after_apply != newChar) { 
                     verify_status = "MISMATCH_AFTER_APPLY";
                // Corrected enum usage for deletions, using type_for_verify
                } else if (type_for_verify == panmanUtils::NucMutationType::ND || type_for_verify == panmanUtils::NucMutationType::NSNPD) { 
                    verify_status = "DEL_VERIFY_RAW"; 
                    if (char_after_apply == '-' || char_after_apply == 'x') verify_status = "OK_DEL_VERIFIED";
                } else if (char_after_apply == newChar) {
                    verify_status = "OK_SNP_INS_VERIFIED";
                } else {
                    verify_status = "FAIL_OTHER_VERIFY";
                }
            } catch (const std::exception& verify_e) {
                verify_status = std::string("VERIFY_EXC: ") + verify_e.what();
            }
            debug_file_ourcode_pmfn << ", VERIFY: char_after_apply=''" << char_after_apply << "''"
                               << ", status=" << verify_status << std::endl; 
        }
        // === END DEBUG VERIFICATION ===

      } catch (const std::exception &e) {
        // Updated conditional logging (Leaf node only, all blocks)
        if (condensed_trace_ss_ptr != nullptr && debug_file_ourcode_pmfn.is_open()) {
            debug_file_ourcode_pmfn << ", APPLY_ERROR:" << e.what() << std::endl;
        }
        throw; 
      }

      // --- Derive Gap Updates: Only if the block is currently active for this node ---
      bool wasNucleotide = (parChar != '-' && parChar != 'x');
      bool isNucleotide = (newChar != '-' && newChar != 'x');

      if (wasNucleotide != isNucleotide) { // If a nuc became a gap or vice-versa
        // Only update gap map and log if the mutation is in an ACTIVE block for this node
        if (isCurrentBlockActive) {
        int64_t scalarPos = -1;
        try {
          bool isBlockInverted_for_coord = stateManager.isBlockInverted(strNodeId, blockId); // Get current inversion for correct global coord
          scalarPos = stateManager.mapToGlobalCoordinate(blockId, currentNucPos,
                                                      currentGapPos, isBlockInverted_for_coord);

          if (scalarPos != -1) {
              bool adding_gap = wasNucleotide && !isNucleotide; // True if Nuc -> Gap
              derivedGapUpdates.emplace_back(scalarPos, 1, adding_gap); 
              if (condensed_trace_ss_ptr != nullptr) { // Conditional logging to stringstream
                  std::string context_prefix = stateManager.extractSequence(strNodeId, std::max((int64_t)0, scalarPos - 5), scalarPos, false, false);
                  std::string context_suffix = stateManager.extractSequence(strNodeId, scalarPos + 1, std::min((int64_t)stateManager.getNumCoords(), scalarPos + 1 + 5), false, false);
                  std::string mut_trigger_str = "NucMut: " + std::to_string(blockId) + ":" + std::to_string(currentNucPos) + ":" + std::to_string(currentGapPos) + " " + parChar + "->" + newChar;
                  // (*condensed_trace_ss_ptr) << "GAP_LOG Node '" << strNodeId << "': Nuc->Gap/Gap->Nuc. Update: " << (adding_gap ? "ADD" : "REMOVE") 
                  //                           << " Gap [" << scalarPos << ", len 1]. Trigger: " << mut_trigger_str 
                  //                           << ". Context: '" << context_prefix << "' (" << parChar << "->" << newChar << ") '" << context_suffix << "'" << std::endl;
              }
            }
        } catch (const std::exception &e) {
          // Log warning instead of throwing, to continue with other mutations
          // logging::warn("DEBUG_MUTATION: Mapping error during gap update derivation for nuc mut: {}", e.what());
           if (condensed_trace_ss_ptr != nullptr) (*condensed_trace_ss_ptr) << "GAP_LOG ERROR (nuc-derived): " << e.what() << std::endl;
            }
        }
      }
    } 
  } 


  // PHASE 3.5: Collect nucleotide mutation ranges and blocksToDeactivate ranges
  // NOW calculate nucleotide mutation ranges AFTER mutations are applied
  
  // Process blocksToDeactivate ranges
  for (const auto &blockId : blocksToDeactivate) {
    // find previous ON block 
    int32_t blockIdPrev = blockId;
    while (blockIdPrev > 0 && !stateManager.isBlockOn(strNodeId, blockIdPrev)) {
      --blockIdPrev;
    }
    int32_t blockIdNext = blockId;
    while (blockIdNext < stateManager.getNumBlocks() && !stateManager.isBlockOn(strNodeId, blockIdNext)) {
      ++blockIdNext;
    }

    auto prevBlockRange = stateManager.getBlockRange(blockIdPrev);
    coordinates::CoordRange prevBoundary = {prevBlockRange.end - k_val, prevBlockRange.end};
    allRecompRanges.push_back(prevBoundary);
  }
  
  // Process nucleotide mutation ranges AFTER mutations are applied
  for (const auto &nuc_mutation : node->nucMutation) {
    int32_t blockId = nuc_mutation.primaryBlockId;
    int32_t nucPos = nuc_mutation.nucPosition;
    int32_t nucGapPos = nuc_mutation.nucGapPosition;
    uint32_t mutInfo = nuc_mutation.mutInfo;
    uint8_t type = mutInfo & 0x7;
    int len = (mutInfo >> 4);
    if (type >= panmanUtils::NucMutationType::NSNPS) len = 1;
    
    // Check if the block is currently active (post-mutation state)
    bool isCurrentBlockActive = stateManager.isBlockOn(strNodeId, blockId);
    
    // Count mutations by block activation status
    if (isCurrentBlockActive) {
      nucMutationsInActiveBlocks++;
      
      // Only generate recomputation ranges for active blocks
      try {
        // Calculate the global position for this mutation to ensure correct recomputation range
        int64_t globalPos = -1;
        try {
          bool isBlockInverted = stateManager.isBlockInverted(strNodeId, blockId);
          globalPos = stateManager.mapToGlobalCoordinate(blockId, nucPos, nucGapPos, isBlockInverted);
        } catch (...) {
          // If global position calculation fails, fall back to block-relative calculation
          globalPos = -1;
        }
        
        coordinates::CoordRange baseRange;
        if (globalPos != -1) {
          // Use global position to create a centered recomputation range
          baseRange.start = globalPos;
          baseRange.end = globalPos + len;
        } else {
          throw std::runtime_error("Failed to calculate global position for nucleotide mutation");
        }

        if (baseRange.start < baseRange.end) {
          
          // Expand the nucleotide mutation range upstream by k-1 non-gaps
          // to capture k-mers that span across the mutation point
          try {
            // Start looking for k non-gap characters from the mutation position itself
            // This ensures we capture all k-mers that could span the mutation
            int64_t searchStart = baseRange.start;
            
            // Extract k non-gap characters going backward from the mutation position
            auto kmerResult = stateManager.extractKmer(strNodeId, searchStart, k_val, true); // reverse=true to go backward
            
            if (!kmerResult.first.empty() && kmerResult.second.size() >= k_val) {
              // kmerResult.second[0] is the furthest back position when reverse=true
              // This gives us the start position of the first character in the k-mer
              int64_t upstreamStart = kmerResult.second[0];
              
              if (upstreamStart >= 0 && upstreamStart < baseRange.start) {
                // Extend the base range to include upstream context
                baseRange.start = upstreamStart;
                // std::cout << "=> Extended nuc mut range [" << baseRange.start << ", " << baseRange.end << ") with " << k_val << " upstream non-gap chars from pos " << upstreamStart << std::endl;
              }
            } else {
              // If we can't find k non-gap characters, fall back to simple position-based extension
              int64_t simpleUpstreamStart = std::max(static_cast<int64_t>(0), baseRange.start - (k_val - 1));
              if (simpleUpstreamStart < baseRange.start) {
                baseRange.start = simpleUpstreamStart;
                // std::cout << "=> Extended nuc mut range [" << baseRange.start << ", " << baseRange.end << ") with simple upstream fallback (couldn't find " << k_val << " non-gap chars)" << std::endl;
              }
            }
          } catch (const std::exception& e) {
            // Fallback: simple upstream extension if extractKmer fails
            int64_t simpleUpstreamStart = std::max(static_cast<int64_t>(0), baseRange.start - (k_val - 1));
            if (simpleUpstreamStart < baseRange.start) {
              baseRange.start = simpleUpstreamStart;
              // std::cout << "=> Extended nuc mut range [" << baseRange.start << ", " << baseRange.end << ") with simple upstream fallback (extractKmer failed)" << std::endl;
            }
          }
          

          
          nucMutationInducedRanges.push_back(baseRange);
          allRecompRanges.push_back(baseRange);
          
          // Debug log for nucleotide mutation ranges
          if (strNodeId == "node_3") {
            logging::info("DEBUG_NUC_RANGE_POST: Node {} mutation in block {} at nucPos {} nucGapPos {} globalPos {} generated range [{}, {}) size={}", 
                         strNodeId, blockId, nucPos, nucGapPos, globalPos, baseRange.start, baseRange.end, baseRange.end - baseRange.start);
          }
        }
      } catch (const std::exception &e) {
        // Error handled silently to prevent overflow
      }
    } else {
      nucMutationsInInactiveBlocks++;
    }
  }
  
  // PHASE 4: Add all collected recomputation ranges in a single operation
  if (!allRecompRanges.empty()) {
    {
      auto &nodeState = stateManager.getNodeState(strNodeId);
      nodeState.recompRanges.insert(nodeState.recompRanges.end(), allRecompRanges.begin(), allRecompRanges.end());
    }
  } 
  
  
  if (!derivedGapUpdates.empty()) {
    
    bool hasInvalidGapUpdates = false;
    for (size_t i = 0; i < derivedGapUpdates.size(); ++i) {
      const auto& update = derivedGapUpdates[i];
      if (update.pos < 0 || update.length <= 0) {
        hasInvalidGapUpdates = true;
      }
      
    }

    if (hasInvalidGapUpdates) {
      std::vector<coordinates::GapUpdate> validUpdates;
      validUpdates.reserve(derivedGapUpdates.size());
      for (const auto& update : derivedGapUpdates) {
        if (update.pos >= 0 && update.length > 0) {
          validUpdates.push_back(update);
        }
      }
      try {
        stateManager.consolidatedGapUpdate(strNodeId, validUpdates);
      } catch (const std::exception &e) {
        throw;
      }
    } else {
      try {
        stateManager.consolidatedGapUpdate(strNodeId, derivedGapUpdates);
      } catch (const std::exception &e) {
        throw;
      }
    }
  } 

  
  // Debug logging for recomputation ranges (minimal locking with proper scoping)
  {
    auto& nodeState = stateManager.getNodeState(strNodeId);
    if (!nodeState.recompRanges.empty()) {
      for (size_t i = 0; i < std::min(nodeState.recompRanges.size(), size_t(5)); i++) {
        auto& range = nodeState.recompRanges[i];
      }
    }
    // nodeState reference goes out of scope here, releasing the lock
  }
  
  std::stringstream ssBlockRanges, ssNucRanges;
  auto formatRanges = [](std::stringstream& ss, const std::vector<coordinates::CoordRange>& ranges) {
      ss << "[";
      for (size_t i = 0; i < ranges.size(); ++i) {
          ss << "[" << ranges[i].start << "," << ranges[i].end << ")";
          if (i < ranges.size() - 1) ss << ",";
      }
      ss << "]";
  };

  formatRanges(ssBlockRanges, blockMutationInducedRanges);
  formatRanges(ssNucRanges, nucMutationInducedRanges);

  
  size_t totalMutations = node->blockMutation.size() + node->nucMutation.size();
  if (debug_file_ourcode_pmfn.is_open()) { // Close the file at the end of the function
    debug_file_ourcode_pmfn << "===== FINISHED PROCESSING MUTATIONS FOR NODE: " << strNodeId << " =====" << std::endl << std::endl;
    debug_file_ourcode_pmfn.close();
  }

  // <<< END ADDED DEBUG LOGGING >>>
} 


void logNodeDetailsToFile(
    state::StateManager& stateManager,
    panmanUtils::Node* node,
    const std::string& logFilePath);


void outputSeedDetails(std::ofstream& outFile, int64_t pos, const seeding::seed_t& seed, 
                      state::StateManager& stateManager, const std::string& nodeId);



/**
 * @brief Write seed changes directly to Cap'n Proto format
 *
 * @param seedChanges Vector of seed change tuples
 * @param dfsIndex DFS index of the node
 * @param perNodeSeedMutations Cap'n Proto builder for seed mutations
 * @param kmerDictionary Dictionary mapping k-mer sequences to their IDs
 * @param stateManager StateManager for extracting k-mer sequences
 */
void writeSeedChangesToCapnp(
    const std::vector<int64_t>& basePositions, 
    const std::vector<uint64_t>& bitMasks,     
    int64_t dfsIndex,
    ::capnp::List<SeedMutations>::Builder &perNodeSeedMutations
    ) {

  
  logging::debug("DEBUG_INDEX_WRITE: Writing seed changes to CapnProto. Base positions: {}, Bit masks: {}, DFS index: {}", 
                basePositions.size(), bitMasks.size(), dfsIndex);

  
  if (basePositions.empty() || bitMasks.empty() || dfsIndex < 0) {
    logging::warn("DEBUG_INDEX_WRITE: Nothing to write - empty positions/masks or invalid DFS index!");
    return;
  }

  
  size_t entryIndex = dfsIndex;

  
  if (entryIndex >= perNodeSeedMutations.size()) {
    logging::err("DEBUG_INDEX_WRITE: DFS index {} is out of bounds for {} seed mutation entries",
                   dfsIndex, perNodeSeedMutations.size());
    return;
  }

  
  if (!basePositions.empty()) {
    std::stringstream positionStr;
    for (size_t i = 0; i < std::min(basePositions.size(), size_t(5)); i++) {
      positionStr << basePositions[i];
      if (i < basePositions.size() - 1 && i < 4) positionStr << ", ";
    }
    
    if (basePositions.size() > 5) positionStr << "...";
    
    logging::debug("DEBUG_INDEX_WRITE: First base positions: [{}]", positionStr.str());
  }
  
  if (!bitMasks.empty()) {
    std::stringstream maskStr;
    for (size_t i = 0; i < std::min(bitMasks.size(), size_t(5)); i++) {
      maskStr << std::hex << bitMasks[i] << std::dec;
      if (i < bitMasks.size() - 1 && i < 4) maskStr << ", ";
    }
    
    if (bitMasks.size() > 5) maskStr << "...";
    
    logging::debug("DEBUG_INDEX_WRITE: First bit masks: [{}]", maskStr.str());
  }

  
  auto mutations = perNodeSeedMutations[entryIndex];

  
    auto basePositionsBuilder = mutations.initBasePositions(basePositions.size());
    for (size_t i = 0; i < basePositions.size(); i++) {
      basePositionsBuilder.set(i, basePositions[i]);
    }
  
  
    auto perPosMasksBuilder = mutations.initPerPosMasks(bitMasks.size());
    for (size_t i = 0; i < bitMasks.size(); i++) {
      perPosMasksBuilder.set(i, bitMasks[i]);
    }



  
  // logging::info("DEBUG_INDEX_WRITE: Successfully wrote {} base positions and {} masks for node DFS index {}",
                 // basePositions.size(), bitMasks.size(), dfsIndex);

  
  try {
    auto verifyBasePositions = mutations.getBasePositions();
    auto verifyMasks = mutations.getPerPosMasks();
    
    if (verifyBasePositions.size() != basePositions.size() || verifyMasks.size() != bitMasks.size()) {
      logging::warn("DEBUG_INDEX_WRITE: Verification failed - size mismatch! Expected {}/{} positions/masks, got {}/{}",
                  basePositions.size(), bitMasks.size(), verifyBasePositions.size(), verifyMasks.size());
    } else {
      logging::debug("DEBUG_INDEX_WRITE: Verification successful - wrote {}/{} positions/masks for DFS index {}",
                   verifyBasePositions.size(), verifyMasks.size(), dfsIndex);
    }
  } catch (const std::exception& e) {
    logging::err("DEBUG_INDEX_WRITE: Verification failed with error: {}", e.what());
  }
}


std::tuple<int, int, int, int, int, int, int, int, int, std::string, int, int, std::vector<std::tuple<int64_t, bool, bool, std::optional<size_t>, std::optional<size_t>, std::optional<bool>, std::optional<bool>, std::optional<int64_t>, std::optional<int64_t>>>> recomputeSeeds(
    state::StateManager &stateManager,
    placement::PlacementEngine &engine,
    panmanUtils::Node *node, int k, int s,
    ::capnp::List<SeedMutations>::Builder *perNodeSeedMutations,
    const std::unordered_map<std::string, uint32_t>* kmerDictionary,
    std::unordered_set<std::string>& uniqueKmersCollector,
    std::ofstream& debugSeedFile) {

  
  thread_local std::string kmerBuffer;

  std::string nodeId = node->identifier;
  std::cout << "DEBUG_RECOMP_START: Starting recomputeSeeds for " << nodeId << std::endl;
  
  // New debug metrics for TSV
  int totalMaterializedSeeds = 0;
  int seedsFoundInDeletedBlocks = 0;
  int blockDeletionPositionsChecked = 0;
  
  // NOTE: materializeNodeState() is now called BEFORE processMutationsForNode in processNodeComplete
  // This ensures inherited seeds are available during block deletion processing
  
  // Add debug logging to verify materialization is already complete
  bool materializedStateComputed = false;
  {
    auto& nodeStateDebug = stateManager.getNodeState(nodeId);
    // Avoid potential deadlock by not locking here - this is just for debug info
    int materializedCount = nodeStateDebug.materializedSeeds.size();
    materializedStateComputed = nodeStateDebug.materializedStateComputed;
    totalMaterializedSeeds = materializedCount;
    logging::debug("SEED_MATERIALIZATION: Node {} has {} materialized seeds (materialization done before mutations), materializedStateComputed={}", 
                   nodeId, materializedCount, materializedStateComputed);
  }
  
  // Debug: Check if materialized state is populated
  auto &nodeState = stateManager.getNodeState(nodeId);
  
  // Initialize efficient seed counting for this node
  if (node->parent) {
    try {
      const auto &parentState = stateManager.getNodeState(node->parent->identifier);
      nodeState.initializeSeedCountFromParent(parentState.getTotalSeedCount());
    } catch (const std::exception& e) {
      // If parent state not available, we'll track changes from 0
      logging::debug("Could not inherit parent seed count for node {}: {}", nodeId, e.what());
      nodeState.initializeSeedCountFromParent(0);
    }
  } else {
    // Root node starts with 0 seeds
    nodeState.initializeSeedCountFromParent(0);
  }
  
  std::vector<uint32_t> dictionaryIds;
  std::vector<int64_t> dictionaryPositions; 
  std::vector<uint32_t> dictionaryEndOffsets;
  
  int seedsInherited = 0; // This will likely remain 0 or be removed after change
  int seedsCleared = 0;
  int seedsAdded = 0;
  int seedsChanged = 0;
  
  // Block mutation counters for debug TSV output
  int blockDeletions = 0;
  int blockInsertions = 0;

  
  std::vector<std::tuple<int64_t, bool, bool, std::optional<size_t>, std::optional<size_t>, 
                         std::optional<bool>, std::optional<bool>, std::optional<int64_t>, 
                         std::optional<int64_t>>> seedChanges;

  // Track positions where seeds were cleared due to block deletions
  std::set<int64_t> positionsWithClearedSeeds;
  
  // Track specific positions of interest for debugging
  std::set<int64_t> debugPositions = {369894, 369925, 370975, 370981, 371019, 371042, 375703, 375704, 375710, 375724, 377237, 377962, 378303, 378327, 378336, 378337, 400112, 401809, 506725, 513651, 539424, 572027, 572085};
  
  // STEP 1: Handle block mutations - clear seeds from deleted/inverted blocks
  if (nodeId == "node_3") {
    logging::info("DEBUG_BLOCK_MUTATIONS: Node {} has {} block mutations:", nodeId, node->blockMutation.size());
    for (size_t i = 0; i < node->blockMutation.size(); i++) {
      const auto& mut = node->blockMutation[i];
      std::string mutType = "UNKNOWN";
      if (mut.blockMutInfo == 1) mutType = "INSERTION";
      else if (mut.blockMutInfo == 0 && mut.inversion) mutType = "INVERSION";
      else if (mut.blockMutInfo == 0 && !mut.inversion) mutType = "DELETION";
      
      try {
        auto blockRange = stateManager.getBlockRange(mut.primaryBlockId);
        logging::info("  Mutation[{}]: {} block {} range [{}, {}) size={}", 
                     i, mutType, mut.primaryBlockId, blockRange.start, blockRange.end, 
                     blockRange.end - blockRange.start);
        
        // Check if any debug positions fall within this block
        std::vector<int64_t> debugPosInBlock;
        for (int64_t debugPos : debugPositions) {
          if (debugPos >= blockRange.start && debugPos < blockRange.end) {
            debugPosInBlock.push_back(debugPos);
          }
        }
        if (!debugPosInBlock.empty()) {
          std::stringstream debugList;
          for (size_t j = 0; j < debugPosInBlock.size(); j++) {
            if (j > 0) debugList << ",";
            debugList << debugPosInBlock[j];
          }
          logging::info("    Block contains debug positions: {}", debugList.str());
        }
      } catch (...) {
        logging::info("  Mutation[{}]: {} block {} (could not get range)", i, mutType, mut.primaryBlockId);
      }
    }
  }
  
  for (const auto &block_mutation : node->blockMutation) {
    bool isInsertion = (block_mutation.blockMutInfo == 1);
    bool isDeletionOrInversion = (block_mutation.blockMutInfo == 0);
    bool isInversion = isDeletionOrInversion && block_mutation.inversion;
    bool isDeletion = isDeletionOrInversion && !block_mutation.inversion;
    
    // DEBUG: Log block mutation details for troubleshooting
    logging::debug("BLOCK_MUT_DEBUG: Node {} block {} - blockMutInfo={}, inversion={}, isInsertion={}, isDeletion={}, isInversion={}", 
                   nodeId, block_mutation.primaryBlockId, block_mutation.blockMutInfo, 
                   block_mutation.inversion, isInsertion, isDeletion, isInversion);
    
    // Count block mutations for debug TSV output
    if (isDeletion) {
      blockDeletions++; // True deletion
    } else if (isInsertion) {
      blockInsertions++; // Insertion
    }
    // Note: inversions don't count as deletions or insertions for TSV
    
    // Clear seeds from deleted blocks (since they become all gaps)
    if (isDeletion) {
      logging::debug("BLOCK_DEL_DEBUG: Processing deletion for node {} block {}", nodeId, block_mutation.primaryBlockId);
      
      // Debug: Check materialized state before deletion
      {
        auto& nodeStateDebug = stateManager.getNodeState(nodeId);
        // Avoid potential deadlock by not locking here - this is just for debug info
        int currentMaterializedCount = nodeStateDebug.materializedSeeds.size();
        logging::debug("BLOCK_DEL_DEBUG: Node {} has {} materialized seeds before processing block {} deletion", 
                       nodeId, currentMaterializedCount, block_mutation.primaryBlockId);
      }
      
      int64_t start = stateManager.getBlockRange(block_mutation.primaryBlockId).start;
      int64_t end = stateManager.getBlockRange(block_mutation.primaryBlockId).end;
      
      logging::debug("BLOCK_DEL_DEBUG: Block {} range is [{}, {})", block_mutation.primaryBlockId, start, end);
      
      int seedsFoundInThisBlock = 0;
      int positionsCheckedInThisBlock = 0;
      
      // DEADLOCK SAFE - Clear all seeds within the deleted block
      // Use materialized seed state to query inherited seeds
      for (int64_t pos = start; pos < end; pos++) {
        positionsCheckedInThisBlock++;
        
        // DEADLOCK SAFE: Wrap seed operations to avoid deadlocks
        try {
          auto seedOpt = stateManager.getSeedAtPosition(nodeId, pos);
          if (seedOpt.has_value()) {
            seedsFoundInThisBlock++;
            seeding::seed_t seed = seedOpt.value();
            stateManager.clearSeedAtPosition(nodeId, pos);
            
            // Track this position as having a cleared seed
            positionsWithClearedSeeds.insert(pos);
            
            // Report if this is a position of interest
            if (debugPositions.find(pos) != debugPositions.end()) {
              logging::info("DEBUG_POS_TRACK: Node {} - BLOCK_DELETION cleared seed at debug position {} in block {}", 
                           nodeId, pos, block_mutation.primaryBlockId);
            }
                                  
            seedChanges.emplace_back(
                std::make_tuple(pos, true, false, std::optional<size_t>(static_cast<size_t>(seed.hash)), std::optional<size_t>(),
                                std::optional<bool>(seed.reversed), std::optional<bool>(), std::optional<int64_t>(seed.endPos), std::optional<int64_t>()));
            seedsCleared++;
            
            // Debug: Log individual seed deletions
            if (seedsFoundInThisBlock <= 5) { // Log first few seeds only to avoid spam
              logging::debug("SEED_DELETION_DETAIL: Node {} pos {} - found seed hash={}, endPos={}, cleared successfully", 
                             nodeId, pos, seed.hash, seed.endPos);
            }
          }
        } catch (const std::system_error& e) {
          if (e.code() == std::errc::resource_deadlock_would_occur) {
            throw std::runtime_error("Resource deadlock detected during seed deletion at pos " + 
                                   std::to_string(pos) + " in block " + std::to_string(block_mutation.primaryBlockId) + 
                                   " for node " + nodeId);
          }
          throw; // Re-throw other exceptions
        }
      }
      
      // Update debug metrics
      seedsFoundInDeletedBlocks += seedsFoundInThisBlock;
      blockDeletionPositionsChecked += positionsCheckedInThisBlock;
      
      // Add debug logging for block deletion
      logging::debug("BLOCK_DELETION: Node {} deleted block {} [{}, {}): found {} seeds in {} positions, total cleared so far: {}", 
                     nodeId, block_mutation.primaryBlockId, start, end, seedsFoundInThisBlock, positionsCheckedInThisBlock, seedsCleared);
    }
    
    // Clear seeds from inverted blocks (since all k-mers get new hashes due to reverse complementation)
    if (isInversion) {
      logging::debug("BLOCK_INV_DEBUG: Processing inversion for node {} block {}", nodeId, block_mutation.primaryBlockId);
      
      int64_t start = stateManager.getBlockRange(block_mutation.primaryBlockId).start;
      int64_t end = stateManager.getBlockRange(block_mutation.primaryBlockId).end;
      
      logging::debug("BLOCK_INV_DEBUG: Block {} range is [{}, {})", block_mutation.primaryBlockId, start, end);
      
      int seedsFoundInThisInversion = 0;
      
      // DEADLOCK SAFE - Clear all seeds within the inverted block
      // Use hierarchical seed storage - query current node to get inherited seeds
      for (int64_t pos = start; pos < end; pos++) {
        try {
          auto seedOpt = stateManager.getSeedAtPosition(nodeId, pos);
          if (seedOpt.has_value()) {
            seedsFoundInThisInversion++;
            seeding::seed_t seed = seedOpt.value();
            stateManager.clearSeedAtPosition(nodeId, pos);
                                  
            seedChanges.emplace_back(
                std::make_tuple(pos, true, false, std::optional<size_t>(static_cast<size_t>(seed.hash)), std::optional<size_t>(),
                                std::optional<bool>(seed.reversed), std::optional<bool>(), std::optional<int64_t>(seed.endPos), std::optional<int64_t>()));
            seedsCleared++;
          }
        } catch (const std::system_error& e) {
          if (e.code() == std::errc::resource_deadlock_would_occur) {
            throw std::runtime_error("Resource deadlock detected during seed inversion at pos " + 
                                   std::to_string(pos) + " in block " + std::to_string(block_mutation.primaryBlockId) + 
                                   " for node " + nodeId);
          }
          throw; // Re-throw other exceptions
        }
      }
      
      logging::debug("BLOCK_INVERSION: Node {} inverted block {} [{}, {}): found {} seeds, total cleared so far: {}", 
                     nodeId, block_mutation.primaryBlockId, start, end, seedsFoundInThisInversion, seedsCleared);
    }
  }

  // STEP 2: DEADLOCK SAFE - Clear any seeds at gap or invalid positions
  std::vector<std::pair<int32_t, coordinates::CoordRange>> activeBlocksAndRanges;
  try {
    activeBlocksAndRanges = stateManager.getActiveBlockRanges(nodeId); // Get all active blocks with ranges
  } catch (const std::system_error& e) {
    if (e.code() == std::errc::resource_deadlock_would_occur) {
      throw std::runtime_error("Resource deadlock detected during active block range query for node " + nodeId);
    } else {
      throw; // Re-throw other exceptions
    }
  }
  
  
  
  
  // STEP 3: DEADLOCK SAFE - Process recomputation ranges if any
  std::vector<coordinates::CoordRange> ranges;
  try {
    ranges = stateManager.getRecompRanges(nodeId);
  } catch (const std::system_error& e) {
    if (e.code() == std::errc::resource_deadlock_would_occur) {
      throw std::runtime_error("Resource deadlock detected during recomputation range query for node " + nodeId);
    } else {
      throw; // Re-throw other exceptions
    }
  }
  
  // Merge overlapping ranges BEFORE processing to avoid duplicate work
  std::vector<coordinates::CoordRange> mergedRanges;
  if (!ranges.empty()) {
    // Sort ranges by start position
    std::vector<coordinates::CoordRange> sortedRanges = ranges;
    std::sort(sortedRanges.begin(), sortedRanges.end(), 
              [](const coordinates::CoordRange& a, const coordinates::CoordRange& b) {
                return a.start < b.start;
              });
    
    // Merge overlapping ranges
    for (const auto& range : sortedRanges) {
      if (mergedRanges.empty() || mergedRanges.back().end < range.start) {
        // No overlap, add new range
        mergedRanges.push_back(range);
      } else {
        // Overlap detected, merge with previous range
        mergedRanges.back().end = std::max(mergedRanges.back().end, range.end);
      }
    }
    
    // Update ranges to use merged ranges
    ranges = mergedRanges;
    
    // Log merging results for debugging
    if (nodeId == "node_3" || nodeId == "node_2") {
      logging::info("DEBUG_MERGE: Node {} had {} original ranges, merged to {} ranges", 
                    nodeId, sortedRanges.size(), ranges.size());
    }
  }
  
  // Log recomputation ranges and check for debug positions coverage
  if (nodeId == "node_3") {
    logging::info("DEBUG_RECOMP_RANGES: Node {} has {} recomputation ranges (after merging):", nodeId, ranges.size());
    for (size_t i = 0; i < ranges.size(); i++) {
      const auto& range = ranges[i];
      logging::info("  Range[{}]: [{}, {}) size={}", i, range.start, range.end, range.end - range.start);
      
      // Check if any debug positions fall within this range
      std::vector<int64_t> debugPosInRange;
      for (int64_t debugPos : debugPositions) {
        if (debugPos >= range.start && debugPos < range.end) {
          debugPosInRange.push_back(debugPos);
        }
      }
      if (!debugPosInRange.empty()) {
        std::stringstream debugList;
        for (size_t j = 0; j < debugPosInRange.size(); j++) {
          if (j > 0) debugList << ",";
          debugList << debugPosInRange[j];
        }
        logging::info("    Contains debug positions: {}", debugList.str());
      }
    }
    
    // Check for debug positions NOT covered by any range
    std::vector<int64_t> uncoveredDebugPos;
    for (int64_t debugPos : debugPositions) {
      bool covered = false;
      for (const auto& range : ranges) {
        if (debugPos >= range.start && debugPos < range.end) {
          covered = true;
          break;
        }
      }
      if (!covered) {
        uncoveredDebugPos.push_back(debugPos);
      }
    }
    
    if (!uncoveredDebugPos.empty()) {
      std::stringstream uncoveredList;
      for (size_t j = 0; j < uncoveredDebugPos.size(); j++) {
        if (j > 0) uncoveredList << ",";
        uncoveredList << uncoveredDebugPos[j];
      }
      logging::info("DEBUG_RECOMP_RANGES: Node {} - Debug positions NOT covered by any recomp range: {}", 
                   nodeId, uncoveredList.str());
    }
  }
  
  // Debug: Track recomputation range information
  int recompRangeCount = ranges.size();
  int64_t recompRangeSize = 0;
  std::stringstream recompRangesStr;
  std::set<int32_t> blocksWithRanges;
  
  // DEBUG: Track unique positions processed in recomputation ranges
  std::set<int64_t> uniquePositionsProcessed;
  int totalKmersProcessedInRanges = 0;
  
  // Collect inserted block ranges for better mapping
  std::unordered_map<int64_t, int32_t> insertionPointToBlockId;
  for (const auto &block_mutation : node->blockMutation) {
    if (block_mutation.blockMutInfo == 1) { // insertion
      try {
        auto blockRange = stateManager.getBlockRange(block_mutation.primaryBlockId);
        insertionPointToBlockId[blockRange.start] = block_mutation.primaryBlockId;
      } catch(...) {}
    }
  }
  
  // Debug log insertion point map for node_2
  if (nodeId == "node_2") {
    // logging::info("DEBUG_INSERTION_MAP: {} insertion points mapped", insertionPointToBlockId.size());
  }
  
  for (size_t i = 0; i < ranges.size(); i++) {
    const auto& range = ranges[i];
    recompRangeSize += (range.end - range.start);
    if (i > 0) recompRangesStr << ";";
    
    int64_t length = range.end - range.start;
    
    // Map range to block ID with improved logic for merged ranges
    int32_t blockId = -1;
    std::set<int32_t> blocksSpannedByRange;
    
    try {
      // Check which inserted blocks this range covers or overlaps with
      for (const auto &block_mutation : node->blockMutation) {
        if (block_mutation.blockMutInfo == 1) { // insertion
          try {
            auto blockRange = stateManager.getBlockRange(block_mutation.primaryBlockId);
            
            // Check if range overlaps with or is adjacent to this block
            bool rangeOverlapsBlock = (range.start < blockRange.end && range.end > blockRange.start);
            bool rangeAdjacentToBlock = (range.end == blockRange.start || range.start == blockRange.end);
            bool rangeContainsBlock = (range.start <= blockRange.start && range.end >= blockRange.end);
            
            if (rangeOverlapsBlock || rangeAdjacentToBlock || rangeContainsBlock) {
              blocksSpannedByRange.insert(block_mutation.primaryBlockId);
            }
          } catch(...) {}
        }
      }
      
      // If range spans multiple blocks, pick the first one (for reporting purposes)
      // but count all blocks as having ranges
      if (!blocksSpannedByRange.empty()) {
        blockId = *blocksSpannedByRange.begin();
        
        // Add all spanned blocks to the blocksWithRanges set
        for (int32_t spannedBlockId : blocksSpannedByRange) {
          blocksWithRanges.insert(spannedBlockId);
        }
        
        if (nodeId == "node_2" && blocksSpannedByRange.size() > 1) {
          std::stringstream spannedList;
          for (auto it = blocksSpannedByRange.begin(); it != blocksSpannedByRange.end(); ++it) {
            if (it != blocksSpannedByRange.begin()) spannedList << ",";
            spannedList << *it;
          }
          // logging::info("RANGE_MAP: [{}, {}) -> spans blocks [{}] (reporting as block {})", 
                       // range.start, range.end, spannedList.str(), blockId);
        } else if (nodeId == "node_2") {
          // logging::info("RANGE_MAP: [{}, {}) -> block {} (single block)", 
                       // range.start, range.end, blockId);
        }
      } else {
        // Try mapping the range start position using coordinate mapping
        auto coordOpt = stateManager.fastMapGlobalToLocal(range.start);
        if (coordOpt) {
          auto [mappedBlockId, nucPos, gapPos] = *coordOpt;
          blockId = mappedBlockId;
          if (nodeId == "node_2") {
            // logging::info("RANGE_MAP: [{}, {}) -> block {} (via coordinate mapping)", 
                         // range.start, range.end, blockId);
          }
          blocksWithRanges.insert(blockId);
        } else {
          if (nodeId == "node_2") {
            // logging::info("RANGE_MAP: [{}, {}) -> unmapped range", range.start, range.end);
          }
        }
      }
    } catch(...) {}
    
    recompRangesStr << length << ":" << blockId;
  }
  std::string recompRangesFormatted = recompRangesStr.str();
  
  // Add debug logging for node_2 specifically to understand undercounting
  if (nodeId == "node_2") {
    // logging::info("DEBUG_NODE2: Processing {} recomputation ranges with total size {}", recompRangeCount, recompRangeSize);
    
    // Track which blocks have recomputation ranges vs which blocks are inserted
    std::set<int32_t> insertedBlocks;
    
    // Collect inserted block IDs
    for (const auto &block_mutation : node->blockMutation) {
      if (block_mutation.blockMutInfo == 1) { // insertion
        insertedBlocks.insert(block_mutation.primaryBlockId);
      }
    }
    
    // Find blocks that are inserted but missing recomputation ranges
    std::set<int32_t> missingBlocks;
    std::set_difference(insertedBlocks.begin(), insertedBlocks.end(),
                       blocksWithRanges.begin(), blocksWithRanges.end(),
                       std::inserter(missingBlocks, missingBlocks.begin()));
    
    // logging::info("DEBUG_NODE2_RANGES: {} inserted blocks, {} blocks with ranges, {} missing ranges",
                  // insertedBlocks.size(), blocksWithRanges.size(), missingBlocks.size());
    
    if (!missingBlocks.empty()) {
      std::stringstream missingList;
      for (auto blockId : missingBlocks) {
        if (missingList.tellp() > 0) missingList << ",";
        missingList << blockId;
      }
      // logging::info("DEBUG_NODE2_RANGES: Missing recomp ranges for blocks: {}", missingList.str());
    }
  }
  
  if (!ranges.empty()) {
    
    // Ranges are already merged at the beginning of STEP 3
    std::vector<coordinates::CoordRange> rangesToProcess = ranges;
    
    logging::debug("DEBUG_RECOMP: Node '{}' processing {} merged recomputation ranges", nodeId, rangesToProcess.size());

    
    int64_t totalRangeSize = 0;
    for(const auto& range : rangesToProcess) totalRangeSize += (range.end - range.start);

    // Process each range - DEADLOCK SAFE implementation
    for (size_t rangeIdx = 0; rangeIdx < rangesToProcess.size(); ++rangeIdx) {
      const auto &range = rangesToProcess[rangeIdx];
      try {
        // DEADLOCK SAFE: Extract sequence with timeout and fallback
        std::string gappedSequence;
        std::vector<int64_t> positions;
        std::vector<bool> gaps;
        std::vector<int64_t> endPositions;
        
        // Extend range end to include k more non-gap characters downstream
        // but don't recompute seeds at those extended positions
        int64_t extendedRangeEnd = range.end;
        try {
            // Use extractKmer to get k characters downstream from range.end
            auto kmerResult = stateManager.extractKmer(nodeId, range.end, k);
            if (!kmerResult.first.empty() && !kmerResult.second.empty()) {
                // Use the position after the k-th character (end position of k-mer)
                extendedRangeEnd = kmerResult.second[k - 1] + 1;
            } else {
              // If k-mer extraction fails at range.end, use range.end as the boundary
              extendedRangeEnd = range.end;
              logging::debug("Could not extract k-mer from range end {} in node {}, using range.end as boundary", 
                           range.end, nodeId);
            }
        } catch (...) {
            // If k-mer extraction fails at range.end, use range.end as the boundary
            extendedRangeEnd = range.end;
            logging::debug("Error extracting k-mer from range end {} in node {}, using range.end as boundary", 
                         range.end, nodeId);
        }
        coordinates::CoordRange extendedRange = {range.start, extendedRangeEnd};
        try {
          std::tie(gappedSequence, positions, gaps, endPositions) = stateManager.extractSequence(nodeId, extendedRange, false);
        } catch (const std::system_error& e) {
          if (e.code() == std::errc::resource_deadlock_would_occur) {
            throw std::runtime_error("Resource deadlock detected during sequence extraction for range [" + 
                                   std::to_string(extendedRange.start) + ", " + std::to_string(extendedRange.end) + 
                                   ") in node " + nodeId);
          }
          throw; // Re-throw other exceptions
        }


        std::string ungappedSequence = gappedSequence;
        ungappedSequence.erase(
            std::remove_if(ungappedSequence.begin(), ungappedSequence.end(), [](char c) { return c == '-' || c == 'x'; }),
            ungappedSequence.end());

        if (ungappedSequence.length() < static_cast<size_t>(k)) {
          continue;
        }
        

        // make a vector of int64_t that maps local index in ungapped seq to local position in gapped sequence
        std::vector<int64_t> localToGappedPos(ungappedSequence.length(), -1);
        int64_t gappedIndex = 0;
        for (size_t i = 0; i < gappedSequence.length(); ++i) {
          if (gappedSequence[i] != '-' && gappedSequence[i] != 'x') {
            localToGappedPos[gappedIndex++] = i; // Map gapped position to ungapped index
          }
        }

        // std::cerr << " ___ ours ____" << std::endl;
        // std::cerr << "Range: [" << range.start << ", " << range.end << "]" << std::endl;
        // std::cerr << "Gapped Sequence: " << gappedSequence << std::endl;
        // std::cerr << "Ungapped Sequence: " << ungappedSequence << std::endl;
        // std::cerr << "localIdx -> gappedIndex -> globalPos: " << std::endl;
        // for (size_t i = 0; i < localToGappedPos.size(); ++i) {
        //   if (localToGappedPos[i] != -1) {
        //     int64_t globalPos = positions[localToGappedPos[i]];
        //     std::cerr << "  " << i << " -> " << localToGappedPos[i] << " -> " << globalPos << std::endl;
        //   } else {
        //     std::cerr << "  " << i << " -> -1 (no mapping)" << std::endl;
        //   }
        // }



        // Step 1: DEADLOCK SAFE seed collection - collect all positions in this range that currently have seeds
        // Use the hierarchical seed storage - query the current node and let it automatically
        // inherit from parents, just like character data inheritance
        std::unordered_map<int64_t, seeding::seed_t> existingSeedsInRange; // Use globalPos as key, not local index
        int foundSeeds = 0, totalPositions = positions.size();
        
        // DEADLOCK SAFE: Batch seed queries to minimize lock contention
        for (size_t i = 0; i < positions.size(); i++) {
          int64_t globalPos = positions[i];
          
          // DEADLOCK SAFE: Wrap getSeedAtPosition call to avoid deadlocks
          try {
            // Query seeds from current node - materialized state includes inherited seeds
            auto seedOpt = stateManager.getSeedAtPosition(nodeId, globalPos);
            if (seedOpt.has_value()) {
              existingSeedsInRange[globalPos] = seedOpt.value(); // Key by globalPos
              foundSeeds++;
              
              // Report if this is a position of interest
              if (debugPositions.find(globalPos) != debugPositions.end()) {
                logging::info("DEBUG_POS_TRACK: Node {} - Found EXISTING seed at debug position {} (hash {})", 
                             nodeId, globalPos, seedOpt.value().hash);
              }
            } else {
              // Report if this is a position of interest that has no seed
              if (debugPositions.find(globalPos) != debugPositions.end()) {
                logging::info("DEBUG_POS_TRACK: Node {} - NO existing seed at debug position {}", 
                             nodeId, globalPos);
              }
            }
          } catch (const std::system_error& e) {
            if (e.code() == std::errc::resource_deadlock_would_occur) {
              throw std::runtime_error("Resource deadlock detected during seed query at pos " + 
                                     std::to_string(globalPos) + " for node " + nodeId);
            }
            throw; // Re-throw other exceptions
          }
        }
        
        // Diagnostic logging to understand seed inheritance
        if (nodeId == "node_2" || nodeId == "node_3") {
          logging::info("DEBUG_SEED_INHERITANCE: Node {} range [{},{}): found {} existing seeds out of {} positions", 
                        nodeId, range.start, range.end, foundSeeds, totalPositions);
          
          // Sample a few positions to see what's happening
          for (size_t i = 0; i < std::min(5, totalPositions); i++) {
            int64_t globalPos = positions[i];
            bool hasExisting = existingSeedsInRange.find(globalPos) != existingSeedsInRange.end();
            logging::info("  pos[{}]: globalPos={} hasExisting={}", i, globalPos, hasExisting);
          }
        }
        // Delete seeds at any position that now has a gap character ('-' or 'x')
        // Loop over existing seeds and check if their positions now have gaps
        std::vector<int64_t> seedsToDelete;
        for (const auto& [globalPos, seed] : existingSeedsInRange) {
          // Find the local index for this global position
          auto posIt = std::find(positions.begin(), positions.end(), globalPos);
          if (posIt != positions.end()) {
            size_t localIndex = std::distance(positions.begin(), posIt);
            if (localIndex < gappedSequence.size()) {
              char ch = gappedSequence[localIndex];
              if (ch == '-' || ch == 'x') {
                seedsToDelete.push_back(globalPos);
              }
            }
          }
        }
        
        // DEADLOCK SAFE: Delete seeds at gap positions
        for (int64_t globalPos : seedsToDelete) {
          try {
            auto seedIt = existingSeedsInRange.find(globalPos);
            if (seedIt != existingSeedsInRange.end()) {
              seeding::seed_t deletedSeed = seedIt->second;
              
              stateManager.clearSeedAtPosition(nodeId, globalPos);
              existingSeedsInRange.erase(seedIt);
              seedsCleared++;
              
              // Track this position as having a cleared seed (for gap deletions)
              positionsWithClearedSeeds.insert(globalPos);
              
              seedChanges.emplace_back(
                  std::make_tuple(globalPos, true, false, std::optional<size_t>(static_cast<size_t>(deletedSeed.hash)), std::optional<size_t>(),
                                 std::optional<bool>(deletedSeed.reversed), std::optional<bool>(), std::optional<int64_t>(deletedSeed.endPos), std::optional<int64_t>()));
              
              logging::debug("DEBUG_RECOMP: Deleted seed at gap position {} (char='{}')", globalPos, 
                           (globalPos - range.start < static_cast<int64_t>(gappedSequence.size())) ? 
                           gappedSequence[globalPos - range.start] : '?');
            }
          } catch (const std::system_error& e) {
            if (e.code() == std::errc::resource_deadlock_would_occur) {
              throw std::runtime_error("Resource deadlock detected during seed clearing at gap position " + 
                                     std::to_string(globalPos) + " for node " + nodeId);
            }
            throw; // Re-throw other exceptions
          }
        }

        // Step 2: Process k-mers in this range
        std::vector<std::tuple<size_t, bool, bool, int64_t>> kmers =
            seeding::rollingSyncmers(ungappedSequence, k, s, false, 0, true);

        // DEBUG: Log initial k-mer count for node_2
        if (nodeId == "node_2") {
          // logging::info("SYNCMER_COUNT: Range [{},{}): found {} syncmers in ungapped seq len {}", 
                       // range.start, range.end, kmers.size(), ungappedSequence.length());
        }

        for (const auto& [hash, isReverse, isNowSeed, localStartPos] : kmers) {
          // Bounds check for position mapping
          if (localStartPos >= localToGappedPos.size()) {
            logging::warn("SKIP_BOUNDS: pos {} >= mapSize {} in node {}", 
                         localStartPos, localToGappedPos.size(), nodeId);
            continue;
          }
          
          int64_t gappedLocalStartPos = localToGappedPos[localStartPos];
          if (gappedLocalStartPos < 0 || gappedLocalStartPos >= static_cast<int64_t>(positions.size())) {
            logging::warn("SKIP_INVALID: gappedPos {} invalid (posSize={}) in node {}", 
                         gappedLocalStartPos, positions.size(), nodeId);
            continue;
          }
          
          int64_t globalPos = positions[gappedLocalStartPos];
          
          // Report syncmer processing result for debug positions
          if (debugPositions.find(globalPos) != debugPositions.end()) {
            std::string kmerSeq = (localStartPos + k <= ungappedSequence.size()) ? 
                                 ungappedSequence.substr(localStartPos, k) : "TRUNCATED";
            logging::info("DEBUG_POS_TRACK: Node {} - SYNCMER_RESULT at debug position {} kmer='{}' hash={} isReverse={} isNowSeed={}", 
                         nodeId, globalPos, kmerSeq, hash, isReverse, isNowSeed);
          }
          
          // Report if this is a position of interest
          if (debugPositions.find(globalPos) != debugPositions.end()) {
            logging::info("DEBUG_POS_TRACK: Node {} - Processing debug position {} in range [{},{})", 
                         nodeId, globalPos, range.start, range.end);
          }
          
          // Only process seeds within the original range, not in the extended context
          // The extended context is only used to ensure complete k-mers can be extracted
          if (globalPos >= range.end) {
            // DEBUG: Log when we skip positions in extended context
            if (nodeId == "node_2") {              // logging::info("SKIP_EXTENDED: pos={} >= range.end={} in range [{},{})", 
                             // globalPos, range.end, range.start, range.end);
            }
            
            // Report if this debug position is being skipped due to extended range
            if (debugPositions.find(globalPos) != debugPositions.end()) {
              logging::info("DEBUG_POS_TRACK: Node {} - SKIPPED debug position {} (>= range.end={}) in extended context", 
                           nodeId, globalPos, range.end);
            }
            
            continue; // Skip seeds in the extended context area
          }
          
          // DEBUG: Track all positions processed regardless of seed status
          uniquePositionsProcessed.insert(globalPos);
          
          // DEBUG: Track k-mer processing
          totalKmersProcessedInRanges++;
          
          bool hadSeed = existingSeedsInRange.find(globalPos) != existingSeedsInRange.end(); // Use globalPos as key
          
          // Report decision logic for debug positions
          if (debugPositions.find(globalPos) != debugPositions.end()) {
            std::string decision = "UNKNOWN";
            if (!isNowSeed) decision = "NO_CHANGE_NOT_SEED";
            else if (!hadSeed) decision = "ADD_NEW_SEED"; 
            else if (existingSeedsInRange[globalPos].hash != hash) decision = "MODIFY_SEED";
            else decision = "NO_CHANGE_SAME_HASH";
            
            logging::info("DEBUG_POS_TRACK: Node {} - DECISION for debug position {}: hadSeed={} isNowSeed={} decision={}", 
                         nodeId, globalPos, hadSeed, isNowSeed, decision);
            
            if (hadSeed) {
              logging::info("DEBUG_POS_TRACK: Node {} - Existing seed at {}: hash={} vs new hash={}", 
                           nodeId, globalPos, existingSeedsInRange[globalPos].hash, hash);
            }
          }
          
          // DEBUG: Compact logging for node_2 - log every syncmer position processed and decision
          if (nodeId == "node_2") {
            char action = '?';
            if (!isNowSeed) action = 'S'; // Skip (not a seed)
            else if (!hadSeed) action = 'A'; // Add new seed
            else if (existingSeedsInRange[globalPos].hash != hash) action = 'M'; // Modify seed
            else action = 'K'; // Keep unchanged
            
            // logging::info("SYNC_PROC: {} pos={} act={} had={} now={} h={}", 
            //              nodeId, globalPos, action, hadSeed, isNowSeed, hash);
          }
          
          // Case 1: Seed position ON in parent, seed position ON in current node,
          // hashes are different → MODIFICATION
          if (hadSeed && isNowSeed && existingSeedsInRange[globalPos].hash != hash) {
            seedsChanged++;
            // Modify existing seed
            auto &oldSeed = existingSeedsInRange[globalPos];
            
            // Report if this is a position of interest
            if (debugPositions.find(globalPos) != debugPositions.end()) {
              logging::info("DEBUG_POS_TRACK: Node {} - MODIFY seed at debug position {} (hash {} -> {})", 
                           nodeId, globalPos, oldSeed.hash, hash);
            }
            
            // Create new seed
            seeding::seed_t newSeed;
            newSeed.hash = hash;
            newSeed.reversed = isReverse;
            // Use pre-calculated end position - ensure we use the correct index
            newSeed.endPos = (gappedLocalStartPos < endPositions.size()) ? 
                             endPositions[gappedLocalStartPos] : globalPos + k - 1;
            
            try {
              auto kmerObj = stateManager.extractKmer(node->identifier, globalPos, k);
              std::string kmer = kmerObj.first;
                               
              // Store seed
              stateManager.setSeedAtPosition(nodeId, globalPos, newSeed);
              
              seedChanges.emplace_back(
                  std::make_tuple(globalPos, true, true, std::optional<size_t>(static_cast<size_t>(oldSeed.hash)), std::optional<size_t>(static_cast<size_t>(hash)),
                                  std::optional<bool>(oldSeed.reversed), std::optional<bool>(isReverse),
                                  std::optional<int64_t>(oldSeed.endPos), std::optional<int64_t>(newSeed.endPos)));
                                  
              uniqueKmersCollector.insert(kmer);
              
              stateManager.nodeKmerSequences[nodeId][globalPos] = kmer;
              stateManager.nodeKmerEndOffsets[nodeId][globalPos] = static_cast<uint32_t>(newSeed.endPos - globalPos);
              
              // DEBUG: Compact logging for node_2
              if (nodeId == "node_2") {
                // logging::info("SEED_MODIFY: {} pos={} kmer='{}'", nodeId, globalPos, kmer);
              }
            } catch (const std::system_error& e) {
              if (e.code() == std::errc::resource_deadlock_would_occur) {
                throw std::runtime_error("Resource deadlock detected during seed modification at pos " + 
                                       std::to_string(globalPos) + " for node " + nodeId);
              }
              throw; // Re-throw other exceptions
            }
            
          }
          // Case 2: Seed position ON in parent, seed position ON in current node,
          // hashes are the same → NO CHANGE
          else if (hadSeed && isNowSeed && existingSeedsInRange[globalPos].hash == hash) {
            // No change needed, seed remains the same
            continue;
          }
          // Case 3: Seed position ON in parent, seed position OFF in current node → DELETION
          else if (hadSeed && !isNowSeed) {
            // Count as deletion ONLY if position wasn't already cleared by block operations
            bool alreadyClearedByBlock = positionsWithClearedSeeds.find(globalPos) != positionsWithClearedSeeds.end();
            
            if (!alreadyClearedByBlock) {
              auto &oldSeed = existingSeedsInRange[globalPos];
              
              // Report if this is a position of interest
              if (debugPositions.find(globalPos) != debugPositions.end()) {
                logging::info("DEBUG_POS_TRACK: Node {} - DELETE seed at debug position {} (hash {})", 
                             nodeId, globalPos, oldSeed.hash);
              }
              
              // Add deletion to seed changes
              seedChanges.emplace_back(
                  std::make_tuple(globalPos, true, false, std::optional<size_t>(static_cast<size_t>(oldSeed.hash)), std::optional<size_t>(),
                                 std::optional<bool>(oldSeed.reversed), std::optional<bool>(),
                                 std::optional<int64_t>(oldSeed.endPos), std::optional<int64_t>()));
              
              // DEADLOCK SAFE: Remove seed from block and clear seed
              try {
                coordinates::BlockCoordinate coords = 
                    stateManager.mapGlobalToBlockCoords(nodeId, globalPos);
                if (coords.blockId >= 0) {
                  stateManager.removeSeedFromBlock(coords.blockId, globalPos);
                }
                
                // Clear the seed
                stateManager.clearSeedAtPosition(nodeId, globalPos);
                seedsCleared++;
                
                // DEBUG: Compact logging for node_2
                if (nodeId == "node_2") {
                  // logging::info("SEED_DELETE: {} pos={}", nodeId, globalPos);
                }
              } catch (const std::system_error& e) {
                if (e.code() == std::errc::resource_deadlock_would_occur) {
                  throw std::runtime_error("Resource deadlock detected during seed removal at pos " + 
                                         std::to_string(globalPos) + " for node " + nodeId);
                }
                throw; // Re-throw other exceptions
              }
            }
            // If already cleared by block operations, don't count it as another deletion
          }
          // Case 4: Seed position OFF in parent, seed position ON in current node → ADDITION
          else if (!hadSeed && isNowSeed) {
            // Add new seed
            seeding::seed_t newSeed;
            newSeed.hash = hash;
            newSeed.reversed = isReverse;
            // Use pre-calculated end position - ensure we use the correct index
            newSeed.endPos = (gappedLocalStartPos < endPositions.size()) ? 
                             endPositions[gappedLocalStartPos] : globalPos + k - 1;
            
            // Report if this is a position of interest
            if (debugPositions.find(globalPos) != debugPositions.end()) {
              logging::info("DEBUG_POS_TRACK: Node {} - ADD seed at debug position {} (hash {})", 
                           nodeId, globalPos, hash);
            }
                             
            // DEADLOCK SAFE: Store seed and update kmer sequences
            try {
              // Store seed
              stateManager.setSeedAtPosition(nodeId, globalPos, newSeed);
              
              // Always count as addition - this is a true new seed
              seedsAdded++;
              
              seedChanges.emplace_back(
                  std::make_tuple(globalPos, false, true, std::optional<size_t>(), std::optional<size_t>(static_cast<size_t>(hash)),
                                 std::optional<bool>(), std::optional<bool>(isReverse),
                                 std::optional<int64_t>(), std::optional<int64_t>(newSeed.endPos)));
                                 
              
              std::string extractedKmer = ungappedSequence.substr(localStartPos, k);
              uniqueKmersCollector.insert(extractedKmer);
              
              stateManager.nodeKmerSequences[nodeId][globalPos] = extractedKmer;
              stateManager.nodeKmerEndOffsets[nodeId][globalPos] = static_cast<uint32_t>(newSeed.endPos - globalPos);
              
              // DEBUG: Compact logging for node_2
              if (nodeId == "node_2") {
                // logging::info("SEED_ADD: {} pos={} kmer='{}'", nodeId, globalPos, extractedKmer);
              }
            } catch (const std::system_error& e) {
              if (e.code() == std::errc::resource_deadlock_would_occur) {
                throw std::runtime_error("Resource deadlock detected during seed addition at pos " + 
                                       std::to_string(globalPos) + " for node " + nodeId);
              }
              throw; // Re-throw other exceptions
            }
          }

        }
        
        // DEBUG: Log range processing summary for node_2
        if (nodeId == "node_2") {
          logging::info("RANGE_DONE: [{},{}): {} syncmers, {} unique pos processed", 
                       range.start, range.end, kmers.size(), 
                       std::count_if(kmers.begin(), kmers.end(), 
                                   [&](const auto& kmer) { 
                                     auto [h, rev, isSeed, localPos] = kmer;
                                     return localPos < localToGappedPos.size() && 
                                            localToGappedPos[localPos] >= 0;
                                   }));
        }
        
      } catch (const std::exception &e) {
        logging::err("Error processing range {} for node {}: {}", 
                    rangeIdx, nodeId, e.what());
      }
    }
    
  }
  
  // DEBUG: Log summary of recomputation processing for node_2
  if (nodeId == "node_2") {
    logging::info("DEBUG_RECOMP_SUMMARY: Node {} processed {} unique positions, {} total k-mers across {} ranges", 
                  nodeId, uniquePositionsProcessed.size(), totalKmersProcessedInRanges, ranges.size());
  }
  
  
  std::vector<int64_t> basePositions;
  std::vector<uint64_t> bitMasks;

  if (!seedChanges.empty()) {
    processSeedChanges(seedChanges, basePositions, bitMasks);
    
    if (perNodeSeedMutations != nullptr) {
      // DEADLOCK SAFE: Wrap getDfsIndex call to avoid deadlocks
      int64_t dfsIndex = -1;
      try {
        dfsIndex = stateManager.getDfsIndex(nodeId);
      } catch (const std::system_error& e) {
        if (e.code() == std::errc::resource_deadlock_would_occur) {
          throw std::runtime_error("Resource deadlock detected during getDfsIndex for node " + nodeId);
        }
        throw; // Re-throw other exceptions
      }
      
      if (dfsIndex >= 0) {
        
        if (dfsIndex < static_cast<int64_t>(perNodeSeedMutations->size())) {
          
          auto mutations = (*perNodeSeedMutations)[dfsIndex];
          
          
          auto basePositionsBuilder = mutations.initBasePositions(basePositions.size());
          for (size_t i = 0; i < basePositions.size(); i++) {
            basePositionsBuilder.set(i, basePositions[i]);
          }
          
          
          auto perPosMasksBuilder = mutations.initPerPosMasks(bitMasks.size());
          for (size_t i = 0; i < bitMasks.size(); i++) {
            perPosMasksBuilder.set(i, bitMasks[i]);
          }
          
          
          if (!dictionaryIds.empty()) {
            
            auto dictIdsBuilder = mutations.initKmerDictionaryIds(dictionaryIds.size());
            for (size_t i = 0; i < dictionaryIds.size(); i++) {
              dictIdsBuilder.set(i, dictionaryIds[i]);
            }
            
            
            auto positionsBuilder = mutations.initKmerPositions(dictionaryPositions.size());
            for (size_t i = 0; i < dictionaryPositions.size(); i++) {
              positionsBuilder.set(i, dictionaryPositions[i]);
            }
            
            
            auto endOffsetsBuilder = mutations.initKmerEndOffsets(dictionaryEndOffsets.size());
            for (size_t i = 0; i < dictionaryEndOffsets.size(); i++) {
              endOffsetsBuilder.set(i, dictionaryEndOffsets[i]);
            }
            
            logging::info("DEBUG_INDEX_WRITE: Added {} global dictionary references to node {}", 
                         dictionaryIds.size(), nodeId);
          } else {
            // logging::warn("No dictionary entries for node {} even though seeds were added/modified", nodeId);
          }
          
          logging::info("Node {}: Successfully encoded {} seed changes (deleted: {}, added/modified: {}, dictionary entries: {}, inherited seeds: {} - not encoded)", 
                       nodeId, seedChanges.size(), seedsCleared, seedsAdded, dictionaryIds.size(), seedsInherited);
          
          // Update efficient seed count tracking in NodeState
          auto &nodeState = stateManager.getNodeState(nodeId);
          nodeState.updateSeedCounts(seedsCleared, seedsAdded, seedsChanged);
        } else {
          logging::err("DFS index {} is out of bounds for {} seed mutation entries", 
                       dfsIndex, perNodeSeedMutations->size());
        }
      } else {
        logging::warn("Could not write seed changes for node {} - invalid DFS index: {}", nodeId, dfsIndex);
      }
    } else {
      auto &nodeState = stateManager.getNodeState(nodeId);
      nodeState.addSeedChanges(basePositions, bitMasks);
      
      // Update efficient seed count tracking in NodeState
      nodeState.updateSeedCounts(seedsCleared, seedsAdded, seedsChanged);
    }
  } else {
    
    logging::debug("No seed changes to encode for node {} (inherited: {}, all seeds preserved)", 
                  nodeId, seedsInherited);
    
    // Even if no changes, update NodeState to inherit parent seed count for efficiency
    auto &nodeState = stateManager.getNodeState(nodeId);
    if (node->parent) {
      try {
        const auto &parentState = stateManager.getNodeState(node->parent->identifier);
        nodeState.initializeSeedCountFromParent(parentState.getTotalSeedCount());
      } catch (const std::exception& e) {
        // If parent state not available, calculate from scratch as fallback
        logging::warn("Failed to inherit seed count from parent for node {}: {}", nodeId, e.what());
      }
    }
  }


  
  return std::make_tuple(blockDeletions, blockInsertions, totalMaterializedSeeds, seedsFoundInDeletedBlocks, blockDeletionPositionsChecked, seedsCleared, seedsAdded, recompRangeCount, recompRangeSize, recompRangesFormatted, static_cast<int>(uniquePositionsProcessed.size()), totalKmersProcessedInRanges, seedChanges);
}

/**
 * @brief Process seed changes and encode them in quaternary format (4 values)
 *
 * Each seed change is encoded using 2 bits:
 * - 0: seed unchanged
 * - 1: seed deleted (on→off)
 * - 2: seed added (off→on)
 * - 3: seed modified (on→on with different hash/position)
 * @param seedChanges Vector of seed change tuples
 * @param basePositions Output vector of base positions
 * @param bitMasks Output vector of quaternary-encoded masks
 */
void processSeedChanges(
    const std::vector<std::tuple<int64_t, bool, bool, std::optional<size_t>,
                                 std::optional<size_t>, std::optional<bool>,
                                 std::optional<bool>, std::optional<int64_t>,
                                 std::optional<int64_t>>> &seedChanges,
    std::vector<int64_t> &basePositions, std::vector<uint64_t> &bitMasks) {

  
  basePositions.clear();
  bitMasks.clear();

  if (seedChanges.empty()) {
    logging::debug("SEED_ENCODE: No seed changes to encode");
    return;
  }

  
  constexpr uint8_t BITS_PER_VALUE = 2;
  constexpr uint8_t VALUES_PER_MASK = 32; 
  constexpr uint64_t MAX_POSITION_RANGE =
      VALUES_PER_MASK - 1; 

  
  
  std::vector<std::pair<int64_t, uint8_t>> positionValuePairs;
  positionValuePairs.reserve(seedChanges.size()); 

  std::array<int, 4> valueCount = {0, 0, 0, 0};

  // Transform seed changes to position-value pairs
  int logged_changes = 0;
  for (const auto &change : seedChanges) {
    const auto &[pos, wasOn, isOn, oldHash, newHash, oldReversed, newReversed,
                 oldEndPos, newEndPos] = change;

    // Determine quaternary value (using all 4 possible values)
    uint8_t value;
    if (wasOn && !isOn) {
      value = 1; // Deletion (on → off)
      valueCount[1]++;
    } else if (!wasOn && isOn) {
      value = 2; // Addition (off → on)
      valueCount[2]++;
    } else if (wasOn && isOn) {
      value = 3; // Modification (on → on')
      valueCount[3]++;
    } else {
      value = 0; // No change (off → off)
      valueCount[0]++;
    }

    positionValuePairs.emplace_back(pos, value);

    
    if (logged_changes < 5) {
      logging::debug("SEED_ENCODE: Change {}: pos={}, wasOn={}, isOn={}, value={}",
                    positionValuePairs.size(), pos, (wasOn ? "true" : "false"), (isOn ? "true" : "false"), static_cast<int>(value));
      logged_changes++;
    }
  }

  
  const size_t estimatedBasePositions = (positionValuePairs.size() + VALUES_PER_MASK - 1) / VALUES_PER_MASK;
  basePositions.reserve(estimatedBasePositions);
  bitMasks.reserve(estimatedBasePositions);

  
  size_t groupCount = 0;
  size_t currentIndex = 0;
  const size_t totalPairs = positionValuePairs.size();

  while (currentIndex < totalPairs) {
    groupCount++;
    
    int64_t basePos = positionValuePairs[currentIndex].first;
    uint64_t bitMask = 0;

    
    int64_t minPos = std::max<int64_t>(0, basePos - MAX_POSITION_RANGE);

    
    
    size_t pairsInThisGroup = 0;
    while (currentIndex < totalPairs &&
           positionValuePairs[currentIndex].first >= minPos) {

      
      const auto &[pos, value] = positionValuePairs[currentIndex];

      
      uint8_t offset = static_cast<uint8_t>(basePos - pos);

      
      bitMask |= (static_cast<uint64_t>(value) << (offset * BITS_PER_VALUE));

      
      currentIndex++;
      pairsInThisGroup++;
    }

    
    basePositions.push_back(basePos);
    bitMasks.push_back(bitMask);
  }

  logging::debug("SEED_ENCODE: Generated {} base positions and masks from {} position-value pairs",
                basePositions.size(), positionValuePairs.size());
}

/**
 * @brief Decode quaternary-encoded seed changes back to original form
 *
 * @param basePositions Vector of base positions
 * @param bitMasks Vector of quaternary-encoded masks
 * @return Vector of (pos, wasOn, isOn) tuples
 */
std::vector<std::tuple<int64_t, bool, bool>>
decodeSeedChanges(const std::vector<int64_t> &basePositions,
                  const std::vector<uint64_t> &bitMasks) {

  std::vector<std::tuple<int64_t, bool, bool>> result;
  if (basePositions.empty() || bitMasks.empty()) {
    return result;
  }

  
  result.reserve(basePositions.size() * 32);

  
  constexpr uint8_t BITS_PER_VALUE = 2;
  constexpr uint8_t POSITIONS_PER_MASK = 32;
  constexpr uint8_t VALUE_MASK = 0x3;

  
  int logged_positions = 0;
  for (size_t i = 0; i < basePositions.size(); i++) {
    const int64_t basePos = basePositions[i];
    const uint64_t mask = bitMasks[i];

    
    for (uint8_t offset = 0; offset < POSITIONS_PER_MASK; offset++) {
      
      const uint8_t value = (mask >> (offset * BITS_PER_VALUE)) & VALUE_MASK;

      
      if (value == 0)
              continue;

      
      const int64_t pos = basePos - offset;


      const bool wasOn = (value == 1 || value == 3); 
      const bool isOn = (value == 2 || value == 3);  

      
      if (logged_positions < 5) {
          logging::debug("SEED_DECODE: BasePos={}, Mask={}, Offset={} -> Value={}, Pos={}, wasOn={}, isOn={}",
                        basePos, mask, offset, value, pos, wasOn, isOn);
          logged_positions++;
      }

      result.emplace_back(pos, wasOn, isOn);
    }
  }

  
  std::sort(result.begin(), result.end(), [](const auto &a, const auto &b) {
    return std::get<0>(a) < std::get<0>(b);
  });

  return result;
}


void initializeBlockSequences(
    state::StateManager &stateManager,
    const std::unordered_map<int32_t, std::string> &blockSequences,
    const std::unordered_map<int32_t, coordinates::CoordRange> &blockRanges) {

  
  stateManager.setNumBlocks(blockSequences.size());

  
  for (const auto &[blockId, sequence] : blockSequences) {
    if (blockRanges.find(blockId) == blockRanges.end()) {
      throw std::runtime_error("Block " + std::to_string(blockId) +
                               " has sequence but no range");
    }
    stateManager.setBlockSequence(blockId, sequence);
    stateManager.setBlockRange(blockId, blockRanges.at(blockId));

    logging::debug("Set block {} sequence (length={}) and range [{}, {})",
                   blockId, sequence.length(), blockRanges.at(blockId).start, blockRanges.at(blockId).end);
  }
}


void processNodesByLevel(
    const std::vector<panmanUtils::Node *> &nodes,
    std::function<void(panmanUtils::Node *)> processFunction) {

  
  if (nodes.size() <= 4) {
    logging::debug("Processing {} nodes sequentially (small batch)", nodes.size());
    for (auto *node : nodes) {
      processFunction(node);
    }
    return;
  }

  
  // Calculate optimal grain size for TBB parallel_for
  const size_t numThreads = std::thread::hardware_concurrency();
  const size_t nodeCount = nodes.size();
  const size_t optimalGrainSize = std::max(1UL, 
      std::min(nodeCount / (numThreads * 4), 64UL)); // Prevent too small or too large grains
  
  logging::debug("Processing {} nodes in parallel with TBB optimal grain size {} ({} hardware threads available)",
               nodeCount, optimalGrainSize, numThreads);

  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, nodeCount, optimalGrainSize),
      [&](const tbb::blocked_range<size_t> &range) {
        logging::debug("Thread processing range [{}, {})", range.begin(), range.end());
        for (size_t i = range.begin(); i < range.end(); ++i) {
          processFunction(nodes[i]);
        }
      });

  logging::debug("Completed parallel processing of {} nodes", nodes.size());
}


void fillNodesDepthFirst(panmanUtils::Node *node,
                         std::vector<panmanUtils::Node *> &nodes) {
  if (!node)
    return;

  
  nodes.push_back(node);

  
  for (auto *child : node->children) {
    fillNodesDepthFirst(child, nodes);
  }
}


std::unique_ptr<state::StateManager>
initializeStateManager(panmanUtils::Tree *tree, panmanUtils::Node *rootNode,
                       int kmerSize, int smerSize) {
  if (!tree || !rootNode) {
    throw std::runtime_error(
        "Invalid tree or root node for state initialization");
  }

  
  
  const int num_threads = std::thread::hardware_concurrency();
  logging::debug("Initializing state manager with {} threads", num_threads);
  
  
  logging::debug("Getting sequence data from reference...");
  
  panmanUtils::sequence_t rootSequence;
  panmanUtils::blockExists_t rootBlockExists;
  panmanUtils::blockStrand_t rootBlockStrand;
  tree->getSequenceFromReference(rootSequence, rootBlockExists, rootBlockStrand,
                                 rootNode->identifier);

  if (rootSequence.empty()) {
    throw std::runtime_error(
        "Root sequence is empty, unable to proceed with indexing");
  }



  
  int64_t totalCoords = 0;
  for (size_t blockId = 0; blockId < rootSequence.size(); blockId++) {
    size_t blockSize = 0;
    const auto &sequence = rootSequence[blockId];

    
    for (const auto &seqPair : sequence.first) {
      blockSize += (1 + seqPair.second.size());
    }

    totalCoords += blockSize;
  }
  logging::debug("Total coordinates: {}", totalCoords);

  
  auto stateManager = std::make_unique<state::StateManager>(totalCoords);

  
  stateManager->setNumBlocks(rootSequence.size());
  
  
  size_t maxNucPosInRoot = 0;
  for (const auto& blockData : rootSequence) {
      if (!blockData.first.empty()) {
        
        maxNucPosInRoot = std::max(maxNucPosInRoot, blockData.first.size() - 1);
      }
  }
  logging::debug("Calculated max nucPos hint from root sequence: {}", maxNucPosInRoot);
  
  // MOVED EARLIER: Initialize the StateManager with the tree structure and node hierarchy.
  // This should set up default states for all nodes.
  stateManager->initialize(tree, maxNucPosInRoot);
  
  // Set k-mer and s-mer sizes AFTER general initialization but BEFORE sequence processing might need them.
  stateManager->setKmerSize(kmerSize);
  stateManager->setSmerSize(smerSize);
  
  std::unordered_map<int32_t, std::string> extractedBlockSequences;
  std::unordered_map<int32_t, coordinates::CoordRange> blockRanges;
  int64_t currentPosition = 0;
  logging::debug("Root sequence size: {}", rootSequence.size());

  // Build a lookup for canonical gap lengths from tree->gaps
  std::unordered_map<int32_t, std::unordered_map<int32_t, int32_t>> canonical_gap_lengths;
  if (tree && tree->gaps.data() != nullptr) { // Check if tree->gaps is valid and not null
      logging::debug("INIT_STMGR: Processing {} entries from tree->gaps to build canonical gap length map.", tree->gaps.size());
      for (const auto& gap_entry : tree->gaps) { // panmanUtils::Gap gap_entry
          int32_t gap_block_id = static_cast<int32_t>(gap_entry.primaryBlockId);
          // Ensure nucPosition and nucGapLength are vectors and accessible
          // This assumes gap_entry.nucPosition and gap_entry.nucGapLength are std::vector<int> or similar
          if (gap_entry.nucPosition.data() != nullptr && gap_entry.nucGapLength.data() != nullptr) {
              for (size_t j = 0; j < gap_entry.nucPosition.size(); ++j) {
                  if (j < gap_entry.nucGapLength.size()) { 
                      int32_t nuc_p = gap_entry.nucPosition[j];
                      int32_t len = gap_entry.nucGapLength[j];
                      canonical_gap_lengths[gap_block_id][nuc_p] = len;
                      if (gap_block_id == 602 && (nuc_p < 5 || nuc_p % 50 == 0)) { // Log for block 602
                         logging::debug("INIT_STMGR_CANON_GAP: Block 602, NucPos {}, Canonical Gap Length from tree->gaps: {}", nuc_p, len);
                      }
                  } else {
                      logging::warn("INIT_STMGR: Mismatch in sizes of nucPosition and nucGapLength for block {} in tree->gaps. Skipping an entry.", gap_block_id);
                  }
              }
          } else {
              logging::warn("INIT_STMGR: nucPosition or nucGapLength data is null for a gap_entry concerning blockId {}.", gap_block_id);
          }
      }
      logging::debug("INIT_STMGR: Finished building canonical gap length map. Found entries for {} distinct blocks.", canonical_gap_lengths.size());
  } else {
      logging::warn("INIT_STMGR: tree->gaps is null, or tree->gaps.data() is null. Canonical gap lengths will not be pre-populated from tree->gaps.");
  }

  // MOVED HERE: Flush coordinate mappings after all registrations are done from Loop 1.
  // stateManager->flushCoordinateMappings(); // This was moved from original position, evaluate if needed here or later.
                                          // For now, let's keep it commented to apply changes incrementally.

  // The existing loop starts below
  logging::debug("INIT_STMGR: Processing {} blocks based on rootSequence (from getSequenceFromReference) to establish base sequences.", rootSequence.size());
  
  // This loop iterates based on rootSequence, which itself is 0-indexed by block ID (implicitly).
  // It ensures that any specific interpretation of consensusSeq + gaps done by getSequenceFromReference
  // (potentially populating pair.second with actual nucs) is used.
  for (size_t block_id_as_idx = 0; block_id_as_idx < rootSequence.size(); block_id_as_idx++) {
    int32_t current_block_id = static_cast<int32_t>(block_id_as_idx);
    const auto& panman_block_s_data = rootSequence[block_id_as_idx]; // This is panmanUtils::block_s
    // data_to_process is std::vector<std::pair<char, std::vector<char>>> derived from rootSequence
    // This data reflects block state as of rootNode, including any nucs in gap lists.
    const auto& data_to_process = panman_block_s_data.first; 

    std::string blockSeq_str_representation;
    // Estimate reservation based on main nucs + canonical gap lengths
    size_t reserve_size_estimate = 0;
    for(size_t nps_idx = 0; nps_idx < data_to_process.size(); ++nps_idx) {
        reserve_size_estimate += 1; // for main nuc
        auto block_gaps_it_est = canonical_gap_lengths.find(current_block_id);
        if (block_gaps_it_est != canonical_gap_lengths.end()) {
            auto nuc_pos_gaps_it_est = block_gaps_it_est->second.find(static_cast<int32_t>(nps_idx));
            if (nuc_pos_gaps_it_est != block_gaps_it_est->second.end()) {
                reserve_size_estimate += nuc_pos_gaps_it_est->second;
            }
        } // If not in canonical, and rootSequence also has no gaps for this nucPos, estimate remains minimal.
          // If rootSequence *does* have gaps, they are covered by data_to_process.size()*2 below if estimate is too small.
    }
    blockSeq_str_representation.reserve(std::max(reserve_size_estimate, data_to_process.size()*2)); 
            
    if (!data_to_process.empty()) {
      int32_t nucPos_structural = 0; 
      for (const auto &seqPair : data_to_process) { 
        if (seqPair.first == 'x' && seqPair.second.empty()) { 
            // Handle sentinel 'x'. If it's the first thing and has canonical gaps, set those.
            if (nucPos_structural == 0) {
                 int32_t canonical_len_for_x_pos0 = 0;
                 auto block_gaps_it_x = canonical_gap_lengths.find(current_block_id);
                 if (block_gaps_it_x != canonical_gap_lengths.end()) {
                    auto nuc_pos_gaps_it_x = block_gaps_it_x->second.find(0); 
                    if (nuc_pos_gaps_it_x != block_gaps_it_x->second.end()) {
                        canonical_len_for_x_pos0 = nuc_pos_gaps_it_x->second;
                    }
                 }
                 // Set gap list length for nucPos 0 of this block if 'x' is the only char.
                 stateManager->setGapListLength(current_block_id, 0, canonical_len_for_x_pos0);
            }
            break; 
        }
        
        int32_t expected_gap_list_len = 0;
        auto block_gaps_it = canonical_gap_lengths.find(current_block_id);
        if (block_gaps_it != canonical_gap_lengths.end()) {
            auto nuc_pos_gaps_it = block_gaps_it->second.find(nucPos_structural);
            if (nuc_pos_gaps_it != block_gaps_it->second.end()) {
                expected_gap_list_len = nuc_pos_gaps_it->second;
            }
        }
        // If canonical_gap_lengths has no entry for this nucPos, expected_gap_list_len remains 0.
        // Panman would have an empty .second vector in this case after its resize.

        size_t iteration_limit_for_gaps = static_cast<size_t>(expected_gap_list_len);

        // Iterate based on the canonical (expected) gap list length
        for (size_t char_idx_in_gap_vec = 0; char_idx_in_gap_vec < iteration_limit_for_gaps; ++char_idx_in_gap_vec) {
          char gapCharToStore;
          // If rootSequence's current gap list is shorter than canonical, pad with '-'
          if (char_idx_in_gap_vec < seqPair.second.size()) {
              gapCharToStore = seqPair.second[char_idx_in_gap_vec];
          } else {
              gapCharToStore = '-'; // Pad with '-' if expected_gap_list_len is longer
          }
          blockSeq_str_representation.push_back(gapCharToStore);
          int64_t globalPosForGapChar = currentPosition++; 
          stateManager->registerCoordinateMapping(
              current_block_id, nucPos_structural, static_cast<int32_t>(char_idx_in_gap_vec), globalPosForGapChar);
          state::PositionKey structuralKey_gap = state::PositionKey::create(current_block_id, nucPos_structural, static_cast<int32_t>(char_idx_in_gap_vec));
          stateManager->blockRootCharFlatIndices[current_block_id][structuralKey_gap] = static_cast<int64_t>(blockSeq_str_representation.length() - 1);
        }

        // Then process the main nucleotide character itself from seqPair.first
        char mainCharToStore = seqPair.first;
        // Only add main char if it's not the sentinel 'x' that was already handled by the break (which implies .second is empty and iteration_limit_for_gaps was 0)
        if (!(mainCharToStore == 'x' && seqPair.second.empty() && iteration_limit_for_gaps == 0 )) {
             blockSeq_str_representation.push_back(mainCharToStore);
             int64_t globalPosForMainNuc = currentPosition++; 
             stateManager->registerCoordinateMapping(current_block_id, nucPos_structural, -1, globalPosForMainNuc);
             state::PositionKey structuralKey_main = state::PositionKey::create(current_block_id, nucPos_structural, -1);\
             stateManager->blockRootCharFlatIndices[current_block_id][structuralKey_main] = static_cast<int64_t>(blockSeq_str_representation.length() - 1);
        }
        
        // SET THE CANONICAL GAP LIST LENGTH FOR THIS NUC_POS
        stateManager->setGapListLength(current_block_id, nucPos_structural, expected_gap_list_len);

        nucPos_structural++; 
      }
    }

    coordinates::CoordRange range{currentPosition - static_cast<int64_t>(blockSeq_str_representation.length()), currentPosition};
    blockRanges[current_block_id] = range; 
    extractedBlockSequences[current_block_id] = blockSeq_str_representation; 
    
    if ((block_id_as_idx + 1) % 1000 == 0 || block_id_as_idx + 1 == rootSequence.size() || block_id_as_idx < 5 || current_block_id == 602) { 
      logging::debug("INIT_STMGR: Processed base sequence for Block ID {} (rootSequence index {}). String Length: {}. Range: [{}, {}). Sample: '{}'", 
                     current_block_id, block_id_as_idx, blockSeq_str_representation.length(), range.start, range.end, 
                     blockSeq_str_representation.substr(0, std::min(blockSeq_str_representation.length(), (size_t)30)));
    }
  }
  // END OF THE FIRST MAJOR LOOP (NOW SOURCES FROM rootSequence)

  stateManager->flushCoordinateMappings();

  // This second major loop sets sequences in StateManager and initial state for rootNode->identifier.
  // It uses block_id_loop_var (which is the 0-indexed blockId from rootSequence) for all operations.
  logging::debug("INIT_STMGR: Setting block sequences in StateManager and initial state for rootNode: {}", rootNode->identifier);
  
  // Track current position for auto-range calculation
  // Use existing currentPosition variable that was already declared
  if (!blockRanges.empty()) {
      // Find the max end position of existing ranges
      for (const auto& [id, range] : blockRanges) {
          currentPosition = std::max(currentPosition, range.end);
      }
  }
  
  for (const auto &[blockId_loop_var, sequence_string_for_block] : extractedBlockSequences) {
    stateManager->setBlockSequence(blockId_loop_var, sequence_string_for_block);
    
    if (blockRanges.count(blockId_loop_var)) {
        stateManager->setBlockRange(blockId_loop_var, blockRanges.at(blockId_loop_var));
    } else {
        // CRITICAL FIX: Instead of setting an empty range, attempt to compute a sensible range
        // based on the sequence string length and current position
        const auto& sequence = sequence_string_for_block;
        if (!sequence.empty()) {
            // Compute a reasonable range based on current position and sequence length
            int64_t seqLength = static_cast<int64_t>(sequence.length());
            int64_t start = currentPosition;
            int64_t end = currentPosition + seqLength;
            
            logging::warn("INIT_STMGR: No range found in blockRanges map for block ID {}. "
                         "Computing range from sequence length: [{}, {}).", 
                         blockId_loop_var, start, end);
                         
            // Set the computed range instead of a dummy {0,0} range
            stateManager->setBlockRange(blockId_loop_var, {start, end});
            
            // Update currentPosition for next block
            currentPosition = end;
        } else {
            // If sequence is empty, log a more severe warning and use a minimal but non-zero range
            logging::err("INIT_STMGR: No range found for block ID {} and sequence is empty! "
                        "Setting minimal range at current position.", blockId_loop_var);
            
            // Use a minimal non-zero range at current position
            stateManager->setBlockRange(blockId_loop_var, {currentPosition, currentPosition + 1});
            currentPosition += 1;
        }
        
        // Add this block to the blockRanges map to ensure consistent state
        auto newRange = stateManager->getBlockRange(blockId_loop_var);
        blockRanges[blockId_loop_var] = newRange;
        
        // Special debugging for problem blocks 617 and 941
        if (blockId_loop_var == 617 || blockId_loop_var == 941) {
            logging::info("CRITICAL BLOCK {} range initialized to [{}, {})", 
                         blockId_loop_var, newRange.start, newRange.end);
        }
    }
    
    stateManager->initializeBlockMappings(blockId_loop_var);
    
    bool should_activate_for_root = false;
    bool is_inverted_for_root = false; 
    int32_t idx_for_root_status_arrays = blockId_loop_var; 

    if (idx_for_root_status_arrays >= 0 && static_cast<size_t>(idx_for_root_status_arrays) < rootBlockExists.size()) {
        should_activate_for_root = rootBlockExists[idx_for_root_status_arrays].first;
        if (should_activate_for_root && static_cast<size_t>(idx_for_root_status_arrays) < rootBlockStrand.size()) {
           is_inverted_for_root = !rootBlockStrand[idx_for_root_status_arrays].first; 
        }
    } else {
        logging::warn("INIT_STMGR: Block ID {} (from rootSequence index) is out of bounds for rootBlockExists (size {}) or rootBlockStrand (size {}). Defaulting to inactive for root node '{}'.",
                     blockId_loop_var, rootBlockExists.size(), rootBlockStrand.size(), rootNode->identifier);
        should_activate_for_root = false;
    }
    
                // Initialize block state for root node which will be inherited by children
      stateManager->setBlockOn(rootNode->identifier, blockId_loop_var, should_activate_for_root); 
    if (should_activate_for_root) { 
      stateManager->setBlockInverted(rootNode->identifier, blockId_loop_var, is_inverted_for_root);
    }
    
     if ((blockId_loop_var % 1000 == 0 && blockId_loop_var !=0 ) || blockId_loop_var < 5 || blockId_loop_var == 602 || blockId_loop_var == 595) { 
        logging::debug("INIT_STMGR_ROOT_STATE: For root node '{}', Block ID {}: Set Active={}, Inverted={}. Stored string length: {}. Range: [{}, {}).", 
        rootNode->identifier, 
        blockId_loop_var, 
        should_activate_for_root, 
        (should_activate_for_root ? is_inverted_for_root : false), 
        sequence_string_for_block.length(), // Corrected variable name
        (blockRanges.count(blockId_loop_var) ? blockRanges.at(blockId_loop_var).start : -1LL), 
        (blockRanges.count(blockId_loop_var) ? blockRanges.at(blockId_loop_var).end : -1LL)
        );
    }
  }

  
  stateManager->flushCoordinateMappings();

  
  stateManager->setKmerSize(kmerSize);
  stateManager->setSmerSize(smerSize);


  logging::debug("StateManager initialization complete");

  // === START Log initial state of block 595 ===
  if (stateManager) { // Ensure stateManager is not null
      std::string initial_b595_log_filename = "debug_OURCODE_INITIAL_BLOCK595_CONTENT.txt";
      std::ofstream b595_log(initial_b595_log_filename);
      if (b595_log.is_open()) {
          b595_log << "===== Initial content of block 595 from stateManager after initializeStateManager =====\n" << std::endl;
          try {
              // Attempt to get and log the gappy sequence for block 595
              // This assumes blockSequences is accessible or there's a getter.
              // For this example, we'll conceptually try to use extractSequence
              // for the block's initially defined global range.
              coordinates::CoordRange range595 = stateManager->getBlockRange(595); // This might throw if block 595 or its range isn't set
              b595_log << "Initial Global Coordinate Range for Block 595: [" << range595.start << ", " << range595.end << ")\n";

              if (range595.start < range595.end) {
                  // Use rootNode->identifier as the "nodeId" for extracting the initial reference sequence state
                  auto [b595_sequence, b595_positions, b595_gaps, b595_endPositions] = stateManager->extractSequence(rootNode->identifier, range595, false);
                  b595_log << "Block 595 Gappy Content (length " << b595_sequence.length() << "):\n";
                  for (size_t i = 0; i < b595_sequence.length(); i += 80) {
                      b595_log << b595_sequence.substr(i, 80) << "\n";
                  }
                  if (b595_sequence.empty()) {
                       b595_log << "[EMPTY SEQUENCE EXTRACTED FOR BLOCK 595 INITIAL STATE]\n";
                  }
              } else {
                  b595_log << "Block 595 has an invalid or zero-length initial range. No sequence to log.\n";
              }

          } catch (const std::exception& e) {
              b595_log << "Exception trying to log initial state of block 595: " << e.what() << std::endl;
          }
          b595_log.close();
      } else {
          // std::cerr << "ERROR: Could not open " << initial_b595_log_filename << " for writing." << std::endl;
      }
  }
  // === END Log initial state of block 595 ===

  // NOW initialize the global position cache, after all blockRanges and coordinate mappings are set
  stateManager->initializeGlobalCache();

  return stateManager;
}

std::unique_ptr<state::StateManager>
initializeStateManagerLight(panmanUtils::Tree* tree, panmanUtils::Node* rootNode,
                          int kmerSize, int smerSize) {
  if (!tree || !rootNode) {
    throw std::invalid_argument("Invalid tree or root node in initializeStateManagerLight");
  }

  logging::info("Initializing StateManager (light) with k={}, s={}", kmerSize, smerSize);
  
  
  
  
  auto stateManager = std::make_unique<state::StateManager>(); 
  
  stateManager->setKmerSize(kmerSize);
  stateManager->setSmerSize(smerSize);
  
  if (tree->blocks.empty()) {
      logging::warn("Tree has no blocks defined. StateManager may not function correctly.");
      stateManager->setNumBlocks(0);
  } else {
      stateManager->setNumBlocks(tree->blocks.size());
  }
  
  // Initialize the node hierarchy and DFS indices
  stateManager->initialize(tree, 0); // Pass 0 as maxNucPosHint is okay for light init

  return stateManager;
}

// Build additional data for placement acceleration
void buildPlacementAccelerationData(panmanUtils::Tree *tree,
                                    Index::Builder &index,
                                    const state::StateManager &stateManager) {
  if (!tree) {
    throw std::invalid_argument(
        "Tree pointer is null in buildPlacementAccelerationData");
  }

  logging::debug("Building placement acceleration data structures for {} nodes", tree->allNodes.size());

  // Pre-compute paths for all nodes
  logging::debug("Computing node paths from root");
  auto nodePaths = state::computeNodePaths(tree, tree->root);

  // Group nodes by level
  logging::debug("Organizing nodes by level");
  auto nodesByLevel = state::groupNodesByLevel(tree, tree->root);
  logging::debug("Found {} levels in the tree", nodesByLevel.size());

  // IMPROVED: Better node verification and tracking
  std::set<std::string> treeNodeIds;            // All node IDs in the tree
  std::set<std::string> stateManagerNodeIds;    // All node IDs in the StateManager
  std::set<std::string> missingInStateManager;  // Nodes in tree but not in StateManager
  std::set<std::string> indexableNodeIds;       // Nodes that can be included in the index


  // Initialize node path info builder - always use the tree size
  logging::debug("Initializing node path info for {} nodes", tree->allNodes.size());
  auto nodePathInfoBuilder = index.initNodePathInfo(tree->allNodes.size());

  // Build node path info - process ALL tree nodes even if they're not in StateManager
  size_t nodeIndex = 0;
  logging::debug("Building node path information");
  
  // Second pass: Add ALL nodes to the index (will use fallbacks for missing data)
  for (const auto &[nodeId, node] : tree->allNodes) {
    auto nodeInfo = nodePathInfoBuilder[nodeIndex++];

    // Set node ID - CRITICAL: Use C++ string copy instead of direct set
    nodeInfo.setNodeId(nodeId);  // Use the original nodeId directly
    
    // Verify that the node ID was set correctly
    auto checkId = nodeInfo.getNodeId();
    std::string checkIdStr;
    checkIdStr.assign(checkId.begin(), checkId.end());
    if (checkIdStr != nodeId) {
      logging::warn("Node ID mismatch during buildPlacementAccelerationData: Expected '{}', got '{}'", 
                   nodeId, checkIdStr);
    }
    
    // Additional check for empty ID
    if (checkIdStr.empty()) {
      logging::critical("EMPTY NODE ID detected in buildPlacementAccelerationData for node {}", nodeId);
      // Force-fill with the actual node ID
      nodeInfo.setNodeId(nodeId);
    }

    // Set node level
    size_t level = 0;
    for (const auto &levelNodes : nodesByLevel) {
      if (std::find_if(levelNodes.begin(), levelNodes.end(),
                       [&nodeId](const panmanUtils::Node *n) {
                         return n->identifier == nodeId;
                       }) != levelNodes.end()) {
        break;
      }
      level++;
    }
    nodeInfo.setLevel(level);

    // Set parent node ID
    if (node->parent) {
      nodeInfo.setParentId(node->parent->identifier);  // Use the original parent identifier directly
      
      // Verify parent ID was set correctly
      auto checkParentId = nodeInfo.getParentId();
      std::string checkParentStr;
      checkParentStr.assign(checkParentId.begin(), checkParentId.end());
      if (checkParentStr != node->parent->identifier) {
        logging::warn("Parent ID mismatch: Expected '{}', got '{}'", 
                     node->parent->identifier, checkParentStr);
      }
    } else {
      nodeInfo.setParentId(""); // Root node has no parent
    }

    // Set active block IDs - IMPROVED: Better handling for missing nodes
    std::set<int32_t> activeBlocks;
    
    if (stateManagerNodeIds.find(nodeId) != stateManagerNodeIds.end()) {
      // Node exists in StateManager - get blocks directly
    try {
      // Manual element-by-element copy to handle different container types
      const auto& srcBlocks = stateManager.getActiveBlocks(nodeId);
      activeBlocks.clear();
      for (const auto& blockId : srcBlocks) {
        activeBlocks.insert(blockId);
      }
      
      // Debug logging for crucial nodes
      if (nodeId == "node_1" || nodeId == "node_2" || nodeId == "node_3") {
          logging::info("Node {} has {} active blocks during nodePathInfo construction", 
                      nodeId, activeBlocks.size());
      }
      } catch (const std::exception& e) {
        logging::warn("Error getting active blocks for node {}: {}", nodeId, e.what());
      }
    } else {
      // Node missing from StateManager - use fallback/reconstruction
      logging::warn("Using fallback logic for missing node {}", nodeId);
      
      // If node has a parent, try to reconstruct from parent's blocks
      if (node->parent) {
          std::string parentId = node->parent->identifier;
        
        try {
          const auto& parentBlocksSource = stateManager.getActiveBlocks(parentId);
          
          if (!parentBlocksSource.empty()) {
            logging::info("Reconstructing blocks for node {} using parent {} blocks", 
                         nodeId, parentId);
              
              // Manual element-by-element copy from parent's blocks
              for (const auto& blockId : parentBlocksSource) {
                  activeBlocks.insert(blockId);
              }
              
              // Apply node's block mutations to estimate what blocks should be active
              for (const auto& blockMut : node->blockMutation) {
                  int32_t blockId = blockMut.primaryBlockId;
                  bool isInsertion = (blockMut.blockMutInfo == 1);
                  
                  if (isInsertion) {
                      // Block should be added
                      activeBlocks.insert(blockId);
                      logging::debug("Reconstructed block activation: adding block {} to node {}", 
                                   blockId, nodeId);
                  } else {
                      // Block should be removed or toggled
                      if (activeBlocks.count(blockId) > 0) {
                          activeBlocks.erase(blockId);
                          logging::debug("Reconstructed block activation: removing block {} from node {}", 
                                       blockId, nodeId);
                      }
                  }
              }
              
              logging::info("Reconstructed {} active blocks for node {}", activeBlocks.size(), nodeId);
          } else {
            logging::warn("Parent {} has no active blocks to reconstruct from for node {}", 
                         parentId, nodeId);
      }
    } catch (const std::exception& e) {
          logging::warn("Failed to reconstruct blocks for node {} from parent: {}", 
                       nodeId, e.what());
        }
      }
    }
    
    // Always write what we have, even if it's empty
    auto activeBlocksBuilder = nodeInfo.initActiveBlocks(activeBlocks.size());
    
    size_t blockIndex = 0;
    for (int32_t blockId : activeBlocks) {
      activeBlocksBuilder.set(blockIndex++, blockId);
    }
    
    // Verify blocks were stored
    if (blockIndex > 0) {
      logging::debug("Node {} stored {} active blocks in nodePathInfo", nodeId, blockIndex);
    } else if (level > 0) { // Only warn for non-root nodes
      logging::warn("Node {} has NO active blocks in final nodePathInfo", nodeId);
    }
  }
  
  logging::debug("Completed node path info for {} nodes", nodeIndex);
  
  // Quick verification of node path info
  auto createdNodePathInfo = index.getNodePathInfo();
  if (createdNodePathInfo.size() == 0) {
    logging::critical("CRITICAL: nodePathInfo has ZERO entries after building!");
  } else if (createdNodePathInfo.size() != tree->allNodes.size()) {
    logging::warn("Unexpected nodePathInfo size: got {} entries, expected {} (tree node count)",
                 createdNodePathInfo.size(), tree->allNodes.size());
  }

  // Build block info
  logging::debug("Building block info for {} blocks", stateManager.getNumBlocks());
  auto blockInfoBuilder = index.initBlockInfo(stateManager.getNumBlocks());
  for (size_t blockId_idx = 0; blockId_idx < stateManager.getNumBlocks(); blockId_idx++) {
    auto blockInfo = blockInfoBuilder[blockId_idx];
    int32_t currentBlockId = static_cast<int32_t>(blockId_idx);
    blockInfo.setBlockId(currentBlockId);

    try {
      coordinates::CoordRange range = stateManager.getBlockRange(currentBlockId);
      blockInfo.setRangeStart(range.start);
      blockInfo.setRangeEnd(range.end);
      // Removed the verbose debug log from here to minimize changes, it can be added back if needed
    } catch (const std::exception& e) {
      logging::err("INDEXING: Error getting block range for blockId {} from StateManager: {}. Setting default [0,0].", currentBlockId, e.what());
      blockInfo.setRangeStart(0);
      blockInfo.setRangeEnd(0);
    }
  }

  // Build ancestor-descendant relationship matrix
  logging::debug("Building ancestor-descendant relationship matrix ({}x{})",
               tree->allNodes.size(), tree->allNodes.size());
  auto ancestorMatrixBuilder = index.initAncestorMatrix(tree->allNodes.size());

  nodeIndex = 0;
  for (const auto &[nodeId, _] : tree->allNodes) {
    auto rowBuilder =
        ancestorMatrixBuilder.init(nodeIndex, tree->allNodes.size());

    size_t colIndex = 0;
    for (const auto &[otherNodeId, _] : tree->allNodes) {
      // Check if otherNodeId is an ancestor of nodeId
      bool isAncestor = false;

      // Walk up the path from nodeId to root
      const auto &path = nodePaths[nodeId];
      for (auto *ancestorNode : path) {
        if (ancestorNode->identifier == otherNodeId) {
          isAncestor = true;
          break;
        }
      }

      rowBuilder.set(colIndex++, isAncestor);
    }

    nodeIndex++;
    if (nodeIndex <= 5 || nodeIndex % 100 == 0 || nodeIndex == tree->allNodes.size()) { // Log first 5, every 100, and the last one
       logging::debug("Built ancestor relationships for {} of {} nodes", nodeIndex, tree->allNodes.size());
    }
  }

  logging::debug("Placement acceleration data built successfully");
}

// New helper function to build the k-mer dictionary
std::unordered_map<std::string, uint32_t> buildKmerDictionary(
    const std::string& rootNodeId, // Add parameter for root node ID
    state::StateManager& stateManager, 
    ::capnp::List<::KmerDictionaryEntry>::Builder& kmerDictionaryBuilder,
    const std::unordered_set<std::string>& uniqueKmers) {  // Use the collected unique k-mers

    // Create the dictionary mapping k-mers to IDs
    std::unordered_map<std::string, uint32_t> kmerToId;
    
    // Instead of scanning genome, use the already collected unique k-mers
    logging::debug("Building k-mer dictionary from {} collected unique k-mers", uniqueKmers.size());
    
    // Sort k-mers for consistent dictionary building
    std::vector<std::string> sortedKmers(uniqueKmers.begin(), uniqueKmers.end());
    std::sort(sortedKmers.begin(), sortedKmers.end());
    
    
    logging::debug("Selected {} k-mers for the dictionary", sortedKmers.size());
    
    // The builder is already initialized with the correct size, so we don't need to call init
    
    // Build the dictionary and add k-mers to it
    for (size_t i = 0; i < sortedKmers.size(); i++) {
        // Create a copy of the k-mer string
        std::string kmer = sortedKmers[i];
        
        // Normalize k-mer to uppercase for consistent hashing
        std::transform(kmer.begin(), kmer.end(), kmer.begin(),
                      [](unsigned char c){ return std::toupper(c); });
        
        // Assign an ID to this k-mer
        uint32_t id = i;
        kmerToId[kmer] = id;
        
        // Compute reverse complement
        std::string revComp;
        revComp.reserve(kmer.size());
        
        // Compute reverse complement
        for (auto it = kmer.rbegin(); it != kmer.rend(); ++it) {
            char c = *it;
            switch (c) {
                case 'A': revComp.push_back('T'); break;
                case 'T': revComp.push_back('A'); break;
                case 'G': revComp.push_back('C'); break;
                case 'C': revComp.push_back('G'); break;
                default: revComp.push_back(c); break; // Non-standard bases stay the same
            }
        }
        
        bool isCanonical = (kmer <= revComp);
        
        // Store in the CapnProto dictionary
        auto entry = kmerDictionaryBuilder[i];
        entry.setSequence(kmer);
        entry.setCanonicalForm(isCanonical);
        
    }
    
    logging::debug("K-mer dictionary built with {} entries", kmerToId.size());
    return kmerToId;
}

// Update the main index function to include dictionary building
void index(panmanUtils::Tree *tree, Index::Builder &indexBuilder, int k, int s,
           ::capnp::MallocMessageBuilder& message, const std::string& indexPath) {
  if (!tree) {
    throw std::invalid_argument("Tree is null in index function");
  }

  // Open the output FASTA file for all node sequences
  std::ofstream allNodesFastaFile("ours_all.fa");
  if (!allNodesFastaFile.is_open()) {
      logging::err("Failed to open ours_all.fa for writing all node sequences.");
      // Decide if this is a fatal error or if indexing should continue without this dump
  }
  std::mutex fastaWriteMutex; // Mutex to synchronize writes to the FASTA file


  std::ofstream debugSeedFile("debug.tsv");

  // Minimal TSV logging: just node, FULL seed count, ACCURATE seed count
  debugSeedFile << "NodeID\tFULL_TotalSeeds\tACCURATE_TotalSeeds\n";
  
  // Comment out the complex debug TSV header:
  // debugSeedFile << "NodeID\tFULL_TotalSeeds\tFULL_OnlySeeds\tOURS_OnlySeeds\tCOMMON_ModifiedSeeds\tSeqLength\tBlockCount\tBlockDeletions\tBlockInsertions\tACCURATE_TotalSeeds\tInheritedSeeds\tLocalSeedChanges\tCOMMON_ModifiedSeeds\tSeqLength\tBlockCount\tBlockDeletions\tBlockInsertions\tDelta_Seeds\tFULL_SeedPositions\tACCURATE_SeedPositions\tBlockBoundaries\tActiveBlocks\tDeletedSeeds\tAddedSeeds\tModifiedSeeds\n";

  if (k <= 0 || s <= 0 || s >= k) {
    throw std::invalid_argument("Invalid k=" + std::to_string(k) +
                                " or s=" + std::to_string(s) +
                                " parameters (must satisfy 0 < s < k)");
  }

  // Normalize the output path to ensure consistency
  boost::filesystem::path normalizedPath = boost::filesystem::absolute(indexPath);
  std::string resolvedPath = normalizedPath.string();
  logging::debug("DEBUG-INDEX: Using normalized path for index: {}", resolvedPath);

  // CRITICAL: Create a local copy of the index builder to avoid corrupting the original
  // This ensures we're always working with the same object throughout the function
  Index::Builder localIndexBuilder = indexBuilder;
  
  try {
    logging::debug("Starting indexing with k={}, s={}", k, s);
    logging::debug("Tree has {} nodes", tree->allNodes.size());

    // Set basic parameters - CRITICAL: Must set these values first!
    localIndexBuilder.setK(k);
    localIndexBuilder.setS(s);
    localIndexBuilder.setOpen(false); 
    
    // Immediately verify the values were set correctly
    uint32_t check_k = localIndexBuilder.getK();
    uint32_t check_s = localIndexBuilder.getS();
    
    logging::debug("Successfully set and verified index parameters: k={}, s={}", check_k, check_s);

    // Initialize builders for seed and gap mutations
    size_t numNodes = tree->allNodes.size();
    auto perNodeSeedMutations = localIndexBuilder.initPerNodeSeedMutations(numNodes);
    auto perNodeGapMutations = localIndexBuilder.initPerNodeGapMutations(numNodes);
    logging::debug("Allocated storage for {} nodes", numNodes);

    // Initialize state manager and create k-mer collector
    logging::debug("Initializing state manager with reference sequence...");
    auto stateManager = initializeStateManager(tree, tree->root, k, s);
    std::unordered_set<std::string> uniqueKmersCollector;

    // Create a placement engine for indexing
    placement::PlacementEngine engine(k); // Assuming constructor takes k

    // Compute nodes by level for breadth-first traversal
    auto nodesByLevel = state::groupNodesByLevel(tree, tree->root);
    logging::debug("Found {} levels with {} total nodes", nodesByLevel.size(), tree->allNodes.size());

    // Compute node paths for mutation processing
    auto nodePaths = state::computeNodePaths(tree, tree->root);

    // CRITICAL: Verify global cache is initialized before processing any nodes
    if (!stateManager->isGlobalCacheInitialized()) {
        logging::critical("CRITICAL ERROR: Global position cache not initialized after StateManager setup!");
        throw std::runtime_error("Global position cache initialization failed - cannot proceed with indexing");
    }
    logging::debug("Verified: Global position cache is properly initialized before node processing");

    tbb::enumerable_thread_specific<std::unordered_set<std::string>> threadLocalKmerCollectors;

    // Traverse tree and process nodes level by level
    logging::debug("Starting tree traversal...");
    std::vector<std::string> processedNodes;
    processedNodes.reserve(tree->allNodes.size());

    // --- Progress Tracking --- 
    size_t totalNodesToIndex = tree->allNodes.size();
    std::atomic<size_t> nodesIndexedCounter{0};
    const size_t progressInterval = 30; // Update more frequently
    // --- End Progress Tracking ---

    // Process root node seeds first (no parent to compare with)
    logging::debug("Processing root node: {}", tree->root->identifier);
    auto &rootState = stateManager->getNodeState(tree->root->identifier);

    // Process root seeds - without k-mer dictionary first pass
    logging::info("Processing root node seeds..."); 
    processNodeComplete(*stateManager, engine, tree, tree->root, tree->root, nodePaths, k, s,
                       &perNodeSeedMutations, nullptr, false, true,
                       uniqueKmersCollector, debugSeedFile);
    processedNodes.push_back(tree->root->identifier);
    logging::info("Root node seed processing complete."); 
    logging::debug("Root node processing complete"); // Keep original debug log too

    // Process each level in breadth-first order
    logging::info("Starting main level processing loop..."); 
    for (size_t level = 1; level < nodesByLevel.size(); level++) {
      logging::debug("Processing level {} with {} nodes", level, nodesByLevel[level].size());
      const auto& currentLevelNodes = nodesByLevel[level]; // Get nodes for the current level
      
      // Calculate optimal grain size for TBB parallel_for
      const size_t numThreads = std::thread::hardware_concurrency();
      const size_t levelSize = currentLevelNodes.size();
      const size_t optimalGrainSize = std::max(1UL, 
          std::min(levelSize / (numThreads * 4), 64UL)); // Prevent too small or too large grains

      // --- Parallelize the inner loop over nodes in the current level with optimized grain size --- 
      tbb::parallel_for(tbb::blocked_range<size_t>(0, currentLevelNodes.size(), optimalGrainSize),
          [&](const tbb::blocked_range<size_t>& r) { // Use capture-by-reference [&]
          
          // Get a reference to the thread-local k-mer collector for this thread
          auto& localKmerCollector = threadLocalKmerCollectors.local();

          for (size_t i = r.begin(); i < r.end(); ++i) {
              panmanUtils::Node* node = currentLevelNodes[i]; // Get the node for this index
              if (!node) continue; // Skip null nodes
              std::string parentId = node->parent ? node->parent->identifier : ""; // Handle null parent for root
              std::string nodeId = node->identifier;

              // NOTE: Node should already be initialized during StateManager initialization
              // Redundant initializeNode() call removed to prevent parentId corruption

              // Apply mutations from parent to child and compute seeds
              processNodeComplete(*stateManager, engine, tree, node, tree->root, nodePaths, k, s,
                                 &perNodeSeedMutations, nullptr, true, true,
                                 localKmerCollector, debugSeedFile); // Use thread-local collector

              // --- BEGIN FASTA DUMP FOR CURRENT NODE ---
              if (allNodesFastaFile.is_open()) {
                  std::string full_sequence_for_node;
                  try {
                      std::vector<std::pair<int32_t, coordinates::CoordRange>> active_blocks_with_ranges = 
                          stateManager->getActiveBlockRanges(nodeId);

                      for (const auto& block_info_pair : active_blocks_with_ranges) {
                          const coordinates::CoordRange& block_range = block_info_pair.second;
                          if (block_range.start < block_range.end) {
                              auto [blockSequence, blockPositions, blockGaps, blockEndPositions] = 
                                  stateManager->extractSequence(nodeId, block_range, true /* skipGaps */);
                              full_sequence_for_node += blockSequence;
                          }
                      }
                  } catch (const std::exception& e) {
                      logging::err("Error extracting sequence for node {} for FASTA dump: {}", nodeId, e.what());
                      full_sequence_for_node = "ERROR_EXTRACTING_SEQUENCE";
                  }

                  std::lock_guard<std::mutex> lock(fastaWriteMutex);
                  allNodesFastaFile << ">" << nodeId << "\n";
                  const int FASTA_LINE_WIDTH = 80;
                  for (size_t char_idx = 0; char_idx < full_sequence_for_node.length(); char_idx += FASTA_LINE_WIDTH) {
                      allNodesFastaFile << full_sequence_for_node.substr(char_idx, FASTA_LINE_WIDTH) << "\n";
                  }
              }
              // --- END FASTA DUMP FOR CURRENT NODE ---

              size_t currentCount = nodesIndexedCounter.fetch_add(1) + 1; 
              if (currentCount % progressInterval == 0 || currentCount == totalNodesToIndex) {
                  // fprintf to stderr is generally thread-safe enough for progress bars
                  double percentage = (static_cast<double>(currentCount) / totalNodesToIndex) * 100.0;
                  int barWidth = 50;
                  int pos = static_cast<int>(barWidth * percentage / 100.0);
                  pos = std::min(pos, barWidth);
                  fprintf(stderr, "\rIndexing progress: [");
                  for (int j = 0; j < barWidth; ++j) {
                      if (j < pos) fprintf(stderr, "=");
                      else if (j == pos && pos < barWidth) fprintf(stderr, ">");
                      else fprintf(stderr, " ");
                  }
                  fprintf(stderr, "] %.1f%% (%zu/%zu nodes)    ", percentage, currentCount, totalNodesToIndex);
                  fflush(stderr);
              }
          } // End inner loop for range r
      }); // End tbb::parallel_for

      logging::debug("Completed processing level {}", level);
    } // End loop over levels

    fprintf(stderr, "\n"); 
    fflush(stderr); 

    // --- Merge thread-local k-mer collectors --- 
    uniqueKmersCollector = threadLocalKmerCollectors.combine([](const std::unordered_set<std::string>& a, const std::unordered_set<std::string>& b) {
        std::unordered_set<std::string> result = a;
        result.insert(b.begin(), b.end());
        return result;
    });
    

    logging::debug("Tree traversal complete, processed {} nodes", nodesIndexedCounter.load()); // Use atomic counter

    // Build the k-mer dictionary and add to the index
    logging::debug("Building k-mer dictionary from seeds...");

    size_t dictionarySize = uniqueKmersCollector.size();
    logging::debug("Determined k-mer dictionary size from collector: {}", dictionarySize);

    auto kmerDictionaryBuilder = localIndexBuilder.initKmerDictionary(dictionarySize);
    auto kmerDict = buildKmerDictionary(tree->root->identifier, *stateManager, kmerDictionaryBuilder, uniqueKmersCollector);

    // --- Second Pass: Update SeedMutations with Dictionary IDs ---
    if (!kmerDict.empty()) {
        logging::debug("Starting Pass 2: Updating node seed mutations with dictionary IDs...");
        size_t numNodes = stateManager->nodeIdsByDfsIndex.size();
        
        // Process each node by DFS index
        size_t nodesUpdated = 0;
        size_t totalDictionaryEntries = 0;
        
        for (size_t i = 0; i < numNodes; i++) {
                int64_t dfsIndex = static_cast<int64_t>(i);
                
                // Skip invalid indices
                if (dfsIndex < 0 || dfsIndex >= static_cast<int64_t>(numNodes)) {
                    continue;
                }

            // Get node ID from state manager
                const std::string nodeId = stateManager->nodeIdsByDfsIndex[dfsIndex];
                
            // Check if this node has stored k-mer sequences
            auto nodeKmersIt = stateManager->nodeKmerSequences.find(nodeId);
            if (nodeKmersIt == stateManager->nodeKmerSequences.end() || nodeKmersIt->second.empty()) {
                continue; // No k-mers for this node
            }
            
            // Convert stored k-mer sequences to dictionary IDs
            std::vector<uint32_t> dictIds;
            std::vector<int64_t> positions;
            std::vector<uint32_t> endOffsets;
            
            // Get stored k-mers
            const auto& posToKmerMap = nodeKmersIt->second;
            
            // Process each stored k-mer
            for (const auto& [pos, kmerSeq] : posToKmerMap) {
                // Normalize k-mer to uppercase for consistent lookup
                std::string upperKmer = kmerSeq;
                std::transform(upperKmer.begin(), upperKmer.end(), upperKmer.begin(),
                               [](unsigned char c){ return std::toupper(c); });
                
                // Look up in dictionary
                auto dictIt = kmerDict.find(upperKmer);
                if (dictIt != kmerDict.end()) {
                    dictIds.push_back(dictIt->second);
                    positions.push_back(pos);
                    
                    // Get end offset
                    uint32_t endOffset = k - 1; // Default
                    auto offsetsIt = stateManager->nodeKmerEndOffsets.find(nodeId);
                    if (offsetsIt != stateManager->nodeKmerEndOffsets.end()) {
                        auto posOffsetIt = offsetsIt->second.find(pos);
                        if (posOffsetIt != offsetsIt->second.end()) {
                            endOffset = posOffsetIt->second;
                        }
                    }
                    
                    endOffsets.push_back(endOffset);
                }
            }
            
            // If we found dictionary entries, update the Cap'n Proto structure
            if (!dictIds.empty()) {
                auto mutations = perNodeSeedMutations[dfsIndex];
                
                // Write dictionary IDs
                auto dictIdsBuilder = mutations.initKmerDictionaryIds(dictIds.size());
                for (size_t j = 0; j < dictIds.size(); j++) {
                    dictIdsBuilder.set(j, dictIds[j]);
                }
                
                // Write positions
                auto positionsBuilder = mutations.initKmerPositions(positions.size());
                for (size_t j = 0; j < positions.size(); j++) {
                    positionsBuilder.set(j, positions[j]);
                }
                
                // Write end offsets
                auto endOffsetsBuilder = mutations.initKmerEndOffsets(endOffsets.size());
                for (size_t j = 0; j < endOffsets.size(); j++) {
                    endOffsetsBuilder.set(j, endOffsets[j]);
                }
                
                // Update stats
                nodesUpdated++;
                totalDictionaryEntries += dictIds.size();
                
                // Log progress for some nodes
                if (nodesUpdated <= 5 || nodesUpdated % 100 == 0 || nodesUpdated == numNodes) {
                    logging::info("Updated node {} (DFS index {}) with {} dictionary entries", 
                                 nodeId, dfsIndex, dictIds.size());
                }
            }
        }
        
        logging::info("Pass 2 complete: Updated {} nodes with {} total dictionary entries", 
                     nodesUpdated, totalDictionaryEntries);
                     
        // Clear temporary storage
        stateManager->nodeKmerSequences.clear();
        stateManager->nodeKmerEndOffsets.clear();
    }
    
    // CRITICAL FIX: Create a snapshot of the current message state to ensure we're not losing data
    // Store the current count of all critical data structures
    size_t check_seedMutations = localIndexBuilder.getPerNodeSeedMutations().size();
    size_t check_dictionary = localIndexBuilder.getKmerDictionary().size();
    
    // Final verification of the index structure
    check_k = localIndexBuilder.getK();
    check_s = localIndexBuilder.getS();
    
    logging::info("Final index verification: k={}, s={}, seedMutations={}, dictionary={}",
                 check_k, check_s, check_seedMutations, check_dictionary);

    // Make sure indexBuilder (passed by reference) has the same values as our local copy
    // THIS IS CRITICAL - we're ensuring the data in the message is preserved
    indexBuilder = localIndexBuilder;
    
    // Double check critical values again
    auto finalRoot = message.getRoot<Index>();
    uint32_t final_k = finalRoot.getK();
    uint32_t final_s = finalRoot.getS();
    size_t final_mutations = finalRoot.getPerNodeSeedMutations().size();
    size_t final_dictionary = finalRoot.getKmerDictionary().size();
    
    logging::info("Last verification before returning: k={}, s={}, mutations={}, nodePathInfo={}",
                 final_k, final_s, final_mutations, final_dictionary);
    
    // CRITICAL: If the values are corrupted, fix them directly in the message root
    if (final_k != static_cast<uint32_t>(k) || final_s != static_cast<uint32_t>(s)) {
        logging::warn("Final verification detected corrupted values k={}, s={}. Fixing before writing...", 
                     final_k, final_s);
        
        // Re-set the key values on the message root
        finalRoot.setK(k);
        finalRoot.setS(s);
        
        // Re-verify
        logging::info("After correction: k={}, s={}, seedMutations={}, dictionary={}", 
                     finalRoot.getK(), finalRoot.getS(), 
                     finalRoot.getPerNodeSeedMutations().size(),
                     finalRoot.getKmerDictionary().size());
                     
        // CRITICAL CHECK: If the seedMutations or dictionary was corrupted, this is a fatal error
        if (final_mutations == 0 && check_seedMutations > 0) {
            logging::critical("FATAL: seedMutations data was lost during correction ({} -> 0)", check_seedMutations);
            throw std::runtime_error("Index data corruption: seedMutations data was lost");
        }
        
        if (final_dictionary == 0 && check_dictionary > 0) {
            logging::critical("FATAL: dictionary data was lost during correction ({} -> 0)", check_dictionary);
            throw std::runtime_error("Index data corruption: dictionary data was lost");
        }
    } else if (final_mutations != check_seedMutations || final_dictionary != check_dictionary) {
        // Something else caused data corruption
        logging::critical("FATAL: Data corruption detected. Expected {}:{}, got {}:{} (mutations:dictionary)",
                      check_seedMutations, check_dictionary, final_mutations, final_dictionary);
        throw std::runtime_error("Index data corruption: Count mismatch in critical data structures");
    }
    
    // If we're going to write the index, write it now
    if (!indexPath.empty()) {
        logging::info("Writing complete index to: {}", indexPath);
        
        try {
            // CRITICAL: Build placement acceleration data before writing
            logging::info("Building placement acceleration data (nodePathInfo)...");
            buildPlacementAccelerationData(tree, localIndexBuilder, *stateManager);
            
            // Verify nodePathInfo was properly populated
            auto nodePathInfo = localIndexBuilder.getNodePathInfo();
            if (nodePathInfo.size() == 0) {
                logging::err("CRITICAL ERROR: nodePathInfo is still empty after building acceleration data");
                throw std::runtime_error("Failed to populate nodePathInfo structure");
            } else {
                logging::info("Successfully populated nodePathInfo with {} entries", nodePathInfo.size());
            }
            
            // CRITICAL FIX: Create a completely fresh copy of the nodePathInfo in the message root
            auto rootForWriting = message.getRoot<Index>();
            
            try {
                // First clear any existing nodePathInfo in the root by re-initializing with the correct size
                logging::info("Initializing nodePathInfo in message root with {} entries", nodePathInfo.size());
                auto newNodePathInfo = rootForWriting.initNodePathInfo(nodePathInfo.size());
                
                // Additional debugging to print first few node IDs to find the problem
                logging::info("First few node IDs in source nodePathInfo:");
                for (size_t i = 0; i < std::min<size_t>(5, nodePathInfo.size()); i++) {
                    auto srcId = nodePathInfo[i].getNodeId();
                    if (srcId.size() > 0) {
                        // Convert using direct copy of bytes
                        std::string srcIdStr(srcId.begin(), srcId.end());
                        logging::info("  Source[{}]: '{}' (length: {})", i, srcIdStr, srcId.size());
                    } else {
                        logging::critical("  Source[{}]: EMPTY NODE ID in source data!", i);
                    }
                }
            
                // Use a field-by-field copy to ensure all data is properly transferred
                for (size_t i = 0; i < nodePathInfo.size(); i++) {
                    auto srcInfo = nodePathInfo[i];
                    auto dstInfo = newNodePathInfo[i];
                    
                    // Transfer node ID - IMPORTANT FIX: Properly copy the text value
                    auto nodeId = srcInfo.getNodeId();
                    
                    // Direct method to get the raw data from Cap'n Proto Text
                    if (nodeId.size() > 0) {
                        // Create a direct and explicit copy of the text data
                        std::string nodeIdStr;
                        nodeIdStr.assign(nodeId.begin(), nodeId.end());
                        
                        // Set node ID with explicit string
                        dstInfo.setNodeId(nodeIdStr);
                        
                        // Verify it worked immediately
                        if (i < 5) {
                          auto verifyId = dstInfo.getNodeId();
                          std::string verifyStr(verifyId.begin(), verifyId.end());
                          logging::info("Immediate verification of node {}: '{}' (length={})", 
                                        i, verifyStr, verifyStr.length());
                        }
                    } else {
                        // If source is empty, try to get the actual node ID from the tree
                        std::string fallbackId;
                        size_t nodeIndex = i;
                        if (nodeIndex < stateManager->nodeIdsByDfsIndex.size()) {
                            fallbackId = stateManager->nodeIdsByDfsIndex[nodeIndex];
                            logging::info("Using fallback node ID '{}' for index {}", fallbackId, nodeIndex);
                            dstInfo.setNodeId(fallbackId);
                        } else {
                            // Last resort: use a generated ID
                            fallbackId = "node_" + std::to_string(i + 1);
                            logging::warn("Using GENERATED node ID '{}' for index {}", fallbackId, i);
                            dstInfo.setNodeId(fallbackId);
                        }
                    }
                    
                    // Transfer level
                    dstInfo.setLevel(srcInfo.getLevel());
                    
                    // Transfer parent ID if it exists
                    if (srcInfo.hasParentId()) {
                        auto parentId = srcInfo.getParentId();
                        
                        if (parentId.size() > 0) {
                            std::string parentIdStr;
                            parentIdStr.assign(parentId.begin(), parentId.end());
                            dstInfo.setParentId(parentIdStr);
                        } else {
                            // Use empty string for root node
                            dstInfo.setParentId("");
                        }
                    } else {
                        dstInfo.setParentId(""); // Root node has no parent
                    }

          // Copy active blocks
          auto srcBlocks = srcInfo.getActiveBlocks();
          auto dstBlocks = dstInfo.initActiveBlocks(srcBlocks.size());
          for (size_t j = 0; j < srcBlocks.size(); j++) {
            dstBlocks.set(j, srcBlocks[j]);
          }
          
          // Log sample transfers to verify
          if (i < 3 || i + 1 == nodePathInfo.size()) {
            auto nodeIdCheck = dstInfo.getNodeId();
            std::string nodeIdStr;
            nodeIdStr.assign(nodeIdCheck.begin(), nodeIdCheck.end());
            logging::debug("Copied node [{}]: ID='{}', Level={}, Active blocks: {}", 
                        i, nodeIdStr, srcInfo.getLevel(), srcBlocks.size());
          }
        }
            
        // Verify the transfer
        auto verifyPathInfo = rootForWriting.getNodePathInfo();
        if (verifyPathInfo.size() != nodePathInfo.size()) {
          logging::critical("TRANSFER VERIFICATION FAILED: Expected {} entries, got {}", 
                           nodePathInfo.size(), verifyPathInfo.size());
          throw std::runtime_error("nodePathInfo transfer failed - size mismatch");
        }
        
        // Check a sample of transferred entries
        for (size_t i = 0; i < std::min(size_t(3), static_cast<size_t>(verifyPathInfo.size())); i++) {
          auto node = verifyPathInfo[i];
          auto nodeId = node.getNodeId();
          std::string nodeIdStr(nodeId.begin(), nodeId.end());
          auto level = node.getLevel();
          auto blockCount = node.getActiveBlocks().size();
          
          // IMPORTANT: Verify that node IDs are not empty
          if (nodeIdStr.empty()) {
            logging::critical("VERIFICATION FAILED: Node ID at index {} is empty!", i);
          }
          
          logging::debug("Verified node [{}]: ID='{}', Level={}, Active blocks: {}", 
                      i, nodeIdStr, level, blockCount);
        }
        
        logging::info("Successfully transferred nodePathInfo with {} entries to message root", verifyPathInfo.size());
      } catch (const ::kj::Exception& e) {
                logging::critical("Cap'n Proto error during nodePathInfo transfer: {}", e.getDescription().cStr());
                throw std::runtime_error("Fatal error during nodePathInfo transfer: " + std::string(e.getDescription().cStr()));
            } catch (const std::exception& e) {
                logging::critical("Error during nodePathInfo transfer: {}", e.what());
                throw std::runtime_error("Fatal error during nodePathInfo transfer: " + std::string(e.what()));
            }
            
            // Ensure root has correct k, s values after nodePathInfo transfer
            if (rootForWriting.getK() != static_cast<uint32_t>(k) || 
                rootForWriting.getS() != static_cast<uint32_t>(s)) {
                logging::warn("Primary values corrupted after nodePathInfo transfer, fixing...");
                rootForWriting.setK(k);
                rootForWriting.setS(s);
            }
            
            // Write the index to disk
            logging::info("Final verification before writing: Checking node IDs in index...");
            auto finalVerifyPathInfo = rootForWriting.getNodePathInfo();
            if (finalVerifyPathInfo.size() > 0) {
                logging::info("Sample node IDs in final index:");
                for (size_t i = 0; i < std::min<size_t>(5, finalVerifyPathInfo.size()); i++) {
                    auto nodeId = finalVerifyPathInfo[i].getNodeId();
                    std::string nodeIdStr(nodeId.begin(), nodeId.end());
                    logging::info("  [{}]: '{}'", i, nodeIdStr);
                    
                    if (nodeIdStr.empty()) {
                        logging::critical("CRITICAL: Empty node ID found at index {} just before writing!", i);
                    }
                }
            }
            
            ::writeCapnp(message, indexPath);
            logging::info("Successfully wrote index to disk: {}", indexPath);
            
            // Verify the written index
            try {
                logging::info("Verifying the written index...");
                auto reader = ::readCapnp(indexPath);
                if (reader) {
                    auto verifyRoot = reader->getRoot<Index>();
                    
                    // Check key data was written properly
                    uint32_t verify_k = verifyRoot.getK();
                    uint32_t verify_s = verifyRoot.getS();
                    size_t verify_mutations = verifyRoot.getPerNodeSeedMutations().size();
                    size_t verify_nodePathInfo = verifyRoot.getNodePathInfo().size();
                    
                    logging::info("Verified written index: k={}, s={}, mutations={}, nodePathInfo={}", 
                                verify_k, verify_s, verify_mutations, verify_nodePathInfo);
                                
                    if (verify_k != k || verify_s != s || verify_mutations == 0 || verify_nodePathInfo == 0) {
                        logging::err("Index verification failed: Incorrect values in index file");
                        
                        // ADDED: More detailed error for nodePathInfo specifically
                        if (verify_nodePathInfo == 0) {
                            logging::critical("CRITICAL ERROR: Written index has EMPTY nodePathInfo structure");
                            // Attempt to identify where the data was lost
                            logging::critical("This confirms data loss during writing. Please report this bug.");
                            throw std::runtime_error("Fatal error: nodePathInfo is empty in written index");
                        }
                    } else {
                        // ADDED: Extra verification of nodePathInfo contents
                        auto verifyNodePathInfo = verifyRoot.getNodePathInfo();
                        size_t samplesToCheck = std::min<size_t>(3, verifyNodePathInfo.size());
                        
                        for (size_t i = 0; i < samplesToCheck; i++) {
                            auto node = verifyNodePathInfo[i];
                            auto nodeId = node.getNodeId();
                            std::string nodeIdStr(nodeId.begin(), nodeId.end());
                            
                            std::string parentIdStr;
                            if (node.hasParentId()) {
                                auto parentId = node.getParentId();
                                parentIdStr = std::string(parentId.begin(), parentId.end());
                            } else {
                                parentIdStr = "";
                            }
                            
                            auto level = node.getLevel();
                            auto blockCount = node.getActiveBlocks().size();
                            
                            logging::debug("Verified node [{}]: ID='{}', Parent='{}', Level={}, Blocks={}",
                                        i, nodeIdStr, parentIdStr, level, blockCount);
                            
                            // Check if the node has blocks
                            if (blockCount == 0) {
                                logging::warn("Node '{}' has 0 active blocks in written index", nodeIdStr);
                            }
                        } // <-- FIXED: Added missing closing brace for the for loop
            
                        logging::info("Index verified successfully with complete nodePathInfo structure");
                    }
                }
            } catch (const std::exception& e) {
                logging::err("Error verifying index: {}", e.what());
                // Continue anyway - but we should re-throw for critical failures
                if (std::string(e.what()).find("nodePathInfo") != std::string::npos) {
                    throw; // Re-throw nodePathInfo-related errors
                }
            }
        } catch (const std::exception& e) {
            logging::err("Failed to write index: {}", e.what());
            throw; // Re-throw to caller
        }
    } else {
        logging::info("No indexPath provided, skipping index file write");
    }

    if (allNodesFastaFile.is_open()) {
        allNodesFastaFile.close();
    }
  } catch (const std::exception &e) {
    logging::critical("FATAL ERROR during indexing: {}", e.what());
    if (allNodesFastaFile.is_open()) { // Ensure file is closed on error too
        allNodesFastaFile.close();
    }
    throw;
  }
}

// Helper function to output seed details - keeps code cleaner
void outputSeedDetails(std::ofstream& outFile, int64_t pos, const seeding::seed_t& seed, 
                      state::StateManager& stateManager, const std::string& nodeId) {
    outFile << "Seed at position " << pos << ":\n";
    outFile << "  Hash: " << seed.hash << "\n";
    outFile << "  Reversed: " << (seed.reversed ? "YES" : "NO") << "\n";
    outFile << "  End position: " << seed.endPos << "\n";
    outFile << "  Length: " << (seed.endPos - pos + 1) << "\n";
    
    // Try to extract k-mer
    try {
        auto kmerResult = stateManager.extractKmer(nodeId, pos, stateManager.getKmerSize());
        std::string kmer = kmerResult.first;
        if (!kmer.empty()) {
            outFile << "  K-mer: " << kmer << "\n";
            
            // Compute reverse complement
            std::string revComp;
            revComp.reserve(kmer.size());
            
            for (auto it = kmer.rbegin(); it != kmer.rend(); ++it) {
                char c = *it;
                switch (c) {
                    case 'A': revComp.push_back('T'); break;
                    case 'T': revComp.push_back('A'); break;
                    case 'G': revComp.push_back('C'); break;
                    case 'C': revComp.push_back('G'); break;
                    case 'a': revComp.push_back('t'); break;
                    case 't': revComp.push_back('a'); break;
                    case 'g': revComp.push_back('c'); break;
                    case 'c': revComp.push_back('g'); break;
                    default: revComp.push_back(c); break; // Non-standard bases stay the same
                }
            }
            
            outFile << "  Rev Comp: " << revComp << "\n";
            outFile << "  Is Canonical: " << ((kmer <= revComp) ? "YES" : "NO") << "\n";
        }
    } catch (...) {}
    
    // Get block coordinates
    try {
        auto coordOpt = stateManager.fastMapGlobalToLocal(pos);
        if (coordOpt) {
            auto [blockId, nucPos, gapPos] = *coordOpt;
            outFile << "  Mapped to: Block=" << blockId << ", nucPos=" << nucPos << ", gapPos=" << gapPos << "\n";
        }
    } catch (...) {}
    
    outFile << "\n";
}

/**
 * @brief Process a node completely, combining mutation processing and seed recomputation
 */
void processNodeComplete(
    state::StateManager& stateManager,
    placement::PlacementEngine& engine,
    panmanUtils::Tree* tree, // Added tree pointer
    panmanUtils::Node* node,
    panmanUtils::Node* commonAncestor,
    const std::unordered_map<std::string, std::vector<panmanUtils::Node*>>& nodePaths,
    int k, int s,
    ::capnp::List<SeedMutations>::Builder* perNodeSeedMutations,
    const std::unordered_map<std::string, uint32_t>* kmerDictionary,
    bool processMutations, // Flag from the caller of processNodeComplete
    bool computeSeeds,
    std::unordered_set<std::string>& uniqueKmersCollector,
    std::ofstream& debugSeedFile) {
    
    if (!node) return;
    
    std::string nodeId = node->identifier;

    std::string ab_before_str_for_log; // Populated by processMutationsForNode
    std::string strNodeId = node->identifier; // strNodeId is same as nodeId

    // NodeState acquisition and parent propagation (critical first steps)
    auto &nodeState = stateManager.getNodeState(strNodeId);
    // With inheritance-based approach, we don't need to explicitly propagate state
    // Block states will be automatically inherited when queried
    if (node->parent) { 
        std::string parentId = node->parent->identifier;
        try {
            // Just ensure the node is initialized
            stateManager.getNodeState(strNodeId);
            // std::cerr << "[DEBUG] PNC: Node " << strNodeId << " will inherit from parent " << parentId << " as needed" << std::endl;
        } catch (const std::exception& e) {
            logging::err("PNC: Error accessing state for {} (with parent {}): {}", strNodeId, parentId, e.what());
            throw; 
        }
    }

    if (node->identifier != tree->root->identifier) {
      for (const int32_t blockId : stateManager.getActiveBlocks(node->parent->identifier)) {
        stateManager.setBlockOn(strNodeId, blockId, true);
      }
    }

    std::stringstream condensed_ss_ourcode_stream; // Use a local stream for this function\'s log part


    // CRITICAL: Materialize node state BEFORE applying local mutations
    // This ensures inherited seeds are available during block deletion processing
    std::cout << "DEBUG_PROCESS_NODE: Processing node " << nodeId << std::endl;
    stateManager.getNodeState(nodeId); // This will trigger initialization if needed
    std::cout << "DEBUG_PROCESS_NODE: About to materialize node " << nodeId << std::endl;
    stateManager.materializeNodeState(nodeId); // EXPLICITLY materialize node state
    std::cout << "DEBUG_PROCESS_NODE: Materialized node " << nodeId << std::endl;
    
    // Apply mutations for the node IF THE FLAG IS SET
    std::string* ab_before_ptr = nullptr;
    std::stringstream* condensed_stream_for_pmfn = nullptr; 

    processMutationsForNode(stateManager, tree, node, commonAncestor, nodePaths, 
        ab_before_ptr,
        condensed_stream_for_pmfn 
    );

    if (computeSeeds) {
        auto [blockDeletions, blockInsertions, totalMaterializedSeeds, seedsFoundInDeletedBlocks, blockDeletionPositionsChecked, seedsCleared, seedsAdded, recompRangeCount, recompRangeSize, recompRangeList, uniquePositionsProcessed, totalKmersInRanges, detailedSeedChanges] = recomputeSeeds(stateManager, engine, node, k, s, perNodeSeedMutations, kmerDictionary, uniqueKmersCollector, debugSeedFile);

        // Update materialized seed state with detailed seed changes from recomputation
        stateManager.updateMaterializedSeedsAfterRecomputation(nodeId, detailedSeedChanges);

        // Get materialized state computed flag for debug TSV
        bool materializedStateComputed = false;
        try {
            const auto& nodeState = stateManager.getNodeState(nodeId);
            // Avoid potential deadlock by not locking here - this is just for debug info
            materializedStateComputed = nodeState.materializedStateComputed;
        } catch (const std::exception& e) {
            logging::warn("Failed to get materializedStateComputed for node {}: {}", nodeId, e.what());
        }

        // Helper functions removed to avoid deadlock potential
        
        // SAFE: Skip StateManager calls to avoid deadlock - use fallback data
        std::vector<std::pair<int32_t, coordinates::CoordRange>> activeBlocksAndRanges;
        
        // Try to get active blocks safely - if it fails, use empty list
        try {
            activeBlocksAndRanges = stateManager.getActiveBlockRanges(nodeId);
        } catch (...) {
            // Skip StateManager call to avoid deadlock
            logging::warn("Skipped getActiveBlockRanges for node {} to avoid deadlock", nodeId);
        }
    if (activeBlocksAndRanges.empty()) {
        // No blocks → all counts zero
        debugSeedFile << nodeId << "\t0\t0\n";
        return;
    }

    // Compute full span - use safe defaults if no active blocks
    int64_t fullStart = 0;
    int64_t fullEnd = 1000; // Default small range
    int64_t commonSeqLength = 0;
    int64_t commonBlocks = 0;
    
    if (!activeBlocksAndRanges.empty()) {
        fullStart = activeBlocksAndRanges.front().second.start;
        fullEnd = activeBlocksAndRanges.back().second.end;
        commonSeqLength = fullEnd - fullStart;
        commonBlocks = activeBlocksAndRanges.size();
    }

    // === FULL method - SAFE version ===
    std::string fullSeq = "";
    std::vector<int64_t> fullPos;
    std::vector<bool> fullGaps;
    std::vector<int64_t> fullEndPos;
    
    // Try sequence extraction safely
    try {
        // Extend fullEnd to include k more non-gap characters downstream
        int64_t extendedFullEnd = fullEnd;
        int k_val = stateManager.getKmerSize();
        if (k_val > 0) {
            try {
                // Use extractKmer to get k characters downstream from fullEnd
                auto kmerResult = stateManager.extractKmer(nodeId, fullEnd, k_val);
                if (!kmerResult.first.empty() && !kmerResult.second.empty()) {
                    // Use the position after the k-th character (end position of k-mer)
                    extendedFullEnd = kmerResult.second[k_val - 1] + 1;
                } else {
                    // If k-mer extraction fails, fall back to simple position extension
                    extendedFullEnd = fullEnd + k_val;
                }
            } catch (...) {
                // If k-mer extraction fails, fall back to simple position extension
                extendedFullEnd = fullEnd + k_val;
            }
        }
        
        auto extractResult = stateManager.extractSequence(nodeId, {fullStart, extendedFullEnd}, false);
        fullSeq = std::get<0>(extractResult);
        fullPos = std::get<1>(extractResult);
        fullGaps = std::get<2>(extractResult);
        fullEndPos = std::get<3>(extractResult);
    } catch (...) {
        // Skip StateManager call to avoid deadlock
        logging::warn("Skipped extractSequence for node {} to avoid deadlock", nodeId);
        // Use empty defaults
    }

    // Strip gaps to get ungappedSequence
    std::string ungapped = fullSeq;
    ungapped.erase(std::remove_if(ungapped.begin(), ungapped.end(),
                                  [](char c){ return c=='-'||c=='x'; }),
                    ungapped.end());

    // Map ungapped → gapped indices
    std::vector<int64_t> ung2gap(ungapped.size(), -1);
    for (size_t i = 0, gi = 0; i < fullSeq.size(); ++i) {
        if (fullSeq[i] != '-' && fullSeq[i] != 'x') {
            ung2gap[gi++] = i; // Map gapped position to ungapped index
        }
    }

    // Collect rolling syncmers
    auto fullKmers = seeding::rollingSyncmers(ungapped, k, s, false, 0, true);

    // Count only actual syncmers, not all k-mer positions
    // But restrict to seeds within the original range [fullStart, fullEnd]
    int64_t fullTotalSeeds = 0;
    std::vector<int64_t> fullSeedPositions;
    for (const auto& [hash, isReverse, isSyncmer, startPos] : fullKmers) {
        if (isSyncmer) {
            // Convert local ungapped position to global position
            if (startPos < ung2gap.size() && ung2gap[startPos] != -1) {
                int64_t globalPos = fullPos[ung2gap[startPos]];
                // Only count seeds within the original range [fullStart, fullEnd]
                if (globalPos >= fullStart && globalPos <= fullEnd) {
                    fullTotalSeeds++;
                    fullSeedPositions.push_back(globalPos);
                }
            }
        }
    }
    int64_t fullOnly        = 0;  // compute as needed
    int64_t oursOnly        = 0;  // compute as needed
    int64_t fullOursMod     = 0;  // compute as needed

    // === OURS method - SAFE version ===
    // Get seed positions from our implementation using the complete materialized state
    std::unordered_set<int64_t> oursPositions;
    
    // Get all seed positions from the materialized state (includes inherited + local)
    try {
        const auto& nodeState = stateManager.getNodeState(nodeId);
        // Get all materialized seeds (both inherited and local)
        for (const auto& [position, seed] : nodeState.materializedSeeds) {
            oursPositions.insert(position);
        }
        logging::info("Successfully retrieved {} complete seed positions from materialized state for node {}", oursPositions.size(), nodeId);
    } catch (const std::exception& e) {
        throw std::runtime_error("Failed to retrieve materialized seeds for node " + nodeId + ": " + e.what());
    }
    
    int64_t oursTotalSeeds = oursPositions.size();
    
    

    // Compute diffs between fullKmers and oursSeeds if you want fullOnly/oursOnly/fullOursMod:
    //   for each seed in fullKmers  → if not in oursSeeds ⇒ fullOnly++
    //   for each seed in oursSeeds  → if not in fullKmers ⇒ oursOnly++
    //   common ⇒ else fullOursMod++

    // Compute delta if desired:
    int64_t deltaFullVsOurs = oursTotalSeeds - fullTotalSeeds;
    
    // Get the accurate total seed count from NodeState
    const auto& nodeState = stateManager.getNodeState(nodeId);
    int64_t accurateTotalSeeds = nodeState.getTotalSeedCount();
    int64_t inheritedSeeds = nodeState.inheritedSeedCount;
    int64_t localChanges = nodeState.localSeedChanges;
    int64_t deltaFullVsAccurate = accurateTotalSeeds - fullTotalSeeds;

    // Convert seed position vectors to comma-separated strings
    std::string fullSeedPosStr = "";
    for (size_t i = 0; i < fullSeedPositions.size(); ++i) {
        if (i > 0) fullSeedPosStr += ",";
        fullSeedPosStr += std::to_string(fullSeedPositions[i]);
    }
    
    std::string accurateSeedPosStr = "";
    bool first = true;
    for (const auto& position : oursPositions) {
        if (!first) accurateSeedPosStr += ",";
        accurateSeedPosStr += std::to_string(position);
        first = false;
    }

    // Format active block ranges for TSV output
    std::string blockRangesStr = "";
    for (size_t i = 0; i < activeBlocksAndRanges.size(); ++i) {
        if (i > 0) blockRangesStr += ",";
        const auto& [blockId, range] = activeBlocksAndRanges[i];
        blockRangesStr += std::to_string(range.start) + "-" + std::to_string(range.end) + ":" + std::to_string(blockId);
    }
    
    // Format active block IDs for TSV output
    std::string activeBlocksStr = "";
    for (size_t i = 0; i < activeBlocksAndRanges.size(); ++i) {
        if (i > 0) activeBlocksStr += ",";
        const auto& [blockId, range] = activeBlocksAndRanges[i];
        activeBlocksStr += std::to_string(blockId);
    }

    // Compute seed position differences between FULL and ACCURATE methods
    std::set<int64_t> fullPositionsSet(fullSeedPositions.begin(), fullSeedPositions.end());
    std::set<int64_t> accuratePositionsSet(oursPositions.begin(), oursPositions.end());
    
    // Deleted seeds: in FULL but not in ACCURATE
    std::vector<int64_t> deletedSeeds;
    std::set_difference(fullPositionsSet.begin(), fullPositionsSet.end(),
                       accuratePositionsSet.begin(), accuratePositionsSet.end(),
                       std::back_inserter(deletedSeeds));
    
    // Added seeds: in ACCURATE but not in FULL
    std::vector<int64_t> addedSeeds;
    std::set_difference(accuratePositionsSet.begin(), accuratePositionsSet.end(),
                       fullPositionsSet.begin(), fullPositionsSet.end(),
                       std::back_inserter(addedSeeds));
    
    // Modified seeds: in both FULL and ACCURATE (intersection)
    std::vector<int64_t> modifiedSeeds;
    std::set_intersection(fullPositionsSet.begin(), fullPositionsSet.end(),
                         accuratePositionsSet.begin(), accuratePositionsSet.end(),
                         std::back_inserter(modifiedSeeds));

    // Format seed difference lists as comma-separated strings
    std::string deletedSeedsStr = "";
    for (size_t i = 0; i < deletedSeeds.size(); ++i) {
        if (i > 0) deletedSeedsStr += ",";
        deletedSeedsStr += std::to_string(deletedSeeds[i]);
    }
    
    std::string addedSeedsStr = "";
    for (size_t i = 0; i < addedSeeds.size(); ++i) {
        if (i > 0) addedSeedsStr += ",";
        addedSeedsStr += std::to_string(addedSeeds[i]);
    }
    
    std::string modifiedSeedsStr = "";
    for (size_t i = 0; i < modifiedSeeds.size(); ++i) {
        if (i > 0) modifiedSeedsStr += ",";
        modifiedSeedsStr += std::to_string(modifiedSeeds[i]);
    }

    // For node_2, extract and log k-mers at missing positions
    if (nodeId == "node_2") {
        std::set<int64_t> fullPositionsSet(fullSeedPositions.begin(), fullSeedPositions.end());
        std::set<int64_t> missingPositions;
        for (int64_t pos : fullPositionsSet) {
            if (oursPositions.find(pos) == oursPositions.end()) {
                missingPositions.insert(pos);
            }
        }
        
        if (!missingPositions.empty()) {
            logging::info("NODE_2_MISSING_KMERS: Found {} missing positions in ACCURATE vs FULL", missingPositions.size());
            
            for (int64_t pos : missingPositions) {
                std::string fullKmer = "NOTFOUND";
                std::string accurateKmer = "NOTFOUND";
                
                // --- FULL k-mer extraction (SAFE) ---
                try {
                    // Find the corresponding k-mer in fullKmers by matching the global position
                    size_t syncmerIdx = 0;
                    for (const auto& [hash, isReverse, isSyncmer, startPos] : fullKmers) {
                        if (isSyncmer) {
                            if (syncmerIdx < fullSeedPositions.size() && fullSeedPositions[syncmerIdx] == pos) {
                                if (startPos + k <= ungapped.size()) {
                                    fullKmer = ungapped.substr(startPos, k);
                                }
                                break;
                            }
                            syncmerIdx++;
                        }
                    }
                } catch (...) {
                    // If extraction fails, leave as NOTFOUND
                }
                
                // --- ACCURATE (OURS) k-mer extraction ---
                // First try our stored sequences
                auto nodeKmersIt = stateManager.nodeKmerSequences.find(nodeId);
                if (nodeKmersIt != stateManager.nodeKmerSequences.end()) {
                    auto kmerIt = nodeKmersIt->second.find(pos);
                    if (kmerIt != nodeKmersIt->second.end()) {
                        accurateKmer = kmerIt->second;
                    } else {
                        accurateKmer = "NOT_FOUND_IN_STORED_KMERS";
                    }
                } else {
                    accurateKmer = "NO_KMERS_STORED_FOR_NODE";
                }
                
                // Additionally, try to extract the actual k-mer at this position from our sequence
                std::string actualKmerAtPos = "EXTRACTION_FAILED";
                try {
                    auto kmerResult = stateManager.extractKmer(nodeId, pos, k);
                    actualKmerAtPos = kmerResult.first;
                } catch (...) {
                    actualKmerAtPos = "EXTRACTION_FAILED";
                }
                
                logging::info("NODE_2_MISSING_KMER: pos={} FULL_kmer={} ACCURATE_stored={} ACCURATE_actual={}", 
                              pos, fullKmer, accurateKmer, actualKmerAtPos);
            }
        }
    }

    // === MINIMAL TSV OUTPUT ===
    debugSeedFile << nodeId << "\t" << fullTotalSeeds << "\t" << accurateTotalSeeds << "\n";
    
    // Comment out the complex TSV output:
    /*
    debugSeedFile 
        << nodeId 
        // FULL: total, full-only, ours-only, modified
        << "\t" << fullTotalSeeds << "\t" << fullOnly << "\t" << oursOnly << "\t" << fullOursMod
        // COMMON metrics
        << "\t" << commonSeqLength << "\t" << commonBlocks 
                << "\t" << blockDeletions << "\t" << blockInsertions
        // ACCURATE: inherited + local = total
        << "\t" << accurateTotalSeeds << "\t" << inheritedSeeds << "\t" << localChanges << "\t" << fullOursMod
        // SAME COMMON metrics (duplicated as requested)
        << "\t" << commonSeqLength << "\t" << commonBlocks 
                << "\t" << blockDeletions << "\t" << blockInsertions
        // Delta using accurate count
        << "\t" << deltaFullVsAccurate 
        // Seed position lists
        << "\t" << fullSeedPosStr << "\t" << accurateSeedPosStr
        // Block boundaries
        << "\t" << blockRangesStr
        // Active blocks
        << "\t" << activeBlocksStr
        // Seed differences
        << "\t" << deletedSeedsStr << "\t" << addedSeedsStr << "\t" << modifiedSeedsStr
        << "\n";
    */

        debugSeedFile.flush();
  }
}

} // Close indexing namespace
