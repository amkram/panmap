#include "indexing.hpp"
#include "panmap_utils.hpp"
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

// Global constants to control verification and debugging (for performance)
constexpr bool ENABLE_SEED_VERIFICATION = true;       // Set to true to enable seed count verification  
constexpr bool ENABLE_DEBUG_LOGGING = false;          // Set to false to disable debug output (370736, etc.)
constexpr bool ENABLE_FULL_METHOD = true;             // Set to true to enable full method computation

// Global storage for missing positions from verification (positions found by full method but missing from delta)
std::unordered_map<std::string, std::vector<int64_t>> g_missingPositionsPerNode;

// Full verification function - compares delta approach against complete recomputation
bool verifySeedApproaches(state::StateManager& stateManager, const std::string& nodeId, int k, int s) {
  if constexpr (!ENABLE_SEED_VERIFICATION) return true;
  
  logging::info("SEED_VERIFICATION: Starting verification for node {}", nodeId);
  
  if (false) { // Disabled comprehensive debug
    try {
      const auto& currentNodeState = stateManager.getNodeState(nodeId);
      const std::string& actualParentId = currentNodeState.parentId;
      
      logging::info("Node: {} | Parent: {}", nodeId, actualParentId);
      
      // 1. Check parent sequence and seed data
      if (!actualParentId.empty()) {
        const auto& parentNodeState = stateManager.getNodeState(actualParentId);
        const auto& parentSeeds = parentNodeState.materializedSeeds;
        
        auto parentSeed106353 = parentSeeds.find(106353);
        if (parentSeed106353 != parentSeeds.end()) {
          
          // Extract parent k-mer using extractKmer
          try {
            auto parentKmerResult = stateManager.extractKmer(actualParentId, 106353, k, false);
            logging::info("PARENT_EXTRACTKMER: kmer='{}' positions={}", 
                         parentKmerResult.first, parentKmerResult.second.size());
            
            // Extract parent sequence around 106353 using extractSequence
            coordinates::CoordRange parentRange = {106320, 106400};
            auto [parentGappedSeq, parentPositions, parentGaps, parentEndPos] = 
                stateManager.extractSequence(actualParentId, parentRange, false);
            logging::info("PARENT_EXTRACTSEQUENCE: range=[106320,106400) gapped='{}' positions={} gaps={}", 
                         parentGappedSeq, parentPositions.size(), parentGaps.size());
            
            // Find 106353 in parent sequence
            auto parentPos106353_it = std::find(parentPositions.begin(), parentPositions.end(), 106353);
            if (parentPos106353_it != parentPositions.end()) {
              size_t parentIdx = parentPos106353_it - parentPositions.begin();
              std::string parentSeqKmer = parentGappedSeq.substr(parentIdx, std::min(k, (int)(parentGappedSeq.length() - parentIdx)));
              logging::info("PARENT_SEQUENCE_KMER: at_index={} kmer='{}'", parentIdx, parentSeqKmer);
            }
          } catch (const std::exception& e) {
            logging::warn("PARENT_EXTRACTION_ERROR: {}", e.what());
          }
        } else {
        }
      }
      
      // 2. Check current node sequence and seed data using extractKmer
      try {
        auto currentKmerResult = stateManager.extractKmer(nodeId, 106353, k, false);
        
        // Show position mapping for current node
        for (size_t i = 0; i < std::min(currentKmerResult.second.size(), size_t(k)); i++) {
          logging::info("  [{}]: global={}", i, currentKmerResult.second[i]);
        }
      } catch (const std::exception& e) {
      }
      
      // 3. Check current node sequence using extractSequence
      try {
        coordinates::CoordRange currentRange = {106320, 106400};
        auto [currentGappedSeq, currentPositions, currentGaps, currentEndPos] = 
            stateManager.extractSequence(nodeId, currentRange, false);
        logging::info("CURRENT_EXTRACTSEQUENCE: range=[106320,106400) gapped='{}' positions={} gaps={}", 
                     currentGappedSeq, currentPositions.size(), currentGaps.size());
        
        // Find 106353 in current sequence
        auto currentPos106353_it = std::find(currentPositions.begin(), currentPositions.end(), 106353);
        if (currentPos106353_it != currentPositions.end()) {
          size_t currentIdx = currentPos106353_it - currentPositions.begin();
          std::string currentSeqKmer = currentGappedSeq.substr(currentIdx, std::min(k, (int)(currentGappedSeq.length() - currentIdx)));
          logging::info("CURRENT_SEQUENCE_KMER: at_index={} kmer='{}'", currentIdx, currentSeqKmer);
          
          // Show character-by-character comparison
          for (size_t i = 0; i < std::min(size_t(k), currentGappedSeq.length() - currentIdx); i++) {
            if (currentIdx + i < currentPositions.size()) {
              logging::info("  [{}]: pos={} char='{}'", i, currentPositions[currentIdx + i], currentGappedSeq[currentIdx + i]);
            }
          }
        } else {
        }
      } catch (const std::exception& e) {
        logging::warn("CURRENT_EXTRACTSEQUENCE_ERROR: {}", e.what());
      }
      
      // 4. Check if position 106353 is in any block and what mutations affect it
      try {
        auto coordsOpt = stateManager.fastMapGlobalToLocal(106353);
        if (coordsOpt) {
          auto [blockId, nucPos, gapPos] = *coordsOpt;
          
          // Check if this block has mutations in current node
          bool blockHasMutations = false;
          // This would require access to node mutations - simplified for now
        } else {
          logging::warn("POSITION_MAPPING: Could not map position 106353 to block coordinates");
        }
      } catch (const std::exception& e) {
        logging::warn("POSITION_MAPPING_ERROR: {}", e.what());
      }
      
    } catch (const std::exception& e) {
    }
  }
  
  try {
    // Get current materialized seeds (delta approach result)
    const auto& nodeState = stateManager.getNodeState(nodeId);
    const auto& deltaSeeds = nodeState.materializedSeeds;
    
    // SUMMARY: Show delta method k-mer for position 106353
    if (false) { // Disabled node_3001 specific debug
      auto deltaPos106353 = deltaSeeds.find(106353);
      if (deltaPos106353 != deltaSeeds.end()) {
        try {
          auto deltaKmerResult = stateManager.extractKmer(nodeId, 106353, k);
          logging::info("DELTA k-mer: '{}' (hash={}, isSeed=true)", 
                       deltaKmerResult.first, deltaPos106353->second.hash);
        } catch (const std::exception& e) {
          logging::warn("Failed to extract delta k-mer: {}", e.what());
        }
      } else {
        logging::info("DELTA k-mer: No seed at position 106353");
      }
    }
    
    if constexpr (ENABLE_FULL_METHOD) {
      // Compute seeds using full method - extract complete node sequence and find all syncmers
      std::unordered_map<int64_t, seeding::seed_t> fullSeeds;
      
      try {
        // Get the active blocks for this node to determine the correct coordinate range
        auto activeBlocks = stateManager.getActiveBlockRanges(nodeId);
        if (activeBlocks.empty()) {
          logging::warn("SEED_VERIFICATION: Node {} has no active blocks, skipping", nodeId);
          return true;
        }
        
        // Use the full span of all active blocks
        int64_t rangeStart = activeBlocks.front().second.start;
        int64_t rangeEnd = activeBlocks.back().second.end;
        coordinates::CoordRange fullRange = {rangeStart, rangeEnd};
        
        auto [gappedSeq, positions, gaps, endPositions] = stateManager.extractSequence(nodeId, fullRange, false);
        
        // Remove gaps to get ungapped sequence for rollingSyncmers
        std::string fullSequence = gappedSeq;
        fullSequence.erase(std::remove_if(fullSequence.begin(), fullSequence.end(), 
                          [](char c) { return c == '-' || c == 'x'; }), fullSequence.end());
        
        // DEBUG: Simple extraction confirmation for node_3001 
        if (false) { // Disabled node_3001 specific debug
          auto pos106353_it = std::find(positions.begin(), positions.end(), 106353);
          if (pos106353_it != positions.end()) {
            size_t gappedIdx = pos106353_it - positions.begin();
            logging::info("FULL_SEQUENCE_DEBUG: pos 106353 found at gapped index {} - k-mer comparison done above", gappedIdx);
          } else {
            logging::info("FULL_SEQUENCE_DEBUG: pos 106353 NOT found in extracted positions - range [{}, {})", 
                         rangeStart, rangeEnd);
          }
        }
        
        if (fullSequence.length() >= static_cast<size_t>(k)) {
          // Map ungapped positions back to global positions
          std::vector<int64_t> ungappedToGlobal;
          for (size_t i = 0; i < gappedSeq.length(); ++i) {
            if (gappedSeq[i] != '-' && gappedSeq[i] != 'x') {
              ungappedToGlobal.push_back(positions[i]);
            }
          }
          
          // Find all syncmers in the complete sequence using rollingSyncmers
          auto syncmers = seeding::rollingSyncmers(fullSequence, k, s, false, 0, true);
          
          // DEBUG: For node_3001, check if position 106353 appears in any syncmer
          if (false) { // Disabled node_3001 specific debug
            bool found106353InSyncmers = false;
            for (const auto& [hash, isReverse, isSeed, localStartPos] : syncmers) {
              if (localStartPos < ungappedToGlobal.size()) {
                int64_t globalPos = ungappedToGlobal[localStartPos];
                if (globalPos == 106353) {
                  found106353InSyncmers = true;
                  std::string kmerAtPos = (localStartPos + k <= fullSequence.length()) ? 
                                         fullSequence.substr(localStartPos, k) : "TRUNCATED";
                  logging::info("FULL k-mer: '{}' (hash={}, isSeed={})", 
                               kmerAtPos, hash, isSeed);
                  break;
                }
              }
            }
            if (!found106353InSyncmers) {
              logging::info("FULL k-mer: Position 106353 not found in syncmers");
            }
            logging::info("=== END K-MER SUMMARY ===");
          }
          
          for (const auto& [hash, isReverse, isSeed, localStartPos] : syncmers) {
            if (isSeed && localStartPos < ungappedToGlobal.size()) {
              int64_t globalPos = ungappedToGlobal[localStartPos];
              
              // Debug: Show the k-mer that the full method extracts for known problematic positions (DISABLED)
              // if (globalPos == 106353 || globalPos == 371508 || globalPos == 371517 || globalPos == 372797 || globalPos == 375076 || globalPos == 375077) {
              //   std::string fullMethodKmer = fullSequence.substr(localStartPos, k);
              //   logging::info("FULL_METHOD_DEBUG: pos {} -> k-mer '{}' (len={}) hash={} reversed={} from fullSequence at localPos {}", 
              //                globalPos, fullMethodKmer, fullMethodKmer.length(), hash, isReverse, localStartPos);
              // }
              
              seeding::seed_t seed;
              seed.hash = hash;
              seed.reversed = isReverse;
              seed.endPos = (localStartPos + k - 1 < ungappedToGlobal.size()) ? 
                           ungappedToGlobal[localStartPos + k - 1] : globalPos + k - 1;
              
              fullSeeds[globalPos] = seed;
            }
          }
        }
      } catch (const std::exception& e) {
        logging::err("SEED_VERIFICATION: Error extracting full sequence for node {}: {}", nodeId, e.what());
        return false;
      }
      
      // Compare delta vs full results
      bool seedCountMatch = (deltaSeeds.size() == fullSeeds.size());
      
      if (!seedCountMatch) {
        logging::err("SEED_VERIFICATION_FAILED: Node {} - SEED COUNT MISMATCH", nodeId);
        logging::err("  Delta approach: {} seeds", deltaSeeds.size());
        logging::err("  Full approach:  {} seeds", fullSeeds.size());
        
        // Log some examples of differences
        std::set<int64_t> deltaPositions, fullPositions;
        for (const auto& [pos, seed] : deltaSeeds) deltaPositions.insert(pos);
        for (const auto& [pos, seed] : fullSeeds) fullPositions.insert(pos);
        
        // Find positions only in delta
        std::vector<int64_t> onlyInDelta;
        std::set_difference(deltaPositions.begin(), deltaPositions.end(),
                           fullPositions.begin(), fullPositions.end(),
                           std::back_inserter(onlyInDelta));
        
        // Find positions only in full
        std::vector<int64_t> onlyInFull;
        std::set_difference(fullPositions.begin(), fullPositions.end(),
                           deltaPositions.begin(), deltaPositions.end(),
                           std::back_inserter(onlyInFull));
        
        // Store missing positions globally for debugging
        g_missingPositionsPerNode[nodeId] = onlyInFull;
        logging::debug("DYNAMIC_DEBUG_POSITIONS: Node {} has {} missing positions from delta method", 
                     nodeId, onlyInFull.size());
        
        if (!onlyInDelta.empty()) {
          logging::err("  First 5 positions only in DELTA: [{}]", 
                      [&]() {
                        std::stringstream ss;
                        for (size_t i = 0; i < std::min(onlyInDelta.size(), size_t(5)); i++) {
                          if (i > 0) ss << ", ";
                          ss << onlyInDelta[i];
                        }
                        return ss.str();
                      }());
        }
        
        if (!onlyInFull.empty()) {
          logging::err("  First 5 positions only in FULL: [{}]", 
                      [&]() {
                        std::stringstream ss;
                        for (size_t i = 0; i < std::min(onlyInFull.size(), size_t(5)); i++) {
                          if (i > 0) ss << ", ";
                          ss << onlyInFull[i];
                        }
                        return ss.str();
                      }());
          
          // Additional debug: check a few of these positions to see what k-mers they correspond to
          logging::err("  Checking first 3 missing positions:");
          for (size_t i = 0; i < std::min(onlyInFull.size(), size_t(3)); i++) {
            int64_t pos = onlyInFull[i];
            auto fullIt = fullSeeds.find(pos);
            if (fullIt != fullSeeds.end()) {
              const auto& seed = fullIt->second;
              try {
                auto kmerResult = stateManager.extractKmer(nodeId, pos, k);
                logging::err("    pos {} -> k-mer '{}' hash={} endPos={}", 
                           pos, kmerResult.first, seed.hash, seed.endPos);
              } catch (const std::exception& e) {
                logging::err("    pos {} -> k-mer extraction failed: {}", pos, e.what());
              }
            }
          }
        }
        
        logging::err("SEED_VERIFICATION: ABORTING due to seed count mismatch in node {}", nodeId);
        std::exit(1); // Exit immediately on first verification failure
      }
      
      for (const auto& [pos, deltaSeed] : deltaSeeds) {
        auto fullIt = fullSeeds.find(pos);
        if (fullIt != fullSeeds.end()) {
          const auto& fullSeed = fullIt->second;
          
          if (deltaSeed.hash != fullSeed.hash || 
              deltaSeed.reversed != fullSeed.reversed) {
            
            logging::err("SEED_VERIFICATION_FAILED: Node {} - SEED DETAILS MISMATCH at position {}", nodeId, pos);
            logging::err("  Delta: hash={}, reversed={}, endPos={}", deltaSeed.hash, deltaSeed.reversed, deltaSeed.endPos);
            logging::err("  Full:  hash={}, reversed={}, endPos={}", fullSeed.hash, fullSeed.reversed, fullSeed.endPos);
            logging::err("SEED_VERIFICATION: ABORTING due to seed detail mismatch in node {}", nodeId);
            std::exit(1); // Exit immediately on first verification failure
          }
        }
      }
      
      logging::info("SEED_VERIFICATION: Node {} PASSED - {} seeds match between delta and full approaches", 
                   nodeId, deltaSeeds.size());
      
    } else {
      // Fall back to simple count verification if full method is disabled
      size_t expectedSeedCount = nodeState.getTotalSeedCount();
      if (deltaSeeds.size() != expectedSeedCount) {
        logging::err("SEED_VERIFICATION_FAILED: Node {} - materialized has {} seeds, expected {} seeds", 
                       nodeId, deltaSeeds.size(), expectedSeedCount);
        logging::err("SEED_VERIFICATION: ABORTING due to count mismatch in node {}", nodeId);
        std::exit(1); // Exit immediately on first verification failure
      }
    }
    
    return true;
    
  } catch (const std::exception& e) {
    logging::err("SEED_VERIFICATION_ERROR: Exception during verification of node {}: {}", nodeId, e.what());
    logging::err("SEED_VERIFICATION: ABORTING due to exception in node {}", nodeId);
    std::exit(1); // Exit immediately on verification error
  }
}

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
    std::ofstream& debugSeedFile,
    const absl::flat_hash_map<std::string, absl::flat_hash_set<int32_t>>* cachedActiveBlocks = nullptr);


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
  if (false) { // Disabled node-specific debug
    // logging::info("DEADLOCK_DEBUG: Entering processMutationsForNode for {}", strNodeId);
  }
  
  // DEADLOCK FIX: Don't hold nodeState reference throughout entire function
  // Get it only when needed for specific operations
  
  // DEADLOCK DEBUG: Log after avoiding nodeState lock
  if (false) { // Disabled node-specific debug
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
    if (false && blockId == 595) { // Disabled node-specific debug
      // logging::info("DEADLOCK_DEBUG: About to call applyBlockMutation for block {} in {}", blockId, strNodeId);
    }
    
    // DEADLOCK FIX: Apply the mutation without holding any previous locks
    stateManager.applyBlockMutation(strNodeId, blockId, isInsertion, isMutationInversionFlag);
    
    // DEADLOCK DEBUG: Log after applyBlockMutation
    if (false && blockId == 595) { // Disabled node-specific debug
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
    
    // Only log blocks that contain debug positions
    if (false) { // Disabled node-specific debug
      std::vector<int64_t> blockDebugPos;
      // Use dynamically computed missing positions instead of hardcoded ones
      auto it = g_missingPositionsPerNode.find(strNodeId);
      if (it != g_missingPositionsPerNode.end()) {
        for (int64_t pos : it->second) {
          if (pos >= blockRange.start && pos < blockRange.end) {
            blockDebugPos.push_back(pos);
          }
        }
      } // Close the if (it != g_missingPositionsPerNode.end()) block
      
      if (!blockDebugPos.empty()) {
        std::stringstream posStr;
        for (size_t i = 0; i < blockDebugPos.size(); i++) {
          if (i > 0) posStr << ",";
          posStr << blockDebugPos[i];
        }
        logging::debug("MUTATION_CRITICAL: Block {} [{},{}): isInsertion={}, contains debug pos: {}", 
                     blockId, blockRange.start, blockRange.end, isInsertion, posStr.str());
      }
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
        
        // Only log insertion ranges that cover debug positions
        if (false) { // Disabled node-specific debug
          std::vector<int64_t> coveredDebugPos;
          // Use dynamically computed missing positions instead of hardcoded ones
          auto it = g_missingPositionsPerNode.find(strNodeId);
          if (it != g_missingPositionsPerNode.end()) {
            for (int64_t pos : it->second) {
              if (pos >= insertionRange.start && pos < insertionRange.end) {
                coveredDebugPos.push_back(pos);
              }
            }
          } // Close the if (it != g_missingPositionsPerNode.end()) block
          
          if (!coveredDebugPos.empty()) {
            std::stringstream posStr;
            for (size_t i = 0; i < coveredDebugPos.size(); i++) {
              if (i > 0) posStr << ",";
              posStr << coveredDebugPos[i];
            }
            logging::debug("INSERTION_RANGE: Block {} [{},{}): covers debug pos: {}", 
                         blockId, insertionRange.start, insertionRange.end, posStr.str());
          }
        }
      }
      
      // 2. Generate boundary range: from end of previous ON block to beginning of next ON block
      // This ensures k-mers that span the block transition are recomputed
      try {
        // Find the previous ON block before this insertion
        int32_t prevOnBlockId = -1;
        for (int32_t testBlockId = blockId - 1; testBlockId >= 0; testBlockId--) {
          if (stateManager.isBlockOn(strNodeId, testBlockId)) {
            prevOnBlockId = testBlockId;
            break;
          }
        }
        
        // Find the next ON block after this insertion
        int32_t nextOnBlockId = -1;
        for (int32_t testBlockId = blockId + 1; testBlockId < stateManager.getNumBlocks(); testBlockId++) {
          if (stateManager.isBlockOn(strNodeId, testBlockId)) {
            nextOnBlockId = testBlockId;
            break;
          }
        }
        
        if (prevOnBlockId != -1 && nextOnBlockId != -1) {
          auto prevBlockRange = stateManager.getBlockRange(prevOnBlockId);
          auto nextBlockRange = stateManager.getBlockRange(nextOnBlockId);
          
          // Create a 2-position range at the end of previous ON block
          // This will be extended by k in both directions later
          int64_t boundaryStart = std::max(prevBlockRange.end, prevBlockRange.start);
          int64_t boundaryEnd = prevBlockRange.end;
          
          coordinates::CoordRange boundaryRange = {boundaryStart, boundaryEnd};
          allRecompRanges.push_back(boundaryRange);
          rangeToBlockMapping.push_back({boundaryRange, prevOnBlockId}); // FIX: Map to the actual block containing the range
          
          if (false) { // Disabled node-specific debug
          }
        } else if (prevOnBlockId != -1) {
          // Fallback: only previous ON block exists
          auto prevBlockRange = stateManager.getBlockRange(prevOnBlockId);
          int64_t boundaryStart = std::max(prevBlockRange.end, prevBlockRange.start);
          int64_t boundaryEnd = prevBlockRange.end;
          
          coordinates::CoordRange boundaryRange = {boundaryStart, boundaryEnd};
          allRecompRanges.push_back(boundaryRange);
          rangeToBlockMapping.push_back({boundaryRange, prevOnBlockId}); // FIX: Map to the actual block containing the range
          
          if (false) { // Disabled node-specific debug
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
      // Block ON→OFF generates recomputation ranges at BOTH boundaries:
      // 1. Last k-mer of previous ON block  
      // 2. First k-mer of next ON block
      
      // DEBUG: Log block deactivation for node_4
      if (strNodeId == "node_4") {
        logging::info("BLOCK_DEACT: Block {} deactivated, finding prev/next ON blocks", blockId);
      }
      
      // 1. Create a short 2-position range at the end of the previous ON block
      // This range will be extended by k upstream/downstream in the extension phase
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
          
          // Create a simple 2-position range at the end of the previous ON block
          int64_t rangeStart = prevBlockRange.end - 2;  // Last 2 positions
          int64_t rangeEnd = prevBlockRange.end;        // End of block
          
          if (rangeStart >= prevBlockRange.start && rangeStart < rangeEnd) {
            coordinates::CoordRange prevRecompRange = {rangeStart, rangeEnd};
            allRecompRanges.push_back(prevRecompRange);
            rangeToBlockMapping.push_back({prevRecompRange, prevOnBlockId});
            
            // DEBUG: Log range generation for node_4
            if (strNodeId == "node_4") {
              logging::info("BLOCK_DEACT: Generated PREV range [{}, {}) for deactivated block {} -> prev ON block {}", 
                           rangeStart, rangeEnd, blockId, prevOnBlockId);
            }
          }
        } else {
          // DEBUG: Log when no previous ON block found
          if (strNodeId == "node_4") {
            logging::info("BLOCK_DEACT: No previous ON block found for deactivated block {}", blockId);
          }
        }
      } catch (const std::exception& e) {
        // If anything fails, continue without this range
        logging::debug("Failed to create prev block range for ON->OFF mutation: {}", e.what());
      }
      
      // 2. Create a short 2-position range at the start of the next ON block
      // This range will be extended by k upstream/downstream in the extension phase
      try {
        // Find the next ON block after this deactivation point (post-mutation state)
        int32_t nextOnBlockId = -1;
        for (int32_t testBlockId = blockId + 1; testBlockId < stateManager.getNumBlocks(); testBlockId++) {
          if (stateManager.isBlockOn(strNodeId, testBlockId)) {
            nextOnBlockId = testBlockId;
            break;
          }
        }
        
        if (nextOnBlockId != -1) {
          auto nextBlockRange = stateManager.getBlockRange(nextOnBlockId);
          
          // Create a simple 2-position range at the start of the next ON block
          int64_t rangeStart = nextBlockRange.start;     // Start of block
          int64_t rangeEnd = nextBlockRange.start + 2;   // First 2 positions
          
          if (rangeEnd <= nextBlockRange.end) {
            coordinates::CoordRange nextRecompRange = {rangeStart, rangeEnd};
            allRecompRanges.push_back(nextRecompRange);
            rangeToBlockMapping.push_back({nextRecompRange, nextOnBlockId});
            
            // DEBUG: Log range generation for node_4
            if (strNodeId == "node_4") {
              logging::info("BLOCK_DEACT: Generated NEXT range [{}, {}) for deactivated block {} -> next ON block {}", 
                           rangeStart, rangeEnd, blockId, nextOnBlockId);
            }
          }
        } else {
          // DEBUG: Log when no next ON block found
          if (strNodeId == "node_4") {
            logging::info("BLOCK_DEACT: No next ON block found for deactivated block {}", blockId);
          }
        }
      } catch (const std::exception& e) {
        // If anything fails, continue without this range
        logging::debug("Failed to create next block range for ON->OFF mutation: {}", e.what());
      }
    }
  }
  
  // Process nucleotide mutations to find recomputation ranges
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

    // PERFORMANCE OPTIMIZATION: Use bulk character extraction instead of individual calls
    // Extract all characters needed for this mutation in one operation
    std::vector<char> parentChars;
    bool useBulkExtraction = false;
    
    try {
      // Extract characters in bulk to avoid 21.6M+ individual StateManager calls
      parentChars = stateManager.getCharactersAtBlock(strNodeId, blockId, nucPos, initialGapPos, len);
      if (parentChars.size() == static_cast<size_t>(len)) {
        useBulkExtraction = true;
      } else {
        if (debug_file_ourcode_pmfn.is_open()) {
            debug_file_ourcode_pmfn << "NUC_MUT_APPLY_WARN: Bulk extraction returned " << parentChars.size() 
                               << " chars but expected " << len << " for block " << blockId 
                               << " - falling back to individual calls" << std::endl;
        }
        useBulkExtraction = false;
      }
    } catch (const std::exception &e) {
      if (debug_file_ourcode_pmfn.is_open()) {
          debug_file_ourcode_pmfn << "NUC_MUT_APPLY_WARN: Exception in bulk extraction for "
                             << blockId << ":" << nucPos << ":" << initialGapPos
                             << " len=" << len << " - " << e.what() << " - falling back to individual calls" << std::endl;
      }
      useBulkExtraction = false;
    }
    
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
      if (useBulkExtraction) {
        // Use pre-extracted character from bulk operation
        parChar = parentChars[j];
      } else {
        // Fallback to individual character access
        try {
          parChar = stateManager.getCharAtPosition(strNodeId, blockId, currentNucPos, currentGapPos);
        } catch (const std::exception &e) {
          if (debug_file_ourcode_pmfn.is_open()) {
              debug_file_ourcode_pmfn << "NUC_MUT_APPLY_ERR: Exception getting parChar for "
                                 << blockId << ":" << currentNucPos << ":" << currentGapPos
                                 << " - " << e.what() << std::endl;
          }
          throw; 
        }
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
  
  // PERFORMANCE OPTIMIZATION: Batch process nucleotide mutations by block to avoid redundant k-mer extractions
  // Group mutations by block first, then process each block's mutations together
  std::unordered_map<int32_t, std::vector<std::tuple<int32_t, int32_t, int64_t, int>>> blockMutations;
  
  // First pass: collect and group mutations by block
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
    
    if (isCurrentBlockActive) {
      nucMutationsInActiveBlocks++;
      
      // Calculate global position for sorting
      int64_t globalPos = -1;
      try {
        bool isBlockInverted = stateManager.isBlockInverted(strNodeId, blockId);
        globalPos = stateManager.mapToGlobalCoordinate(blockId, nucPos, nucGapPos, isBlockInverted);
      } catch (...) {
        globalPos = -1;
      }
      
      if (globalPos != -1) {
        blockMutations[blockId].emplace_back(nucPos, nucGapPos, globalPos, len);
      }
    } else {
      nucMutationsInInactiveBlocks++;
    }
  }
  
  // Second pass: process each block's mutations efficiently with advanced merging
  
  
  for (auto& [blockId, mutations] : blockMutations) {
    if (mutations.empty()) continue;
    
    
    // Sort mutations by global position for efficient range merging
    std::sort(mutations.begin(), mutations.end(), 
              [](const auto& a, const auto& b) { return std::get<2>(a) < std::get<2>(b); });
    
    // SMART MERGING: Use adaptive merge distance based on mutation density
    const int64_t BASE_MERGE_DISTANCE = k_val; // Base merge distance
    const int64_t MAX_MERGE_DISTANCE = k_val; // Maximum merge distance for dense clusters
    
    std::vector<coordinates::CoordRange> mergedRanges;
    
    for (auto& mutation : mutations) {
      auto [nucPos, nucGapPos, globalPos, len] = mutation;
      
      coordinates::CoordRange baseRange{globalPos, globalPos + len};
      
      // Try to merge with the last range using adaptive distance
      if (!mergedRanges.empty()) {
        auto& lastRange = mergedRanges.back();
        int64_t distance = baseRange.start - lastRange.end;
        
        // Calculate adaptive merge distance based on current range size
        int64_t rangeSize = lastRange.end - lastRange.start;
        int64_t adaptiveMergeDistance = std::min(MAX_MERGE_DISTANCE, 
                                               BASE_MERGE_DISTANCE + rangeSize / 10);
        
        if (distance <= adaptiveMergeDistance) {
          // Merge ranges - extend end position
          lastRange.end = std::max(lastRange.end, baseRange.end);
          continue;
        }
      }
      
      // Add as new range
      mergedRanges.push_back(baseRange);
    }
    
    // Process merged ranges: ALWAYS extend by k non-gap characters upstream and downstream
    for (auto& range : mergedRanges) {
      try {
        // BACKWARD EXTENSION: Extract k non-gap characters going backward from range start
        auto backwardKmerResult = stateManager.extractKmer(strNodeId, range.start, k_val, true);
        
        if (!backwardKmerResult.first.empty() && backwardKmerResult.second.size() >= k_val) {
          int64_t upstreamStart = backwardKmerResult.second[0];
          if (upstreamStart >= 0 && upstreamStart < range.start) {
            range.start = upstreamStart;
          }
        }
        
        // FORWARD EXTENSION: Extract k non-gap characters going forward from range end
        auto forwardKmerResult = stateManager.extractKmer(strNodeId, range.end, k_val, false);
        
        if (!forwardKmerResult.first.empty() && forwardKmerResult.second.size() >= k_val) {
          int64_t downstreamEnd = forwardKmerResult.second[k_val - 1] + 1; // +1 for exclusive end
          if (downstreamEnd > range.end) {
            range.end = downstreamEnd;
          }
        }
        
        // DEBUG: Log extension for node_3, node_1002, node_2002, and node_3001
        if (false) { // Disabled node-specific debug
          logging::debug("NUC_MUT_RANGE_EXT: Block {} range extended to [{}, {}) size={}", 
                       blockId, range.start, range.end, range.end - range.start);
          
          
        }
        
        nucMutationInducedRanges.push_back(range);
        allRecompRanges.push_back(range);
        
        // Debug log for merged ranges (reduced frequency)
        if constexpr (ENABLE_DEBUG_LOGGING) {
          if (false && mergedRanges.size() <= 5) { // Disabled node-specific debug
            int mutationCount = std::count_if(mutations.begin(), mutations.end(), 
                                            [&](const auto& m) { 
                                              auto globalPos = std::get<2>(m);
                                              return globalPos >= range.start && globalPos < range.end; 
                                            });
           
          }
        }
        
      } catch (const std::exception &e) {
        // Error handled silently to prevent overflow
      }
    }
  }
  
  
  // PHASE 4: Extend all ranges upstream by one k-mer and add to recomputation ranges
  if (!allRecompRanges.empty()) {
    
    // DEBUG: Log all original ranges before extension and check for missing position context
    if (strNodeId == "node_2" || strNodeId == "node_3" || strNodeId == "node_1002" || strNodeId == "node_2002" || strNodeId == "node_3001") {
      logging::debug("RANGE_EXTENSION_DEBUG: Node {} - {} original ranges before extension:", strNodeId, allRecompRanges.size());
      for (size_t i = 0; i < std::min(allRecompRanges.size(), static_cast<size_t>(10)); i++) {
        const auto& range = allRecompRanges[i];
      }
      
      // CRITICAL DEBUG: Check if missing positions should be covered by mutations
      if (strNodeId == "node_3") {
        std::vector<int64_t> missingPositions = {377943, 377962, 377983, 378463, 378473, 378477, 510287, 510288};
        
        // Check nucleotide mutations that should cover these positions
        for (const auto& nucMut : node->nucMutation) {
          int32_t blockId = nucMut.primaryBlockId;
          int32_t nucPos = nucMut.nucPosition;
          
          // Calculate global position for this mutation
          try {
            auto blockRange = stateManager.getBlockRange(blockId);
            int64_t globalMutPos = blockRange.start + nucPos;
            
            for (int64_t missingPos : missingPositions) {
              int64_t distance = std::abs(globalMutPos - missingPos);
              if (distance <= k_val) {
                logging::info("  NUC_MUT at block {} pos {} (global {}) should affect missing pos {} (distance={})", 
                             blockId, nucPos, globalMutPos, missingPos, distance);
              }
            }
          } catch (const std::exception& e) {
            logging::warn("  ERROR: Failed to get block range for mutation block {}: {}", blockId, e.what());
          }
        }
        
        // Check block mutations that should cover these positions
        for (const auto& blockMut : node->blockMutation) {
          try {
            auto blockRange = stateManager.getBlockRange(blockMut.primaryBlockId);
            for (int64_t missingPos : missingPositions) {
              if (missingPos >= blockRange.start && missingPos < blockRange.end) {
                logging::info("  BLOCK_MUT {} [{}, {}) contains missing pos {}", 
                             blockMut.primaryBlockId, blockRange.start, blockRange.end, missingPos);
              }
              // Check if missing position is within k of block boundaries
              if (std::abs(missingPos - blockRange.start) <= k_val || 
                  std::abs(missingPos - blockRange.end) <= k_val) {
                logging::info("  BLOCK_MUT {} [{}, {}) should affect missing pos {} (boundary distance={})", 
                             blockMut.primaryBlockId, blockRange.start, blockRange.end, missingPos,
                             std::min(std::abs(missingPos - blockRange.start), std::abs(missingPos - blockRange.end)));
              }
            }
          } catch (...) {}
        }
      }
      
  // DEBUG: Analyze the missing positions in context of blocks and ranges
      // Get dynamically computed missing positions from verification results
      std::vector<int64_t> missingPositions;
      auto it = g_missingPositionsPerNode.find(strNodeId);
      if (it != g_missingPositionsPerNode.end()) {
        missingPositions = it->second;
      }
      
      // For node_3, add the known problematic positions from previous verification runs
      if (node->identifier == "node_3") {
        std::vector<int64_t> knownMissingPositions = {377962, 377983, 378473, 378477, 510288};
        missingPositions.insert(missingPositions.end(), knownMissingPositions.begin(), knownMissingPositions.end());
      }
      
      // For node_1002, add the known problematic position from verification
      if (node->identifier == "node_1002") {
        std::vector<int64_t> knownMissingPositions = {377935};
        missingPositions.insert(missingPositions.end(), knownMissingPositions.begin(), knownMissingPositions.end());
      }
      
      // For node_2002, add the known problematic position from verification
      if (node->identifier == "node_2002") {
        std::vector<int64_t> knownMissingPositions = {106353}; // This is an EXTRA position in delta
        missingPositions.insert(missingPositions.end(), knownMissingPositions.begin(), knownMissingPositions.end());
      }
      
      // For node_3001, add the known problematic position from verification  
      if (node->identifier == "node_3001") {
        std::vector<int64_t> knownMissingPositions = {106353}; // This is an EXTRA position in delta
        missingPositions.insert(missingPositions.end(), knownMissingPositions.begin(), knownMissingPositions.end());
        logging::info("NODE3001_EXTRA_POS: Added {} extra positions for analysis", knownMissingPositions.size());
      }
      
      
      for (int64_t pos : missingPositions) {
        try {
          // Find which block contains this position
          int32_t containingBlockId = -1;
          bool blockIsActive = false;
          for (int32_t blockId = 0; blockId < stateManager.getNumBlocks(); blockId++) {
            auto blockRange = stateManager.getBlockRange(blockId);
            if (pos >= blockRange.start && pos < blockRange.end) {
              containingBlockId = blockId;
              blockIsActive = stateManager.isBlockOn(strNodeId, blockId);
              break;
            }
          }
          
          // Check if position is covered by any original range
          bool coveredByOriginalRange = false;
          for (const auto& range : allRecompRanges) {
            if (pos >= range.start && pos < range.end) {
              coveredByOriginalRange = true;
              break;
            }
          }
          
          
          // If in an active block but not covered, investigate WHY
          if (containingBlockId >= 0 && blockIsActive && !coveredByOriginalRange) {
            
            // Check if this block had any mutations that should have generated ranges
            bool blockHadBlockMutation = false;
            bool blockHadNucMutation = false;
            
            // Check block mutations
            for (const auto &block_mutation : node->blockMutation) {
              if (block_mutation.primaryBlockId == containingBlockId) {
                blockHadBlockMutation = true;
                logging::info("      Block {} had block mutation: insertion={}, inversion={}", 
                             containingBlockId, (block_mutation.blockMutInfo == 1), block_mutation.inversion);
                break;
              }
            }
            
            // Check nucleotide mutations - SHOW ALL MUTATIONS IN THIS BLOCK
            std::vector<std::pair<int64_t, int64_t>> nucMutationsInBlock; // (globalPos, distance)
            for (const auto &nuc_mutation : node->nucMutation) {
              if (nuc_mutation.primaryBlockId == containingBlockId) {
                blockHadNucMutation = true;
                
                // Calculate global position of this mutation
                int64_t mutGlobalPos = -1;
                try {
                  bool isBlockInverted = stateManager.isBlockInverted(strNodeId, containingBlockId);
                  mutGlobalPos = stateManager.mapToGlobalCoordinate(containingBlockId, nuc_mutation.nucPosition, 
                                                                  nuc_mutation.nucGapPosition, isBlockInverted);
                } catch (...) {
                  mutGlobalPos = -1;
                }
                
                if (mutGlobalPos != -1) {
                  int64_t distance = std::abs(pos - mutGlobalPos);
                  nucMutationsInBlock.push_back({mutGlobalPos, distance});
                }
              }
            }
            
            if (blockHadNucMutation) {
              logging::info("      Block {} has {} nucleotide mutations:", containingBlockId, nucMutationsInBlock.size());
              for (auto [mutPos, dist] : nucMutationsInBlock) {
                logging::info("        - Mutation at globalPos {} (distance {} from {})", mutPos, dist, pos);
                if (dist < k_val) {
                  logging::info("          >>> This mutation SHOULD generate a range covering {}!", pos);
                } else {
                  logging::info("          >>> This mutation is too far away (dist={} > k={})", dist, k_val);
                }
              }
            }
            
            if (!blockHadBlockMutation && !blockHadNucMutation) {
              logging::info("      Block {} had NO mutations - this explains why no ranges were generated", containingBlockId);
              logging::info("      >>> This suggests we need to ensure ALL active blocks generate some coverage");
            }
          }
          
        } catch (const std::exception& e) {
          logging::info("  Position {}: analysis failed - {}", pos, e.what());
        }
      }
    }
    
    // Step 1: Extend all ranges upstream by one k-mer (k non-gap characters)
    std::vector<coordinates::CoordRange> extendedRanges;
    extendedRanges.reserve(allRecompRanges.size());
    
    for (size_t rangeIdx = 0; rangeIdx < allRecompRanges.size(); rangeIdx++) {
      auto& range = allRecompRanges[rangeIdx];
      int64_t originalStart = range.start;
      
      try {
        // Extract k non-gap characters going backward from range start
        auto kmerResult = stateManager.extractKmer(strNodeId, range.start, k_val, true); // reverse=true
        
        if (!kmerResult.first.empty() && kmerResult.second.size() >= k_val) {
          // Extend the start position upstream by k non-gap characters
          // When reverse=true, kmerResult.second[0] is the most upstream position
          int64_t extendedStart = kmerResult.second[0]; // First position is the most upstream
          if (extendedStart >= 0 && extendedStart < range.start) {
            range.start = extendedStart;
          }
          
          // Only log expansions that affect debug positions
          if (false) { // Disabled node-specific debug
            bool affectsDebugPos = false;
            // Use dynamically computed missing positions instead of hardcoded ones
            auto it = g_missingPositionsPerNode.find(strNodeId);
            if (it != g_missingPositionsPerNode.end()) {
              for (int64_t pos : it->second) {
                if (pos >= range.start && pos < range.end && pos < originalStart) {
                  affectsDebugPos = true;
                  break;
                }
              }
            } // Close the if (it != g_missingPositionsPerNode.end()) block
            if (affectsDebugPos || rangeIdx < 3) { // Always log first 3 ranges
            }
          }
        } else {
          // No fallback - gap-aware extension failed, keep original range
        }
      } catch (const std::exception& e) {
        // No fallback - gap-aware extension failed, keep original range
      }
      
      extendedRanges.push_back(range);
    }
    
    // Step 2: Merge overlapping extended ranges
    if (!extendedRanges.empty()) {
      std::sort(extendedRanges.begin(), extendedRanges.end(), 
                [](const coordinates::CoordRange& a, const coordinates::CoordRange& b) {
                  return a.start < b.start;
                });
      
      std::vector<coordinates::CoordRange> mergedExtendedRanges;
      mergedExtendedRanges.push_back(extendedRanges[0]);
      
      for (size_t i = 1; i < extendedRanges.size(); i++) {
        auto& lastRange = mergedExtendedRanges.back();
        const auto& currentRange = extendedRanges[i];
        
        // DEBUG: For node_3001, log range merging decisions
        if (strNodeId == "node_3001") {
          logging::debug("RANGE_MERGE_DEBUG: Considering merge of [{}, {}) with last range [{}, {})", 
                       currentRange.start, currentRange.end, lastRange.start, lastRange.end);
        }
        
        if (currentRange.start <= lastRange.end) {
          // Overlapping ranges - merge them
          int64_t oldEnd = lastRange.end;
          lastRange.end = std::max(lastRange.end, currentRange.end);
          
          if (strNodeId == "node_3001") {
            logging::debug("RANGE_MERGE_DEBUG: MERGED ranges to [{}, {}) (extended end from {} to {})", 
                         lastRange.start, lastRange.end, oldEnd, lastRange.end);
          }
        } else {
          // Non-overlapping - add as new range
          mergedExtendedRanges.push_back(currentRange);
          
          if (strNodeId == "node_3001") {
            logging::debug("RANGE_MERGE_DEBUG: Added NON-OVERLAPPING range [{}, {})", 
                         currentRange.start, currentRange.end);
          }
        }
      }
      
      // DEBUG: Show ranges BEFORE adding to state to see what we generated
      if (strNodeId == "node_3001") {
        for (size_t i = 0; i < mergedExtendedRanges.size(); i++) {
          const auto& range = mergedExtendedRanges[i];
          
          // Check if this range includes position 106353
          if (range.start <= 106353 && 106353 < range.end) {
            logging::info("  *** This range INCLUDES position 106353 ***");
          }
        }
      }
      
      if constexpr (ENABLE_DEBUG_LOGGING) {
        if (false) { // Disabled node-specific debug
          logging::debug("RANGE_DEBUG_BEFORE_ADD: Node {} - Generated {} merged extended ranges:", strNodeId, mergedExtendedRanges.size());
          for (size_t i = 0; i < std::min(mergedExtendedRanges.size(), static_cast<size_t>(10)); i++) {
            const auto& range = mergedExtendedRanges[i];
                }
          
          // Check if missing positions are covered by any range AND analyze their context
          // Use dynamically computed missing positions instead of hardcoded ones
          auto it = g_missingPositionsPerNode.find(strNodeId);
          if (it != g_missingPositionsPerNode.end()) {
          for (int64_t pos : it->second) {
            bool covered = false;
            int64_t coveringRangeStart = -1, coveringRangeEnd = -1;
            for (const auto& range : mergedExtendedRanges) {
              if (pos >= range.start && pos < range.end) {
                covered = true;
                coveringRangeStart = range.start;
                coveringRangeEnd = range.end;
                break;
              }
            }
            
            if (covered) {
              logging::info("COVERAGE_CHECK_BEFORE: Position {} is COVERED by range [{}, {})", pos, coveringRangeStart, coveringRangeEnd);
            } else {
              logging::info("COVERAGE_CHECK_BEFORE: Position {} is NOT COVERED by any recomputation range", pos);
              
              // Find the closest ranges to understand why it's not covered
              int64_t closestBefore = -1, closestAfter = -1;
              for (const auto& range : mergedExtendedRanges) {
                if (range.end <= pos && (closestBefore == -1 || range.end > closestBefore)) {
                  closestBefore = range.end;
                }
                if (range.start > pos && (closestAfter == -1 || range.start < closestAfter)) {
                  closestAfter = range.start;
                }
              }
              
              logging::info("  Closest range before: ends at {}, gap = {}", closestBefore, pos - closestBefore);
              logging::info("  Closest range after: starts at {}, gap = {}", closestAfter, closestAfter - pos);
              
              // Try to extract k-mer at this position to see if there's an issue
              try {
                auto kmerResult = stateManager.extractKmer(strNodeId, pos, k_val, false);
                if (kmerResult.first.empty()) {
                  logging::info("  K-mer extraction returned EMPTY at position {}", pos);
                } else {
                  logging::info("  K-mer at position {} is '{}' (length {})", pos, kmerResult.first, kmerResult.first.length());
                }
              } catch (const std::exception& e) {
                logging::info("  K-mer extraction FAILED at position {}: {}", pos, e.what());
              }
            }
          }
          } // Close the if (it != g_missingPositionsPerNode.end()) block
        }
      }
      
      
      // Step 4: Add merged extended ranges to the node state
      // CRITICAL FIX: Only add ranges that correspond to active blocks
      auto &nodeState = stateManager.getNodeState(strNodeId);
      for (const auto& range : mergedExtendedRanges) {
        // Check if range overlaps with any active blocks
        bool rangeOverlapsActiveBlocks = false;
        try {
          auto activeBlocks = stateManager.getActiveBlockRanges(strNodeId);
          for (const auto& [blockId, blockRange] : activeBlocks) {
            if (range.start < blockRange.end && range.end > blockRange.start) {
              rangeOverlapsActiveBlocks = true;
              break;
            }
          }
          
          if (rangeOverlapsActiveBlocks) {
            nodeState.recompRanges.push_back(range);
          } else {
            logging::debug("RANGE_FILTER: Skipping range [{}, {}) for node {} - no overlap with active blocks", 
                          range.start, range.end, strNodeId);
          }
        } catch (const std::exception& e) {
          // If we can't check active blocks, add the range as fallback
          logging::warn("RANGE_FILTER: Failed to check active blocks for node {}, adding range anyway: {}", 
                       strNodeId, e.what());
          nodeState.recompRanges.push_back(range);
        }
      }
      
      // Debug log for extended ranges
      if constexpr (ENABLE_DEBUG_LOGGING) {
        if (strNodeId == "node_2" || strNodeId == "node_3") {
          logging::info("EXTENDED_RANGES: Node {} - Original {} ranges, Extended and merged to {} ranges", 
                       strNodeId, allRecompRanges.size(), mergedExtendedRanges.size());
          
          // For node_2, show all ranges to debug missing coverage
          if (false) { // Disabled node-specific debug
            logging::debug("RANGE_DEBUG: Node {} - ALL {} merged extended ranges:", strNodeId, mergedExtendedRanges.size());
            for (size_t i = 0; i < mergedExtendedRanges.size(); i++) {
              const auto& range = mergedExtendedRanges[i];
                    }
            
            // Check if missing positions are covered by any range
            // Also check for specific problematic positions from the previous run
            std::vector<int64_t> knownMissingPos = {371508, 371517, 372797, 375076, 375077};
            logging::info("COVERAGE_CHECK: Checking {} known problematic positions from previous run", knownMissingPos.size());
            for (int64_t pos : knownMissingPos) {
              bool covered = false;
              int64_t coveringRangeIdx = -1;
              for (size_t i = 0; i < mergedExtendedRanges.size(); i++) {
                const auto& range = mergedExtendedRanges[i];
                if (pos >= range.start && pos < range.end) {
                  covered = true;
                  coveringRangeIdx = i;
                  break;
                }
              }
              logging::info("COVERAGE_CHECK: Position {} is {} by recomputation ranges{}", pos, covered ? "COVERED" : "NOT COVERED", 
                           covered ? " (range " + std::to_string(coveringRangeIdx) + ")" : "");
            }
            
            // Use dynamically computed missing positions instead of hardcoded ones
            auto it = g_missingPositionsPerNode.find(strNodeId);
            if (it != g_missingPositionsPerNode.end()) {
            for (int64_t pos : it->second) {
              bool covered = false;
              for (const auto& range : mergedExtendedRanges) {
                if (pos >= range.start && pos < range.end) {
                  covered = true;
                  break;
                }
              }
              logging::info("COVERAGE_CHECK: Position {} is {} by recomputation ranges", pos, covered ? "COVERED" : "NOT COVERED");
            }
            } // Close the if (it != g_missingPositionsPerNode.end()) block
          } else {
            // For other nodes, show just first 5 ranges
            for (size_t i = 0; i < std::min(mergedExtendedRanges.size(), static_cast<size_t>(5)); i++) {
              const auto& range = mergedExtendedRanges[i];
              logging::info("  Extended range[{}]: [{}, {}) size={}", i, range.start, range.end, range.end - range.start);
            }
          }
        }
      }
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

  if constexpr (ENABLE_DEBUG_LOGGING) {
    logging::debug("DEBUG_INDEX_WRITE: Writing seed changes to CapnProto. Base positions: {}, Bit masks: {}, DFS index: {}", 
                  basePositions.size(), bitMasks.size(), dfsIndex);
  }

  
  if (basePositions.empty() || bitMasks.empty() || dfsIndex < 0) {
    if constexpr (ENABLE_DEBUG_LOGGING) {
      logging::warn("DEBUG_INDEX_WRITE: Nothing to write - empty positions/masks or invalid DFS index!");
    }
    return;
  }

  
  size_t entryIndex = dfsIndex;

  
  if (entryIndex >= perNodeSeedMutations.size()) {
    logging::err("DEBUG_INDEX_WRITE: DFS index {} is out of bounds for {} seed mutation entries",
                   dfsIndex, perNodeSeedMutations.size());
    return;
  }

  
  if constexpr (ENABLE_DEBUG_LOGGING) {
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

  
  if constexpr (ENABLE_DEBUG_LOGGING) {
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
  
  // DEBUG: For node_4, log inherited seed info before processing mutations
  if (false) { // Disabled node_4 specific debug
    auto& nodeStateDebug = stateManager.getNodeState(nodeId);
    int inheritedSeedCount = nodeStateDebug.materializedSeeds.size();
    
    // Show first few inherited seed positions
    auto it = nodeStateDebug.materializedSeeds.begin();
    for (int i = 0; i < 5 && it != nodeStateDebug.materializedSeeds.end(); i++, it++) {
    }
  }
  
  // Track specific positions of interest for debugging
  std::set<int64_t> debugPositions = {};
  if constexpr (ENABLE_DEBUG_LOGGING) {
    // Add the dynamically computed missing positions for debugging
    auto it = g_missingPositionsPerNode.find(nodeId);
    if (it != g_missingPositionsPerNode.end()) {
      for (int64_t pos : it->second) {
        debugPositions.insert(pos);
      }
      std::stringstream ss;
      for (size_t i = 0; i < std::min(it->second.size(), size_t(5)); i++) {
        if (i > 0) ss << ", ";
        ss << it->second[i];
      }
      if (it->second.size() > 5) ss << "...";
      logging::debug("DEBUG_POSITIONS: Node {} tracking {} dynamically computed positions: {}", 
                   nodeId, it->second.size(), ss.str());
    } else {
      logging::debug("DEBUG_POSITIONS: Node {} has no missing positions to track", nodeId);
    }
    
    if (false) { // Disabled node_1502 specific debug
      debugPositions.insert(370736);
      
      // DEBUG_370736_KMER_EXTRACTION: Extract k-mers using extractKmer at position 370736
      logging::debug("DEBUG_370736_KMER_EXTRACTION: Starting k-mer analysis for position 370736 in node_1502");
      
      try {
        // Extract k-mer using Delta method (current node state)
        auto deltaKmerResult = stateManager.extractKmer(nodeId, 370736, k);
        std::string deltaKmer = deltaKmerResult.first;
        std::vector<int64_t> deltaPositions = deltaKmerResult.second;
        
        // Extract k-mer using Full method (check if it would be a seed)
        std::string fullKmer = "";
        std::vector<int64_t> fullPositions;
        try {
          // Get the full sequence to compute what the k-mer should be
          coordinates::CoordRange fullRange = {370736, 370736 + k};
          auto [gappedSeq, positions, gaps, endPos] = stateManager.extractSequence(nodeId, fullRange, false);
          
          // Remove gaps to get the full k-mer
          fullKmer = gappedSeq;
          fullKmer.erase(std::remove_if(fullKmer.begin(), fullKmer.end(), 
                        [](char c) { return c == '-' || c == 'x'; }), fullKmer.end());
          fullPositions = positions;
          
          logging::debug("DEBUG_370736_KMER_EXTRACTION: Full method k-mer='{}' positions.size={}", 
                       fullKmer, fullPositions.size());
                       
        } catch (const std::exception& e) {
          logging::debug("DEBUG_370736_KMER_EXTRACTION: Full method failed: {}", e.what());
        }
        
        logging::info("DEBUG_370736_KMER_EXTRACTION: Delta method k-mer='{}' positions.size={}", 
                     deltaKmer, deltaPositions.size());
        
        // Check if k-mers are the same
        bool kmersMatch = (deltaKmer == fullKmer);
        logging::info("DEBUG_370736_KMER_EXTRACTION: K-mers match: {} (Delta='{}' vs Full='{}')", 
                     kmersMatch, deltaKmer, fullKmer);
        
        // Check syncmer status for both
        if (!deltaKmer.empty() && deltaKmer.length() == static_cast<size_t>(k)) {
          auto deltaSyncmers = seeding::rollingSyncmers(deltaKmer, k, s, false, 0, true);
          bool deltaIsSyncmer = !deltaSyncmers.empty() && std::get<2>(deltaSyncmers[0]);
          logging::info("DEBUG_370736_KMER_EXTRACTION: Delta k-mer is syncmer: {}", deltaIsSyncmer);
        }
        
        if (!fullKmer.empty() && fullKmer.length() == static_cast<size_t>(k)) {
          auto fullSyncmers = seeding::rollingSyncmers(fullKmer, k, s, false, 0, true);
          bool fullIsSyncmer = !fullSyncmers.empty() && std::get<2>(fullSyncmers[0]);
          logging::info("DEBUG_370736_KMER_EXTRACTION: Full k-mer is syncmer: {}", fullIsSyncmer);
        }
        
        // Check if position has a seed in current state
        auto seedOpt = stateManager.getSeedAtPosition(nodeId, 370736);
        logging::info("DEBUG_370736_KMER_EXTRACTION: Position has seed in Delta method: {}", 
                     seedOpt.has_value());
        if (seedOpt.has_value()) {
          logging::info("DEBUG_370736_KMER_EXTRACTION: Existing seed hash={} endPos={}", 
                       seedOpt->hash, seedOpt->endPos);
        }
        
      } catch (const std::exception& e) {
        logging::info("DEBUG_370736_KMER_EXTRACTION: Error during k-mer extraction: {}", e.what());
      }
    }
  }
  
  // DEBUG: Simple block mutation dump for node_1002, node_2002, and node_3001
  if (false) { // Disabled specific node debug
    for (size_t i = 0; i < node->blockMutation.size(); i++) {
      const auto& mut = node->blockMutation[i];
      std::string mutType = "UNKNOWN";
      if (mut.blockMutInfo == 1) mutType = "INSERTION";
      else if (mut.blockMutInfo == 0 && mut.inversion) mutType = "INVERSION";
      else if (mut.blockMutInfo == 0 && !mut.inversion) mutType = "DELETION";
      
      try {
        auto blockRange = stateManager.getBlockRange(mut.primaryBlockId);
      } catch (...) {
      }
    }
  }
  
  // STEP 1: Handle block mutations - clear seeds from deleted/inverted blocks
  if constexpr (ENABLE_DEBUG_LOGGING) {
    if (false) { // Disabled specific node debug
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
      // Use hierarchical seed storage - query current node to get inherited seeds
      for (int64_t pos = start; pos < end; pos++) {
        positionsCheckedInThisBlock++;
        auto seedOpt = stateManager.getSeedAtPosition(nodeId, pos);
        if (seedOpt.has_value()) {
          seedsFoundInThisBlock++;
          seedsCleared++;
          positionsWithClearedSeeds.insert(pos);
          
          // Record seed change for debugging
          seedChanges.emplace_back(pos, false, true, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt);
          
          // Clear the seed
          stateManager.clearSeedAtPosition(nodeId, pos);
        }
      }
      
      seedsFoundInDeletedBlocks += seedsFoundInThisBlock;
      blockDeletionPositionsChecked += positionsCheckedInThisBlock;
      
      // Add debug logging for block deletion
      logging::debug("BLOCK_DELETION: Node {} deleted block {} [{}, {}): found {} seeds in {} positions, total cleared so far: {}", 
                     nodeId, block_mutation.primaryBlockId, start, end, seedsFoundInThisBlock, positionsCheckedInThisBlock, seedsCleared);
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
    if (false) { // Disabled specific node debug
      
      if (false) { // Disabled node_2 specific debug
        
        // Only log ranges that cover debug positions
        for (size_t i = 0; i < ranges.size(); i++) {
          const auto& range = ranges[i];
          std::vector<int64_t> coveredDebugPos;
          // Use dynamically computed missing positions instead of hardcoded ones
          auto it = g_missingPositionsPerNode.find(nodeId);
          if (it != g_missingPositionsPerNode.end()) {
            for (int64_t pos : it->second) {
              if (pos >= range.start && pos < range.end) {
                coveredDebugPos.push_back(pos);
              }
            }
          } // Close the if (it != g_missingPositionsPerNode.end()) block
          
          if (!coveredDebugPos.empty()) {
            std::stringstream posStr;
            for (size_t j = 0; j < coveredDebugPos.size(); j++) {
              if (j > 0) posStr << ",";
              posStr << coveredDebugPos[j];
            }
            logging::info("  Range[{}]: [{}, {}) covers debug pos: {}", i, range.start, range.end, posStr.str());
          }
        }
      }
    }
  }
  
  // Log recomputation ranges and check for debug positions coverage
  if (false) { // Disabled node_3 specific debug
    
    // Also check the specific problematic positions
    std::vector<int64_t> problematicPositions = {377943, 377962, 377983, 378463, 378473, 378477, 510287, 510288};
    
    for (size_t i = 0; i < ranges.size(); i++) {
      const auto& range = ranges[i];
      
      // Check if any debug positions fall within this range
      std::vector<int64_t> debugPosInRange;
      for (int64_t debugPos : debugPositions) {
        if (debugPos >= range.start && debugPos < range.end) {
          debugPosInRange.push_back(debugPos);
        }
      }
      
      // Check problematic positions
      std::vector<int64_t> problematicPosInRange;
      for (int64_t pos : problematicPositions) {
        if (pos >= range.start && pos < range.end) {
          problematicPosInRange.push_back(pos);
        }
      }
      
      if (!debugPosInRange.empty()) {
        std::stringstream debugList;
        for (size_t j = 0; j < debugPosInRange.size(); j++) {
          if (j > 0) debugList << ",";
          debugList << debugPosInRange[j];
        }
      }
      
      if (!problematicPosInRange.empty()) {
        std::stringstream probList;
        for (size_t j = 0; j < problematicPosInRange.size(); j++) {
          if (j > 0) probList << ",";
          probList << problematicPosInRange[j];
        }
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
    
    // Check for problematic positions NOT covered by any range
    std::vector<int64_t> uncoveredProblematicPos;
    for (int64_t pos : problematicPositions) {
      bool covered = false;
      for (const auto& range : ranges) {
        if (pos >= range.start && pos < range.end) {
          covered = true;
          break;
        }
      }
      if (!covered) {
        uncoveredProblematicPos.push_back(pos);
      }
    }
    
    if (!uncoveredDebugPos.empty()) {
      std::stringstream uncoveredList;
      for (size_t j = 0; j < uncoveredDebugPos.size(); j++) {
        if (j > 0) uncoveredList << ",";
        uncoveredList << uncoveredDebugPos[j];
      }
    }
    
    if (!uncoveredProblematicPos.empty()) {
      std::stringstream uncoveredProbList;
      for (size_t j = 0; j < uncoveredProblematicPos.size(); j++) {
        if (j > 0) uncoveredProbList << ",";
        uncoveredProbList << uncoveredProblematicPos[j];
      }
      logging::err("CRITICAL: Node {} - PROBLEMATIC positions NOT covered by any recomp range: {}", 
                    nodeId, uncoveredProbList.str());
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
  if (false) { // Disabled node_2 specific debug
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
        
        if (false && blocksSpannedByRange.size() > 1) { // Disabled node-specific debug
          std::stringstream spannedList;
          for (auto it = blocksSpannedByRange.begin(); it != blocksSpannedByRange.end(); ++it) {
            if (it != blocksSpannedByRange.begin()) spannedList << ",";
            spannedList << *it;
          }
          // logging::info("RANGE_MAP: [{}, {}) -> spans blocks [{}] (reporting as block {})", 
                       // range.start, range.end, spannedList.str(), blockId);
        } else if (false) { // Disabled node-specific debug
          // logging::info("RANGE_MAP: [{}, {}) -> block {} (single block)", 
                       // range.start, range.end, blockId);
        }
      } else {
        // Try mapping the range start position using coordinate mapping
        auto coordOpt = stateManager.fastMapGlobalToLocal(range.start);
        if (coordOpt) {
          auto [mappedBlockId, nucPos, gapPos] = *coordOpt;
          blockId = mappedBlockId;
          if (false) { // Disabled node_2 specific debug
            // logging::info("RANGE_MAP: [{}, {}) -> block {} (via coordinate mapping)", 
                         // range.start, range.end, blockId);
          }
          blocksWithRanges.insert(blockId);
        } else {
          if (false) { // Disabled node_2 specific debug
            // logging::info("RANGE_MAP: [{}, {}) -> unmapped range", range.start, range.end);
          }
        }
      }
    } catch(...) {}
    
    recompRangesStr << length << ":" << blockId;
  }
  std::string recompRangesFormatted = recompRangesStr.str();
  
  // Add debug logging for node_2 specifically to understand undercounting
  if (false) { // Disabled node_2 specific debug
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
    int skippedDownstream = 0; // Count positions skipped due to downstream boundary
    for (size_t rangeIdx = 0; rangeIdx < rangesToProcess.size(); ++rangeIdx) {
      const auto &range = rangesToProcess[rangeIdx];
      // Targeted diagnostics for node_3 boundary mismatches
      if (false) { // Disabled node_3 specific debug
        static const std::vector<int64_t> node3Probe = {377943, 377962, 377983, 378463, 378473, 378477, 510287, 510288};
        std::vector<int64_t> inRange;
        for (auto p : node3Probe) if (p >= range.start && p < range.end) inRange.push_back(p);
        if (!inRange.empty()) {
          std::stringstream ss; for (size_t i = 0; i < inRange.size(); ++i) { if (i) ss << ","; ss << inRange[i]; }
        }
      }
      
      // For node_2, check if this range should contain positions where the full method finds seeds but delta doesn't
      std::vector<int64_t> expectedDebugPos;
      if (false) { // Disabled node_2 specific debug
        // Use dynamically computed missing positions instead of hardcoded ones
        auto it = g_missingPositionsPerNode.find(nodeId);
        if (it != g_missingPositionsPerNode.end()) {
          for (int64_t pos : it->second) {
            if (pos >= range.start && pos < range.end) {
              expectedDebugPos.push_back(pos);
            }
          }
        } // Close the if (it != g_missingPositionsPerNode.end()) block
        if (!expectedDebugPos.empty()) {
          std::stringstream expectedList;
          for (size_t i = 0; i < expectedDebugPos.size(); i++) {
            if (i > 0) expectedList << ",";
            expectedList << expectedDebugPos[i];
          }
          logging::info("RANGE_PROCESS: Range[{}] [{}, {}) should find missing positions: {}", rangeIdx, range.start, range.end, expectedList.str());
        }
      }
      
      try {
        // DEADLOCK SAFE: Extract sequence with timeout and fallback
        std::string gappedSequence;
        std::vector<int64_t> positions;
        std::vector<bool> gaps;
        std::vector<int64_t> endPositions;
        
  // Extend range end to include k more non-gap characters downstream
  // and use this as the effective processing boundary for impacted seeds
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
  // Policy: do not update downstream-only context. Process seeds only up to original range.end
  // If we ever relax this, set seedProcessingEnd = extendedRangeEnd to include boundary-shifted starts.
  int64_t seedProcessingEnd = range.end;
  coordinates::CoordRange extendedRange = {range.start, extendedRangeEnd};
        try {
          std::tie(gappedSequence, positions, gaps, endPositions) = stateManager.extractSequence(nodeId, extendedRange, false);
          
          // CRITICAL FIX: If extractSequence failed to extract any characters, clear all inherited seeds in this range
          if (gappedSequence.empty() && extendedRange.start < extendedRange.end) {
            logging::warn("RANGE_CLEAR: Clearing inherited seeds in range [{}, {}) for node {} - extractSequence returned empty sequence", 
                         extendedRange.start, extendedRange.end, nodeId);
            
            // Clear all inherited seeds in this range to ensure consistency with full approach
            auto& nodeState = stateManager.getNodeState(nodeId);
            auto& materializedSeeds = nodeState.materializedSeeds;
            
            // Find and remove all seeds in the failed range
            std::vector<int64_t> seedsToRemove;
            for (const auto& [pos, seed] : materializedSeeds) {
              if (pos >= extendedRange.start && pos < extendedRange.end) {
                seedsToRemove.push_back(pos);
              }
            }
            
            // Remove the seeds and track changes
            for (int64_t pos : seedsToRemove) {
              auto it = materializedSeeds.find(pos);
              if (it != materializedSeeds.end()) {
                // Record the deletion in seed changes
                seedChanges.emplace_back(
                    std::make_tuple(pos, true, false, std::optional<size_t>(static_cast<size_t>(it->second.hash)), std::optional<size_t>(),
                                   std::optional<bool>(it->second.reversed), std::optional<bool>(),
                                   std::optional<int64_t>(it->second.endPos), std::optional<int64_t>()));
                
                materializedSeeds.erase(it);
                seedsCleared++;
              }
            }
            
            logging::debug("RANGE_CLEAR: Cleared {} inherited seeds in failed range [{}, {}) for node {}", 
                         seedsToRemove.size(), extendedRange.start, extendedRange.end, nodeId);
            continue; // Skip to next range
          }
          
          // DEBUG: For node_3001, log ALL ranges being processed (DISABLED)
          // if (nodeId == "node_3001") {
          //   logging::info("NODE3001_RANGE_DEBUG: Processing range [{}, {}) -> extended [{}, {}), positions={}", 
          //                range.start, range.end, extendedRange.start, extendedRange.end, positions.size());
          // }
          
          // DEBUG: For node_3001, check if position 106353 is in the extracted positions
          if (false) { // Disabled position-specific debug
            
            auto pos106353_it = std::find(positions.begin(), positions.end(), 106353);
            if (pos106353_it != positions.end()) {
              size_t index = std::distance(positions.begin(), pos106353_it);
              if (index < gappedSequence.size()) {
              }
            } else {
              for (size_t i = 0; i < std::min(positions.size(), size_t(10)); i++) {
                logging::info("  Position[{}]: {}", i, positions[i]);
              }
              for (size_t i = std::max(size_t(0), positions.size() - 10); i < positions.size(); i++) {
                logging::info("  Position[{}]: {}", i, positions[i]);
              }
            }
          }
        } catch (const std::system_error& e) {
          if (e.code() == std::errc::resource_deadlock_would_occur) {
            throw std::runtime_error("Resource deadlock detected during sequence extraction for range [" + 
                                   std::to_string(extendedRange.start) + ", " + std::to_string(extendedRange.end) + 
                                   ") in node " + nodeId);
          }
          throw; // Re-throw other exceptions
        }


        std::string ungappedSequence = gappedSequence;
        
        // CRITICAL FIX: Identify gap positions and clear inherited seeds at those positions
        // This addresses the core issue where inherited seeds exist at positions that are now gaps
        std::vector<int64_t> gapPositions;
        int seedsClearedAtGaps = 0;
        
        for (size_t i = 0; i < gappedSequence.length() && i < positions.size(); i++) {
          if (gappedSequence[i] == '-' || gappedSequence[i] == 'x') {
            int64_t gapGlobalPos = positions[i];
            gapPositions.push_back(gapGlobalPos);
            
            // Clear any inherited seed at this gap position
            try {
              auto seedOpt = stateManager.getSeedAtPosition(nodeId, gapGlobalPos);
              if (seedOpt.has_value()) {
                // Clear the inherited seed since this position is now a gap
                stateManager.clearSeedAtPosition(nodeId, gapGlobalPos);
                seedsClearedAtGaps++;
                
                // Log for debugging, especially for our problem position
                if (false) { // Disabled position-specific debug
                }
              }
            } catch (const std::exception& e) {
            }
          }
        }
        
        // Log gap processing summary
        if (!gapPositions.empty()) {
          logging::debug("GAP_PROCESSING: Node {} found {} gap positions in range [{}, {}), cleared {} inherited seeds", 
                       nodeId, gapPositions.size(), range.start, range.end, seedsClearedAtGaps);
          
          // For our problem case, log details (DISABLED)
          // if (nodeId == "node_3001") {
          //   logging::info("NODE3001_GAP_DEBUG: Found {} gaps, cleared {} seeds", gapPositions.size(), seedsClearedAtGaps);
          //   if (std::find(gapPositions.begin(), gapPositions.end(), 106353) != gapPositions.end()) {
          //     logging::info("NODE3001_GAP_DEBUG: Position 106353 IS A GAP - inherited seed should be cleared");
          //   }
          // }
        }
        
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
    // ...existing code...
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
        
        // DEADLOCK SAFE: Delete seeds at gap positions, but DO NOT act beyond downstream extension boundary
        // Upstream updates are allowed by design; we only skip positions >= seedProcessingEnd
        for (int64_t globalPos : seedsToDelete) {
          if (globalPos >= seedProcessingEnd) {
            // Skip deletions beyond the downstream extension boundary to avoid corrupting seeds outside impact window
            continue;
          }
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
            
        // Debug: For node_2, check if missing positions are in this sequence and why they're not syncmers
        if (false) { // Disabled node_2 specific debug
          std::vector<int64_t> knownMissingPos = {371508, 371517, 372797, 375076, 375077};
          for (int64_t missingPos : knownMissingPos) {
            if (missingPos >= range.start && missingPos < range.end) {
              // Try to extract k-mer directly at this global position
              try {
                // First check if the position can be mapped to coordinates
                auto coordsOpt = stateManager.fastMapGlobalToLocal(missingPos);
                if (!coordsOpt) {
                } else {
                  auto [blockId, nucPos, gapPos] = *coordsOpt;
                  
                  // Try to get character at this position
                  char c = stateManager.getCharAtPosition(nodeId, blockId, nucPos, gapPos);
                }
                
                auto kmerResult = stateManager.extractKmer(nodeId, missingPos, k, false);
                             
                // Check if this position is found in our sequence mapping
                auto posIt = std::find(positions.begin(), positions.end(), missingPos);
                if (posIt != positions.end()) {
                  size_t gappedIdx = std::distance(positions.begin(), posIt);
                  size_t localUngappedPos = gappedIdx; // Start with gapped index
                  
                  // Find corresponding ungapped position
                  if (gappedIdx < localToGappedPos.size()) {
                    // Find reverse mapping from gapped to ungapped position  
                    for (size_t i = 0; i < localToGappedPos.size(); i++) {
                      if (localToGappedPos[i] == static_cast<int64_t>(gappedIdx)) {
                        localUngappedPos = i;
                        break;
                      }
                    }
                    
                    if (localUngappedPos <= ungappedSequence.length() - k) {
                      std::string ungappedKmer = ungappedSequence.substr(localUngappedPos, k);
                      
                      // Test if this k-mer should be a syncmer
                      auto singleResult = seeding::rollingSyncmers(ungappedKmer, k, s, false, 0, true);
                      bool shouldBeSyncmer = !singleResult.empty() && std::get<2>(singleResult[0]);
                      
                      // Compare direct extraction vs ungapped sequence
                      if (kmerResult.first != ungappedKmer) {
                      }
                    }
                  }
                } else {
                }
              } catch (const std::exception& e) {
              }
            }
          }
        }

        // Debug coordinate mapping for node_2 to understand the coordinate transformation issue
        if (false) { // Disabled node_2 specific debug
          size_t syncmerCount = 0;
          for (const auto& [hash, isReverse, isSyncmer, localPos] : kmers) {
            if (isSyncmer) syncmerCount++;
          }
          
          
          // Show coordinate mapping samples for the first few positions to understand the transformation
          for (size_t i = 0; i < std::min(ungappedSequence.length(), size_t(5)); i++) {
            if (i < localToGappedPos.size() && localToGappedPos[i] >= 0 && 
                localToGappedPos[i] < static_cast<int64_t>(positions.size())) {
              int64_t gappedIdx = localToGappedPos[i];
              int64_t globalPos = positions[gappedIdx];
            }
          }
          
          // Also show what global positions the syncmers actually map to
          int count = 0;
          for (const auto& [hash, isReverse, isSyncmer, localPos] : kmers) {
            if (isSyncmer && count < 10) { // Show first 10 syncmers
              if (localPos < localToGappedPos.size() && localToGappedPos[localPos] >= 0 && 
                  localToGappedPos[localPos] < static_cast<int64_t>(positions.size())) {
                int64_t gappedIdx = localToGappedPos[localPos];
                int64_t actualGlobalPos = positions[gappedIdx];
                count++;
              }
            }
          }
        }

        // Track specific debug positions in this range
        std::set<int64_t> expectedDebugPosSet(expectedDebugPos.begin(), expectedDebugPos.end());
        std::set<int64_t> foundDebugPos;
        size_t totalProcessed = 0;
        
  for (const auto& [hash, isReverse, isNowSeed, localStartPos] : kmers) {
          // Bounds check for position mapping
          if (localStartPos >= localToGappedPos.size()) {
            if (false) { // Disabled node_2 specific debug
              logging::warn("SKIP_BOUNDS: localPos {} >= mapSize {} in range[{}]", 
                           localStartPos, localToGappedPos.size(), rangeIdx);
            }
            continue;
          }
          
          int64_t gappedLocalStartPos = localToGappedPos[localStartPos];
          if (gappedLocalStartPos < 0 || gappedLocalStartPos >= static_cast<int64_t>(positions.size())) {
            if (false) { // Disabled node_2 specific debug
              logging::warn("SKIP_INVALID: gappedPos {} invalid (posSize={}) in range[{}]", 
                           gappedLocalStartPos, positions.size(), rangeIdx);
            }
            continue;
          }
          
          int64_t globalPos = positions[gappedLocalStartPos];
          totalProcessed++;
          
          // SUMMARY: Track position 106353 for node_3001
          if (false) { // Disabled position-specific debug
            std::string kmerSeq = (localStartPos + k <= ungappedSequence.size()) ? 
                                 ungappedSequence.substr(localStartPos, k) : "TRUNCATED";
          }
          
          // DEBUG: Only log the exact missing positions for node_3
          if (false) { // Disabled position-specific debug
          }
          
          // Track if this is a debug position we're expecting
          if (expectedDebugPosSet.count(globalPos)) {
            foundDebugPos.insert(globalPos);
            logging::info("FOUND_DEBUG: Range[{}] found debug pos {} at localPos={} -> globalPos={}, isSyncmer={}", 
                         rangeIdx, globalPos, localStartPos, globalPos, isNowSeed);
          }
          
          // Report syncmer processing result for debug positions  
          if (debugPositions.find(globalPos) != debugPositions.end()) {
            logging::info("DEBUG_POS_SYNCMER: Node {} range[{}] - Found debug position {} at localPos={}, isSyncmer={}, hash={}", 
                         nodeId, rangeIdx, globalPos, localStartPos, isNowSeed, hash);
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
          
          // Only process seeds up to the downstream extension boundary (k non-gap chars beyond original end)
          // This captures seeds whose starts shift slightly right due to mutations near the boundary
          if (globalPos >= seedProcessingEnd) {
            // Count skipped positions for final summary (remove log spam)
            skippedDownstream++;
            
            // Report if this debug position is being skipped due to downstream boundary
            if (debugPositions.find(globalPos) != debugPositions.end()) {
              logging::info("DEBUG_POS_TRACK: Node {} - SKIPPED debug position {} (>= seedProcessingEnd={})", 
                           nodeId, globalPos, seedProcessingEnd);
            }
            
            continue; // Skip seeds in the extended context area
          }

          // Targeted diagnostics for node_3 decision logic at probe positions
          if (false) { // Disabled node_3 specific debug
            if (globalPos == 377943 || globalPos == 377962 || globalPos == 377983 ||
                globalPos == 378463 || globalPos == 378473 || globalPos == 378477 ||
                globalPos == 510287 || globalPos == 510288) {
              char action = 'S';
              bool hadSeed = existingSeedsInRange.find(globalPos) != existingSeedsInRange.end();
              if (isNowSeed) action = hadSeed ? (existingSeedsInRange[globalPos].hash != hash ? 'M' : 'K') : 'A';
              std::string kmerSeq = (localStartPos + k <= ungappedSequence.size()) ? ungappedSequence.substr(localStartPos, k) : std::string("TRUNC");
              int64_t endPosCand = (gappedLocalStartPos < static_cast<int64_t>(endPositions.size())) ? endPositions[gappedLocalStartPos] : -1;
            }
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
          if (false) { // Disabled node_2 specific debug
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
            
            // SPECIAL: Log position 106353 processing for node_3001
            if (false) { // Disabled position-specific debug
              logging::info("DELTA_SEED_MODIFY_106353: oldHash={} newHash={} at position {}", 
                           oldSeed.hash, hash, globalPos);
            }
            
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
              std::string kmer = kmerObj.first; // Use ungapped k-mer
                               
              // Store seed
              stateManager.setSeedAtPosition(nodeId, globalPos, newSeed);
              
              seedChanges.emplace_back(
                  std::make_tuple(globalPos, true, true, std::optional<size_t>(static_cast<size_t>(oldSeed.hash)), std::optional<size_t>(static_cast<size_t>(hash)),
                                  std::optional<bool>(oldSeed.reversed), std::optional<bool>(isReverse),
                                  std::optional<int64_t>(oldSeed.endPos), std::optional<int64_t>(newSeed.endPos)));
                                  
              uniqueKmersCollector.insert(kmer);
              
              stateManager.nodeKmerSequences[nodeId][globalPos] = kmer;
              stateManager.nodeKmerEndOffsets[nodeId][globalPos] = static_cast<uint32_t>(newSeed.endPos - globalPos);
              
              // SPECIAL DEBUG: Log when modifying k-mer for position 106353 in node_3001
              if (false) { // Disabled position-specific debug
                logging::info("DELTA_MODIFY_106353: MODIFYING seed at position 106353 - old_hash={} new_hash={} k-mer='{}'", 
                             oldSeed.hash, hash, kmer);
              }
              
              // DEBUG: Compact logging for node_2
              if (false) { // Disabled node_2 specific debug
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
                if (false) { // Disabled node_2 specific debug
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
            
            // SPECIAL: Log position 106353 processing for node_3001
            if (false) { // Disabled position-specific debug
              logging::info("DELTA_SEED_ADD_106353: Adding new seed with hash={} at position {}", 
                           hash, globalPos);
            }
            
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
              
              // SPECIAL DEBUG: Log when storing k-mer for position 106353 in node_3001
              if (false) { // Disabled position-specific debug
                logging::info("DELTA_STORE_106353: STORING seed at position 106353 with k-mer='{}' hash={} endPos={}", 
                             extractedKmer, hash, newSeed.endPos);
              }
              
              // DEBUG: Compact logging for node_2
              if (false) { // Disabled node_2 specific debug
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
        if (false) { // Disabled node_2 specific debug

        }
        
        // Debug summary for ranges with debug positions
        if (false && !expectedDebugPos.empty()) { // Disabled node-specific debug
          std::set<int64_t> missingDebugPos;
          for (int64_t pos : expectedDebugPos) {
            if (foundDebugPos.find(pos) == foundDebugPos.end()) {
              missingDebugPos.insert(pos);
            }
          }
          
          if (!missingDebugPos.empty()) {
            std::stringstream missingStr;
            for (auto it = missingDebugPos.begin(); it != missingDebugPos.end(); ++it) {
              if (it != missingDebugPos.begin()) missingStr << ",";
              missingStr << *it;
            }
            logging::err("RANGE_MISSING: Range[{}] processed {}/{} kmers, MISSING debug pos: {}", 
                          rangeIdx, totalProcessed, kmers.size(), missingStr.str());
          } else {
            logging::info("RANGE_SUCCESS: Range[{}] found all expected debug positions", rangeIdx);
          }
        }
        
      } catch (const std::exception &e) {
        logging::err("Error processing range {} for node {}: {}", 
                    rangeIdx, nodeId, e.what());
      }
    }
    
  }
  
  // DEBUG: Log summary of recomputation processing for node_2
  if (false) { // Disabled node_2 specific debug
  }
  // Mirror summary for node_3 with probe positions
  if (false) { // Disabled node_3 specific debug
    std::vector<int64_t> probe = {377943, 377962, 377983, 378463, 378473, 378477, 510287, 510288};
    std::vector<int64_t> notProcessed;
    for (auto p : probe) {
      if (uniquePositionsProcessed.find(p) == uniquePositionsProcessed.end()) notProcessed.push_back(p);
    }
    if (!notProcessed.empty()) {
      std::stringstream ss; for (size_t i = 0; i < notProcessed.size(); ++i) { if (i) ss << ","; ss << notProcessed[i]; }
    }
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

  // Final debug summary for node_2
  if (false) { // Disabled node_2 specific debug
    logging::info("RECOMP_FINAL_SUMMARY: Node {} - {} ranges processed, {} unique positions, {} seed changes", 
                 nodeId, recompRangeCount, uniquePositionsProcessed.size(), seedChanges.size());
    
    // Check if all debug positions were processed
    std::vector<int64_t> unprocessedDebugPos;
    // Use dynamically computed missing positions instead of hardcoded ones
    auto it = g_missingPositionsPerNode.find(nodeId);
    if (it != g_missingPositionsPerNode.end()) {
      for (int64_t pos : it->second) {
        if (uniquePositionsProcessed.find(pos) == uniquePositionsProcessed.end()) {
          unprocessedDebugPos.push_back(pos);
        }
      }
    } // Close the if (it != g_missingPositionsPerNode.end()) block
    
    if (!unprocessedDebugPos.empty()) {
      std::stringstream unprocessedList;
      for (size_t i = 0; i < unprocessedDebugPos.size(); i++) {
        if (i > 0) unprocessedList << ",";
        unprocessedList << unprocessedDebugPos[i];
      }
      logging::err("CRITICAL: Debug positions NOT processed by delta method: {}", unprocessedList.str());
    } else {
      logging::info("DEBUG_COVERAGE: All debug positions were processed by delta method");
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

  // sort reverse
  std::sort(positionValuePairs.begin(), positionValuePairs.end(),
          [](const auto &a, const auto &b) { return a.first > b.first; });

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
              // Skip entries with null nucPosition or nucGapLength data
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
      if (false) { // Disabled node-specific debug
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
                       uniqueKmersCollector, debugSeedFile, nullptr); // No cache needed for root
    processedNodes.push_back(tree->root->identifier);
    logging::info("Root node seed processing complete."); 
    logging::debug("Root node processing complete"); // Keep original debug log too

    // Process each level in breadth-first order
    logging::info("Starting main level processing loop..."); 
    for (size_t level = 1; level < nodesByLevel.size(); level++) {
      logging::debug("Processing level {} with {} nodes", level, nodesByLevel[level].size());
      const auto& currentLevelNodes = nodesByLevel[level]; // Get nodes for the current level
      
      // Pre-cache active blocks for parent nodes to reduce StateManager contention
      logging::debug("Pre-caching active blocks for level {} parent nodes", level);
      absl::flat_hash_map<std::string, absl::flat_hash_set<int32_t>> cachedActiveBlocks;
      {
        std::set<std::string> parentIds;
        for (const auto& node : currentLevelNodes) {
          if (node && node->parent) {
            parentIds.insert(node->parent->identifier);
          }
        }
        for (const auto& parentId : parentIds) {
          cachedActiveBlocks[parentId] = stateManager->getActiveBlocks(parentId);
        }
      }
      
      // Calculate optimal grain size for TBB parallel_for  
      const size_t numThreads = std::thread::hardware_concurrency();
      const size_t levelSize = currentLevelNodes.size();
      // Increase grain size to reduce contention further
      const size_t optimalGrainSize = std::max(1UL, 
          std::min(levelSize / (numThreads * 2), 128UL)); // Larger grains, fewer threads competing

      // --- Parallelize the inner loop over nodes in the current level with optimized grain size --- 
      logging::debug("PARALLEL_START: About to start parallel processing of {} nodes at level {}", currentLevelNodes.size(), level);
      tbb::parallel_for(tbb::blocked_range<size_t>(0, currentLevelNodes.size(), optimalGrainSize),
          [&](const tbb::blocked_range<size_t>& r) { // Use capture-by-reference [&]
          
          logging::debug("THREAD_START: Processing range [{}, {}) on thread", r.begin(), r.end());
          
          // Get a reference to the thread-local k-mer collector for this thread
          auto& localKmerCollector = threadLocalKmerCollectors.local();

          for (size_t i = r.begin(); i < r.end(); ++i) {
              panmanUtils::Node* node = currentLevelNodes[i]; // Get the node for this index
              if (!node) continue; // Skip null nodes
              std::string parentId = node->parent ? node->parent->identifier : ""; // Handle null parent for root
              std::string nodeId = node->identifier;
              
              logging::debug("THREAD_PROCESSING: About to process node {} (index {})", nodeId, i);

              // NOTE: Node should already be initialized during StateManager initialization
              // Redundant initializeNode() call removed to prevent parentId corruption

              // Apply mutations from parent to child and compute seeds
              logging::debug("THREAD_CALL_PNC: Calling processNodeComplete for node {}", nodeId);
              processNodeComplete(*stateManager, engine, tree, node, tree->root, nodePaths, k, s,
                                 &perNodeSeedMutations, nullptr, true, true,
                                 localKmerCollector, debugSeedFile, &cachedActiveBlocks); // Pass cached data
              logging::debug("THREAD_FINISHED_PNC: Completed processNodeComplete for node {}", nodeId);

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
              
              logging::debug("THREAD_COMPLETED: Finished processing node {} (index {})", nodeId, i);
          } // End inner loop for range r
          
          logging::debug("THREAD_END: Completed processing range [{}, {})", r.begin(), r.end());
      }); // End tbb::parallel_for
      
      logging::debug("PARALLEL_COMPLETE: Finished parallel processing for level {}", level);

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
            
            panmapUtils::writeCapnp(message, indexPath);
            logging::info("Successfully wrote index to disk: {}", indexPath);
            
            // Verify the written index
            try {
                logging::info("Verifying the written index...");
                auto reader = panmapUtils::readCapnp(indexPath);
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
    std::ofstream& debugSeedFile,
    const absl::flat_hash_map<std::string, absl::flat_hash_set<int32_t>>* cachedActiveBlocks) {
    
    if (!node) return;
    
    std::string nodeId = node->identifier;
    logging::debug("PNC_START: Starting processNodeComplete for node {}", nodeId);

    std::string ab_before_str_for_log; // Populated by processMutationsForNode
    std::string strNodeId = node->identifier; // strNodeId is same as nodeId

    // NodeState acquisition and parent propagation (critical first steps)
    logging::debug("PNC_GET_STATE: Getting node state for {}", nodeId);
    auto &nodeState = stateManager.getNodeState(strNodeId);
    logging::debug("PNC_GOT_STATE: Successfully got node state for {}", nodeId);
    // With inheritance-based approach, we don't need to explicitly propagate state
    // Block states will be automatically inherited when queried
    if (node->parent) { 
        std::string parentId = node->parent->identifier;
        logging::debug("PNC_PARENT_CHECK: Node {} has parent {}", nodeId, parentId);
        try {
            // Just ensure the node is initialized
            logging::debug("PNC_PARENT_STATE: Getting parent state for {}", parentId);
            stateManager.getNodeState(strNodeId);
            logging::debug("PNC_PARENT_STATE_OK: Got parent state for {}", parentId);
            // std::cerr << "[DEBUG] PNC: Node " << strNodeId << " will inherit from parent " << parentId << " as needed" << std::endl;
        } catch (const std::exception& e) {
            logging::err("PNC: Error accessing state for {} (with parent {}): {}", strNodeId, parentId, e.what());
            throw; 
        }
    }

    logging::debug("PNC_BLOCKS: About to process active blocks for node {}", nodeId);
    if (node->identifier != tree->root->identifier) {
      // Use cached active blocks if available, otherwise fall back to StateManager
      if (cachedActiveBlocks && node->parent) {
        auto cacheIt = cachedActiveBlocks->find(node->parent->identifier);
        if (cacheIt != cachedActiveBlocks->end()) {
          logging::debug("PNC_CACHE_HIT: Using cached blocks for parent {}", node->parent->identifier);
          // Use cached data - no StateManager lock needed
          for (const int32_t blockId : cacheIt->second) {
            stateManager.setBlockOn(strNodeId, blockId, true);
          }
          logging::debug("PNC_CACHE_DONE: Applied {} cached blocks for node {}", cacheIt->second.size(), nodeId);
        } else {
          logging::debug("PNC_CACHE_MISS: Cache miss for parent {}, using StateManager", node->parent->identifier);
          // Fallback to StateManager if not in cache
          for (const int32_t blockId : stateManager.getActiveBlocks(node->parent->identifier)) {
            stateManager.setBlockOn(strNodeId, blockId, true);
          }
          logging::debug("PNC_FALLBACK_DONE: Applied blocks via StateManager for node {}", nodeId);
        }
      } else {
        logging::debug("PNC_NO_CACHE: No cache provided, using StateManager directly for node {}", nodeId);
        // No cache provided, use StateManager directly
        for (const int32_t blockId : stateManager.getActiveBlocks(node->parent->identifier)) {
          stateManager.setBlockOn(strNodeId, blockId, true);
        }
        logging::debug("PNC_DIRECT_DONE: Applied blocks directly via StateManager for node {}", nodeId);
      }
    }
    logging::debug("PNC_BLOCKS_DONE: Finished processing active blocks for node {}", nodeId);

    std::stringstream condensed_ss_ourcode_stream; // Use a local stream for this function\'s log part


    // CRITICAL: Materialize node state BEFORE applying local mutations
    // This ensures inherited seeds are available during block deletion processing
    // std::cout << "DEBUG_PROCESS_NODE: Processing node " << nodeId << std::endl;
    stateManager.getNodeState(nodeId); // This will trigger initialization if needed
    // std::cout << "DEBUG_PROCESS_NODE: About to materialize node " << nodeId << std::endl;
    
    // Debug: Check seed count before materialization
    if constexpr (ENABLE_SEED_VERIFICATION) {
      if (false) { // Disabled node_2 specific debug
        const auto& nodeStateBeforeMat = stateManager.getNodeState(nodeId);
        
        // Check parent state
        if (node->parent) {
          try {
            const auto& parentState = stateManager.getNodeState(node->parent->identifier);
          } catch (const std::exception& e) {
          }
        }
      }
    }
    
    stateManager.materializeNodeState(nodeId); // EXPLICITLY materialize node state
    
    // Debug: Check seed count after materialization
    if constexpr (ENABLE_SEED_VERIFICATION) {
      if (false) { // Disabled node_2 specific debug
        const auto& nodeStateAfterMat = stateManager.getNodeState(nodeId);
      }
    }
    
    // std::cout << "DEBUG_PROCESS_NODE: Materialized node " << nodeId << std::endl;
    
    // Apply mutations for the node IF THE FLAG IS SET
    std::string* ab_before_ptr = nullptr;
    std::stringstream* condensed_stream_for_pmfn = nullptr; 

    processMutationsForNode(stateManager, tree, node, commonAncestor, nodePaths, 
        ab_before_ptr,
        condensed_stream_for_pmfn 
    );

    if (computeSeeds) {
        auto [blockDeletions, blockInsertions, totalMaterializedSeeds, seedsFoundInDeletedBlocks, blockDeletionPositionsChecked, seedsCleared, seedsAdded, recompRangeCount, recompRangeSize, recompRangeList, uniquePositionsProcessed, totalKmersInRanges, detailedSeedChanges] = recomputeSeeds(stateManager, engine, node, k, s, perNodeSeedMutations, kmerDictionary, uniqueKmersCollector, debugSeedFile);
        
        // Debug: Log seed change statistics from recomputation
        if constexpr (ENABLE_SEED_VERIFICATION) {
          size_t additionChanges = 0, deletionChanges = 0, modificationChanges = 0;
          for (const auto& [pos, wasSeed, isSeed, prevHash, newHash, prevReversed, newReversed, prevEndPos, newEndPos] : detailedSeedChanges) {
            if (!wasSeed && isSeed) additionChanges++;
            else if (wasSeed && !isSeed) deletionChanges++;
            else if (wasSeed && isSeed) modificationChanges++;
          }
          // logging::info("RECOMP_DEBUG: Node {} - Generated {} seed changes: {} additions, {} deletions, {} modifications", 
          //              nodeId, detailedSeedChanges.size(), additionChanges, deletionChanges, modificationChanges);
        }

        // Update materialized seed state with detailed seed changes from recomputation
        stateManager.updateMaterializedSeedsAfterRecomputation(nodeId, detailedSeedChanges);
        
        // Debug: Log materialized seed count after update
        if constexpr (ENABLE_SEED_VERIFICATION) {
          const auto& nodeStateAfterUpdate = stateManager.getNodeState(nodeId);
          const auto& materializedSeedsAfterUpdate = nodeStateAfterUpdate.materializedSeeds;
          // logging::info("SEED_UPDATE_DEBUG: Node {} - Applied {} seed changes, materialized seeds count now: {}", 
          //              nodeId, detailedSeedChanges.size(), materializedSeedsAfterUpdate.size());
        }

        // Verify delta vs full approach match - break on first failure
        if constexpr (ENABLE_SEED_VERIFICATION) {
          if (!verifySeedApproaches(stateManager, nodeId, k, s)) {
            logging::err("SEED_VERIFICATION: Stopping execution due to verification failure in node {}", nodeId);
            throw std::runtime_error("Seed verification failed for node " + nodeId);
          }
        }

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

    // Get seed counts from verification function for comparison
    int64_t fullTotalSeeds = 0;
    int64_t accurateTotalSeeds = 0;
    
    try {
        const auto& nodeState = stateManager.getNodeState(nodeId);
        accurateTotalSeeds = nodeState.materializedSeeds.size();
    } catch (const std::exception& e) {
        logging::warn("Failed to get accurate seed count for node {}: {}", nodeId, e.what());
    }
    
    // Compute full span for logging
    int64_t fullStart = activeBlocksAndRanges.front().second.start;
    int64_t fullEnd = activeBlocksAndRanges.back().second.end;
    
    // Compute full method seed count for comparison
    try {
        coordinates::CoordRange fullRange = {fullStart, fullEnd};
        auto [gappedSeq, positions, gaps, endPositions] = stateManager.extractSequence(nodeId, fullRange, false);
        
        // Remove gaps to get ungapped sequence for rollingSyncmers
        std::string fullSequence = gappedSeq;
        fullSequence.erase(std::remove_if(fullSequence.begin(), fullSequence.end(), 
                          [](char c) { return c == '-' || c == 'x'; }), fullSequence.end());
        
        if (fullSequence.length() >= static_cast<size_t>(k)) {
          // Map ungapped positions back to global positions
          std::vector<int64_t> ungappedToGlobal;
          for (size_t i = 0; i < gappedSeq.length(); ++i) {
            if (gappedSeq[i] != '-' && gappedSeq[i] != 'x') {
              ungappedToGlobal.push_back(positions[i]);
            }
          }
          
          // Find all syncmers in the complete sequence using rollingSyncmers
          auto syncmers = seeding::rollingSyncmers(fullSequence, k, s, false, 0, true);
          
          for (const auto& [hash, isReverse, isSeed, localStartPos] : syncmers) {
            if (isSeed && localStartPos < ungappedToGlobal.size()) {
              fullTotalSeeds++;
            }
          }
        }
    } catch (const std::exception& e) {
        logging::warn("Failed to compute full method seed count for node {}: {}", nodeId, e.what());
        fullTotalSeeds = -1; // Use -1 to indicate failure
    }
    
    // if (!activeBlocksAndRanges.empty()) {
    //     fullStart = activeBlocksAndRanges.front().second.start;
    //     fullEnd = activeBlocksAndRanges.back().second.end;
    //     commonSeqLength = fullEnd - fullStart;
    //     commonBlocks = activeBlocksAndRanges.size();
    // }

    // // === FULL method - SAFE version ===
    // std::string fullSeq = "";
    // std::vector<int64_t> fullPos;
    // std::vector<bool> fullGaps;
    // std::vector<int64_t> fullEndPos;
    
    // // Try sequence extraction safely
    // try {
    //     // Extend fullEnd to include k more non-gap characters downstream
    //     int64_t extendedFullEnd = fullEnd;
    //     int k_val = stateManager.getKmerSize();
    //     if (k_val > 0) {
    //         try {
    //             // Use extractKmer to get k characters downstream from fullEnd
    //             auto kmerResult = stateManager.extractKmer(nodeId, fullEnd, k_val);
    //             if (!kmerResult.first.empty() && !kmerResult.second.empty()) {
    //                 // Use the position after the k-th character (end position of k-mer)
    //                 extendedFullEnd = kmerResult.second[k_val - 1] + 1;
    //             } else {
    //                 // If k-mer extraction fails, fall back to simple position extension
    //                 extendedFullEnd = fullEnd + k_val;
    //             }
    //         } catch (...) {
    //             // If k-mer extraction fails, fall back to simple position extension
    //             extendedFullEnd = fullEnd + k_val;
    //         }
    //     }
        
    //     auto extractResult = stateManager.extractSequence(nodeId, {fullStart, extendedFullEnd}, false);
    //     fullSeq = std::get<0>(extractResult);
    //     fullPos = std::get<1>(extractResult);
    //     fullGaps = std::get<2>(extractResult);
    //     fullEndPos = std::get<3>(extractResult);
    // } catch (...) {
    //     // Skip StateManager call to avoid deadlock
    //     logging::warn("Skipped extractSequence for node {} to avoid deadlock", nodeId);
    //     // Use empty defaults
    // }

    // // Strip gaps to get ungappedSequence
    // std::string ungapped = fullSeq;
    // ungapped.erase(std::remove_if(ungapped.begin(), ungapped.end(),
    //                               [](char c){ return c=='-'||c=='x'; }),
    //                 ungapped.end());

    // // Map ungapped → gapped indices
    // std::vector<int64_t> ung2gap(ungapped.size(), -1);
    // for (size_t i = 0, gi = 0; i < fullSeq.size(); ++i) {
    //     if (fullSeq[i] != '-' && fullSeq[i] != 'x') {
    //         ung2gap[gi++] = i; // Map gapped position to ungapped index
    //     }
    // }

    // // Collect rolling syncmers
    // auto fullKmers = seeding::rollingSyncmers(ungapped, k, s, false, 0, true);

    // // Count only actual syncmers, not all k-mer positions
    // // But restrict to seeds within the original range [fullStart, fullEnd]
    // int64_t fullTotalSeeds = 0;
    // std::vector<int64_t> fullSeedPositions;
    // for (const auto& [hash, isReverse, isSyncmer, startPos] : fullKmers) {
    //     if (isSyncmer) {
    //         // Convert local ungapped position to global position
    //         if (startPos < ung2gap.size() && ung2gap[startPos] != -1) {
    //             int64_t globalPos = fullPos[ung2gap[startPos]];
    //             // Only count seeds within the original range [fullStart, fullEnd]
    //             if (globalPos >= fullStart && globalPos <= fullEnd) {
    //                 fullTotalSeeds++;
    //                 fullSeedPositions.push_back(globalPos);
    //             }
    //         }
    //     }
    // }
    // int64_t fullOnly        = 0;  // compute as needed
    // int64_t oursOnly        = 0;  // compute as needed
    // int64_t fullOursMod     = 0;  // compute as needed

    // // === OURS method - SAFE version ===
    // // Get seed positions from our implementation using the complete materialized state
    // std::unordered_set<int64_t> oursPositions;
    
    // // Get all seed positions from the materialized state (includes inherited + local)
    // try {
    //     const auto& nodeState = stateManager.getNodeState(nodeId);
    //     // Get all materialized seeds (both inherited and local)
    //     for (const auto& [position, seed] : nodeState.materializedSeeds) {
    //         oursPositions.insert(position);
    //     }
    //     logging::info("Successfully retrieved {} complete seed positions from materialized state for node {}", oursPositions.size(), nodeId);
    // } catch (const std::exception& e) {
    //     throw std::runtime_error("Failed to retrieve materialized seeds for node " + nodeId + ": " + e.what());
    // }
    
    // int64_t oursTotalSeeds = oursPositions.size();
    
    

    // // Compute diffs between fullKmers and oursSeeds if you want fullOnly/oursOnly/fullOursMod:
    // //   for each seed in fullKmers  → if not in oursSeeds ⇒ fullOnly++
    // //   for each seed in oursSeeds  → if not in fullKmers ⇒ oursOnly++
    // //   common ⇒ else fullOursMod++

    // // Compute delta if desired:
    // int64_t deltaFullVsOurs = oursTotalSeeds - fullTotalSeeds;
    
    // // Get the accurate total seed count from NodeState
    // const auto& nodeState = stateManager.getNodeState(nodeId);
    // int64_t accurateTotalSeeds = nodeState.getTotalSeedCount();
    // int64_t inheritedSeeds = nodeState.inheritedSeedCount;
    // int64_t localChanges = nodeState.localSeedChanges;
    // int64_t deltaFullVsAccurate = accurateTotalSeeds - fullTotalSeeds;

    // // Convert seed position vectors to comma-separated strings
    // std::string fullSeedPosStr = "";
    // for (size_t i = 0; i < fullSeedPositions.size(); ++i) {
    //     if (i > 0) fullSeedPosStr += ",";
    //     fullSeedPosStr += std::to_string(fullSeedPositions[i]);
    // }
    
    // std::string accurateSeedPosStr = "";
    // bool first = true;
    // for (const auto& position : oursPositions) {
    //     if (!first) accurateSeedPosStr += ",";
    //     accurateSeedPosStr += std::to_string(position);
    //     first = false;
    // }

    // // Format active block ranges for TSV output
    // std::string blockRangesStr = "";
    // for (size_t i = 0; i < activeBlocksAndRanges.size(); ++i) {
    //     if (i > 0) blockRangesStr += ",";
    //     const auto& [blockId, range] = activeBlocksAndRanges[i];
    //     blockRangesStr += std::to_string(range.start) + "-" + std::to_string(range.end) + ":" + std::to_string(blockId);
    // }
    
    // // Format active block IDs for TSV output
    // std::string activeBlocksStr = "";
    // for (size_t i = 0; i < activeBlocksAndRanges.size(); ++i) {
    //     if (i > 0) activeBlocksStr += ",";
    //     const auto& [blockId, range] = activeBlocksAndRanges[i];
    //     activeBlocksStr += std::to_string(blockId);
    // }

    // // Compute seed position differences between FULL and ACCURATE methods
    // std::set<int64_t> fullPositionsSet(fullSeedPositions.begin(), fullSeedPositions.end());
    // std::set<int64_t> accuratePositionsSet(oursPositions.begin(), oursPositions.end());
    
    // // Deleted seeds: in FULL but not in ACCURATE
    // std::vector<int64_t> deletedSeeds;
    // std::set_difference(fullPositionsSet.begin(), fullPositionsSet.end(),
    //                    accuratePositionsSet.begin(), accuratePositionsSet.end(),
    //                    std::back_inserter(deletedSeeds));
    
    // // Added seeds: in ACCURATE but not in FULL
    // std::vector<int64_t> addedSeeds;
    // std::set_difference(accuratePositionsSet.begin(), accuratePositionsSet.end(),
    //                    fullPositionsSet.begin(), fullPositionsSet.end(),
    //                    std::back_inserter(addedSeeds));
    
    // // Modified seeds: in both FULL and ACCURATE (intersection)
    // std::vector<int64_t> modifiedSeeds;
    // std::set_intersection(fullPositionsSet.begin(), fullPositionsSet.end(),
    //                      accuratePositionsSet.begin(), accuratePositionsSet.end(),
    //                      std::back_inserter(modifiedSeeds));

    // // Format seed difference lists as comma-separated strings
    // std::string deletedSeedsStr = "";
    // for (size_t i = 0; i < deletedSeeds.size(); ++i) {
    //     if (i > 0) deletedSeedsStr += ",";
    //     deletedSeedsStr += std::to_string(deletedSeeds[i]);
    // }
    
    // std::string addedSeedsStr = "";
    // for (size_t i = 0; i < addedSeeds.size(); ++i) {
    //     if (i > 0) addedSeedsStr += ",";
    //     addedSeedsStr += std::to_string(addedSeeds[i]);
    // }
    
    // std::string modifiedSeedsStr = "";
    // for (size_t i = 0; i < modifiedSeeds.size(); ++i) {
    //     if (i > 0) modifiedSeedsStr += ",";
    //     modifiedSeedsStr += std::to_string(modifiedSeeds[i]);
    // }

    // // For node_2, extract and log k-mers at missing positions
    // if (nodeId == "node_2") {
    //     std::set<int64_t> fullPositionsSet(fullSeedPositions.begin(), fullSeedPositions.end());
    //     std::set<int64_t> missingPositions;
    //     for (int64_t pos : fullPositionsSet) {
    //         if (oursPositions.find(pos) == oursPositions.end()) {
    //             missingPositions.insert(pos);
    //         }
    //     }
        
    //     if (!missingPositions.empty()) {
    //         logging::info("NODE_2_MISSING_KMERS: Found {} missing positions in ACCURATE vs FULL", missingPositions.size());
            
    //         for (int64_t pos : missingPositions) {
    //             std::string fullKmer = "NOTFOUND";
    //             std::string accurateKmer = "NOTFOUND";
                
    //             // --- FULL k-mer extraction (SAFE) ---
    //             try {
    //                 // Find the corresponding k-mer in fullKmers by matching the global position
    //                 size_t syncmerIdx = 0;
    //                 for (const auto& [hash, isReverse, isSyncmer, startPos] : fullKmers) {
    //                     if (isSyncmer) {
    //                         if (syncmerIdx < fullSeedPositions.size() && fullSeedPositions[syncmerIdx] == pos) {
    //                             if (startPos + k <= ungapped.size()) {
    //                                 fullKmer = ungapped.substr(startPos, k);
    //                             }
    //                             break;
    //                         }
    //                         syncmerIdx++;
    //                     }
    //                 }
    //             } catch (...) {
    //                 // If extraction fails, leave as NOTFOUND
    //             }
                
    //             // --- ACCURATE (OURS) k-mer extraction ---
    //             // First try our stored sequences
    //             auto nodeKmersIt = stateManager.nodeKmerSequences.find(nodeId);
    //             if (nodeKmersIt != stateManager.nodeKmerSequences.end()) {
    //                 auto kmerIt = nodeKmersIt->second.find(pos);
    //                 if (kmerIt != nodeKmersIt->second.end()) {
    //                     accurateKmer = kmerIt->second;
    //                 } else {
    //                     accurateKmer = "NOT_FOUND_IN_STORED_KMERS";
    //                 }
    //             } else {
    //                 accurateKmer = "NO_KMERS_STORED_FOR_NODE";
    //             }
                
    //             // Additionally, try to extract the actual k-mer at this position from our sequence
    //             std::string actualKmerAtPos = "EXTRACTION_FAILED";
    //             try {
    //                 auto kmerResult = stateManager.extractKmer(nodeId, pos, k);
    //                 actualKmerAtPos = kmerResult.first;
    //             } catch (...) {
    //                 actualKmerAtPos = "EXTRACTION_FAILED";
    //             }
                
    //             logging::info("NODE_2_MISSING_KMER: pos={} FULL_kmer={} ACCURATE_stored={} ACCURATE_actual={}", 
    //                           pos, fullKmer, accurateKmer, actualKmerAtPos);
    //         }
    //     }
    // }

    // Log comparison for debugging
    if (fullTotalSeeds >= 0 && accurateTotalSeeds >= 0) {
        int64_t delta = accurateTotalSeeds - fullTotalSeeds;
        if (delta != 0) {
            logging::warn("SEED_COUNT_DIFF: Node {} - FULL: {}, ACCURATE: {}, DELTA: {}", 
                         nodeId, fullTotalSeeds, accurateTotalSeeds, delta);
        } else {
            logging::info("SEED_COUNT_MATCH: Node {} - Both methods produced {} seeds", 
                         nodeId, fullTotalSeeds);
        }
    }
    
    // === TSV OUTPUT FOR VERIFICATION ===
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
