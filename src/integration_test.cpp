#include "coordinates.hpp"
#include "indexing.hpp"
#include "placement.hpp"
#include "seed_annotated_tree.hpp"
#include "logging.hpp"
#include "gap_map.hpp"
#include "index.capnp.h"
#include <capnp/message.h>
#include <capnp/serialize-packed.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <stdexcept>
#include <string>
#include <vector>
#include <filesystem>

using namespace std;
using namespace logging;
using namespace coordinates;
using namespace placement;
using namespace seed_annotated_tree;
using namespace indexing;
namespace fs = std::filesystem;

// Forward declaration of loadPanMAN from main.cpp
panmanUtils::TreeGroup* loadPanMAN(const std::string& panman_file);

// Test function to deliberately introduce invalid gap ranges
void introduceInvalidGapRanges(CoordinateManager& manager) {
    msg("Deliberately introducing invalid gap ranges for testing...");
    
    // Get the gap map for modification
    auto& gap_map = const_cast<coordinates::GapMap&>(manager.getGapMap());
    int64_t totalCoords = manager.size();
    
    // Add some invalid gap ranges:
    
    // 1. Negative position
    gap_map[-5] = 10;
    msg("Added invalid gap with negative position: -5");
    
    // 2. Zero length
    gap_map[100] = 0;
    msg("Added invalid gap with zero length at position 100");
    
    // 3. Negative length
    gap_map[200] = -10;
    msg("Added invalid gap with negative length at position 200");
    
    // 4. Out of bounds (if we have bounds)
    if (totalCoords > 0) {
        gap_map[totalCoords - 10] = 20;
        msg("Added invalid gap extending beyond bounds: [{}, {}]", 
            totalCoords - 10, totalCoords + 10);
    }
    
    // 5. Overlapping gaps
    gap_map[50] = 20;
    gap_map[60] = 20;
    msg("Added overlapping gaps: [50, 69] and [60, 79]");
    
    msg("Successfully introduced invalid gap ranges");
}

// Test function to validate gap ranges
bool testGapRangeValidation(panmanUtils::Tree* tree) {
    msg("Testing gap range validation...");
    
    try {
        // Create test data structures
        sequence_t sequence;
        blockExists_t blockExists;
        blockStrand_t blockStrand;
        
        // Setup the coordinate system
        setupIndexing(sequence, blockExists, blockStrand, tree);
        
        // Verify sequence was created properly
        if (sequence.empty()) {
            throw std::runtime_error("Sequence is empty after setup");
        }
        
        // Create coordinate manager
        CoordinateManager coordManager(&sequence, &blockExists, &blockStrand);
        
        // Initialize gap map
        CoordinateTraverser traverser(coordManager.getLeftmostCoordinate(), &coordManager);
        initializeGapMapFromSequence(traverser, coordManager);
        
        // Verify gap map has entries
        int initial_gap_count = coordManager.getGapMap().size();
        msg("Gap map initialized with {} entries", initial_gap_count);
        if (initial_gap_count == 0) {
            throw std::runtime_error("No gaps found in sequence");
        }
        
        // First validate the initial gap map
        bool initialValid = seed_annotated_tree::validateAndFixGapMap(coordManager, "initial validation");
        msg("Initial gap map valid: {}", initialValid ? "YES" : "NO");
        
        // Introduce invalid gap ranges for testing
        introduceInvalidGapRanges(coordManager);
        
        // Count gaps after introducing invalid entries
        int after_invalid_count = coordManager.getGapMap().size();
        msg("Gap map now has {} entries after introducing invalid gaps", after_invalid_count);
        
        // Run the validation and fixing function
        msg("Running gap map validation to fix issues...");
        bool afterFixValid = seed_annotated_tree::validateAndFixGapMap(coordManager, "after introducing invalid gaps");
        
        // Count gaps after fixing
        int after_fix_count = coordManager.getGapMap().size();
        msg("Gap map now has {} entries after validation and fixing", after_fix_count);
        
        // Verify all gaps are now valid
        bool allValid = true;
        const auto& gap_map = coordManager.getGapMap();
        int64_t totalCoords = coordManager.size();
        
        for (const auto& [pos, length] : gap_map) {
            if (pos < 0 || length <= 0 || pos + length > totalCoords) {
                msg("Found invalid gap after fixing: position={}, length={}", pos, length);
                allValid = false;
            }
        }
        
        // Check for overlaps
        auto it = gap_map.begin();
        while (it != gap_map.end()) {
            auto next_it = std::next(it);
            if (next_it != gap_map.end()) {
                int64_t end = it->first + it->second - 1;
                if (end >= next_it->first) {
                    msg("Found overlapping gaps after fixing: [{},{}] and [{},{}]", 
                        it->first, end, next_it->first, next_it->first + next_it->second - 1);
                    allValid = false;
                }
            }
            it = next_it;
        }
        
        msg("Gap map validation test result: {}", allValid ? "PASSED" : "FAILED");
        return allValid;
    }
    catch (const std::exception& e) {
        err("Gap range validation test failed: {}", e.what());
        return false;
    }
}

// Function to create a minimal FASTQ file for testing
bool createTestFastq(const std::string& filename, panmanUtils::Tree* tree) {
    msg("Creating test FASTQ file {}...", filename);
    
    try {
        // Get sequence from a random node in the tree
        std::vector<std::string> nodeIds;
        for (const auto& [id, node] : tree->allNodes) {
            nodeIds.push_back(id);
        }
        
        if (nodeIds.empty()) {
            throw std::runtime_error("No nodes found in tree");
        }
        
        // Pick a random node
        const std::string& randomNodeId = nodeIds[0];
        std::string sequence = tree->getStringFromReference(randomNodeId, false, true);
        
        if (sequence.empty()) {
            throw std::runtime_error("Failed to extract sequence from node");
        }
        
        // Create a few reads from the sequence
        std::ofstream outfile(filename);
        if (!outfile.is_open()) {
            throw std::runtime_error("Failed to create test FASTQ file");
        }
        
        // Write a few reads from different parts of the sequence
        const int READ_LENGTH = 150;
        const int NUM_READS = 5;
        
        for (int i = 0; i < NUM_READS; i++) {
            // Position to extract read (evenly distributed across sequence)
            size_t pos = 0;
            if (sequence.length() > READ_LENGTH) {
                pos = (i * (sequence.length() - READ_LENGTH)) / (NUM_READS - 1);
            }
            
            // Extract read
            std::string read = sequence.substr(pos, READ_LENGTH);
            if (read.empty()) continue;
            
            // Write FASTQ format
            outfile << "@test_read_" << i << " from:" << randomNodeId << " pos:" << pos << "\n";
            outfile << read << "\n";
            outfile << "+\n";
            // Add quality scores (all high quality 'I')
            outfile << std::string(read.length(), 'I') << "\n";
        }
        
        outfile.close();
        
        msg("Created test FASTQ with {} reads from node {}", NUM_READS, randomNodeId);
        return true;
    }
    catch (const std::exception& e) {
        err("Failed to create test FASTQ file: {}", e.what());
        return false;
    }
}

int main(int argc, char* argv[]) {
    // Set up paths
    std::string panman_file = "dev/panmans/rsv_4000.panman";
    std::string test_fastq = "test_reads.fastq";
    
    // Allow command-line override of panman file
    if (argc > 1) {
        panman_file = argv[1];
    }
    
    msg("Starting integration tests with panman file: {}", panman_file);
    
    // Load the panman file
    panmanUtils::TreeGroup* treeGroup = loadPanMAN(panman_file);
    if (!treeGroup) {
        err("Failed to load panman file");
        return 1;
    }
    
    msg("Successfully loaded pangenome with {} trees", treeGroup->trees.size());
    
    if (treeGroup->trees.empty()) {
        err("No trees found in the pangenome");
        return 1;
    }
    
    // Use the first tree for testing
    panmanUtils::Tree* tree = &treeGroup->trees[0];
    
    msg("Using tree with {} nodes and {} blocks for testing", 
        tree->allNodes.size(), tree->blocks.size());
    
    // Create test FASTQ file
    if (!createTestFastq(test_fastq, tree)) {
        err("Failed to create test FASTQ file");
        return 1;
    }
    
    // Run gap range validation test
    if (!testGapRangeValidation(tree)) {
        err("Gap range validation test failed");
        return 1;
    }
    
    // All tests passed
    msg("All integration tests passed!");
    
    // Clean up
    delete treeGroup;
    
    return 0;
} 