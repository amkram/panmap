#ifndef VISUALIZATION_HPP
#define VISUALIZATION_HPP

#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <map>
#include <set>

// Forward declarations for coordinates types to avoid circular dependencies
namespace coordinates {
    typedef std::vector<std::pair<bool, std::vector<bool>>> blockExists_t;
    typedef std::vector<std::pair<bool, std::vector<bool>>> blockStrand_t;
    typedef std::map<int64_t, int64_t> GapMap;
}

namespace visualization {

/**
 * @brief Create a symbolic representation of block structure
 * @param blockExists Whether blocks exist (on/off)
 * @param blockStrand Strand orientation (forward/inverted)
 * @return ASCII representation of blocks
 */
inline std::string visualizeBlocks(const coordinates::blockExists_t &blockExists,
                           const coordinates::blockStrand_t &blockStrand) {
    std::stringstream ss;
    
    for (size_t i = 0; i < blockExists.size(); ++i) {
        if (!blockExists[i].first) {
            ss << "_"; // Block is off
        } else {
            // Block is on - check strand
            if (blockStrand[i].first) {
                ss << ">"; // Forward strand
            } else {
                ss << "<"; // Inverted strand
            }
        }
    }
    
    return ss.str();
}

/**
 * @brief Detect overlapping gap runs in a gap map
 * @param gapMap The gap map to check
 * @return Vector of overlapping gap pairs
 */
inline std::vector<std::pair<std::pair<int64_t, int64_t>, std::pair<int64_t, int64_t>>> 
detectOverlappingGaps(const coordinates::GapMap &gapMap) {
    std::vector<std::pair<std::pair<int64_t, int64_t>, std::pair<int64_t, int64_t>>> overlaps;
    
    // Convert gaps to (start, end) format for easier overlap checking
    std::vector<std::pair<int64_t, int64_t>> gapRanges;
    for (const auto &[start, length] : gapMap) {
        gapRanges.emplace_back(start, start + length - 1);
    }
    
    // Sort by start position
    std::sort(gapRanges.begin(), gapRanges.end());
    
    // Check for overlaps
    for (size_t i = 0; i < gapRanges.size(); i++) {
        for (size_t j = i + 1; j < gapRanges.size(); j++) {
            // If gap i ends after gap j starts, they overlap
            if (gapRanges[i].second >= gapRanges[j].first) {
                overlaps.push_back({gapRanges[i], gapRanges[j]});
            }
            // Once we find a gap that starts after i ends, we can stop checking
            if (gapRanges[j].first > gapRanges[i].second) {
                break;
            }
        }
    }
    
    return overlaps;
}

/**
 * @brief Calculate the true coverage of gaps accounting for overlaps
 * @param gapMap The gap map
 * @param totalLength The total sequence length
 * @return Pair of (coverage length, percentage)
 */
inline std::pair<int64_t, double> calculateTrueGapCoverage(
    const coordinates::GapMap &gapMap, int64_t totalLength) {
    
    // Create a set of all positions covered by gaps
    std::set<int64_t> coveredPositions;
    
    for (const auto &[start, length] : gapMap) {
        for (int64_t pos = start; pos < start + length; pos++) {
            coveredPositions.insert(pos);
        }
    }
    
    int64_t trueCoverage = coveredPositions.size();
    double percentage = (totalLength > 0) ? 
        (100.0 * trueCoverage / totalLength) : 0.0;
    
    return {trueCoverage, percentage};
}

/**
 * @brief Calculate the true coverage of gaps accounting for overlaps, 
 * with percentage relative to active blocks only
 * @param gapMap The gap map
 * @param totalLength The total sequence length
 * @param blockExists Block existence info
 * @param blockStrand Block strand info
 * @return Pair of (coverage length, percentage relative to active blocks)
 * 
 * Note: This function estimates the active block length based on the proportion
 * of active blocks to total blocks, assuming blocks are roughly equal in size.
 * If blocks have significantly different sizes, this estimate may not be precise.
 * 
 * The function only counts gaps that fall within active blocks for the percentage
 * calculation, providing a more accurate representation of gap density in active regions.
 */
inline std::pair<int64_t, double> calculateTrueGapCoverageActiveOnly(
    const coordinates::GapMap &gapMap, 
    int64_t totalLength,
    const coordinates::blockExists_t &blockExists,
    const coordinates::blockStrand_t &blockStrand) {
    
    // Count active blocks (on blocks)
    int activeBlockCount = 0;
    int totalBlockCount = blockExists.size();
    
    for (size_t i = 0; i < blockExists.size(); ++i) {
        if (blockExists[i].first) {
            // Block exists (is on)
            activeBlockCount++;
        }
    }
    
    // Calculate active block percentage of the total sequence
    // We assume blocks are roughly equal in size for this calculation
    double activeBlockRatio = (totalBlockCount > 0) ? 
        static_cast<double>(activeBlockCount) / totalBlockCount : 1.0;
    
    // Estimate active sequence length based on the proportion of active blocks
    int64_t activeLength = static_cast<int64_t>(totalLength * activeBlockRatio);
    
    // Ensure we don't divide by zero
    if (activeLength <= 0) {
        activeLength = 1; // Avoid division by zero
    }
    
    // Calculate block boundaries (approximate)
    // This assumes blocks are evenly distributed across the sequence
    std::vector<std::pair<int64_t, int64_t>> activeBlockRanges;
    if (totalBlockCount > 0) {
        // Handle the case where there are no active blocks
        if (activeBlockCount == 0) {
            // No active blocks, so no gaps in active blocks
            return {0, 0.0};
        }
        
        // Calculate approximate block size
        int64_t blockSize = totalLength / totalBlockCount;
        if (blockSize == 0) blockSize = 1; // Ensure minimum block size
        
        // Calculate boundaries for each active block
        for (size_t i = 0; i < blockExists.size(); ++i) {
            if (blockExists[i].first) {
                // This is an active block
                int64_t blockStart = i * blockSize;
                int64_t blockEnd = (i + 1) * blockSize - 1;
                
                // Ensure the last block extends to the end of the sequence
                if (i == blockExists.size() - 1) {
                    blockEnd = totalLength - 1;
                }
                
                // Ensure we don't exceed sequence bounds
                blockEnd = std::min(blockEnd, totalLength - 1);
                
                // Only add valid ranges
                if (blockStart <= blockEnd && blockStart < totalLength) {
                    activeBlockRanges.emplace_back(blockStart, blockEnd);
                }
            }
        }
    }
    
    // If we couldn't determine any active block ranges, return zero
    if (activeBlockRanges.empty()) {
        return {0, 0.0};
    }
    
    // Create a set of all positions covered by gaps WITHIN ACTIVE BLOCKS
    std::set<int64_t> coveredPositions;
    
    for (const auto &[start, length] : gapMap) {
        for (int64_t pos = start; pos < start + length; pos++) {
            // Check if this position falls within any active block
            bool inActiveBlock = false;
            for (const auto &[blockStart, blockEnd] : activeBlockRanges) {
                if (pos >= blockStart && pos <= blockEnd) {
                    inActiveBlock = true;
                    break;
                }
            }
            
            if (inActiveBlock) {
                coveredPositions.insert(pos);
            }
        }
    }
    
    int64_t trueCoverage = coveredPositions.size();
    
    // Calculate percentage relative to active blocks only
    double percentage = 100.0 * trueCoverage / activeLength;
    
    return {trueCoverage, percentage};
}

/**
 * @brief Create a symbolic representation of a gap map
 * @param gapMap The gap map to visualize
 * @param totalLength Total sequence length for scaling
 * @param width Desired width of the visualization
 * @return ASCII representation of the gap map
 */
inline std::string visualizeGapMap(const coordinates::GapMap &gapMap, 
                          int64_t totalLength,
                          int width = 80) {
    // Prepare visualization components
    std::stringstream ss;
    std::string viz(width, '-'); // Base visualization with '-' characters
    std::vector<std::pair<int, std::string>> annotations; // To store positions and labels
    
    // Check for overlapping gaps
    auto overlaps = detectOverlappingGaps(gapMap);
    bool hasOverlaps = !overlaps.empty();
    
    // Calculate the true coverage accounting for overlaps
    auto [trueCoverage, truePercentage] = calculateTrueGapCoverage(gapMap, totalLength);
    
    // Calculate naïve total (this may exceed totalLength if gaps overlap)
    int64_t rawTotalGapLength = 0;
    for (const auto &[_, length] : gapMap) {
        rawTotalGapLength += length;
    }
    double rawPercentage = (totalLength > 0) ? 
        (100.0 * rawTotalGapLength / totalLength) : 0.0;
    
    // Scale factor to map coordinates to visualization width
    double scale = static_cast<double>(width) / totalLength;
    
    // Maintain a count of overlaps at each position for visualization
    std::vector<int> overlapCounts(width, 0);
    
    // Process each gap run
    for (const auto &[start, length] : gapMap) {
        int64_t end = start + length - 1;
        int vizStart = static_cast<int>(start * scale);
        int vizEnd = static_cast<int>(end * scale);
        
        // Ensure at least one character for small gaps
        vizEnd = std::max(vizStart, vizEnd);
        
        // Mark the gap region
        for (int i = vizStart; i <= vizEnd && i < width; ++i) {
            viz[i] = 'G';
            overlapCounts[i]++;
        }
        
        // Add annotation for significant gap runs (longer than a threshold)
        if (length > 10) {
            std::string label = std::to_string(length);
            int labelPos = (vizStart + vizEnd) / 2; // Center the label
            annotations.push_back({labelPos, label});
        }
    }
    
    // Mark positions with overlapping gaps using a different character
    if (hasOverlaps) {
        for (int i = 0; i < width; ++i) {
            if (overlapCounts[i] > 1) {
                viz[i] = 'O'; // Use 'O' to indicate overlapping gaps
            }
        }
    }
    
    // Render the base visualization
    ss << viz << std::endl;
    
    // Add scale markers
    ss << "0" << std::string(width-10, ' ') << totalLength << std::endl;
    
    // Add annotations for gap lengths (avoiding overlaps)
    if (!annotations.empty()) {
        // Sort annotations by position
        std::sort(annotations.begin(), annotations.end());
        
        // Create annotation line
        std::string annLine(width, ' ');
        for (const auto& [pos, label] : annotations) {
            // Ensure we don't write outside bounds
            if (pos >= 0 && pos < width) {
                // Try to center the label, but ensure we stay within bounds
                int startPos = std::max(0, pos - (int)label.length()/2);
                int endPos = std::min(width-1, startPos + (int)label.length()-1);
                
                // Check if there's already text here (from another annotation)
                bool hasOverlap = false;
                for (int i = startPos; i <= endPos; i++) {
                    if (annLine[i] != ' ') {
                        hasOverlap = true;
                        break;
                    }
                }
                
                // If no overlap, add the label
                if (!hasOverlap) {
                    for (size_t i = 0; i < label.length() && (startPos + i) < width; i++) {
                        annLine[startPos + i] = label[i];
                    }
                }
            }
        }
        ss << annLine << std::endl;
    }
    
    // Add summary line
    int64_t nucleotideCount = totalLength - trueCoverage;
    
    // Forced to be non-negative
    if (nucleotideCount < 0) {
        ss << "WARNING: INVALID NUCLEOTIDE COUNT - Gap coverage exceeds sequence length!" << std::endl;
        nucleotideCount = 0;
    }
    
    // Show detailed gap information
    if (hasOverlaps) {
        ss << gapMap.size() << " gaps with " << overlaps.size() << " overlaps!" << std::endl;
        ss << "Raw gap length: " << rawTotalGapLength << " (" << std::fixed << std::setprecision(2) 
           << rawPercentage << "% of sequence)" << std::endl;
        ss << "True gap coverage: " << trueCoverage << " (" << std::fixed << std::setprecision(2) 
           << truePercentage << "% of sequence)" << std::endl;
        
        // Show first few overlaps for debugging
        int overlapsToPrint = std::min(5, (int)overlaps.size());
        ss << "First " << overlapsToPrint << " overlaps:" << std::endl;
        for (int i = 0; i < overlapsToPrint; i++) {
            auto &overlap = overlaps[i];
            ss << "  [" << overlap.first.first << "," << overlap.first.second 
               << "] and [" << overlap.second.first << "," << overlap.second.second << "]" << std::endl;
        }
    } else {
        ss << gapMap.size() << " gaps covering " 
           << std::fixed << std::setprecision(2) << truePercentage << "% of sequence";
    }
    
    return ss.str();
}

/**
 * @brief Create a symbolic representation of a gap map with active block percentage
 * @param gapMap The gap map to visualize
 * @param totalLength Total sequence length for scaling
 * @param blockExists Block existence info
 * @param blockStrand Block strand info
 * @param width Desired width of the visualization
 * @return ASCII representation of the gap map
 */
inline std::string visualizeGapMap(const coordinates::GapMap &gapMap, 
                          int64_t totalLength,
                          const coordinates::blockExists_t &blockExists,
                          const coordinates::blockStrand_t &blockStrand,
                          int width = 80) {
    // Prepare visualization components
    std::stringstream ss;
    std::string viz(width, '-'); // Base visualization with '-' characters
    std::vector<std::pair<int, std::string>> annotations; // To store positions and labels
    
    // Check for overlapping gaps
    auto overlaps = detectOverlappingGaps(gapMap);
    bool hasOverlaps = !overlaps.empty();
    
    // Calculate the true coverage accounting for overlaps
    auto [trueCoverage, truePercentage] = calculateTrueGapCoverage(gapMap, totalLength);
    
    // Calculate the active-blocks-only percentage
    auto [_, activeBlockPercentage] = calculateTrueGapCoverageActiveOnly(
        gapMap, totalLength, blockExists, blockStrand);
    
    // Note: The active block percentage is an estimate based on the assumption
    // that blocks are roughly equal in size. The actual percentage may vary
    // if blocks have significantly different sizes.
    
    // Calculate naïve total (this may exceed totalLength if gaps overlap)
    int64_t rawTotalGapLength = 0;
    for (const auto &[_, length] : gapMap) {
        rawTotalGapLength += length;
    }
    double rawPercentage = (totalLength > 0) ? 
        (100.0 * rawTotalGapLength / totalLength) : 0.0;
    
    // Scale factor to map coordinates to visualization width
    double scale = static_cast<double>(width) / totalLength;
    
    // Maintain a count of overlaps at each position for visualization
    std::vector<int> overlapCounts(width, 0);
    
    // Process each gap run
    for (const auto &[start, length] : gapMap) {
        int64_t end = start + length - 1;
        int vizStart = static_cast<int>(start * scale);
        int vizEnd = static_cast<int>(end * scale);
        
        // Ensure at least one character for small gaps
        vizEnd = std::max(vizStart, vizEnd);
        
        // Mark the gap region
        for (int i = vizStart; i <= vizEnd && i < width; ++i) {
            viz[i] = 'G';
            overlapCounts[i]++;
        }
        
        // Add annotation for significant gap runs (longer than a threshold)
        if (length > 10) {
            std::string label = std::to_string(length);
            int labelPos = (vizStart + vizEnd) / 2; // Center the label
            annotations.push_back({labelPos, label});
        }
    }
    
    // Mark positions with overlapping gaps using a different character
    if (hasOverlaps) {
        for (int i = 0; i < width; ++i) {
            if (overlapCounts[i] > 1) {
                viz[i] = 'O'; // Use 'O' to indicate overlapping gaps
            }
        }
    }
    
    // Render the base visualization
    ss << viz << std::endl;
    
    // Add scale markers
    ss << "0" << std::string(width-10, ' ') << totalLength << std::endl;
    
    // Add annotations for gap lengths (avoiding overlaps)
    if (!annotations.empty()) {
        // Sort annotations by position
        std::sort(annotations.begin(), annotations.end());
        
        // Create annotation line
        std::string annLine(width, ' ');
        for (const auto& [pos, label] : annotations) {
            // Ensure we don't write outside bounds
            if (pos >= 0 && pos < width) {
                // Try to center the label, but ensure we stay within bounds
                int startPos = std::max(0, pos - (int)label.length()/2);
                int endPos = std::min(width-1, startPos + (int)label.length()-1);
                
                // Check if there's already text here (from another annotation)
                bool hasOverlap = false;
                for (int i = startPos; i <= endPos; i++) {
                    if (annLine[i] != ' ') {
                        hasOverlap = true;
                        break;
                    }
                }
                
                // If no overlap, add the label
                if (!hasOverlap) {
                    for (size_t i = 0; i < label.length() && (startPos + i) < width; i++) {
                        annLine[startPos + i] = label[i];
                    }
                }
            }
        }
        ss << annLine << std::endl;
    }
    
    // Add summary line
    int64_t nucleotideCount = totalLength - trueCoverage;
    
    // Forced to be non-negative
    if (nucleotideCount < 0) {
        ss << "WARNING: INVALID NUCLEOTIDE COUNT - Gap coverage exceeds sequence length!" << std::endl;
        nucleotideCount = 0;
    }
    
    // Count active blocks
    int onBlocks = 0, invertedBlocks = 0, offBlocks = 0;
    for (size_t i = 0; i < blockExists.size(); ++i) {
        if (!blockExists[i].first) {
            offBlocks++;
        } else if (!blockStrand[i].first) {
            invertedBlocks++;
        } else {
            onBlocks++;
        }
    }
    
    // Show detailed gap information
    if (hasOverlaps) {
        ss << gapMap.size() << " gaps with " << overlaps.size() << " overlaps!" << std::endl;
        ss << "Raw gap length: " << rawTotalGapLength << " (" << std::fixed << std::setprecision(2) 
           << rawPercentage << "% of sequence)" << std::endl;
        ss << "True gap coverage: " << trueCoverage << " (" << std::fixed << std::setprecision(2) 
           << truePercentage << "% of sequence, " << activeBlockPercentage << "% of active blocks - gaps in active blocks only)" << std::endl;
        ss << "Block structure: " << onBlocks << " on, " << invertedBlocks << " inverted, " 
           << offBlocks << " off (total: " << blockExists.size() << ")" << std::endl;
        
        // Show first few overlaps for debugging
        int overlapsToPrint = std::min(5, (int)overlaps.size());
        ss << "First " << overlapsToPrint << " overlaps:" << std::endl;
        for (int i = 0; i < overlapsToPrint; i++) {
            auto &overlap = overlaps[i];
            ss << "  [" << overlap.first.first << "," << overlap.first.second 
               << "] and [" << overlap.second.first << "," << overlap.second.second << "]" << std::endl;
        }
    } else {
        ss << gapMap.size() << " gaps covering " 
           << std::fixed << std::setprecision(2) << truePercentage << "% of sequence, "
           << activeBlockPercentage << "% of active blocks (gaps in active blocks only)";
    }
    
    return ss.str();
}

/**
 * @brief Print node visualization to stdout
 * @param nodeName Name of the node
 * @param blockExists Block existence info
 * @param blockStrand Block strand info
 * @param gapMap Gap map for the node
 * @param totalLength Total sequence length
 */
inline void printNodeVisualization(const std::string &nodeName,
                          const coordinates::blockExists_t &blockExists,
                          const coordinates::blockStrand_t &blockStrand,
                          const coordinates::GapMap &gapMap,
                          int64_t totalLength) {
    std::cout << "Node: " << nodeName << std::endl;
    
    // Count block types
    int onBlocks = 0, invertedBlocks = 0, offBlocks = 0;
    for (size_t i = 0; i < blockExists.size(); ++i) {
        if (!blockExists[i].first) {
            offBlocks++;
        } else if (!blockStrand[i].first) {
            invertedBlocks++;
        } else {
            onBlocks++;
        }
    }
    
    // Calculate true nucleotide count accounting for overlapping gaps
    // Use the active-blocks-only calculation for percentage
    auto [trueCoverage, activeBlockPercentage] = calculateTrueGapCoverageActiveOnly(
        gapMap, totalLength, blockExists, blockStrand);
    
    // Calculate the total number of non-gap nucleotides
    // This is the total length minus the number of positions covered by gaps
    auto [totalGapCoverage, _] = calculateTrueGapCoverage(gapMap, totalLength);
    int64_t nucleotideCount = totalLength - totalGapCoverage;
    
    // Also calculate the regular percentage for reference
    double regularPercentage = (totalLength > 0) ? 
        (100.0 * trueCoverage / totalLength) : 0.0;
    
    // Calculate raw gap length (could exceed totalLength if gaps overlap)
    int64_t rawGapLength = 0;
    for (const auto &[_, length] : gapMap) {
        rawGapLength += length;
    }
    
    // Detect gap overlaps
    auto overlaps = detectOverlappingGaps(gapMap);
    bool hasOverlaps = !overlaps.empty();
    
    // Print block structure
    std::cout << "Blocks: " << visualizeBlocks(blockExists, blockStrand) << std::endl;
    std::cout << "       " << onBlocks << " on, " 
              << invertedBlocks << " inverted, " 
              << offBlocks << " off (total: " << blockExists.size() << ")" << std::endl;
    
    // Print nucleotide stats with validation
    if (nucleotideCount < 0) {
        std::cout << "Nucs:  ⚠️ INVALID - GAP OVERLAPS DETECTED! ⚠️" << std::endl;
        std::cout << "       Raw data shows " << nucleotideCount << " nucleotides, " 
                  << rawGapLength << " gap length (exceeds total: " << totalLength << ")" << std::endl;
        std::cout << "       Corrected: 0 nucleotides, " << trueCoverage << " unique gap positions" << std::endl;
    } else {
        if (hasOverlaps) {
            std::cout << "Nucs:  " << nucleotideCount << " nucleotides, " 
                      << trueCoverage << " gap positions (total: " << totalLength << ")" << std::endl;
            std::cout << "       Gaps: " << std::fixed << std::setprecision(2) 
                      << regularPercentage << "% of total sequence, " 
                      << activeBlockPercentage << "% of active blocks (gaps in active blocks only)" << std::endl;
            std::cout << "       ⚠️ WARNING: " << overlaps.size() << " overlapping gap runs detected! ⚠️" << std::endl;
            std::cout << "       Raw gap length sum is " << rawGapLength << " which " 
                      << (rawGapLength > totalLength ? "exceeds" : "equals") 
                      << " total length" << std::endl;
        } else {
            std::cout << "Nucs:  " << nucleotideCount << " nucleotides, " 
                      << trueCoverage << " gaps (total length: " << totalLength << ")" << std::endl;
            std::cout << "       Gaps: " << std::fixed << std::setprecision(2) 
                      << regularPercentage << "% of total sequence, " 
                      << activeBlockPercentage << "% of active blocks (gaps in active blocks only)" << std::endl;
        }
    }
    
    // Print gap map
    std::cout << "Gaps:  " << visualizeGapMap(gapMap, totalLength, blockExists, blockStrand) << std::endl;
    std::cout << std::endl;
}

/**
 * @brief Check if a gap map is valid (no overlaps, gaps within range)
 * @param gapMap The gap map to validate
 * @param totalLength The total sequence length
 * @return True if valid, false if issues found
 */
inline bool validateGapMap(const coordinates::GapMap &gapMap, int64_t totalLength) {
    bool valid = true;
    std::stringstream issues;
    
    // Check for negative or zero lengths
    for (const auto &[start, length] : gapMap) {
        if (length <= 0) {
            issues << "Gap at position " << start << " has invalid length " << length << std::endl;
            valid = false;
        }
    }
    
    // Check for gaps starting outside sequence range
    for (const auto &[start, length] : gapMap) {
        if (start < 0) {
            issues << "Gap starts at negative position " << start << std::endl;
            valid = false;
        }
        if (start >= totalLength) {
            issues << "Gap starts at " << start << " which is beyond sequence length " 
                   << totalLength << std::endl;
            valid = false;
        }
    }
    
    // Check for gaps extending beyond sequence
    for (const auto &[start, length] : gapMap) {
        if (start + length > totalLength) {
            issues << "Gap at position " << start << " with length " << length 
                   << " extends beyond sequence length " << totalLength << std::endl;
            valid = false;
        }
    }
    
    // Check for overlapping gaps
    auto overlaps = detectOverlappingGaps(gapMap);
    if (!overlaps.empty()) {
        issues << "Found " << overlaps.size() << " overlapping gap runs:" << std::endl;
        // Show first few overlaps
        int overlapsToPrint = std::min(5, (int)overlaps.size());
        for (int i = 0; i < overlapsToPrint; i++) {
            auto &overlap = overlaps[i];
            issues << "  [" << overlap.first.first << "," << overlap.first.second 
                   << "] and [" << overlap.second.first << "," << overlap.second.second << "]" << std::endl;
        }
        valid = false;
    }
    
    // Print issues if validation failed
    if (!valid) {
        std::cerr << "Gap map validation failed:" << std::endl;
        std::cerr << issues.str();
    }
    
    return valid;
}

} // namespace visualization

#endif // VISUALIZATION_HPP 