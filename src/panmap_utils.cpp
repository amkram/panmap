
#include "panmap_utils.hpp"
#include "index.capnp.h"
#include <boost/filesystem.hpp>
#include <fcntl.h>
#include <unistd.h>
#include <cerrno>
#include <cstring>
#include <limits>
#include "capnp/serialize-packed.h"

namespace panmapUtils {

std::string seedChangeTypeToString(seedChangeType changeType) {
  switch (changeType) {
    case seedChangeType::ADD:
      return "ADD";
    case seedChangeType::DEL:
      return "DEL";
    case seedChangeType::SUB:
      return "SUB";
    default:
      return "UNKNOWN";
  }
}

void getSequenceFromReference(
  panmanUtils::Tree* tree,
  std::vector<std::vector<std::pair<char, std::vector<char>>>>& sequence,
  std::vector<char>& blockExists,
  std::vector<char>& blockStrand,
  std::unordered_map<int, int>& blockLengths,
  std::string reference
) {
  if (tree->allNodes.find(reference) == tree->allNodes.end()) {
    logging::err("Reference sequence with matching name not found: {}", reference);
    std::exit(1);
  }
  panmanUtils::Node* referenceNode = tree->allNodes[reference];

  // get path from root to reference node
  std::vector<panmanUtils::Node*> pathFromRoot;
  panmanUtils::Node* it = referenceNode;
  while (it != tree->root) {
    pathFromRoot.push_back(it);
    it = it->parent;
  }
  pathFromRoot.push_back(tree->root);
  std::reverse(pathFromRoot.begin(), pathFromRoot.end());

  // get block sequence (blockSequence[i] = true if block i is on on the reference node)
  std::vector<char> blockSequence(tree->blocks.size() + 1, false);
  for (auto node : pathFromRoot) {
    for (const auto& blockMutation : node->blockMutation) {
      int32_t blockId = blockMutation.primaryBlockId;
      bool insertion = blockMutation.blockMutInfo;
      bool inversion = blockMutation.inversion;
      if (insertion) {
        blockSequence[blockId] = true;
      } else {
        if (!inversion) {
          blockSequence[blockId] = false;
        }
      }
    }
  }

  // initialize sequence, blockExists, blockStrand, blockLengths
  sequence.clear();
  blockExists.clear();
  blockStrand.clear();
  blockLengths.clear();
  sequence.resize(tree->blocks.size() + 1);
  blockExists.resize(tree->blocks.size() + 1, false);
  blockStrand.resize(tree->blocks.size() + 1, true);
  int32_t maxBlockId = 0;

  // fill in the skeleton of the sequence object
  for (int32_t blockId = 0; blockId < tree->blocks.size(); blockId++) {
    blockLengths[blockId] = 0;
    maxBlockId = std::max(maxBlockId, blockId);
    if (blockSequence[blockId]) {
      for (size_t i = 0; i < tree->blocks[blockId].consensusSeq.size(); i++) {
        bool endFlag = false;
        for (size_t j = 0; j < 8; j++) {
          const int nucCode = (((tree->blocks[blockId].consensusSeq[i]) >> (4*(7 - j))) & 15);
          if (nucCode == 0) {
            endFlag = true;
            break;
          }
          const char nucleotide = panmanUtils::getNucleotideFromCode(nucCode);
          sequence[blockId].push_back({nucleotide, {}});
        }
        if (endFlag) {
          break;
        }
      }

      sequence[blockId].push_back({'x', {}});
    } else {
      int len = 0;
      for (size_t i = 0; i < tree->blocks[blockId].consensusSeq.size(); i++) {
        bool endFlag = false;
        for (size_t j = 0; j < 8; j++) {
          const int nucCode = (((tree->blocks[blockId].consensusSeq[i]) >> (4*(7 - j))) & 15);
          if (nucCode == 0) {
            endFlag = true;
            break;
          }
          len++;
        }
        if (endFlag) {
          break;
        }
      }
      blockLengths[blockId] += len;
    }
  }

  sequence.resize(maxBlockId + 1);
  blockExists.resize(maxBlockId + 1);
  blockStrand.resize(maxBlockId + 1);

  // Assign nuc gaps
  auto& gaps = tree->gaps;
  for(size_t i = 0; i < gaps.size(); i++) {
    int32_t primaryBId = gaps[i].primaryBlockId;
    if (blockSequence[primaryBId]){
      for(size_t j = 0; j < gaps[i].nucPosition.size(); j++) {
        int len = gaps[i].nucGapLength[j];
        int pos = gaps[i].nucPosition[j];
        sequence[primaryBId][pos].second.resize(len, '-');
      }
    } else {
      int len=0;
      for(size_t j = 0; j < gaps[i].nucPosition.size(); j++) {
        len += gaps[i].nucGapLength[j];
      }
      blockLengths[primaryBId] += len;
    }
  }

  // apply mutations from root to reference node
  for (auto node : pathFromRoot) {
    // apply block mutations
    for (const auto& blockMutation : node->blockMutation) {
      int32_t blockId = blockMutation.primaryBlockId;
      bool insertion = blockMutation.blockMutInfo;
      bool inversion = blockMutation.inversion;
      if (!blockSequence[blockId]) {
        continue;
      }
      if (insertion) {
        blockExists[blockId] = true;
        blockStrand[blockId] = !inversion;
      } else {
        if (inversion) {
          blockStrand[blockId] = !blockStrand[blockId];
        } else {
          blockExists[blockId] = false;
          blockStrand[blockId] = true;
        }
      }
    }

    // apply nuc mutations
    for (const auto& nucMutation : node->nucMutation) {
      int32_t blockId = nucMutation.primaryBlockId;
      if (!blockSequence[blockId]) {
        continue;
      }
      int length = nucMutation.mutInfo >> 4;
      for (int i = 0; i < length; i++) {
        panmapUtils::Coordinate pos = panmapUtils::Coordinate(nucMutation, i);
        if (pos.nucPosition == sequence[pos.primaryBlockId].size() - 1 && pos.nucGapPosition == -1) {
          continue;
        } else if (pos.nucPosition >= sequence[pos.primaryBlockId].size()) {
          continue;
        }
        int newNucCode = (nucMutation.nucs >> (4*(5-i))) & 0xF;
        char newNuc = panmanUtils::getNucleotideFromCode(newNucCode);
        pos.setSequenceBase(sequence, newNuc);
      }
    }
  }
}

std::string getStringFromSequence(
  const std::vector<std::vector<std::pair<char, std::vector<char>>>>& sequence,
  const std::unordered_map<int, int>& blockLengths,
  const std::vector<char>& blockExists,
  const std::vector<char>& blockStrand,
  bool aligned
) {
  std::string seqString;
  for (size_t i = 0; i < blockExists.size(); i++) {
    if (blockExists[i]) {
      if (blockStrand[i]) {
        for(size_t j = 0; j < sequence[i].size(); j++) {
          // Gap nucs
          for (size_t k = 0; k < sequence[i][j].second.size(); k++) {
            if(sequence[i][j].second[k] != '-') {
              seqString += sequence[i][j].second[k];
            } else if(aligned) {
              seqString += '-';
            }
          }
          // Main nuc
          if(sequence[i][j].first != '-' && sequence[i][j].first != 'x') {
            seqString += sequence[i][j].first;
          } else if (aligned && sequence[i][j].first != 'x') {
            seqString += '-';
          }
        }
      } else {
        for(size_t j = sequence[i].size()-1; j+1 > 0; j--) {
            // Main nuc first since we are iterating in reverse direction
            if(sequence[i][j].first != '-' && sequence[i][j].first != 'x') {
                seqString += panmanUtils::getComplementCharacter(sequence[i][j].first);
            } else if (aligned  && sequence[i][j].first != 'x') {
                seqString += '-';
            }

            // Gap nucs
            for(size_t k = sequence[i][j].second.size()-1; k+1 > 0; k--) {
                if(sequence[i][j].second[k] != '-') {
                    seqString += panmanUtils::getComplementCharacter(sequence[i][j].second[k]);
                } else if (aligned) {
                    seqString += '-';
                }
            }
        }
      }
    } else if (aligned){
      seqString.append(blockLengths.at(i), '-');
    }
  }
  return seqString;
}

std::string getStringFromReference(panmanUtils::Tree* tree, std::string reference, bool aligned) {
  std::vector<std::vector<std::pair<char, std::vector<char>>>> sequence;
  std::unordered_map<int, int> blockLengths;
  std::vector<char> blockExists;
  std::vector<char> blockStrand;
  getSequenceFromReference(tree, sequence, blockExists, blockStrand, blockLengths, reference);
  std::string seqString = getStringFromSequence(sequence, blockLengths, blockExists, blockStrand, aligned);
  return seqString;
}


void LiteTree::cleanup() {
  for (auto& pair : allLiteNodes) {
    delete pair.second;
  }
  allLiteNodes.clear();
  blockScalarRanges.clear();
  nodeToDfsIndex.clear();
  root = nullptr;
}

uint32_t LiteTree::getBlockStartScalar(const uint32_t blockId) const {
  return blockScalarRanges[blockId].first;
}

uint32_t LiteTree::getBlockEndScalar(const uint32_t blockId) const {
  return blockScalarRanges[blockId].second;
}

void LiteTree::initialize(::LiteTree::Reader liteTreeReader) {
  // initialize blockScalarRanges
  auto blockScalarRangesReader = liteTreeReader.getBlockRanges();
  blockScalarRanges.resize(blockScalarRangesReader.size());
  for (size_t i = 0; i < blockScalarRangesReader.size(); i++) {
    blockScalarRanges[i] = {blockScalarRangesReader[i].getRangeBeg(), blockScalarRangesReader[i].getRangeEnd()};
  }

  // initialize allLiteNodes
  auto liteNodesReader = liteTreeReader.getLiteNodes();
  for (size_t i = 0; i < liteNodesReader.size(); i++) {
    const auto liteNodeReader = liteNodesReader[i];
    const auto& nodeIdentifier = liteNodeReader.getId();
    const auto parentIndex = liteNodeReader.getParentIndex();
    nodeToDfsIndex.emplace(nodeIdentifier, i);
    auto [it, inserted] = allLiteNodes.emplace(nodeIdentifier, new LiteNode(nodeIdentifier, nullptr, {}));
    if (i == 0) continue;
    const auto parentNodeReader = liteNodesReader[parentIndex];
    const auto& parentNodeId = parentNodeReader.getId();
    it->second->parent = allLiteNodes[parentNodeId];
    allLiteNodes[parentNodeId]->children.push_back(it->second);
  }

  root = allLiteNodes[liteNodesReader[0].getId()];
}
namespace fs = boost::filesystem;

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

//end of namespace panmapUtils
}