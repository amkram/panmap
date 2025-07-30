#include "alignment.hpp"
#include "capnp/list.h"
#include "conversion.hpp"
#include "htslib/sam.h"
#include "index.capnp.h"
#include "panman.hpp"
#include "seeding.hpp"
#include <algorithm>
#include <bits/getopt_core.h>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <tuple>
#include <unistd.h>
#include <unordered_map>
#include <utility>
#include <vector>

extern "C" {
#include "../src/3rdparty/bwa/bwa.h"
}

using namespace panmanUtils;

int extractPosition(char *line) {
  int fieldCount = 0;
  char *ptr = line;
  char *fieldStart = ptr;

  while (*ptr != '\0') {
    if (*ptr == '\t') {
      fieldCount++;
      if (fieldCount == 4) {
        // We've found the end of the 4th field
        *ptr = '\0'; // Temporarily terminate the 4th field string

        // Convert the 4th field to an integer
        int position = atoi(fieldStart);

        *ptr = '\t'; // Restore the original character
        return position;
      }
      // Move to the start of the next field
      fieldStart = ptr + 1;
    }
    ptr++;
  }

  // Check if the line ends exactly after the 4th field
  if (fieldCount == 3) {
    // Line ends after the 4th field
    int position = atoi(fieldStart);
    return position;
  }

  // If we reach here, the line doesn't have at least 4 fields
  throw std::runtime_error("Line does not have at least 4 fields.");
}

void prepareAndRunBwa(const std::vector<std::string> &idx_args,
                      const std::vector<std::string> &aln_args1,
                      const std::vector<std::string> &aln_args2,
                      const std::vector<std::string> &samaln_args,
                      const std::string &reads1Path,
                      const std::string &reads2Path,
                      std::vector<std::pair<int, char *>> &samAlignmentPairs,
                      std::vector<std::string> &samHeaders) {
  std::vector<std::vector<char>> idx_argv_buf;
  std::vector<char *> idx_argv;

  std::cerr << "Constructing idx_argv..." << std::endl;
  for (const auto &arg : idx_args) {
    std::cerr << "Processing argument: " << arg << std::endl;
    idx_argv_buf.push_back(std::vector<char>(arg.begin(), arg.end()));
    idx_argv_buf.back().push_back('\0');            // Ensure null termination
    idx_argv.push_back(idx_argv_buf.back().data()); // Push the C-string pointer
  }
  idx_argv.push_back(nullptr); // Null-terminate the argv array

  std::cerr << "Constructed idx_argv:" << std::endl;
  for (size_t i = 0; i < idx_argv.size(); ++i) {
    if (idx_argv[i] != nullptr)
      std::cerr << "idx_argv[" << i << "]: " << idx_argv[i] << std::endl;
    else
      std::cerr << "idx_argv[" << i << "]: (null)" << std::endl;
  }

  std::vector<std::vector<char>> aln_argv_buf1;
  std::vector<char *> aln_argv1;
  std::cerr << "Constructing aln_argv1..." << std::endl;
  for (const auto &arg : aln_args1) {
    std::cerr << "Processing argument: " << arg << std::endl;
    aln_argv_buf1.push_back(std::vector<char>(arg.begin(), arg.end()));
    aln_argv_buf1.back().push_back('\0');
    aln_argv1.push_back(aln_argv_buf1.back().data());
  }
  aln_argv1.push_back(nullptr);

  std::cerr << "Constructed aln_argv1:" << std::endl;
  for (size_t i = 0; i < aln_argv1.size(); ++i) {
    if (aln_argv1[i] != nullptr)
      std::cerr << "aln_argv1[" << i << "]: " << aln_argv1[i] << std::endl;
    else
      std::cerr << "aln_argv1[" << i << "]: (null)" << std::endl;
  }

  std::vector<std::vector<char>> aln_argv_buf2;
  std::vector<char *> aln_argv2;
  if (aln_args2.size() > 0) {
    std::cerr << "Constructing aln_argv2..." << std::endl;
    for (const auto &arg : aln_args2) {
      std::cerr << "Processing argument: " << arg << std::endl;
      aln_argv_buf2.push_back(std::vector<char>(arg.begin(), arg.end()));
      aln_argv_buf2.back().push_back('\0');
      aln_argv2.push_back(aln_argv_buf2.back().data());
    }
    aln_argv2.push_back(nullptr);

    std::cerr << "Constructed aln_argv2:" << std::endl;
    for (size_t i = 0; i < aln_argv2.size(); ++i) {
      if (aln_argv2[i] != nullptr)
        std::cerr << "aln_argv2[" << i << "]: " << aln_argv2[i] << std::endl;
      else
        std::cerr << "aln_argv2[" << i << "]: (null)" << std::endl;
    }
  }

  std::vector<std::vector<char>> samaln_argv_buf;
  std::vector<char *> samaln_argv;
  std::cerr << "Constructing samaln_argv..." << std::endl;
  for (const auto &arg : samaln_args) {
    std::cerr << "Processing argument: " << arg << std::endl;
    samaln_argv_buf.push_back(std::vector<char>(arg.begin(), arg.end()));
    samaln_argv_buf.back().push_back('\0');
    samaln_argv.push_back(samaln_argv_buf.back().data());
  }
  samaln_argv.push_back(nullptr);

  std::cerr << "Constructed samaln_argv:" << std::endl;
  for (size_t i = 0; i < samaln_argv.size(); ++i) {
    if (samaln_argv[i] != nullptr)
      std::cerr << "samaln_argv[" << i << "]: " << samaln_argv[i] << std::endl;
    else
      std::cerr << "samaln_argv[" << i << "]: (null)" << std::endl;
  }

  std::cerr << "About to call run_bwa with idx_argv..." << std::endl;

  try {
    std::cerr << "Running run_bwa with idx_argv..." << std::endl;
    optind = 1;
    if (run_bwa(idx_argv.size() - 1, idx_argv.data()) != 0) {
      throw std::runtime_error("BWA index failed");
    }
    std::cerr << "Finished run_bwa with idx_argv." << std::endl;

    std::cerr << "Running run_bwa with aln_argv1..." << std::endl;
    optind = 1;
    int stdoutFd = dup(fileno(stdout));
    if (run_bwa(aln_argv1.size() - 1, aln_argv1.data()) != 0) {
      throw std::runtime_error("BWA aln failed");
    }
    fflush(stdout);
    dup2(stdoutFd, fileno(stdout));
    close(stdoutFd);
    std::cerr << "Finished run_bwa with aln_argv1." << std::endl;

    if (aln_args2.size() > 0) {
      std::cerr << "Running run_bwa with aln_argv2..." << std::endl;
      optind = 1;
      int stdoutFd2 = dup(fileno(stdout));
      if (run_bwa(aln_argv2.size() - 1, aln_argv2.data()) != 0) {
        throw std::runtime_error("BWA aln failed");
      }
      fflush(stdout);
      dup2(stdoutFd2, fileno(stdout));
      close(stdoutFd2);
      std::cerr << "Finished run_bwa with aln_argv2." << std::endl;
    }

    std::cerr << "Running run_bwa with samaln_argv..." << std::endl;
    optind = 1;
    FILE *tempFile = std::tmpfile();
    if (!tempFile) {
      std::cerr
          << "Failed to create temporary file for capturing bwa SAM alignments."
          << std::endl;
      return;
    }
    int stdoutFd3 = dup(fileno(stdout));
    if (stdoutFd3 == -1) {
      std::cerr << "Failed to duplicate stdout file descriptor for capturing "
                   "bwa SAM alignments."
                << std::endl;
      return;
    }
    if (dup2(fileno(tempFile), fileno(stdout)) == -1) {
      std::cerr << "Failed to redirect stdout to temporary file for capturing "
                   "bwa SAM alignments."
                << std::endl;
      return;
    }
    if (run_bwa(samaln_argv.size() - 1, samaln_argv.data()) != 0) {
      throw std::runtime_error("BWA samaln failed");
    }
    fflush(stdout);
    dup2(stdoutFd3, fileno(stdout));
    close(stdoutFd3);
    rewind(tempFile);
    char buffer[512];
    size_t i = 0;
    while (fgets(buffer, sizeof(buffer), tempFile)) {
      char *line = new char[strlen(buffer) + 1];
      std::strcpy(line, buffer);
      if (aln_args2.size() > 0) {
        if (i < 4) {
          samHeaders.push_back(line);
        } else {
          size_t len = strlen(line);
          if (len > 0 && line[len - 1] == '\n')
            line[len - 1] = '\0';
          int pos = extractPosition(line);
          samAlignmentPairs.push_back(std::make_pair(pos, line));
        }
      } else {
        if (i < 3) {
          samHeaders.push_back(line);
        } else {
          size_t len = strlen(line);
          if (len > 0 && line[len - 1] == '\n')
            line[len - 1] = '\0';
          int pos = extractPosition(line);
          samAlignmentPairs.push_back(std::make_pair(pos, line));
        }
      }
      ++i;
    }
    fclose(tempFile);
    std::cerr << "Finished run_bwa with samaln_argv." << std::endl;

    // Clean up allocated char arrays
    for (auto &pair : samAlignmentPairs) {
      delete[] pair.second; // Clean up allocated char arrays
    }
  } catch (...) {
    std::cerr << "run_bwa() caused an exception!" << std::endl;
  }
}
void alignment::getAnchors(
    std::vector<std::tuple<int64_t, int32_t, int>> &anchors,
    const std::vector<std::vector<seeding::seed_t>> &readSeeds,
    const std::vector<std::string> &readSequences,
    const std::unordered_map<
        size_t, std::pair<std::vector<uint32_t>, std::vector<uint32_t>>>
        &seedToRefPositions,
    int k) {

  for (size_t i = 0; i < readSequences.size(); ++i) {
    const auto &curReadSeeds = readSeeds[i];

    for (const auto &seed : curReadSeeds) {
      if (seedToRefPositions.find(seed.hash) == seedToRefPositions.end())
        continue;
      auto it = seedToRefPositions.find(seed.hash);
      if (it == seedToRefPositions.end())
        continue;
      const auto &[forwardSeedToRefPositions, reverseSeedToRefPositions] =
          it->second;

      // Add forward anchors
      for (const auto &refPos : forwardSeedToRefPositions) {
        anchors.push_back(std::make_tuple(refPos + k - 1, seed.pos + k - 1, k));
      }

      // Add reverse anchors
      for (const auto &refPos : reverseSeedToRefPositions) {
        anchors.push_back(std::make_tuple(refPos + k - 1, seed.pos + k - 1, k));
      }
    }
  }
}

float alignment::align(
    std::vector<std::optional<seeding::seed_t>> &bestNodeSeeds,
    std::unordered_map<size_t,
                       std::pair<std::vector<uint32_t>, std::vector<uint32_t>>>
        &seedToRefPositions,
    std::string &nodeSequence, panmanUtils::Node *node, panmanUtils::Tree *T,
    int32_t k, int32_t s, int32_t t, bool open, int32_t l,
    ::capnp::List<SeedMutations>::Reader &perNodeSeedMutations_Reader,
    ::capnp::List<GapMutations>::Reader &perNodeGapMutations_Reader,
    const std::string &reads1Path, const std::string &reads2Path,
    std::vector<std::vector<seeding::seed_t>> &readSeeds,
    std::vector<std::string> &readSequences,
    std::vector<std::string> &readQuals, std::vector<std::string> &readNames,
    std::string &samFileName, std::string &bamFileName, bool pairedEndReads,
    std::string &refFileName, std::string aligner) {

  // Get the full sequence from the reference node
  std::string bestMatchSequence =
      T->getStringFromReference(node->identifier, false, true);

  if (refFileName.size() > 0) {

    std::ofstream outFile{refFileName};

    if (outFile.is_open()) {

      // outFile << ">" << node->identifier << "|Panmap.PLACED_SEQUENCE\n";
      outFile << ">ref"
              << "\n";
      outFile << bestMatchSequence << "\n";

      std::cerr << "Wrote reference fasta to " << refFileName << std::endl;
    } else {
      std::cerr << "Error: failed to write to file " << refFileName
                << std::endl;
    }
  }

  double alignmentScore = 0;

  if (aligner == "minimap2") {

    // Align reads to node sequence
    std::vector<char *> samAlignments;
    std::string samHeader;
    createSam(readSeeds, readSequences, readQuals, readNames, nodeSequence,
              seedToRefPositions, samFileName, k, false, pairedEndReads, samAlignments,
              samHeader);

    // Compute alignment metrics
    // metrics = computeAlignmentMetrics(samAlignments, samHeader, nodeSequence,
    // readSequences);

    // Compute final score using default weights
    // metrics.computeScore();
    // alignmentScore = metrics.finalScore;

    // Clean up
    for (char *sam : samAlignments) {
      free(sam);
    }

    if (samFileName.size() == 0) {
      return alignmentScore;
    }

    // Convert to BAM
    sam_hdr_t *header;
    bam1_t **bamRecords;

    createBam(samAlignments, samHeader, bamFileName,

              header, bamRecords);
  } else {
    std::cout << "Aligning with bwa..." << std::endl;

    std::vector<std::string> idx_args = {"bwa", "index", refFileName};
    std::vector<std::string> aln_args1;
    std::vector<std::string> aln_args2;
    std::vector<std::string> samaln_args;
    std::vector<std::pair<int, char *>> samAlignmentPairs;
    std::vector<std::string> samHeaders;
    if (reads2Path.empty()) {
      aln_args1 = {
          "bwa",       "aln",     "-l", "1024", "-n",
          "0.01",      "-o",      "2",  "-f",   reads1Path + ".tmp.sai",
          refFileName, reads1Path};
      samaln_args = {"bwa", "samse", refFileName, reads1Path + ".tmp.sai",
                     reads1Path};
    } else {
      aln_args1 = {
          "bwa",       "aln",     "-l", "1024", "-n",
          "0.01",      "-o",      "2",  "-f",   reads1Path + ".tmp.sai",
          refFileName, reads1Path};
      aln_args2 = {
          "bwa",       "aln",     "-l", "1024", "-n",
          "0.01",      "-o",      "2",  "-f",   reads2Path + ".tmp.sai",
          refFileName, reads2Path};
      samaln_args = {"bwa",
                     "sampe",
                     refFileName,
                     reads1Path + ".tmp.sai",
                     reads2Path + ".tmp.sai",
                     reads1Path,
                     reads2Path};
    }

    prepareAndRunBwa(idx_args, aln_args1, aln_args2, samaln_args, reads1Path,
                     reads2Path, samAlignmentPairs, samHeaders);

    std::sort(
        samAlignmentPairs.begin(), samAlignmentPairs.end(),
        [](const std::pair<int, char *> &a, const std::pair<int, char *> &b) {
          return a.first < b.first;
        });

    std::vector<char *> samAlignments(samAlignmentPairs.size());
    for (size_t i = 0; i < samAlignmentPairs.size(); ++i) {
      samAlignments[i] = samAlignmentPairs[i].second;
    }

    if (samFileName.size() == 0) {
      return alignmentScore; // TODO implement for bwa
    }

    std::ofstream samOut{samFileName};
    for (const auto &header : samHeaders) {
      samOut << header;
    }
    for (const auto &line : samAlignments) {
      samOut << line << "\n";
    }
    samOut.close();
    std::cout << "Wrote sam data to " << samFileName << std::endl;
  }

  return alignmentScore;
}

// Parse cs tag to count variants
int count_variants_from_cs(const std::string &cs_str, const int64_t mask_left,
                           const int64_t mask_right) {
  int var_count = 0;
  size_t pos = 0;
  int64_t seq_pos = 0;       // Track position in sequence
  int64_t total_seq_len = 0; // Track total sequence length

  // First pass to calculate total sequence length
  size_t tmp_pos = 0;
  while (tmp_pos < cs_str.length()) {
    char op = cs_str[tmp_pos++];
    switch (op) {
    case '=': // Identical sequence
      while (tmp_pos < cs_str.length() && isalpha(cs_str[tmp_pos])) {
        total_seq_len++;
        tmp_pos++;
      }
      break;
    case '*': // Substitution
      total_seq_len++;
      tmp_pos += 2; // Skip the two bases
      break;
    case '-': // Deletion
      while (tmp_pos < cs_str.length() && isalpha(cs_str[tmp_pos])) {
        total_seq_len++;
        tmp_pos++;
      }
      break;
    case '+': // Insertion
      while (tmp_pos < cs_str.length() && isalpha(cs_str[tmp_pos]))
        tmp_pos++;
      break;
    case '~': // Intron
      while (tmp_pos < cs_str.length() &&
             (isalnum(cs_str[tmp_pos]) || isalpha(cs_str[tmp_pos]))) {
        if (isalpha(cs_str[tmp_pos]))
          total_seq_len++;
        tmp_pos++;
      }
      break;
    }
  }

  // Second pass to count variants with proper masking
  while (pos < cs_str.length()) {
    char op = cs_str[pos];
    pos++; // Move past operation character

    switch (op) {
    case '=': // Identical sequence
      while (pos < cs_str.length() && isalpha(cs_str[pos])) {
        seq_pos++;
        pos++;
      }
      break;

    case '*': // Substitution
      // Only count if not in masked regions
      if (seq_pos >= mask_left &&
          (mask_right == 0 || seq_pos < total_seq_len - mask_right)) {
        var_count++;
      }
      seq_pos++; // Advance sequence position
      pos += 2;  // Skip the two bases
      break;

    case '+': // Insertion
    {
      // Count bases in insertion
      int64_t start_pos = seq_pos;
      while (pos < cs_str.length() && isalpha(cs_str[pos])) {
        pos++;
      }
      // Only count if not in masked regions
      if (start_pos >= mask_left &&
          (mask_right == 0 || start_pos < total_seq_len - mask_right)) {
        var_count++;
      }
    } break;

    case '-': // Deletion
    {
      // Count bases in deletion
      int64_t start_pos = seq_pos;
      while (pos < cs_str.length() && isalpha(cs_str[pos])) {
        seq_pos++;
        pos++;
      }
      // Only count if not in masked regions
      if (start_pos >= mask_left &&
          (mask_right == 0 || start_pos < total_seq_len - mask_right)) {
        var_count++;
      }
    } break;

    case '~': // Intron - format: ~gt12ag
    {
      int64_t start_pos = seq_pos;
      // Skip splice signal and length
      while (pos < cs_str.length() &&
             (isalnum(cs_str[pos]) || isalpha(cs_str[pos]))) {
        if (isalpha(cs_str[pos]))
          seq_pos++;
        pos++;
      }
      // Only count if not in masked regions
      if (start_pos >= mask_left &&
          (mask_right == 0 || start_pos < total_seq_len - mask_right)) {
        var_count++;
      }
    } break;
    }
  }
  return var_count;
}
