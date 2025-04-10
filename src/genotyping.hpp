#pragma once

#include "htslib/sam.h"
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include "conversion.hpp"
#include "panman.hpp"
#include "state.hpp"
#include <iostream>

namespace genotyping {


struct VariationSite {
  VariationSite(size_t sid, char ref, size_t position, int variation_types,
                const std::string &nucs, const std::vector<std::string> &insertion_seqs,
                const std::vector<std::string> &deletion_seqs, const std::string &errors,
                const mutationMatrices &mutMat);

  size_t site_id;
  size_t ref_position;
  int8_t site_info; // 2 bit -> reference nuc, 3 bit -> varaition types
  char ref_nuc;

  // substitution
  std::vector<std::vector<double>> read_errs;

  // deletion
  std::map<std::string, std::vector<double>> deletions;

  // insertion
  std::map<std::string, std::vector<double>> insertions;

  size_t most_probable_idx;
  std::vector<double> likelihoods;
  std::vector<double> posteriors;
  std::vector<size_t> read_depth;
};

// Ensure correct usage of std::getline
  void stringSplit(const std::string& str, char delimiter, std::vector<std::string>& out);

  // Ensure correct initialization of vectors and maps
  struct mutationMatrices {
    std::vector<std::vector<double>> submat;
    std::unordered_map<int64_t, double> insmat;
    std::unordered_map<int64_t, double> delmat;
    double maxInsLogProb;
    double maxDelLogProb;
    bool filled;

    mutationMatrices() : submat(4, std::vector<double>(4, 0.0)), maxInsLogProb(100.0), maxDelLogProb(100.0), filled(false) {}
  };

  /**
   * Fill mutation matrices from a file
   * 
   * @param mutMat Mutation matrices to fill
   * @param inf Input file stream to read from
   */
  void fillMutationMatricesFromFile(mutationMatrices &mutMat, std::ifstream &inf);

  /**
   * Build mutation matrices helper function
   * 
   * @param mutMat Mutation matrices to fill
   * @param tree Tree to process
   * @param node Current node being processed
   * @param stateManager StateManager containing sequence data
   * @param parentBaseCounts Parent base counts
   * @param totalBaseCounts Total base counts
   * @param subCount Substitution counts
   * @param insCount Insertion counts
   * @param delCount Deletion counts
   */
  void buildMutationMatricesHelper(
      mutationMatrices &mutMat,
      panmanUtils::Tree *tree,
      panmanUtils::Node *node,
      state::StateManager &stateManager,
      std::vector<int64_t> &parentBaseCounts,
      std::vector<int64_t> &totalBaseCounts,
      std::vector<std::vector<int64_t>> &subCount,
      std::unordered_map<int64_t, int64_t> &insCount,
      std::unordered_map<int64_t, int64_t> &delCount);

  /**
   * Fill mutation matrices from a tree
   * 
   * @param mutMat Mutation matrices to fill
   * @param tree Tree to process
   * @param path File path to save matrices to
   */
  void fillMutationMatricesFromTree_test(
      mutationMatrices &mutMat, 
      panmanUtils::Tree *tree, 
      const std::string& path);

std::vector<std::vector<double>>
scaleMutationSpectrum(const mutationMatrices &mutMat, double mutationRate);
std::string
applyMutationSpectrum(const std::string &line,
                      const std::vector<std::vector<double>> &scaled_submat);
std::pair<std::vector<VariationSite>, std::pair<size_t, size_t>>
getVariantSites(std::istream &fin, const mutationMatrices &mutMat);
void printSamplePlacementVCF(std::istream &fin, const mutationMatrices &mutMat,
                             bool variantOnly, size_t maskSize,
                             std::ofstream &fout);

/**
 * @brief Genotype variants from alignment data
 * @param prefix Output file prefix
 * @param refFileName Reference FASTA file
 * @param bestMatchSequence Best matching sequence from placement
 * @param bamFileName BAM alignment file
 * @param mpileupFileName MPileup output file
 * @param vcfFileName VCF output file
 * @param mutMat Mutation matrices
 */
void genotype(std::string prefix, std::string refFileName,
              std::string bestMatchSequence, std::string bamFileName,
              std::string mpileupFileName, std::string vcfFileName,
              mutationMatrices &mutMat);

} // namespace genotyping

inline void extractSam(const std::string &samFileName,
                       std::vector<char *> &samAlignments,
                       std::string &samHeader) {
  std::ifstream samFile{samFileName};
  if (!samFile.is_open()) {
    std::cerr << "Error: failed to open file " << samFileName << std::endl;
    return;
  }
  std::string line;
  while (std::getline(samFile, line)) {
    if (line.rfind("@", 0) == 0) {
      if (line.rfind("@SQ", 0) == 0) {
        samHeader += line;
      }
    } else {
      char *cstr = new char[line.length() + 1];
      std::strcpy(cstr, line.c_str());
      samAlignments.push_back(cstr);
    }
  }
  samFile.close();
}

inline void extractRef(const std::string &refFileName,
                       std::string &bestMatchSequence) {
  std::ifstream refFile{refFileName};
  if (!refFile.is_open()) {
    std::cerr << "Error: failed to open file " << refFileName << std::endl;
    return;
  }
  std::string line;
  while (std::getline(refFile, line)) {
    if (line[0] != '>') {
      if (!line.empty() && line.back() == '\n') {
        line.pop_back();
      }
      bestMatchSequence += line;
    }
  }
  refFile.close();
}

inline void callVariantsFromSAM(const std::string &samFileName,
                                const std::string &refFileName,
                                genotyping::mutationMatrices &mutMat,
                                const std::string &prefix) {
  std::vector<char *> samAlignments;
  std::string samHeader;
  std::string bestMatchSequence;

  extractSam(samFileName, samAlignments, samHeader);
  extractRef(refFileName, bestMatchSequence);
  std::string current_prefix = prefix;
  std::string current_refFileName = refFileName;

  // Convert to BAM
  sam_hdr_t *header;
  bam1_t **bamRecords;
  std::string bamFileName = current_prefix + ".bam";
  createBam(samAlignments, samHeader, bamFileName, header, bamRecords);

  std::string mpileupFileName = current_prefix + ".mpileup";
  // //Convert to Mplp
  // char *mplpString;
  // createMplp(
  //     bestMatchSequence,
  //     header,
  //     bamRecords,
  //     samAlignments.size(),
  //     mpileupFileName,

  //     mplpString
  // );
  createMplpBcf(current_prefix, current_refFileName, bestMatchSequence,
                bamFileName, mpileupFileName);

  std::string vcfFileName = current_prefix + ".vcf";
  // //Convert to VCF
  // createVcf(
  //     mplpString,
  //     mutMat,
  //     vcfFileName,
  //     true
  // );

  createVcfWithMutationMatrices(current_prefix, mpileupFileName, mutMat,
                                vcfFileName, 0.0011);
} // namespace genotyping