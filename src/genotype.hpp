#ifndef __GENOTYPE_HPP
#define __GENOTYPE_HPP

#pragma once
#include <iostream>
#include "seed_annotated_tree.hpp"
#include "conversion.hpp"

namespace genotype {
    using namespace std;
    using namespace seed_annotated_tree;

    struct VariationSite {
        VariationSite(
            size_t sid, char ref, size_t position, int variation_types, const string& nucs,
            const vector<string>& insertion_seqs, const vector<string>& deletion_seqs, const string& errors,
            const mutationMatrices& mutMat
        );

        size_t site_id;
        size_t ref_position;
        int8_t site_info; // 2 bit -> reference nuc, 3 bit -> varaition types
        char ref_nuc;
        
        // substitution
        vector< vector<double> > read_errs;

        // deletion
        // map<size_t, size_t> deletions;
        map<string, vector<double> > deletions;
        
        // insertion
        // map<string, size_t> insertions;
        map<string, vector<double> > insertions;

        size_t most_probable_idx;
        vector<double> likelihoods;
        vector<double> posteriors;
        vector<size_t> read_depth;
    };
  
    vector<std::vector<double>> scaleMutationSpectrum(const mutationMatrices& mutMat, double mutationRate);
    std::string applyMutationSpectrum(const std::string& line, const std::vector<std::vector<double>>& scaled_submat);
    pair< vector<VariationSite>, pair<size_t, size_t> > getVariantSites(std::istream& fin, const mutationMatrices& mutMat);
    void printSamplePlacementVCF(std::istream& fin, const mutationMatrices& mutMat, bool variantOnly, size_t maskSize, std::ofstream& fout);
}


inline void extractSam(const std::string& samFileName, std::vector<char *>& samAlignments, std::string& samHeader){
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

inline void extractRef(const std::string& refFileName, std::string& bestMatchSequence){
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

inline void callVariantsFromSAM(
  const std::string& samFileName, const std::string& refFileName,
  seed_annotated_tree::mutationMatrices& mutMat, const std::string& prefix
) {
  std::vector<char *> samAlignments;
  std::string samHeader;
  std::string bestMatchSequence;

  extractSam(samFileName, samAlignments, samHeader);
  extractRef(refFileName, bestMatchSequence);
  

  //Convert to BAM
  sam_hdr_t *header;
  bam1_t **bamRecords;
  std::string bamFileName = prefix + ".bam";
  createBam(
      samAlignments,
      samHeader,
      bamFileName,
      header,
      bamRecords
  );

  std::string mpileupFileName = prefix + ".mplp";
  //Convert to Mplp
  char *mplpString;
  createMplp(
      bestMatchSequence,
      header,
      bamRecords,
      samAlignments.size(),
      mpileupFileName,

      mplpString
  );

  std::string vcfFileName = prefix + ".vcf";
  //Convert to VCF
  createVcf(
      mplpString,
      mutMat,
      vcfFileName,
      true
  );
}

#endif