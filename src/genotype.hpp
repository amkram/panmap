#pragma once
#include <iostream>
#include "tree.hpp"
#include "PangenomeMAT.hpp"

namespace genotype {
    using namespace std;
    using namespace tree;

    struct VariationSite {
        VariationSite(
            size_t sid, char ref, size_t position, int variation_types, const string& nucs,
            const vector<string>& insertion_seqs, const vector<string>& deletion_seqs, const string& errors,
            const mutationMatrices& mutMat
        );

        size_t site_id;
        size_t ref_position;
        int8_t site_info; // 2 bit -> reference nuc, 3 bit -> varaition types
        
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

    pair< vector<VariationSite>, pair<size_t, size_t> > getVariantSites(std::ifstream& fin, const mutationMatrices& mutMat);
    void printSamplePlacementVCF(std::ifstream& fin, mutationMatrices& mutMat, bool variantOnly, size_t maskSize);
}