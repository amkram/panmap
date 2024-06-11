#include <algorithm>
#include <cmath>
#include <htslib/sam.h>
#include <sys/_types/_int64_t.h>
#include "place.hpp"
#include "pmi.hpp"
#include "util.hpp"
#include "tree.hpp"
#include "genotype.hpp"
#include "conversion.hpp"

using namespace PangenomeMAT;
using namespace tree;
using namespace seeding;
using namespace util;
using namespace genotype; //Remove


void perfect_shuffle(vector<std::string>& v) {
    int n = v.size();

    vector<std::string> canvas(n);

    for (int i = 0; i < n / 2; i++) {
        canvas[i*2] = v[i];
        canvas[i*2+1] = v[i + n/2];
    }

    v = canvas;
}


std::string reverseComplement(std::string dna_sequence) {
    std::string complement = "";
    for (char c : dna_sequence) {
        switch (c) {
            case 'A': complement += 'T'; break;
            case 'T': complement += 'A'; break;
            case 'C': complement += 'G'; break;
            case 'G': complement += 'C'; break;
            default: complement += c; break;
        }
    }
    std::reverse(complement.begin(), complement.end());
    return complement;
}


// second value in pair in seedMap is the number of seed hits at the current node at the genomic position of the key
void placeHelper(std::map<int64_t, std::pair<std::string, int32_t>> &seedMap, Tree *T, Node *node, std::map<std::string, int32_t> &scores, SeedmerIndex &index, int64_t &pb_idx, std::unordered_map<std::string, int32_t> &readSeedCounts) {
    /// This node's indexed seed mutations
    NodeSeedmerMutations pb_node_mutations = index.per_node_mutations(pb_idx);
    pb_idx++;

    std::vector<std::pair<int32_t, std::pair<std::string, int32_t>>> backtrack;
    std::vector<int64_t> delSeeds;
    std::vector<std::pair<int64_t, std::string>> addSeeds;

    // Apply seed mutations
    for (int mut_i = 0; mut_i < pb_node_mutations.mutations_size(); mut_i++) {
        const auto &mut = pb_node_mutations.mutations(mut_i);
        int64_t pos = mut.pos();
        if (mut.is_deletion()) {
        backtrack.push_back({pos, {seedMap[pos]}});
        delSeeds.push_back(pos);
        } else {
        if (seedMap.find(pos) != seedMap.end()) {
            backtrack.push_back({pos, seedMap[pos]});
        } else {
            backtrack.push_back({pos, {"",0}});
        }
        addSeeds.push_back({pos, mut.seq()});
        
        }
    }
    for (const auto &p : delSeeds) {
        seedMap.erase(p);
    }
    for (const auto &p : addSeeds) {
        seedMap[p.first] = {p.second, readSeedCounts[p.second]};
    }
    
    int32_t thisNodeScore = 0;
    for (const auto &p : seedMap) {
        thisNodeScore += p.second.second;
    }
    std::cout << "Node: " << node->identifier << "\tScore: " << thisNodeScore << std::endl;

    // Recursive step
    for (Node *child : node->children) {
        placeHelper(seedMap, T, child, scores, index, pb_idx, readSeedCounts);
    }

    // Backtrack
      for (const auto &p : backtrack) {
    if (p.second.first.size() > 0) {
      seedMap[p.first].first = p.second.first;
      seedMap[p.first].second = p.second.second;

    } else {
      seedMap.erase(p.first);
    }
  }
}
void seedsFromFastq(SeedmerIndex &index, std::unordered_map<std::string, int32_t> &readSeedCounts, std::vector<std::string> &readSequences, std::vector<std::string> &readQuals, std::vector<std::string> &readNames, std::vector<std::vector<seed>> &readSeedsFwd, std::vector<std::vector<seed>> &readSeedsBwd,  const std::string &fastqPath1, const std::string &fastqPath2) {
    FILE *fp;
    kseq_t *seq;
    fp = fopen(fastqPath1.c_str(), "r");
    if(!fp){
        std::cerr << "Error: File " << fastqPath1 << " not found" << std::endl;
        exit(0);
    }
    seq = kseq_init(fileno(fp));
    int line;
    while ((line = kseq_read(seq)) >= 0) {
        readSequences.push_back(seq->seq.s);
        readNames.push_back(seq->name.s);
        readQuals.push_back(seq->qual.s);
    }
    if (fastqPath2.size() > 0) {
        fp = fopen(fastqPath2.c_str(), "r");
        if(!fp){
            std::cerr << "Error: File " << fastqPath2 << " not found" << std::endl;
            exit(0);
        }
        seq = kseq_init(fileno(fp));

        line = 0;
        int forwardReads = readSequences.size();
        while ((line = kseq_read(seq)) >= 0) {
            readSequences.push_back(reverseComplement(seq->seq.s));
            readNames.push_back(seq->name.s);
            readQuals.push_back(seq->qual.s);
        }

        if (readSequences.size() != forwardReads*2){
            std::cerr << "Error: File " << fastqPath2 << " does not contain the same number of reads as " << fastqPath1 << std::endl;
            exit(0);
        }
        
        //Shuffle reads together, so that pairs are next to eatch other
        perfect_shuffle(readSequences);
        perfect_shuffle(readNames);
        perfect_shuffle(readQuals);
    }

    for (int i = 0; i < readSequences.size(); i++) {
        std::string seq = readSequences[i];
        std::string name = readNames[i];
        std::string rc = reverseComplement(seq);

        std::vector<seed> syncmers = syncmerize(seq, index.k(), index.s(), false, true, 0);
        readSeedsFwd.push_back(syncmers);

        std::vector<seed> syncmersReverse = syncmerize(rc, index.k(), index.s(), false, true, 0);
        readSeedsBwd.push_back(syncmersReverse);

        std::vector<seedmer> seedmers = seedmerize(syncmers, index.j());
        std::vector<seedmer> seedmersReverse = seedmerize(syncmersReverse, index.j());
        for (const auto &m : seedmers) {
            if (readSeedCounts.find(m.seq) == readSeedCounts.end()) {
                readSeedCounts[m.seq] = 1;
            } else {
                readSeedCounts[m.seq] += 1;
            }
        }
        for (const auto &m : seedmersReverse) {
            if (readSeedCounts.find(m.seq) == readSeedCounts.end()) {
                readSeedCounts[m.seq] = 1;
            } else {
                readSeedCounts[m.seq] += 1;
            }
        }
    }
}

void place::placeIsolate(SeedmerIndex &index, const tree::mutationMatrices& mutMat, const std::string &reads1Path, const std::string &reads2Path, std::string &samFileName, std::string &bamFileName, std::string &mpileupFileName, std::string &vcfFileName, std::string &refFileName, Tree *T, bool use_root) {
    tree::mutableTreeData data;
    tree::globalCoords_t globalCoords;
    tree::setup(data, globalCoords, T);

    /* Read processing */
    std::cout << "\n◠◡ Processing reads ... " << std::flush;
    std::vector<std::string> readSequences;
    std::vector<std::string> readQuals;
    std::vector<std::string> readNames;
    std::vector<std::vector<seed>> readSeedsFwd;
    std::vector<std::vector<seed>> readSeedsBwd;
    std::unordered_map<std::string, int32_t> readSeedCounts;
    seedsFromFastq(index, readSeedCounts, readSequences, readQuals, readNames, readSeedsFwd, readSeedsBwd, reads1Path, reads2Path);
    
    /* Collecting forward and backward seeds into one vector */
    std::vector<std::vector<seed>> readSeeds;
    readSeeds.reserve( readSeedsFwd.size());
    for(int i = 0; i < readSeedsFwd.size(); i++) {
        std::vector<seed> thisReadsSeeds;
        thisReadsSeeds.reserve(readSeedsFwd[i].size() + readSeedsBwd[i].size());

        for(int j = 0; j < readSeedsFwd[i].size(); j++) {
            readSeedsFwd[i][j].reversed = false;
            readSeedsFwd[i][j].pos = readSeedsFwd[i][j].pos + index.k() - 1; // Minimap standard
            thisReadsSeeds.push_back(readSeedsFwd[i][j]);
        }
        for(int j = 0; j < readSeedsBwd[i].size(); j++) {
            readSeedsBwd[i][j].reversed = true;
            readSeedsBwd[i][j].pos = readSeedsBwd[i][j].pos + index.k() - 1; // Minimap standard
            thisReadsSeeds.push_back(readSeedsBwd[i][j]);
        }
        readSeeds.push_back(thisReadsSeeds);
    }

    bool pairedEndReads = reads2Path.size();

    /* Sample placement */
    std::cout << "⋌⋋ Placing sample ... " << std::flush;
    std::map<std::string, int32_t> scores;
    seedmerIndex_t seedmerIndex;
    enum class ScoringMethod {
        NUM_SEED_HITS,
        JACCARD
    };
    enum class NormalizationMethod {
        NONE,
        BY_GENOME_LEN,
        BY_PHYLO_FREQ
    };

    ScoringMethod scoringMethod = ScoringMethod::NUM_SEED_HITS;
    NormalizationMethod normalizationMethod = NormalizationMethod::NONE;
    
    std::string bestMatch;
    if ( !use_root ) {
        int64_t pb_idx = 0;
        std::map<int64_t, std::pair<std::string, int32_t>> seedMap;
        
        placeHelper(seedMap, T, T->root, scores, index, pb_idx, readSeedCounts);

        /* Setup target sequence and seeds */
        std::vector<std::pair<std::string, float>> targetNodes;
        std::copy(scores.begin(), scores.end(), back_inserter<std::vector<std::pair<std::string, float>>>(targetNodes));
        std::sort(targetNodes.begin(), targetNodes.end(), [] (auto &left, auto &right) { return left.second > right.second; });
        
        bestMatch = targetNodes[0].first;
    } else{
        bestMatch = T->root->identifier;
    }
    
    std::string bestMatchSequence = "";
    std::string gappedSeq = T->getStringFromReference(bestMatch, true);
    std::vector<int32_t> degapVec;
    for (int32_t i = 0; i < gappedSeq.size(); i ++) {
        char &c = gappedSeq[i];
        degapVec.push_back(bestMatchSequence.size());
        if (c != '-') {
            bestMatchSequence += c;
        }
    }

    // path format {target}.*.fastq
    std::string targetId = reads1Path.substr(0, reads1Path.find_first_of('.')) + ".1";
    std::cerr << "\n" << targetId << "\t" << bestMatch << std::endl;
    std::unordered_map<std::string, std::vector<int32_t>> seedToRefPositions;
    std::unordered_map<int32_t, std::string> targetSeedmers;

    for (int i = 0; i < index.per_node_mutations_size(); i++) {
        const NodeSeedmerMutations muts = index.per_node_mutations(i);
        seedmerIndex[muts.node_id()] = {};
        for (int j = 0; j < muts.mutations_size(); j++) {
            const SeedmerMutation mut = muts.mutations(j);
            seedmerIndex[muts.node_id()].insert({mut.pos(), (mut.is_deletion() ? '@' + mut.seq() : mut.seq())});
        }
    }
    

    //Create SAM
    std::vector<char *> samAlignments;
    std::string samHeader;

    createSam(
        readSeeds,
        readSequences,
        readQuals,
        readNames,
        bestMatchSequence,
        seedToRefPositions,
        samFileName,
        index.k(),
        pairedEndReads,

        samAlignments,
        samHeader
    );


    //Convert to BAM
    sam_hdr_t *header;
    bam1_t **bamRecords;

    createBam(
        samAlignments,
        samHeader,
        bamFileName,

        header,
        bamRecords
    );


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


    //Convert to VCF
    createVcf(
        mplpString,
        mutMat,
        vcfFileName
    );

}