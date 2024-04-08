#include <algorithm>
#include "place.hpp"
#include "pmi.hpp"
#include "util.hpp"
#include "tree.hpp"
#include "genotype.hpp"
#include <cmath>
#include <htslib/sam.h>

#include "conversion.hpp"

using namespace PangenomeMAT;
using namespace tree;
using namespace seeding;
using namespace util;
using namespace genotype; //Remove

void mutateSeedmerMap(std::unordered_map<int32_t, std::string> &seedmers, const std::string &nid, seedmerIndex_t &index, std::vector<std::pair<int32_t, std::string>> &seedmersToRevert) {
    // std::cout << "Forward => " << nid << "\n";
    for (const auto &op : index[nid]) {
        std::string seedmer = op.second;
        int32_t pos = op.first;
        if (seedmer[0] == '@') { // deletion
            if (seedmers.find(pos) != seedmers.end()) {
                seedmers.erase(seedmers.find(pos));
            }
        } else { // insertion
            seedmersToRevert.push_back({pos, seedmers[pos]});
            seedmers[pos] = seedmer;
        }
    }
}
void revertSeedmerMap(std::unordered_map<int32_t, std::string> &seedmers, const std::string &nid, seedmerIndex_t &index, std::vector<std::pair<int32_t, std::string>> &seedmersToRevert) {
    for (auto it = index[nid].rbegin(); it != index[nid].rend(); it++) {
        std::string seedmer = it->second;
        int32_t pos = it->first;
        if (seedmer[0] == '@') { // deletion
            seedmers[pos] = seedmer.substr(1);
        } else {
            if (seedmers.find(pos) != seedmers.end()) {
                seedmers.erase(seedmers.find(pos));
            }
        }
    }
    for (const auto &p : seedmersToRevert) {
        if (p.second.size() > 1) {
            seedmers[p.first] = p.second;
        }
    }
}

void getPhyloCounts(seedmerIndex_t &index, std::unordered_map<std::string, int32_t> &counts, std::unordered_map<int32_t, std::string> &seedmers, Tree *T, Node *node) {
    std::vector<std::pair<int32_t, std::string>> seedmersToRevert;
    mutateSeedmerMap(seedmers, node->identifier, index, seedmersToRevert);
    std::unordered_map<std::string, bool> seen;
    for (const auto &seed : seedmers) {
        if (seen.find(seed.second) != seen.end()) {
            continue;
        }
        seen[seed.second] = true;
        counts[seed.second] += 1;
    }
    for (Node *child : node->children) {
        getPhyloCounts(index, counts, seedmers, T, child);
    }
    revertSeedmerMap(seedmers, node->identifier, index, seedmersToRevert);
}

void placeDFS(std::ofstream *out, std::unordered_map<std::string, std::set<std::pair<int32_t, std::string>, decltype(seed_cmp)>> &index, std::unordered_map<int32_t, std::string> &seedmers, std::map<std::string, float> &scores, std::unordered_map<std::string, int32_t> &seedmerCounts, std::unordered_map<std::string, int32_t> &phyloCounts, const Node *node, Tree *T, const std::string &optionalTarget, std::unordered_map<int32_t, std::string> *optionalOutputSeedmers) {
    std::vector<std::pair<int32_t, std::string>> seedmersToRevert;

    mutateSeedmerMap(seedmers, node->identifier, index, seedmersToRevert);

    if (optionalTarget != "" && optionalTarget == node->identifier) {
        *optionalOutputSeedmers = seedmers;
        return;
    }
    double score = 0;
    double score_inf = 0;
    double score_count = 0;
    double score_other = 0;
 //   if (node->identifier != T->root->identifier) {
    if(1) {
        std::unordered_map<std::string, int32_t> nodeSeedmerCounts;
        for (const auto &seed : seedmers) {
            if (seed.second == "") {
                continue;
            }
            if (nodeSeedmerCounts.find(seed.second) == nodeSeedmerCounts.end()) {
                nodeSeedmerCounts[seed.second] = 1;
            } else {
                nodeSeedmerCounts[seed.second] += 1;
            }
        }
        int32_t max = 0;
        for (const auto &seed : seedmers) {
            double informativity = 1.0 - phyloCounts[seed.second] / T->allNodes.size(); 
            if (seedmerCounts.find(seed.second) != seedmerCounts.end() && seedmers.size() > 10) {
                if (seedmerCounts[seed.second] < 2) {
                    continue;
                }
                // it's in the reads
                int32_t num = seedmerCounts[seed.second]; 
                if (num > max) {
                    max = num;
                }
                score += informativity * std::min(1.0, std::log(num)) / (double) nodeSeedmerCounts[seed.second];
                score_inf += informativity;
                score_count += 1;
                score_other += std::min(1.0, std::log(num)) / (double) nodeSeedmerCounts[seed.second];
            }
        }
        score /= (double) seedmers.size();
        score_other /= (double) seedmers.size();
        score_count /= (double) seedmers.size();
        score_inf /= (double) seedmers.size();
        scores[node->identifier] = score_inf;
    } 
    if (out != nullptr) {
        *out << node->identifier << "\t" << score << "\t" << score_inf << "\t" << score_count << "\t" << score_other << "\n";
    }
    for (Node *child : node->children) {
        placeDFS(out, index, seedmers, scores, seedmerCounts, phyloCounts, child, T, optionalTarget, optionalOutputSeedmers);
    }
    revertSeedmerMap(seedmers, node->identifier, index, seedmersToRevert);
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
void placeHelper(std::unordered_map<std::string, std::set<std::pair<int32_t, std::string>, decltype(seed_cmp)>> &index, std::unordered_map<int32_t, std::string> &seedmers, std::map<std::string, float> &scores, std::unordered_map<std::string, int32_t> &seedmerCounts, std::unordered_map<std::string, int32_t> &phyloCounts, const Node *node, Tree *T, const std::string optionalTarget) {
    std::ofstream out("placement.out");
    placeDFS(&out, index, seedmers, scores, seedmerCounts, phyloCounts, T->root, T, "", nullptr);
}

seedmerIndex_t seedsFromFastq(std::ifstream &indexFile, int32_t *k, int32_t *s, int32_t *j, std::unordered_map<std::string, int32_t> &seedmerCounts, std::vector<std::string> &readSequences, std::vector<std::string> &readQuals, std::vector<std::string> &readNames, std::vector<std::vector<seed>> &readSeedsFwd, std::vector<std::vector<seed>> &readSeedsBwd,  const std::string &fastqPath) {
    util::scopedTimer();
    seedmerIndex_t seedmerIndex;
    std::string line0;
    std::getline(indexFile, line0);
    std::vector<std::string> spltTop;
    stringSplit(line0, ' ', spltTop);
    int32_t tempK = std::stoi(spltTop[0]);
    int32_t tempS = std::stoi(spltTop[1]);
    int32_t tempJ = std::stoi(spltTop[2]);

    FILE *fp;
    kseq_t *seq;
    fp = fopen(fastqPath.c_str(), "r");
    if(!fp){
        std::cerr << "Error: File " << fastqPath << " not found" << std::endl;
        exit(0);
    }

    seq = kseq_init(fileno(fp));
    
    int line;
    while ((line = kseq_read(seq)) >= 0) {
        readSequences.push_back(seq->seq.s);
        readNames.push_back(seq->name.s);
        readQuals.push_back(seq->qual.s);
    }
    for (int i = 0; i < readSequences.size(); i++) {
        std::string seq = readSequences[i];
        std::string name = readNames[i];
        std::string rc = reverseComplement(seq);

        std::vector<seed> syncmers = syncmerize(seq, tempK, tempS, false, true, 0);
        readSeedsFwd.push_back(syncmers);

        std::vector<seed> syncmersReverse = syncmerize(rc, tempK, tempS, false, true, 0);
        readSeedsBwd.push_back(syncmersReverse);

        std::vector<seedmer> seedmers = seedmerize(syncmers, tempJ);
        std::vector<seedmer> seedmersReverse = seedmerize(syncmersReverse, tempJ);
        for (const auto &m : seedmers) {
            if (seedmerCounts.find(m.seq) == seedmerCounts.end()) {
                seedmerCounts[m.seq] = 1;
            } else {
                seedmerCounts[m.seq] += 1;
            }
        }
        for (const auto &m : seedmersReverse) {
            if (seedmerCounts.find(m.seq) == seedmerCounts.end()) {
                seedmerCounts[m.seq] = 1;
            } else {
                seedmerCounts[m.seq] += 1;
            }
        }
    }

    while (std::getline(indexFile, line0)) {
        std::vector<std::string> splt;
        stringSplit(line0, ' ', splt);
        std::string nid = splt[0];
        seedmerIndex[nid] = {};
        for (int32_t i = 1; i < splt.size(); i++) {
            std::vector<std::string> splt2;
            stringSplit(splt[i], ':', splt2);
            int32_t pos = std::stoi(splt2[0]);
            std::string seedmer = splt2[1];
            seedmerIndex[nid].insert(std::make_pair(pos, seedmer));
        }
    }
    *k = tempK;
    *s = tempS;
    *j = tempJ;
    return seedmerIndex;
}


void place::placeIsolate(std::ifstream &indexFile, const tree::mutationMatrices& mutMat, const std::string &reads1Path, const std::string &reads2Path, std::string &samFileName, std::string &bamFileName, std::string &mpileupFileName, std::string &vcfFileName, std::string &refFileName, Tree *T) {
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
    std::unordered_map<std::string, int32_t> seedmerCounts;
    int32_t k, s, j;
    seedmerIndex_t seedmerIndex = seedsFromFastq(indexFile, &k, &s, &j, seedmerCounts, readSequences, readQuals, readNames, readSeedsFwd, readSeedsBwd, reads1Path);
    
    /* Collecting forward and backward seeds into one vector */
    std::vector<std::vector<seed>> readSeeds;
    readSeeds.reserve(readSeedsFwd.size());
    for(int i = 0; i < readSeedsFwd.size(); i++) {
        std::vector<seed> thisReadsSeeds;
        thisReadsSeeds.reserve(readSeedsFwd[i].size() + readSeedsBwd[i].size());

        for(int j = 0; j < readSeedsFwd[i].size(); j++) {
            readSeedsFwd[i][j].reversed = false;
            readSeedsFwd[i][j].pos = readSeedsFwd[i][j].pos + k - 1; // Minimap standard
            thisReadsSeeds.push_back(readSeedsFwd[i][j]);
        }
        for(int j = 0; j < readSeedsBwd[i].size(); j++) {
            readSeedsBwd[i][j].reversed = true;
            readSeedsBwd[i][j].pos = readSeedsBwd[i][j].pos + k - 1; // Minimap standard
            thisReadsSeeds.push_back(readSeedsBwd[i][j]);
        }
        readSeeds.push_back(thisReadsSeeds);
    }

    /* Sample placement */
    std::cout << "⋌⋋ Placing sample ... " << std::flush;
    std::map<std::string, float> scores;
    std::unordered_map<std::string, int32_t> phyloCounts;
    std::unordered_map<int32_t, std::string> dynamicSeedmersPhylo;
    std::unordered_map<int32_t, std::string> dynamicSeedmersPlace;

    getPhyloCounts(seedmerIndex, phyloCounts, dynamicSeedmersPhylo, T, T->root);

    placeHelper(seedmerIndex, dynamicSeedmersPlace, scores, seedmerCounts, phyloCounts, T->root, T, "");
    
    /* Setup target sequence and seeds */
    std::vector<std::pair<std::string, float>> targetNodes;
    std::copy(scores.begin(), scores.end(), back_inserter<std::vector<std::pair<std::string, float>>>(targetNodes));
    std::sort(targetNodes.begin(), targetNodes.end(), [] (auto &left, auto &right) { return left.second > right.second; });
    // redo DFS with target node -> returns early with dynamicSeedmers in target node's state

    std::string bestMatch = targetNodes[0].first;
    std::string bestMatchSequence = "";
    std::string gappedSeq = T->getStringFromReference(bestMatch, true);
    std::vector<int32_t> degap;
    for (int32_t i = 0; i < gappedSeq.size(); i ++) {
        char &c = gappedSeq[i];
        degap.push_back(bestMatchSequence.size());
        if (c != '-') {
            bestMatchSequence += c;
        }
    }
    
    // path format {target}.*.fastq
    std::string targetId = reads1Path.substr(0, reads1Path.find_first_of('.')) + ".1";
    std::cerr << "\n" << targetId << "\t" << bestMatch << std::endl;
    std::unordered_map<std::string, std::vector<int32_t>> seedToRefPositions;
    std::unordered_map<int32_t, std::string> targetSeedmers;
    std::unordered_map<int32_t, std::string> dynamicSeedmersTarget;
    
    placeDFS(nullptr, seedmerIndex, dynamicSeedmersTarget, scores, seedmerCounts, phyloCounts, T->root, T, bestMatch, &targetSeedmers);

    /* Debug Print statements */
    for (const auto &seedmer : targetSeedmers) {
        if (seedmer.second == "" ) {
            continue;
        }
        int32_t refPos = seedmer.first;
        std::string seed = seedmer.second.substr(0, k);
        if (seedToRefPositions.find(seed) == seedToRefPositions.end()) {
            seedToRefPositions[seed] = {};
        }
        seedToRefPositions[seed].push_back(degap[refPos]);
    }



    //Print out reference
    if(refFileName.size() > 0){
        std::ofstream outFile{refFileName};

        if (outFile.is_open()) {
            
            outFile << ">ref\n";
            outFile << bestMatchSequence << "\n";

            std::cout << "Wrote reference fasta to " << refFileName << std::endl;
        } else {
            std::cerr << "Error: failed to write to file " << refFileName << std::endl;
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
        k,
        
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