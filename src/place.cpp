#include <algorithm>
#include "place.hpp"
#include "pmi.hpp"
#include "util.hpp"
#include "tree.hpp"
#include "mgsr.hpp"
#include "seeding.hpp"
#include "genotype.hpp"
#include <cmath>
#include <htslib/sam.h>

#include "conversion.hpp"
#include <tbb/parallel_for.h>
#include <tbb/concurrent_vector.h>


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

void perfect_shuffle(vector<std::string>& v) {
    int n = v.size();

    vector<std::string> canvas(n);

    for (int i = 0; i < n / 2; i++) {
        canvas[i*2] = v[i];
        canvas[i*2+1] = v[i + n/2];
    }

    v = canvas;
}

seedmerIndex_t seedsFromFastq(std::ifstream &indexFile, int32_t *k, int32_t *s, int32_t *j, std::unordered_map<std::string, int32_t> &seedmerCounts, std::vector<std::string> &readSequences, std::vector<std::string> &readQuals, std::vector<std::string> &readNames, std::vector<std::vector<seed>> &readSeedsFwd, std::vector<std::vector<seed>> &readSeedsBwd,  const std::string &fastqPath1, const std::string &fastqPath2) {
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


void place::placeIsolate(std::ifstream &indexFile, const tree::mutationMatrices& mutMat, const std::string &reads1Path, const std::string &reads2Path, const std::string& prefix, const bool& makeSam, const bool& makeBam, const bool& makeMPileup, const bool& makeVCF, const bool& makeRef, Tree *T, bool use_root) {
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
    

    seedmerIndex_t seedmerIndex = seedsFromFastq(indexFile, &k, &s, &j, seedmerCounts, readSequences, readQuals, readNames, readSeedsFwd, readSeedsBwd, reads1Path, reads2Path);
    
    /* Collecting forward and backward seeds into one vector */
    std::vector<std::vector<seed>> readSeeds;
    readSeeds.reserve( readSeedsFwd.size());
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


    bool pairedEndReads = reads2Path.size();


    /* Sample placement */
    std::cout << "⋌⋋ Placing sample ... " << std::flush;
    std::map<std::string, float> scores;
    std::unordered_map<std::string, int32_t> phyloCounts;
    std::unordered_map<int32_t, std::string> dynamicSeedmersPhylo;
    std::unordered_map<int32_t, std::string> dynamicSeedmersPlace;
    
    std::string bestMatch;
    if ( !use_root ) {
        
        getPhyloCounts(seedmerIndex, phyloCounts, dynamicSeedmersPhylo, T, T->root);
        
        placeHelper(seedmerIndex, dynamicSeedmersPlace, scores, seedmerCounts, phyloCounts, T->root, T, "");
        
        /* Setup target sequence and seeds */
        std::vector<std::pair<std::string, float>> targetNodes;
        std::copy(scores.begin(), scores.end(), back_inserter<std::vector<std::pair<std::string, float>>>(targetNodes));
        std::sort(targetNodes.begin(), targetNodes.end(), [] (auto &left, auto &right) { return left.second > right.second; });
        // redo DFS with target node -> returns early with dynamicSeedmers in target node's state
        
        bestMatch = targetNodes[0].first;
    }else{
        bestMatch = T->root->identifier;
    }
    
    
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
    exit(0);
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
    std::string refFileName = prefix + ".fa";
    if(makeRef){
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

    std::string samFileName = prefix + ".sam";
    createSam(
        readSeeds,
        readSequences,
        readQuals,
        readNames,
        bestMatchSequence,
        seedToRefPositions,
        makeSam,
        samFileName,
        k,
        pairedEndReads,
        
        samAlignments,
        samHeader
    );


    //Convert to BAM
    sam_hdr_t *header;
    bam1_t **bamRecords;

    std::string bamFileName = prefix + ".bam";
    createBam(
        samAlignments,
        samHeader,
        makeBam,
        bamFileName,

        header,
        bamRecords
    );


    //Convert to Mplp
    char *mplpString;

    std::string mpileupFileName = prefix + ".mpileup";
    createMplp(
        bestMatchSequence,
        header,
        bamRecords,
        samAlignments.size(),
        makeMPileup,
        mpileupFileName,

        mplpString
    );


    //Convert to VCF

    std::string vcfFileName = prefix + ".vcf";
    createVcf(
        mplpString,
        mutMat,
        makeVCF,
        vcfFileName
    );

}

void place::placeAccio(
    PangenomeMAT::Tree *T, const tree::mutationMatrices& mutMat, const int32_t& accioK, const int32_t& accioS, const int32_t& accioL,
    const std::string& defaultKmiPath, const std::string& reads1File, const std::string& reads2File, const std::string& prefix,
    const bool& makeSam, const bool& makeBam, const bool& makeMPileup, const bool& makeVCF, const bool& makeRef, const int& maximumGap,
    const int& minimumCount, const int& minimumScore, const double& errorRate, const bool& confidence, const int32_t& roundsRemove, const double& removeThreshold
    ) {
    std::unordered_map<std::string, tbb::concurrent_vector<std::pair<int32_t, double>>> allScores;
    std::vector<std::pair<int32_t, std::vector<size_t>>> numReadDuplicates;
    std::unordered_map<std::string, std::string> leastRecentIdenticalAncestor;
    std::unordered_map<std::string, std::unordered_set<std::string>> identicalSets;
    std::vector<std::string> readSequences;
    std::vector<std::string> readQuals;
    std::vector<std::string> readNames;
    std::vector<std::vector<seeding::seed>> readSeeds;
    std::ifstream seedmersIndex(defaultKmiPath);
    int32_t numReads;
    mgsr::scorePseudo(seedmersIndex, reads1File, reads2File, allScores, numReadDuplicates, leastRecentIdenticalAncestor, identicalSets, numReads, T, readSeeds, readSequences, readQuals, readNames, maximumGap, minimumCount, minimumScore, errorRate);       
    Eigen::MatrixXd probs;
    Eigen::VectorXd props;
    std::vector<std::string> nodes;
    double llh;
    mgsr::squaremHelper(T, allScores, numReadDuplicates, numReads, leastRecentIdenticalAncestor, identicalSets, probs, nodes, props, llh, roundsRemove, removeThreshold, "");

    std::vector<std::pair<std::string, double>> sortedOut(nodes.size());
    for (size_t i = 0; i < nodes.size(); ++i) {
        sortedOut.at(i) = {nodes[i], props(i)};
    }
    std::sort(sortedOut.begin(), sortedOut.end(), [](const std::pair<std::string, double>& a, const std::pair<std::string, double>& b) {
        return a.second > b.second;
    });

    std::vector<double> llhdiffs(nodes.size());
    std::vector<double> exclllhs(nodes.size());
    if (confidence) {
        tbb::parallel_for(tbb::blocked_range<size_t>(0, nodes.size()), [&](const tbb::blocked_range<size_t>& range) {
            for (size_t i = range.begin(); i < range.end(); ++i) {
                Eigen::MatrixXd curprobs;
                Eigen::VectorXd curprops;
                std::vector<std::string> curnodes;
                double curllh;
                mgsr::squaremHelper(T, allScores, numReadDuplicates, numReads, leastRecentIdenticalAncestor, identicalSets, curprobs, curnodes, curprops, curllh, roundsRemove, removeThreshold, nodes[i]);
                exclllhs.at(i) = curllh;
                llhdiffs.at(i) = llh - curllh;
            }
        });
    }

    std::unordered_map<std::string, std::unordered_set<size_t>> assignedReads;
    std::cerr << "\nStarting reads assignment" << std::endl;
    mgsr::accio(allScores, nodes, probs, props, numReadDuplicates, assignedReads);
    std::unordered_map<std::string, std::pair<double, double>> readAssignmentAccuracy = mgsr::getReadAssignmentAccuracy(assignedReads, nodes, readNames, leastRecentIdenticalAncestor);
    std::cerr << "Finished reads assignment" << std::endl;

    std::string abundanceOutFile = prefix + ".abundance";
    std::ofstream abundanceOut(abundanceOutFile);
    abundanceOut << "@likelihood: " << llh << "\n";
    for (size_t i = 0; i < sortedOut.size(); ++i) {
        const auto& node = sortedOut[i];
        abundanceOut << node.first;
        if (identicalSets.find(node.first) != identicalSets.end()) {
            for (const auto& identicalNode : identicalSets.at(node.first)) {
                abundanceOut << "," << identicalNode;
            }
        }
        abundanceOut << "\t"
                  << node.second << "\t"
                  << assignedReads[node.first].size() << "\t"
                  << readAssignmentAccuracy[node.first].first << "\t"
                  << readAssignmentAccuracy[node.first].second;
        // std::cout << "\n";

        /* for assessing read assignment accuracy*/
        abundanceOut << "\t" << readAssignmentAccuracy[node.first].first / readAssignmentAccuracy[node.first].second << "\n";
    }

    std::cerr << "Wrote abundance estimates and reads assignment data to " << abundanceOutFile << std::endl;

    bool pairedEndReads = reads2File.size();
    for (const auto& node : sortedOut) {
        const std::string& nodeName = node.first;
        std::string nodeSeq = T->getStringFromReference(nodeName, false);
        std::vector<std::vector<seeding::seed>> assignedReadSeeds;
        std::vector<std::string> assignedReadSequences;
        std::vector<std::string> assignedReadQuals;
        std::vector<std::string> assignedReadNames;
        assignedReadSeeds.reserve(assignedReads[nodeName].size());
        for (const size_t& readIdx : assignedReads[nodeName]) {
            assignedReadSeeds.push_back(readSeeds[readIdx]);
            assignedReadSequences.push_back(readSequences[readIdx]);
            assignedReadQuals.push_back(readQuals[readIdx]);
            assignedReadNames.push_back(readNames[readIdx]);
        }
        std::unordered_map<std::string, std::vector<int32_t>> seedToRefPositions;
        for (size_t i = 0; i < nodeSeq.size() - accioK + 1; ++i) {
            std::string kmer = nodeSeq.substr(i, accioK);
            if (seeding::is_syncmer(kmer, accioS, false)) {
                seedToRefPositions[kmer].push_back(i);
            }
        }

        if (!(makeSam || makeBam || makeMPileup || makeVCF)) {
            return;
        }
        //Create SAM
        std::vector<char *> samAlignments;
        std::string samHeader;

        std::string samFileName = prefix + "_" + nodeName + ".sam";
        createSam(
            assignedReadSeeds,
            assignedReadSequences,
            assignedReadQuals,
            assignedReadNames,
            nodeSeq,
            seedToRefPositions,
            makeSam,
            samFileName,
            accioK,
            pairedEndReads,
            samAlignments,
            samHeader
        );


        //Convert to BAM
        sam_hdr_t *header;
        bam1_t **bamRecords;

        std::string bamFileName = prefix + "_" + nodeName + ".bam";
        createBam(
            samAlignments,
            samHeader,
            makeBam,
            bamFileName,

            header,
            bamRecords
        );

        //Convert to Mplp
        char *mplpString;

        std::string mpileupFileName = prefix + "_" + nodeName + ".mpileup";
        createMplp(
            nodeSeq,
            header,
            bamRecords,
            samAlignments.size(),
            makeMPileup,
            mpileupFileName,

            mplpString
        );


        //Convert to VCF
        std::string vcfFileName = prefix + "_" + nodeName + ".vcf";
        createVcf(
            mplpString,
            mutMat,
            makeVCF,
            vcfFileName
        );
    }
}