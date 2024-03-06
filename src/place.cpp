#include <algorithm>
#include "place.hpp"
#include "pmi.hpp"
#include "util.hpp"
#include <cmath>
#include <htslib/sam.h>

using namespace PangenomeMAT;
using namespace tree;
using namespace seeding;
using namespace util;

void mutateSeedmerMap(std::unordered_map<int32_t, std::string> &seedmers, std::string nid, seedmerIndex_t &index, std::vector<std::pair<int32_t, std::string>> &seedmersToRevert) {
    // std::cout << "Forward => " << nid << "\n";
    for (const auto &op : index[nid]) {
        std::string seedmer = op.second;
        if (seedmer == "" || seedmer == "@") {
            continue;
        }
        int32_t pos = op.first;
        if (seedmer[0] == '@') { // deletion
            if (seedmers.find(pos) != seedmers.end()) {
                seedmers.erase(seedmers.find(pos));
            }
        } else {
            seedmersToRevert.push_back({pos, seedmers[pos]});
            seedmers[pos] = seedmer;
        }
    }
}
void revertSeedmerMap(std::unordered_map<int32_t, std::string> &seedmers, std::string nid, seedmerIndex_t &index, std::vector<std::pair<int32_t, std::string>> &seedmersToRevert) {
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
        if (seed.second == "") {
            continue;
        }
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

void placeDFS(std::ofstream *out, std::unordered_map<std::string, std::set<std::pair<int32_t, std::string>, decltype(seed_cmp)>> &index, std::unordered_map<int32_t, std::string> &seedmers, std::map<std::string, float> &scores, std::unordered_map<std::string, int32_t> &seedmerCounts, std::unordered_map<std::string, int32_t> &phyloCounts, const Node *node, Tree *T, const std::string optionalTarget, std::unordered_map<int32_t, std::string> *optionalOutputSeedmers) {
    std::vector<std::pair<int32_t, std::string>> seedmersToRevert;

    mutateSeedmerMap(seedmers, node->identifier, index, seedmersToRevert);

    if (optionalTarget != "" && optionalTarget == node->identifier) {
        *optionalOutputSeedmers = seedmers;
        return;
    }
    float score = 0;
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
  //  std::cout << node->identifier << "\t" << seedmers.size() << std::endl;
      double covered = 0;
      int32_t occ = 0;
      std::unordered_map<std::string, bool> seen;
    int32_t max = 0;
    for (const auto &seed : seedmers) {
         if (seed.second == "" || seed.first < 100) {
            continue;
        }
        if (seen.find(seed.second) != seen.end()) {
            continue;
        }
        if (seedmerCounts.find(seed.second) != seedmerCounts.end() && seedmers.size() > 10) {
            int32_t num = seedmerCounts[seed.second]; 
            if (num > max) {
                max = num;
            }
            occ += std::min(1.0, std::log(seedmerCounts[seed.second])) / (double) nodeSeedmerCounts[seed.second];
        }
    }
    score = occ / (double) seedmers.size();
    scores[node->identifier] = score;
    if (out != nullptr) {
        *out << node->identifier << "\t" << score << "\n";
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
    util::scopedTimer();
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

extern "C" {
    void align_reads(const char *reference, int n_reads, const char **reads, const char **quality, const char **read_names, int *r_lens, int *seed_counts, uint8_t **reversed, int **ref_positions, int **qry_positions, char** sam_alignments, int syncmer_k);
}


//TODO FIXME DELETE
void findSyncmers(const std::string& sequence, int k, int s, std::unordered_map<std::string, std::vector<int32_t>> &refSeeds) {
        std::unordered_set<std::string> syncmers;
        int len = sequence.length();

        for (int i = 0; i <= len - k; ++i) {
            std::string kmer = sequence.substr(i, k);
            std::string start_smer = kmer.substr(0, s);
            std::string end_smer = kmer.substr(k - s, s);

            bool start_is_min = true, end_is_min = true;

            for (int j = 1; j < k - 1; ++j) {
                std::string smer = kmer.substr(j, s);
                if (start_smer > smer) {
                    start_is_min = false;
                }
                if (end_smer > smer) {
                    end_is_min = false;
                }
            }

            if (start_is_min || end_is_min) {
                refSeeds[kmer] = {i};
            }
        }
}



void place::placeIsolate(std::ifstream &indexFile, const std::string &reads1Path, const std::string &reads2Path, Tree *T) {
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

    for (const auto &seedmer : targetSeedmers) {
        if (seedmer.second == "" ) {
            continue;
        }
        int32_t refPos = seedmer.first;
        std::string seed = seedmer.second.substr(0, k);
        //std::cerr << seed << " R\n";
        if (seedToRefPositions.find(seed) == seedToRefPositions.end()) {
            seedToRefPositions[seed] = {};
        }
        seedToRefPositions[seed].push_back(degap[refPos]);
    }
    

    /* Alignment to target */

    /*  @nico: at this point,
    **    - bestMatch should contain the target node id
    **    - bestMatchSequence has the target node's sequence without gaps
    **    - readSequences, readQuals, readNames contain the relevant info
    **    - seedToRefPositions maps a seed sequence (just seeds not seed-mers) to all matching positions in bestMatchSequence
    */
    std::cerr << "\n" << bestMatchSequence << std::endl;

    //TODO FIXME quickhacks                                        //
    seedToRefPositions = {};                                       //
    findSyncmers(bestMatchSequence, k, s, seedToRefPositions);     //
    
    std::vector<seed> refSeeds;
    for(auto kv : seedToRefPositions) {
        seed thisSeed;
        thisSeed.seq = kv.first;
        thisSeed.pos = kv.second[0];

        refSeeds.push_back(thisSeed);
    }
    
    //Sorting ref seeds
    sort(refSeeds.begin(), refSeeds.end(), []( const seed& lhs, const seed& rhs )
    {
        return lhs.seq < rhs.seq ;
    });

    //Sorting read seeds
    for(int i = 0; i < readSeeds.size(); i++){
        sort(readSeeds[i].begin(), readSeeds[i].end(), []( const seed& lhs, const seed& rhs )
        {
            return lhs.seq < rhs.seq;
        });
    }
    /* Debug Print statements */
    
    //std::cerr << bestMatchSequence << "\n\n\n\n\n";
    

    /*

    std::cerr << "READSEEDS\n";
    for(int j = 0; j < refSeeds.size() ; j++){
        std::cerr << "B\n";
        std::cerr << refSeeds[j].seq << " ";
        std::cerr << refSeeds[j].pos << "\n";
    }
    /*

    /*
    
    for(int i = 0; i < readSequences.size(); i++) {
        std::cerr << "\n\n" << i <<"\n";
        std::cerr << readSequences[i] << "\n";
        std::cerr << readQuals[i] << "\n";
        std::cerr << readSeeds[i].size() << "\n";
        for(int j = 0; j < readSeeds[i].size() ; j++){
            std::cerr << "X\n";
            std::cerr << readSeeds[i][j].seq.size() << "\n";
            std::cerr << readSeeds[i][j].seq << "\n";
            std::cerr << readSeeds[i][j].pos << "\n";
            std::cerr << readSeeds[i][j].reversed << "\n";
        }
    }

    
    for(int i = 0; i < readSequences.size(); i++) {
        std::cerr << "\n\n" << i << "\n";
        for(int j = 0; j < readSeeds[i].size() ; j++){
            if(!readSeeds[i][j].reversed){
                std::cerr << readSeeds[i][j].seq << "\n";
            }
        }
    }
    */
    


    //Finding syncmer matches
    for(int r = 0; r < readSequences.size() ; r++){

        int refSeedIndex = 0;
        int readSeedIndex = 0;

        std::vector<seed> matchingSeeds;
        while(readSeedIndex < readSeeds[r].size() && refSeedIndex < refSeeds.size()) {

            if (readSeeds[r][readSeedIndex].seq < refSeeds[refSeedIndex].seq) {
                
                readSeedIndex++;
            } else if (readSeeds[r][readSeedIndex].seq > refSeeds[refSeedIndex].seq) {
                
                refSeedIndex++;
            } else {

                matchingSeeds.push_back(readSeeds[r][readSeedIndex]);
                readSeedIndex++;
                refSeedIndex++;
            }
        }
        readSeeds[r] = matchingSeeds;
    }
    
    /*
    std::cerr << "MATCHED--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n";
    

    for(int i = 0; i < readSequences.size(); i++) {
        std::cerr << "\n" << i <<"\n";
        std::cerr << readSequences[i] << "\n";
        std::cerr << readQuals[i] << "\n";
        std::cerr << readSeeds[i].size() << "\n";
        for(int j = 0; j < readSeeds[i].size() ; j++){
            std::cerr << "XXXXXXXXXXXXXXXXXXXXXXXX\n";
            std::cerr << readSeeds[i][j].seq.size() << "\n";
            std::cerr << readSeeds[i][j].seq << "\n";
            std::cerr << readSeeds[i][j].pos << "\n";
            std::cerr << readSeeds[i][j].reversed << "\n";
        }
    }
    */
    



    //TODO SPLIT THIS FUNCTION APART
    
    //C++ to C interface will have to be cleaned up
    const char *reference = bestMatchSequence.c_str();
    int n_reads = readSequences.size();
    const char **read_strings = (const char **)malloc(n_reads*sizeof(char *));
    const char **qual_strings = (const char **)malloc(n_reads*sizeof(char *));
    const char **read_names = (const char **)malloc(n_reads*sizeof(char *));

    int *r_lens         = (int *)malloc(n_reads*sizeof(int));
    int *seed_counts    = (int *)malloc(n_reads*sizeof(int));

    uint8_t **reversed  = (uint8_t **)malloc(n_reads*sizeof(uint8_t *));
    int **ref_positions = (int **)malloc(n_reads*sizeof(int *));
    int **qry_positions = (int **)malloc(n_reads*sizeof(int *));
    

    for(int i = 0; i < n_reads; i++) {
        int n_seeds = readSeeds[i].size();

        seed_counts[i] = n_seeds;
        read_strings[i] = readSequences[i].c_str();
        qual_strings[i] = readQuals[i].c_str();
        read_names[i] = readNames[i].c_str();

        r_lens[i] = readSequences[i].length();

        uint8_t *reversed_array = (uint8_t *)malloc(n_seeds*sizeof(uint8_t));
        int *ref_pos_array = (int *)malloc(n_seeds*sizeof(int));
        int *qry_pos_array = (int *)malloc(n_seeds*sizeof(int));


        for(int j = 0; j < n_seeds; j++){
            reversed_array[j] = readSeeds[i][j].reversed;
            qry_pos_array[j] = readSeeds[i][j].pos;
            ref_pos_array[j] = seedToRefPositions[readSeeds[i][j].seq][0] + k - 1;
        }

        reversed[i]      = reversed_array;
        ref_positions[i] = ref_pos_array;
        qry_positions[i] = qry_pos_array;
    }
    
    std::string sam_header = "@SQ\tSN:reference\tLN:";
    sam_header += std::to_string(bestMatchSequence.length());

    std::cout << "\n" << sam_header << "\n";

    char *sam_alignments[n_reads]; //constituants must be freed


    align_reads(reference, n_reads, read_strings,qual_strings, read_names, r_lens, seed_counts, reversed, ref_positions, qry_positions, sam_alignments, k);
    
    
    /* Sorting the alignments by their reference position */
    std::pair<int, char*> sam_lines[n_reads];
    for(int i = 0; i < n_reads; i++) {
        sam_lines[i] = std::make_pair(r_lens[i],sam_alignments[i]);
    }
    sort(sam_lines, sam_lines + n_reads, [](const std::pair<int, char*>& a, const std::pair<int, char*>& b) {
        return a.first < b.first;
    });


    for(int i = 0; i < n_reads; i++) {
        sam_alignments[i] = sam_lines[i].second;
    }

    for(int i = 0; i < n_reads; i++) {
        if(sam_alignments[i]) {
            std::cout << sam_alignments[i] << std::endl;
        }
    }

    for(int i = 0; i < n_reads; i++) {
        free(reversed[i]);
        free(ref_positions[i]);
        free(qry_positions[i]);
    }
    free(qry_positions);
    free(ref_positions);
    free(reversed);
    free(seed_counts);
    free(read_strings);
    free(qual_strings);
    free(read_names);










    //Convert to BAM and beyond
    sam_hdr_t *header = NULL;
    bam1_t *record = bam_init1();

    
    // Parse SAM header
    header = sam_hdr_parse(sam_header.length(), sam_header.c_str());

    htsFile *bam_file = hts_open("output.bam", "wb");
    if (!bam_file) {
        fprintf(stderr, "Error: Failed to open output BAM file.\n");
        bam_hdr_destroy(header);
        bam_destroy1(record);
        return;
    }


    // Write BAM header
    if (sam_hdr_write(bam_file, header) < 0) {
        fprintf(stderr, "Error: Failed to write BAM header.\n");
        bam_hdr_destroy(header);
        bam_destroy1(record);
        hts_close(bam_file);
        return;
    }


    for (int i = 0; i < n_reads; i++) {
        if(sam_alignments[i]){

            kstring_t line = KS_INITIALIZE;
            kputs(sam_alignments[i], &line);
            
            sam_parse1(&line, header, record);

            //free line
            if (bam_write1(bam_file->fp.bgzf, record) < 0) {
                fprintf(stderr, "Error: Failed to write BAM record.\n");
                bam_hdr_destroy(header);
                bam_destroy1(record);
                hts_close(bam_file);
                return;
            }
        }
    }




    bam_hdr_destroy(header);
    bam_destroy1(record);
    hts_close(bam_file);
    /*
    */

    







    for(int i = 0; i < n_reads; i++) {
        free(sam_alignments[i]);
    }
    free(r_lens);
}