#include <iostream>
#include "conversion.hpp"
#include "genotype.hpp"

extern "C" {
    #include "pileup.h"
    #include "mm_align.h"
}

//samAlignment is sorted at the end

void addSeeds(std::vector<seeding::seed>& fwdmatchingSeeds,std::vector<seeding::seed>& bwdmatchingSeeds, const std::vector<uint32_t>& positions, 
              const seeding::seed& curSeed, bool reverseCondition, int k, int readlen) {
    for (uint32_t rpos : positions) {
        seed thisSeed;
        thisSeed.reversed = curSeed.reversed == reverseCondition;
        thisSeed.hash = curSeed.hash;
        thisSeed.pos = (thisSeed.reversed) ? readlen - curSeed.pos + k - 2 :curSeed.pos;
        
        thisSeed.rpos = rpos;

        if(thisSeed.reversed){
            bwdmatchingSeeds.push_back(thisSeed);
        }else{
            fwdmatchingSeeds.push_back(thisSeed);
        }
        
    }
}

//Elements of samAlignments must be freed
void createSam(
    std::vector<std::vector<seeding::seed>> &readSeeds,
    std::vector<std::string> &readSequences,
    std::vector<std::string> &readQuals,
    std::vector<std::string> &readNames,
    std::string &bestMatchSequence,
    std::unordered_map<size_t, std::pair<std::vector<uint32_t>, std::vector<uint32_t>>> &seedToRefPositions,
    std::string &samFileName,
    int k,
    bool pairedEndReads,
    
    std::vector<char *> &samAlignments,
    std::string &samHeader
)   {

  /*
  for(int i = 0; i < readSeeds.size(); i++){
    std::cerr << "Read " <<i<< " " << readSequences[i] << "\n";
    for(int j = 0; j < readSeeds[i].size(); j++){
        std::cerr << readSeeds[i][j].hash << " " << readSeeds[i][j].reversed << " " << readSeeds[i][j].pos << "\n";
    }
  }*/

  
  for (size_t i = 0; i < readSequences.size(); ++i) {
    std::vector<seeding::seed>& curReadSeeds = readSeeds[i];
    std::vector<seeding::seed> fwdMatchingSeeds;
    std::vector<seeding::seed> bwdMatchingSeeds;
    

    for (size_t j = 0; j < curReadSeeds.size(); ++j) {
      if (seedToRefPositions.find(curReadSeeds[j].hash) == seedToRefPositions.end()) continue;
      const auto& [forwardSeedToRefPositions, reverseSeedToRefPositions] = seedToRefPositions[curReadSeeds[j].hash];

      
      addSeeds(fwdMatchingSeeds, bwdMatchingSeeds, forwardSeedToRefPositions, curReadSeeds[j], true, k , readSequences[i].size());
      addSeeds(fwdMatchingSeeds, bwdMatchingSeeds, reverseSeedToRefPositions, curReadSeeds[j], false, k, readSequences[i].size());
    }

    std::reverse(bwdMatchingSeeds.begin(), bwdMatchingSeeds.end());
    fwdMatchingSeeds.insert(fwdMatchingSeeds.end(), bwdMatchingSeeds.begin(), bwdMatchingSeeds.end());

    readSeeds[i] = fwdMatchingSeeds;
    
  }

    
  

    
    /*
    for(int r = 0; r < readSequences.size() ; r++){

        std::vector<seed> matchingSeeds;
        for(int i = 0; i < readSeeds[r].size(); i++){

            for (int32_t rpos : seedToRefPositions[readSeeds[r][i].hash]) {
                //if(readSeeds[r][i].reversed == ){ //TODO remove this
                seed thisSeed;
                thisSeed.reversed = readSeeds[r][i].reversed;
                thisSeed.hash = readSeeds[r][i].hash;
                thisSeed.pos = readSeeds[r][i].pos;
                thisSeed.rpos = rpos;
                matchingSeeds.push_back(thisSeed);
                //}
            }

        }

        readSeeds[r] = matchingSeeds;
    }*/

   /*
    for(int r = 0; r < readSequences.size() ; r++){
        std::cout << "READ " << r << ": " << readSequences[r] << "\n";
        for(int i = 0; i < readSeeds[r].size(); i++){
            std::cout << "seed: " << readSeeds[r][i].hash << " " << readSeeds[r][i].reversed 
                << " " <<  readSequences[r].substr(readSeeds[r][i].pos, k) << " " << bestMatchSequence.substr(readSeeds[r][i].rpos, k) << "\n";
            if(readSequences[r].substr(readSeeds[r][i].pos, k) !=  bestMatchSequence.substr(readSeeds[r][i].rpos, k)){
                std::cout << "\n\n\nSCREAM \n\n\n\n";
            }
        }
    }*/

    
    //Preparing C structures for minimap
    const char *reference = bestMatchSequence.c_str();
    int n_reads = readSequences.size();
    const char **read_strings = (const char **)malloc(n_reads*sizeof(char *));
    const char **qual_strings = (const char **)malloc(n_reads*sizeof(char *));
    const char **read_names = (const char **)malloc(n_reads*sizeof(char *));

    int *r_lens         = (int *)malloc(n_reads*sizeof(int));
    int *seed_counts    = (int *)malloc(n_reads*sizeof(int));

    for(int i = 0; i < n_reads; i++) {
        int n_seeds = readSeeds[i].size();

        seed_counts[i] = n_seeds;
        read_strings[i] = readSequences[i].c_str();
        qual_strings[i] = readQuals[i].c_str();
        read_names[i] = readNames[i].c_str();
        r_lens[i] = readSequences[i].length();
        
    }

    uint8_t **reversed;
    int **ref_positions;
    int **qry_positions;

    if(pairedEndReads){
        reversed  = (uint8_t **)malloc((n_reads/2)*sizeof(uint8_t *));
        ref_positions = (int **)malloc((n_reads/2)*sizeof(int *));
        qry_positions = (int **)malloc((n_reads/2)*sizeof(int *));

        for(int i = 0; i < n_reads/2; i++) {
            
            int n_seeds = seed_counts[i*2] + seed_counts[i*2+1];

            uint8_t *reversed_array = (uint8_t *)malloc(n_seeds*sizeof(uint8_t));
            int *ref_pos_array = (int *)malloc(n_seeds*sizeof(int));
            int *qry_pos_array = (int *)malloc(n_seeds*sizeof(int));

            for(int j = 0; j < seed_counts[i*2]; j++){
                reversed_array[j] = readSeeds[i*2][j].reversed;
                qry_pos_array[j] = readSeeds[i*2][j].pos;
                ref_pos_array[j] = readSeeds[i*2][j].rpos + k - 1;
            }
            for(int j = 0; j < seed_counts[i*2 + 1]; j++){
                reversed_array[j + seed_counts[i*2]] = readSeeds[i*2 + 1][j].reversed;
                qry_pos_array[j + seed_counts[i*2]] = readSeeds[i*2 + 1][j].pos;
                ref_pos_array[j + seed_counts[i*2]] = readSeeds[i*2 + 1][j].rpos + k - 1;
            }

            reversed[i]      = reversed_array;
            ref_positions[i] = ref_pos_array;
            qry_positions[i] = qry_pos_array;
        }

    }else{

        reversed  = (uint8_t **)malloc(n_reads*sizeof(uint8_t *));
        ref_positions = (int **)malloc(n_reads*sizeof(int *));
        qry_positions = (int **)malloc(n_reads*sizeof(int *));

        for(int i = 0; i < n_reads; i++) {

            int n_seeds = seed_counts[i];

            uint8_t *reversed_array = (uint8_t *)malloc(n_seeds*sizeof(uint8_t));
            int *ref_pos_array = (int *)malloc(n_seeds*sizeof(int));
            int *qry_pos_array = (int *)malloc(n_seeds*sizeof(int));

            for(int j = 0; j < n_seeds; j++){
                reversed_array[j] = readSeeds[i][j].reversed;
                qry_pos_array[j] = readSeeds[i][j].pos;
                ref_pos_array[j] = readSeeds[i][j].rpos + k - 1;
            }

            reversed[i]      = reversed_array;
            ref_positions[i] = ref_pos_array;
            qry_positions[i] = qry_pos_array;
        }
    }

    
    samHeader = "@SQ\tSN:ref\tLN:";
    samHeader += std::to_string(bestMatchSequence.length());

    std::vector<char *> samAlignmentsBuffer(n_reads);

    char **sam_alignments = &samAlignmentsBuffer[0];



    align_reads(reference,n_reads,read_strings,qual_strings, read_names, r_lens, seed_counts, reversed, ref_positions, qry_positions, sam_alignments, k, pairedEndReads);


    
    //Print out sam
    if(samFileName.size() > 0){
        std::ofstream outFile{samFileName};

        if (outFile.is_open()) {

            outFile << samHeader << std::endl;
            
            for(int i = 0; i < n_reads; i++) {
                if(sam_alignments[i]) {
                    outFile << sam_alignments[i] << std::endl;
                }
            }

            std::cerr << "Wrote sam data to " << samFileName << std::endl;
        } else {
            std::cerr << "Error: failed to write to file " << samFileName << std::endl;
        }
    }
    

    
    /* Sorting the alignments by their reference position */
    std::vector<std::pair<int, char*>> sam_lines(n_reads);

    for(int i = 0; i < n_reads; i++) {
        sam_lines[i] = std::make_pair(r_lens[i],sam_alignments[i]);
    }
    sort(sam_lines.begin(), sam_lines.end(), [](const std::pair<int, char*>& a, const std::pair<int, char*>& b) {
        return a.first < b.first;
    });

    
    
    int numAlignedReads = n_reads;
    
    for(int i = 0; i < n_reads; i++) {
        sam_alignments[i] = sam_lines[i].second;
        if(!sam_alignments[i]){
            numAlignedReads = i;         //Some reads failed
            break;
        }
    }

    
    samAlignments.resize(numAlignedReads);
    for(int i = 0; i < numAlignedReads; i++){
        samAlignments[i] = sam_alignments[i];
    }


    if( pairedEndReads )
        n_reads /= 2;

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
    free(r_lens);

}

//samAlignments elements are freed
//
//bam_records must be freed and elements must be freed. 
void createBam(
    std::vector<char *> &samAlignments,
    std::string &samHeader,
    std::string &bamFileName,

    sam_hdr_t * &header,
    bam1_t ** &bamRecords
)   {

    // Parse SAM header
    header = sam_hdr_parse(samHeader.length(), samHeader.c_str());
    htsFile *bam_file = NULL;

    if (bamFileName.size() > 0) {
        bam_file = hts_open(bamFileName.c_str(), "wb");
        if (!bam_file) {
            fprintf(stderr, "Error: Failed to open output BAM file.\n");
            hts_close(bam_file);
        }
        // Write BAM header
        else if (sam_hdr_write(bam_file, header) < 0) {
            fprintf(stderr, "Error: Failed to write BAM header.\n");
        }
    }
    
    //Prepare list of bam1_t
    bamRecords = (bam1_t **)malloc(sizeof(bam1_t *) * samAlignments.size());

    for (int i = 0; i < samAlignments.size(); i++) {
        if(samAlignments[i]){
            
            bamRecords[i] = bam_init1();

            kstring_t line = KS_INITIALIZE;
            kputs(samAlignments[i], &line);
            
            sam_parse1(&line, header, bamRecords[i]);

            //Write to bam file
            if (bam_file && bam_write1(bam_file->fp.bgzf, bamRecords[i]) < 0) {
                fprintf(stderr, "Error: Failed to write BAM record.\n");
                bam_hdr_destroy(header);
                hts_close(bam_file);
            }
        }
    }
    
    if(bam_file){
        std::cerr << "Wrote bam files to " << bamFileName << "\n";
    }
    /// Converted to Bam
    hts_close(bam_file);
    for(int i = 0; i < samAlignments.size(); i++) {
        free(samAlignments[i]);
    }

    return;
}





//Destroys header and bamRecords
// 
//mplpString must be freed
void createMplp(
    std::string &bestMatchSequence,
    sam_hdr_t *header,
    bam1_t **bamRecords,
    int numBams,
    std::string &mpileupFileName,

    char * &mplpString
){

    char* ref_string = new char[bestMatchSequence.length() + 1];
    std::strcpy(ref_string, bestMatchSequence.c_str());

    kstring_t mplp_string = KS_INITIALIZE;
    bam_and_ref_to_mplp(header, bamRecords, numBams, ref_string, bestMatchSequence.size(), &mplp_string);

    //Print out mpileup
    if(mpileupFileName.size() > 0){
        std::ofstream outFile{mpileupFileName};

        if (outFile.is_open()) {
            
            outFile << mplp_string.s;

            std::cerr << "Wrote mpileup data to " << mpileupFileName << std::endl;
        } else {
            std::cerr << "Error: failed to write to file " << mpileupFileName << std::endl;
        }
    }

    mplpString = mplp_string.s;

    for(int i = 0; i < numBams; i++){
        bam_destroy1(bamRecords[i]);
    }
    
    free(bamRecords);
}


//destroys mplpString
void createVcf(
    char *mplpString,
    const seed_annotated_tree::mutationMatrices& mutMat,
    std::string &vcfFileName,
    bool keep_alts
)   {
    // Convert c string of mpileup to ifstream
    std::istringstream mpileipStream(mplpString);

    std::ofstream vcfOutFile;
    if(vcfFileName.size() > 0) {
        vcfOutFile.open(vcfFileName);
        if (vcfOutFile.is_open()) {

            genotype::printSamplePlacementVCF(mpileipStream, mutMat, keep_alts, 0, vcfOutFile);

            std::cerr << "Wrote vcf data to " << vcfFileName << std::endl;
        }else{

            std::cerr << "Error: failed to write to file " << vcfFileName << std::endl;
        }
    }

    free(mplpString);
}