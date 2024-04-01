#include <algorithm>
#include "place.hpp"
#include "pmi.hpp"
#include "util.hpp"
#include "tree.hpp"
#include "genotype.hpp"
#include <cmath>
#include <htslib/sam.h>

#include "conversion.hpp"

extern "C" {
    void align_reads(const char *reference, int n_reads, const char **reads, const char **quality, const char **read_names, int *r_lens, int *seed_counts, uint8_t **reversed, int **ref_positions, int **qry_positions, char** sam_alignments, int syncmer_k);
    void bam_and_ref_to_mplp(sam_hdr_t *header, bam1_t **bam_lines, int nbams, char *ref_string, int lref, kstring_t *mplp_string);
}


//We assume refSeeds and readSeeds are sorted.

//samAlignment is sorted at the end

//Elements of samAlignments must be freed
void createSam(
    std::vector<seeding::seed> &refSeeds,
    std::vector<std::vector<seeding::seed>> &readSeeds,
    std::vector<std::string> &readSequences,
    std::vector<std::string> &readQuals,
    std::vector<std::string> &readNames,
    std::string &bestMatchSequence,
    std::unordered_map<std::string, std::vector<int32_t>> &seedToRefPositions,
    std::string &samFileName,
    int k,
    
    std::vector<char *> &samAlignments,
    std::string &samHeader
)   {

    /* Finding syncmer matches between reads and ref */
    // replaces readSeeds elements with only matching syncmers
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
    
    //Preparing C structures for minimap
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
    
    samHeader = "@SQ\tSN:ref\tLN:";
    samHeader += std::to_string(bestMatchSequence.length());
    
    char *sam_alignments[n_reads];

    
    align_reads(reference, n_reads, read_strings,qual_strings, read_names, r_lens, seed_counts, reversed, ref_positions, qry_positions, sam_alignments, k);
    
    
    /* Sorting the alignments by their reference position */
    std::pair<int, char*> sam_lines[n_reads];
    for(int i = 0; i < n_reads; i++) {
        sam_lines[i] = std::make_pair(r_lens[i],sam_alignments[i]);
    }
    sort(sam_lines, sam_lines + n_reads, [](const std::pair<int, char*>& a, const std::pair<int, char*>& b) {
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

            std::cout << "Wrote sam data to " << samFileName << std::endl;
        } else {
            std::cerr << "Error: failed to write to file " << samFileName << std::endl;
        }
    }

    samAlignments.resize(numAlignedReads);
    for(int i = 0; i < numAlignedReads; i++){
        samAlignments[i] = sam_alignments[i];
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

            std::cout << "Wrote mpileup data to " << mpileupFileName << std::endl;
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
    const tree::mutationMatrices& mutMat,
    std::string &vcfFileName
)   {
    // Convert c string of mpileup to ifstream
    std::istringstream mpileipStream(mplpString);

    std::ofstream vcfOutFile;
    if(vcfFileName.size() > 0) {
        vcfOutFile.open(vcfFileName);
        if (vcfOutFile.is_open()) {

            genotype::printSamplePlacementVCF(mpileipStream, mutMat, true, 0, vcfOutFile);

            std::cout << "Wrote vcf data to " << vcfFileName << std::endl;
        }else{

            std::cerr << "Error: failed to write to file " << vcfFileName << std::endl;
        }
    }

    free(mplpString);
}