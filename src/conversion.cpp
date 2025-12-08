#include "conversion.hpp"
#include "alignment.hpp"
#include "genotyping.hpp"
#include <iostream>
extern "C" {
#include "mm_align.h"
#include "pileup.h"
#include <bcftools/bcftools.h>
}

// samAlignment is sorted at the end

void addSeeds(std::vector<seeding::seed_t> &fwdmatchingSeeds,
              std::vector<seeding::seed_t> &bwdmatchingSeeds,
              const std::vector<uint32_t> &positions,
              const seeding::seed_t &curSeed,
              bool reverseCondition,
              int k,
              bool shortenSyncmers,
              int readlen) {
  for (uint32_t rpos : positions) {
    seeding::seed_t thisSeed = curSeed;
    thisSeed.reversed = thisSeed.reversed == reverseCondition;
    // endPos is the last position of the seed in the query (0-indexed)
    // This is always pos + k - 1, regardless of strand. The reversed flag
    // tells minimap2 which strand the seed is on.
    thisSeed.endPos = thisSeed.pos + k - 1;
    thisSeed.rpos = rpos;

    if (thisSeed.reversed) {
      bwdmatchingSeeds.push_back(thisSeed);
    } else {
      fwdmatchingSeeds.push_back(thisSeed);
    }
  }
}


void createSam(
  std::vector<std::vector<seeding::seed_t>> &readSeeds,
  std::vector<std::string> &readSequences,
  std::vector<std::string> &readQuals,
  std::vector<std::string> &readNames,
  std::string &bestMatchSequence,
  std::unordered_map<size_t, std::pair<std::vector<uint32_t>, std::vector<uint32_t>>> &seedToRefPositions,
  std::string &samFileName,
  int k,
  bool shortenSyncmers,
  bool pairedEndReads,
  std::vector<char *> &samAlignments,
  std::string &samHeader
)   {
  for (size_t i = 0; i < readSequences.size(); ++i) {
    std::vector<seeding::seed_t>& curReadSeeds = readSeeds[i];
    std::vector<seeding::seed_t> fwdMatchingSeeds;
    std::vector<seeding::seed_t> bwdMatchingSeeds;
    
    for (size_t j = 0; j < curReadSeeds.size(); ++j) {
      if (seedToRefPositions.find(curReadSeeds[j].hash) == seedToRefPositions.end()) continue;
      const auto& [forwardSeedToRefPositions, reverseSeedToRefPositions] = seedToRefPositions[curReadSeeds[j].hash];
      addSeeds(fwdMatchingSeeds, bwdMatchingSeeds, forwardSeedToRefPositions, curReadSeeds[j], true, k, shortenSyncmers, readSequences[i].size());
      addSeeds(fwdMatchingSeeds, bwdMatchingSeeds, reverseSeedToRefPositions, curReadSeeds[j], false, k, shortenSyncmers, readSequences[i].size());
    }

    std::reverse(bwdMatchingSeeds.begin(), bwdMatchingSeeds.end());
    fwdMatchingSeeds.insert(fwdMatchingSeeds.end(), bwdMatchingSeeds.begin(), bwdMatchingSeeds.end());

    readSeeds[i] = fwdMatchingSeeds;
  }

    
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
                qry_pos_array[j] = readSeeds[i*2][j].endPos;
                ref_pos_array[j] = readSeeds[i*2][j].rpos + k - 1;
            }
            for(int j = 0; j < seed_counts[i*2 + 1]; j++){
                reversed_array[j + seed_counts[i*2]] = readSeeds[i*2 + 1][j].reversed;
                qry_pos_array[j + seed_counts[i*2]] = readSeeds[i*2 + 1][j].endPos;
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
                qry_pos_array[j] = readSeeds[i][j].endPos;
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

// samAlignments elements are freed
//
// bam_records must be freed and elements must be freed.
void createBam(
  std::vector<char *> &samAlignments,
  std::string &samHeader,
  std::string &bamFileName,
  sam_hdr_t *&header,
  bam1_t **&bamRecords) {

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

  // Prepare list of bam1_t
  bamRecords = (bam1_t **)malloc(sizeof(bam1_t *) * samAlignments.size());

  for (int i = 0; i < samAlignments.size(); i++) {
    if (samAlignments[i]) {

      bamRecords[i] = bam_init1();

      kstring_t line = KS_INITIALIZE;
      kputs(samAlignments[i], &line);

      sam_parse1(&line, header, bamRecords[i]);

      // Write to bam file
      if (bam_file && bam_write1(bam_file->fp.bgzf, bamRecords[i]) < 0) {
        fprintf(stderr, "Error: Failed to write BAM record.\n");
        bam_hdr_destroy(header);
        hts_close(bam_file);
      }
    }
  }

  if (bam_file) {
    std::cerr << "Wrote bam files to " << bamFileName << "\n";
  }
  /// Converted to Bam
  hts_close(bam_file);
  for (int i = 0; i < samAlignments.size(); i++) {
    free(samAlignments[i]);
  }

  return;
}

// Destroys header and bamRecords
//
// mplpString must be freed
void createMplp(
  std::string &bestMatchSequence,
  sam_hdr_t *header,
  bam1_t **bamRecords,
  int numBams,
  std::string &mpileupFileName,
  char *&mplpString) {

  char *ref_string = new char[bestMatchSequence.length() + 1];
  std::strcpy(ref_string, bestMatchSequence.c_str());

  kstring_t mplp_string = KS_INITIALIZE;
  bam_and_ref_to_mplp(header, bamRecords, numBams, ref_string,
                      bestMatchSequence.size(), &mplp_string);

  // Print out mpileup
  if (mpileupFileName.size() > 0) {
    std::ofstream outFile{mpileupFileName};

    if (outFile.is_open()) {

      outFile << mplp_string.s;

      std::cerr << "Wrote mpileup data to " << mpileupFileName << std::endl;
    } else {
      std::cerr << "Error: failed to write to file " << mpileupFileName
                << std::endl;
    }
  }

  mplpString = mplp_string.s;

  for (int i = 0; i < numBams; i++) {
    bam_destroy1(bamRecords[i]);
  }

  free(bamRecords);
}

void createMplpBcf(
  std::string &prefix,
  std::string &refFileName,
  std::string &bestMatchSequence,
  std::string &bamFileName,
  std::string &mpileupFileName
) {
  std::string outRefFileName = "";
  if (refFileName.size() == 0) {
    outRefFileName = prefix + ".tmp.reference.fa";
    std::ofstream outRefFile{outRefFileName};
    outRefFile << ">ref\n";
    outRefFile << bestMatchSequence << "\n";
    outRefFile.close();
  } else {
    outRefFileName = refFileName;
  }

  if (mpileupFileName.size() == 0) {
    mpileupFileName = prefix + ".mpileup";
  }

  optind = 1;
  const char *mpileup_args[] = {"mpileup",
                                "-Ou",
                                "-f",
                                outRefFileName.c_str(),
                                bamFileName.c_str(),
                                "-o",
                                mpileupFileName.c_str()};
  main_mpileup(7, const_cast<char **>(mpileup_args));
}

void createVcfWithMutationMatrices(
    std::string &prefix,
    std::string &mpileupFileName,
    const genotyping::mutationMatrices &mutMat,
    std::string &vcfFileName,
    double mutationRate
) {
  if (mpileupFileName.size() == 0) {
    mpileupFileName = prefix + ".mpileup";
  }

  FILE *tempFile = std::tmpfile();
  int stdoutFd = dup(fileno(stdout));
  dup2(fileno(tempFile), fileno(stdout));
  optind = 1;
  const char *mpileup_args[] = {
      "call", "--ploidy", "1", "-c", "-A", "-O", "v", mpileupFileName.c_str()};
  main_vcfcall(8, const_cast<char **>(mpileup_args));

  fflush(stdout);
  dup2(stdoutFd, fileno(stdout));
  close(stdoutFd);
  rewind(tempFile);

  std::vector<std::vector<double>> scaled_submat =
      genotyping::scaleMutationSpectrum(mutMat, mutationRate);
  std::ofstream vcfOutFile{vcfFileName};
  std::vector<std::string> vcfLines;
  char buffer[512];
  while (fgets(buffer, sizeof(buffer), tempFile)) {
    std::string line(buffer);
    if (line.size() > 0 && line[line.size() - 1] == '\n')
      line.pop_back();
    std::string spectrum_applied_line =
        genotyping::applyMutationSpectrum(line, scaled_submat);
    if (spectrum_applied_line.size() > 0)
      vcfOutFile << spectrum_applied_line << "\n";
  }
  fclose(tempFile);
}

// destroys mplpString
void createVcf(char *mplpString,
               const genotyping::mutationMatrices &mutMat,
               std::string &vcfFileName, bool keep_alts) {
  // Convert c string of mpileup to ifstream
  std::istringstream mpileipStream(mplpString);

  std::ofstream vcfOutFile;
  if (vcfFileName.size() > 0) {
    vcfOutFile.open(vcfFileName);
    if (vcfOutFile.is_open()) {

      genotyping::printSamplePlacementVCF(mpileipStream, mutMat, keep_alts, 0,
                                          vcfOutFile);

      std::cerr << "Wrote vcf data to " << vcfFileName << std::endl;
    } else {

      std::cerr << "Error: failed to write to file " << vcfFileName
                << std::endl;
    }
  }

  free(mplpString);
}