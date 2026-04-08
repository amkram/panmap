#include "conversion.hpp"
#include "genotyping.hpp"
#include "logging.hpp"
#include <algorithm>
#include <sys/wait.h>
extern "C" {
#include "mm_align.h"
#include "pileup.h"
#include <bcftools/bcftools.h>
}

// Run a bcftools function in a forked child process to avoid global state conflicts.
// bcftools main_mpileup/main_vcfcall use global optind and other process-wide state,
// so forking gives each call its own address space for safe parallelism.
static int run_bcftools_in_fork(int (*func)(int, char**), int argc, char** argv) {
    pid_t pid = fork();
    if (pid == 0) {
        // Child: run bcftools function and exit
        optind = 1;
        func(argc, argv);
        _exit(0);
    } else if (pid > 0) {
        int status;
        waitpid(pid, &status, 0);
        if (WIFEXITED(status) && WEXITSTATUS(status) == 0) return 0;
        return 1;
    } else {
        logging::err("fork() failed for bcftools call");
        return 1;
    }
}

// samAlignment is sorted at the end

void createSam(
  std::vector<std::string> &readSequences,
  std::vector<std::string> &readQuals,
  std::vector<std::string> &readNames,
  std::string &bestMatchSequence,
  std::string &samFileName,
  bool pairedEndReads,
  std::vector<char *> &samAlignments,
  std::string &samHeader
)   {
    //Preparing C structures for minimap
    const char *reference = bestMatchSequence.c_str();
    int n_reads = readSequences.size();
    const char **read_strings = (const char **)malloc(n_reads*sizeof(char *));
    const char **qual_strings = (const char **)malloc(n_reads*sizeof(char *));
    const char **read_names = (const char **)malloc(n_reads*sizeof(char *));
    int *r_lens = (int *)malloc(n_reads*sizeof(int));

    for(int i = 0; i < n_reads; i++) {
        read_strings[i] = readSequences[i].c_str();
        qual_strings[i] = readQuals[i].c_str();
        read_names[i] = readNames[i].c_str();
        r_lens[i] = readSequences[i].length();
    }

    samHeader = "@SQ\tSN:ref\tLN:";
    samHeader += std::to_string(bestMatchSequence.length());

    std::vector<char *> samAlignmentsBuffer(n_reads);
    char **sam_alignments = &samAlignmentsBuffer[0];

    align_reads(reference, n_reads, read_strings, qual_strings, read_names,
                r_lens, sam_alignments, pairedEndReads);

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
            logging::info("Wrote sam data to {}", samFileName);
        } else {
            logging::err("Failed to write to file {}", samFileName);
        }
    }

    /* Sorting the alignments by their reference position */
    std::vector<std::pair<int, char*>> sam_lines(n_reads);
    for(int i = 0; i < n_reads; i++) {
        sam_lines[i] = std::make_pair(r_lens[i], sam_alignments[i]);
    }
    sort(sam_lines.begin(), sam_lines.end(), [](const std::pair<int, char*>& a, const std::pair<int, char*>& b) {
        return a.first < b.first;
    });

    int numAlignedReads = n_reads;
    for(int i = 0; i < n_reads; i++) {
        sam_alignments[i] = sam_lines[i].second;
        if(!sam_alignments[i]){
            numAlignedReads = i;
            break;
        }
    }

    samAlignments.resize(numAlignedReads);
    for(int i = 0; i < numAlignedReads; i++){
        samAlignments[i] = sam_alignments[i];
    }

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
      output::error("Failed to open output BAM file: {}", bamFileName);
      return;
    }
    // Write BAM header
    if (sam_hdr_write(bam_file, header) < 0) {
      output::error("Failed to write BAM header");
      hts_close(bam_file);
      return;
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
        output::error("Failed to write BAM record");
      }
    }
  }

  if (bam_file) {
    logging::info("Wrote bam files to {}", bamFileName);
    hts_close(bam_file);
  }
  for (int i = 0; i < samAlignments.size(); i++) {
    free(samAlignments[i]);
  }

  return;
}

void createMplpBcf(
  std::string &prefix,
  std::string &refFileName,
  std::string &bestMatchSequence,
  std::string &bamFileName,
  std::string &mpileupFileName,
  bool baq
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

  {
    if (baq) {
      const char *mpileup_args[] = {"mpileup",
                                    "-Ou",
                                    "-f",
                                    outRefFileName.c_str(),
                                    "-o",
                                    mpileupFileName.c_str(),
                                    bamFileName.c_str()};
      run_bcftools_in_fork(main_mpileup, 7, const_cast<char **>(mpileup_args));
    } else {
      const char *mpileup_args[] = {"mpileup",
                                    "-Ou",
                                    "-B",
                                    "-f",
                                    outRefFileName.c_str(),
                                    "-o",
                                    mpileupFileName.c_str(),
                                    bamFileName.c_str()};
      run_bcftools_in_fork(main_mpileup, 8, const_cast<char **>(mpileup_args));
    }
  }
}

void createVcfWithMutationMatrices(
    std::string &prefix,
    std::string &mpileupFileName,
    std::string &vcfFileName,
    const std::vector<std::vector<double>> &substMatrixPhred
) {
  if (mpileupFileName.size() == 0) {
    mpileupFileName = prefix + ".mpileup";
  }

  // Use a temp file for raw VCF output (thread-safe via fork isolation)
  std::string rawVcfFile = prefix + ".vcf.raw";
  {
    const char *call_args[] = {
        "call", "--ploidy", "1", "-c", "-A", "-O", "v",
        "-o", rawVcfFile.c_str(), mpileupFileName.c_str()};
    run_bcftools_in_fork(main_vcfcall, 10, const_cast<char **>(call_args));
  }

  // Read raw VCF and apply mutation spectrum or filter non-ref variants
  std::ifstream rawVcfIn(rawVcfFile);
  bool hasSpectrum = !substMatrixPhred.empty();
  std::ofstream vcfOutFile{vcfFileName};
  std::string line;
  while (std::getline(rawVcfIn, line)) {
    if (hasSpectrum) {
      std::string spectrum_applied_line =
          genotyping::applyMutationSpectrum(line, substMatrixPhred);
      if (spectrum_applied_line.size() > 0)
        vcfOutFile << spectrum_applied_line << "\n";
    } else {
      // No spectrum: write non-ref variants only
      if (line.empty()) continue;
      if (line[0] == '#') {
        vcfOutFile << line << "\n";
        continue;
      }
      std::istringstream iss(line);
      std::string f;
      std::vector<std::string> fields;
      while (std::getline(iss, f, '\t')) fields.push_back(f);
      if (fields.size() >= 10 && fields[4] != "." && fields[9][0] != '0') {
        vcfOutFile << line << "\n";
      }
    }
  }
  rawVcfIn.close();
  std::remove(rawVcfFile.c_str());
}

// ============================================================================
// Direct alignment-to-BAM (opts 2+3: parallel align, skip SAM text)
// ============================================================================

static uint16_t compute_sam_flags(bool is_paired, bool is_read1,
                                   uint8_t rev, uint8_t mate_rev,
                                   uint8_t proper_frag, bool mate_unmapped) {
    uint16_t flag = 0;
    if (is_paired) {
        flag |= BAM_FPAIRED;
        if (proper_frag) flag |= BAM_FPROPER_PAIR;
        if (rev) flag |= BAM_FREVERSE;
        if (mate_rev) flag |= BAM_FMREVERSE;
        if (mate_unmapped) flag |= BAM_FMUNMAP;
        if (is_read1) flag |= BAM_FREAD1;
        else flag |= BAM_FREAD2;
    } else {
        if (rev) flag |= BAM_FREVERSE;
    }
    return flag;
}

static int32_t compute_tlen(int32_t this_rs, int32_t this_re, uint8_t this_rev,
                             int32_t mate_rs, int32_t mate_re, uint8_t mate_rev) {
    int this_pos5 = this_rev ? this_re - 1 : this_rs;
    int mate_pos5 = mate_rev ? mate_re - 1 : mate_rs;
    int tlen = mate_pos5 - this_pos5;
    if (tlen > 0) tlen++;
    else if (tlen < 0) tlen--;
    return tlen;
}

static bam1_t *build_bam_from_result(
    const std::string &qname_full, const std::string &seq, const std::string &qual,
    const read_align_t *aln, int read_len,
    bool is_paired, bool is_read1, uint8_t mate_rev, int32_t mate_pos,
    int32_t this_rs, int32_t this_re, int32_t mate_rs, int32_t mate_re,
    uint8_t proper_frag, bool mate_unmapped,
    sam_hdr_t *header)
{
    bam1_t *b = bam_init1();

    // Strip /1 or /2 suffix from read name (matching minimap2 convention)
    std::string qname = qname_full;
    if (qname.size() >= 2 && qname[qname.size()-2] == '/' &&
        (qname.back() == '1' || qname.back() == '2')) {
        qname.resize(qname.size() - 2);
    }

    uint16_t flag = compute_sam_flags(is_paired, is_read1, aln->rev, mate_rev,
                                       proper_frag, mate_unmapped);

    // Build full CIGAR with soft clips for unaligned query ends
    uint32_t clip5 = aln->rev ? (read_len - aln->qe) : aln->qs;
    uint32_t clip3 = aln->rev ? aln->qs : (read_len - aln->qe);

    int full_n_cigar = aln->n_cigar + (clip5 > 0 ? 1 : 0) + (clip3 > 0 ? 1 : 0);
    std::vector<uint32_t> full_cigar(full_n_cigar);
    int ci = 0;
    if (clip5 > 0) full_cigar[ci++] = (clip5 << 4) | BAM_CSOFT_CLIP;
    for (int j = 0; j < aln->n_cigar; j++) full_cigar[ci++] = aln->cigar[j];
    if (clip3 > 0) full_cigar[ci++] = (clip3 << 4) | BAM_CSOFT_CLIP;

    // Prepare sequence and quality in BAM orientation
    // BAM stores reverse-strand reads as reverse-complemented seq with reversed qual
    std::string bam_seq(read_len, 'N');
    std::string bam_qual(read_len, '\0');
    if (aln->rev) {
        for (int i = 0; i < read_len; i++) {
            char c = seq[read_len - 1 - i];
            switch (c) {
                case 'A': case 'a': bam_seq[i] = 'T'; break;
                case 'T': case 't': bam_seq[i] = 'A'; break;
                case 'C': case 'c': bam_seq[i] = 'G'; break;
                case 'G': case 'g': bam_seq[i] = 'C'; break;
                default: bam_seq[i] = 'N'; break;
            }
        }
        for (int i = 0; i < read_len; i++) {
            bam_qual[i] = qual[read_len - 1 - i] - 33;
        }
    } else {
        bam_seq = seq;
        for (int i = 0; i < read_len; i++) {
            bam_qual[i] = qual[i] - 33;
        }
    }

    // Compute TLEN and mate info
    int32_t tlen = 0;
    int32_t mtid = -1;
    hts_pos_t mpos_val = -1;
    if (is_paired) {
        tlen = compute_tlen(this_rs, this_re, aln->rev, mate_rs, mate_re, mate_rev);
        mtid = 0;
        mpos_val = mate_pos;
    }

    bam_set1(b,
             qname.length(), qname.c_str(),
             flag,
             0,                  // tid = 0
             (hts_pos_t)aln->rs, // 0-based position
             aln->mapq,
             full_n_cigar, full_cigar.data(),
             mtid, mpos_val, (hts_pos_t)tlen,
             read_len, bam_seq.c_str(), bam_qual.c_str(),
             0);

    return b;
}

void alignAndWriteBam(
    std::vector<std::string> &readSequences,
    std::vector<std::string> &readQuals,
    std::vector<std::string> &readNames,
    std::string &reference,
    const std::string &bamFileName,
    bool pairedEndReads,
    int n_threads)
{
    int n_reads = (int)readSequences.size();

    // Prepare C arrays
    std::vector<const char *> read_ptrs(n_reads);
    std::vector<const char *> qual_ptrs(n_reads);
    std::vector<const char *> name_ptrs(n_reads);
    std::vector<int> r_lens(n_reads);
    for (int i = 0; i < n_reads; i++) {
        read_ptrs[i] = readSequences[i].c_str();
        qual_ptrs[i] = readQuals[i].c_str();
        name_ptrs[i] = readNames[i].c_str();
        r_lens[i] = (int)readSequences[i].length();
    }

    int n_results = pairedEndReads ? n_reads / 2 : n_reads;
    std::vector<align_pair_result_t> results(n_results);

    align_reads_direct(reference.c_str(), n_reads,
                       read_ptrs.data(), qual_ptrs.data(),
                       name_ptrs.data(), r_lens.data(),
                       results.data(), pairedEndReads, n_threads);

    // Build SAM header
    std::string samHeaderStr = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:ref\tLN:"
                               + std::to_string(reference.length());
    sam_hdr_t *header = sam_hdr_parse(samHeaderStr.length(), samHeaderStr.c_str());

    // Build bam1_t records with (sort_pos, bam1_t*) for position sorting
    std::vector<std::pair<int32_t, bam1_t*>> bam_entries;
    bam_entries.reserve(pairedEndReads ? n_results * 2 : n_results);

    for (int k = 0; k < n_results; k++) {
        align_pair_result_t *res = &results[k];
        if (!res->mapped) continue;

        if (pairedEndReads) {
            int r1_idx = k * 2;
            int r2_idx = k * 2 + 1;

            bam1_t *b1 = build_bam_from_result(
                readNames[r1_idx], readSequences[r1_idx], readQuals[r1_idx],
                &res->r1, r_lens[r1_idx],
                true, true, res->r2.rev, res->r2.rs,
                res->r1.rs, res->r1.re, res->r2.rs, res->r2.re,
                res->r1.proper_frag, false, header);
            bam_entries.push_back({res->r1.pos, b1});

            bam1_t *b2 = build_bam_from_result(
                readNames[r2_idx], readSequences[r2_idx], readQuals[r2_idx],
                &res->r2, r_lens[r2_idx],
                true, false, res->r1.rev, res->r1.rs,
                res->r2.rs, res->r2.re, res->r1.rs, res->r1.re,
                res->r2.proper_frag, false, header);
            bam_entries.push_back({res->r2.pos, b2});
        } else {
            bam1_t *b = build_bam_from_result(
                readNames[k], readSequences[k], readQuals[k],
                &res->r1, r_lens[k],
                false, false, 0, -1,
                0, 0, 0, 0, 0, false, header);
            bam_entries.push_back({res->r1.pos, b});
        }
    }

    // Sort by reference position
    std::sort(bam_entries.begin(), bam_entries.end(),
              [](const auto &a, const auto &b) { return a.first < b.first; });

    // Write BAM file
    if (!bamFileName.empty()) {
        htsFile *bam_file = hts_open(bamFileName.c_str(), "wb");
        if (bam_file) {
            sam_hdr_write(bam_file, header);
            for (auto &[pos, rec] : bam_entries) {
                bam_write1(bam_file->fp.bgzf, rec);
            }
            hts_close(bam_file);
            logging::info("Wrote bam files to {}", bamFileName);
        } else {
            logging::err("Failed to open output BAM file: {}", bamFileName);
        }
    }

    // Cleanup
    for (auto &[pos, rec] : bam_entries) {
        bam_destroy1(rec);
    }
    sam_hdr_destroy(header);
    for (int k = 0; k < n_results; k++) {
        free(results[k].r1.cigar);
        free(results[k].r2.cigar);
    }
}