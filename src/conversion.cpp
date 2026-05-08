#include "conversion.hpp"
#include "genotyping.hpp"
#include "logging.hpp"
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <unistd.h>
#include <sys/wait.h>
extern "C" {
#include "mm_align.h"
#include "pileup.h"
#include <bcftools/bcftools.h>
#include <htslib/bgzf.h>
#include <htslib/sam.h>
#include <htslib/tbx.h>
}

// Run a bcftools function in a forked child process to avoid global state conflicts.
// bcftools main_mpileup/main_vcfcall use global optind and other process-wide state,
// so forking gives each call its own address space for safe parallelism.
static int run_bcftools_in_fork(int (*func)(int, char**), int argc, char** argv, bool silenceStderr = true) {
    bool verbose = output::config().verbose;
    bool fullSilence = silenceStderr && !verbose;
    bool filterNote = !silenceStderr && !verbose;  // keep meaningful chatter, drop "Note: ..." lines

    int pipeFds[2] = {-1, -1};
    if (filterNote && pipe(pipeFds) != 0) {
        filterNote = false;  // fall back to passthrough on failure
    }

    pid_t pid = fork();
    if (pid == 0) {
        if (fullSilence) {
            FILE* devnull = std::fopen("/dev/null", "w");
            if (devnull) dup2(fileno(devnull), STDERR_FILENO);
        } else if (filterNote) {
            close(pipeFds[0]);
            dup2(pipeFds[1], STDERR_FILENO);
            close(pipeFds[1]);
        }
        optind = 1;
        int ret = func(argc, argv);
        _exit(ret);
    } else if (pid > 0) {
        if (filterNote) {
            close(pipeFds[1]);
            FILE* in = fdopen(pipeFds[0], "r");
            if (in) {
                char* line = nullptr;
                size_t cap = 0;
                while (getline(&line, &cap, in) != -1) {
                    if (std::strncmp(line, "Note:", 5) == 0) continue;
                    fputs(line, stderr);
                }
                free(line);
                fclose(in);
            } else {
                close(pipeFds[0]);
            }
        }
        int status;
        waitpid(pid, &status, 0);
        if (WIFEXITED(status) && WEXITSTATUS(status) == 0) return 0;
        return 1;
    } else {
        logging::err("fork() failed for bcftools call");
        if (filterNote) {
            close(pipeFds[0]);
            close(pipeFds[1]);
        }
        return 1;
    }
}

void createMplpBcf(const std::string& prefix,
                   const std::string& refFileName,
                   const std::string& bestMatchSequence,
                   const std::string& bamFileName,
                   std::string& mpileupFileName,
                   bool baq) {
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
            const char* mpileup_args[] = {
                "mpileup", "-Ou", "-f", outRefFileName.c_str(), "-o", mpileupFileName.c_str(), bamFileName.c_str()};
            run_bcftools_in_fork(main_mpileup, 7, const_cast<char**>(mpileup_args));
        } else {
            const char* mpileup_args[] = {"mpileup",
                                          "-Ou",
                                          "-B",
                                          "-f",
                                          outRefFileName.c_str(),
                                          "-o",
                                          mpileupFileName.c_str(),
                                          bamFileName.c_str()};
            run_bcftools_in_fork(main_mpileup, 8, const_cast<char**>(mpileup_args));
        }
    }
}

void createVcfWithMutationMatrices(std::string& prefix,
                                   std::string& mpileupFileName,
                                   std::string& vcfFileName,
                                   const std::vector<std::vector<double>>& substMatrixPhred) {
    if (mpileupFileName.size() == 0) {
        mpileupFileName = prefix + ".mpileup";
    }

    // Use a temp file for raw VCF output (thread-safe via fork isolation)
    std::string rawVcfFile = prefix + ".vcf.raw";
    {
        const char* call_args[] = {
            "call", "--ploidy", "1", "-c", "-A", "-O", "v", "-o", rawVcfFile.c_str(), mpileupFileName.c_str()};
        run_bcftools_in_fork(main_vcfcall, 10, const_cast<char**>(call_args));
    }

    // Read raw VCF and apply mutation spectrum or filter non-ref variants
    std::ifstream rawVcfIn(rawVcfFile);
    bool hasSpectrum = !substMatrixPhred.empty();
    std::ofstream vcfOutFile{vcfFileName};
    std::string line;
    while (std::getline(rawVcfIn, line)) {
        if (hasSpectrum) {
            std::string spectrum_applied_line = genotyping::applyMutationSpectrum(line, substMatrixPhred);
            if (spectrum_applied_line.size() > 0) vcfOutFile << spectrum_applied_line << "\n";
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

int createConsensus(const std::string& vcfFileName,
                    const std::string& refFileName,
                    const std::string& consensusFileName) {
    // bcftools consensus requires a bgzipped + tabix-indexed VCF.
    // Transparently produce a .vcf.gz next to the plain .vcf for the consensus call.
    std::string bgzVcf = vcfFileName + ".gz";

    {
        std::ifstream in(vcfFileName, std::ios::binary);
        if (!in) {
            logging::err("Cannot open VCF for bgzip: {}", vcfFileName);
            return 1;
        }
        BGZF* out = bgzf_open(bgzVcf.c_str(), "w");
        if (!out) {
            logging::err("Cannot create bgzipped VCF: {}", bgzVcf);
            return 1;
        }
        constexpr size_t bufSize = 64 * 1024;
        std::vector<char> buf(bufSize);
        while (in) {
            in.read(buf.data(), bufSize);
            std::streamsize n = in.gcount();
            if (n > 0 && bgzf_write(out, buf.data(), static_cast<size_t>(n)) < 0) {
                logging::err("bgzf_write failed for {}", bgzVcf);
                bgzf_close(out);
                std::remove(bgzVcf.c_str());
                return 1;
            }
        }
        if (bgzf_close(out) < 0) {
            logging::err("bgzf_close failed for {}", bgzVcf);
            std::remove(bgzVcf.c_str());
            return 1;
        }
    }

    if (tbx_index_build(bgzVcf.c_str(), 0, &tbx_conf_vcf) != 0) {
        logging::err("Failed to build tabix index for {}", bgzVcf);
        std::remove(bgzVcf.c_str());
        return 1;
    }

    const char* args[] = {
        "consensus", "-f", refFileName.c_str(), "-o", consensusFileName.c_str(), bgzVcf.c_str()};
    return run_bcftools_in_fork(main_consensus, 6, const_cast<char**>(args), /*silenceStderr=*/false);
}

static uint16_t compute_sam_flags(
    bool is_paired, bool is_read1, uint8_t rev, uint8_t mate_rev, uint8_t proper_frag, bool mate_unmapped) {
    uint16_t flag = 0;
    if (is_paired) {
        flag |= BAM_FPAIRED;
        if (proper_frag) flag |= BAM_FPROPER_PAIR;
        if (rev) flag |= BAM_FREVERSE;
        if (mate_rev) flag |= BAM_FMREVERSE;
        if (mate_unmapped) flag |= BAM_FMUNMAP;
        if (is_read1)
            flag |= BAM_FREAD1;
        else
            flag |= BAM_FREAD2;
    } else {
        if (rev) flag |= BAM_FREVERSE;
    }
    return flag;
}

static int32_t
compute_tlen(int32_t this_rs, int32_t this_re, uint8_t this_rev, int32_t mate_rs, int32_t mate_re, uint8_t mate_rev) {
    int this_pos5 = this_rev ? this_re - 1 : this_rs;
    int mate_pos5 = mate_rev ? mate_re - 1 : mate_rs;
    int tlen = mate_pos5 - this_pos5;
    if (tlen > 0)
        tlen++;
    else if (tlen < 0)
        tlen--;
    return tlen;
}

static bam1_t* build_bam_from_result(const std::string& qname_full,
                                     const std::string& seq,
                                     const std::string& qual,
                                     const read_align_t* aln,
                                     int read_len,
                                     bool is_paired,
                                     bool is_read1,
                                     uint8_t mate_rev,
                                     int32_t mate_pos,
                                     int32_t this_rs,
                                     int32_t this_re,
                                     int32_t mate_rs,
                                     int32_t mate_re,
                                     uint8_t proper_frag,
                                     bool mate_unmapped,
                                     sam_hdr_t* header) {
    bam1_t* b = bam_init1();

    // Strip /1 or /2 suffix from read name (matching minimap2 convention)
    std::string qname = qname_full;
    if (qname.size() >= 2 && qname[qname.size() - 2] == '/' && (qname.back() == '1' || qname.back() == '2')) {
        qname.resize(qname.size() - 2);
    }

    uint16_t flag = compute_sam_flags(is_paired, is_read1, aln->rev, mate_rev, proper_frag, mate_unmapped);

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
                case 'A':
                case 'a': bam_seq[i] = 'T'; break;
                case 'T':
                case 't': bam_seq[i] = 'A'; break;
                case 'C':
                case 'c': bam_seq[i] = 'G'; break;
                case 'G':
                case 'g': bam_seq[i] = 'C'; break;
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
             qname.length(),
             qname.c_str(),
             flag,
             0,                                // tid = 0
             static_cast<hts_pos_t>(aln->rs),  // 0-based position
             aln->mapq,
             full_n_cigar,
             full_cigar.data(),
             mtid,
             mpos_val,
             static_cast<hts_pos_t>(tlen),
             read_len,
             bam_seq.c_str(),
             bam_qual.c_str(),
             0);

    return b;
}

void alignAndWriteBam(std::vector<std::string>& readSequences,
                      std::vector<std::string>& readQuals,
                      std::vector<std::string>& readNames,
                      std::string& reference,
                      const std::string& bamFileName,
                      bool pairedEndReads,
                      int n_threads) {
    int n_reads = static_cast<int>(readSequences.size());

    // Prepare C arrays
    std::vector<const char*> read_ptrs(n_reads);
    std::vector<const char*> qual_ptrs(n_reads);
    std::vector<const char*> name_ptrs(n_reads);
    std::vector<int> r_lens(n_reads);
    for (int i = 0; i < n_reads; i++) {
        read_ptrs[i] = readSequences[i].c_str();
        qual_ptrs[i] = readQuals[i].c_str();
        name_ptrs[i] = readNames[i].c_str();
        r_lens[i] = static_cast<int>(readSequences[i].length());
    }

    int n_results = pairedEndReads ? n_reads / 2 : n_reads;
    std::vector<align_pair_result_t> results(n_results);

    align_reads_direct(reference.c_str(),
                       n_reads,
                       read_ptrs.data(),
                       qual_ptrs.data(),
                       name_ptrs.data(),
                       r_lens.data(),
                       results.data(),
                       pairedEndReads,
                       n_threads);

    std::string samHeaderStr = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:ref\tLN:" + std::to_string(reference.length());
    sam_hdr_t* header = sam_hdr_parse(samHeaderStr.length(), samHeaderStr.c_str());

    std::vector<std::pair<int32_t, bam1_t*>> bam_entries;
    bam_entries.reserve(pairedEndReads ? n_results * 2 : n_results);

    for (int k = 0; k < n_results; k++) {
        align_pair_result_t* res = &results[k];
        if (!res->mapped) continue;

        if (pairedEndReads) {
            int r1_idx = k * 2;
            int r2_idx = k * 2 + 1;

            bam1_t* b1 = build_bam_from_result(readNames[r1_idx],
                                               readSequences[r1_idx],
                                               readQuals[r1_idx],
                                               &res->r1,
                                               r_lens[r1_idx],
                                               true,
                                               true,
                                               res->r2.rev,
                                               res->r2.rs,
                                               res->r1.rs,
                                               res->r1.re,
                                               res->r2.rs,
                                               res->r2.re,
                                               res->r1.proper_frag,
                                               false,
                                               header);
            bam_entries.push_back({res->r1.pos, b1});

            bam1_t* b2 = build_bam_from_result(readNames[r2_idx],
                                               readSequences[r2_idx],
                                               readQuals[r2_idx],
                                               &res->r2,
                                               r_lens[r2_idx],
                                               true,
                                               false,
                                               res->r1.rev,
                                               res->r1.rs,
                                               res->r2.rs,
                                               res->r2.re,
                                               res->r1.rs,
                                               res->r1.re,
                                               res->r2.proper_frag,
                                               false,
                                               header);
            bam_entries.push_back({res->r2.pos, b2});
        } else {
            bam1_t* b = build_bam_from_result(readNames[k],
                                              readSequences[k],
                                              readQuals[k],
                                              &res->r1,
                                              r_lens[k],
                                              false,
                                              false,
                                              0,
                                              -1,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              false,
                                              header);
            bam_entries.push_back({res->r1.pos, b});
        }
    }

    std::sort(bam_entries.begin(), bam_entries.end(), [](const auto& a, const auto& b) { return a.first < b.first; });

    if (!bamFileName.empty()) {
        htsFile* bam_file = hts_open(bamFileName.c_str(), "wb");
        if (bam_file) {
            sam_hdr_write(bam_file, header);
            for (auto& [pos, rec] : bam_entries) {
                bam_write1(bam_file->fp.bgzf, rec);
            }
            hts_close(bam_file);
            logging::info("Wrote bam files to {}", bamFileName);
            if (sam_index_build(bamFileName.c_str(), 0) != 0) {
                logging::warn("Failed to index BAM file: {}", bamFileName);
            }
        } else {
            logging::err("Failed to open output BAM file: {}", bamFileName);
        }
    }

    // Cleanup
    for (auto& [pos, rec] : bam_entries) {
        bam_destroy1(rec);
    }
    sam_hdr_destroy(header);
    for (int k = 0; k < n_results; k++) {
        free(results[k].r1.cigar);
        free(results[k].r2.cigar);
    }
}
