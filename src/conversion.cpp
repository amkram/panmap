#include "conversion.hpp"
#include "genotyping.hpp"
#include "logging.hpp"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <mutex>
#include <unistd.h>
#include <sys/wait.h>
extern "C" {
#include "mm_align.h"
#include "pileup.h"
#include "alignment_options.h"
#include <bcftools/bcftools.h>
#include <htslib/bgzf.h>
#include <htslib/sam.h>
#include <htslib/tbx.h>
}

// Run a bcftools function in a forked child. bcftools main_* use global optind and other
// process-wide state; forking gives each call its own address space for safe parallelism.
static int run_bcftools_in_fork(int (*func)(int, char**), int argc, char** argv, bool silenceStderr = true) {
    // Serialize forks: fork() from a TBB worker clones only the calling thread, so a lock
    // another worker holds (stdio/htslib) can deadlock the child. This mutex bounds it to one
    // bcftools fork in flight (malloc is already fork-safe via atfork).
    static std::mutex forkMutex;
    std::lock_guard<std::mutex> forkLock(forkMutex);

    bool verbose = output::config().verbose;
    bool fullSilence = silenceStderr && !verbose;
    bool filterNote = !silenceStderr && !verbose;  // keep stderr, drop "Note: ..." lines

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

int createMplpBcf(const std::string& prefix,
                  const std::string& refFileName,
                  const std::string& bestMatchSequence,
                  const std::string& bamFileName,
                  std::string& mpileupFileName,
                  bool baq,
                  const std::string& refName) {
    std::string outRefFileName = "";
    if (refFileName.size() == 0) {
        outRefFileName = prefix + ".tmp.reference.fa";
        std::ofstream outRefFile{outRefFileName};
        outRefFile << ">" << refName << "\n";
        outRefFile << bestMatchSequence << "\n";
        outRefFile.close();
    } else {
        outRefFileName = refFileName;
    }

    if (mpileupFileName.size() == 0) {
        mpileupFileName = prefix + ".mpileup";
    }

    int rc = 0;
    {
        if (baq) {
            const char* mpileup_args[] = {
                "mpileup", "-Ou", "-f", outRefFileName.c_str(), "-o", mpileupFileName.c_str(), bamFileName.c_str()};
            rc = run_bcftools_in_fork(main_mpileup, 7, const_cast<char**>(mpileup_args));
        } else {
            const char* mpileup_args[] = {"mpileup",
                                          "-Ou",
                                          "-B",
                                          "-f",
                                          outRefFileName.c_str(),
                                          "-o",
                                          mpileupFileName.c_str(),
                                          bamFileName.c_str()};
            rc = run_bcftools_in_fork(main_mpileup, 8, const_cast<char**>(mpileup_args));
        }
    }
    if (rc != 0) {
        output::error("bcftools mpileup failed (exit code {}) for {}", rc, bamFileName);
        return 1;
    }
    return 0;
}

int createVcfWithMutationMatrices(std::string& prefix,
                                  std::string& mpileupFileName,
                                  std::string& vcfFileName,
                                  const std::vector<std::vector<double>>& substMatrixPhred,
                                  int minDepth,
                                  double minQual) {
    if (mpileupFileName.size() == 0) {
        mpileupFileName = prefix + ".mpileup";
    }

    // Temp file for raw VCF output (thread-safe via fork isolation)
    std::string rawVcfFile = prefix + ".vcf.raw";
    {
        // -m = multiallelic caller (current bcftools default; -c is the legacy
        // consensus caller). --ploidy 1 for a haploid consensus; -A keeps ALT alleles.
        const char* call_args[] = {
            "call", "--ploidy", "1", "-m", "-A", "-O", "v", "-o", rawVcfFile.c_str(), mpileupFileName.c_str()};
        int rc = run_bcftools_in_fork(main_vcfcall, 10, const_cast<char**>(call_args));
        if (rc != 0) {
            output::error("bcftools call failed (exit code {})", rc);
            return 1;
        }
    }

    std::ifstream rawVcfIn(rawVcfFile);
    bool hasSpectrum = !substMatrixPhred.empty();
    std::ofstream vcfOutFile{vcfFileName};
    std::string line;
    while (std::getline(rawVcfIn, line)) {
        if (hasSpectrum) {
            std::string spectrum_applied_line =
                genotyping::applyMutationSpectrum(line, substMatrixPhred, minDepth, minQual);
            if (spectrum_applied_line.size() > 0) vcfOutFile << spectrum_applied_line << "\n";
        } else {
            // No spectrum: write header plus non-ref variants only
            if (line.empty()) continue;
            if (line[0] == '#') {
                vcfOutFile << line << "\n";
                continue;
            }
            std::istringstream iss(line);
            std::string f;
            std::vector<std::string> fields;
            while (std::getline(iss, f, '\t')) fields.push_back(f);
            bool qualOk = fields.size() >= 6 && (fields[5] == "." || std::stod(fields[5]) >= minQual);
            if (fields.size() >= 10 && fields[4] != "." && fields[9][0] != '0' && qualOk &&
                genotyping::passesConsensusGate(fields[9], minDepth)) {
                vcfOutFile << line << "\n";
            }
        }
    }
    rawVcfIn.close();
    std::remove(rawVcfFile.c_str());
    return 0;
}

int createConsensus(const std::string& vcfFileName,
                    const std::string& refFileName,
                    const std::string& consensusFileName,
                    const std::string& consensusHeader) {
    // bcftools consensus requires a bgzipped + tabix-indexed VCF.
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

    const char* args[] = {"consensus", "-f", refFileName.c_str(), "-o", consensusFileName.c_str(), bgzVcf.c_str()};
    int rc = run_bcftools_in_fork(main_consensus, 6, const_cast<char**>(args), /*silenceStderr=*/true);
    if (rc != 0) return rc;

    // bcftools consensus names the output after the reference contig; rename it.
    if (!consensusHeader.empty()) {
        std::string tmpName = consensusFileName + ".tmp";
        std::ifstream in(consensusFileName);
        std::ofstream out(tmpName);
        if (in && out) {
            std::string line;
            bool renamed = false;
            while (std::getline(in, line)) {
                if (!renamed && !line.empty() && line[0] == '>') {
                    out << ">" << consensusHeader << "\n";
                    renamed = true;
                } else {
                    out << line << "\n";
                }
            }
            in.close();
            out.close();
            std::rename(tmpName.c_str(), consensusFileName.c_str());
        }
    }
    return 0;
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
                                     [[maybe_unused]] sam_hdr_t* header) {
    bam1_t* b = bam_init1();

    // Strip /1 or /2 suffix from read name (matching minimap2 convention)
    std::string qname = qname_full;
    if (qname.size() >= 2 && qname[qname.size() - 2] == '/' && (qname.back() == '1' || qname.back() == '2')) {
        qname.resize(qname.size() - 2);
    }

    // R2 was reverse-complemented upstream, so the aligner's rev bit is inverted vs the original
    // fragment. Report true strand in FLAG/TLEN; SEQ/CIGAR stay forward-reference, as BAM expects.
    const uint8_t effective_rev = (is_paired && !is_read1) ? static_cast<uint8_t>(!aln->rev) : aln->rev;
    uint16_t flag = compute_sam_flags(is_paired, is_read1, effective_rev, mate_rev, proper_frag, mate_unmapped);

    // Build full CIGAR with soft clips for unaligned query ends
    uint32_t clip5 = aln->rev ? (read_len - aln->qe) : aln->qs;
    uint32_t clip3 = aln->rev ? aln->qs : (read_len - aln->qe);
    int full_n_cigar = aln->n_cigar + (clip5 > 0 ? 1 : 0) + (clip3 > 0 ? 1 : 0);

    std::vector<uint32_t> full_cigar(full_n_cigar);
    int ci = 0;
    if (clip5 > 0) full_cigar[ci++] = (clip5 << 4) | BAM_CSOFT_CLIP;
    for (int j = 0; j < aln->n_cigar; j++) full_cigar[ci++] = aln->cigar[j];
    if (clip3 > 0) full_cigar[ci++] = (clip3 << 4) | BAM_CSOFT_CLIP;
    


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

    int32_t tlen = 0;
    int32_t mtid = -1;
    hts_pos_t mpos_val = -1;
    if (is_paired) {
        tlen = compute_tlen(this_rs, this_re, effective_rev, mate_rs, mate_re, mate_rev);
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

    bam_aux_update_int(b, "NM", aln->nm);

    if (aln->md) { 
        bam_aux_append(b, "MD", 'Z', static_cast<int>(strlen(aln->md)) + 1,
                       reinterpret_cast<const uint8_t*>(aln->md));
    }
    
    return b;
}

int alignAndWriteBam(std::vector<std::string>& readSequences,
                     std::vector<std::string>& readQuals,
                     std::vector<std::string>& readNames,
                     std::string& reference,
                     const std::string& bamFileName,
                     bool pairedEndReads,
                     int n_threads,
                     const aligner_options_t* alignerOptions,
                     std::string aligner,
                     const std::string& refName) {
    int n_reads = static_cast<int>(readSequences.size());

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

    // "bwa" is bwa aln (aDNA settings); mem is intentionally not offered.
    auto aligner_fn = (aligner == "bwa") ? bwa_aln_align_reads_direct : align_reads_direct;
    aligner_fn(reference.c_str(),
               refName.c_str(),
               n_reads,
               read_ptrs.data(),
               qual_ptrs.data(),
               name_ptrs.data(),
               r_lens.data(),
               results.data(),
               pairedEndReads,
               n_threads,
               alignerOptions);

    std::string samHeaderStr =
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:" + refName + "\tLN:" + std::to_string(reference.length());
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
                                               static_cast<uint8_t>(!res->r2.rev),  // mate (R2) true strand
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

    int rc = 0;
    if (!bamFileName.empty()) {
        htsFile* bam_file = hts_open(bamFileName.c_str(), "wb");
        if (bam_file) {
            if (sam_hdr_write(bam_file, header) < 0) rc = 1;
            for (auto& [pos, rec] : bam_entries) {
                if (bam_write1(bam_file->fp.bgzf, rec) < 0) {
                    rc = 1;
                    break;
                }
            }
            if (hts_close(bam_file) < 0) rc = 1;  // flushes buffered records; a full disk surfaces here
            if (rc != 0) {
                logging::err("Failed to write BAM file: {}", bamFileName);
            } else {
                logging::info("Wrote bam files to {}", bamFileName);
                if (sam_index_build(bamFileName.c_str(), 0) != 0) {
                    logging::warn("Failed to index BAM file: {}", bamFileName);
                }
            }
        } else {
            logging::err("Failed to open output BAM file: {}", bamFileName);
            rc = 1;
        }
    }

    for (auto& [pos, rec] : bam_entries) {
        bam_destroy1(rec);
    }
    sam_hdr_destroy(header);
    for (int k = 0; k < n_results; k++) {
        free(results[k].r1.cigar);
        free(results[k].r2.cigar);
        free(results[k].r1.md);
        free(results[k].r2.md);
    }
    return rc;
}

struct CigarOp {
  uint32_t op;
  uint32_t len;
};


static int consume_query_from_front(std::vector<CigarOp>& ops, int delta) {
  int refConsumed = 0;
  int queryConsumed = 0;
  size_t idx = 0;
  while (queryConsumed < delta && idx < ops.size()) {
    uint32_t op = ops[idx].op;
    bool consumesQuery = op == BAM_CMATCH || op == BAM_CINS || op == BAM_CEQUAL || op == BAM_CDIFF;
    bool consumesRef = op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CEQUAL || op == BAM_CDIFF;
    if (consumesQuery) {
      int take = std::min<int>(ops[idx].len, delta - queryConsumed);
      queryConsumed += take;
      if (consumesRef) refConsumed += take;
      if (take == static_cast<int>(ops[idx].len)) {
        idx++;
      } else {
        ops[idx].len -= take;
        break;
      }
    } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
      refConsumed += ops[idx].len;
      idx++;
    } else {
      idx++;
    }
  }
  ops.erase(ops.begin(), ops.begin() + idx);
  return refConsumed;
}

static void consume_query_from_back(std::vector<CigarOp>& ops, int delta) {
  std::reverse(ops.begin(), ops.end());
  consume_query_from_front(ops, delta);
  std::reverse(ops.begin(), ops.end());
}


static bool compute_clipped_cigar(const uint32_t* origCigar, int nOrigCigar, bool reverse, int mask5, int mask3,
                                  std::vector<CigarOp>& outOps, int& posDelta) {
  std::vector<CigarOp> ops(nOrigCigar);
  for (int i = 0; i < nOrigCigar; i++) {
    ops[i].op = bam_cigar_op(origCigar[i]);
    ops[i].len = bam_cigar_oplen(origCigar[i]);
  }

  int leftTarget = reverse ? mask3 : mask5;
  int rightTarget = reverse ? mask5 : mask3;

  int existingLeft = 0, existingRight = 0;
  if (!ops.empty() && ops.front().op == BAM_CSOFT_CLIP) {
    existingLeft = ops.front().len;
    ops.erase(ops.begin());
  }
  if (!ops.empty() && ops.back().op == BAM_CSOFT_CLIP) {
    existingRight = ops.back().len;
    ops.pop_back();
  }

  posDelta = 0;
  if (leftTarget > existingLeft) {
    posDelta = consume_query_from_front(ops, leftTarget - existingLeft);
    existingLeft = leftTarget;
  }
  if (rightTarget > existingRight) {
    consume_query_from_back(ops, rightTarget - existingRight);
    existingRight = rightTarget;
  }

  bool hasAlignment = std::any_of(ops.begin(), ops.end(), [](const CigarOp& o) {
    return o.op == BAM_CMATCH || o.op == BAM_CEQUAL || o.op == BAM_CDIFF || o.op == BAM_CINS;
  });
  if (!hasAlignment) return false;

  outOps.clear();
  if (existingLeft > 0) outOps.push_back({BAM_CSOFT_CLIP, static_cast<uint32_t>(existingLeft)});
  outOps.insert(outOps.end(), ops.begin(), ops.end());
  if (existingRight > 0) outOps.push_back({BAM_CSOFT_CLIP, static_cast<uint32_t>(existingRight)});
  return true;
}

static bool cigar_unchanged(const uint32_t* orig, int nOrig, const std::vector<CigarOp>& newOps) {
  if (static_cast<size_t>(nOrig) != newOps.size()) return false;
  for (int i = 0; i < nOrig; i++) {
    if (bam_cigar_op(orig[i]) != newOps[i].op || bam_cigar_oplen(orig[i]) != newOps[i].len) return false;
  }
  return true;
}

static void recompute_nm_md(bam1_t* b, const char* ref, int64_t refLen, std::string& mdBuf) {
  uint32_t* cigar = bam_get_cigar(b);
  int nCigar = static_cast<int>(b->core.n_cigar);
  uint8_t* seq = bam_get_seq(b);

  mdBuf.clear();
  int64_t refPos = b->core.pos;
  int qPos = 0;
  int matchRun = 0;
  int nm = 0;

  auto refBase = [&](int64_t pos) -> char { return (pos < refLen) ? ref[pos] : 'N'; };
  // toupper wants an int in [0, UCHAR_MAX] or EOF; a plain (signed) char above 127
  // is UB, hence the unsigned char cast.
  auto refBaseUpper = [&](int64_t pos) -> char {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(refBase(pos))));
  };

  for (int i = 0; i < nCigar; i++) {
    uint32_t op = bam_cigar_op(cigar[i]);
    uint32_t len = bam_cigar_oplen(cigar[i]);
    switch (op) {
      case BAM_CMATCH:
      case BAM_CEQUAL:
      case BAM_CDIFF:
        for (uint32_t j = 0; j < len; j++, qPos++, refPos++) {
          int c1 = bam_seqi(seq, qPos);
          int c2 = seq_nt16_table[static_cast<uint8_t>(refBase(refPos))];
          // Mirrors samtools calmd: c1==c2 alone isn't enough, since 15==15
          // would count two N's as a match; c1==0 covers literal '=' in SEQ.
          bool isMatch = (c1 == c2 && c1 != 15 && c2 != 15) || c1 == 0;
          if (isMatch) {
            matchRun++;
          } else {
            nm++;
            mdBuf += std::to_string(matchRun);
            mdBuf += refBaseUpper(refPos);
            matchRun = 0;
          }
        }
        break;
      case BAM_CINS:
        nm += len;
        qPos += len;
        break;
      case BAM_CDEL:
        nm += len;
        mdBuf += std::to_string(matchRun);
        mdBuf += '^';
        for (uint32_t j = 0; j < len; j++, refPos++) mdBuf += refBaseUpper(refPos);
        matchRun = 0;
        break;
      case BAM_CREF_SKIP:
        refPos += len;
        break;
      case BAM_CSOFT_CLIP:
        qPos += len;
        break;
      default:
        break;
    }
  }
  mdBuf += std::to_string(matchRun);

  bam_aux_update_int(b, "NM", nm);
  bam_aux_update_str(b, "MD", static_cast<int>(mdBuf.size()) + 1, mdBuf.c_str());
}

static void soft_clip_record(bam1_t* b, int mask5, int mask3, std::string& seqBuf, std::string& qualBuf,
                             std::vector<uint32_t>& cigarBuf, const char* ref, int64_t refLen, std::string& mdBuf) {
  int posDelta = 0;
  std::vector<CigarOp> newOps;
  bool stillAligned = compute_clipped_cigar(bam_get_cigar(b), static_cast<int>(b->core.n_cigar), bam_is_rev(b), mask5, mask3, newOps, posDelta);

  uint16_t flag = b->core.flag;
  bool unmap = false;
  if (!stillAligned) {
    newOps.clear();
    posDelta = 0;
    flag |= BAM_FUNMAP;
    unmap = true;
  } else if (posDelta == 0 && cigar_unchanged(bam_get_cigar(b), static_cast<int>(b->core.n_cigar), newOps)) {
    return;
  }

  int l = b->core.l_qseq;
  seqBuf.resize(l);
  uint8_t* seq = bam_get_seq(b);
  for (int i = 0; i < l; i++) seqBuf[i] = seq_nt16_str[bam_seqi(seq, i)];
  qualBuf.assign(reinterpret_cast<char*>(bam_get_qual(b)), l);

  cigarBuf.resize(newOps.size());
  for (size_t i = 0; i < newOps.size(); i++) cigarBuf[i] = bam_cigar_gen(newOps[i].len, newOps[i].op);

  std::string qname(bam_get_qname(b));
  hts_pos_t newPos = unmap ? b->core.pos : b->core.pos + posDelta;

  bam_set1(b, qname.size(), qname.c_str(), flag, b->core.tid, newPos, unmap ? 0 : b->core.qual, cigarBuf.size(),
          cigarBuf.data(), b->core.mtid, b->core.mpos, b->core.isize, l, seqBuf.c_str(), qualBuf.c_str(), 0);

  if (!unmap) recompute_nm_md(b, ref, refLen, mdBuf);
}

static void hard_mask_ends(bam1_t* b, int mask5, int mask3, const char* ref, int64_t refLen, std::string& mdBuf) {
  int len = b->core.l_qseq;
  if (len == 0) return;

  bool rev = bam_is_rev(b);
  int leftCount = rev ? mask3 : mask5;
  int rightCount = rev ? mask5 : mask3;

  if (leftCount < 0) leftCount = 0;
  if (rightCount < 0) rightCount = 0;

  if (leftCount > len) leftCount = len;
  if (rightCount > len - leftCount) rightCount = len - leftCount;

  uint8_t* seq = bam_get_seq(b);
  uint8_t* qual = bam_get_qual(b);
  bool hasQual = qual[0] != 0xff;

  auto setN = [&](int i) {
    // Two 4-bit bases per byte; even i is the high nibble, odd i the low nibble.
    if (i & 1) seq[i >> 1] = (seq[i >> 1] & 0xf0) | 0x0f;
    else       seq[i >> 1] = (seq[i >> 1] & 0x0f) | 0xf0;
    if (hasQual) qual[i] = 0;
  };
  for (int i = 0; i < leftCount; i++) setN(i);
  for (int i = len - rightCount; i < len; i++) setN(i);

  // recompute new NM and MD tags
  recompute_nm_md(b, ref, refLen, mdBuf);
}

int maskBam(const std::string& inputBamFileName,
           const std::string& outputBamFileName,
           int mask5,
           int mask3,
           bool clipMode,
           int nThreads,
           const std::string& reference) {
  samFile* in = sam_open(inputBamFileName.c_str(), "r");
  if (!in) {
    logging::err("Cannot open BAM for masking: {}", inputBamFileName);
    return 1;
  }
  if (nThreads > 1) hts_set_threads(in, nThreads);

  sam_hdr_t* hdr = sam_hdr_read(in);
  if (!hdr) {
    logging::err("Cannot read header from: {}", inputBamFileName);
    sam_close(in);
    return 1;
  }

  std::string seqBuf, qualBuf, mdBuf;
  std::vector<uint32_t> cigarBuf;
  std::vector<std::pair<hts_pos_t, bam1_t*>> entries;
  const char* ref = reference.data();
  int64_t refLen = static_cast<int64_t>(reference.size());

  int rc = 0;
  int ret = 0;
  while (true) {
    bam1_t* b = bam_init1();
    ret = sam_read1(in, hdr, b);
    if (ret < 0) {
      bam_destroy1(b);
      break;
    }
    if (!(b->core.flag & BAM_FUNMAP)) {
      if (clipMode) {
        soft_clip_record(b, mask5, mask3, seqBuf, qualBuf, cigarBuf, ref, refLen, mdBuf);
      } else {
        hard_mask_ends(b, mask5, mask3, ref, refLen, mdBuf);
      }
    }
    entries.push_back({b->core.pos, b});
  }
  if (ret < -1) {
    logging::err("Error reading record while masking: {}", inputBamFileName);
    rc = 1;
  }
  sam_close(in);

  if (rc == 0) {
    // input order can't be assumed sorted anymore.
    std::sort(entries.begin(), entries.end(),
             [](const auto& a, const auto& b) { return a.first < b.first; });

    samFile* out = sam_open(outputBamFileName.c_str(), "wb");
    if (!out) {
      logging::err("Cannot open output BAM for masking: {}", outputBamFileName);
      rc = 1;
    } else {
      if (nThreads > 1) hts_set_threads(out, nThreads);
      if (sam_hdr_write(out, hdr) < 0) {
        logging::err("Failed to write header to: {}", outputBamFileName);
        rc = 1;
      }
      for (auto& [pos, rec] : entries) {
        if (rc == 0 && sam_write1(out, hdr, rec) < 0) {
          logging::err("Failed writing masked record to: {}", outputBamFileName);
          rc = 1;
        }
      }
      if (sam_close(out) < 0) rc = 1;
    }
  }

  for (auto& [pos, rec] : entries) bam_destroy1(rec);
  sam_hdr_destroy(hdr);

  if (rc == 0 && sam_index_build(outputBamFileName.c_str(), 0) != 0) {
    logging::warn("Failed to index masked BAM file: {}", outputBamFileName);
  }
  return rc;
}