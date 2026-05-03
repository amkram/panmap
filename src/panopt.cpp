/**
 * panopt — field-standard parsimony scoring for single-tree PanMANs.
 *
 * Computes:
 *   - SNV parsimony  (Fitch 1971, Syst Zool 20:406)
 *   - Consistency Index  (Kluge & Farris 1969, Syst Zool 18:1)
 *   - Retention Index    (Farris 1989, Cladistics 5:417)
 *   - Indel & block-level SV event counts
 *
 * For SNV likelihood under GTR+Γ (Felsenstein 1981, J Mol Evol 17:368),
 * pass --export-fasta --export-newick and run IQ-TREE on the fixed
 * topology.
 */

#include <boost/iostreams/filter/lzma.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/program_options.hpp>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "panman.hpp"

namespace po = boost::program_options;

namespace {

constexpr uint8_t kBitA = 1, kBitC = 2, kBitG = 4, kBitT = 8;

inline uint8_t baseBits(char b) {
    switch (b) {
        case 'A': case 'a': return kBitA;
        case 'C': case 'c': return kBitC;
        case 'G': case 'g': return kBitG;
        case 'T': case 't': case 'U': case 'u': return kBitT;
        default: return 0;
    }
}

struct MutCounts {
    long subs = 0;
    long ins_events = 0;
    long ins_bases = 0;
    long del_events = 0;
    long del_bases = 0;
    long block_ins = 0;
    long block_del = 0;
    long block_inv = 0;
    long n_branches = 0;
    long n_leaves = 0;
};

void countMutations(const panmanUtils::Node* node, MutCounts& c, bool isRoot) {
    if (!isRoot) c.n_branches++;
    if (node->children.empty()) c.n_leaves++;

    for (const auto& nm : node->nucMutation) {
        int len = nm.mutInfo >> 4;
        int type = nm.mutInfo & 0xF;
        switch (type) {
            case panmanUtils::NucMutationType::NS:    c.subs += len; break;
            case panmanUtils::NucMutationType::NSNPS: c.subs += 1; break;
            case panmanUtils::NucMutationType::NI:    c.ins_events++; c.ins_bases += len; break;
            case panmanUtils::NucMutationType::NSNPI: c.ins_events++; c.ins_bases += 1; break;
            case panmanUtils::NucMutationType::ND:    c.del_events++; c.del_bases += len; break;
            case panmanUtils::NucMutationType::NSNPD: c.del_events++; c.del_bases += 1; break;
        }
    }
    // Block mutation semantics (per panman.hpp BlockMut + panmap_utils.cpp:131-141):
    //   blockMutInfo=true                  -> insertion (strand normal or inverted at insert)
    //   blockMutInfo=false, inversion=true -> pure inversion (toggle existing strand)
    //   blockMutInfo=false, inversion=false-> deletion
    for (const auto& bm : node->blockMutation) {
        if (bm.blockMutInfo)      c.block_ins++;
        else if (bm.inversion)    c.block_inv++;
        else                       c.block_del++;
    }
    for (auto* child : node->children) countMutations(child, c, false);
}

// Hartigan (1973, J Am Stat Assoc 68:340) postorder — Fitch generalized to
// multifurcations. For binary trees this reduces to Fitch (1971).
// At an internal node: count, across children with observed state, how many
// children's state-set contains each base; the parent's state is the bases
// with maximum count, parsimony increments by (k_observed - max_count).
// Missing children (state=0) are skipped (PAUP*/phangorn convention).
uint8_t fitchPostorder(const panmanUtils::Node* node,
                       const std::unordered_map<std::string, size_t>& leafToRow,
                       const std::vector<std::string>& msa,
                       size_t col,
                       int& parsimony) {
    if (node->children.empty()) {
        auto it = leafToRow.find(node->identifier);
        if (it == leafToRow.end()) return 0;
        return baseBits(msa[it->second][col]);
    }
    int counts[4] = {0, 0, 0, 0};
    int k = 0;
    for (auto* child : node->children) {
        uint8_t s = fitchPostorder(child, leafToRow, msa, col, parsimony);
        if (s == 0) continue;
        k++;
        if (s & kBitA) counts[0]++;
        if (s & kBitC) counts[1]++;
        if (s & kBitG) counts[2]++;
        if (s & kBitT) counts[3]++;
    }
    if (k == 0) return 0;
    int maxc = counts[0];
    for (int i = 1; i < 4; i++) if (counts[i] > maxc) maxc = counts[i];
    uint8_t state = 0;
    if (counts[0] == maxc) state |= kBitA;
    if (counts[1] == maxc) state |= kBitC;
    if (counts[2] == maxc) state |= kBitG;
    if (counts[3] == maxc) state |= kBitT;
    parsimony += (k - maxc);
    return state;
}

std::unique_ptr<panmanUtils::TreeGroup> loadPanMAN(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file) throw std::runtime_error("Cannot open PanMAN file: " + path);
    auto buffer = std::make_unique<boost::iostreams::filtering_streambuf<boost::iostreams::input>>();
    buffer->push(boost::iostreams::lzma_decompressor());
    buffer->push(file);
    std::istream stream(buffer.get());
    return std::make_unique<panmanUtils::TreeGroup>(stream);
}

void parseAlignedFasta(std::istream& in,
                       std::vector<std::string>& seqs,
                       std::unordered_map<std::string, size_t>& nameToRow) {
    std::string line, name, seq;
    auto flush = [&]() {
        if (!name.empty()) {
            nameToRow[name] = seqs.size();
            seqs.push_back(std::move(seq));
            seq.clear();
            name.clear();
        }
    };
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line[0] == '>') {
            flush();
            size_t b = 1;
            while (b < line.size() && (line[b] == ' ' || line[b] == '\t')) b++;
            size_t e = line.find_first_of(" \t", b);
            if (e == std::string::npos) e = line.size();
            name = line.substr(b, e - b);
        } else {
            seq += line;
        }
    }
    flush();
}

}  // namespace

int main(int argc, char** argv) {
    std::string panmanPath, outPrefix;
    int treeIndex = 0;
    bool exportFasta = false, exportNewick = false, skipCi = false;
    long maxLeavesForCi = 100000;

    po::options_description desc("panopt — field-standard parsimony scoring for single-tree PanMANs");
    desc.add_options()
        ("help,h", "Show help")
        ("input", po::value<std::string>(&panmanPath), "Input .panman file (positional)")
        ("output,o", po::value<std::string>(&outPrefix), "Output prefix (default: derived from input)")
        ("tree-index", po::value<int>(&treeIndex)->default_value(0), "Which PanMAT in the PanMAN")
        ("export-fasta", po::bool_switch(&exportFasta), "Export aligned leaf MSA (.leaves.fasta)")
        ("export-newick", po::bool_switch(&exportNewick), "Export topology (.tree.nwk)")
        ("skip-ci", po::bool_switch(&skipCi), "Skip CI/RI (no MSA reconstruction)")
        ("max-leaves-for-ci", po::value<long>(&maxLeavesForCi)->default_value(100000),
            "Auto-skip CI/RI if leaf count exceeds this (memory guard)");

    po::positional_options_description pos;
    pos.add("input", 1);
    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
        po::notify(vm);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n\n" << desc << "\n";
        return 2;
    }

    if (vm.count("help") || panmanPath.empty()) {
        std::cout << "Usage: panopt <input.panman> [options]\n\n" << desc << "\n";
        return panmanPath.empty() ? 2 : 0;
    }

    if (outPrefix.empty()) {
        outPrefix = panmanPath;
        const std::string ext = ".panman";
        if (outPrefix.size() > ext.size() &&
            outPrefix.compare(outPrefix.size() - ext.size(), ext.size(), ext) == 0) {
            outPrefix.resize(outPrefix.size() - ext.size());
        }
    }

    std::cerr << "Loading " << panmanPath << "...\n";
    auto t0 = std::chrono::steady_clock::now();
    std::unique_ptr<panmanUtils::TreeGroup> tg;
    try {
        tg = loadPanMAN(panmanPath);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    if (!tg || tg->trees.empty()) {
        std::cerr << "Error: PanMAN has no trees\n";
        return 1;
    }
    if (treeIndex < 0 || (size_t)treeIndex >= tg->trees.size()) {
        std::cerr << "Error: --tree-index " << treeIndex << " out of range (PanMAN has "
                  << tg->trees.size() << " tree(s))\n";
        return 1;
    }
    if (tg->trees.size() > 1) {
        std::cerr << "Note: PanMAN has " << tg->trees.size()
                  << " trees; scoring tree " << treeIndex << " only.\n";
    }
    if (!tg->complexMutations.empty()) {
        std::cerr << "Note: ignoring " << tg->complexMutations.size()
                  << " complex (network) mutations — panopt is single-tree.\n";
    }

    panmanUtils::Tree& tree = tg->trees[treeIndex];
    auto t1 = std::chrono::steady_clock::now();
    std::cerr << "Loaded in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << " ms\n";

    MutCounts c;
    countMutations(tree.root, c, /*isRoot=*/true);

    // Reference length: try root, but in many PanMATs the root has no blocks
    // present (blocks are inserted on branches), so we also report the median
    // ungapped leaf length from the MSA when available (computed below).
    long rootLen = 0;
    try {
        std::string rootSeq = tree.getStringFromReference(tree.root->identifier, false, true);
        for (char ch : rootSeq) {
            if (ch != '-' && ch != 'N' && ch != 'n') rootLen++;
        }
    } catch (...) {
    }
    long medianLeafLen = 0;

    long P_SNV = c.subs;
    long P_indel_events = c.ins_events + c.del_events;
    long P_indel_bases = c.ins_bases + c.del_bases;
    long P_SV = c.block_ins + c.block_del + c.block_inv;

    bool computedCi = false;
    long P_min_total = 0, P_max_total = 0, P_fitch_total = 0, n_variable = 0;
    long alignmentLen = 0;

    if (skipCi) {
        std::cerr << "Skipping CI/RI (--skip-ci).\n";
    } else if (c.n_leaves > maxLeavesForCi) {
        std::cerr << "Skipping CI/RI: " << c.n_leaves << " leaves > --max-leaves-for-ci "
                  << maxLeavesForCi << ".\n";
    } else if (c.n_leaves == 0) {
        std::cerr << "Skipping CI/RI: tree has no leaves.\n";
    } else {
        std::cerr << "Generating aligned MSA (" << c.n_leaves << " leaves)...\n";
        std::stringstream msaStream;
        bool msaOk = true;
        try {
            tree.printFASTA(msaStream, /*aligned=*/true, /*rootSeq=*/false);
        } catch (const std::exception& e) {
            std::cerr << "Warning: aligned FASTA generation failed (" << e.what()
                      << ") — skipping CI/RI.\n";
            msaOk = false;
        }
        if (msaOk) {
            std::vector<std::string> msa;
            std::unordered_map<std::string, size_t> name2row;
            parseAlignedFasta(msaStream, msa, name2row);
            if (msa.empty()) {
                std::cerr << "Warning: empty MSA — skipping CI/RI.\n";
            } else {
                alignmentLen = (long)msa[0].size();
                bool ragged = false;
                for (const auto& s : msa) {
                    if ((long)s.size() != alignmentLen) { ragged = true; break; }
                }
                if (ragged) {
                    std::cerr << "Warning: MSA sequences differ in length — CI/RI may be unreliable.\n";
                }
                std::cerr << "MSA: " << msa.size() << " × " << alignmentLen
                          << ". Computing per-column Fitch parsimony...\n";

                std::vector<long> leafLens;
                leafLens.reserve(msa.size());
                for (const auto& s : msa) {
                    long n = 0;
                    for (char ch : s) if (baseBits(ch)) n++;
                    leafLens.push_back(n);
                }
                if (!leafLens.empty()) {
                    std::nth_element(leafLens.begin(),
                                     leafLens.begin() + leafLens.size() / 2,
                                     leafLens.end());
                    medianLeafLen = leafLens[leafLens.size() / 2];
                }

                for (long col = 0; col < alignmentLen; col++) {
                    int counts[4] = {0, 0, 0, 0};
                    for (const auto& s : msa) {
                        if ((long)s.size() <= col) continue;
                        uint8_t b = baseBits(s[col]);
                        if (b == kBitA) counts[0]++;
                        else if (b == kBitC) counts[1]++;
                        else if (b == kBitG) counts[2]++;
                        else if (b == kBitT) counts[3]++;
                    }
                    int distinct = 0, max_freq = 0, n_obs = 0;
                    for (int i = 0; i < 4; i++) {
                        if (counts[i] > 0) distinct++;
                        if (counts[i] > max_freq) max_freq = counts[i];
                        n_obs += counts[i];
                    }
                    if (distinct < 2) continue;  // invariant — contributes nothing
                    n_variable++;
                    int colP = 0;
                    fitchPostorder(tree.root, name2row, msa, (size_t)col, colP);
                    P_fitch_total += colP;
                    P_min_total += (distinct - 1);
                    P_max_total += (n_obs - max_freq);
                }
                computedCi = true;
            }
        }
    }

    auto t2 = std::chrono::steady_clock::now();
    long ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t0).count();

    std::cout << "\n=== panopt ===\n";
    std::cout << "Input:               " << panmanPath << "\n";
    std::cout << "Tree index:          " << treeIndex << " of " << tg->trees.size() << "\n";
    std::cout << "Leaves:              " << c.n_leaves << "\n";
    std::cout << "Internal nodes:      " << (c.n_branches + 1 - c.n_leaves) << "\n";
    std::cout << "Branches:            " << c.n_branches << "\n";
    std::cout << "Root ungapped len:   " << rootLen << " bp\n";
    if (medianLeafLen > 0) {
        std::cout << "Median leaf length:  " << medianLeafLen << " bp (ungapped)\n";
    }
    if (computedCi) std::cout << "Alignment length:    " << alignmentLen << " columns\n";

    std::cout << "\n--- Parsimony (from PanMAT branch annotations) ---\n";
    std::cout << "SNV substitutions:   " << P_SNV << "\n";
    std::cout << "Insertions:          " << c.ins_events << " events ("
              << c.ins_bases << " bases)\n";
    std::cout << "Deletions:           " << c.del_events << " events ("
              << c.del_bases << " bases)\n";
    std::cout << "Indel total:         " << P_indel_events << " events ("
              << P_indel_bases << " bases)\n";
    std::cout << "Block insertions:    " << c.block_ins << "\n";
    std::cout << "Block deletions:     " << c.block_del << "\n";
    std::cout << "Block inversions:    " << c.block_inv << "\n";
    std::cout << "SV total (block):    " << P_SV << "\n";

    long denomLen = medianLeafLen > 0 ? medianLeafLen : rootLen;
    const char* denomLabel = medianLeafLen > 0 ? "median-leaf-bp" : "root-bp";
    if (denomLen > 0) {
        std::cout << "\n--- Per-site rates (denominator: " << denomLabel << ") ---\n";
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "SNV / site:          " << (double)P_SNV / denomLen << "\n";
        std::cout << "Indel-events / site: " << (double)P_indel_events / denomLen << "\n";
    }

    if (computedCi) {
        std::cout << "\n--- SNV homoplasy (Fitch on aligned MSA) ---\n";
        std::cout << "Variable sites:      " << n_variable << "\n";
        std::cout << "P_fitch (MSA):       " << P_fitch_total << "\n";
        std::cout << "P_min (Σ(k−1)):      " << P_min_total << "\n";
        std::cout << "P_max (Σ(n−maxfq)):  " << P_max_total << "\n";
        if (P_fitch_total > 0) {
            double CI = (double)P_min_total / (double)P_fitch_total;
            double denom = (double)(P_max_total - P_min_total);
            double RI = denom > 0 ? (double)(P_max_total - P_fitch_total) / denom : 1.0;
            std::cout << std::setprecision(4);
            std::cout << "CI (Kluge-Farris):   " << CI << "\n";
            std::cout << "RI (Farris):         " << RI << "\n";
            std::cout << "RC (CI·RI):          " << (CI * RI) << "\n";
        }
        if (P_fitch_total != P_SNV) {
            std::cout << "\nNote: branch-annotated SNVs (" << P_SNV
                      << ") ≠ MSA-recomputed Fitch parsimony (" << P_fitch_total
                      << "); branch labels may not be Fitch-optimal, or indels in the\n"
                      << "MSA are absorbing some SNV columns. Use P_fitch for CI/RI.\n";
        }
    }

    std::cout << "\nWall time: " << ms << " ms\n";

    std::string tsvPath = outPrefix + ".scores.tsv";
    {
        std::ofstream tsv(tsvPath);
        tsv << "metric\tvalue\n";
        tsv << "n_leaves\t" << c.n_leaves << "\n";
        tsv << "n_branches\t" << c.n_branches << "\n";
        tsv << "root_length_ungapped\t" << rootLen << "\n";
        tsv << "median_leaf_length\t" << medianLeafLen << "\n";
        tsv << "P_SNV_branch\t" << P_SNV << "\n";
        tsv << "P_ins_events\t" << c.ins_events << "\n";
        tsv << "P_ins_bases\t" << c.ins_bases << "\n";
        tsv << "P_del_events\t" << c.del_events << "\n";
        tsv << "P_del_bases\t" << c.del_bases << "\n";
        tsv << "P_indel_events\t" << P_indel_events << "\n";
        tsv << "P_indel_bases\t" << P_indel_bases << "\n";
        tsv << "P_block_ins\t" << c.block_ins << "\n";
        tsv << "P_block_del\t" << c.block_del << "\n";
        tsv << "P_block_inv\t" << c.block_inv << "\n";
        tsv << "P_SV\t" << P_SV << "\n";
        if (computedCi) {
            tsv << "alignment_length\t" << alignmentLen << "\n";
            tsv << "n_variable_sites\t" << n_variable << "\n";
            tsv << "P_fitch_msa\t" << P_fitch_total << "\n";
            tsv << "P_min\t" << P_min_total << "\n";
            tsv << "P_max\t" << P_max_total << "\n";
            if (P_fitch_total > 0) {
                double CI = (double)P_min_total / (double)P_fitch_total;
                double denom = (double)(P_max_total - P_min_total);
                double RI = denom > 0 ? (double)(P_max_total - P_fitch_total) / denom : 1.0;
                tsv << std::setprecision(8);
                tsv << "CI\t" << CI << "\n";
                tsv << "RI\t" << RI << "\n";
                tsv << "RC\t" << (CI * RI) << "\n";
            }
        }
    }
    std::cerr << "Wrote " << tsvPath << "\n";

    if (exportFasta) {
        std::string p = outPrefix + ".leaves.fasta";
        std::ofstream f(p);
        try {
            tree.printFASTA(f, /*aligned=*/true, /*rootSeq=*/false);
            std::cerr << "Wrote " << p << "\n";
        } catch (const std::exception& e) {
            std::cerr << "FASTA export failed: " << e.what() << "\n";
        }
    }
    if (exportNewick) {
        std::string p = outPrefix + ".tree.nwk";
        std::ofstream f(p);
        f << tree.getNewickString(tree.root) << "\n";
        std::cerr << "Wrote " << p << "\n";
    }
    if (exportFasta && exportNewick) {
        std::cerr << "\nFor SNV likelihood (GTR+Γ on fixed topology), run:\n"
                  << "  iqtree -s " << outPrefix << ".leaves.fasta -te "
                  << outPrefix << ".tree.nwk -m GTR+G -n 0 -pre " << outPrefix << ".iqtree\n";
    }

    return 0;
}
