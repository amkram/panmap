/**
 * panopt_spr — placement-driven SPR move proposer for single-tree PanMANs.
 *
 * For each leaf L:
 *   1. Reconstruct L's sequence from the PanMAN.
 *   2. Run panmap placement with skipNodeIndex = L (leave-one-out).
 *   3. Compare best non-self placement vs L's current parent score.
 *   4. Emit TSV row: leaf, parent, best_node, scores, delta, would_move.
 *
 * Phase-1 scope: this binary does NOT mutate the PanMAN. The output is the
 * proposal set; actual prune/regraft + Fitch update is a follow-up.
 */

#include <boost/iostreams/filter/lzma.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/program_options.hpp>
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <climits>
#include <atomic>

#include "capnp/message.h"
#include "capnp/serialize.h"
#include "kj/io.h"
#include "kj/std/iostream.h"
#include <boost/iostreams/filter/lzma.hpp>
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include "index_lite.capnp.h"
#include "panman.hpp"
#include "panmap_utils.hpp"
#include "placement.hpp"
#include "zstd_compression.hpp"

#include <spdlog/spdlog.h>

namespace po = boost::program_options;

namespace {

// Mirrors panmap's IndexReader (zstd + capnp lite-index).
class IndexReader : public ::capnp::MessageReader {
   public:
    std::vector<uint8_t> data;
    std::unique_ptr<::capnp::FlatArrayMessageReader> reader;

    explicit IndexReader(const std::string& path, int numThreads = 0)
        : ::capnp::MessageReader(makeOptions()) {
        if (!panmap_zstd::decompressFromFile(path, data, numThreads)) {
            throw std::runtime_error("Failed to decompress index: " + path);
        }
        reader = std::make_unique<::capnp::FlatArrayMessageReader>(
            kj::ArrayPtr<const capnp::word>(reinterpret_cast<const capnp::word*>(data.data()),
                                            data.size() / sizeof(capnp::word)),
            makeOptions());
    }
    kj::ArrayPtr<const capnp::word> getSegment(uint id) override { return reader->getSegment(id); }

   private:
    static ::capnp::ReaderOptions makeOptions() {
        ::capnp::ReaderOptions opts;
        opts.traversalLimitInWords = kj::maxValue;
        opts.nestingLimit = 1024;
        return opts;
    }
};

std::unique_ptr<panmanUtils::TreeGroup> loadPanMAN(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file) throw std::runtime_error("Cannot open PanMAN: " + path);
    auto buffer = std::make_unique<boost::iostreams::filtering_streambuf<boost::iostreams::input>>();
    buffer->push(boost::iostreams::lzma_decompressor());
    buffer->push(file);
    std::istream stream(buffer.get());
    return std::make_unique<panmanUtils::TreeGroup>(stream);
}

void writeFasta(const std::string& path, const std::string& name, const std::string& seq) {
    std::ofstream f(path);
    f << ">" << name << "\n";
    for (size_t i = 0; i < seq.size(); i += 80) {
        f << seq.substr(i, std::min<size_t>(80, seq.size() - i)) << "\n";
    }
}

// Scope guard that swallows std::cout / std::cerr during its lifetime.
// Used around placeLite to silence placement's own log spam (per-call seed
// stats, "Tree traversal:", "Best LogRaw score:", etc).
class SilenceStdio {
    std::streambuf* oldCout;
    std::streambuf* oldCerr;
    std::ofstream sink;
   public:
    SilenceStdio() : sink("/dev/null") {
        oldCout = std::cout.rdbuf(sink.rdbuf());
        oldCerr = std::cerr.rdbuf(sink.rdbuf());
    }
    ~SilenceStdio() {
        std::cout.rdbuf(oldCout);
        std::cerr.rdbuf(oldCerr);
    }
};

// ---- Fitch / CI / RI on (topology, MSA) — copied from panopt.cpp ----

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

// Hartigan 1973 — generalizes Fitch 1971 to multifurcations.
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

struct ParsimonyResult {
    long n_variable = 0;
    long P_fitch = 0;
    long P_min = 0;
    long P_max = 0;
    double CI = 0.0;
    double RI = 0.0;
};

// Compute Fitch parsimony + CI + RI on the current tree topology + MSA.
ParsimonyResult computeParsimony(const panmanUtils::Tree& tree,
                                  const std::vector<std::string>& msa,
                                  const std::unordered_map<std::string, size_t>& name2row) {
    ParsimonyResult r;
    if (msa.empty()) return r;
    long L = (long)msa[0].size();
    for (long col = 0; col < L; col++) {
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
        if (distinct < 2) continue;
        r.n_variable++;
        int colP = 0;
        fitchPostorder(tree.root, name2row, msa, (size_t)col, colP);
        r.P_fitch += colP;
        r.P_min += (distinct - 1);
        r.P_max += (n_obs - max_freq);
    }
    if (r.P_fitch > 0) {
        r.CI = (double)r.P_min / (double)r.P_fitch;
        double denom = (double)(r.P_max - r.P_min);
        r.RI = denom > 0 ? (double)(r.P_max - r.P_fitch) / denom : 1.0;
    }
    return r;
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

// ---- Incremental Fitch cache ---------------------------------------------
// Stores per-(node, variable-column) Fitch state and per-node parsimony
// contribution. After an SPR topology change, recomputes states only on the
// affected ancestor chains (depth-deepest-first). Layout is column-major,
// flat: states_[col * stride + nodeIdx].
class FitchCache {
   public:
    FitchCache(const std::vector<std::string>& msa,
               const std::unordered_map<std::string, size_t>& name2row,
               size_t maxExtraNodes)
        : msa_(msa), name2row_(name2row), maxExtraNodes_(maxExtraNodes) {}

    // Identify variable columns, allocate cache, run full Fitch.
    void initialize(panmanUtils::Tree& tree) {
        // Index every existing node first.
        for (auto& [_, n] : tree.allNodes) registerNode(n);
        size_t nNow = nNodes_;
        stride_ = nNow + maxExtraNodes_;

        // Identify variable columns from MSA.
        if (msa_.empty()) return;
        long L = (long)msa_[0].size();
        variableCols_.reserve(L / 16);
        for (long col = 0; col < L; col++) {
            int counts[4] = {0, 0, 0, 0};
            for (const auto& s : msa_) {
                if ((long)s.size() <= col) continue;
                uint8_t b = baseBits(s[col]);
                if (b == kBitA) counts[0]++;
                else if (b == kBitC) counts[1]++;
                else if (b == kBitG) counts[2]++;
                else if (b == kBitT) counts[3]++;
            }
            int distinct = 0;
            for (int i = 0; i < 4; i++) if (counts[i] > 0) distinct++;
            if (distinct >= 2) variableCols_.push_back(col);
        }
        nVariable_ = variableCols_.size();

        states_.assign(nVariable_ * stride_, 0);
        contrib_.assign(stride_, 0);

        // Full postorder Fitch from root, populating cache.
        std::vector<panmanUtils::Node*> postorder;
        postorderTraversal(tree.root, postorder);
        for (auto* n : postorder) {
            recomputeNode(n);  // populates states + contrib, accumulates totalP_
        }
        // recomputeNode adds delta; on first build, oldContrib is 0, so totalP_ ends up correct.
    }

    long totalParsimony() const { return totalP_; }
    size_t nVariableCols() const { return nVariable_; }

    void registerNode(panmanUtils::Node* node) {
        if (nodeIdx_.count(node)) return;
        nodeIdx_[node] = nNodes_++;
    }

    // Subtract this node's contribution from totalP_ and forget it. Use for nodes
    // about to be deleted (e.g., newInternal during undoSprMove).
    void forgetNode(panmanUtils::Node* node) {
        auto it = nodeIdx_.find(node);
        if (it == nodeIdx_.end()) return;
        size_t idx = it->second;
        totalP_ -= contrib_[idx];
        contrib_[idx] = 0;
        nodeIdx_.erase(it);
    }

    // Update cache and totalP_ after a topology change. Pass the bottom-most
    // affected nodes; we walk parents and recompute deepest-first.
    // Also accepts "orphan" nodes (suppressed parents): we recompute them too
    // so their stale contribution gets subtracted from totalP_.
    //
    // NOTE: we compute true depth-from-root rather than reading node->level,
    // because earlier SPR moves' spr_node insertions invalidate the level
    // field on descendants (we don't propagate). True depth via parent-walk
    // is O(depth) per node — chains are small so this is cheap.
    long updateAfterChange(const std::vector<panmanUtils::Node*>& bottoms) {
        std::unordered_set<panmanUtils::Node*> seen;
        std::vector<panmanUtils::Node*> chain;
        for (auto* b : bottoms) {
            for (auto* n = b; n != nullptr; n = n->parent) {
                if (!seen.insert(n).second) break;
                chain.push_back(n);
            }
        }
        std::unordered_map<panmanUtils::Node*, int> depthCache;
        auto trueDepth = [&](panmanUtils::Node* n) -> int {
            auto it = depthCache.find(n);
            if (it != depthCache.end()) return it->second;
            int d = 0;
            for (auto* p = n; p != nullptr; p = p->parent) d++;
            depthCache[n] = d;
            return d;
        };
        std::sort(chain.begin(), chain.end(),
                  [&](panmanUtils::Node* a, panmanUtils::Node* b) {
                      return trueDepth(a) > trueDepth(b);
                  });
        for (auto* n : chain) recomputeNode(n);
        return totalP_;
    }

   private:
    const std::vector<std::string>& msa_;
    const std::unordered_map<std::string, size_t>& name2row_;
    size_t maxExtraNodes_;

    std::vector<size_t> variableCols_;
    size_t nVariable_ = 0;

    std::unordered_map<panmanUtils::Node*, size_t> nodeIdx_;
    size_t nNodes_ = 0;
    size_t stride_ = 0;
    std::vector<uint8_t> states_;
    std::vector<long> contrib_;
    long totalP_ = 0;

    void postorderTraversal(panmanUtils::Node* node, std::vector<panmanUtils::Node*>& out) {
        for (auto* c : node->children) postorderTraversal(c, out);
        out.push_back(node);
    }

    // Recompute states + contribution at this node across all variable cols.
    // Reads children's cached states. Updates totalP_ by the contribution delta.
    long recomputeNode(panmanUtils::Node* node) {
        if (!nodeIdx_.count(node)) registerNode(node);
        size_t idx = nodeIdx_[node];
        if (idx >= stride_) {
            // Should never happen if maxExtraNodes_ was sized correctly.
            std::cerr << "FitchCache: node index " << idx << " >= stride " << stride_
                      << " (increase maxExtraNodes)\n";
            std::abort();
        }
        long oldContrib = contrib_[idx];
        long newContrib = 0;

        bool isLeaf = node->children.empty();
        size_t leafRow = (size_t)-1;
        if (isLeaf) {
            auto it = name2row_.find(node->identifier);
            if (it != name2row_.end()) leafRow = it->second;
        }

        for (size_t ci = 0; ci < nVariable_; ci++) {
            size_t col = variableCols_[ci];
            uint8_t state;
            int transition = 0;
            if (isLeaf) {
                if (leafRow == (size_t)-1) {
                    state = 0;
                } else {
                    state = baseBits(msa_[leafRow][col]);
                }
            } else {
                int counts[4] = {0, 0, 0, 0};
                int k = 0;
                for (auto* child : node->children) {
                    auto cit = nodeIdx_.find(child);
                    if (cit == nodeIdx_.end()) continue;  // unregistered child
                    uint8_t s = states_[ci * stride_ + cit->second];
                    if (s == 0) continue;
                    k++;
                    if (s & kBitA) counts[0]++;
                    if (s & kBitC) counts[1]++;
                    if (s & kBitG) counts[2]++;
                    if (s & kBitT) counts[3]++;
                }
                if (k == 0) {
                    state = 0;
                } else {
                    int maxc = counts[0];
                    for (int i = 1; i < 4; i++) if (counts[i] > maxc) maxc = counts[i];
                    state = 0;
                    if (counts[0] == maxc) state |= kBitA;
                    if (counts[1] == maxc) state |= kBitC;
                    if (counts[2] == maxc) state |= kBitG;
                    if (counts[3] == maxc) state |= kBitT;
                    transition = (k - maxc);
                }
            }
            states_[ci * stride_ + idx] = state;
            newContrib += transition;
        }

        contrib_[idx] = newContrib;
        long delta = newContrib - oldContrib;
        totalP_ += delta;
        return delta;
    }
};

// ---- SPR move primitives (topology-only) ----

// Records everything needed to undo a single SPR move. When --reject-on-parsimony-increase
// is set, the move is staged; if parsimony got worse, undoSprMove restores the prior tree.
//
// Note: when oldParent is "suppressed" (degree-1 + not root), we DON'T delete it — we just
// orphan it (remove from allNodes, leave heap-allocated). undo re-inserts. If the move is
// kept, the orphaned Node* leaks until process exit; for RSV-4K scale that's <few MB.
struct UndoState {
    panmanUtils::Node* leaf = nullptr;
    panmanUtils::Node* leafOldParent = nullptr;
    int leafIdxInOldParent = -1;

    // Suppression bookkeeping
    bool oldParentSuppressed = false;
    panmanUtils::Node* sibling = nullptr;     // the sibling that took oldParent's slot
    panmanUtils::Node* gp = nullptr;          // oldParent's parent (== sibling's new parent)
    int oldParentIdxInGp = -1;                // where to re-insert oldParent in gp's children
    std::vector<panmanUtils::NucMut> savedSiblingNucMut;
    std::vector<panmanUtils::BlockMut> savedSiblingBlockMut;

    // Insertion bookkeeping
    panmanUtils::Node* newInternal = nullptr;
    panmanUtils::Node* targetParent = nullptr;
    panmanUtils::Node* target = nullptr;
    int targetIdxInTargetParent = -1;
};

// Recursively set node levels from a fixed parent level. Used after SPR moves
// to keep the level field consistent with actual tree depth (panman's
// getNewickString uses n->level to drive paren depth, so stale levels produce
// malformed Newick).
static void resetSubtreeLevels(panmanUtils::Node* n, size_t parentLevel) {
    n->level = parentLevel + 1;
    for (auto* c : n->children) resetSubtreeLevels(c, n->level);
}

// Apply an SPR move; record undo info. Returns false (with `err` set) if rejected.
bool applySprMove(panmanUtils::Tree& tree,
                  panmanUtils::Node* leaf,
                  panmanUtils::Node* target,
                  long& internalCounter,
                  UndoState& u,
                  std::string& err) {
    auto* oldParent = leaf->parent;
    if (!oldParent) { err = "leaf has no parent"; return false; }
    if (target == leaf) { err = "target == leaf"; return false; }
    auto* targetParent = target->parent;
    if (!targetParent) { err = "target is root"; return false; }
    if (target == oldParent) { err = "target == current parent (no-op)"; return false; }

    // Reject if target is in oldParent's subtree (move would attach into own clade)
    {
        std::vector<panmanUtils::Node*> stk{oldParent};
        while (!stk.empty()) {
            auto* n = stk.back(); stk.pop_back();
            if (n == target) { err = "target in oldParent's subtree"; return false; }
            for (auto* c : n->children) stk.push_back(c);
        }
    }

    u.leaf = leaf;
    u.leafOldParent = oldParent;

    // 1. Detach leaf from oldParent (record original index)
    auto& opC = oldParent->children;
    auto it = std::find(opC.begin(), opC.end(), leaf);
    u.leafIdxInOldParent = (int)std::distance(opC.begin(), it);
    if (it != opC.end()) opC.erase(it);

    // 2. Suppress oldParent if degree-1 and not root
    if (opC.size() == 1 && oldParent != tree.root) {
        auto* sibling = opC[0];
        auto* gp = oldParent->parent;
        if (gp) {
            // Save sibling's branch mutations for restore
            u.savedSiblingNucMut = sibling->nucMutation;
            u.savedSiblingBlockMut = sibling->blockMutation;
            u.sibling = sibling;
            u.gp = gp;

            // Find oldParent's index in gp's children, replace with sibling
            auto& gpC = gp->children;
            auto git = std::find(gpC.begin(), gpC.end(), oldParent);
            u.oldParentIdxInGp = (int)std::distance(gpC.begin(), git);
            if (git != gpC.end()) *git = sibling;
            sibling->parent = gp;

            // Concatenate oldParent's branch mutations onto sibling's
            std::vector<panmanUtils::NucMut> merged;
            merged.reserve(oldParent->nucMutation.size() + sibling->nucMutation.size());
            merged.insert(merged.end(), oldParent->nucMutation.begin(), oldParent->nucMutation.end());
            merged.insert(merged.end(), sibling->nucMutation.begin(), sibling->nucMutation.end());
            sibling->nucMutation = std::move(merged);
            std::vector<panmanUtils::BlockMut> mergedB;
            mergedB.reserve(oldParent->blockMutation.size() + sibling->blockMutation.size());
            mergedB.insert(mergedB.end(), oldParent->blockMutation.begin(), oldParent->blockMutation.end());
            mergedB.insert(mergedB.end(), sibling->blockMutation.begin(), sibling->blockMutation.end());
            sibling->blockMutation = std::move(mergedB);

            tree.allNodes.erase(oldParent->identifier);
            // Do NOT delete oldParent — keep it heap-allocated for potential undo.
            u.oldParentSuppressed = true;
            // Sibling moved up one level (was gp's grandchild, now gp's child).
            // Cascade through sibling's subtree.
            resetSubtreeLevels(sibling, gp->level);
        }
    }

    // 3. Insert new internal node N between targetParent and target
    std::string newId = "spr_node_" + std::to_string(++internalCounter);
    auto* N = new panmanUtils::Node(newId, 0.0f);
    N->parent = targetParent;
    N->level = targetParent->level + 1;

    auto& tpC = targetParent->children;
    auto tit = std::find(tpC.begin(), tpC.end(), target);
    u.targetIdxInTargetParent = (int)std::distance(tpC.begin(), tit);
    if (tit != tpC.end()) *tit = N;

    target->parent = N;
    leaf->parent = N;
    N->children.push_back(target);
    N->children.push_back(leaf);
    tree.allNodes[N->identifier] = N;
    // target moved one level deeper (now under N). Cascade through its subtree.
    resetSubtreeLevels(target, N->level);
    leaf->level = N->level + 1;

    u.newInternal = N;
    u.targetParent = targetParent;
    u.target = target;
    return true;
}

// ---- Column → (blockId, nucPos, gapPos) coordinate map -------------------
// The aligned MSA's column ordering follows panman's block→gaps→main scheme
// (see panmap_utils.cpp:getStringFromSequence).
// gapPos = -1 means "main nucleotide at (blockId, nucPos)".
// gapPos = 0..L-1 means "k-th gap-region char before that main nuc".
struct ColCoord {
    int32_t blockId;
    int32_t nucPos;
    int32_t gapPos;
};

// Walk tree.blocks/gaps to enumerate column → coord. n_cols_out is the total
// alignment length (must match MSA[0].size()).
// Mirrors panman's fasta.cpp:printSequenceLines emit order:
//   For each primaryBlockId i: emit secondary blocks i.second[0..k] in order, then primary i.first.
//   For each block: for each nucPos 0..nucCount inclusive, emit gap chars at that nucPos
//                   then the main char (only if nucPos < nucCount).
std::vector<ColCoord> buildColMap(const panmanUtils::Tree& tree, long& n_cols_out) {
    // (primaryBlockId, secondaryBlockId+1) -> nucleotide count from consensusSeq.
    auto bkey = [](int primary, int secondary) -> long long {
        return ((long long)primary << 40) | ((long long)(secondary + 1) << 20);
    };
    auto gkey = [](int primary, int secondary, int nucPos) -> long long {
        return ((long long)primary << 40) | ((long long)(secondary + 1) << 20) | (long long)nucPos;
    };

    std::unordered_map<long long, int> blockNucCount;
    int maxPrimary = -1;
    for (const auto& blk : tree.blocks) {
        int nucCount = 0;
        for (size_t i = 0; i < blk.consensusSeq.size(); i++) {
            bool stop = false;
            for (size_t j = 0; j < 8; j++) {
                int code = (blk.consensusSeq[i] >> (4 * (7 - j))) & 15;
                if (code == 0) { stop = true; break; }
                nucCount++;
            }
            if (stop) break;
        }
        blockNucCount[bkey(blk.primaryBlockId, blk.secondaryBlockId)] = nucCount;
        if (blk.primaryBlockId > maxPrimary) maxPrimary = blk.primaryBlockId;
    }

    std::unordered_map<long long, int> gapLen;
    for (const auto& gl : tree.gaps) {
        for (size_t j = 0; j < gl.nucPosition.size(); j++) {
            gapLen[gkey(gl.primaryBlockId, gl.secondaryBlockId, gl.nucPosition[j])] = gl.nucGapLength[j];
        }
    }

    // # of secondary slots per primary position (from blockGaps).
    std::unordered_map<int, int> nSecPerPrimary;
    for (size_t i = 0; i < tree.blockGaps.blockPosition.size(); i++) {
        nSecPerPrimary[tree.blockGaps.blockPosition[i]] = tree.blockGaps.blockGapLength[i];
    }

    std::vector<ColCoord> colMap;
    long blocksWithCount = 0, secondariesEmitted = 0;
    auto emitBlock = [&](int primaryBlockId, int secondaryBlockId) {
        long long k = bkey(primaryBlockId, secondaryBlockId);
        auto it = blockNucCount.find(k);
        if (it == blockNucCount.end()) return;  // empty slot — no cols
        int nucCount = it->second;
        blocksWithCount++;
        if (secondaryBlockId >= 0) secondariesEmitted++;
        // Iterate nucPos 0..nucCount INCLUSIVE. The 'x' sentinel at nucPos==nucCount
        // emits a single '-' column in aligned mode (fasta.cpp:78-82), so we DO push
        // a main entry at every nucPos, including the sentinel.
        for (int32_t nucPos = 0; nucPos <= nucCount; nucPos++) {
            long long gk = gkey(primaryBlockId, secondaryBlockId, nucPos);
            int gl = gapLen.count(gk) ? gapLen[gk] : 0;
            for (int gp = 0; gp < gl; gp++) {
                colMap.push_back({(int32_t)primaryBlockId, nucPos, gp});
            }
            colMap.push_back({(int32_t)primaryBlockId, nucPos, -1});
        }
    };

    for (int blockId = 0; blockId <= maxPrimary; blockId++) {
        // Secondary slots first
        int nSec = nSecPerPrimary.count(blockId) ? nSecPerPrimary[blockId] : 0;
        for (int s = 0; s < nSec; s++) emitBlock(blockId, s);
        // Then primary
        emitBlock(blockId, -1);
    }

    std::cerr << "  buildColMap: " << blocksWithCount << " blocks emitted ("
              << secondariesEmitted << " secondary), " << colMap.size() << " cols\n";
    n_cols_out = (long)colMap.size();
    return colMap;
}

// 5-state encoding: gap as a 5th bit. For Fitch on the regeneration pass.
inline uint8_t baseBits5(char b) {
    switch (b) {
        case 'A': case 'a': return 0x01;
        case 'C': case 'c': return 0x02;
        case 'G': case 'g': return 0x04;
        case 'T': case 't': case 'U': case 'u': return 0x08;
        case '-':           return 0x10;  // gap as 5th state
        default:            return 0;
    }
}
inline char charFromBits5(uint8_t s) {
    if (s & 0x01) return 'A';
    if (s & 0x02) return 'C';
    if (s & 0x04) return 'G';
    if (s & 0x08) return 'T';
    if (s & 0x10) return '-';
    return 'N';
}
inline int panmanCodeFromChar(char c) {
    switch (c) {
        case 'A': return 1;
        case 'C': return 2;
        case 'G': return 4;
        case 'T': case 'U': return 8;
        default: return 0;
    }
}

// BFS edge-distance from `a` to `b`, capped at `maxR`. Returns -1 if > maxR.
// Used by --local-radius to constrain SPR moves to nearby branches.
int bfsDistance(panmanUtils::Node* a, panmanUtils::Node* b, int maxR) {
    if (a == b) return 0;
    std::unordered_map<panmanUtils::Node*, int> dist;
    dist[a] = 0;
    std::vector<panmanUtils::Node*> frontier{a};
    while (!frontier.empty()) {
        std::vector<panmanUtils::Node*> next;
        for (auto* n : frontier) {
            int d = dist[n];
            if (d >= maxR) continue;
            std::vector<panmanUtils::Node*> neigh;
            if (n->parent) neigh.push_back(n->parent);
            for (auto* c : n->children) neigh.push_back(c);
            for (auto* nb : neigh) {
                if (dist.count(nb)) continue;
                dist[nb] = d + 1;
                if (nb == b) return d + 1;
                next.push_back(nb);
            }
        }
        frontier = std::move(next);
    }
    return -1;
}

// Regenerate per-branch nucMutation lists for the entire tree from the MSA + topology.
// Uses 5-state Fitch (gap as a 5th state). After this runs, the tree's NucMut chain is
// consistent: getStringFromReference(leaf) reproduces the MSA sequence, and writeToFile
// produces a valid PanMAN.
//
// Block mutations are preserved as-is — we only regenerate nucleotide-level mutations.
//
// Implementation: column-streaming. For each MSA column we run Fitch postorder + preorder
// to pick concrete states and emit any NucMut at that position. Memory O(n_nodes); avoids
// the O(n_nodes × n_cols) state matrix that would be ~8 GB on RSV-4K.
void regenerateBranchMutations(panmanUtils::Tree& tree,
                                const std::vector<std::string>& msa,
                                const std::unordered_map<std::string, size_t>& name2row,
                                const std::vector<ColCoord>& colMap) {
    long n_cols = (long)(msa.empty() ? 0 : msa[0].size());
    if (n_cols == 0 || (long)colMap.size() != n_cols) {
        std::cerr << "regenerateBranchMutations: MSA/colMap size mismatch ("
                  << n_cols << " vs " << colMap.size() << ")\n";
        return;
    }

    // Output: per-node mutation lists, accumulated across columns.
    using NodeMutsMap = std::unordered_map<panmanUtils::Node*, std::vector<panmanUtils::NucMut>>;
    using NodeStateMap = std::unordered_map<panmanUtils::Node*, uint8_t>;
    NodeMutsMap nodeMuts;

    auto pickLowBit = [](uint8_t st) -> uint8_t {
        for (int b = 0; b < 5; b++) if (st & (1 << b)) return (uint8_t)(1 << b);
        return 0;
    };

    // Process one MSA column on a private state-map and per-thread muts-map.
    auto processCol = [&](long c, NodeStateMap& states, NodeMutsMap& localMuts) {
        states.clear();
        // Postorder Fitch (5-state). Stores state per node in `states`.
        std::function<uint8_t(panmanUtils::Node*)> post = [&](panmanUtils::Node* n) -> uint8_t {
            uint8_t s;
            if (n->children.empty()) {
                auto it = name2row.find(n->identifier);
                s = (it != name2row.end() && c < (long)msa[it->second].size())
                    ? baseBits5(msa[it->second][c]) : 0;
            } else {
                int counts[5] = {0,0,0,0,0};
                int k = 0;
                for (auto* ch : n->children) {
                    uint8_t cs = post(ch);
                    if (cs == 0) continue;
                    k++;
                    if (cs & 0x01) counts[0]++;
                    if (cs & 0x02) counts[1]++;
                    if (cs & 0x04) counts[2]++;
                    if (cs & 0x08) counts[3]++;
                    if (cs & 0x10) counts[4]++;
                }
                if (k == 0) s = 0;
                else {
                    int maxc = counts[0];
                    for (int i = 1; i < 5; i++) if (counts[i] > maxc) maxc = counts[i];
                    s = 0;
                    if (counts[0] == maxc) s |= 0x01;
                    if (counts[1] == maxc) s |= 0x02;
                    if (counts[2] == maxc) s |= 0x04;
                    if (counts[3] == maxc) s |= 0x08;
                    if (counts[4] == maxc) s |= 0x10;
                }
            }
            states[n] = s;
            return s;
        };
        uint8_t rootState = post(tree.root);
        uint8_t rootPick = rootState ? pickLowBit(rootState) : 0;

        // Preorder pick + emit.
        std::function<void(panmanUtils::Node*, uint8_t)> pre =
            [&](panmanUtils::Node* n, uint8_t parentPick) {
                uint8_t mySet = states[n];
                uint8_t pick;
                if (mySet == 0) pick = parentPick ? parentPick : 0x10;
                else if (mySet & parentPick) pick = parentPick;
                else pick = pickLowBit(mySet);
                if (n != tree.root && pick != parentPick) {
                    const ColCoord& cc = colMap[c];
                    char pickChar = charFromBits5(pick);
                    char parentChar = charFromBits5(parentPick);
                    int type = -1, code = 0;
                    if (parentChar == '-' && pickChar != '-') {
                        type = panmanUtils::NucMutationType::NSNPI;
                        code = panmanCodeFromChar(pickChar);
                    } else if (parentChar != '-' && pickChar == '-') {
                        type = panmanUtils::NucMutationType::NSNPD;
                    } else if (parentChar != '-' && pickChar != '-' && pickChar != 'N') {
                        type = panmanUtils::NucMutationType::NSNPS;
                        code = panmanCodeFromChar(pickChar);
                    }
                    if (type >= 0) {
                        localMuts[n].push_back(panmanUtils::NucMut(
                            std::make_tuple((int)cc.blockId, -1,
                                            (int)cc.nucPos, (int)cc.gapPos,
                                            type, code)));
                    }
                }
                for (auto* ch : n->children) pre(ch, pick);
            };
        for (auto* child : tree.root->children) pre(child, rootPick);
    };

    auto t0 = std::chrono::steady_clock::now();

    // ---- Parallel column processing ----
    // Each thread keeps its own colStates + nodeMuts. Cols within a thread are
    // processed in increasing order; after parallel_for we merge per-node lists
    // and sort by (block, sec, nucPos, gapPos) for consolidation.
    tbb::enumerable_thread_specific<NodeStateMap> tlsStates;
    tbb::enumerable_thread_specific<NodeMutsMap> tlsMuts;
    std::atomic<long> colsDone{0};

    tbb::parallel_for(
        tbb::blocked_range<long>(0, n_cols, 4096),
        [&](const tbb::blocked_range<long>& r) {
            auto& states = tlsStates.local();
            auto& muts = tlsMuts.local();
            for (long c = r.begin(); c < r.end(); c++) {
                processCol(c, states, muts);
            }
            long done = colsDone.fetch_add(r.end() - r.begin()) + (r.end() - r.begin());
            if ((done / 100000) != ((done - (r.end() - r.begin())) / 100000)) {
                auto now = std::chrono::steady_clock::now();
                long ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - t0).count();
                std::cerr << "  regenerate: ~" << done << "/" << n_cols
                          << " cols (" << (ms / 1000) << "s)\n";
            }
        });

    // Merge thread-local mut maps into nodeMuts.
    for (auto& tlm : tlsMuts) {
        for (auto& kv : tlm) {
            auto& dest = nodeMuts[kv.first];
            dest.insert(dest.end(),
                        std::make_move_iterator(kv.second.begin()),
                        std::make_move_iterator(kv.second.end()));
        }
        tlm.clear();
    }
    // Sort each node's muts by (block, sec, nucPos, gapPos with -1 last) — required
    // for bulk consolidation to find consecutive same-type runs.
    for (auto& kv : nodeMuts) {
        std::sort(kv.second.begin(), kv.second.end(),
                  [](const panmanUtils::NucMut& a, const panmanUtils::NucMut& b) {
                      if (a.primaryBlockId != b.primaryBlockId) return a.primaryBlockId < b.primaryBlockId;
                      if (a.secondaryBlockId != b.secondaryBlockId) return a.secondaryBlockId < b.secondaryBlockId;
                      if (a.nucPosition != b.nucPosition) return a.nucPosition < b.nucPosition;
                      int ag = a.nucGapPosition < 0 ? INT_MAX : a.nucGapPosition;
                      int bg = b.nucGapPosition < 0 ? INT_MAX : b.nucGapPosition;
                      return ag < bg;
                  });
    }

    // ---- Consolidate single-position NucMuts into bulk (NS/NI/ND, up to length 6) ----
    // panman stores up to 6 consecutive same-type same-block changes in one NucMut.
    // Two entries are position-consecutive (in moveForward sense) iff:
    //   both gapPos==-1 and nucPos differs by 1 (and same block); OR
    //   same nucPos and same block, gap-consecutive (gapPos differs by 1).
    auto singleToBulkType = [](int t) -> int {
        switch (t) {
            case panmanUtils::NucMutationType::NSNPS: return panmanUtils::NucMutationType::NS;
            case panmanUtils::NucMutationType::NSNPI: return panmanUtils::NucMutationType::NI;
            case panmanUtils::NucMutationType::NSNPD: return panmanUtils::NucMutationType::ND;
            default: return -1;
        }
    };

    auto consolidate = [&](std::vector<panmanUtils::NucMut>& muts) {
        if (muts.size() < 2) return;
        std::vector<panmanUtils::NucMut> out;
        out.reserve(muts.size());

        size_t i = 0;
        while (i < muts.size()) {
            const auto& start = muts[i];
            int startType = start.mutInfo & 0xF;
            int bulkType = singleToBulkType(startType);
            // Length-cap groups (up to 6 entries).
            size_t j = i + 1;
            int prevPos = start.nucPosition;
            int prevGap = start.nucGapPosition;
            std::vector<int> codes;
            codes.push_back((start.nucs >> 20) & 0xF);
            while (j < muts.size() && (j - i) < 6 && bulkType != -1) {
                const auto& cur = muts[j];
                int curType = cur.mutInfo & 0xF;
                int curBulkType = singleToBulkType(curType);
                if (curBulkType != bulkType) break;
                if (cur.primaryBlockId != start.primaryBlockId) break;
                if (cur.secondaryBlockId != start.secondaryBlockId) break;
                bool consec = false;
                if (start.nucGapPosition < 0 && cur.nucGapPosition < 0
                    && cur.nucPosition == prevPos + 1) consec = true;
                else if (start.nucGapPosition >= 0 && cur.nucPosition == prevPos
                         && cur.nucGapPosition == prevGap + 1) consec = true;
                if (!consec) break;
                codes.push_back((cur.nucs >> 20) & 0xF);
                prevPos = cur.nucPosition;
                prevGap = cur.nucGapPosition;
                j++;
            }
            int L = (int)codes.size();
            if (L == 1 || bulkType == -1) {
                out.push_back(start);
            } else {
                std::vector<std::tuple<int, int, int, int, int, int>> arr;
                arr.reserve(L);
                for (int k = 0; k < L; k++) {
                    arr.emplace_back(start.primaryBlockId, start.secondaryBlockId,
                                     start.nucPosition, start.nucGapPosition,
                                     startType, codes[k]);
                }
                out.emplace_back(arr, 0, L);
            }
            i = j;
        }
        muts = std::move(out);
    };

    // Replace nucMutation lists. Nodes not in nodeMuts get cleared.
    long nodes_with_muts = 0, total_muts = 0, total_consolidated = 0;
    for (auto& kv : tree.allNodes) {
        auto* n = kv.second;
        if (n == tree.root) {
            n->nucMutation.clear();
            continue;
        }
        auto it = nodeMuts.find(n);
        if (it != nodeMuts.end()) {
            long pre = (long)it->second.size();
            consolidate(it->second);
            n->nucMutation = std::move(it->second);
            nodes_with_muts++;
            total_muts += pre;
            total_consolidated += (long)n->nucMutation.size();
        } else {
            n->nucMutation.clear();
        }
    }
    auto t1 = std::chrono::steady_clock::now();
    std::cerr << "  regenerate: " << total_muts << " single-pos NucMuts → "
              << total_consolidated << " bulk-consolidated, on "
              << nodes_with_muts << " branches  ("
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
              << " ms)\n";
}

// Reverse a successful applySprMove.
void undoSprMove(panmanUtils::Tree& tree, UndoState& u) {
    // Reverse step 3: remove newInternal from targetParent, restore target in its slot
    auto& tpC = u.targetParent->children;
    tpC.erase(std::remove(tpC.begin(), tpC.end(), u.newInternal), tpC.end());
    if (u.targetIdxInTargetParent >= 0 && (size_t)u.targetIdxInTargetParent <= tpC.size()) {
        tpC.insert(tpC.begin() + u.targetIdxInTargetParent, u.target);
    } else {
        tpC.push_back(u.target);
    }
    u.target->parent = u.targetParent;
    // target moved one level shallower again.
    resetSubtreeLevels(u.target, u.targetParent->level);
    tree.allNodes.erase(u.newInternal->identifier);
    delete u.newInternal;
    u.newInternal = nullptr;

    // Reverse step 2: if oldParent was suppressed, un-suppress
    if (u.oldParentSuppressed) {
        // Restore sibling's branch mutations
        u.sibling->nucMutation = std::move(u.savedSiblingNucMut);
        u.sibling->blockMutation = std::move(u.savedSiblingBlockMut);
        // Reattach sibling under oldParent (sibling moves one level deeper again).
        u.sibling->parent = u.leafOldParent;
        resetSubtreeLevels(u.sibling, u.leafOldParent->level);
        // oldParent's children should be just {sibling} (we left it that way after detach)
        u.leafOldParent->children.clear();
        u.leafOldParent->children.push_back(u.sibling);
        // Replace sibling in gp's children with oldParent (at its original index)
        auto& gpC = u.gp->children;
        auto sit = std::find(gpC.begin(), gpC.end(), u.sibling);
        if (sit != gpC.end()) *sit = u.leafOldParent;
        // Restore oldParent in allNodes
        tree.allNodes[u.leafOldParent->identifier] = u.leafOldParent;
        u.leafOldParent->parent = u.gp;
    }

    // Reverse step 1: re-attach leaf to oldParent at its original index
    auto& opC = u.leafOldParent->children;
    if (u.leafIdxInOldParent >= 0 && (size_t)u.leafIdxInOldParent <= opC.size()) {
        opC.insert(opC.begin() + u.leafIdxInOldParent, u.leaf);
    } else {
        opC.push_back(u.leaf);
    }
    u.leaf->parent = u.leafOldParent;
}

}  // namespace

int main(int argc, char** argv) {
    std::string panmanPath, indexPath, outPrefix;
    int treeIndex = 0, threads = 1, maxLeaves = -1, progressEvery = 100;
    double margin = 0.0;

    po::options_description desc("panopt_spr — placement-driven SPR move proposer");
    desc.add_options()
        ("help,h", "Show help")
        ("panman,p", po::value<std::string>(&panmanPath), "Input .panman file")
        ("index,i", po::value<std::string>(&indexPath), "Pre-built panmap index (.idx)")
        ("output,o", po::value<std::string>(&outPrefix)->default_value("spr"), "Output prefix")
        ("tree-index", po::value<int>(&treeIndex)->default_value(0), "Which PanMAT in PanMAN")
        ("threads,t", po::value<int>(&threads)->default_value(1), "Threads (passed to placement)")
        ("max-leaves,n", po::value<int>(&maxLeaves)->default_value(-1),
            "Cap leaves processed (-1 = all). Useful for testing.")
        ("margin,m", po::value<double>(&margin)->default_value(0.0),
            "Min log-containment delta to flag would_move (default 0)")
        ("progress-every", po::value<int>(&progressEvery)->default_value(500),
            "Print mid-pass progress every N leaves (0 = only at end of each pass)")
        ("placement-verbose", po::bool_switch(),
            "Show placement's internal log output (suppressed by default)")
        ("apply-moves", po::bool_switch(),
            "Actually perform SPR moves on the in-memory tree (Phase 2)")
        ("score-every", po::value<int>(),
            "Compute panopt parsimony every N moves and append to .trace.tsv")
        ("reject-on-parsimony-increase", po::bool_switch(),
            "After each apply, recompute parsimony; revert if it increased.")
        ("max-passes", po::value<int>()->default_value(10),
            "Max optimization passes per phase (default 10)")
        ("local-radius", po::value<int>()->default_value(0),
            "After global converges, do local passes restricting each move to "
            "branches within R BFS hops of leaf's old parent. 0 = skip local phase.")
        ("refresh-index", po::bool_switch(),
            "Rebuild .panman + .idx between passes (subprocesses panmanUtils + panmap). "
            "Required for multi-pass to find new improvements; without it, later passes "
            "see stale placement scores.")
        ("panmap-bin", po::value<std::string>()->default_value("panmap"),
            "panmap executable for index refresh")
        ("panman-utils-bin", po::value<std::string>()->default_value("panmanUtils"),
            "panmanUtils executable for index refresh")
        ("write-output-panman", po::bool_switch(),
            "After all SPR passes: regenerate per-branch NucMuts via Fitch on (topology, MSA), "
            "then write a valid optimized .panman to <output>.optimized.panman");

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
        if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
        po::notify(vm);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n\n" << desc << "\n";
        return 2;
    }
    if (panmanPath.empty() || indexPath.empty()) {
        std::cerr << "Error: --panman and --index are required\n\n" << desc << "\n";
        return 2;
    }

    // Silence placement's per-call info logging.
    spdlog::set_pattern("%v");
    spdlog::set_level(spdlog::level::err);

    // Set TBB max parallelism — placement, regenerate, and any TBB-based code below
    // honors this. Default 1 means serial; user passes -t 12 etc.
    tbb::global_control tbbCtl(tbb::global_control::max_allowed_parallelism,
                                std::max(1, threads));

    // ---- Flags ----
    bool applyMoves = vm["apply-moves"].as<bool>();
    int scoreEvery = vm.count("score-every") ? vm["score-every"].as<int>() : 0;
    bool rejectOnIncrease = vm["reject-on-parsimony-increase"].as<bool>();
    int maxPasses = vm["max-passes"].as<int>();
    int localRadius = vm["local-radius"].as<int>();
    bool refreshIndex = vm["refresh-index"].as<bool>();
    bool placementVerbose = vm["placement-verbose"].as<bool>();
    std::string panmapBin = vm["panmap-bin"].as<std::string>();
    std::string panmanUtilsBin = vm["panman-utils-bin"].as<std::string>();
    if (rejectOnIncrease && !applyMoves) {
        std::cerr << "--reject-on-parsimony-increase requires --apply-moves\n";
        return 2;
    }

    // ---- State that lives across passes (rebuilt from in-memory tree state) ----
    std::unique_ptr<panmanUtils::TreeGroup> tg;
    panmanUtils::Tree* tree = nullptr;
    std::unique_ptr<IndexReader> reader;
    panmapUtils::LiteTree liteTree;
    std::vector<panmanUtils::Node*> leaves;
    int idxK = 0, idxS = 0, idxL = 0, idxT = 0;
    bool idxOpen = false;
    long internalCounter = 0;
    long currentPFitch = -1;
    std::unique_ptr<FitchCache> fcache;

    // Persistent (computed once)
    std::vector<std::string> msa;
    std::unordered_map<std::string, size_t> name2row;
    long P_min_const = -1, P_max_const = -1;
    int n_var_const = 0;
    std::string msaFastaPath;  // file path; reused across refreshes

    // Lambda: (re)load tg/tree/reader/liteTree/leaves from given .panman + .idx paths.
    auto loadAll = [&](const std::string& pmPath, const std::string& idxPath) {
        fcache.reset();
        leaves.clear();
        liteTree = panmapUtils::LiteTree();
        reader.reset();
        tree = nullptr;
        tg.reset();

        tg = loadPanMAN(pmPath);
        if (!tg || tg->trees.empty() || treeIndex < 0 || (size_t)treeIndex >= tg->trees.size()) {
            throw std::runtime_error("invalid tree-index or empty PanMAN: " + pmPath);
        }
        tree = &tg->trees[treeIndex];

        reader = std::make_unique<IndexReader>(idxPath, threads);
        auto idx = reader->getRoot<LiteIndex>();
        idxK = idx.getK(); idxS = idx.getS(); idxL = idx.getL();
        idxT = idx.getT(); idxOpen = idx.getOpen();
        liteTree.initialize(idx.getLiteTree());

        leaves.reserve(tree->allNodes.size() / 2);
        for (auto& kv : tree->allNodes) {
            if (kv.second->children.empty()) leaves.push_back(kv.second);
        }
        std::sort(leaves.begin(), leaves.end(),
                  [](panmanUtils::Node* a, panmanUtils::Node* b) {
                      return a->identifier < b->identifier;
                  });
        if (maxLeaves > 0 && (size_t)maxLeaves < leaves.size()) leaves.resize(maxLeaves);
    };

    // Lambda: rebuild PanMAN+index from current Newick + cached MSA, then loadAll.
    auto refresh = [&](int passN) -> bool {
        std::string nwkPath = outPrefix + ".pass" + std::to_string(passN) + ".nwk";
        {
            std::ofstream f(nwkPath);
            f << tree->getNewickString(tree->root) << "\n";
        }
        std::string newPrefix = outPrefix + ".pass" + std::to_string(passN);
        std::string newPanman = newPrefix + ".panman";
        std::string newIdx = newPrefix + ".idx";

        std::string cmd1 = panmanUtilsBin + " -M " + msaFastaPath + " -N " + nwkPath
                           + " -o " + newPrefix + " 1>&2";
        std::cerr << "  $ " << cmd1 << "\n";
        if (system(cmd1.c_str()) != 0) {
            std::cerr << "  panmanUtils failed\n";
            return false;
        }
        std::string cmd2 = panmapBin + " " + newPanman + " --stop index -o " + newPrefix
                           + " -t " + std::to_string(threads) + " 1>&2";
        std::cerr << "  $ " << cmd2 << "\n";
        if (system(cmd2.c_str()) != 0) {
            std::cerr << "  panmap index build failed\n";
            return false;
        }
        try { loadAll(newPanman, newIdx); }
        catch (const std::exception& e) {
            std::cerr << "  reload failed: " << e.what() << "\n";
            return false;
        }
        return true;
    };

    // ---- Initial load ----
    std::cerr << "Loading initial PanMAN/index...\n";
    auto t0 = std::chrono::steady_clock::now();
    try { loadAll(panmanPath, indexPath); }
    catch (const std::exception& e) { std::cerr << "Error: " << e.what() << "\n"; return 1; }
    auto tLoaded = std::chrono::steady_clock::now();
    std::cerr << "Loaded "
              << std::chrono::duration_cast<std::chrono::milliseconds>(tLoaded - t0).count()
              << " ms (" << leaves.size() << " leaves, k=" << idxK << " s=" << idxS << ")\n";

    // ---- Build MSA only if needed (reject-mode or write-output-panman) ----
    // Without these flags, the placement-driven loop doesn't need leaf states.
    bool needMsa = applyMoves && (rejectOnIncrease || vm["write-output-panman"].as<bool>());
    std::ofstream traceTsv;
    if (needMsa) {
        std::cerr << "Generating aligned MSA (needed for "
                  << (rejectOnIncrease ? "parsimony-reject" : "")
                  << ((rejectOnIncrease && vm["write-output-panman"].as<bool>()) ? " + " : "")
                  << (vm["write-output-panman"].as<bool>() ? "panman regeneration" : "")
                  << ")...\n";
        std::stringstream msaStream;
        try { tree->printFASTA(msaStream, true, false); }
        catch (const std::exception& e) {
            std::cerr << "Error: aligned FASTA: " << e.what() << "\n";
            return 1;
        }
        parseAlignedFasta(msaStream, msa, name2row);
        std::cerr << "MSA: " << msa.size() << " × "
                  << (msa.empty() ? 0 : (long)msa[0].size()) << "\n";

        // Constants needed for CI/RI from cache.totalParsimony() across passes.
        std::cerr << "Computing P_min, P_max constants from MSA...\n";
        auto pBase = computeParsimony(*tree, msa, name2row);
        P_min_const = pBase.P_min;
        P_max_const = pBase.P_max;
        n_var_const = pBase.n_variable;

        if (refreshIndex) {
            msaFastaPath = outPrefix + ".initial_msa.fa";
            std::cerr << "Writing initial MSA to " << msaFastaPath << "...\n";
            std::ofstream f(msaFastaPath);
            for (auto& kv : name2row) {
                f << ">" << kv.first << "\n" << msa[kv.second] << "\n";
            }
            std::cerr << "  Done.\n";
        }

        traceTsv.open(outPrefix + ".trace.tsv");
        traceTsv << "pass\tphase\tradius\tcum_applied\tn_variable\tP_fitch\tP_min\tP_max\tCI\tRI\n";

        // Baseline trace row
        double CI = pBase.P_fitch > 0 ? (double)pBase.P_min / pBase.P_fitch : 0;
        double RI = (pBase.P_max - pBase.P_min) > 0
                  ? (double)(pBase.P_max - pBase.P_fitch) / (pBase.P_max - pBase.P_min) : 1.0;
        traceTsv << "0\tbaseline\t0\t0\t" << pBase.n_variable << "\t" << pBase.P_fitch
                 << "\t" << pBase.P_min << "\t" << pBase.P_max << "\t"
                 << std::setprecision(8) << CI << "\t" << RI << "\n";
        traceTsv.flush();
        std::cerr << "Baseline: P_fitch=" << pBase.P_fitch
                  << " CI=" << std::setprecision(4) << CI << " RI=" << RI << "\n";
    } else if (applyMoves) {
        // No MSA needed; still open trace.tsv so we record per-pass move counts.
        traceTsv.open(outPrefix + ".trace.tsv");
        traceTsv << "pass\tphase\tradius\tcum_applied\tn_variable\tP_fitch\tP_min\tP_max\tCI\tRI\n";
    }

    std::string tmpFasta = outPrefix + ".tmp_leaf.fa";
    std::string tmpOut = outPrefix + ".tmp_placement.tsv";
    std::string r2Empty;

    std::ofstream tsv(outPrefix + ".moves.tsv");
    tsv << "pass\tphase\tleaf_id\tparent_id\tbest_node_id\tparent_score\tbest_score\tdelta\twould_move\tapplied\n";

    // Cumulative counters across all passes
    int total_better = 0, total_would_move = 0, total_skipped = 0;
    int total_applied = 0, total_failed = 0, total_rejected = 0;
    double max_delta = 0.0;
    auto tStart = std::chrono::steady_clock::now();

    int phase = 0;          // 0 = global, 1 = local
    int currentRadius = 0;
    int passN = 0;          // global pass counter (resets to 0 when phase switches)
    int globalPassN = 0;    // reported in trace
    int totalPassN = 0;     // for output tmp filenames (always increments)

    bool needCacheRebuild = true;  // set to true after refresh or on pass 0
    while (passN < maxPasses) {
        // (Re)build FitchCache only when needed: pass 0, or after a refresh.
        // The cache is incrementally maintained DURING a pass, so its end-state
        // matches the start-state of the next pass.
        if (applyMoves && rejectOnIncrease && needCacheRebuild) {
            std::cerr << "Pass " << passN << " (phase=" << (phase==0?"global":"local")
                      << ", radius=" << currentRadius << "): building Fitch cache...\n";
            auto tFc0 = std::chrono::steady_clock::now();
            fcache = std::make_unique<FitchCache>(msa, name2row, leaves.size() + 16);
            fcache->initialize(*tree);
            currentPFitch = fcache->totalParsimony();
            auto tFc1 = std::chrono::steady_clock::now();
            std::cerr << "  Cache built in "
                      << std::chrono::duration_cast<std::chrono::milliseconds>(tFc1 - tFc0).count()
                      << " ms. P_fitch=" << currentPFitch << "\n";
            needCacheRebuild = false;
        } else if (applyMoves && rejectOnIncrease) {
            std::cerr << "Pass " << passN << " (phase=" << (phase==0?"global":"local")
                      << ", radius=" << currentRadius << "): reusing cache from previous pass.\n";
        }

        int pass_better = 0, pass_would_move = 0, pass_skipped = 0;
        int pass_applied = 0, pass_failed = 0, pass_rejected = 0;
        double pass_max_delta = 0.0;
        const char* phaseStr = (phase == 0 ? "global" : "local");

        for (size_t i = 0; i < leaves.size(); i++) {
            auto* leaf = leaves[i];
            auto* parent = leaf->parent;
            const std::string& leafId = leaf->identifier;
            std::string parentId = parent ? parent->identifier : "";
            if (!parent) { pass_skipped++; continue; }

            auto itL = liteTree.allLiteNodes.find(leafId);
            if (itL == liteTree.allLiteNodes.end()) { pass_skipped++; continue; }
            uint32_t leafIdx = itL->second->nodeIndex;

            std::string seq;
            try { seq = tree->getStringFromReference(leafId, false, true); }
            catch (...) { pass_skipped++; continue; }
            if (seq.empty()) { pass_skipped++; continue; }
            writeFasta(tmpFasta, leafId, seq);

            placement::PlacementResult result;
            placement::TraversalParams params;
            params.k = idxK; params.s = idxS; params.l = idxL;
            params.t = idxT; params.open = idxOpen;
            params.skipNodeIndex = leafIdx;
            params.dedupReads = false; params.refineEnabled = false;
            try {
                if (placementVerbose) {
                    placement::placeLite(result, &liteTree, *reader, tmpFasta, r2Empty,
                                         tmpOut, params, nullptr);
                } else {
                    SilenceStdio s;
                    placement::placeLite(result, &liteTree, *reader, tmpFasta, r2Empty,
                                         tmpOut, params, nullptr);
                }
            } catch (const std::exception& e) {
                std::cerr << "  placement failed for " << leafId << ": " << e.what() << "\n";
                pass_skipped++;
                continue;
            }

            double parentScore = 0.0;
            auto itP = liteTree.allLiteNodes.find(parentId);
            if (itP != liteTree.allLiteNodes.end()) {
                parentScore = (double)itP->second->logContainmentScore;
            }
            double bestScore = result.bestLogContainmentScore;
            std::string bestId = result.bestLogContainmentNodeId;
            bool different = !bestId.empty() && (bestId != parentId);
            double delta = different ? (bestScore - parentScore) : 0.0;
            bool wouldMove = different && (delta > margin);

            // Local-radius constraint
            if (wouldMove && currentRadius > 0) {
                auto itT = tree->allNodes.find(bestId);
                if (itT == tree->allNodes.end() ||
                    bfsDistance(parent, itT->second, currentRadius) < 0) {
                    wouldMove = false;
                }
            }

            if (different) {
                pass_better++;
                if (delta > pass_max_delta) pass_max_delta = delta;
            }
            if (wouldMove) pass_would_move++;

            bool applied = false;
            if (applyMoves && wouldMove) {
                auto itT = tree->allNodes.find(bestId);
                if (itT == tree->allNodes.end()) {
                    pass_failed++;
                } else {
                    std::string err;
                    UndoState undoState;
                    bool ok = applySprMove(*tree, leaf, itT->second, internalCounter,
                                           undoState, err);
                    if (!ok) {
                        pass_failed++;
                    } else if (rejectOnIncrease) {
                        fcache->registerNode(undoState.newInternal);
                        std::vector<panmanUtils::Node*> applyBottoms{undoState.newInternal};
                        if (undoState.oldParentSuppressed) {
                            applyBottoms.push_back(undoState.leafOldParent);
                            applyBottoms.push_back(undoState.gp);
                        } else {
                            applyBottoms.push_back(undoState.leafOldParent);
                        }
                        long newP = fcache->updateAfterChange(applyBottoms);
                        if (newP > currentPFitch) {
                            fcache->forgetNode(undoState.newInternal);
                            undoSprMove(*tree, undoState);
                            std::vector<panmanUtils::Node*> undoBottoms{
                                undoState.targetParent, undoState.leafOldParent};
                            fcache->updateAfterChange(undoBottoms);
                            pass_rejected++;
                        } else {
                            currentPFitch = newP;
                            pass_applied++;
                            applied = true;
                        }
                    } else {
                        pass_applied++;
                        applied = true;
                    }
                }
            }

            tsv << passN << "\t" << phaseStr << "\t"
                << leafId << "\t" << parentId << "\t" << bestId << "\t"
                << std::fixed << std::setprecision(8)
                << parentScore << "\t" << bestScore << "\t" << delta << "\t"
                << (wouldMove ? "1" : "0") << "\t" << (applied ? "1" : "0") << "\n";

            if (progressEvery > 0 && ((i + 1) % progressEvery) == 0) {
                auto now = std::chrono::steady_clock::now();
                long ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - tStart).count();
                std::cerr << "  pass " << passN << " (" << phaseStr << ", R=" << currentRadius
                          << ") " << (i + 1) << "/" << leaves.size()
                          << "  applied=" << pass_applied
                          << " rejected=" << pass_rejected
                          << " failed=" << pass_failed;
                if (rejectOnIncrease) std::cerr << "  P=" << currentPFitch;
                std::cerr << "  (" << ms / 1000 << "s)\n";
            }
        }

        // Pass summary
        total_better += pass_better;
        total_would_move += pass_would_move;
        total_applied += pass_applied;
        total_failed += pass_failed;
        total_rejected += pass_rejected;
        total_skipped += pass_skipped;
        if (pass_max_delta > max_delta) max_delta = pass_max_delta;

        std::cerr << "\nPass " << passN << " (" << phaseStr << ", R=" << currentRadius
                  << ") done: applied=" << pass_applied << " rejected=" << pass_rejected
                  << " failed=" << pass_failed << " P_fitch=" << currentPFitch << "\n";

        // Per-pass trace row
        if (applyMoves && rejectOnIncrease && currentPFitch >= 0) {
            double CI = currentPFitch > 0 ? (double)P_min_const / currentPFitch : 0;
            double RI = (P_max_const - P_min_const) > 0
                      ? (double)(P_max_const - currentPFitch) / (P_max_const - P_min_const) : 1.0;
            traceTsv << (passN + 1) << "\t" << phaseStr << "\t" << currentRadius
                     << "\t" << total_applied << "\t" << n_var_const << "\t"
                     << currentPFitch << "\t" << P_min_const << "\t" << P_max_const << "\t"
                     << std::setprecision(8) << CI << "\t" << RI << "\n";
            traceTsv.flush();
        }

        totalPassN++;

        // Convergence: pass with 0 accepted moves means done with this phase.
        if (pass_applied == 0) {
            if (phase == 0 && localRadius > 0) {
                std::cerr << "Global converged after " << (passN + 1)
                          << " pass(es). Switching to local (radius=" << localRadius << ")\n";
                phase = 1;
                currentRadius = localRadius;
                passN = 0;  // reset per-phase counter
                continue;
            }
            std::cerr << "Phase " << phaseStr << " converged after " << (passN + 1) << " pass(es).\n";
            break;
        }

        passN++;
        if (passN >= maxPasses) {
            std::cerr << "Hit max-passes (" << maxPasses << ") in phase " << phaseStr << "\n";
            break;
        }

        // Refresh index for next pass
        if (refreshIndex) {
            std::cerr << "Refreshing index for pass " << totalPassN << "...\n";
            if (!refresh(totalPassN)) {
                std::cerr << "Refresh failed; stopping.\n";
                break;
            }
            needCacheRebuild = true;  // tree was reloaded; cache must rebuild
        }
    }  // end pass loop

    auto tEnd = std::chrono::steady_clock::now();
    long totalMs = std::chrono::duration_cast<std::chrono::milliseconds>(tEnd - tStart).count();

    std::cout << "\n=== panopt_spr ===\n";
    std::cout << "Total passes:           " << passN << "\n";
    std::cout << "Final phase:            " << (phase == 0 ? "global" : "local")
              << " (radius=" << currentRadius << ")\n";
    std::cout << "Cumulative better:      " << total_better << "\n";
    std::cout << "Cumulative would-move:  " << total_would_move << "\n";
    std::cout << "Cumulative applied:     " << total_applied << "\n";
    std::cout << "Cumulative failed:      " << total_failed << "\n";
    std::cout << "Cumulative rejected:    " << total_rejected << "\n";
    std::cout << "Max delta seen:         " << max_delta << "\n";
    std::cout << "Wall time:              " << totalMs << " ms\n";
    std::cout << "Wrote " << outPrefix << ".moves.tsv\n";

    if (applyMoves) {
        if (needMsa) {
            std::cerr << "Computing final parsimony...\n";
            auto pFinal = computeParsimony(*tree, msa, name2row);
            if (traceTsv.is_open()) {
                traceTsv << "final\tdone\t" << currentRadius << "\t" << total_applied
                         << "\t" << pFinal.n_variable << "\t" << pFinal.P_fitch
                         << "\t" << pFinal.P_min << "\t" << pFinal.P_max
                         << "\t" << std::setprecision(8) << pFinal.CI << "\t" << pFinal.RI << "\n";
            }
            std::cout << "Final P_fitch:          " << pFinal.P_fitch << "\n";
            std::cout << "Final CI:               " << std::setprecision(4) << pFinal.CI << "\n";
            std::cout << "Final RI:               " << pFinal.RI << "\n";
            std::cout << "Wrote " << outPrefix << ".trace.tsv\n";
        }

        std::ofstream nwk(outPrefix + ".optimized.nwk");
        nwk << tree->getNewickString(tree->root) << "\n";
        std::cout << "Wrote " << outPrefix << ".optimized.nwk\n";

        // Regenerate NucMuts and write a valid output PanMAN.
        if (vm["write-output-panman"].as<bool>()) {
            std::cerr << "Building column → (block,pos,gap) coordinate map...\n";
            long n_cols_check = 0;
            auto colMap = buildColMap(*tree, n_cols_check);
            std::cerr << "  " << n_cols_check << " columns mapped (MSA has "
                      << (msa.empty() ? 0 : msa[0].size()) << ").\n";
            if (n_cols_check != (long)(msa.empty() ? 0 : msa[0].size())) {
                std::cerr << "  Warning: colMap size != MSA cols. Output may be inconsistent.\n";
            }

            std::cerr << "Regenerating per-branch NucMuts via 5-state Fitch on (topology, MSA)...\n";
            regenerateBranchMutations(*tree, msa, name2row, colMap);

            std::string outPanman = outPrefix + ".optimized.panman";
            std::cerr << "Writing " << outPanman << " (lzma-compressed)...\n";
            try {
                std::ofstream rawOut(outPanman, std::ios::binary);
                boost::iostreams::filtering_streambuf<boost::iostreams::output> outBuf;
                boost::iostreams::lzma_params lz; lz.level = 9;
                outBuf.push(boost::iostreams::lzma_compressor(lz));
                outBuf.push(rawOut);
                std::ostream outStream(&outBuf);
                kj::std::StdOutputStream kjOut(outStream);
                tg->writeToFile(kjOut);  // TreeGroup format (matches panmanUtils' input expectation)
                boost::iostreams::close(outBuf);
            } catch (const std::exception& e) {
                std::cerr << "writeToFile failed: " << e.what() << "\n";
            }
            std::cout << "Wrote " << outPanman << "\n";
        }
    }

    std::remove(tmpFasta.c_str());
    std::remove(tmpOut.c_str());
    return 0;
}
