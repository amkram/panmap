#include <boost/test/unit_test.hpp>

#include "seeding.hpp"

#include <random>
#include <string>

namespace {
std::string randomDna(std::mt19937_64& gen, size_t len) {
    static const char bases[] = {'A', 'C', 'G', 'T'};
    std::uniform_int_distribution<int> d(0, 3);
    std::string s(len, 'A');
    for (auto& c : s) c = bases[d(gen)];
    return s;
}
}  // namespace

BOOST_AUTO_TEST_SUITE(seeding_tests)

BOOST_AUTO_TEST_CASE(hashseq_determinism_and_canonical) {
    std::mt19937_64 gen(1234);
    for (int i = 0; i < 50; ++i) {
        std::string s = randomDna(gen, 21);

        // Determinism: same input -> same hashes.
        auto a = seeding::hashSeq(s);
        auto b = seeding::hashSeq(s);
        BOOST_TEST((a.first == b.first && a.second == b.second));

        // Canonical hash min(fwd,rev) matches for a sequence and its reverse
        // complement (orientation-invariant).
        auto rc = seeding::hashSeq(seeding::reverseComplement(s));
        BOOST_TEST(std::min(a.first, a.second) == std::min(rc.first, rc.second));
    }
}

BOOST_AUTO_TEST_CASE(rollingsyncmers_contract) {
    std::mt19937_64 gen(99);
    for (int rep = 0; rep < 10; ++rep) {
        std::string seq = randomDna(gen, 200);
        for (int k : {15, 19, 31}) {
            for (int s : {6, 8}) {
                auto all = seeding::rollingSyncmers(seq, k, s, /*open=*/false, /*t=*/0, /*returnAll=*/true);
                // One entry per k-mer window.
                BOOST_TEST(all.size() == seq.size() - k + 1);
                // Determinism.
                auto all2 = seeding::rollingSyncmers(seq, k, s, false, 0, true);
                BOOST_REQUIRE(all.size() == all2.size());
                for (size_t i = 0; i < all.size(); ++i) {
                    BOOST_TEST(std::get<0>(all[i]) == std::get<0>(all2[i]));
                }
                // Non-sentinel entries store the canonical k-mer hash; returnAll=true yields
                // UINT64_MAX sentinels at non-syncmer positions.
                for (const auto& [hash, isRev, isSync, pos] : all) {
                    BOOST_REQUIRE(pos >= 0);
                    BOOST_REQUIRE(static_cast<size_t>(pos) + k <= seq.size());
                    if (hash != UINT64_MAX) {
                        auto h = seeding::hashSeq(seq.substr(pos, k));
                        BOOST_TEST(hash == std::min(h.first, h.second));
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(rollingsyncmers_syncmer_subset_of_all) {
    std::mt19937_64 gen(7);
    std::string seq = randomDna(gen, 300);
    const int k = 19, s = 8;
    auto all = seeding::rollingSyncmers(seq, k, s, false, 0, /*returnAll=*/true);
    auto only = seeding::rollingSyncmers(seq, k, s, false, 0, /*returnAll=*/false);

    size_t flagged = 0;
    for (const auto& e : all)
        if (std::get<2>(e)) flagged++;
    // The syncmer-only call returns exactly the flagged subset.
    BOOST_TEST(only.size() == flagged);
    BOOST_TEST(flagged > 0);
    BOOST_TEST(flagged <= all.size());
    for (const auto& e : only) BOOST_TEST(std::get<2>(e) == true);
}

BOOST_AUTO_TEST_CASE(reversecomplement_and_comp) {
    // Involution over random sequences.
    std::mt19937_64 gen(2024);
    for (int i = 0; i < 100; ++i) {
        std::string s = randomDna(gen, 30);
        BOOST_TEST(seeding::reverseComplement(seeding::reverseComplement(s)) == s);
    }
    BOOST_TEST(seeding::reverseComplement("ACGT") == "ACGT");  // palindrome
    BOOST_TEST(seeding::reverseComplement("AAAA") == "TTTT");
    BOOST_TEST(seeding::reverseComplement("GCGC") == "GCGC");  // palindrome
    BOOST_TEST(seeding::reverseComplement("A") == "T");
    // comp() base complement table: involution, case preserved.
    BOOST_TEST(seeding::comp(seeding::comp('A')) == 'A');
    BOOST_TEST(seeding::comp(seeding::comp('c')) == 'c');
    BOOST_TEST(seeding::comp('A') == 'T');
    BOOST_TEST(seeding::comp('c') == 'g');
    BOOST_TEST(seeding::comp('N') == 'N');
}

BOOST_AUTO_TEST_CASE(hpc_compress_and_mapping) {
    BOOST_TEST(seeding::hpcCompress("") == "");
    BOOST_TEST(seeding::hpcCompress("AAAA") == "A");
    BOOST_TEST(seeding::hpcCompress("ACGT") == "ACGT");
    BOOST_TEST(seeding::hpcCompress("AAACCCGGG") == "ACG");

    auto [compressed, mapping] = seeding::hpcCompressWithMapping("AAACCCGGGT");
    BOOST_TEST(compressed == "ACGT");
    BOOST_REQUIRE(mapping.size() == compressed.size());
    std::string orig = "AAACCCGGGT";
    for (size_t i = 0; i < compressed.size(); ++i) {
        BOOST_REQUIRE(mapping[i] < orig.size());
        BOOST_TEST(static_cast<char>(::toupper(orig[mapping[i]])) == compressed[i]);
        if (i > 0) BOOST_TEST(mapping[i] > mapping[i - 1]);
    }
}

BOOST_AUTO_TEST_SUITE_END()
