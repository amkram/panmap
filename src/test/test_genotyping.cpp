// Genotyping unit tests: mutation-matrix (.mm) parsing.
//
// NOTE: the orphaned src/test/data/genotype_test_data/*.mm fixtures were in an OLD
// format (plain numbers on the indel lines); the current parser (genotyping.cpp:63-92)
// requires "size:prob" pairs. Rather than commit stale fixtures, these tests write a
// valid .mm in the CURRENT format inline and assert the real parser's behavior.
#include <boost/test/unit_test.hpp>

#include "genotyping.hpp"

#include <atomic>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unistd.h>

namespace {
std::string writeTempMm(const std::string& contents) {
    static std::atomic<uint64_t> ctr{0};
    auto name = "panmap_test_" + std::to_string(::getpid()) + "_" + std::to_string(ctr.fetch_add(1)) + ".mm";
    std::string path = (std::filesystem::temp_directory_path() / name).string();
    std::ofstream out(path);
    out << contents;
    out.close();
    return path;
}

genotyping::mutationMatrices parseMm(const std::string& contents) {
    std::string path = writeTempMm(contents);
    std::ifstream inf(path);
    genotyping::mutationMatrices mm;
    try {
        genotyping::fillMutationMatricesFromFile(mm, inf);
    } catch (...) {
        std::error_code ec;
        std::filesystem::remove(path, ec);
        throw;
    }
    std::error_code ec;
    std::filesystem::remove(path, ec);
    return mm;
}

// 4x4 substitution matrix + ins/del lines in size:prob format.
const std::string kValidMm =
    "1 20 25 34\n"
    "21 1 24 12\n"
    "20 22 1 23\n"
    "20 21 19 1\n"
    "1:0.05 2:0.40 3:0.50\n"
    "1:0.05 2:0.45\n";
}  // namespace

BOOST_AUTO_TEST_SUITE(genotyping_tests)

BOOST_AUTO_TEST_CASE(mm_parse_valid) {
    auto mm = parseMm(kValidMm);
    BOOST_TEST(mm.filled == true);

    BOOST_REQUIRE(mm.submat.size() == 4);
    BOOST_REQUIRE(mm.submat[0].size() == 4);
    BOOST_TEST(mm.submat[0][0] == 1.0);
    BOOST_TEST(mm.submat[0][3] == 34.0);
    BOOST_TEST(mm.submat[3][0] == 20.0);

    BOOST_TEST(mm.insmat.size() == 3u);
    BOOST_TEST(mm.delmat.size() == 2u);
    BOOST_TEST(mm.insmat.at(1) == 0.05);
    BOOST_TEST(mm.insmat.at(3) == 0.50);
    BOOST_TEST(mm.delmat.at(2) == 0.45);

    // maxInsLogProb is the largest insertion log-prob seen.
    BOOST_TEST(mm.maxInsLogProb == 0.50);
    BOOST_TEST(mm.maxDelLogProb == 0.45);
}

BOOST_AUTO_TEST_CASE(mm_parse_invalid_throws) {
    // Truncated: only the 4 submat rows, no indel lines (idx != 6).
    BOOST_CHECK_THROW(parseMm("1 20 25 34\n21 1 24 12\n20 22 1 23\n20 21 19 1\n"),
                      std::invalid_argument);

    // Submat row with the wrong field count.
    BOOST_CHECK_THROW(parseMm("1 20 25\n21 1 24 12\n20 22 1 23\n20 21 19 1\n1:0.05\n1:0.05\n"),
                      std::invalid_argument);

    // Indel line not in size:prob format (no colon).
    BOOST_CHECK_THROW(parseMm("1 20 25 34\n21 1 24 12\n20 22 1 23\n20 21 19 1\n40 50\n1:0.05\n"),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
