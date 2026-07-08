// Conversion/alignment-IO tests: exercise the real alignAndWriteBam (minimap2 +
// direct bam1_t construction) and verify the output BAM structurally by reading it
// back with htslib. The full align->bam->vcf->consensus path is covered by the e2e
// test; the stale samtools-generated conversion_test_data goldens are not revived.
#include <boost/test/unit_test.hpp>

#include "conversion.hpp"
#include "seeding.hpp"
#include "helpers/paths.hpp"
#include "helpers/seed_helpers.hpp"

extern "C" {
#include <htslib/sam.h>
}

#include <atomic>
#include <filesystem>
#include <string>
#include <unistd.h>
#include <vector>

namespace {
std::string tmpBam() {
    static std::atomic<uint64_t> ctr{0};
    auto name = "panmap_test_" + std::to_string(::getpid()) + "_" + std::to_string(ctr.fetch_add(1)) + ".bam";
    return (std::filesystem::temp_directory_path() / name).string();
}
}  // namespace

BOOST_AUTO_TEST_SUITE(conversion_tests)

BOOST_AUTO_TEST_CASE(align_and_write_bam_structure) {
    std::string reference = ts::readFasta(ts::dataPath("MZ515733.1.fa"));
    BOOST_REQUIRE(!reference.empty());

    // Real reads for the same genome.
    std::vector<std::string> readSequences, readQuals, readNames;
    seeding::readFastqPaired(readSequences, readQuals, readNames, ts::dataPath("MZ515733.1.fastq"), "");
    BOOST_REQUIRE(!readSequences.empty());

    std::string bamPath = tmpBam();
    alignAndWriteBam(readSequences,
                     readQuals,
                     readNames,
                     reference,
                     bamPath,
                     /*pairedEndReads=*/false,
                     /*n_threads=*/1);

    BOOST_REQUIRE(std::filesystem::exists(bamPath));
    BOOST_REQUIRE(std::filesystem::file_size(bamPath) > 0);

    // Read the BAM back with htslib and assert structure.
    samFile* in = sam_open(bamPath.c_str(), "r");
    BOOST_REQUIRE(in != nullptr);
    sam_hdr_t* hdr = sam_hdr_read(in);
    BOOST_REQUIRE(hdr != nullptr);

    BOOST_TEST(sam_hdr_nref(hdr) == 1);
    BOOST_TEST(static_cast<size_t>(sam_hdr_tid2len(hdr, 0)) == reference.size());

    bam1_t* rec = bam_init1();
    int total = 0, mapped = 0;
    int64_t lastPos = -1;
    bool sortedByPos = true;
    while (sam_read1(in, hdr, rec) >= 0) {
        total++;
        if (!(rec->core.flag & BAM_FUNMAP)) {
            mapped++;
            if (rec->core.pos < lastPos) sortedByPos = false;
            lastPos = rec->core.pos;
        }
    }
    bam_destroy1(rec);
    sam_hdr_destroy(hdr);
    sam_close(in);

    BOOST_TEST(total > 0);
    // Reads come from this genome, so the large majority must map.
    BOOST_TEST(mapped > total / 2);
    // Output is coordinate-sorted (mapped records appear in non-decreasing position).
    BOOST_TEST(sortedByPos);

    std::error_code ec;
    std::filesystem::remove(bamPath, ec);
    std::filesystem::remove(bamPath + ".bai", ec);
}

BOOST_AUTO_TEST_SUITE_END()
