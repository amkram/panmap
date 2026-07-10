// MGSR unit tests. Only getDust (free function, mgsr.cpp:1541) is unit-testable;
// the EM/abundance path is behind private squareEM members, covered by the
// metagenomic e2e scenario (TESTING.md §4 / §7).
#include <boost/test/unit_test.hpp>

#include "mgsr.hpp"

#include <string>

BOOST_AUTO_TEST_SUITE(mgsr_tests)

BOOST_AUTO_TEST_CASE(dust_score_deterministic) {
    const std::string lowComplexity(200, 'A');
    std::string highComplexity;
    const char bases[] = {'A', 'C', 'G', 'T'};
    for (int i = 0; i < 200; ++i) highComplexity += bases[(i * 7 + (i / 4)) % 4];

    // Deterministic: same input -> same score.
    BOOST_TEST(mgsr::getDust(lowComplexity) == mgsr::getDust(lowComplexity));

    // Low-complexity sequence scores strictly higher than a varied one.
    double lo = mgsr::getDust(lowComplexity);
    double hi = mgsr::getDust(highComplexity);
    BOOST_TEST(lo > hi);

    // Edge cases must not crash.
    BOOST_CHECK_NO_THROW((void)mgsr::getDust(""));
    BOOST_CHECK_NO_THROW((void)mgsr::getDust("ACG"));
    BOOST_CHECK_NO_THROW((void)mgsr::getDust(std::string(50, 'N')));
}

BOOST_AUTO_TEST_SUITE_END()
