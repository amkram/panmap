# Testing

Tests are off by default. Turn them on with `-DOPTION_BUILD_TESTS=ON`. CTest runs three
tests: the Boost.Test binary (`panmap_tests`, the `unit` test), the end-to-end script
(`e2e`), and the README-demo regression check (`examples`, which diffs CLI output against
`examples/expected/`).

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DOPTION_BUILD_TESTS=ON
cmake --build build -j
ctest --test-dir build --output-on-failure        # unit + e2e

./build/bin/panmap_tests --run_test=placement_tests  # one suite
ctest --test-dir build -R e2e                         # just the script
```

## Layout

```
src/test/
  helpers/              shared helpers (testsupport lib)
    test_index.*        TestIndex / loadIndex / IndexData
    seed_helpers.*      extractSeeds, readFasta, generateReads
    tree_helpers.*      findNodeIndex, pathToRoot, RSVPanmanFixture
    traversal.*         reconstructGenomeSeeds
    metrics_oracle.*    GroundTruthMetrics
    paths.hpp
  test_seeding.cpp      hashSeq, syncmers, reverseComplement, HPC
  test_core_utils.cpp   zstd round-trip, LiteTree invariants
  test_index.cpp        index round-trip; delta reconstruction vs direct extraction
  test_placement.cpp    the 5 metrics vs the ground-truth oracle
  test_mgsr.cpp         getDust
  test_genotyping.cpp   .mm parsing
  test_conversion.cpp   alignAndWriteBam (BAM read back and checked)
  test_main.cpp         Boost.Test entry point
  e2e/run_e2e.sh        full CLI pipelines: index/place/genotype/consensus,
                        variant-call correctness, paired-end, --meta mixture
  data/                 fixtures
examples/
  check_examples.sh     the three README demos diffed against examples/expected/
```