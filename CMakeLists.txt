cmake_minimum_required (VERSION 3.20)

project(panmap)

set(CMAKE_CXX_STANDARD 20)

cmake_policy(SET CMP0135 NEW)

include(FetchContent)

add_executable(panmap
    src/main.cpp
    src/place.cpp
    src/place.hpp
    src/genotype.cpp
    src/genotype.hpp
    src/pmi.cpp
    src/pmi.hpp
    src/tree.cpp
    src/tree.hpp
    src/util.hpp
    src/seed.hpp
    src/align.cpp
    src/align.hpp
)

# panmat
FetchContent_Declare(
    panmat
    GIT_REPOSITORY https://github.com/TurakhiaLab/pangenome-mat
    GIT_TAG 1a05d22f4e0922e0cf98f0589d0e7bf06c6de92f # testing-cmake branch
    GIT_SHALLOW TRUE
    )

FetchContent_MakeAvailable(panmat)

target_link_libraries(panmap PUBLIC samtools htslib panmat minimap2)
target_compile_options(panmap PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)

# Unit testing data
configure_file(src/test/data/test.aligned.fa test.aligned.fa COPYONLY)
configure_file(src/test/data/test.nwk test.nwk COPYONLY)
configure_file(src/test/data/test.json test.json COPYONLY)
configure_file(src/test/data/test.pmat test.pmat COPYONLY)
configure_file(src/test/data/test.fastq test.fastq COPYONLY)