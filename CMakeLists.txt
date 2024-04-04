cmake_minimum_required (VERSION 3.28) 

project(panmap)

set(CMAKE_CXX_STANDARD 20)

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
    src/align_from_seeds.c
    src/bam_to_mpileup.c
    src/kseq.h)

# panmat
FetchContent_Declare(
    panmat
    GIT_REPOSITORY https://github.com/TurakhiaLab/pangenome-mat
    GIT_TAG aef68be76d7a2b6f4cb26a69fca337f260906a9c # testing-cmake branch
    GIT_SHALLOW TRUE
    FIND_PACKAGE_ARGS CONFIG
    )

FetchContent_MakeAvailable(panmat)
include_directories(${protobuf_INCLUDE_DIRS})

include(ExternalProject)



# minimap2

cmake_host_system_information(RESULT has_sse2 QUERY HAS_SSE2)

if (APPLE AND CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "arm64")
    set(APPLE_ARM_FLAGS "aarch64=1 arm_neon=1")
endif()

ExternalProject_Add(
    minimap2-remote
    GIT_REPOSITORY https://github.com/lh3/minimap2
    GIT_TAG b6762503a9c9af34ae02f56784035ff852ed04a5
    GIT_SHALLOW TRUE
    CONFIGURE_COMMAND ""
    BUILD_COMMAND cd ../minimap2-remote && make ${APPLE_ARM_FLAGS} -j4
    INSTALL_COMMAND "")


add_library(minimap2 STATIC IMPORTED)
set_property(TARGET minimap2
    PROPERTY IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/minimap2-remote-prefix/src/minimap2-remote/libminimap2.a)

# samtools + htslib
ExternalProject_Add(
    samtools-remote
    URL https://github.com/samtools/samtools/releases/download/1.19.2/samtools-1.19.2.tar.bz2
    URL_HASH SHA256=71f60499668e4c08e7d745fbff24c15cc8a0977abab1acd5d2bb419bdb065e96
    CONFIGURE_COMMAND ../samtools-remote/configure
    BUILD_COMMAND cd ../samtools-remote && make -j4
    INSTALL_COMMAND "")

add_library(samtools STATIC IMPORTED)
set_property(TARGET samtools
    PROPERTY IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/samtools-remote-prefix/src/samtools-remote/libst.a)

target_link_libraries(panmap PUBLIC minimap2 samtools panmat z)
target_compile_options(panmap PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)

# Unit testing data
configure_file(src/test/data/test.aligned.fa test.aligned.fa COPYONLY)
configure_file(src/test/data/test.nwk test.nwk COPYONLY)
configure_file(src/test/data/test.json test.json COPYONLY)
configure_file(src/test/data/test.pmat test.pmat COPYONLY)
configure_file(src/test/data/test.fastq test.fastq COPYONLY)