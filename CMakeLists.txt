cmake_minimum_required (VERSION 3.20)

project(panmap)

set(CMAKE_CXX_STANDARD 20)

# Project options
option(OPTION_DEBUG "Compile with debug symbols and no optimizations." ON)
option(OPTION_BUILD_TESTS "Build tests." ON)
option(OPTION_BUILD_SIMULATE "Build simulate." ON)

cmake_policy(SET CMP0135 NEW)

include(FetchContent)

add_executable(panmap
    src/main.cpp
    src/place.cpp
    src/genotype.cpp
    src/pmi.cpp
    src/tree.cpp
    src/conversion.cpp
    src/mm_align.c
    src/pileup.c
)

# system
include(cmake/detect-arch.cmake)
detect_host_arch()
message("neon: ${HAS_NEON}")
message("Third party dependencies...")


add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/libminimap2.a
    COMMAND make -C src/3rdparty/minimap2 $<$<BOOL:${HAS_NEON}>:aarch64=1>
    COMMAND mv src/3rdparty/minimap2/libminimap2.a ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND rm src/3rdparty/minimap2/*.o
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/libhts.a
    COMMAND src/3rdparty/samtools/htslib-1.20/configure --disable-bz2
        --disable-largefile --disable-lzma --disable-libcurl
    COMMAND make -C src/3rdparty/samtools/htslib-1.20 lib-static
    COMMAND mv src/3rdparty/samtools/htslib-1.20/libhts.a ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND rm src/3rdparty/samtools/htslib-1.20/*.o
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/libst.a
    COMMAND src/3rdparty/samtools/configure
    COMMAND make -C src/3rdparty/samtools libst.a # specific .o files here
    COMMAND mv src/3rdparty/samtools/libst.a ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND rm -f src/3rdparty/samtools/*.o src/3rdparty/samtools/htslib-1.20/libhts.a
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/libhts.a
)


add_custom_target(build-deps DEPENDS
    libminimap2.a
    libst.a
    libhts.a
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_library(minimap2 STATIC IMPORTED)
set_target_properties(minimap2 PROPERTIES IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/libminimap2.a")
add_dependencies(minimap2 build-deps)

add_library(samtools STATIC IMPORTED)
set_target_properties(samtools PROPERTIES IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/libst.a")
add_dependencies(samtools build-deps)

add_library(htslib STATIC IMPORTED)
set_target_properties(htslib PROPERTIES IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/libhts.a")
add_dependencies(htslib build-deps)

# panmat
FetchContent_Declare(
    panmat
    GIT_REPOSITORY https://github.com/TurakhiaLab/pangenome-mat
    GIT_TAG 1a05d22f4e0922e0cf98f0589d0e7bf06c6de92f # testing-cmake branch
    GIT_SHALLOW TRUE
    )

FetchContent_MakeAvailable(panmat)

target_include_directories(panmap PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/_deps/panmat-build)
target_include_directories(panmap PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty)
target_include_directories(panmap PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/samtools/htslib-1.20)
# Lets us do: 
#    #include <minimap2/kseq.h>
#    #include <htslib/sam.h>    
#    ...

target_link_libraries(panmap PUBLIC panmat minimap2 samtools htslib)
target_compile_options(panmap PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
target_link_options(panmap PRIVATE -lm -lz -lpthread)

# Unit testing data
configure_file(src/test/data/test.aligned.fa test.aligned.fa COPYONLY)
configure_file(src/test/data/test.nwk test.nwk COPYONLY)
configure_file(src/test/data/test.json test.json COPYONLY)
configure_file(src/test/data/test.pmat test.pmat COPYONLY)
configure_file(src/test/data/test.fastq test.fastq COPYONLY)
