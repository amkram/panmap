cmake_minimum_required (VERSION 3.26) 

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
    src/align.cpp
    src/align.hpp
)

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

# minimap2
include(ExternalProject)
FetchContent_Declare(
    minimap2-remote
    GIT_REPOSITORY https://github.com/nanoporetech/ont-minimap2/
    GIT_TAG 41b2c01bba3b1a89304f8220694385fcdbe9ed67
)
FetchContent_GetProperties(minimap2-remote)
if(NOT minimap2-remote_POPULATED)
  FetchContent_Populate(minimap2-remote)
endif()

set(minimap2_ROOT ${CMAKE_BINARY_DIR}/_deps/minimap2-remote-src)

ExternalProject_Add(minimap2-remote
    SOURCE_DIR ${minimap2_ROOT}
    INSTALL_COMMAND ""
    BUILD_COMMAND cd ${minimap2_ROOT} && cmake -S . -B cmake-build && cmake --build cmake-build --config Release -j
    # remove 'static' keyword from map.c functions
    PATCH_COMMAND cmake
        -Din_file:FILEPATH=${minimap2_ROOT}/src/3rdparty/minimap2/map.c
        -Dpatch_file:FILEPATH=${CMAKE_CURRENT_SOURCE_DIR}/cmake/map.c.patch
        -Dout_file:FILEPATH=${minimap2_ROOT}/src/3rdparty/minimap2/map.c
        -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/patch-file.cmake
)


add_library(minimap2 STATIC IMPORTED GLOBAL)

set_target_properties(minimap2 PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${minimap2_ROOT}/src/3rdparty/minimap2
    IMPORTED_LOCATION ${minimap2_ROOT}/cmake-build/libminimap2.a
)

add_dependencies(minimap2 minimap2-remote)

# samtools
set(samtools_ROOT ${CMAKE_CURRENT_BINARY_DIR}/samtools-remote-prefix/src/samtools-remote)
ExternalProject_Add(
    samtools-remote
    URL https://github.com/samtools/samtools/releases/download/1.19.2/samtools-1.19.2.tar.bz2
    URL_HASH SHA256=71f60499668e4c08e7d745fbff24c15cc8a0977abab1acd5d2bb419bdb065e96
    CONFIGURE_COMMAND ../samtools-remote/configure
    BUILD_COMMAND cd ../samtools-remote && make
    INSTALL_COMMAND ""
    PATCH_COMMAND cmake
    -Din_file:FILEPATH=${samtools_ROOT}/bam_plcmd.c
    -Dpatch_file:FILEPATH=${CMAKE_CURRENT_SOURCE_DIR}/cmake/bam_plcmd.c.patch
    -Dout_file:FILEPATH=${samtools_ROOT}/bam_plcmd.c
    -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/patch-file.cmake
)


add_library(samtools STATIC IMPORTED)
target_sources(samtools INTERFACE ${samtools_ROOT}/bam_plcmd.o ${samtools_ROOT}/sample.o)
set_property(TARGET samtools
    PROPERTY IMPORTED_LOCATION ${samtools_ROOT}/libst.a)
set_property(TARGET samtools
    PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${samtools_ROOT})

add_library(htslib STATIC IMPORTED)
set_property(TARGET htslib
    PROPERTY IMPORTED_LOCATION ${samtools_ROOT}/htslib-1.19.1/libhts.a)

# curl
FetchContent_Declare(curl
    URL https://curl.se/download/curl-8.5.0.tar.gz
)
FetchContent_MakeAvailable(curl)
find_package(curl REQUIRED)
target_link_libraries(minimap2 INTERFACE $<TARGET_NAME_IF_EXISTS:CURL::libcurl> $<TARGET_NAME_IF_EXISTS:CURL::libcurl_shared>)

target_link_libraries(panmap PUBLIC samtools htslib panmat minimap2)
target_compile_options(panmap PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
set_target_properties(panmap PROPERTIES LINK_FLAGS "-lbz2 -llzma")

# Unit testing data
configure_file(src/test/data/test.aligned.fa test.aligned.fa COPYONLY)
configure_file(src/test/data/test.nwk test.nwk COPYONLY)
configure_file(src/test/data/test.json test.json COPYONLY)
configure_file(src/test/data/test.pmat test.pmat COPYONLY)
configure_file(src/test/data/test.fastq test.fastq COPYONLY)