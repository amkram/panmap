cmake_minimum_required(VERSION 3.20)

project(panmap VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(FetchContent)
set(FETCHCONTENT_QUIET OFF)

# Project options
option(OPTION_DEBUG "Compile with debug symbols and no optimizations." OFF)
option(OPTION_BUILD_TESTS "Build tests." OFF)
option(OPTION_BUILD_SIMULATE "Build simulate." ON)
option(OPTION_PROFILE "Enable profiling." ON)
set(TBB_DIR ${CMAKE_CURRENT_BINARY_DIR}/oneTBB-2019_U9)
set(VCPKG_DIR ${CMAKE_CURRENT_BINARY_DIR}/vcpkg)
set(CMAKE_TOOLCHAIN_FILE ${VCPKG_DIR}/scripts/buildsystems/vcpkg.cmake)
set(CMAKE_PREFIX_PATH ${TBB_DIR}/cmake)
# Fetch panman
FetchContent_Declare(
    panman
    GIT_REPOSITORY https://github.com/TurakhiaLab/panman
    GIT_TAG 3bc6f2c8cf29a8a106d58586c3bfb6c63dec6df1
    PATCH_COMMAND rm -rf ${TBB_DIR} ${VCPKG_DIR} && cd ${CMAKE_CURRENT_BINARY_DIR} && wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz && tar -xvzf 2019_U9.tar.gz && git clone https://github.com/microsoft/vcpkg.git ${VCPKG_DIR} && ${VCPKG_DIR}/bootstrap-vcpkg.sh && ${VCPKG_DIR}/vcpkg install jsoncpp
)

FetchContent_MakeAvailable(panman)


add_executable(panmap
    src/main.cpp
    src/conversion.cpp
    src/place.cpp
    src/genotype.cpp
    src/pmi.cpp
    src/tree.cpp
    src/mm_align.c
    src/pileup.c
)
if (OPTION_BUILD_TESTS)
add_executable(tests
   src/test/main.test.cpp
   src/test/pmi.test.cpp
   src/conversion.cpp
   src/place.cpp
   src/genotype.cpp
   src/pmi.cpp
   src/tree.cpp
   src/mm_align.c
   src/pileup.c
   )
endif()

target_include_directories(panmap PRIVATE  ${CMAKE_CURRENT_BINARY_DIR}/_deps/panman-build ${CMAKE_CURRENT_BINARY_DIR}/_deps/panman-src/src ${CMAKE_CURRENT_BINARY_DIR}/_deps/panman-src)


include(cmake/detect-arch.cmake)
detect_host_arch()
message("neon: ${HAS_NEON}")
message("Third party dependencies...")

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/libminimap2.a
    COMMAND make -C src/3rdparty/minimap2 $<$<BOOL:${HAS_NEON}>:aarch64=1>
    COMMAND mv src/3rdparty/minimap2/libminimap2.a ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND mv src/3rdparty/minimap2/minimap2 ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND rm src/3rdparty/minimap2/*.o
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
# FetchContent for bzip2
FetchContent_Declare(
    bzip2
    GIT_REPOSITORY https://github.com/libarchive/bzip2.git
    GIT_TAG master
    GIT_SHALLOW TRUE
)

FetchContent_Declare(
    bm
    GIT_REPOSITORY https://github.com/tlk00/BitMagic.git
    GIT_TAG master
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(bm)

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/libst.a ${CMAKE_CURRENT_BINARY_DIR}/libhts.a
    COMMAND cd src/3rdparty/samtools && autoheader && autoconf && ./configure --disable-lzma --disable-bz2 --disable-libcurl --enable-configure-htslib
    COMMAND make -C src/3rdparty/samtools
    COMMAND mv src/3rdparty/samtools/libst.a ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND mv src/3rdparty/samtools/htslib-1.20/libhts.a ${CMAKE_CURRENT_BINARY_DIR}
    #COMMAND rm -f src/3rdparty/samtools/*.o src/3rdparty/samtools/htslib-1.20/libhts.a
    #COMMAND rm -f src/3rdparty/samtools/*/*.o src/3rdparty/samtools/*/*/*.o
    BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/libhts.a
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_custom_target(build-deps DEPENDS
    libminimap2.a
    libst.a
    # libpanman.a
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

add_dependencies(panmap build-deps)
# Find dependencies
find_package(Boost COMPONENTS system iostreams filesystem unit_test_framework REQUIRED)
find_package(ZLIB REQUIRED)

target_include_directories(panmap PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/_deps/bm-src/src)
target_include_directories(panmap PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(panmap PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty)
target_include_directories(panmap PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/samtools/htslib-1.20)

if (OPTION_BUILD_TESTS)
target_include_directories(tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty)
target_include_directories(tests PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/_deps/bm-build)
target_include_directories(tests PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/samtools/htslib-1.20)
endif()

# Lets us do: 
#    #include <minimap2/kseq.h>
#    #include <htslib/sam.h>
#    ...


add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/capnp ${CMAKE_CURRENT_BINARY_DIR}/capnpc-c++ ${CMAKE_CURRENT_BINARY_DIR}/.libs/libcapnp.a
    COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/capnproto/c++ && autoupdate && autoheader && autoconf && ./configure
    COMMAND make -C ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/capnproto/c++
    COMMAND mv ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/capnproto/c++/capnp ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND mv ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/capnproto/c++/capnpc-c++ ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND cp -r ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/capnproto/c++/.libs/ ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Compiling capnproto"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_library(capnp STATIC IMPORTED)
set_target_properties(capnp PROPERTIES IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/.libs/libcapnp.a")
set_target_properties(capnp PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/capnproto/c++/src")

add_library(kj STATIC IMPORTED)
set_target_properties(kj PROPERTIES IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/.libs/libkj.a")
set_target_properties(kj PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/capnproto/c++/src")

set(CAPNP_SRCS ${CMAKE_CURRENT_BINARY_DIR}/index.capnp.c++)
set(CAPNP_HDRS ${CMAKE_CURRENT_BINARY_DIR}/index.capnp.h)


add_custom_command(
    OUTPUT ${CAPNP_SRCS} ${CAPNP_HDRS}
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/capnp compile -o ${CMAKE_CURRENT_BINARY_DIR}/capnpc-c++ ${CMAKE_CURRENT_SOURCE_DIR}/index.capnp --src-prefix ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND mv ${CMAKE_CURRENT_SOURCE_DIR}/index.capnp.c++ ${CMAKE_CURRENT_BINARY_DIR}/index.capnp.c++
    COMMAND mv ${CMAKE_CURRENT_SOURCE_DIR}/index.capnp.h ${CMAKE_CURRENT_BINARY_DIR}/index.capnp.h
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/index.capnp ${CMAKE_CURRENT_BINARY_DIR}/capnp ${CMAKE_CURRENT_BINARY_DIR}/capnpc-c++
    COMMENT "Compiling index.capnp"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

# Add the generated capnp sources to the target
target_sources(panmap PRIVATE ${CAPNP_SRCS})

target_include_directories(panmap PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

FetchContent_Declare(
  docopt
  GIT_REPOSITORY https://github.com/docopt/docopt.cpp
  GIT_TAG master
)

FetchContent_MakeAvailable(docopt)


# Ensure Boost iostreams, zlib, and protobuf are linked
target_link_libraries(panmap PUBLIC minimap2 samtools htslib docopt capnp kj Boost::iostreams Boost::filesystem Boost::system ZLIB::ZLIB panman)
target_link_options(panmap PUBLIC -lm -lpthread -lz -Wl,--allow-multiple-definition)


if (OPTION_BUILD_TESTS)
target_sources(tests PUBLIC ${CAPNP_SRCS})
target_include_directories(tests PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(tests minimap2 samtools htslib
                     ${Boost_FILESYSTEM_LIBRARY}
                     ${Boost_SYSTEM_LIBRARY}
                     ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
                     ${Boost_IOSTREAMS_LIBRARY}
                     ZLIB::ZLIB
                     docopt
                     capnp
                     kj
                     panman
                     )
target_compile_options(tests PUBLIC -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
target_link_options(tests PUBLIC -lm -lpthread -lz)
target_compile_options(tests PUBLIC -g -O3 -ggdb)
endif()

if (OPTION_DEBUG)
    target_compile_options(panmap PUBLIC -g -O1 -ggdb)
else()
    target_compile_options(panmap PUBLIC -g -O3)
endif()

# Unit testing data
configure_file(src/test/data/test.aligned.fa test.aligned.fa COPYONLY)
configure_file(src/test/data/test.nwk test.nwk COPYONLY)
configure_file(src/test/data/test.json test.json COPYONLY)
configure_file(src/test/data/test.pmat test.pmat COPYONLY)
configure_file(src/test/data/test.fastq test.fastq COPYONLY)
