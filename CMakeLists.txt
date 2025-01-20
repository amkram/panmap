cmake_minimum_required(VERSION 3.12)

project(panmap VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 20)


include(FetchContent)
set(FETCHCONTENT_QUIET OFF)
# set(BUILD_SHARED_LIBS ON)
set(POSITION_INDEPENDENT_CODE ON)
set(ABSL_PROPAGATE_CXX_STD ON)


# Project options
option(OPTION_DEBUG "Compile with debug symbols and no optimizations." OFF)
option(OPTION_BUILD_TESTS "Build tests." OFF)
option(OPTION_BUILD_SIMULATE "Build simulate." ON)
option(OPTION_PROFILE "Enable profiling." ON)

set(TBB_DIR ${CMAKE_CURRENT_BINARY_DIR}/oneTBB-2019_U9)
# set(VCPKG_DIR ${CMAKE_CURRENT_BINARY_DIR}/vcpkg)
# set(CMAKE_TOOLCHAIN_FILE ${VCPKG_DIR}/scripts/buildsystems/vcpkg.cmake)

# Fetch panman
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
FetchContent_Declare(
    panman
    GIT_REPOSITORY https://github.com/TurakhiaLab/panman
    GIT_TAG c7ca48da772ae8038dfd6e0bc6cc410319c9b9ab
    PATCH_COMMAND rm -rf ${TBB_DIR} && cd ${CMAKE_CURRENT_BINARY_DIR} && wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz && tar -xvzf 2019_U9.tar.gz
)
FetchContent_MakeAvailable(panman)

include_directories(${TBB_DIR}/include)

add_executable(panmap
    src/main.cpp
    src/conversion.cpp
    src/place.cpp
    src/genotype.cpp
    src/pmi.cpp
    src/seed_annotated_tree.cpp
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
   src/seed_annotated_tree.cpp
   src/mm_align.c
   src/pileup.c
   )
endif()

target_include_directories(panmap PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/_deps/panman-build ${CMAKE_CURRENT_BINARY_DIR}/_deps/panman-src/src ${CMAKE_CURRENT_BINARY_DIR}/_deps/panman-src)


include(cmake/detect-arch.cmake)
detect_host_arch()
message("neon: ${HAS_NEON}")
message("Third party dependencies...")

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/lib/libminimap2.a
    COMMAND cd src/3rdparty/minimap2 && make clean
    COMMAND make -C src/3rdparty/minimap2 $<$<BOOL:${HAS_NEON}>:aarch64=1>
    COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/lib ${CMAKE_CURRENT_BINARY_DIR}/bin
    COMMAND cp src/3rdparty/minimap2/libminimap2.a ${CMAKE_CURRENT_BINARY_DIR}/lib/libminimap2.a
    COMMAND cp src/3rdparty/minimap2/minimap2 ${CMAKE_CURRENT_BINARY_DIR}/bin/minimap2
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
FetchContent_Declare(deflate
    GIT_REPOSITORY https://github.com/ebiggers/libdeflate
    GIT_SHALLOW TRUE
    )
FetchContent_MakeAvailable(deflate)
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/lib/libst.a ${CMAKE_CURRENT_BINARY_DIR}/lib/libhts.a
    COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/samtools/htslib-1.20 && autoreconf -i -Wno-syntax && ./configure --prefix=${CMAKE_CURRENT_BINARY_DIR} --disable-lzma --disable-bz2 --disable-libcurl && make -j && make install
    COMMAND cd src/3rdparty/samtools && autoheader && autoconf -Wno-syntax && ./configure --with-htslib=${CMAKE_CURRENT_BINARY_DIR} --prefix=${CMAKE_CURRENT_BINARY_DIR} --disable-lzma --disable-bz2 --disable-libcurl --disable-shared --with-pic && make -j && cp libst.a ${CMAKE_CURRENT_BINARY_DIR}/lib/libst.a
    BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/lib/libhts.a ${CMAKE_CURRENT_BINARY_DIR}/lib/libst.a
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/bin/capnp ${CMAKE_CURRENT_BINARY_DIR}/bin/capnpc-c++ ${CMAKE_CURRENT_BINARY_DIR}/lib/libcapnp.a ${CMAKE_CURRENT_BINARY_DIR}/lib/libkj.a ${CMAKE_CURRENT_BINARY_DIR}/include/capnp
    COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/capnproto/c++ && autoupdate && autoheader && autoconf && ./configure --prefix=${CMAKE_CURRENT_BINARY_DIR} --exec-prefix=${CMAKE_CURRENT_BINARY_DIR} && make clean && make -j && make install
    COMMENT "Compiling capnproto"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

set(BWA_SOURCES
    src/3rdparty/bwa/bwa.c
    src/3rdparty/bwa/run.c
    src/3rdparty/bwa/bwase.c
    src/3rdparty/bwa/bwaseqio.c
    src/3rdparty/bwa/bwt.c
    src/3rdparty/bwa/bwtaln.c
    src/3rdparty/bwa/bwtindex.c
    src/3rdparty/bwa/bwt_gen.c
    src/3rdparty/bwa/bwtsw2_core.c
    src/3rdparty/bwa/bwtsw2_main.c
    src/3rdparty/bwa/bwtsw2_aux.c
    src/3rdparty/bwa/bwt_lite.c
    src/3rdparty/bwa/bwtsw2_chain.c
    src/3rdparty/bwa/fastmap.c
    src/3rdparty/bwa/bwtsw2_pair.c
    # src/3rdparty/bwa/main.c
    src/3rdparty/bwa/utils.c
    src/3rdparty/bwa/bwape.c
    src/3rdparty/bwa/kopen.c
    src/3rdparty/bwa/pemerge.c
    src/3rdparty/bwa/maxk.c
    src/3rdparty/bwa/bwashm.c
    src/3rdparty/bwa/bntseq.c
    src/3rdparty/bwa/is.c
    src/3rdparty/bwa/bwamem.c
    src/3rdparty/bwa/bwamem_pair.c
    src/3rdparty/bwa/bwamem_extra.c
    src/3rdparty/bwa/malloc_wrap.c
    src/3rdparty/bwa/QSufSort.c
    src/3rdparty/bwa/bwt_gen.c
    src/3rdparty/bwa/rope.c
    src/3rdparty/bwa/utils.c
    src/3rdparty/bwa/kstring.c
    src/3rdparty/bwa/ksw.c
    src/3rdparty/bwa/bwt.c
    src/3rdparty/bwa/bamlite.c
    src/3rdparty/bwa/bntseq.c
    src/3rdparty/bwa/bwashm.c
    src/3rdparty/bwa/bwase.c
    src/3rdparty/bwa/bwaseqio.c
    src/3rdparty/bwa/bwt_lite.c
    src/3rdparty/bwa/bwtaln.c
    src/3rdparty/bwa/bwtindex.c
    src/3rdparty/bwa/is.c
    src/3rdparty/bwa/kopen.c
    src/3rdparty/bwa/kstring.c
    src/3rdparty/bwa/ksw.c
    src/3rdparty/bwa/utils.c
    src/3rdparty/bwa/bwt_gen.c
    src/3rdparty/bwa/rope.c
    src/3rdparty/bwa/rle.c
    src/3rdparty/bwa/bwtsw2_aux.c
    src/3rdparty/bwa/bwtsw2_chain.c
    src/3rdparty/bwa/bwtsw2_main.c
    src/3rdparty/bwa/bwtsw2_pair.c
    src/3rdparty/bwa/bwt_gen.c
    src/3rdparty/bwa/bwtgap.c
    src/3rdparty/bwa/bwtaln.c
    src/3rdparty/bwa/bwtindex.c
    src/3rdparty/bwa/is.c
    src/3rdparty/bwa/bwase.c
    src/3rdparty/bwa/bwaseqio.c
    src/3rdparty/bwa/bwashm.c
    src/3rdparty/bwa/bwt_lite.c
    src/3rdparty/bwa/bwtsw2_aux.c
    src/3rdparty/bwa/bwtsw2_chain.c
    src/3rdparty/bwa/bwtsw2_main.c
    src/3rdparty/bwa/bwtsw2_pair.c
)

add_library(bwa STATIC ${BWA_SOURCES})
target_include_directories(bwa PUBLIC src/3rdparty/bwa)

file(GLOB BCFTOOLS_SOURCES src/3rdparty/bcftools/*.c)
add_library(bcftools STATIC ${BCFTOOLS_SOURCES})

target_include_directories(bcftools PUBLIC src/3rdparty/bcftools)


add_custom_target(build-deps DEPENDS
    ${CMAKE_CURRENT_BINARY_DIR}/lib/libst.a
    ${CMAKE_CURRENT_BINARY_DIR}/lib/libhts.a
    ${CMAKE_CURRENT_BINARY_DIR}/lib/libminimap2.a
    ${CMAKE_CURRENT_BINARY_DIR}/lib/libcapnp.a
    ${CMAKE_CURRENT_BINARY_DIR}/lib/libkj.a
    ${CMAKE_CURRENT_BINARY_DIR}/bin/capnp
    ${CMAKE_CURRENT_BINARY_DIR}/include/capnp
)
add_library(kj STATIC IMPORTED)
set_target_properties(kj PROPERTIES IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/lib/libkj.a")

add_library(capnp STATIC IMPORTED)
set_target_properties(capnp PROPERTIES IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/lib/libcapnp.a")

add_dependencies(kj build-deps)
add_dependencies(capnp kj)

add_library(minimap2 STATIC IMPORTED)
set_target_properties(minimap2 PROPERTIES IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/lib/libminimap2.a")
add_dependencies(minimap2 build-deps)

add_library(samtools STATIC IMPORTED)
set_target_properties(samtools PROPERTIES IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/lib/libst.a")
add_dependencies(samtools build-deps)

add_library(htslib STATIC IMPORTED)
set_target_properties(htslib PROPERTIES IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/lib/libhts.a")
add_dependencies(htslib build-deps)

add_dependencies(panmap build-deps)



# Find dependencies
find_package(Boost COMPONENTS system iostreams filesystem unit_test_framework REQUIRED)
find_package(ZLIB REQUIRED)
target_include_directories(panmap PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(panmap PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty)
target_include_directories(panmap PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/samtools/htslib-1.20)

if (OPTION_BUILD_TESTS)
target_include_directories(tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty)
target_include_directories(tests PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/samtools/htslib-1.20)
endif()

if(OPTION_BUILD_SIMULATE)
add_executable(simulate
    src/test/simulate.cpp
    src/conversion.cpp
    src/place.cpp
    src/genotype.cpp
    src/pmi.cpp
    src/seed_annotated_tree.cpp
    src/mm_align.c
    src/pileup.c
)
target_link_libraries(simulate PRIVATE bwa bcftools minimap2 samtools htslib docopt capnp kj Boost::iostreams Boost::filesystem Boost::system ZLIB::ZLIB panman tbb tbbmalloc tbbmalloc_proxy deflate m)
target_compile_options(simulate PUBLIC -g -O3)
target_include_directories(simulate PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(simulate PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty)
target_include_directories(simulate PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/samtools/htslib-1.20)
target_include_directories(simulate PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/_deps/panman-build ${CMAKE_CURRENT_BINARY_DIR}/_deps/panman-src/src ${CMAKE_CURRENT_BINARY_DIR}/_deps/panman-src)
target_include_directories(simulate PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(simulate PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/include)

endif()

set(CAPNP_SRCS ${CMAKE_CURRENT_BINARY_DIR}/index.capnp.c++)
set(CAPNP_HDRS ${CMAKE_CURRENT_BINARY_DIR}/index.capnp.h)

add_custom_command(
    OUTPUT ${CAPNP_SRCS} ${CAPNP_HDRS}
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/bin/capnp compile -o ${CMAKE_CURRENT_BINARY_DIR}/bin/capnpc-c++ ${CMAKE_CURRENT_SOURCE_DIR}/index.capnp --src-prefix ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND mv ${CMAKE_CURRENT_SOURCE_DIR}/index.capnp.c++ ${CMAKE_CURRENT_BINARY_DIR}/index.capnp.c++
    COMMAND mv ${CMAKE_CURRENT_SOURCE_DIR}/index.capnp.h ${CMAKE_CURRENT_BINARY_DIR}/index.capnp.h
    COMMENT "Compiling index.capnp"
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/bin/capnp
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

# Add the generated capnp sources to the target
target_sources(panmap PRIVATE ${CAPNP_SRCS})

target_include_directories(panmap PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(panmap PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/include)

FetchContent_Declare(
  docopt
  GIT_REPOSITORY https://github.com/docopt/docopt.cpp
  GIT_TAG master
)

FetchContent_MakeAvailable(docopt)

target_compile_options(bwa PRIVATE -fPIC)
target_compile_options(bcftools PRIVATE -fPIC)

target_link_libraries(panmap PRIVATE bwa bcftools minimap2 samtools htslib docopt capnp kj Boost::iostreams Boost::filesystem Boost::system ZLIB::ZLIB panman tbb tbbmalloc tbbmalloc_proxy deflate m)
target_link_options(panmap PRIVATE -lm -lpthread)
target_compile_options(panmap PRIVATE -Wall)
target_compile_features(panmap PRIVATE cxx_std_20)


if (OPTION_BUILD_TESTS)
target_sources(tests PUBLIC ${CAPNP_SRCS})
target_include_directories(tests PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(tests bwa bcftools minimap2 samtools htslib
                     ${Boost_FILESYSTEM_LIBRARY}
                     ${Boost_SYSTEM_LIBRARY}
                     ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
                     ${Boost_IOSTREAMS_LIBRARY}
                     ZLIB::ZLIB
                     docopt
                     capnp
                     kj
                     panman
                     tbb_static tbbmalloc_static tbbmalloc_proxy_static 
                     deflate
                     m
                     )
target_compile_options(tests PUBLIC -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
target_link_options(tests PUBLIC -lm -lpthread )
target_compile_options(tests PUBLIC -g -O3 -ggdb)
endif()


if (OPTION_DEBUG)
    target_compile_options(panmap PUBLIC -g -O3 -ggdb)
else()
    target_compile_options(panmap PUBLIC -g -O3)
endif()

# Unit testing data
configure_file(src/test/data/test.aligned.fa test.aligned.fa COPYONLY)
configure_file(src/test/data/test.nwk test.nwk COPYONLY)
configure_file(src/test/data/test.json test.json COPYONLY)
configure_file(src/test/data/test.pmat test.pmat COPYONLY)
configure_file(src/test/data/test.fastq test.fastq COPYONLY)


install(TARGETS panmap RUNTIME DESTINATION /usr/local/bin)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib/ DESTINATION /usr/local/lib FILES_MATCHING PATTERN "*.a" PATTERN "*.so" PATTERN "*.dylib" PATTERN "*.la" PATTERN "*.so.*")
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin/ DESTINATION /usr/local/bin FILES_MATCHING PATTERN "*")
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/ DESTINATION /usr/local/include)
