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
option(OPTION_ENABLE_FIXED_KMER "Enable optimized k-mer size implementations (8, 16, 32)." ON)

FetchContent_Declare(ftxui
    GIT_REPOSITORY https://github.com/ArthurSonzogni/ftxui
    GIT_TAG cdf28903a7781f97ba94d30b79c3a4b0c97ccce7
)
FetchContent_MakeAvailable(ftxui)

# Fetch and configure spdlog
FetchContent_Declare(
    spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog.git
    GIT_TAG v1.12.0
)
FetchContent_MakeAvailable(spdlog)

# Get spdlog properties and ensure it's available
FetchContent_GetProperties(spdlog)
if(NOT spdlog_POPULATED)
    FetchContent_Populate(spdlog)
    add_subdirectory(${spdlog_SOURCE_DIR} ${spdlog_BINARY_DIR})
endif()

# Add include directories for spdlog globally
include_directories(${spdlog_SOURCE_DIR}/include)

# Add profiling flags if enabled
if(OPTION_PROFILE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
endif()



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
    src/genotyping.cpp
    src/coordinates.cpp
    src/gap_map.cpp
    src/indexing.cpp
    src/alignment.cpp
    src/placement.cpp
    src/seed_annotated_tree.cpp
    src/seeding.cpp
    src/mm_align.c
    src/pileup.c
    src/fixed_kmer.hpp
    src/fixed_kmer.cpp
    src/performance.cpp
)
if(OPTION_BUILD_TESTS)
    add_executable(tests
        src/test/main.test.cpp
        src/test/pmi.test.cpp
        src/conversion.cpp
        src/indexing.cpp
        src/alignment.cpp
        src/coordinates.cpp
        src/gap_map.cpp
        src/placement.cpp
        src/seed_annotated_tree.cpp
        src/seeding.cpp
        src/mm_align.c
        src/pileup.c
        src/fixed_kmer.hpp
        src/fixed_kmer.cpp
        src/performance.cpp
    )
endif()

FetchContent_GetProperties(ftxui)
if(NOT ftxui_POPULATED)
    FetchContent_Populate(ftxui)
    add_subdirectory(${ftxui_SOURCE_DIR} ${ftxui_BINARY_DIR})
endif()

# Add include directories for spdlog
target_include_directories(panmap PRIVATE ${spdlog_SOURCE_DIR}/include)
target_include_directories(panmap PRIVATE ${ftxui_SOURCE_DIR}/include)

target_include_directories(panmap PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/_deps/panman-build ${CMAKE_CURRENT_BINARY_DIR}/_deps/panman-src/src ${CMAKE_CURRENT_BINARY_DIR}/_deps/panman-src)


include(cmake/detect-arch.cmake)
detect_host_arch()
message("neon: ${HAS_NEON}")
message("Third party dependencies...")

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/lib/libminimap2.a
    COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/lib ${CMAKE_CURRENT_BINARY_DIR}/bin
    COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/minimap2 && make clean
    COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/minimap2 && make $<$<BOOL:${HAS_NEON}>:aarch64=1>
    COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/minimap2/libminimap2.a ${CMAKE_CURRENT_BINARY_DIR}/lib/libminimap2.a
    COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/minimap2/minimap2 ${CMAKE_CURRENT_BINARY_DIR}/bin/minimap2
    COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/minimap2 && rm -f *.o
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

# Define source directory for samtools build to ensure consistency
set(SAMTOOLS_SOURCE_DIR ${CMAKE_SOURCE_DIR})

# Verify SAMTOOLS_SOURCE_DIR is not empty
if("${SAMTOOLS_SOURCE_DIR}" STREQUAL "")
    message(FATAL_ERROR "SAMTOOLS_SOURCE_DIR is empty. This will cause configure scripts to fail.")
endif()

# Try to find autoreconf in common locations
find_program(AUTORECONF_PATH autoreconf PATHS /usr/bin /usr/local/bin /opt/local/bin /opt/homebrew/bin)
if(NOT AUTORECONF_PATH)
    message(FATAL_ERROR "autoreconf not found! Please install autoconf or specify AUTORECONF_PATH manually.")
endif()

# Verify CMAKE_CURRENT_BINARY_DIR is not empty
if("${CMAKE_CURRENT_BINARY_DIR}" STREQUAL "")
    message(FATAL_ERROR "CMAKE_CURRENT_BINARY_DIR is empty. This will cause configure scripts to fail.")
endif()

# Define absolute paths for scripts
set(HTSLIB_CONFIG_SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/htslib_config.sh)
set(SAMTOOLS_CONFIG_SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/samtools_config.sh)

# Generate htslib_config.sh at configure time with proper Unix line endings
file(WRITE ${HTSLIB_CONFIG_SCRIPT} "#!/bin/bash\n")
file(APPEND ${HTSLIB_CONFIG_SCRIPT} "cd \"${SAMTOOLS_SOURCE_DIR}/src/3rdparty/samtools/htslib-1.20\"\n")
file(APPEND ${HTSLIB_CONFIG_SCRIPT} "autoreconf -i\n")
file(APPEND ${HTSLIB_CONFIG_SCRIPT} "autoconf\n")
file(APPEND ${HTSLIB_CONFIG_SCRIPT} "./configure --prefix=\"${CMAKE_CURRENT_BINARY_DIR}\" --disable-lzma --disable-bz2 --disable-libcurl\n")

# Ensure the script is executable
file(CHMOD ${HTSLIB_CONFIG_SCRIPT} PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

add_custom_command(
    OUTPUT ${SAMTOOLS_CONFIG_SCRIPT}
    COMMAND ${CMAKE_COMMAND} -E echo "#!/bin/bash" > ${SAMTOOLS_CONFIG_SCRIPT}
    COMMAND ${CMAKE_COMMAND} -E echo "cd ${SAMTOOLS_SOURCE_DIR}/src/3rdparty/samtools" >> ${SAMTOOLS_CONFIG_SCRIPT}
    COMMAND ${CMAKE_COMMAND} -E echo "autoreconf -i" >> ${SAMTOOLS_CONFIG_SCRIPT}
    COMMAND ${CMAKE_COMMAND} -E echo "autoconf" >> ${SAMTOOLS_CONFIG_SCRIPT}
    COMMAND ${CMAKE_COMMAND} -E echo "./configure --with-htslib=${CMAKE_CURRENT_BINARY_DIR} --prefix=${CMAKE_CURRENT_BINARY_DIR} --disable-lzma --disable-bz2 --disable-libcurl --disable-shared --with-pic" >> ${SAMTOOLS_CONFIG_SCRIPT}
    COMMAND chmod +x ${SAMTOOLS_CONFIG_SCRIPT}
)

# Ensure scripts are generated before running samtools build
add_custom_target(generate_samtools_scripts DEPENDS ${HTSLIB_CONFIG_SCRIPT} ${SAMTOOLS_CONFIG_SCRIPT})

# Update your existing custom command to depend on these scripts
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/lib/libst.a ${CMAKE_CURRENT_BINARY_DIR}/lib/libhts.a
    DEPENDS generate_samtools_scripts
    COMMAND echo "[samtools build] Step 1: Checking for tools"
    COMMAND echo "Autoreconf path: ${AUTORECONF_PATH}"
    COMMAND echo "Source dir: ${SAMTOOLS_SOURCE_DIR}"
    COMMAND echo "Current binary dir: ${CMAKE_CURRENT_BINARY_DIR}"
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/lib ${CMAKE_CURRENT_BINARY_DIR}/bin

    COMMAND echo "[samtools build] Step 2: Running autoreconf"
    COMMAND cd ${SAMTOOLS_SOURCE_DIR}/src/3rdparty/samtools/htslib-1.20 && ${AUTORECONF_PATH} -i || echo "Warning: autoreconf failed, configure script may be missing"

    COMMAND echo "[samtools build] Step 3: Configuring htslib"
    COMMAND bash ${HTSLIB_CONFIG_SCRIPT}

    COMMAND echo "[samtools build] Step 4: Building htslib"
    COMMAND cd ${SAMTOOLS_SOURCE_DIR}/src/3rdparty/samtools/htslib-1.20 && make -j

    COMMAND echo "[samtools build] Step 5: Installing htslib"
    COMMAND cd ${SAMTOOLS_SOURCE_DIR}/src/3rdparty/samtools/htslib-1.20 && make install

    COMMAND echo "[samtools build] Step 6: Running autoheader for samtools"
    COMMAND cd ${SAMTOOLS_SOURCE_DIR}/src/3rdparty/samtools && ${AUTORECONF_PATH} -i

    COMMAND echo "[samtools build] Step 7: Running autoconf for samtools"
    COMMAND cd ${SAMTOOLS_SOURCE_DIR}/src/3rdparty/samtools && ${AUTORECONF_PATH}

    COMMAND echo "[samtools build] Step 8: Configuring samtools"
    COMMAND bash ${SAMTOOLS_CONFIG_SCRIPT}

    COMMAND echo "[samtools build] Step 9: Building samtools"
    COMMAND cd ${SAMTOOLS_SOURCE_DIR}/src/3rdparty/samtools && make -j

    COMMAND echo "[samtools build] Step 10: Copying samtools library"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${SAMTOOLS_SOURCE_DIR}/src/3rdparty/samtools/libst.a
    ${CMAKE_CURRENT_BINARY_DIR}/lib/libst.a
    || echo "Error: Could not find samtools library"

    COMMAND echo "[samtools build] Step 11: Completed samtools build"
    BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/lib/libhts.a ${CMAKE_CURRENT_BINARY_DIR}/lib/libst.a
    WORKING_DIRECTORY ${SAMTOOLS_SOURCE_DIR}
)

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/bin/capnp ${CMAKE_CURRENT_BINARY_DIR}/bin/capnpc-c++ ${CMAKE_CURRENT_BINARY_DIR}/lib/libcapnp.a ${CMAKE_CURRENT_BINARY_DIR}/lib/libkj.a ${CMAKE_CURRENT_BINARY_DIR}/include/capnp
    COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/lib ${CMAKE_CURRENT_BINARY_DIR}/bin ${CMAKE_CURRENT_BINARY_DIR}/include
    # Skip autoconf/autoheader and try using cmake directly for capnproto
    COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/capnproto/c++ && mkdir -p build && cd build && cmake -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR} -DBUILD_TESTING=OFF -DBUILD_SHARED_LIBS=OFF ..
    COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/capnproto/c++/build && make -j
    COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/capnproto/c++/build && make install
    COMMENT "Compiling capnproto using CMake"
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
    ${CMAKE_CURRENT_BINARY_DIR}/bin/capnp
    ${CMAKE_CURRENT_BINARY_DIR}/lib/libcapnp.a
    ${CMAKE_CURRENT_BINARY_DIR}/lib/libkj.a
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

if(OPTION_BUILD_TESTS)
    target_include_directories(tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty)
    target_include_directories(tests PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
    target_include_directories(tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/samtools/htslib-1.20)
endif()

if(OPTION_BUILD_SIMULATE)
    add_executable(simulate
        src/test/simulate.cpp
        src/conversion.cpp
        src/indexing.cpp
        src/placement.cpp
        src/alignment.cpp
        src/coordinates.cpp
        src/gap_map.cpp
        src/genotyping.cpp
        src/seed_annotated_tree.cpp
        src/seeding.cpp
        src/mm_align.c
        src/pileup.c
        src/fixed_kmer.hpp
        src/fixed_kmer.cpp
        src/performance.cpp
    )

    target_link_libraries(simulate PRIVATE
        bwa
        bcftools
        minimap2
        samtools
        htslib
        spdlog::spdlog
        ftxui::screen
        ftxui::dom
        ftxui::component
        docopt
        capnp
        kj
        Boost::iostreams
        Boost::filesystem
        Boost::system
        ZLIB::ZLIB
        panman
        tbb
        tbbmalloc
        tbbmalloc_proxy
        deflate
        m
    )
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

target_link_libraries(panmap PRIVATE
    bwa
    bcftools
    minimap2
    samtools
    htslib
    docopt
    capnp
    spdlog::spdlog
    ftxui::screen
    ftxui::dom
    ftxui::component
    kj
    Boost::iostreams
    Boost::filesystem
    Boost::system
    ZLIB::ZLIB
    panman
    tbb
    tbbmalloc
    tbbmalloc_proxy
    deflate
    m)
target_link_options(panmap PRIVATE -lm -lpthread)
target_compile_options(panmap PRIVATE -Wall)
target_compile_features(panmap PRIVATE cxx_std_20)


if(OPTION_BUILD_TESTS)
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
    target_link_options(tests PUBLIC -lm -lpthread)
    target_compile_options(tests PUBLIC -g -O3 -ggdb)
endif()


if(OPTION_DEBUG)
    target_compile_options(panmap PUBLIC -g -O0 -ggdb)
else()
    # Use higher optimization levels and native CPU features for better k-mer performance
    target_compile_options(panmap PUBLIC -g -pg -O3 -march=native -mtune=native)
endif()

# If tests are enabled, add optimization flags
if(OPTION_BUILD_TESTS)
    target_compile_options(tests PUBLIC -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
    target_link_options(tests PUBLIC -lm -lpthread)

    if(OPTION_DEBUG)
        target_compile_options(tests PUBLIC -g -O0 -ggdb)
    else()
        # Use higher optimization levels and native CPU features for better k-mer performance
        target_compile_options(tests PUBLIC -g -O3 -ggdb -march=native -mtune=native)
    endif()
endif()

# If simulate is enabled, add optimization flags
if(OPTION_BUILD_SIMULATE)
    if(OPTION_DEBUG)
        target_compile_options(simulate PUBLIC -g -O0 -ggdb)
    else()
        # Use higher optimization levels and native CPU features for better k-mer performance
        target_compile_options(simulate PUBLIC -g -O3 -march=native -mtune=native)
    endif()
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

# Print build configuration summary
message(STATUS "")
message(STATUS "Build Configuration Summary:")
message(STATUS "----------------------------")
message(STATUS "CMake Version: ${CMAKE_VERSION}")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "Architecture: ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "")
message(STATUS "Optimizations:")
message(STATUS "  Fixed k-mer optimization: ${OPTION_ENABLE_FIXED_KMER}")
if(OPTION_ENABLE_FIXED_KMER)
    message(STATUS "  TBB memory allocator: Enabled")
    message(STATUS "  TBB memory pool size: 1GB")
    message(STATUS "  TBB thread cache size: 16MB")
else()
    message(STATUS "  TBB memory allocator: Default")
endif()
if(DEFINED ARCH_SIMD_FLAGS)
    message(STATUS "  SIMD flags: ${ARCH_SIMD_FLAGS}")
else()
    message(STATUS "  SIMD flags: None")
endif()
message(STATUS "  Debug mode: ${OPTION_DEBUG}")
message(STATUS "  Profiling: ${OPTION_PROFILE}")
message(STATUS "")
message(STATUS "Components:")
message(STATUS "  Build tests: ${OPTION_BUILD_TESTS}")
message(STATUS "  Build simulate: ${OPTION_BUILD_SIMULATE}")
message(STATUS "----------------------------")
message(STATUS "")

# Configure compiler flags for SIMD support and k-mer optimizations
include(CheckCXXCompilerFlag)

# Detect architecture
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64)")
    set(IS_ARM TRUE)
    message(STATUS "Detected ARM architecture: ${CMAKE_SYSTEM_PROCESSOR}")
else()
    set(IS_ARM FALSE)
    message(STATUS "Detected x86 architecture: ${CMAKE_SYSTEM_PROCESSOR}")
endif()

# Check for various SIMD instruction sets across different architectures
if(IS_ARM)
    # ARM-specific optimizations
    check_cxx_compiler_flag("-mcpu=apple-m1" COMPILER_SUPPORTS_APPLE_M1)
    check_cxx_compiler_flag("-mcpu=apple-m2" COMPILER_SUPPORTS_APPLE_M2)
    check_cxx_compiler_flag("-mfpu=neon" COMPILER_SUPPORTS_NEON)

    if(COMPILER_SUPPORTS_APPLE_M1)
        message(STATUS "Compiler supports Apple M1 - enabling optimized ARM processing")
        set(ARCH_SIMD_FLAGS "${ARCH_SIMD_FLAGS} -mcpu=apple-m1")
    elseif(COMPILER_SUPPORTS_APPLE_M2)
        message(STATUS "Compiler supports Apple M2 - enabling optimized ARM processing")
        set(ARCH_SIMD_FLAGS "${ARCH_SIMD_FLAGS} -mcpu=apple-m2")
    elseif(COMPILER_SUPPORTS_NEON)
        message(STATUS "Compiler supports NEON - enabling ARM-specific optimizations")
        set(ARCH_SIMD_FLAGS "${ARCH_SIMD_FLAGS} -mfpu=neon")
    endif()
else()
    # x86 optimizations
    check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)
    check_cxx_compiler_flag("-msse4.1" COMPILER_SUPPORTS_SSE4_1)

    if(COMPILER_SUPPORTS_AVX2)
        message(STATUS "Compiler supports AVX2 - enabling optimized k-mer processing")
        set(ARCH_SIMD_FLAGS "${ARCH_SIMD_FLAGS} -mavx2")
    elseif(COMPILER_SUPPORTS_SSE4_1)
        message(STATUS "Compiler supports SSE4.1 - enabling optimized k-mer processing")
        set(ARCH_SIMD_FLAGS "${ARCH_SIMD_FLAGS} -msse4.1")
    endif()
endif()

# Add SIMD flags to the compilation if supported
if(DEFINED ARCH_SIMD_FLAGS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ARCH_SIMD_FLAGS}")
endif()

# Add fixed k-mer optimization definition if enabled
if(OPTION_ENABLE_FIXED_KMER)
    message(STATUS "Enabling optimized fixed k-mer implementations (8, 16, 32)")
    add_compile_definitions(ENABLE_FIXED_KMER_OPTIMIZATION)

    # Enable tbbmalloc as the global memory allocator for better performance
    add_compile_definitions(TBB_USE_THREADING_TOOLS=1)
    add_compile_definitions(USE_TBBMALLOC=1)

    # Configure TBB memory pools for optimized k-mer processing
    add_compile_definitions(TBB_MALLOC_USE_HUGE_PAGES=1) # Use huge pages where available
    add_compile_definitions(TBB_MALLOC_SET_THREAD_CACHE_SIZE=16384) # Larger thread cache (16MB)
    add_compile_definitions(TBB_MALLOC_POOL_SIZE=1073741824) # 1GB memory pool for better performance
endif()
