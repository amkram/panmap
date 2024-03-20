cmake_minimum_required (VERSION 3.26 FATAL_ERROR) 

project(panmap)

# Project options
option(OPTION_DEBUG "Compile with debug symbols and no optimizations." ON)
option(OPTION_BUILD_TESTS "Build tests." ON)

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 ") 
if(OPTION_DEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wno-unused-function -Wno-deprecated-declarations -Wno-deprecated-volatile -g -ggdb3") 
else(OPTION_DEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3") 
endif(OPTION_DEBUG)

# Cmake options
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(FetchContent)
FetchContent_Declare(docopt
    GIT_REPOSITORY https://github.com/docopt/docopt.cpp.git
    GIT_TAG v0.6.3)
FetchContent_Populate(docopt)
include(CMakePrintHelpers)
cmake_print_variables(docopt_SOURCE_DIR)
cmake_print_variables(docopt_BINARY_DIR)

add_library(remote-deps ${docopt_SOURCE_DIR}/docopt.cpp)

# Dependency: panmat-utils
find_package(libpanmat REQUIRED)
include_directories(${libpanmat_DIR}/include ${libpanmat_DIR}/oneTBB-2019_U9/include)


# Dependency: Boost
find_package(Boost COMPONENTS system filesystem program_options unit_test_framework iostreams date_time REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# Dependency: minimap2
add_custom_target(libminimap2.a
   COMMAND make arm_neon=1 aarch64=1 libminimap2.a
   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/minimap2_src/
)
add_custom_target(clean-minimap
    COMMAND make clean
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/minimap2_src/
)

# Dependency: samtools
add_custom_target(sample.o
   COMMAND make arm_neon=1 aarch64=1 sample.o bedidx.o sam_opts.o bam_plbuf.o sam_utils.o version.h
   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/samtools_src/
)

# Dependency: htslib
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/build/htslib_install/include)


file(GLOB PANMAP_SRCS "src/*.c*" "src/*.h*" "src/minimap2_src/libminimap2.a" "src/samtools_src/sample.c" "src/samtools_src/bedidx.c" "src/samtools_src/sam_opts.c" "src/samtools_src/bam_plbuf.c" "src/samtools_src/sam_utils.c" "${CMAKE_CURRENT_SOURCE_DIR}/build/htslib_install/lib/libhts.a")
add_executable(panmap ${PANMAP_SRCS})
add_dependencies(panmap libminimap2.a)
add_dependencies(panmap sample.o)



if(OPTION_BUILD_TESTS)
    get_filename_component(main_src src/panmap.cpp ABSOLUTE)
    list(REMOVE_ITEM PANMAP_SRCS "${main_src}")
    add_executable(tests ${PANMAP_SRCS} "src/test/main.test.cpp" "src/test/pmi.test.cpp" "src/test/seed.test.cpp" "src/test/tree.test.cpp") 
    add_dependencies(tests libminimap2.a)
    add_dependencies(tests sample.o)
endif(OPTION_BUILD_TESTS)


target_link_libraries(panmap PRIVATE stdc++ ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} remote-deps libpanmat ${CMAKE_SOURCE_DIR}/src/minimap2_src/libminimap2.a z hts)
target_link_libraries(tests PRIVATE stdc++ ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} libpanmat ${CMAKE_SOURCE_DIR}/src/minimap2_src/libminimap2.a z hts)
target_compile_options(panmap PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
target_compile_options(tests PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_BINARY_DIR}/clean-minimap)

# Unit testing data
configure_file(src/test/data/test.aligned.fa test.aligned.fa COPYONLY)
configure_file(src/test/data/test.nwk test.nwk COPYONLY)
configure_file(src/test/data/test.json test.json COPYONLY)
configure_file(src/test/data/test.pmat test.pmat COPYONLY)
configure_file(src/test/data/test.fastq test.fastq COPYONLY)