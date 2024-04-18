cmake_minimum_required (VERSION 3.20)

project(panmap)

<<<<<<< HEAD
set(CMAKE_CXX_STANDARD 20)
=======
# Project options
option(OPTION_DEBUG "Compile with debug symbols and no optimizations." ON)
option(OPTION_BUILD_TESTS "Build tests." ON)
option(OPTION_BUILD_SIMULATE "Build simulate." ON)
>>>>>>> be62fc15e57015f876d9a33861fa2b1ea98da3fa

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

<<<<<<< HEAD
target_link_libraries(panmap PUBLIC samtools htslib panmat minimap2)
target_compile_options(panmap PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
=======

file(GLOB PANMAP_SRCS "src/*.c*" "src/*.h*" "src/minimap2_src/libminimap2.a" "src/samtools_src/sample.c" "src/samtools_src/bedidx.c" "src/samtools_src/sam_opts.c" "src/samtools_src/bam_plbuf.c" "src/samtools_src/sam_utils.c" "${CMAKE_CURRENT_SOURCE_DIR}/build/htslib_install/lib/libhts.a")
add_executable(panmap ${PANMAP_SRCS})
add_dependencies(panmap libminimap2.a)
add_dependencies(panmap sample.o)


if(OPTION_BUILD_TESTS)
    get_filename_component(main_src src/panmap.cpp ABSOLUTE)
    list(REMOVE_ITEM PANMAP_SRCS "${main_src}")
    add_executable(tests ${PANMAP_SRCS} "src/test/main.test.cpp" "src/test/pmi.test.cpp" "src/test/seed.test.cpp" "src/test/tree.test.cpp" "src/test/genotype.test.cpp" "src/test/conversions.test.cpp") 
    add_dependencies(tests libminimap2.a)
    add_dependencies(tests sample.o)
endif(OPTION_BUILD_TESTS)

if(OPTION_BUILD_SIMULATE)
    get_filename_component(main_src src/panmap.cpp ABSOLUTE)
    list(REMOVE_ITEM PANMAP_SRCS "${main_src}")
    add_executable(simulate ${PANMAP_SRCS} "src/test/simulate.cpp")
    add_dependencies(simulate libminimap2.a)
    add_dependencies(simulate sample.o)
endif(OPTION_BUILD_SIMULATE)



target_link_libraries(panmap PRIVATE stdc++ ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} remote-deps libpanmat ${CMAKE_SOURCE_DIR}/src/minimap2_src/libminimap2.a z ${CMAKE_CURRENT_SOURCE_DIR}/build/htslib_install/lib/libhts.dylib)
target_link_libraries(simulate PRIVATE stdc++ ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} remote-deps libpanmat ${CMAKE_SOURCE_DIR}/src/minimap2_src/libminimap2.a z ${CMAKE_CURRENT_SOURCE_DIR}/build/htslib_install/lib/libhts.dylib)
target_link_libraries(tests PRIVATE stdc++ ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} libpanmat ${CMAKE_SOURCE_DIR}/src/minimap2_src/libminimap2.a z ${CMAKE_CURRENT_SOURCE_DIR}/build/htslib_install/lib/libhts.dylib)
target_compile_options(panmap PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
target_compile_options(simulate PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
target_compile_options(tests PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_BINARY_DIR}/clean-minimap)
>>>>>>> be62fc15e57015f876d9a33861fa2b1ea98da3fa

# Unit testing data
configure_file(src/test/data/test.aligned.fa test.aligned.fa COPYONLY)
configure_file(src/test/data/test.nwk test.nwk COPYONLY)
configure_file(src/test/data/test.json test.json COPYONLY)
configure_file(src/test/data/test.pmat test.pmat COPYONLY)
configure_file(src/test/data/test.fastq test.fastq COPYONLY)
