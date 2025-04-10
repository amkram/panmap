# Path to the file we want to modify
# Try multiple potential locations since the path is ambiguous
set(POTENTIAL_PATHS
  "${CMAKE_BINARY_DIR}/_deps/panman-src/CMakeLists.txt"
  "${CMAKE_CURRENT_BINARY_DIR}/_deps/panman-src/CMakeLists.txt"
  "${CMAKE_CURRENT_BINARY_DIR}/panman-src/CMakeLists.txt"
  "${CMAKE_SOURCE_DIR}/build/_deps/panman-src/CMakeLists.txt"
)

set(FILE_TO_PATCH "")
foreach(PATH ${POTENTIAL_PATHS})
  message(STATUS "Checking for panman CMakeLists.txt at: ${PATH}")
  if(EXISTS "${PATH}")
    set(FILE_TO_PATCH "${PATH}")
    message(STATUS "Found panman CMakeLists.txt at: ${FILE_TO_PATCH}")
    break()
  endif()
endforeach()

if("${FILE_TO_PATCH}" STREQUAL "")
  message(FATAL_ERROR "Unable to find panman CMakeLists.txt in any of the potential locations")
endif()

# Read the file
file(READ "${FILE_TO_PATCH}" CONTENT)

# Add comment markers to disable TBB related lines
string(REPLACE "include(${TBB_DIR}/cmake/TBBBuild.cmake)" 
       "# include(${TBB_DIR}/cmake/TBBBuild.cmake) # Commented out by panmap" 
       CONTENT "${CONTENT}")
       
string(REPLACE "tbb_build(TBB_ROOT ${TBB_DIR} CONFIG_DIR TBB_DIR MAKE_ARGS tbb_cpf=1)" 
       "# tbb_build(TBB_ROOT ${TBB_DIR} CONFIG_DIR TBB_DIR MAKE_ARGS tbb_cpf=1) # Commented out by panmap" 
       CONTENT "${CONTENT}")
       
string(REPLACE "find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview)" 
       "# find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview) # Commented out by panmap" 
       CONTENT "${CONTENT}")

# Add comment markers to disable JsonCpp related lines
string(REPLACE "include(${CMAKE_TOOLCHAIN_FILE})" 
       "# include(${CMAKE_TOOLCHAIN_FILE}) # Commented out by panmap" 
       CONTENT "${CONTENT}")
       
string(REPLACE "find_package(jsoncpp CONFIG REQUIRED)" 
       "# find_package(jsoncpp CONFIG REQUIRED) # Commented out by panmap" 
       CONTENT "${CONTENT}")

# Write the modified content back to the file
message(STATUS "Writing patched panman CMakeLists.txt")
file(WRITE "${FILE_TO_PATCH}" "${CONTENT}")

message(STATUS "Panman CMakeLists.txt patched successfully") 