# Dummy TBBConfig.cmake file to satisfy panman's dependency
# This file defines the necessary TBB targets without actually building TBB

if(NOT TARGET TBB::tbb)
  add_library(TBB::tbb INTERFACE IMPORTED)
endif()

if(NOT TARGET TBB::tbbmalloc)
  add_library(TBB::tbbmalloc INTERFACE IMPORTED)
endif()

if(NOT TARGET TBB::tbbmalloc_proxy)
  add_library(TBB::tbbmalloc_proxy INTERFACE IMPORTED)
endif()

if(NOT TARGET TBB::tbb_preview)
  add_library(TBB::tbb_preview INTERFACE IMPORTED)
endif()

# Set found flag for find_package
set(TBB_FOUND TRUE)

# Define TBB imported targets variable
set(TBB_IMPORTED_TARGETS TBB::tbb TBB::tbbmalloc TBB::tbbmalloc_proxy TBB::tbb_preview)

message(STATUS "Using dummy TBB configuration (TBB handled by parent project)") 