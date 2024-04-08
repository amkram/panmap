if(WIN32)
  find_package(Git)
  if(Git_FOUND)
    get_filename_component(GIT_DIR ${GIT_EXECUTABLE} DIRECTORY)
    get_filename_component(GIT_DIR ${GIT_DIR} DIRECTORY)
  endif()
endif()

find_program(PATCH
NAMES patch
HINTS ${GIT_DIR}
PATH_SUFFIXES usr/bin
)

if(NOT PATCH)
  message(FATAL_ERROR "Did not find GNU Patch")
endif()

execute_process(COMMAND yes no COMMAND ${PATCH} ${in_file} --input=${patch_file} --output=${out_file} --ignore-whitespace
    TIMEOUT 15
    COMMAND_ECHO STDOUT
    RESULT_VARIABLE ret
)
