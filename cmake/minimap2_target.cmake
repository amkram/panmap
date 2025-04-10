add_custom_target(minimap2_target 
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/lib
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/include
    COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/minimap2 && make -j$(nproc) libminimap2.a
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/minimap2/libminimap2.a ${CMAKE_CURRENT_BINARY_DIR}/lib/
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/minimap2/*.h ${CMAKE_CURRENT_BINARY_DIR}/include/
    COMMENT "Building minimap2 library from src/3rdparty/minimap2"
    SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/minimap2/Makefile
)
