set(VERSION_MAJOR 3)
set(VERSION_MINOR 5)
set(VERSION_PATCH 6)
set(BUILD_NUMBER 49)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/version.rc.in
               ${CMAKE_BINARY_DIR}/version.rc @ONLY)
