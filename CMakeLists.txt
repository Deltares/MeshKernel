cmake_minimum_required(VERSION 3.23)

set(MESHKERNEL_VERSION 4.3.8)

project(
  MeshKernel
  VERSION ${MESHKERNEL_VERSION}
  DESCRIPTION "Library for creating and editing meshes."
  LANGUAGES CXX C
)

# get options
include(cmake/options.cmake)

# configure the compiler
include(cmake/compiler_config.cmake)

# fetch dependencies, must appear after options.cmake
include(cmake/fetch_content.cmake)

# organize targets into folders, can be remove in version > 3.26 (ON by default)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Docs only available if this is the main app
add_subdirectory(docs)

# Tools
# Note: tools/Version sets the semantic version
add_subdirectory(tools)

# Run packaging scripts: requires semantic version
add_subdirectory(package)

# Add third party dependencies
add_subdirectory(extern)

# Libraries
add_subdirectory(libs)
