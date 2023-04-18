cmake_minimum_required(VERSION 3.23)

include(cmake/get_version_from_git.cmake)

# get version from latest tag
get_version_from_git(
  LOG OFF
  TIMESTAMP "%Y%m%d%H%M%S"
)

project(
  MeshKernel
  VERSION ${VERSION_STRING_BASE}
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

# Run packaging scripts
add_subdirectory(package)

# Tools
add_subdirectory(tools)

# Add third party directory
add_subdirectory(extern)

# Libraries
add_subdirectory(libs)
