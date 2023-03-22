cmake_minimum_required(VERSION 3.23)

project(
  MeshKernelLibrary
  VERSION 1.0.0
  DESCRIPTION "Library for creating and editing meshes."
  LANGUAGES CXX C
)

include(cmake/options.cmake)
include(cmake/compiler_config.cmake)
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
