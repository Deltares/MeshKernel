# Works with 3.14 and tested through 3.18
cmake_minimum_required(VERSION 3.14...3.18)

# Project name and a few useful settings. Other commands can pick up the results
project(
  MeshKernel
  VERSION 1.0.0
  DESCRIPTION "Library for creating and editing meshes."
  LANGUAGES CXX C)
  
option(PRODUCE_CODE_COVERAGE "Produce code coverage files under GNU compilers" OFF)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_compile_options("-std=c++17;-Werror;-Wall;-Wextra;-pedantic;-Wno-unused-function")
    add_compile_options("$<$<CONFIG:RELEASE>:-O2>")
    add_compile_options("$<$<CONFIG:DEBUG>:-g>")
  else()
    add_compile_options("/EHsc;/MP;/std:c++17;/W3;/WX")
    add_compile_options("$<$<CONFIG:RELEASE>:/O2>")
    add_compile_options("$<$<CONFIG:DEBUG>:/Od;/DEBUG>")
    add_compile_definitions("_USE_MATH_DEFINES")
  endif()

  # Create position independent binaries
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)

  # Disable compiler specific extensions
  set(CMAKE_CXX_EXTENSIONS OFF)

  # Let's nicely support folders in IDEs
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)

  # Note this needs to be done in the main CMakeLists since it calls
  # enable_testing, which must be in the main CMakeLists.
  include(CTest)

  # Docs only available if this is the main app
  add_subdirectory(docs)
endif()

# Determine version suffix from environment variable
if(DEFINED ENV{VERSION_SUFFIX})
  set(VERSION_SUFFIX $ENV{VERSION_SUFFIX})
else()
  set(VERSION_SUFFIX .0)
endif()

# Run packaging scripts
add_subdirectory(package)

# Add third party directory
add_subdirectory(extern)

# FetchContent added in CMake 3.11, downloads during the configure step
include(FetchContent)

# Boost
if(WIN32)
  # In windows use static libraries
  set(Boost_USE_STATIC_LIBS ON)
endif()
find_package(Boost REQUIRED COMPONENTS system filesystem)
include_directories(SYSTEM ${Boost_INCLUDE_DIR})

# NetCDF
if(WIN32)
  find_path(
    NETCDF_BIN_PATH
    NAMES "${CMAKE_SHARED_LIBRARY_PREFIX}netcdf${CMAKE_SHARED_LIBRARY_SUFFIX}"
    HINTS "${CMAKE_CURRENT_SOURCE_DIR}/extern/netcdf/netCDF 4.6.1/bin")
  message(STATUS "netcdf bin dir is ${NETCDF_BIN_PATH}")
endif()

# OpenMP
find_package(OpenMP REQUIRED)

# The version file
add_subdirectory(include/Version)

# The static library
add_subdirectory(src/MeshKernel)

# The dynamic library
add_subdirectory(src/MeshKernelApi)

# Testing only available if this is the main app.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  add_subdirectory(tests)
endif()
