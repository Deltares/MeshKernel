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
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -O2 -std=c++17")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -std=c++17")
  else()
    set(CMAKE_CXX_FLAGS_RELEASE
        "${CMAKE_CXX_FLAGS_RELEASE} /EHsc /MP /std:c++17")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /EHsc /MP /std:c++17")
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
include_directories(${Boost_INCLUDE_DIR})

# NetCDF:
# Even though only the dlls are needed, use find_package to cleanly determine:
# - the install prefix, which allows constructing and setting the path to the bin directory (win only)
# - the include directory
# Invoke find_package quietly to allow issuing a descriptive configuration error
find_package(netCDF QUIET)
if (netCDF_FOUND)
  message(STATUS "Found NetCDF: in ${netCDF_INSTALL_PREFIX} (found version \"${NetCDFVersion}\")")
  if(WIN32)
    set (netCDF_BIN_DIR "${netCDF_INSTALL_PREFIX}/bin")
  endif()
  include_directories(${netCDF_INCLUDE_DIR})
else()
  message(FATAL_ERROR "Could not find NetCDF. Please make sure it is installed and available in the system path.")
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
