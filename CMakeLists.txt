cmake_minimum_required(VERSION 3.14...3.18)
project(MeshKernel)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

if(WIN32)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /EHsc /MP")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /EHsc /MP")
else()
  set(CMAKE_CXX_FLAGS_RELEASE
      "-g -O2"
      CACHE STRING "List of C++ compiler flags for a Release build")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -g")
endif(WIN32)

# Note this needs to be done in the main CMakeLists since it calls
# enable_testing, which must be in the main CMakeLists.
include(CTest)

# Docs only available if this is the main app
find_package(Doxygen)
if(Doxygen_FOUND)
  add_subdirectory(docs)
else()
  message(STATUS "Doxygen not found, not building docs")
endif()

# Add third party directory
add_subdirectory(thirdParty)

# FetchContent added in CMake 3.11, downloads during the configure step
include(FetchContent)

# Boost C++
find_package(Boost)
include_directories(${Boost_INCLUDE_DIR})

# The compiled library code is here
add_subdirectory(src)

# Add target link dependency on OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(MeshKernel PUBLIC OpenMP::OpenMP_CXX)
endif()

# Testing only available if this is the main app Emergency override
# MODERN_CMAKE_BUILD_TESTING provided as well if((CMAKE_PROJECT_NAME STREQUAL
# PROJECT_NAME OR MODERN_CMAKE_BUILD_TESTING) AND BUILD_TESTING)
add_subdirectory(tests)
# endif()
