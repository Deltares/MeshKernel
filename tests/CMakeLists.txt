# Fetch google test
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG release-1.11.0)
  
FetchContent_GetProperties(googletest)
if(NOT googletest_POPULATED)
  FetchContent_Populate(googletest)
  add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

set(TEST_FOLDER "${CMAKE_CURRENT_SOURCE_DIR}")
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/utils/test_utils/include/TestUtils/Definitions.hpp.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/utils/test_utils/include/TestUtils/Definitions.hpp"
)

# handle benchmarking options
option(
  ENABLE_BENCHMARKS
  "Enables building and running the performance benchmarking executable."
  OFF
)

include(CMakeDependentOption)
cmake_dependent_option(
  ENABLE_MEM_REPORT
  "When performance benchmarking is enabled, enables reporting memory usage statistics."
  ON 
  "ENABLE_BENCHMARKS"
  OFF
)

if(ENABLE_BENCHMARKS)
  # check if it's a release or release with deug information build
  set(VALID_BUILD_TYPES "RELEASE" "RELWITHDEBINFO")
  # upcase the build type to make the option case-insensitive
  string(TOUPPER "${CMAKE_BUILD_TYPE}" UPCASED_CMAKE_BUILD_TYPE)
  if(UPCASED_CMAKE_BUILD_TYPE IN_LIST VALID_BUILD_TYPES)
    # Fetch google benchmark    
    FetchContent_Declare(
      googlebenchmark
      GIT_REPOSITORY https://github.com/google/benchmark.git
      GIT_TAG v1.7.1
    )

    FetchContent_GetProperties(benchmark)

    # Prevent the google benchmark's own tests from running
    set(BENCHMARK_ENABLE_TESTING OFF)

    if(NOT googlebenchmark_POPULATED)
      FetchContent_Populate(googlebenchmark)
      add_subdirectory(
        ${googlebenchmark_SOURCE_DIR}
        ${googlebenchmark_BINARY_DIR}
        EXCLUDE_FROM_ALL
      )
    endif()
  else()
    # disable all benchmarking options
    set(ENABLE_BENCHMARKS OFF)
    set(ENABLE_MEM_REPORT OFF)
    message(
      WARNING
      "The benchmarks and their depenedencies can be built only if the build is configured "
      "with CMAKE_BUILD_TYPE set to Release or RelWithDebInfo. "
      "The current build is configured with CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}."
      "All benchmarking configuration options are ignored."
    )
  endif()
endif()

add_subdirectory(utils)
add_subdirectory(unit)
add_subdirectory(api)
add_subdirectory(benchmark)
