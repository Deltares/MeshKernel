include(FetchContent)

if(ENABLE_UNIT_TESTING)
  # Fetch google test
  FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.13.0)

  FetchContent_GetProperties(googletest)
  if(NOT googletest_POPULATED)
    FetchContent_Populate(googletest)
    add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif()

  # Note this needs to be done in the main CMakeLists since it calls
  # enable_testing, which must be in the main CMakeLists.
  include(CTest)
endif()

if(ENABLE_BENCHMARKING)
  # check if it's a release or release with debug information build
  set(VALID_BUILD_TYPES "RELEASE" "RELWITHDEBINFO")
  # upcase the build type to make the option case-insensitive
  string(TOUPPER "${CMAKE_BUILD_TYPE}" UPCASED_CMAKE_BUILD_TYPE)
  if(UPCASED_CMAKE_BUILD_TYPE IN_LIST VALID_BUILD_TYPES)
    # Fetch google benchmark    
    FetchContent_Declare(
      googlebenchmark
      GIT_REPOSITORY https://github.com/google/benchmark.git
      GIT_TAG v1.8.2
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
    set(ENABLE_BENCHMARKING OFF)
    set(ENABLE_BENCHMARKING_MEM_REPORT OFF)
    message(
      WARNING
      "The benchmarks and their depenedencies can be built only if the build is configured "
      "with CMAKE_BUILD_TYPE set to Release or RelWithDebInfo. "
      "The current build is configured with CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}."
      "All benchmarking configuration options are ignored."
    )
  endif()
endif()


if(${USE_LIBFMT})
  set(LIBFMT_VERSION 10.0.0)

  message(
    STATUS 
    "${CMAKE_CXX_COMPILER_ID} v.${CMAKE_CXX_COMPILER_VERSION} does not support std::format, "
    "libfmt v.${LIBFMT_VERSION} will be used instead."
  )

  
  FetchContent_Declare(
    fmt
    GIT_REPOSITORY https://github.com/fmtlib/fmt.git
    GIT_TAG ${LIBFMT_VERSION}
  )

  FetchContent_GetProperties(fmt)

  if(NOT fmt_POPULATED)
    FetchContent_Populate(fmt)
    add_subdirectory(${fmt_SOURCE_DIR} ${fmt_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif()

  endif()
  
# Eigen
# Note: v3.4.0 seems to have a problem detecting c++11 when MSVC is used, so master will be used here.
FetchContent_Declare(
  Eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG        master
)

FetchContent_GetProperties(Eigen)

set(BUILD_TESTING OFF)

if(NOT eigen_POPULATED)
  FetchContent_Populate(Eigen)
  add_subdirectory(${eigen_SOURCE_DIR} ${eigen_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()