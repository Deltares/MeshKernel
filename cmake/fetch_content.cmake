include(FetchContent)

if(ENABLE_UNIT_TESTING)
  # Fetch google test
  FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.17.0
  )

  FetchContent_MakeAvailable(googletest)

  # Note this needs to be done in the main CMakeLists since it calls
  # enable_testing, which must be in the main CMakeLists.
  include(CTest)
endif()

if(ENABLE_BENCHMARKING)
  # check if it's a release or release with debug information build
  set(VALID_BUILD_TYPES "RELEASE" "RELWITHDEBINFO" "DEBUG")
  # upcase the build type to make the option case-insensitive
  string(TOUPPER "${CMAKE_BUILD_TYPE}" UPCASED_CMAKE_BUILD_TYPE)
  if(UPCASED_CMAKE_BUILD_TYPE IN_LIST VALID_BUILD_TYPES)
    # Fetch google benchmark
    FetchContent_Declare(
      googlebenchmark
      GIT_REPOSITORY https://github.com/google/benchmark.git
      GIT_TAG v1.8.2
    )
    # Prevent the google benchmark's own tests from running
    set(BENCHMARK_ENABLE_TESTING OFF)
    # Modern replacement for deprecated GetProperties/Populate sequence
    FetchContent_MakeAvailable(googlebenchmark)
  else()
    # disable all benchmarking options
    set(ENABLE_BENCHMARKING OFF)
    set(ENABLE_BENCHMARKING_MEM_REPORT OFF)
    message(
      WARNING
      "The benchmarks and their dependencies can be built only if the build is configured "
      "with CMAKE_BUILD_TYPE set to Release or RelWithDebInfo. "
      "The current build is configured with CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}."
      "All benchmarking configuration options are ignored."
    )
  endif()
endif()


if(${USE_LIBFMT})
  set(LIBFMT_VERSION 11.0.0)

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

  FetchContent_MakeAvailable(fmt)

  endif()

# Eigen
# Note: v3.4.0 seems to have a problem detecting c++11 when MSVC is used, so the head master will be used here.
FetchContent_Declare(
  Eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG        21cd3fe20990a5ac1d683806f605110962aac3f1
)

FetchContent_MakeAvailable(Eigen)

set(BUILD_TESTING OFF)

