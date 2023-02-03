if(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  # Fetch google benchmark dependency
  FetchContent_Declare(
    googlebenchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.7.1
  )

  FetchContent_GetProperties(benchmark)
  if(NOT googlebenchmark_POPULATED)
    FetchContent_Populate(googlebenchmark)
    add_subdirectory(
      ${googlebenchmark_SOURCE_DIR}
      ${googlebenchmark_BINARY_DIR}
      EXCLUDE_FROM_ALL
    )
  endif()

  set(BENCHMARK_ENABLE_TESTING OFF)

  # add unit and performance tests
  add_subdirectory(unit)
  add_subdirectory(performance)
endif()
