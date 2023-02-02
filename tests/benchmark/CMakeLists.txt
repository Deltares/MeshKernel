if(CMAKE_BUILD_TYPE STREQUAL "Debug")
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

  # add unit and performance tests
  add_subdirectory(unit)
  add_subdirectory(performance)
endif()
