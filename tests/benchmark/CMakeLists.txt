if(ENABLE_BENCHMARKING)
  # upcase the build type to make the option case-insensitive
  string(TOUPPER "${CMAKE_BUILD_TYPE}" UPCASED_CMAKE_BUILD_TYPE)

  if(UPCASED_CMAKE_BUILD_TYPE STREQUAL "RELEASE" 
    OR UPCASED_CMAKE_BUILD_TYPE STREQUAL "RELWITHDEBINFO"
  )
    # Fetch google benchmark dependency
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

    # Add unit performance tests
    add_subdirectory(performance)

  else()
    message(
      AUTHOR_WARNING
      "The benchmark unit and performance tests can be built only if the "
      "build is configured with CMAKE_BUILD_TYPE set to Release or RelWithDebInfo. "
      "The current build is configured with CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}."
    )
  endif()
endif()
