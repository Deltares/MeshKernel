# NOTE: ENABLE_BENCHMARKING configuration option is set to ON for the time being.
#       Therefore, the unit and performance executables will be always built using the 
#       current build server configuration.
# TODO:
#       (1) Have a custom build configuration that sets the configuration option
#           -DENABLE_BENCHMARKING=ON. This allows max automation, especially if
#           the benchmark builds are to be scheduled, say on a weekly basis.
#           Otherwise check if TC allows building with parameters to enable
#           the build and execution of the benchamrking executables on demand.
#       (2) Make the default value of ENABLE_BENCHMARKING OFF in this script.
#       (3) Decide if this option should be move to the root CMakeLists.txt
#       (4) Change the mode of the message in the else body of if(ENABLE_BENCHMARKING)
#           from AUTHOR_WARNING (Warning (dev), continue processing)
#           to FATAL_ERROR (stop processing and generation)
#       (5) Document this in README.md (how to manually configure a build with benchmarking)
#       (6) Remove this commnent.
option(
  ENABLE_BENCHMARKING
  "Enable building and running the benchmark unit and performance executables."
  ON
)

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

    # Add custom_memory_manager
    add_subdirectory(custom_memory_manager)

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
