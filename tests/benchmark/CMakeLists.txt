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

    # Syncs the python tools of the fetched google benchmark with the repository
    # Line endings are ignored
    # This is not good practice because CMake is a build tool, it should not modify the source
    function (SyncPyScriptsFromFetchedContent)
      set(DO_SYNC FALSE)
      set(SOURCE_DIR ${googlebenchmark_SOURCE_DIR}/tools)
      file(GLOB_RECURSE SOURCE_LIST ${SOURCE_DIR}/*.py)
      list(SORT SOURCE_LIST COMPARE STRING)
      set(DESTINATION_DIR ${CMAKE_CURRENT_SOURCE_DIR}/performance/scripts/googlebechmark/tools)
      file(GLOB_RECURSE DESTINATION_LIST ${DESTINATION_DIR}/*.py)
      list(LENGTH DESTINATION_LIST DESTINATION_LIST_LENGTH)
      if( ${DESTINATION_LIST_LENGTH} EQUAL 0)
       # destination is empty, prepare to sync
       set(DO_SYNC TRUE)
      else()
        list(SORT DESTINATION_LIST COMPARE STRING)
        set(COUNTER 0)
        foreach(SOURCE_PY_SCRIPT IN LISTS SOURCE_LIST)
          list(GET DESTINATION_LIST ${COUNTER} DESTINATION_PY_SCRIPT)
          # The process below does not compare timestamps, meaning that the files will remain
          # identical even if the build dir is deleted and the content is refetched in the next build
          # (the fetched content will have a time stamp newer than the one of the file in the destination)
          execute_process(
            COMMAND ${CMAKE_COMMAND} -E compare_files --ignore-eol
            ${SOURCE_PY_SCRIPT}
            ${DESTINATION_PY_SCRIPT}
            RESULT_VARIABLE COMPARE_FILES_RETURN_VALUE
          )
          # Use brute force: if anything is incoherent with a source, 
          # assume all other source files need to be updated
          # remove the destination directory all together and prepare to sync
          if(NOT(COMPARE_FILES_RETURN_VALUE EQUAL 0))
            file(REMOVE_RECURSE ${DESTINATION_DIR})
            set(DO_SYNC TRUE)
            break()
          endif()
          # increment the counter
          math(EXPR COUNTER "${COUNTER} + 1")
        endforeach()
      endif()

      if(DO_SYNC)
        # copy the py scripts from the source dir and remove the write permissions
        # note: do not remove the trailing / from the source directory after COPY
        file(
          COPY ${SOURCE_DIR}/
          DESTINATION ${DESTINATION_DIR}
          NO_SOURCE_PERMISSIONS
          FILE_PERMISSIONS WORLD_READ WORLD_EXECUTE
          FILES_MATCHING
            PATTERN "*.py"
            PATTERN "Inputs" EXCLUDE
        )
      endif()
    endfunction()

    SyncPyScriptsFromFetchedContent()

  else()
    message(
      AUTHOR_WARNING
      "The benchmark unit and performance tests can be built only if the "
      "build is configured with CMAKE_BUILD_TYPE set to Release or RelWithDebInfo. "
      "The current build is configured with CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}."
    )
  endif()
endif()
