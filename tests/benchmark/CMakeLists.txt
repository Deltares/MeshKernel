# this target is built only if benchmarking is enabled

if(ENABLE_BENCHMARKS)
  # target name
  set(target_name MeshKernelBenchmarkTests)

  # executable
  add_executable(${target_name})

  # list of target sources
  set(
    TARGET_SRC_LIST
    src/main.cpp
    src/perf_rtree.cpp
    src/perf_mesh_refinement.cpp
    src/perf_orthogonalization.cpp
    src/perf_curvilinear_uniform.cpp
  )

  # list of tagret headers
  set(
    TARGET_HEADER_LIST
  )

  # add sources to target
  target_sources(
    ${target_name}
    PRIVATE
      ${TARGET_SRC_LIST}
    PRIVATE
      FILE_SET HEADERS
        BASE_DIRS include
        FILES ${TARGET_HEADER_LIST}
  )

  # specify libraries to be linked to the executable
  target_link_libraries(
    ${target_name}
    PRIVATE
    MeshKernel
    MeshKernelTestUtils
    benchmark::benchmark
  )

  if(ENABLE_MEM_REPORT)
    # link memory manager
    target_link_libraries(
      ${target_name}
      PRIVATE
      MeshKernelBenchmarkMemoryManager
    )
    # add compiler definition for reporting the memory statistics in src
    target_compile_definitions(
      ${target_name}
      PRIVATE
      ENABLE_MEM_REPORT
    )
  endif()

  # group the sources and headers
  source_group("Source Files" FILES ${TARGET_SRC_LIST})
  source_group("Header Files" FILES ${TARGET_HEADER_LIST})

  # add the test to ctest
  add_test(
    NAME ${target_name}
    COMMAND ${target_name}
  )

  # Generate script copy_py_scripts.cmake which copies python scripts from a given
  # source dir to a desired destination dir while preserving the structure of the source dir
  file(
    WRITE ${CMAKE_CURRENT_BINARY_DIR}/copy_py_scripts.cmake
    "message(STATUS \"Copying Python scripts from \$\{SOURCE_DIR\} to \$\{DESTINATION_DIR\}\")
    file(
      COPY \$\{SOURCE_DIR\} 
      DESTINATION \$\{DESTINATION_DIR\}
      FILES_MATCHING PATTERN \"*.py\"
    )"
  )
  
  # custom commans that executes copy_py_scripts.cmake to copy python
  # scripts from the project source dir to the current build dir
  add_custom_command(
    TARGET ${target_name}
    POST_BUILD
    COMMAND ${CMAKE_COMMAND}
      -DSOURCE_DIR=${googlebenchmark_SOURCE_DIR}/tools/
      -DDESTINATION_DIR=$<TARGET_FILE_DIR:${target_name}>/scripts/googlebenchmark
      -DPATTERN_EXCLUDE=Inputs
      -P ${CMAKE_CURRENT_BINARY_DIR}/copy_py_scripts.cmake
    COMMAND ${CMAKE_COMMAND}
      -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/../scripts/
      -DDESTINATION_DIR=$<TARGET_FILE_DIR:${target_name}>/scripts/meshkernel
      -P ${CMAKE_CURRENT_BINARY_DIR}/copy_py_scripts.cmake
    COMMENT "Copying python scripts..."
  )
endif()
