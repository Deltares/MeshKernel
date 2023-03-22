# target name
string(CONCAT TARGET_NAME ${PROJECT_NAME} "Benchmarks")

# executable
add_executable(${TARGET_NAME})

# source directory
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# list of target sources
set(
  SRC_LIST
  ${SRC_DIR}/main.cpp
  ${SRC_DIR}/perf_rtree.cpp
  ${SRC_DIR}/perf_mesh_refinement.cpp
  ${SRC_DIR}/perf_orthogonalization.cpp
  ${SRC_DIR}/perf_curvilinear_uniform.cpp
)

# add sources to target
target_sources(${TARGET_NAME} PRIVATE ${SRC_LIST})

# specify libraries to be linked to the executable
target_link_libraries(
  ${TARGET_NAME}
  PRIVATE
  MeshKernel
  MeshKernelTestUtils
  benchmark::benchmark
)

if(ENABLE_BENCHMARKING_MEM_REPORT)
  # link memory manager
  target_link_libraries(
    ${TARGET_NAME}
    PRIVATE
    MeshKernelBenchmarkMemoryManager
  )
  # add compiler definition for reporting the memory statistics in src
  target_compile_definitions(
    ${TARGET_NAME}
    PRIVATE
    ENABLE_BENCHMARKING_MEM_REPORT
  )
endif()

# group the sources
source_group("Source Files" FILES ${SRC_LIST})

# add the test to ctest
add_test(
  NAME ${TARGET_NAME}
  COMMAND ${TARGET_NAME}
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
  TARGET ${TARGET_NAME}
  POST_BUILD
  COMMAND ${CMAKE_COMMAND}
    -DSOURCE_DIR=${googlebenchmark_SOURCE_DIR}/tools/
    -DDESTINATION_DIR=$<TARGET_FILE_DIR:${TARGET_NAME}>/scripts/googlebenchmark
    -DPATTERN_EXCLUDE=Inputs
    -P ${CMAKE_CURRENT_BINARY_DIR}/copy_py_scripts.cmake
  COMMAND ${CMAKE_COMMAND}
    -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/../scripts/
    -DDESTINATION_DIR=$<TARGET_FILE_DIR:${TARGET_NAME}>/scripts/meshkernel
    -P ${CMAKE_CURRENT_BINARY_DIR}/copy_py_scripts.cmake
  COMMENT "Copying python scripts..."
)
