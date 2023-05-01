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

