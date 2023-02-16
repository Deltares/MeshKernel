# target name
set(target_name MeshKernelBenchmarkPerformanceTests)

# executable
add_executable(${target_name})

# list of target sources
set(
  TARGET_SRC_LIST
  src/main.cpp
  src/perf_rtree.cpp
  src/perf_mesh_refinement.cpp
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
  UtilsStatic
  benchmark::benchmark
  MeshKernelBenchmarkMemoryManager
)

if (UNIX)
   target_link_libraries(${target_name} PUBLIC ${CMAKE_DL_LIBS})
endif (UNIX)

# group the sources and headers
source_group("Source Files" FILES ${TARGET_SRC_LIST})
source_group("Header Files" FILES ${TARGET_HEADER_LIST})

# add the test to ctest
add_test(
  NAME ${target_name}
  COMMAND ${target_name}
)

#message("${CMAKE_CURRENT_BINARY_DIR}")

add_custom_command(
  TARGET ${target_name} 
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory 
     ${googlebenchmark_SOURCE_DIR}/tools/
    $<TARGET_FILE_DIR:${target_name}>/tools
  COMMENT "Copying Google Benchmarl python tools"
)
