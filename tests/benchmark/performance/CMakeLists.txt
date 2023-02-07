# target name
set(target_name MeshKernelBenchmarkPerformanceTests)

# executable
add_executable(${target_name})

# list of target sources
set(
  TARGET_SRC_LIST
  src/main.cpp
  src/dummy_memory_test.cpp
  src/dummy_test.cpp
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
    ${CUSTOM_MEMORY_MANAGER_SRC_LIST}
  PRIVATE
    FILE_SET HEADERS
      BASE_DIRS
        include
        ${CUSTOM_MEMORY_MANAGER_INCLUDE_DIR}
      FILES
        ${TARGET_HEADER_LIST}
        ${CUSTOM_MEMORY_MANAGER_HEADER_LIST}
)

# specify libraries to be linked to the executable
target_link_libraries(
  ${target_name}
  PRIVATE
  MeshKernel
  UtilsStatic
  benchmark::benchmark
)

# group the sources and headers:
# target
source_group("Source Files" FILES ${TARGET_SRC_LIST})
source_group("Header Files" FILES ${TARGET_HEADER_LIST})
# custom meory manager
source_group("Source Files/Custom Memory Manager" FILES ${CUSTOM_MEMORY_MANAGER_SRC_LIST})
source_group("Header Files/Custom Memory Manager" FILES ${CUSTOM_MEMORY_MANAGER_HEADER_LIST})

# add the test to ctest
add_test(
  NAME ${target_name}
  COMMAND ${target_name}
)

