# target name
set(target_name MeshKernelBenchmarkUnitTests)

# executable
add_executable(${target_name})

# list of target sources
set(
  TARGET_SRC_LIST
  src/test_new_delete.cpp
  src/test_containers.cpp
  src/test_smart_pointers.cpp
)

# list of tagret headers
set(
  TARGET_HEADER_LIST
  include/test_config.hpp
  include/test_types.hpp
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
  gtest
  gtest_main
  benchmark::benchmark
  MeshKernelBenchmarkMemoryManager
)

# group the sources and headers
source_group("Source Files" FILES ${TARGET_SRC_LIST})
source_group("Header Files" FILES ${TARGET_HEADER_LIST})

# add the test to ctest
add_test(
  NAME ${target_name}
  COMMAND ${target_name}
)
