# target name
set(target_name MeshKernelBenchmarkMemoryManager)

# executable
add_library(${target_name} STATIC)

# list of target sources
set(
  TARGET_SRC_LIST
  src/custom_memory_manager.cpp
  src/new_delete_global_replacements.cpp
)

# list of tagret headers
set(
  TARGET_HEADER_LIST
  include/custom_memory_manager.hpp
  include/platform.hpp
)

# add sources to target
target_sources(
  ${target_name}
  PRIVATE
    ${TARGET_SRC_LIST}
  PUBLIC
    FILE_SET HEADERS
      BASE_DIRS include
      FILES ${TARGET_HEADER_LIST}
)

# specify libraries to be linked to the executable
target_link_libraries(
  ${target_name}
  PRIVATE
  benchmark::benchmark
)

# group the sources and headers
source_group("Source Files" FILES ${TARGET_SRC_LIST})
source_group("Header Files" FILES ${TARGET_HEADER_LIST})
