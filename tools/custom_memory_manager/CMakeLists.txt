# project name
project(
  MeshKernelBenchmarkMemoryManager
  VERSION 1.0.0
  DESCRIPTION "MeshKernel benchamrk memory manager"
  LANGUAGES CXX
)

# target name
set(TARGET_NAME ${PROJECT_NAME})

# executable
add_library(${TARGET_NAME} STATIC)

# source directories
set(SRC_DIR ${PROJECT_SOURCE_DIR}/src)

# include directories
set(INC_DIR ${PROJECT_SOURCE_DIR}/include)

# list of target sources
set(
  SRC_LIST
  ${SRC_DIR}/custom_memory_manager.cpp
  ${SRC_DIR}/new_delete_global_replacements.cpp
)

# list of tagret headers
set(
  INC_LIST
  ${INC_DIR}/custom_memory_manager.hpp
  ${INC_DIR}/platform.hpp
)

# add sources to target
target_sources(
  ${TARGET_NAME}
  PRIVATE
    ${SRC_LIST}
  PUBLIC
    FILE_SET HEADERS
      BASE_DIRS ${INC_DIR}
      FILES ${INC_LIST}
)

# specify libraries to be linked to the executable
target_link_libraries(
  ${TARGET_NAME}
  PRIVATE
    benchmark::benchmark
)

# group the sources and headers
source_group("Source Files" FILES ${SRC_LIST})
source_group("Header Files" FILES ${INC_LIST})

# add unit tests subdirectory
if(ENABLE_UNIT_TESTING)
  add_subdirectory(tests)
endif()
