# project name
project(
  MeshKernelBenchmarkMemorySystemQuery
  VERSION ${CMAKE_PROJECT_VERSION}
  DESCRIPTION "MeshKernel benchmark memory system query"
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
  ${SRC_DIR}/memory_system_query.cpp
)

# list of tagret headers
set(
  INC_LIST
  ${INC_DIR}/memory_system_query.hpp
  ${INC_DIR}/memory_monitor.hpp
  ${INC_DIR}/memory_monitor_linux.hpp
  ${INC_DIR}/memory_monitor_windows.hpp
  ${INC_DIR}/memory_monitor_factory.hpp
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
