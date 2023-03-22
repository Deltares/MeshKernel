# target name
string(CONCAT TARGET_NAME ${PROJECT_NAME} "UnitTests")

# executable
add_executable(${TARGET_NAME})

# source directories
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# include directories
set(INC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

# list of target sources
set(
  SRC_LIST
  ${SRC_DIR}/test_new_delete.cpp
  ${SRC_DIR}/test_containers.cpp
  ${SRC_DIR}/test_smart_pointers.cpp
  ${SRC_DIR}/test_multi_threading.cpp
)

# list of tagret headers
set(
  INC_LIST
  ${INC_DIR}/test_config.hpp
  ${INC_DIR}/test_types.hpp
)

# add sources to target
target_sources(
  ${TARGET_NAME}
  PRIVATE
    ${SRC_LIST}
  PRIVATE
    FILE_SET HEADERS
      BASE_DIRS ${INC_DIR}
      FILES ${INC_LIST}
)

# specify libraries to be linked to the executable
target_link_libraries(
  ${TARGET_NAME}
  PRIVATE
    gtest_main
    benchmark::benchmark
    MeshKernelBenchmarkMemoryManager
)

# group the sources and headers
source_group("Source Files" FILES ${SRC_LIST})
source_group("Header Files" FILES ${INC_LIST})

# add the test to ctest
add_test(
  NAME ${TARGET_NAME}
  COMMAND ${TARGET_NAME}
)
