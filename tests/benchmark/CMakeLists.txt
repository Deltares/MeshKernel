string(TOUPPER "${CMAKE_BUILD_TYPE}" UPCASED_CMAKE_BUILD_TYPE)

if(UPCASED_CMAKE_BUILD_TYPE STREQUAL "RELEASE" 
  OR UPCASED_CMAKE_BUILD_TYPE STREQUAL "RELWITHDEBINFO"
)
  # Fetch google benchmark dependency
  FetchContent_Declare(
    googlebenchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.7.1
  )

  FetchContent_GetProperties(benchmark)
  if(NOT googlebenchmark_POPULATED)
    FetchContent_Populate(googlebenchmark)
    add_subdirectory(
      ${googlebenchmark_SOURCE_DIR}
      ${googlebenchmark_BINARY_DIR}
      EXCLUDE_FROM_ALL
    )
  endif()

  # Prevent the google benchmark's own tests from running
  set(BENCHMARK_ENABLE_TESTING OFF)

  # set the root, src and include paths to custom memory manager
  set(CUSTOM_MEMORY_MANAGER_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/custom_memory_manager)
  set(CUSTOM_MEMORY_MANAGER_SRC_DIR ${CUSTOM_MEMORY_MANAGER_ROOT_DIR}/src)
  set(CUSTOM_MEMORY_MANAGER_INCLUDE_DIR ${CUSTOM_MEMORY_MANAGER_ROOT_DIR}/include)

  # set list of custom memory manager sources
  set(CUSTOM_MEMORY_MANAGER_SRC_LIST
    ${CUSTOM_MEMORY_MANAGER_SRC_DIR}/custom_memory_management.cpp
    ${CUSTOM_MEMORY_MANAGER_SRC_DIR}/new_delete_global_replacements.cpp
  )

  # set list of custom memory manager headers
  set(CUSTOM_MEMORY_MANAGER_HEADER_LIST
    ${CUSTOM_MEMORY_MANAGER_INCLUDE_DIR}/custom_memory_manager.hpp
    ${CUSTOM_MEMORY_MANAGER_INCLUDE_DIR}/memory_management.hpp
    ${CUSTOM_MEMORY_MANAGER_INCLUDE_DIR}/platform.hpp
  )

  # Add unit tests
  add_subdirectory(unit)

  # Add unit performance tests
  add_subdirectory(performance)
endif()
