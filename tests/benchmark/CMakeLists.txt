# Tests need to be added as executables first
file(GLOB BENCHMARK_TEST_LIST CONFIGURE_DEPENDS "src/*.cpp")

set(target_name MeshKernelBenchmarkTests)

add_executable(${target_name} ${BENCHMARK_TEST_LIST})

target_include_directories(${target_name} PRIVATE include)

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

# Should be linked to the main library, as well as the google test library
target_link_libraries(
  ${target_name}
  PRIVATE
  MeshKernel
  UtilsStatic
  gtest
  benchmark::benchmark
  benchmark::benchmark_main
)

if (UNIX)
   target_link_libraries(${target_name} PUBLIC ${CMAKE_DL_LIBS})
endif (UNIX)

# Make sure that coverage information is produced when using gcc
if(PRODUCE_CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_target_properties( ${target_name} PROPERTIES COMPILE_FLAGS "-fprofile-arcs -ftest-coverage -static-libgcc -static-libstdc++")
  set_target_properties(${target_name} PROPERTIES LINK_FLAGS "-ldl -lgcov --coverage")
endif()

# If you register a test, then ctest and make test will run it. You can also run
# examples and check the output, as well. Command can be a target.
add_test(NAME ${target_name}  COMMAND ${target_name} )

# On Windows we bundle the dlls and then copy them to the target file directory
if(WIN32)
  file(
    GLOB
    NETCDF_FILES_LIST
    CONFIGURE_DEPENDS
    "${NETCDF_BIN_PATH}/${CMAKE_SHARED_LIBRARY_PREFIX}*${CMAKE_SHARED_LIBRARY_SUFFIX}"
  )
  foreach(CurrentNetcdfLibFile IN LISTS NETCDF_FILES_LIST)
    add_custom_command(
      TARGET ${target_name} 
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CurrentNetcdfLibFile}
              "$<TARGET_FILE_DIR:${target_name}>"
      COMMENT "Copying netcdf lib file: ${CurrentNetcdfLibFile}")
  endforeach()
endif()