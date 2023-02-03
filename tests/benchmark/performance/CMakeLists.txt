# Tests need to be added as executables first
file(GLOB BENCHMARK_PERFORMANCE_TEST_LIST CONFIGURE_DEPENDS
  "src/*.cpp"
  "../custom_memory_manager/src/*.cpp"
)

set(target_name MeshKernelBenchmarkPerformanceTests)

add_executable(${target_name} ${BENCHMARK_PERFORMANCE_TEST_LIST})

target_include_directories(
  ${target_name}
  PRIVATE
  include
  ../custom_memory_manager/include
)

# Should be linked to the main library, as well as the google test library
target_link_libraries(
  ${target_name}
  PRIVATE
  MeshKernel
  UtilsStatic
  benchmark::benchmark
)

if (UNIX)
   target_link_libraries(${target_name} PUBLIC ${CMAKE_DL_LIBS})
endif (UNIX)

add_test(
  NAME ${target_name}
  COMMAND ${target_name}
)

