# Tests need to be added as executables first
file(GLOB BENCHMARK_UNIT_TEST_LIST CONFIGURE_DEPENDS
  "src/*.cpp"
  "../custom_memory_manager/src/*.cpp"
)

set(target_name MeshKernelBenchmarkUnitTests)

add_executable(${target_name} ${BENCHMARK_UNIT_TEST_LIST})

target_include_directories(
  ${target_name}
  PRIVATE
  include
  ../custom_memory_manager/include
)

target_link_libraries(
  ${target_name}
  PRIVATE
  gtest
  gtest_main
  benchmark::benchmark
)

if (UNIX)
   target_link_libraries(${target_name} PUBLIC ${CMAKE_DL_LIBS})
endif (UNIX)

add_test(
  NAME ${target_name}
  COMMAND ${target_name}
)
