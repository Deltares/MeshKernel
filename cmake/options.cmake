include(CMakeDependentOption)

# unit testing option
option(
  ENABLE_UNIT_TESTING
  "Enables building the unit test executables"
  ON
)

# benchmarking option
option(
  ENABLE_BENCHMARKING
  "Enables building the benchmark executable."
  OFF
)

# memory usgae reporting option
cmake_dependent_option(
  ENABLE_BENCHMARKING_MEM_REPORT
  "When benchmarking is enabled, enables reporting memory usage statistics."
  ON 
  "ENABLE_BENCHMARKING"
  OFF
)

if(ENABLE_BENCHMARKING_MEM_REPORT)
  set(
    MEM_COLLECTION_METHOD "" 
    CACHE STRING 
    "When memory reporting is enabled, specifies how memory metrics are collected."
  )
  if( NOT (MEM_COLLECTION_METHOD STREQUAL "QUERY_SYSTEM" OR MEM_COLLECTION_METHOD STREQUAL "COUNT_BYTES"))
    message(WARNING "Option MEM_COLLECTION_METHOD not set. QUERY_SYSTEM will be used.")
    set(MEM_COLLECTION_METHOD "QUERY_SYSTEM")
  endif()
endif()

# code coverage option
if(LINUX AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  option(
    ENABLE_CODE_COVERAGE
    "Generates code coverage files under GNU compilers"
    OFF
  )
endif()