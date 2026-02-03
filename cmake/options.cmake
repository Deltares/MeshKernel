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
include(CMakeDependentOption)
cmake_dependent_option(
  ENABLE_BENCHMARKING_MEM_REPORT
  "When benchmarking is enabled, enables reporting memory usage statistics."
  ON
  "ENABLE_BENCHMARKING"
  OFF
)

# code coverage option
if(LINUX AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  option(
    ENABLE_CODE_COVERAGE
    "Generates code coverage files under GNU compilers."
    OFF
  )
endif()

# error messages
option(
  HAVE_SRC_LOC_IN_ERR_MSGS
  "Includes source location information in customized exceptions."
  OFF
)

# suppress unused parameter warnings (sometimes necessary for platform-specific builds)
option(
  SUPPRESS_UNUSED_PARAMETER_WARNING
  "Suppress -Wno-unused-parameter warning flag (used for platform-specific compiler issues)."
  OFF
)
