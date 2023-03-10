# Download google test
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG release-1.11.0)
  
FetchContent_GetProperties(googletest)
if(NOT googletest_POPULATED)
  FetchContent_Populate(googletest)
  add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

set(TEST_FOLDER "${CMAKE_CURRENT_SOURCE_DIR}")
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/utils/test_utils/include/TestUtils/Definitions.hpp.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/utils/test_utils/include/TestUtils/Definitions.hpp"
)

# NOTE: ENABLE_BENCHMARKING configuration option is set to ON for the time being.
#       Therefore, the unit and performance executables will be always built using the 
#       current build server configuration.
# TODO:
#       (1) Have a custom build configuration that sets the configuration option
#           -DENABLE_BENCHMARKING=ON. This allows max automation, especially if
#           the benchmark builds are to be scheduled, say on a weekly basis.
#           Otherwise check if TC allows building with parameters to enable
#           the build and execution of the benchamrking executables on demand.
#       (2) Make the default value of ENABLE_BENCHMARKING OFF in this script.
#       (3) Decide if this option should be move to the root CMakeLists.txt
#       (4) Change the mode of the message in the else body of if(ENABLE_BENCHMARKING)
#           from AUTHOR_WARNING (Warning (dev), continue processing)
#           to FATAL_ERROR (stop processing and generation)
#       (5) Document this in README.md (how to manually configure a build with benchmarking)
#       (6) Remove this commnent.
option(
  ENABLE_BENCHMARKING
  "Enable building and running the benchmark unit and performance executables."
  ON
)

add_subdirectory(utils)
add_subdirectory(unit)
add_subdirectory(api)
add_subdirectory(benchmark)
