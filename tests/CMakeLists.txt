# Download google test
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG release-1.10.0)
FetchContent_MakeAvailable(googletest)

set(TEST_FOLDER "${CMAKE_CURRENT_SOURCE_DIR}")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/utils/include/TestUtils/Definitions.hpp.in" "${CMAKE_CURRENT_SOURCE_DIR}/utils/include/TestUtils/Definitions.hpp")

add_subdirectory(utils)
add_subdirectory(unit)
add_subdirectory(api)
