# Download google test
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG release-1.10.0)
FetchContent_MakeAvailable(googletest)

file(GLOB UTILS_TEST_LIST CONFIGURE_DEPENDS "utils/*.cpp")

add_subdirectory(utils)
add_subdirectory(unit)
