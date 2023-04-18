# Determine version suffix (current build number) from environment variable
if(DEFINED ENV{VERSION_SUFFIX})
  set(VERSION_SUFFIX $ENV{VERSION_SUFFIX})
else()
  set(VERSION_SUFFIX 0)
endif()

set(VERSION_STRING ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}.${VERSION_SUFFIX})
message(STATUS "Version is set to ${VERSION_STRING}")

# Configure version file
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Version.hpp.in"
               "${CMAKE_CURRENT_SOURCE_DIR}/Version.hpp" @ONLY)
