# Determine version suffix from environment variable
if(DEFINED ENV{VERSION_SUFFIX})
  set(VERSION_SUFFIX $ENV{VERSION_SUFFIX})
else()
  set(VERSION_SUFFIX .0)
endif()

# Configure version file
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Version.hpp.in"
               "${CMAKE_CURRENT_SOURCE_DIR}/Version.hpp" @ONLY)
