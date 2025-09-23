# zlib
if(WIN32)
  # Set the path to the shared library
  if(DEFINED ENV{ZLIB_ROOT})
    message(VERBOSE "ZLIB_ROOT is an env var")
    set(ZLIB_LIBRARY_PATH "$ENV{ZLIB_ROOT}")
  elseif(DEFINED ZLIB_ROOT)
    message(VERBOSE "ZLIB_ROOT is a config option")
    set(ZLIB_LIBRARY_PATH "${ZLIB_ROOT}")
  else()
    message(FATAL_ERROR "ZLIB_ROOT is undefined")
  endif()
  set(ZLIB_LIBRARY "${ZLIB_LIBRARY_PATH}/bin/zlib.dll")
endif()
find_package(ZLIB REQUIRED)
if (ZLIB_FOUND)
  message(STATUS "Found ZLIB ${ZLIB_VERSION}")
else()
  message(FATAL_ERROR "Could not find ZLIB")
endif()

check_runtime_dependency(ZLIB::ZLIB UNKNOWN_LIBRARY)

set(
  THIRD_PARTY_RUNTIME_DEPS
  $<TARGET_FILE:ZLIB::ZLIB> # not an element of TARGET_RUNTIME_DLLS, helaas speculaas
  CACHE STRING "Third-party runtime dependencies" FORCE
)

# cache all runtime dependencies
set(
  ALL_RUNTIME_DEPS
  ${THIRD_PARTY_RUNTIME_DEPS}
  $<TARGET_FILE:UGridCSharpWrapper>
  CACHE STRING "All runtime dependencies" FORCE
)
