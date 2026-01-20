# Set the c++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Disable compiler specific extensions
set(CMAKE_CXX_EXTENSIONS OFF)

# Create position-independent executables and shared libraries
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Add compiler-specific options and definitions per supported platform
if(APPLE)
  # macOS: support AppleClang / Clang (system toolchain). GCC may also be supported but Clang is preferred.
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(STATUS "Configuring build for macOS with ${CMAKE_CXX_COMPILER_ID} (${CMAKE_CXX_COMPILER_VERSION}).")
    # Common warning and visibility flags
    add_compile_options("-fvisibility=hidden;-Wall;-Wextra;-pedantic;-Werror;-Wno-unused-function")
    # Conditionally suppress unused parameter warnings (used for platform-specific compiler issues)
    if(SUPPRESS_UNUSED_PARAMETER_WARNING)
      add_compile_options("-Wno-unused-parameter")
      message(STATUS "Suppressing -Wno-unused-parameter warnings for macOS Clang build on macOS15-intel")
    endif()
    # Be lenient for Eigen (deprecated enum conversion) on all macOS Clang builds
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-Wno-deprecated-enum-enum-conversion>)
    # Optimization / debug flags
    add_compile_options("$<$<CONFIG:RELEASE>:-O2>")
    add_compile_options("$<$<CONFIG:DEBUG>:-g>")
  else()
    message(FATAL_ERROR "Unsupported compiler on macOS. Supported: AppleClang/Clang. Found ${CMAKE_CXX_COMPILER_ID}.")
  endif()
elseif(UNIX)
  # Linux: require GNU
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(STATUS "Configuring build for Linux with GNU ${CMAKE_CXX_COMPILER_VERSION}.")
    add_compile_options("-fvisibility=hidden;-Werror;-Wall;-Wextra;-pedantic;-Wno-unused-function")
    add_compile_options("$<$<CONFIG:RELEASE>:-O2>")
    add_compile_options("$<$<CONFIG:DEBUG>:-g>")
  else()
    message(FATAL_ERROR "Unsupported compiler on Linux. Only GNU is supported. Found ${CMAKE_CXX_COMPILER_ID}.")
  endif()
elseif(WIN32)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    add_compile_options("/EHsc;/MP;/W3;/WX")
    add_compile_options("$<$<CONFIG:RELEASE>:/O2>")
    add_compile_options("$<$<CONFIG:DEBUG>:/Od;/DEBUG>")
    add_compile_definitions("_USE_MATH_DEFINES")
    add_compile_definitions("_CRT_SECURE_NO_WARNINGS")
  else()
    message(FATAL_ERROR "Unsupported compiler. Only MSVC is supported under Windows. Found ${CMAKE_CXX_COMPILER_ID}.")
  endif()
else()
    message(FATAL_ERROR "Unsupported platform. Only macOS, Linux and Windows are supported.")
endif()

# CMAKE_SOURCE_DIR is passed to the src in order to strip it out of the path of srcs where exceptions may occur
add_compile_definitions(CMAKE_SRC_DIR=${CMAKE_SOURCE_DIR})

# Show the source location in the exception message?
add_compile_definitions(HAVE_SRC_LOC_IN_ERR_MSGS=$<BOOL:${HAVE_SRC_LOC_IN_ERR_MSGS}>)

# format library: from the standard lib or third-party?
# When supported, std::format is preferred. Otherwise, fmtlib should be used.
set(USE_LIBFMT 0)
if(
  (CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
    AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 13.1)
  OR
  (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC"
    AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 16.11.14)
)
  set(USE_LIBFMT 1)
endif()
add_compile_definitions(USE_LIBFMT=${USE_LIBFMT})
