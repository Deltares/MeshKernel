# Set the c++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Disable compiler specific extensions
set(CMAKE_CXX_EXTENSIONS OFF)

# Create position-independent executables and shared libraries
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Add compiler-specific options and definitions per supported platform
if (UNIX)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_compile_options("-fvisibility=hidden;-Werror;-Wall;-Wextra;-pedantic;-Wno-unused-function")
    if(APPLE AND (CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "arm64"))
      # CMake automatically sets -Xarch_arm64 (for clang) but gcc doesn't support it
      unset(_CMAKE_APPLE_ARCHS_DEFAULT)
      # Be lenient on macos with arm64 toolchain to prevent Eigen -Werror=deprecated-enum-enum-conversion error
      add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-Wno-deprecated-enum-enum-conversion>)
      # Suppress notes related to ABI changes
      add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-Wno-psabi>)
    endif()
    add_compile_options("$<$<CONFIG:RELEASE>:-O2>")
    add_compile_options("$<$<CONFIG:DEBUG>:-g>")
  else()
    message(FATAL_ERROR "Unsupported compiler. Only GNU is supported under Linux. Found ${CMAKE_CXX_COMPILER_ID}.")
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
    message(FATAL_ERROR "Unsupported platform. Only Linux and Windows are supported.")
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
