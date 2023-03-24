# Set the c++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Disable compiler specific extensions
set(CMAKE_CXX_EXTENSIONS OFF)

# Create position-independent executables and shared libraries
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Add compiler-specific options and definitions per supported platform
if (UNIX)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_compile_options("-Werror;-Wall;-Wextra;-pedantic;-Wno-unused-function")
    add_compile_options("$<$<CONFIG:RELEASE>:-O2>")
    add_compile_options("$<$<CONFIG:DEBUG>:-g>")
  else()
    message(FATAL_ERROR "Unsupported compiler. Only GNU is supported under Linux.")
  endif()
elseif(WIN32)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    add_compile_options("/EHsc;/MP;/W3;/WX")
    add_compile_options("$<$<CONFIG:RELEASE>:/O2>")
    add_compile_options("$<$<CONFIG:DEBUG>:/Od;/DEBUG>")
    add_compile_definitions("_USE_MATH_DEFINES")
  else()
    message(FATAL_ERROR "Unsupported compiler. Only MSVC is supported under Windows.")
  endif()
else()
    message(FATAL_ERROR "Unsupported platform. Only Linux and Windows are supported.")
endif()




