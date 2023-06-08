# target name
set(target_name Triangle)

# library
add_library(${target_name} STATIC)

# list of target sources
set(
  TARGET_SRC_LIST
  src/triangle.c
  src/TriangleFacade.c
)

# list of target headers
set(
  TARGET_HEADER_LIST
  src/triangle.h
)

# set target sources
target_sources(
  ${target_name}
  PRIVATE
    ${TARGET_SRC_LIST}
  PRIVATE
    FILE_SET HEADERS
      BASE_DIRS src
      FILES ${TARGET_HEADER_LIST}
)

# target definitions
target_compile_definitions(
  ${target_name}
  PRIVATE
    -DANSI_DECLARATORS
    -DTRILIBRARY
    -DNO_TIMER
)

# platform-specific compiler options and definitions for disabling warnings
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  target_compile_options(${target_name} PRIVATE "-w")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  target_compile_options(${target_name} PRIVATE "/w")
endif()

 # group the sources and headers in IDE project generation
 source_group("Source Files" FILES ${TARGET_SRC_LIST})
 source_group("Header Files" FILES ${TARGET_HEADER_LIST})