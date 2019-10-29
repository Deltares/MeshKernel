project( gridgeomStateful_dll )
cmake_minimum_required(VERSION 2.8)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

#add all headers and cpp files
FILE( GLOB all_headers src/*.hpp )
source_group(include FILES ${all_headers})

FILE( GLOB all_srcs src/*.cpp)
source_group(source FILES ${all_srcs})

add_library( gridgeomStateful_dll STATIC ${all_srcs} ${all_headers})

target_link_libraries(gridgeomStateful_dll LINK_PUBLIC)

#-----------------------------------------------------------------------------------------------------------------#
# Google tests 
# https://crascit.com/2015/07/25/cmake-gtest/
# Download and unpack googletest at configure time
#-----------------------------------------------------------------------------------------------------------------#

configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()
# Prevent overriding the parent project's compiler/linker settings on Windows
if (MSVC)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
endif()

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)
				 
# Add tests
add_executable( GridGeomTest src/tests/MeshTest.cpp src/tests/OrthogonalizationTest.cpp)
target_link_libraries(GridGeomTest LINK_PUBLIC gtest_main )
add_test( NAME GridGeomTest COMMAND GridGeomTest )

if(WIN32)
# Add main executable
add_executable( GridGeomMainTest src/main.cpp )
target_link_libraries(GridGeomMainTest LINK_PUBLIC gridgeomStateful_dll )
endif()