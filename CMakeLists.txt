project( GridGeom )
cmake_minimum_required(VERSION 3.1)

FIND_PACKAGE( Boost )
FIND_PACKAGE( Cgal REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

FILE( GLOB all_headers src/*.h )
source_group(include FILES ${all_headers})

FILE( GLOB all_srcs src/*.cpp)
source_group(source FILES ${all_srcs})

add_executable( GridGeom
  ${all_srcs}
  ${all_headers}
)

target_link_libraries(GridGeom LINK_PUBLIC ${Boost_LIBRARIES} CGAL::CGAL)

# Google tests 
# Download and unpack googletest at configure time
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

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

# Add tests
add_executable( GridGeomTest src/tests/GridGeomTest.cpp)
target_link_libraries(GridGeomTest LINK_PUBLIC gtest_main ${Boost_LIBRARIES} CGAL::CGAL)
add_test( NAME GridGeomTest COMMAND GridGeomTest )
