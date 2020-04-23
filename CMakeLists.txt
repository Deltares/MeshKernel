project( gridgeomStateful_dll )
cmake_minimum_required(VERSION 2.8)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

#Boost
FIND_PACKAGE( Boost )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )
 
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /EHsc") 

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

# third party directory
add_subdirectory(thirdParty)

#add all headers and cpp files
FILE( GLOB all_cpp_headers src/*.hpp )
source_group(include FILES ${all_cpp_headers})

FILE( GLOB all_cpp_srcs src/*.cpp)
source_group(source FILES ${all_cpp_srcs})

# definitation for triangle
add_library( gridgeomStateful_dll SHARED ${all_cpp_headers}  ${all_cpp_srcs})
target_link_libraries(gridgeomStateful_dll LINK_PUBLIC ${Boost_LIBRARIES} triangle)

set(CMAKE_CXX_FLAGS_RELEASE "-g -O2" CACHE STRING "List of C++ compiler flags for a Release build")

# Add tests
FILE( GLOB all_test_srcs src/tests/*.cpp)
source_group(source FILES ${all_test_srcs})

add_executable( GridGeomTest ${all_test_srcs})
target_link_libraries(GridGeomTest LINK_PUBLIC ${Boost_LIBRARIES} gtest_main triangle )
add_test( NAME GridGeomTest COMMAND GridGeomTest )

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(gridgeomStateful_dll PUBLIC OpenMP::OpenMP_CXX )
endif()

if(WIN32)
   # Add main executable
   add_executable( GridGeomMainTest src/main.cpp )
   target_link_libraries(GridGeomMainTest LINK_PUBLIC gridgeomStateful_dll ${Boost_LIBRARIES})
   if(OpenMP_CXX_FOUND)
      target_link_libraries(GridGeomMainTest PUBLIC OpenMP::OpenMP_CXX)
   endif()
endif()