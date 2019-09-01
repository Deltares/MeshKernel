project( GridGeom )
cmake_minimum_required(VERSION 3.1)

#Boost
FIND_PACKAGE( Boost )

#Vc SIMD classes
#set(Vc_DIR "thirdParty/Vc/installIntel/lib/cmake/Vc")
#FIND_PACKAGE( Vc )

#CGAL geometric libraries
FIND_PACKAGE( Cgal REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

SET( EIGEN3_INCLUDE_DIR thirdParty/eigen )
INCLUDE_DIRECTORIES ( "${EIGEN3_INCLUDE_DIR}" )

FILE( GLOB all_headers src/*.hpp )
source_group(include FILES ${all_headers})

FILE( GLOB all_srcs src/*.cpp)
source_group(source FILES ${all_srcs})

add_executable( GridGeom
  ${all_srcs}
  ${all_headers}
)

if (UNIX)
	set(cxxVersion "-std=c++11")
	MESSAGE(STATUS "CXX version is set to ${cxxVersion}")
	
	#Release
    set(CMAKE_C_FLAGS_RELEASE "-g -O2" CACHE STRING "List of C compiler flags for a Release build")
    set(CMAKE_CXX_FLAGS_RELEASE "-g -O2 ${cxxVersion}" CACHE STRING "List of C++ compiler flags for a Release build")
	MESSAGE(STATUS "CMAKE_CXX_FLAGS_RELEASE is set to ${CMAKE_CXX_FLAGS_RELEASE}")

    #Debug
    set(CMAKE_C_FLAGS_DEBUG"-g -O0 ${preciseFpModel}" CACHE STRING "List of C compiler flags for a Debug build")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0  ${preciseFpModel} ${cxxVersion}" CACHE STRING "List of C++ compiler flags for a Debug build")

endif()

target_link_libraries(GridGeom LINK_PUBLIC) #  ${Boost_LIBRARIES} CGAL::CGAL

# Google tests 
# https://crascit.com/2015/07/25/cmake-gtest/
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
add_executable( GridGeomTest src/tests/GridGeomTest.cpp src/tests/OrthogonalizationTest.cpp)
target_link_libraries(GridGeomTest LINK_PUBLIC gtest_main ) # ${Boost_LIBRARIES} CGAL::CGAL
add_test( NAME GridGeomTest COMMAND GridGeomTest )