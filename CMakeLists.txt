project( gridgeomStateful_dll )
cmake_minimum_required(version 2.8)
set_property(global property use_folders on)

# Find Boost C++ libraries
find_package( Boost )
include_directories( ${Boost_INCLUDE_DIR} )
 
if(win32)
    # Use parallel compilation
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /EHsc /MP")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /EHsc /MP")
else()
    set(CMAKE_CXX_FLAGS_RELEASE "-g -O2" cache string "List of C++ compiler flags for a Release build")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -g") 
endif()

# third party libraries
add_subdirectory(thirdParty)

# add all headers and cpp files
FILE( glob all_cpp_headers src/*.hpp )
source_group(include FILES ${all_cpp_headers})
FILE( glob all_cpp_srcs src/*.cpp)
list(remove_item all_cpp_srcs "src/main.cpp")
source_group(source FILES ${all_cpp_srcs})
set(header_version_file "src/gridgeomStateful_dll_version.h")
set(rc_version_file "src/gridgeomStateful_dll_version.rc")

# add gridgeomStateful_dll stateful library
add_library( gridgeomStateful_dll SHARED ${all_cpp_headers}  ${all_cpp_srcs} ${header_version_file} ${rc_version_file})

# add triangle link target
target_link_libraries(gridgeomStateful_dll link_public ${Boost_LIBRARIES} triangle)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    # add OpenMP link target
    target_link_libraries(gridgeomStateful_dll PUBLIC OpenMP::OpenMP_CXX )
endif()

# Add svn versioning
if(win32)
    set(update_version_script "${CMAKE_CURRENT_SOURCE_DIR}/thirdParty/scripts_lgpl/win32/update_version.cmd")
endif(win32)

set(ini_version_file "${CMAKE_CURRENT_SOURCE_DIR}/src/version_number.ini")
set(header_version_file_path "${CMAKE_CURRENT_SOURCE_DIR}/${header_version_file}")
set(rc_version_file_path "${CMAKE_CURRENT_SOURCE_DIR}/${rc_version_file}")
set(version_directory "${CMAKE_CURRENT_SOURCE_DIR}/src")

add_custom_command( target gridgeomStateful_dll
                    pre_build
                    comand  call ${update_version_script} ${header_version_file_path} ${version_directory} ${ini_version_file})

add_custom_command( target gridgeomStateful_dll
                    pre_build
                    comand  call ${update_version_script} ${rc_version_file_path} ${version_directory} ${ini_version_file})

# Add all unit tests tests
# Download and compile GTests
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(comand ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
                result_variable result
                working_directory ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(fatal_error "CMake step for googletest failed: ${result}")
endif()
execute_process(comand ${CMAKE_COMMAND} --build .
                result_variable result
                working_directory ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(fatal_error "Build step for googletest failed: ${result}")
endif()
if (MSVC)
   set(gtest_force_shared_crt on cache BOOL "" FORCE)
endif()
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                 exclude_from_all)

FILE( glob all_test_srcs src/tests/*.cpp)
source_group(source FILES ${all_test_srcs})

add_executable( GridGeomTest ${all_test_srcs})
target_link_libraries(GridGeomTest link_public ${Boost_LIBRARIES} gtest_main triangle )
add_test( NAME GridGeomTest comand GridGeomTest )