cmake_minimum_required(VERSION 3.14...3.18)
project( MeshKernel )
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

if(WIN32)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /EHsc /MP" )
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /EHsc /MP")
else()
    set(CMAKE_CXX_FLAGS_RELEASE "-g -O2" CACHE STRING "List of C++ compiler flags for a Release build")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -g") 
endif(WIN32)

# Add third party directory
add_subdirectory(thirdParty)

# Boost C++
FIND_PACKAGE( Boost )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

# Add all headers and cpp files
# TODO: Take care of Operations and Constants and TriangulationWrapper
FILE( GLOB all_cpp_headers include/MeshKernel/*.hpp  include/MeshKernel/Operations.cpp include/MeshKernel/Constants.cpp include/MeshKernel/TriangulationWrapper.cpp)
source_group(include FILES ${all_cpp_headers})

FILE( GLOB all_cpp_srcs src/*.cpp)
list(REMOVE_ITEM all_cpp_srcs "src/main.cpp")
source_group(source FILES ${all_cpp_srcs})

# Add MeshKernel library
add_library( MeshKernel SHARED ${all_cpp_headers}  ${all_cpp_srcs} ${header_version_file} ${rc_version_file})

# We need this directory, and users of our library will need it too
target_include_directories(MeshKernel PUBLIC ./include)

# Add target link dependency on triangle
target_link_libraries(MeshKernel LINK_PUBLIC ${Boost_LIBRARIES} triangle)

# Add target link dependency on OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(MeshKernel PUBLIC OpenMP::OpenMP_CXX )
endif()

# Add dll versioning
set(header_version_file "src/MeshKernel_version.h")
set(rc_version_file "src/MeshKernel_version.rc")

if(WIN32)
    set(update_version_script "${CMAKE_CURRENT_SOURCE_DIR}/thirdParty/scripts_lgpl/win32/update_version.cmd")
endif(WIN32)
set(ini_version_file "${CMAKE_CURRENT_SOURCE_DIR}/src/version_number.ini")
set(header_version_file_path "${CMAKE_CURRENT_SOURCE_DIR}/${header_version_file}")
set(rc_version_file_path "${CMAKE_CURRENT_SOURCE_DIR}/${rc_version_file}")
set(version_directory "${CMAKE_CURRENT_SOURCE_DIR}/src")

add_custom_command( TARGET MeshKernel
                    PRE_BUILD
                    COMMAND  call ${update_version_script} ${header_version_file_path} ${version_directory} ${ini_version_file})

add_custom_command( TARGET MeshKernel
                    PRE_BUILD
                    COMMAND  call ${update_version_script} ${rc_version_file_path} ${version_directory} ${ini_version_file})
                    
set(build_config $<CONFIG>) 
set_target_properties( MeshKernel PROPERTIES 
                       LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${build_config}/lib"
                       ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${build_config}/lib"
                       RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${build_config}/lib")
					  
					   
# # Add Google tests
# configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
# execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
#                 RESULT_VARIABLE result
#                 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
# if(result)
#    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
# endif()
# execute_process(COMMAND ${CMAKE_COMMAND} --build .
#                 RESULT_VARIABLE result
#                 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
# if(result)
#   message(FATAL_ERROR "Build step for googletest failed: ${result}")
# endif()
# if (MSVC)
#   set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
# endif()

# # Add googletest directly to our build.
# add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
#                  ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
#                  EXCLUDE_FROM_ALL)
# # Add tests
# FILE( GLOB all_test_srcs src/tests/*.cpp)
# source_group(source FILES ${all_test_srcs})

# add_executable( MeshKernelTests ${all_test_srcs})
# target_link_libraries(MeshKernelTests LINK_PUBLIC ${Boost_LIBRARIES} gtest_main triangle )
# add_test( NAME MeshKernelTests COMMAND MeshKernelTests )

# file(GLOB CurrentNetcdfLibFiles "${CMAKE_CURRENT_SOURCE_DIR}/thirdParty/third_party_open/netcdf/netCDF 4.6.1/bin/*.dll")
# foreach(CurrentNetcdfLibFile IN LISTS CurrentNetcdfLibFiles)
# add_custom_command( TARGET MeshKernelTests 
#                     POST_BUILD
#                     COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CurrentNetcdfLibFile} "${CMAKE_BINARY_DIR}/${build_config}"
#                     COMMENT "Copying netcdf lib file: ${CurrentNetcdfLibFile}")
# endforeach()
