# Download google test
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG release-1.10.0)
FetchContent_MakeAvailable(googletest)

# Tests need to be added as executables first
file(GLOB TEST_LIST CONFIGURE_DEPENDS "*.cpp")
add_executable(meshkernel_tests ${TEST_LIST})

# I'm using C++11 in the test
target_compile_features(meshkernel_tests PRIVATE cxx_std_11)

# Should be linked to the main library, as well as the Catch2 testing library
target_link_libraries(meshkernel_tests PRIVATE meshkernel ${Boost_LIBRARIES} gtest_main triangle)

# If you register a test, then ctest and make test will run it.
# You can also run examples and check the output, as well.
add_test(NAME meshkernel_tests COMMAND meshkernel_tests) # Command can be a target

file(GLOB CurrentNetcdfLibFiles "${PROJECT_SOURCE_DIR}/thirdParty/third_party_open/netcdf/netCDF 4.6.1/bin/*.dll")
foreach(CurrentNetcdfLibFile IN LISTS CurrentNetcdfLibFiles)
add_custom_command( TARGET meshkernel_tests 
                    POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CurrentNetcdfLibFile} "${CMAKE_BINARY_DIR}/${build_config}"
                    COMMENT "Copying netcdf lib file: ${CurrentNetcdfLibFile}")
endforeach()
