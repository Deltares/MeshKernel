# Download google test
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG release-1.10.0)
FetchContent_MakeAvailable(googletest)

# Tests need to be added as executables first
file(GLOB TEST_LIST CONFIGURE_DEPENDS "*.cpp")
add_executable(MeshKernelTests ${TEST_LIST})

# Using C++11 in the test
target_compile_features(MeshKernelTests PRIVATE cxx_std_11)

# Should be linked to the main library, as well as the Catch2 testing library
target_link_libraries(MeshKernelTests PRIVATE MeshKernel ${Boost_LIBRARIES}
                                              gtest_main triangle)

# If you register a test, then ctest and make test will run it. You can also run
# examples and check the output, as well.
add_test(NAME MeshKernelTests COMMAND MeshKernelTests) # Command can be a target

file(GLOB NETCDF_FILES_LIST CONFIGURE_DEPENDS
     "${PROJECT_SOURCE_DIR}/extern/netcdf/netCDF 4.6.1/bin/*.dll")
foreach(CurrentNetcdfLibFile IN LISTS NETCDF_FILES_LIST)
  add_custom_command(
    TARGET MeshKernelTests
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CurrentNetcdfLibFile}
            "$<TARGET_FILE_DIR:MeshKernelTests>"
    COMMENT "Copying netcdf lib file: ${CurrentNetcdfLibFile}")
endforeach()
