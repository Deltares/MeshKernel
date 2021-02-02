# Tests need to be added as executables first
file(GLOB API_TEST_LIST CONFIGURE_DEPENDS "*.cpp")

# Make a test executable
add_executable(MeshKernelApiTests ${API_TEST_LIST})

# Add a dependency to ensure the correct build order
add_dependencies(MeshKernelApiTests MeshKernelStatic)

# Should be linked to the main library, as well as the google test library
target_link_libraries(
  MeshKernelApiTests PRIVATE ${Boost_LIBRARIES} MeshKernelStatic UtilsStatic
                             gtest_main)

# If you register a test, then ctest and make test will run it. You can also run
# examples and check the output, as well. Command can be a target.
add_test(NAME MeshKernelApiTests COMMAND MeshKernelApiTests)

file(GLOB NETCDF_FILES_LIST CONFIGURE_DEPENDS
     "${PROJECT_SOURCE_DIR}/extern/netcdf/netCDF 4.6.1/bin/*.dll")
foreach(CurrentNetcdfLibFile IN LISTS NETCDF_FILES_LIST)
  add_custom_command(
    TARGET MeshKernelApiTests
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CurrentNetcdfLibFile}
            "$<TARGET_FILE_DIR:MeshKernelApiTests>"
    COMMENT "Copying netcdf lib file: ${CurrentNetcdfLibFile}")
endforeach()
