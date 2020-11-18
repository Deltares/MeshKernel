# Tests need to be added as executables first
file(GLOB UNIT_TEST_LIST CONFIGURE_DEPENDS "*.cpp")
add_executable(MeshKernelUnitTests ${UNIT_TEST_LIST})

# Should be linked to the main library, as well as the google test library
target_link_libraries(
  MeshKernelUnitTests PRIVATE MeshKernelStatic UtilsStatic ${Boost_LIBRARIES}
                              triangle gtest_main)

# If you register a test, then ctest and make test will run it. You can also run
# examples and check the output, as well. Command can be a target.
add_test(NAME MeshKernelTests COMMAND MeshKernelUnitTests)

file(GLOB NETCDF_FILES_LIST CONFIGURE_DEPENDS
     "${PROJECT_SOURCE_DIR}/extern/netcdf/netCDF 4.6.1/bin/*.dll")
foreach(CurrentNetcdfLibFile IN LISTS NETCDF_FILES_LIST)
  add_custom_command(
    TARGET MeshKernelUnitTests
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CurrentNetcdfLibFile}
            "$<TARGET_FILE_DIR:MeshKernelUnitTests>"
    COMMENT "Copying netcdf lib file: ${CurrentNetcdfLibFile}")
endforeach()
