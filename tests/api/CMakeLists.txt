# Tests need to be added as executables first
file(GLOB UNIT_TEST_LIST CONFIGURE_DEPENDS "*.cpp")
add_executable(MeshKernelApiUnitTests ${UNIT_TEST_LIST})

# Should be linked to the main library, as well as the google test library
target_link_libraries(MeshKernelApiUnitTests PRIVATE ${Boost_LIBRARIES} MeshKernelStatic UtilsStatic gtest_main)

# If you register a test, then ctest and make test will run it. You can also run
# examples and check the output, as well. Command can be a target.
add_test(NAME MeshKernelApiTests COMMAND MeshKernelApiUnitTests)
