# Tests need to be added as executables first
file(GLOB UNIT_TEST_LIST CONFIGURE_DEPENDS "*.cpp")
add_executable(MeshKernelUnitTests ${UNIT_TEST_LIST})

# Should be linked to the main library, as well as the google test library
target_link_libraries(MeshKernelUnitTests PRIVATE MeshKernelStatic UtilsStatic
                                                  gtest_main)

# Make sure that coverage information is produced when using gcc
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_target_properties(
    MeshKernelUnitTests
    PROPERTIES COMPILE_FLAGS
               "-fprofile-arcs -ftest-coverage -static-libgcc -static-libstdc++"
  )
  set_target_properties(MeshKernelUnitTests PROPERTIES LINK_FLAGS
                                                       "-ldl -lgcov --coverage")
endif()

# If you register a test, then ctest and make test will run it. You can also run
# examples and check the output, as well. Command can be a target.
add_test(NAME MeshKernelUnitTests COMMAND MeshKernelUnitTests)

# On Windows we bundle the dlls and then copy them to the target file directory
if(WIN32)
  file(
    GLOB
    NETCDF_FILES_LIST
    CONFIGURE_DEPENDS
    "${NETCDF_BIN_PATH}/${CMAKE_SHARED_LIBRARY_PREFIX}*${CMAKE_SHARED_LIBRARY_SUFFIX}"
  )
  foreach(CurrentNetcdfLibFile IN LISTS NETCDF_FILES_LIST)
    add_custom_command(
      TARGET MeshKernelUnitTests
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CurrentNetcdfLibFile}
              "$<TARGET_FILE_DIR:MeshKernelUnitTests>"
      COMMENT "Copying netcdf lib file: ${CurrentNetcdfLibFile}")
  endforeach()
endif()
