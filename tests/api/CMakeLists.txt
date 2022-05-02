# Tests need to be added as executables first
file(GLOB API_TEST_LIST CONFIGURE_DEPENDS "*.cpp")

# Make a test executable
add_executable(MeshKernelApiTests ${API_TEST_LIST})

# Add a dependency to ensure the correct build order
add_dependencies(MeshKernelApiTests MeshKernelApi)

# Should be linked to the main library, as well as the google test library
target_link_libraries(MeshKernelApiTests PRIVATE MeshKernelApi UtilsStatic
                                                 gtest_main gmock)

if (UNIX)
   target_link_libraries(MeshKernelApiTests PUBLIC ${CMAKE_DL_LIBS})
endif (UNIX)

# Add precompiled header to speed up compilation
target_precompile_headers(MeshKernelApiTests PRIVATE <exception> <memory> <numeric>)

# Make sure that coverage information is produced when using gcc
if(PRODUCE_CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_target_properties( MeshKernelApiTests PROPERTIES COMPILE_FLAGS "-fprofile-arcs -ftest-coverage -static-libgcc -static-libstdc++")
  set_target_properties(MeshKernelApiTests PROPERTIES LINK_FLAGS "-ldl -lgcov --coverage")
endif()

# If you register a test, then ctest and make test will run it. You can also run
# examples and check the output, as well. Command can be a target.
add_test(NAME MeshKernelApiTests COMMAND MeshKernelApiTests)

# On Windows we bundle the netcdf dlls and then copy them to the target file
# directory
if(WIN32)
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
endif()

# Copy the MeshKernel shared library to the target directory
add_custom_command(
  TARGET MeshKernelApiTests
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:MeshKernelApi>"
          "$<TARGET_FILE_DIR:MeshKernelApiTests>"
  COMMENT "Copying MeshKernel shared library")
