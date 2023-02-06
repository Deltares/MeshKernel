# Tests need to be added as executables first
file(GLOB UNIT_TEST_LIST CONFIGURE_DEPENDS "*.cpp")

set(target_name MeshKernelUnitTests )
add_executable(${target_name} ${UNIT_TEST_LIST})

# Should be linked to the main library, as well as the google test library
target_link_libraries(${target_name}  PRIVATE MeshKernel UtilsStatic
                                                  gtest_main gmock)

if (UNIX)
   target_link_libraries(${target_name} PUBLIC ${CMAKE_DL_LIBS})
endif (UNIX)

# Make sure that coverage information is produced when using gcc
if(PRODUCE_CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_target_properties( ${target_name} PROPERTIES COMPILE_FLAGS "-fprofile-arcs -ftest-coverage -static-libgcc -static-libstdc++")
  set_target_properties(${target_name} PROPERTIES LINK_FLAGS "-ldl -lgcov --coverage")
endif()

# If you register a test, then ctest and make test will run it. You can also run
# examples and check the output, as well. Command can be a target.
add_test(
  NAME ${target_name}
  COMMAND ${target_name}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_BUILD_TYPE}
)

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
      TARGET ${target_name} 
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CurrentNetcdfLibFile}
              "$<TARGET_FILE_DIR:${target_name}>"
      COMMENT "Copying netcdf lib file: ${CurrentNetcdfLibFile}")
  endforeach()
endif()