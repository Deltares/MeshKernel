if(DEFINED ENV{BUILD_COUNTER})
  set(BUILD_NUMBER "-dev.$ENV{BUILD_COUNTER}+$ENV{BUILD_NUMBER}")
else()
  set(BUILD_NUMBER 0)
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/version.rc.in
               ${CMAKE_BINARY_DIR}/version.rc @ONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Deltares.MeshKernel.nuspec.in
               ${CMAKE_BINARY_DIR}/Deltares.MeshKernel.nuspec @ONLY)
