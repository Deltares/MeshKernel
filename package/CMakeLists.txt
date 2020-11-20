if(DEFINED ENV{VERSION_SUFFIX})
  set(VERSION_SUFFIX $ENV{VERSION_SUFFIX})
else()
  set(VERSION_SUFFIX .0)
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/version.rc.in
               ${CMAKE_BINARY_DIR}/version.rc @ONLY)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Deltares.MeshKernel.Devel.nuspec.in
  ${CMAKE_BINARY_DIR}/package/Deltares.MeshKernel.Devel.nuspec @ONLY)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Deltares.MeshKernel.Release.nuspec.in
  ${CMAKE_BINARY_DIR}/package/Deltares.MeshKernel.Release.nuspec @ONLY)

configure_file(Deltares.MeshKernel.targets
               ${CMAKE_BINARY_DIR}/package/Deltares.MeshKernel.targets COPYONLY)
