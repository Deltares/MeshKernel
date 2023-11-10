configure_file(${CMAKE_CURRENT_SOURCE_DIR}/version.rc.in
               ${CMAKE_BINARY_DIR}/version.rc @ONLY)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Deltares.MeshKernel.Development.nuspec.in
  ${CMAKE_BINARY_DIR}/package/Deltares.MeshKernel.Development.nuspec @ONLY)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Deltares.MeshKernel.Release.nuspec.in
  ${CMAKE_BINARY_DIR}/package/Deltares.MeshKernel.Release.nuspec @ONLY)

configure_file(Deltares.MeshKernel.Development.targets
               ${CMAKE_BINARY_DIR}/package/Deltares.MeshKernel.Development.targets COPYONLY)

configure_file(Deltares.MeshKernel.Release.targets
               ${CMAKE_BINARY_DIR}/package/Deltares.MeshKernel.Release.targets COPYONLY)