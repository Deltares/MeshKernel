# Set api files
file(GLOB HEADER_LIST CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/include/MeshKernelApi/*.hpp")
file(GLOB SOURCE_LIST CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/src/MeshKernelApi/*.cpp")

# Make a shared library
add_library(MeshKernel SHARED ${SOURCE_LIST} ${HEADER_LIST}
                              ${CMAKE_BINARY_DIR}/version.rc)

# Add a dependency to ensure correct build order
add_dependencies(MeshKernel MeshKernelStatic)

# Link shared lib to static lib, but don't expose it
target_link_libraries(MeshKernel PRIVATE MeshKernelStatic)

# IDEs should put the headers in a nice place
source_group(
  TREE "${PROJECT_SOURCE_DIR}/include/MeshKernelApi"
  PREFIX "Header Files"
  FILES ${HEADER_LIST})
