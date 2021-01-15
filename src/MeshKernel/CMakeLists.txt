# Set library files
file(GLOB HEADER_LIST CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/include/MeshKernel/*.hpp")
file(GLOB SOURCE_LIST CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/src/MeshKernel/*.cpp")

# Set api files
set(API_HEADER "${PROJECT_SOURCE_DIR}/include/MeshKernelApi/MeshKernel.hpp")
set(API_SOURCE "${PROJECT_SOURCE_DIR}/src/MeshKernelApi/MeshKernel.cpp")

# Create the static lib
add_library(MeshKernelStatic STATIC ${SOURCE_LIST} ${HEADER_LIST} ${API_HEADER} ${API_SOURCE})

# Expose the interface of the static lib
target_include_directories(MeshKernelStatic PUBLIC "${PROJECT_SOURCE_DIR}/include")

# Add target link dependency on boost and triangle
target_link_libraries(MeshKernelStatic LINK_PUBLIC ${Boost_LIBRARIES} triangle)

# All users of this library will need at least C++11
target_compile_features(MeshKernelStatic PUBLIC cxx_std_17)

# IDEs should put the headers in a nice place
source_group(
  TREE "${PROJECT_SOURCE_DIR}/include/MeshKernel"
  PREFIX "Header Files"
  FILES ${HEADER_LIST})
