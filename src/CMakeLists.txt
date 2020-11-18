# Note that headers are optional, and do not affect add_library, but they will
# not show up in IDEs unless they are listed in add_library.
file(GLOB HEADER_LIST CONFIGURE_DEPENDS
     "${PROJECT_SOURCE_DIR}/include/MeshKernel/*.hpp")
file(GLOB SOURCE_LIST CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/src/*.cpp")

# Set API variable
set(API_HEADER "${PROJECT_SOURCE_DIR}/include/MeshKernel/MeshKernel.hpp")
set(API_SOURCE "${PROJECT_SOURCE_DIR}/src/MeshKernel.cpp")

# Remove API files from static lib
list(REMOVE_ITEM HEADER_LIST API_HEADER)
list(REMOVE_ITEM SOURCE_LIST API_SOURCE)

# TODO: Get rid of this "bonus" list
file(GLOB BONUS_LIST CONFIGURE_DEPENDS
     "${PROJECT_SOURCE_DIR}/include/MeshKernel/*.cpp")

# Create the static lib
add_library(MeshKernelStatic STATIC ${SOURCE_LIST} ${HEADER_LIST} ${BONUS_LIST})

# Expose the interface of the static lib
target_include_directories(MeshKernelStatic PUBLIC ../include)

# Add target link dependency on boost and triangle
target_link_libraries(MeshKernelStatic LINK_PUBLIC ${Boost_LIBRARIES} triangle)

# All users of this library will need at least C++11
target_compile_features(MeshKernelStatic PUBLIC cxx_std_11)

# Make an automatic library
add_library(MeshKernel SHARED ${API_HEADER} ${API_SOURCE}
                              ${CMAKE_BINARY_DIR}/version.rc)

# Link shared lib to static lib, but don't expose it
target_link_libraries(MeshKernel PRIVATE MeshKernelStatic)

# IDEs should put the headers in a nice place
source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  PREFIX "Header Files"
  FILES ${HEADER_LIST})
