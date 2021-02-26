# Set library files
file(GLOB HEADER_LIST CONFIGURE_DEPENDS
     "${PROJECT_SOURCE_DIR}/include/MeshKernel/*.hpp")
file(GLOB SOURCE_LIST CONFIGURE_DEPENDS
     "${PROJECT_SOURCE_DIR}/src/MeshKernel/*.cpp")

# Set api files
set(API_HEADER "${PROJECT_SOURCE_DIR}/include/MeshKernelApi/MeshKernel.hpp")
set(API_SOURCE "${PROJECT_SOURCE_DIR}/src/MeshKernelApi/MeshKernel.cpp")

# Create the static lib
add_library(MeshKernelStatic STATIC ${SOURCE_LIST} ${HEADER_LIST} ${API_HEADER}
                                    ${API_SOURCE})

# Make sure that coverage information is produced when using gcc
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_target_properties(
    MeshKernelStatic
    PROPERTIES COMPILE_FLAGS
               "-fprofile-arcs -ftest-coverage -static-libgcc -static-libstdc++"
  )
  set_target_properties(MeshKernelStatic PROPERTIES LINK_FLAGS
                                                    "-ldl -lgcov --coverage")
endif()

# Expose the interface of the static lib
target_include_directories(MeshKernelStatic
                           PUBLIC "${PROJECT_SOURCE_DIR}/include")

# Add target link dependencies
target_link_libraries(MeshKernelStatic LINK_PUBLIC ${Boost_LIBRARIES} triangle
                      OpenMP::OpenMP_CXX)

# IDEs should put the headers in a nice place
source_group(
  TREE "${PROJECT_SOURCE_DIR}/include/MeshKernel"
  PREFIX "Header Files"
  FILES ${HEADER_LIST})
