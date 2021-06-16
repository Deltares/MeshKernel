# Configure version file
configure_file("${PROJECT_SOURCE_DIR}/include/MeshKernelApi/Version.hpp.in"
               "${PROJECT_SOURCE_DIR}/include/MeshKernelApi/Version.hpp" @ONLY)

# Set api files
file(GLOB HEADER_LIST CONFIGURE_DEPENDS
     "${PROJECT_SOURCE_DIR}/include/MeshKernelApi/*.hpp")
file(GLOB SOURCE_LIST CONFIGURE_DEPENDS
     "${PROJECT_SOURCE_DIR}/src/MeshKernelApi/*.cpp")

# Make a shared library
add_library(MeshKernelApi SHARED ${SOURCE_LIST} ${HEADER_LIST}
                                 ${CMAKE_BINARY_DIR}/version.rc)

# Expose the interface of the shared lib
target_include_directories(MeshKernelApi PUBLIC "${PROJECT_SOURCE_DIR}/include")

# Make sure that coverage information is produced when using gcc
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_target_properties(
    MeshKernelApi
    PROPERTIES COMPILE_FLAGS
               "-fprofile-arcs -ftest-coverage -static-libgcc -static-libstdc++"
  )
  set_target_properties(MeshKernelApi PROPERTIES LINK_FLAGS
                                                 "-ldl -lgcov --coverage")
endif()

# Link shared lib to static lib, but don't expose it
target_link_libraries(MeshKernelApi PRIVATE MeshKernel)

# IDEs should put the headers in a nice place
source_group(
  TREE "${PROJECT_SOURCE_DIR}/include/MeshKernelApi"
  PREFIX "Header Files"
  FILES ${HEADER_LIST})
