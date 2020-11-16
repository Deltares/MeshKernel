# Note that headers are optional, and do not affect add_library, but they will
# not show up in IDEs unless they are listed in add_library.
file(GLOB HEADER_LIST CONFIGURE_DEPENDS
     "${PROJECT_SOURCE_DIR}/include/MeshKernel/*.hpp")
file(GLOB SOURCE_LIST CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/src/*.cpp")

# TODO: Get rid of this "bonus" list
file(GLOB BONUS_LIST CONFIGURE_DEPENDS
     "${PROJECT_SOURCE_DIR}/include/MeshKernel/*.cpp")

# Make an automatic library
add_library(MeshKernel SHARED ${SOURCE_LIST} ${HEADER_LIST} ${BONUS_LIST})

# We need this directory, and users of our library will need it too
target_include_directories(MeshKernel PUBLIC ../include)

# Add target link dependency on boost and triangle
target_link_libraries(MeshKernel LINK_PUBLIC ${Boost_LIBRARIES} triangle)

# All users of this library will need at least C++11
target_compile_features(MeshKernel PUBLIC cxx_std_11)

# IDEs should put the headers in a nice place
source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  PREFIX "Header Files"
  FILES ${HEADER_LIST})
