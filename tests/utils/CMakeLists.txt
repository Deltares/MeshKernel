# Note that headers are optional, and do not affect add_library, but they will
# not show up in IDEs unless they are listed in add_library.
file(GLOB UTILS_HEADER_LIST CONFIGURE_DEPENDS
     "${CMAKE_CURRENT_SOURCE_DIR}/include/TestUtils/*.hpp")
file(GLOB UTILS_SOURCE_LIST CONFIGURE_DEPENDS
     "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

# Create the utils library
add_library(UtilsStatic STATIC ${UTILS_HEADER_LIST} ${UTILS_SOURCE_LIST})

# Should be linked to the main (static) library
target_link_libraries(UtilsStatic PRIVATE MeshKernel)

# Expose the interface of the utils library
target_include_directories(UtilsStatic PUBLIC include)
