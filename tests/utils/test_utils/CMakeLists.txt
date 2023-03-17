# Note that headers are optional, and do not affect add_library, but they will
# not show up in IDEs unless they are listed in add_library.
file(GLOB UTILS_HEADER_LIST CONFIGURE_DEPENDS
     "${CMAKE_CURRENT_SOURCE_DIR}/include/TestUtils/*.hpp")
file(GLOB UTILS_SOURCE_LIST CONFIGURE_DEPENDS
     "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

set(target_name MeshKernelTestUtils)

# Create the utils library
add_library(${target_name} STATIC ${UTILS_HEADER_LIST} ${UTILS_SOURCE_LIST})

# NetCDF
find_package(netCDF REQUIRED COMPONENTS C)
if (netCDF_FOUND)
  message(STATUS "Found NetCDF: in ${netCDF_INSTALL_PREFIX} (found version \"${NetCDFVersion}\")")
endif()

# Should be linked to the main (static) library
target_link_libraries(
  ${target_name}
  PRIVATE
  MeshKernel
  PUBLIC
  netCDF::netcdf
)

# Expose the interface of the utils library
target_include_directories(${target_name} PUBLIC include)
