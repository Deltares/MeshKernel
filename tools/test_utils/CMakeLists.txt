# project name
project(
  MeshKernelTestUtils
  VERSION 1.0.0
  DESCRIPTION "MeshKernel unit testing utilities"
  LANGUAGES CXX C
)

# target name
set(TARGET_NAME ${PROJECT_NAME})

# static library
add_library(${TARGET_NAME} STATIC)

# source directories
set(SRC_DIR ${PROJECT_SOURCE_DIR}/src)

# include directories
set(INC_DIR ${PROJECT_SOURCE_DIR}/include)
set(DOMAIN_INC_DIR ${INC_DIR}/TestUtils)

set(TEST_FOLDER ${CMAKE_SOURCE_DIR}/data/test)
configure_file(
  ${DOMAIN_INC_DIR}/Definitions.hpp.in
  ${DOMAIN_INC_DIR}/Definitions.hpp
)

# list of target sources
set(
  SRC_LIST
  ${SRC_DIR}/MakeCurvilinearGrids.cpp
  ${SRC_DIR}/MakeMeshes.cpp
  ${SRC_DIR}/SampleFileReader.cpp
)

# list of target headers
set(
  INC_LIST
  ${DOMAIN_INC_DIR}/Definitions.hpp
  ${DOMAIN_INC_DIR}/MakeCurvilinearGrids.hpp
  ${DOMAIN_INC_DIR}/MakeMeshes.hpp
  ${DOMAIN_INC_DIR}/SampleFileReader.hpp
)

# add sources to target
target_sources(
  ${TARGET_NAME}
  PRIVATE
    ${SRC_LIST}
  PUBLIC
    FILE_SET HEADERS
      BASE_DIRS
        ${DOMAIN_INC_DIR}
      FILES
        ${INC_LIST}
)

# NetCDF
find_package(netCDF REQUIRED COMPONENTS C)
message(STATUS "Found NetCDF: in ${netCDF_INSTALL_PREFIX} (found version \"${NetCDFVersion}\")")

# Should be linked to the main (static) library
target_link_libraries(
  ${TARGET_NAME}
  PUBLIC
    MeshKernel
  PRIVATE
    netCDF::netcdf
)

# Expose the interface of the utils library
target_include_directories(${TARGET_NAME} PUBLIC ${INC_DIR})

# group the sources in IDE tree
source_group("Source Files" FILES ${SRC_LIST})
# group the headers in IDE tree
source_group("Header Files" FILES ${INC_LIST})

