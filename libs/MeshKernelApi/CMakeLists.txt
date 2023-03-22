project(
  MeshKernelApi
  VERSION 1.0.0
  DESCRIPTION "MeshKernel API"
  LANGUAGES CXX C
)

# target name
set(TARGET_NAME ${PROJECT_NAME})

# shared library
add_library(${TARGET_NAME} SHARED ${CMAKE_BINARY_DIR}/version.rc)

# source directory
set(SRC_DIR ${PROJECT_SOURCE_DIR}/src)

# include directory
set(INC_DIR ${PROJECT_SOURCE_DIR}/include)
set(DOMAIN_INC_DIR ${INC_DIR}/MeshKernelApi)

# list of target sources
set(SRC_LIST ${SRC_DIR}/MeshKernel.cpp)

# list of target headers
set(
  INC_LIST
  ${DOMAIN_INC_DIR}/Contacts.hpp
  ${DOMAIN_INC_DIR}/CurvilinearGrid.hpp
  ${DOMAIN_INC_DIR}/CurvilinearParameters.hpp
  ${DOMAIN_INC_DIR}/GeometryList.hpp
  ${DOMAIN_INC_DIR}/MakeGridParameters.hpp
  ${DOMAIN_INC_DIR}/Mesh1D.hpp
  ${DOMAIN_INC_DIR}/Mesh2D.hpp
  ${DOMAIN_INC_DIR}/MeshKernel.hpp
  ${DOMAIN_INC_DIR}/MeshRefinementParameters.hpp
  ${DOMAIN_INC_DIR}/OrthogonalizationParameters.hpp
  ${DOMAIN_INC_DIR}/SplinesToCurvilinearParameters.hpp
  ${DOMAIN_INC_DIR}/State.hpp
  ${DOMAIN_INC_DIR}/Utils.hpp
  ${CMAKE_SOURCE_DIR}/tools/Version/Version.hpp
)

# add sources to target
target_sources(
  ${TARGET_NAME}
  PRIVATE
    ${SRC_LIST}
  PUBLIC
    FILE_SET HEADERS
      BASE_DIRS
        ${INC_DIR} ${CMAKE_SOURCE_DIR}/tools #/Version
      FILES
        ${INC_LIST}
)

# set(VERSION_RESOURCE_FILE ${CMAKE_BINARY_DIR}/version.rc)

# set_target_properties(
#   ${TARGET_NAME} 
#   PROPERTIES 
#   RESOURCE "${VERSION_RESOURCE_FILE}"
# )

# define the public header
set_target_properties(${TARGET_NAME} PROPERTIES PUBLIC_HEADER "${INC_LIST}")

# Expose the interface of the shared lib
target_include_directories(${TARGET_NAME} PUBLIC ${INC_DIR})

# link to static MeshKernel library
target_link_libraries(${TARGET_NAME} PRIVATE MeshKernel)

# Make sure that coverage information is produced when using gcc
if(ENABLE_CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  target_compile_options(
    ${TARGET_NAME}
    PRIVATE
      --coverage # equivalent to -fprofile-arcs -ftest-coverage
      -static-libgcc
      -static-libstdc++
  )
  target_link_libraries(${TARGET_NAME} PRIVATE gcov)
endif()

# Add precompiled header to speed up compilation
target_precompile_headers(
  ${TARGET_NAME}
  PRIVATE
    <cstring>
    <map>
    <memory> 
    <stdexcept>
    <unordered_set>
    <vector> 
)

install(TARGETS ${TARGET_NAME} PUBLIC_HEADER DESTINATION "include/MeshKernelApi")

# group the sources in IDE tree
source_group("Source Files" FILES ${SRC_LIST})

# group the headers in IDE tree
source_group("Header Files" FILES ${INC_LIST})

# Add unit tests
if(ENABLE_UNIT_TESTING)
  add_subdirectory(tests)
endif()