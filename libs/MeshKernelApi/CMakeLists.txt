project(
  MeshKernelApi
  VERSION ${CMAKE_PROJECT_VERSION}
  DESCRIPTION "MeshKernel API"
  LANGUAGES CXX C
)

# target name
set(TARGET_NAME ${PROJECT_NAME})

# shared library
add_library(${TARGET_NAME} SHARED ${CMAKE_BINARY_DIR}/version.rc)

# source directory
set(SRC_DIR ${PROJECT_SOURCE_DIR}/src)
set(CACHE_SRC_DIR ${SRC_DIR}/ApiCache)

# include directory
set(INC_DIR ${PROJECT_SOURCE_DIR}/include)
set(DOMAIN_INC_DIR ${INC_DIR}/MeshKernelApi)
set(CACHE_INC_DIR ${DOMAIN_INC_DIR}/ApiCache)
set(VERSION_INC_DIR ${CMAKE_SOURCE_DIR}/tools)

# list of target sources
set(SRC_LIST
  ${SRC_DIR}/MKStateUndoAction.cpp
  ${SRC_DIR}/MeshKernel.cpp
  ${SRC_DIR}/PropertyCalculator.cpp
  ${SRC_DIR}/EdgeLengthPropertyCalculator.cpp
  ${SRC_DIR}/FaceCircumcenterPropertyCalculator.cpp
  ${SRC_DIR}/OrthogonalityPropertyCalculator.cpp
  ${SRC_DIR}/InterpolatedSamplePropertyCalculator.cpp
)

set(CACHE_SRC_LIST
  ${CACHE_SRC_DIR}/CachedIntegerValues.cpp
  ${CACHE_SRC_DIR}/CachedPointValues.cpp
  ${CACHE_SRC_DIR}/CurvilinearBoundariesAsPolygonCache.cpp
  ${CACHE_SRC_DIR}/FacePolygonPropertyCache.cpp
  ${CACHE_SRC_DIR}/HangingEdgeCache.cpp
  ${CACHE_SRC_DIR}/NodeInPolygonCache.cpp
  ${CACHE_SRC_DIR}/MeshBoundariesAsPolygonCache.cpp
  ${CACHE_SRC_DIR}/ObtuseTriangleCentreCache.cpp
  ${CACHE_SRC_DIR}/PolygonRefinementCache.cpp
  ${CACHE_SRC_DIR}/SmallFlowEdgeCentreCache.cpp
)


# list of target headers
set(
  INC_LIST
  ${DOMAIN_INC_DIR}/BoundingBox.hpp
  ${DOMAIN_INC_DIR}/Contacts.hpp
  ${DOMAIN_INC_DIR}/CurvilinearFrozenLinesAddUndoAction.hpp
  ${DOMAIN_INC_DIR}/CurvilinearFrozenLinesDeleteUndoAction.hpp
  ${DOMAIN_INC_DIR}/CurvilinearGrid.hpp
  ${DOMAIN_INC_DIR}/GeometryList.hpp
  ${DOMAIN_INC_DIR}/GriddedSamples.hpp
  ${DOMAIN_INC_DIR}/MKStateUndoAction.hpp
  ${DOMAIN_INC_DIR}/Mesh1D.hpp
  ${DOMAIN_INC_DIR}/Mesh2D.hpp
  ${DOMAIN_INC_DIR}/MeshKernel.hpp
  ${DOMAIN_INC_DIR}/PropertyCalculator.hpp
  ${DOMAIN_INC_DIR}/EdgeLengthPropertyCalculator.hpp
  ${DOMAIN_INC_DIR}/FaceCircumcenterPropertyCalculator.hpp
  ${DOMAIN_INC_DIR}/OrthogonalityPropertyCalculator.hpp
  ${DOMAIN_INC_DIR}/InterpolatedSamplePropertyCalculator.hpp
  ${DOMAIN_INC_DIR}/State.hpp
  ${DOMAIN_INC_DIR}/Utils.hpp
  ${VERSION_INC_DIR}/Version/Version.hpp
)

set(
  CACHE_INC_LIST
  ${CACHE_INC_DIR}/CachedIntegerValues.hpp
  ${CACHE_INC_DIR}/CachedPointValues.hpp
  ${CACHE_INC_DIR}/CurvilinearBoundariesAsPolygonCache.hpp
  ${CACHE_INC_DIR}/FacePolygonPropertyCache.hpp
  ${CACHE_INC_DIR}/HangingEdgeCache.hpp
  ${CACHE_INC_DIR}/NodeInPolygonCache.hpp
  ${CACHE_INC_DIR}/MeshBoundariesAsPolygonCache.hpp
  ${CACHE_INC_DIR}/ObtuseTriangleCentreCache.hpp
  ${CACHE_INC_DIR}/PolygonRefinementCache.hpp
  ${CACHE_INC_DIR}/SmallFlowEdgeCentreCache.hpp
)


# add sources to target
target_sources(
  ${TARGET_NAME}
  PRIVATE
  ${SRC_LIST}
  ${CACHE_SRC_LIST}
  PUBLIC
    FILE_SET HEADERS
      BASE_DIRS
        ${INC_DIR} ${VERSION_INC_DIR}
      FILES
      ${INC_LIST}
      ${CACHE_INC_LIST}
)

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
    <memory>
    <stdexcept>
    <unordered_map>
    <unordered_set>
    <vector>
)

install(TARGETS ${TARGET_NAME} FILE_SET HEADERS DESTINATION "include")

# group the sources in IDE tree
source_group("Source Files" FILES ${SRC_LIST} ${CACHE_SRC_LIST})

# group the headers in IDE tree
source_group("Header Files" FILES ${INC_LIST} ${CACHE_INC_LIST})

# Add unit tests
if(ENABLE_UNIT_TESTING)
  add_subdirectory(tests)
endif()
