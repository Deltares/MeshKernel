# project name
project(
  MeshKernel
  VERSION ${CMAKE_PROJECT_VERSION}
  DESCRIPTION "MeshKernel static library"
  LANGUAGES CXX C
)

# target name
set(TARGET_NAME ${PROJECT_NAME})

# static library
add_library(${TARGET_NAME} STATIC)

# source directories
set(SRC_DIR ${PROJECT_SOURCE_DIR}/src)
set(AVERAGING_STRATEGIES_SRC_DIR ${SRC_DIR}/AveragingStrategies)
set(CURVILINEAR_GRID_SRC_DIR ${SRC_DIR}/CurvilinearGrid)
set(UTILITIES_SRC_DIR ${SRC_DIR}/Utilities)

# include directories
set(INC_DIR ${PROJECT_SOURCE_DIR}/include)
set(DOMAIN_INC_DIR ${INC_DIR}/MeshKernel)
set(AVERAGING_STRATEGIES_INC_DIR ${DOMAIN_INC_DIR}/AveragingStrategies)
set(CURVILINEAR_GRID_INC_DIR ${DOMAIN_INC_DIR}/CurvilinearGrid)
set(UTILITIES_INC_DIR ${DOMAIN_INC_DIR}/Utilities)

# list of target sources
set(
  SRC_LIST
  ${SRC_DIR}/AveragingInterpolation.cpp
  ${SRC_DIR}/BilinearInterpolationOnGriddedSamples.cpp
  ${SRC_DIR}/ConnectMeshes.cpp
  ${SRC_DIR}/Contacts.cpp
  ${SRC_DIR}/Definitions.cpp
  ${SRC_DIR}/Entities.cpp
  ${SRC_DIR}/FlipEdges.cpp
  ${SRC_DIR}/LandBoundaries.cpp
  ${SRC_DIR}/LandBoundary.cpp
  ${SRC_DIR}/Mesh.cpp
  ${SRC_DIR}/Mesh1D.cpp
  ${SRC_DIR}/Mesh2D.cpp
  ${SRC_DIR}/MeshRefinement.cpp
  ${SRC_DIR}/Network1D.cpp
  ${SRC_DIR}/Operations.cpp
  ${SRC_DIR}/OrthogonalizationAndSmoothing.cpp
  ${SRC_DIR}/Orthogonalizer.cpp
  ${SRC_DIR}/Polygon.cpp
  ${SRC_DIR}/PolygonalEnclosure.cpp
  ${SRC_DIR}/Polygons.cpp
  ${SRC_DIR}/RemoveDisconnectedRegions.cpp
  ${SRC_DIR}/Smoother.cpp
  ${SRC_DIR}/SplineAlgorithms.cpp
  ${SRC_DIR}/Splines.cpp
  ${SRC_DIR}/TriangulationInterpolation.cpp
  ${SRC_DIR}/TriangulationWrapper.cpp
)

set(
  AVERAGING_STRATEGIES_SRC_LIST
  ${AVERAGING_STRATEGIES_SRC_DIR}/AveragingStrategyFactory.cpp
  ${AVERAGING_STRATEGIES_SRC_DIR}/ClosestAveragingStrategy.cpp
  ${AVERAGING_STRATEGIES_SRC_DIR}/InverseWeightedAveragingStrategy.cpp
  ${AVERAGING_STRATEGIES_SRC_DIR}/MaxAveragingStrategy.cpp
  ${AVERAGING_STRATEGIES_SRC_DIR}/MinAbsAveragingStrategy.cpp
  ${AVERAGING_STRATEGIES_SRC_DIR}/MinAveragingStrategy.cpp
  ${AVERAGING_STRATEGIES_SRC_DIR}/SimpleAveragingStrategy.cpp
)

set(
  CURVILINEAR_GRID_SRC_LIST
  ${CURVILINEAR_GRID_SRC_DIR}/CurvilinearGrid.cpp
  ${CURVILINEAR_GRID_SRC_DIR}/CurvilinearGridAlgorithm.cpp
  ${CURVILINEAR_GRID_SRC_DIR}/CurvilinearGridDeRefinement.cpp
  ${CURVILINEAR_GRID_SRC_DIR}/CurvilinearGridFromPolygon.cpp
  ${CURVILINEAR_GRID_SRC_DIR}/CurvilinearGridFromSplines.cpp
  ${CURVILINEAR_GRID_SRC_DIR}/CurvilinearGridFromSplinesTransfinite.cpp
  ${CURVILINEAR_GRID_SRC_DIR}/CurvilinearGridLine.cpp
  ${CURVILINEAR_GRID_SRC_DIR}/CurvilinearGridLineAttractionRepulsion.cpp
  ${CURVILINEAR_GRID_SRC_DIR}/CurvilinearGridLineMirror.cpp
  ${CURVILINEAR_GRID_SRC_DIR}/CurvilinearGridLineShift.cpp
  ${CURVILINEAR_GRID_SRC_DIR}/CurvilinearGridOrthogonalization.cpp
  ${CURVILINEAR_GRID_SRC_DIR}/CurvilinearGridRectangular.cpp
  ${CURVILINEAR_GRID_SRC_DIR}/CurvilinearGridRefinement.cpp
  ${CURVILINEAR_GRID_SRC_DIR}/CurvilinearGridSmoothing.cpp
  ${CURVILINEAR_GRID_SRC_DIR}/CurvilinearGridSnapping.cpp
)

set(
  UTILITIES_SRC_LIST
  ${UTILITIES_SRC_DIR}/RTree.cpp
)

# list of target headers
set(
  INC_LIST
  ${DOMAIN_INC_DIR}/AveragingInterpolation.hpp
  ${DOMAIN_INC_DIR}/BilinearInterpolationOnGriddedSamples.hpp
  ${DOMAIN_INC_DIR}/BoundingBox.hpp
  ${DOMAIN_INC_DIR}/ConnectMeshes.hpp
  ${DOMAIN_INC_DIR}/Constants.hpp
  ${DOMAIN_INC_DIR}/Contacts.hpp
  ${DOMAIN_INC_DIR}/Definitions.hpp
  ${DOMAIN_INC_DIR}/Entities.hpp
  ${DOMAIN_INC_DIR}/Exceptions.hpp
  ${DOMAIN_INC_DIR}/FlipEdges.hpp
  ${DOMAIN_INC_DIR}/Formatting.hpp
  ${DOMAIN_INC_DIR}/LandBoundaries.hpp
  ${DOMAIN_INC_DIR}/LandBoundary.hpp
  ${DOMAIN_INC_DIR}/Mesh.hpp
  ${DOMAIN_INC_DIR}/Mesh1D.hpp
  ${DOMAIN_INC_DIR}/Mesh2D.hpp
  ${DOMAIN_INC_DIR}/MeshInterpolation.hpp
  ${DOMAIN_INC_DIR}/MeshRefinement.hpp
  ${DOMAIN_INC_DIR}/Network1D.hpp
  ${DOMAIN_INC_DIR}/Operations.hpp
  ${DOMAIN_INC_DIR}/OrthogonalizationAndSmoothing.hpp
  ${DOMAIN_INC_DIR}/Orthogonalizer.hpp
  ${DOMAIN_INC_DIR}/Parameters.hpp
  ${DOMAIN_INC_DIR}/Point.hpp
  ${DOMAIN_INC_DIR}/Polygon.hpp
  ${DOMAIN_INC_DIR}/PolygonalEnclosure.hpp
  ${DOMAIN_INC_DIR}/Polygons.hpp
  ${DOMAIN_INC_DIR}/RangeCheck.hpp
  ${DOMAIN_INC_DIR}/RemoveDisconnectedRegions.hpp
  ${DOMAIN_INC_DIR}/Smoother.hpp
  ${DOMAIN_INC_DIR}/SplineAlgorithms.hpp
  ${DOMAIN_INC_DIR}/Splines.hpp
  ${DOMAIN_INC_DIR}/TriangulationInterpolation.hpp
  ${DOMAIN_INC_DIR}/TriangulationWrapper.hpp
  ${DOMAIN_INC_DIR}/Vector.hpp
)

set(
  AVERAGING_STRATEGIES_INC_LIST
  ${AVERAGING_STRATEGIES_INC_DIR}/AveragingStrategy.hpp
  ${AVERAGING_STRATEGIES_INC_DIR}/AveragingStrategyFactory.hpp
  ${AVERAGING_STRATEGIES_INC_DIR}/ClosestAveragingStrategy.hpp
  ${AVERAGING_STRATEGIES_INC_DIR}/InverseWeightedAveragingStrategy.hpp
  ${AVERAGING_STRATEGIES_INC_DIR}/MaxAveragingStrategy.hpp
  ${AVERAGING_STRATEGIES_INC_DIR}/MinAbsAveragingStrategy.hpp
  ${AVERAGING_STRATEGIES_INC_DIR}/MinAveragingStrategy.hpp
  ${AVERAGING_STRATEGIES_INC_DIR}/SimpleAveragingStrategy.hpp
)

set(
  CURVILINEAR_GRID_INC_LIST
  ${CURVILINEAR_GRID_INC_DIR}/CurvilinearGrid.hpp
  ${CURVILINEAR_GRID_INC_DIR}/CurvilinearGridAlgorithm.hpp
  ${CURVILINEAR_GRID_INC_DIR}/CurvilinearGridDeRefinement.hpp
  ${CURVILINEAR_GRID_INC_DIR}/CurvilinearGridFromPolygon.hpp
  ${CURVILINEAR_GRID_INC_DIR}/CurvilinearGridFromSplines.hpp
  ${CURVILINEAR_GRID_INC_DIR}/CurvilinearGridFromSplinesTransfinite.hpp
  ${CURVILINEAR_GRID_INC_DIR}/CurvilinearGridLine.hpp
  ${CURVILINEAR_GRID_INC_DIR}/CurvilinearGridLineAttractionRepulsion.hpp
  ${CURVILINEAR_GRID_INC_DIR}/CurvilinearGridLineMirror.hpp
  ${CURVILINEAR_GRID_INC_DIR}/CurvilinearGridLineShift.hpp
  ${CURVILINEAR_GRID_INC_DIR}/CurvilinearGridNodeIndices.hpp
  ${CURVILINEAR_GRID_INC_DIR}/CurvilinearGridOrthogonalization.hpp
  ${CURVILINEAR_GRID_INC_DIR}/CurvilinearGridRectangular.hpp
  ${CURVILINEAR_GRID_INC_DIR}/CurvilinearGridRefinement.hpp
  ${CURVILINEAR_GRID_INC_DIR}/CurvilinearGridSmoothing.hpp
  ${CURVILINEAR_GRID_INC_DIR}/CurvilinearGridSnapping.hpp
)

set(
  UTILITIES_INC_LIST
  ${UTILITIES_INC_DIR}/LinearAlgebra.hpp
  ${UTILITIES_INC_DIR}/NumericFunctions.hpp
  ${UTILITIES_INC_DIR}/RTree.hpp
)

set(
  HEADERS
  ${INC_LIST}
  ${AVERAGING_STRATEGIES_INC_LIST}
  ${CURVILINEAR_GRID_INC_LIST}
  ${UTILITIES_INC_LIST}
)

# add sources to target
target_sources(
  ${TARGET_NAME}
  PRIVATE
    ${SRC_LIST}
    ${AVERAGING_STRATEGIES_SRC_LIST}
    ${CURVILINEAR_GRID_SRC_LIST}
    ${UTILITIES_SRC_LIST}
  PUBLIC
    FILE_SET HEADERS
      BASE_DIRS
        ${INC_DIR}
      FILES
        ${HEADERS}
)

# set the public headers
set_target_properties(${TARGET_NAME} PROPERTIES PUBLIC_HEADER "${HEADERS}")

# Expose the interface of the static lib
target_include_directories(  ${TARGET_NAME} PUBLIC ${INC_DIR})

# Add target link dependencies

# Boost
find_package(Boost REQUIRED)
target_include_directories(${TARGET_NAME} SYSTEM PUBLIC ${Boost_INCLUDE_DIR})

# OpenMP
find_package(OpenMP REQUIRED)

# link required libraries
target_link_libraries(
  ${TARGET_NAME}
  LINK_PUBLIC
    Triangle
    OpenMP::OpenMP_CXX
    Eigen3::Eigen
)

if(${USE_LIBFMT})
  target_link_libraries(${TARGET_NAME} PUBLIC fmt::fmt)
endif()

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

# Add precompiled headers to speed up compilation
target_precompile_headers(
  ${TARGET_NAME}
  PRIVATE
    <algorithm>
    <array>
    <cassert>
    <cmath>
    <initializer_list>
    <limits>
    <map>
    <memory>
    <numeric>
    <source_location>
    <stdexcept>
    <sstream>
    <string>
    <string_view>
    <tuple>
    <type_traits>
    <unordered_map>
    <utility>
    <vector>
)

install(TARGETS ${TARGET_NAME} PUBLIC_HEADER DESTINATION "include/MeshKernel")

# group the sources in IDE tree
source_group("Source Files" FILES ${SRC_LIST})
source_group("Source Files/AveragingStrategies" FILES ${AVERAGING_STRATEGIES_SRC_LIST})
source_group("Source Files/CurvilinearGrid" FILES ${CURVILINEAR_GRID_SRC_LIST})
source_group("Source Files/Utilities" FILES ${UTILITIES_SRC_LIST})

# group the headers in IDE tree
source_group("Header Files" FILES ${INC_LIST})
source_group("Header Files/AveragingStrategies" FILES ${AVERAGING_STRATEGIES_INC_LIST})
source_group("Header Files/CurvilinearGrid" FILES ${CURVILINEAR_GRID_INC_LIST})
source_group("Header Files/Utilities" FILES ${UTILITIES_INC_LIST})

# Add unit tests
if(ENABLE_UNIT_TESTING)
  add_subdirectory(tests)
endif()

# Add benchmarks
if(ENABLE_BENCHMARKING)
  add_subdirectory(benchmark)
endif()
