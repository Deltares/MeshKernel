# C++20 translation of the SEPRAN mesh-generation library.
# Replaces the original Fortran sources.

set(target_name Sepran)


# source directories
set(SRC_DIR src)
set(INC_DIR include)

set(
  TARGET_SRC_LIST
  ${SRC_DIR}/SepranSort.cpp
  ${SRC_DIR}/SepranGeometry.cpp
  ${SRC_DIR}/SepranCurveIntersection.cpp
  ${SRC_DIR}/SepranBoundary.cpp
  ${SRC_DIR}/SepranTransform.cpp
  ${SRC_DIR}/SepranTopology.cpp
  ${SRC_DIR}/SepranFront.cpp
  ${SRC_DIR}/SepranQuadratic.cpp
  ${SRC_DIR}/Msho2d.cpp
  ${SRC_DIR}/Mshoce.cpp
)

set(
  TARGET_HEADER_LIST
  ${INC_DIR}/Msho2d.hpp
  ${INC_DIR}/Mshoce.hpp
  ${INC_DIR}/SepranBoundary.hpp
  ${INC_DIR}/SepranContext.hpp
  ${INC_DIR}/SepranCurveIntersection.hpp
  ${INC_DIR}/SepranFront.hpp
  ${INC_DIR}/SepranGeometry.hpp
  ${INC_DIR}/SepranQuadratic.hpp
  ${INC_DIR}/SepranSort.hpp
  ${INC_DIR}/SepranTopology.hpp
  ${INC_DIR}/SepranTransform.hpp
)

add_library(${target_name} STATIC ${TARGET_SRC_LIST} ${TARGET_HEADER_LIST})

# IMPROVEMENT: Explicitly include the include directory for cleaner downstream usage
target_include_directories(${target_name}
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${INC_DIR}>
    $<INSTALL_INTERFACE:include>
)

target_include_directories(${target_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_compile_features(${target_name} PUBLIC cxx_std_20)

# Add precompiled headers to speed up compilation
target_precompile_headers(${target_name}
  PRIVATE
    <algorithm>
    <array>
    <cassert>
    <cmath>
    <iostream>
    <limits>
    <map>
    <memory>
    <numeric>
    <span>
    <stdexcept>
    <string>
    <string_view>
    <utility>
    <vector>
)

set_target_properties(
  ${target_name}
  PROPERTIES
    POSITION_INDEPENDENT_CODE ON
)

source_group("Source Files" FILES ${TARGET_HEADER_LIST}  ${TARGET_SRC_LIST})
