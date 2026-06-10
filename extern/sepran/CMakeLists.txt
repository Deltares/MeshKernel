# C++20 translation of the SEPRAN mesh-generation library.
# Replaces the original Fortran sources.

set(target_name Sepran)

add_library(${target_name} STATIC)

set(
  TARGET_SRC_LIST
  SepranSort.cpp
  SepranGeometry.cpp
  SepranCurveIntersection.cpp
  SepranBoundary.cpp
  SepranTransform.cpp
  SepranTopology.cpp
  SepranFront.cpp
  SepranQuadratic.cpp
  Msho2d.cpp
  Mshoce.cpp
)

target_sources(${target_name} PRIVATE ${TARGET_SRC_LIST})

target_include_directories(${target_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_compile_features(${target_name} PUBLIC cxx_std_20)

set_target_properties(
  ${target_name}
  PROPERTIES
    POSITION_INDEPENDENT_CODE ON
)

source_group("Source Files" FILES ${TARGET_SRC_LIST})