# target name
string(CONCAT TARGET_NAME ${PROJECT_NAME} "UnitTests")

# executable
add_executable(${TARGET_NAME})

# source directory
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# list of target sources
set(
  SRC_LIST
  ${SRC_DIR}/AveragingStrategyTests.cpp
  ${SRC_DIR}/AveragingTests.cpp
  ${SRC_DIR}/BoundingBoxTests.cpp
  ${SRC_DIR}/CompoundUndoTests.cpp
  ${SRC_DIR}/ContactsTests.cpp
  ${SRC_DIR}/CurvilinearGridBasicTests.cpp
  ${SRC_DIR}/CurvilinearGridBoundaryPolygonTests.cpp
  ${SRC_DIR}/CurvilinearGridCurvatureTests.cpp
  ${SRC_DIR}/CurvilinearGridDeRefinementTests.cpp
  ${SRC_DIR}/CurvilinearGridFromPolygonTests.cpp
  ${SRC_DIR}/CurvilinearGridFromSplinesTests.cpp
  ${SRC_DIR}/CurvilinearGridFromSplinesTransfiniteTests.cpp
  ${SRC_DIR}/CurvilinearGridGenerationTests.cpp
  ${SRC_DIR}/CurvilinearGridLineAttractionRepulsionTest.cpp
  ${SRC_DIR}/CurvilinearGridLineMirrorTests.cpp
  ${SRC_DIR}/CurvilinearGridLineShiftTests.cpp
  ${SRC_DIR}/CurvilinearGridOrthogonalizationTests.cpp
  ${SRC_DIR}/CurvilinearGridRectangularTests.cpp
  ${SRC_DIR}/CurvilinearGridRefinementTests.cpp
  ${SRC_DIR}/CurvilinearGridSmoothingTests.cpp
  ${SRC_DIR}/CurvilinearGridSmoothnessTests.cpp
  ${SRC_DIR}/CurvilinearGridSnappingTests.cpp
  ${SRC_DIR}/FlipEdgesTests.cpp
  ${SRC_DIR}/FormattingTests.cpp
  ${SRC_DIR}/FunctionsTests.cpp
  ${SRC_DIR}/LandBoundaryTests.cpp
  ${SRC_DIR}/LinearAlgebraTests.cpp
  ${SRC_DIR}/Mesh1DTest.cpp
  ${SRC_DIR}/Mesh2DConnectDDTest.cpp
  ${SRC_DIR}/Mesh2DGlobalGridTests.cpp
  ${SRC_DIR}/Mesh2DTest.cpp
  ${SRC_DIR}/MeshConversionTests.cpp
  ${SRC_DIR}/MeshConversionTests.cpp
  ${SRC_DIR}/MeshRefinementTests.cpp
  ${SRC_DIR}/MeshTests.cpp
  ${SRC_DIR}/MeshTransformationTest.cpp
  ${SRC_DIR}/OrthogonalizationTests.cpp
  ${SRC_DIR}/ParametersTests.cpp
  ${SRC_DIR}/PointTests.cpp
  ${SRC_DIR}/PolygonTests.cpp
  ${SRC_DIR}/PolygonalEnclosureTests.cpp
  ${SRC_DIR}/PolygonsTests.cpp
  ${SRC_DIR}/RangeCheckTests.cpp
  ${SRC_DIR}/SpatialTreesTests.cpp
  ${SRC_DIR}/SplineAlgorithmTests.cpp
  ${SRC_DIR}/SplineTests.cpp
  ${SRC_DIR}/TriangleInterpolationTests.cpp
  ${SRC_DIR}/UndoActionConstructionTests.cpp
  ${SRC_DIR}/UndoStackTests.cpp
  ${SRC_DIR}/UndoTests.cpp
)

 # add sources to target
target_sources(${TARGET_NAME} PRIVATE ${SRC_LIST})

# Should be linked to the main library, as well as the google test library
target_link_libraries(
  ${TARGET_NAME}
  PRIVATE
    MeshKernel
    MeshKernelTestUtils
    gmock
    gtest_main
)

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

# If you register a test, then ctest and make test will run it. You can also run
# examples and check the output, as well. Command can be a target.
add_test(
  NAME ${TARGET_NAME}
  COMMAND ${TARGET_NAME}
)
