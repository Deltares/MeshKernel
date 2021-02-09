#include <memory>

#include <gtest/gtest.h>

#include <MeshKernelApi/GeometryList.hpp>
#include <MeshKernelApi/MakeMeshParameters.hpp>
#include <MeshKernelApi/Mesh1d.hpp>
#include <MeshKernelApi/MeshGeometry.hpp>
#include <MeshKernelApi/MeshGeometryDimensions.hpp>
#include <MeshKernelApi/MeshKernel.hpp>
#include <MeshKernelApi/Utils.hpp>
#include <TestUtils/MakeMeshes.hpp>

class ApiTests : public ::testing::Test
{
public:
    void AllocateNewMesh(int& meshKernelId) const
    {
        // Allocate new mesh
        const auto successful = meshkernelapi::mkernel_new_mesh(meshKernelId);
        ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, successful);
    }

    void MakeMesh(int n = 4, int m = 3, double delta = 1.0) const
    {
        // Allocate new mesh
        int meshKernelId;
        AllocateNewMesh(meshKernelId);

        // Set-up new mesh
        auto meshData = MakeRectangularMeshForApiTesting(n, m, delta);
        auto errorCode = mkernel_set_mesh2d(meshKernelId, std::get<1>(meshData), std::get<0>(meshData), false);
        ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
        // Clean up client-allocated memory
        DeleteRectangularMeshForApiTesting(std::get<0>(meshData));
    }

protected:
    void TearDown() override
    {
        //De-allocate meshkernel
        const auto errorCode = meshkernelapi::mkernel_deallocate_state(0);
        ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    }
};

TEST_F(ApiTests, DeleteNodeThroughApi)
{
    // Prepare
    MakeMesh();

    // Execute
    auto errorCode = meshkernelapi::mkernel_delete_node(0, 0);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::MeshGeometryDimensions meshGeometryDimensions{};
    meshkernelapi::MeshGeometry meshGeometry{};
    errorCode = mkernel_get_mesh2d(0, meshGeometryDimensions, meshGeometry);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(11, meshGeometryDimensions.numnode);
    ASSERT_EQ(15, meshGeometryDimensions.numedge);
}

TEST_F(ApiTests, FlipEdgesThroughApi)
{
    // Prepare
    MakeMesh();

    // Execute
    const int isTriangulationRequired = 1;
    const int projectToLandBoundaryOption = 1;
    auto errorCode = meshkernelapi::mkernel_flip_edges(0, isTriangulationRequired, projectToLandBoundaryOption);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::MeshGeometryDimensions meshGeometryDimensions{};
    meshkernelapi::MeshGeometry meshGeometry{};
    errorCode = mkernel_get_mesh2d(0, meshGeometryDimensions, meshGeometry);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(12, meshGeometryDimensions.numnode);
    ASSERT_EQ(23, meshGeometryDimensions.numedge);
}

TEST_F(ApiTests, InsertEdgeThroughApi)
{
    // Prepare
    MakeMesh();

    // Execute
    int newEdgeIndex;
    auto errorCode = meshkernelapi::mkernel_insert_edge(0, 0, 4, newEdgeIndex);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(17, newEdgeIndex);

    meshkernelapi::MeshGeometryDimensions meshGeometryDimensions{};
    meshkernelapi::MeshGeometry meshGeometry{};
    errorCode = mkernel_get_mesh2d(0, meshGeometryDimensions, meshGeometry);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(12, meshGeometryDimensions.numnode);
    ASSERT_EQ(18, meshGeometryDimensions.numedge);
}

TEST_F(ApiTests, MergeTwoNodesThroughApi)
{
    // Prepare
    MakeMesh();

    // Execute
    auto errorCode = meshkernelapi::mkernel_merge_two_nodes(0, 0, 4);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::MeshGeometryDimensions meshGeometryDimensions{};
    meshkernelapi::MeshGeometry meshGeometry{};
    errorCode = mkernel_get_mesh2d(0, meshGeometryDimensions, meshGeometry);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(11, meshGeometryDimensions.numnode);
    ASSERT_EQ(15, meshGeometryDimensions.numedge);
}

TEST_F(ApiTests, MergeNodesThroughApi)
{
    // Prepare
    MakeMesh();
    meshkernelapi::GeometryList geometry_list{};

    // Execute
    auto errorCode = mkernel_merge_nodes(0, geometry_list);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::MeshGeometryDimensions meshGeometryDimensions{};
    meshkernelapi::MeshGeometry meshGeometry{};
    errorCode = mkernel_get_mesh2d(0, meshGeometryDimensions, meshGeometry);

    // Assert (nothing is done, just check that the api communication works)
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(12, meshGeometryDimensions.numnode);
    ASSERT_EQ(17, meshGeometryDimensions.numedge);
}

TEST_F(ApiTests, OrthogonalizationThroughApi)
{
    // Set a new mesh in mesh
    MakeMesh();

    // Prepare
    meshkernelapi::OrthogonalizationParameters orthogonalizationParameters{};
    orthogonalizationParameters.OuterIterations = 1;
    orthogonalizationParameters.BoundaryIterations = 25;
    orthogonalizationParameters.InnerIterations = 25;
    orthogonalizationParameters.OrthogonalizationToSmoothingFactor = 0.975;

    meshkernelapi::GeometryList geometryList{};
    meshkernelapi::GeometryList landBoundaries{};

    // Execute
    auto errorCode = mkernel_orthogonalize_initialize(0,
                                                      1,
                                                      orthogonalizationParameters,
                                                      landBoundaries,
                                                      geometryList);

    // Assert (nothing is done, just check that the api communication works)
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    errorCode = meshkernelapi::mkernel_orthogonalize_prepare_outer_iteration(0);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    errorCode = meshkernelapi::mkernel_orthogonalize_inner_iteration(0);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    errorCode = meshkernelapi::mkernel_orthogonalize_finalize_outer_iteration(0);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    errorCode = meshkernelapi::mkernel_orthogonalize_delete(0);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::MeshGeometryDimensions meshGeometryDimensions{};
    meshkernelapi::MeshGeometry meshGeometry{};
    errorCode = mkernel_get_mesh2d(0, meshGeometryDimensions, meshGeometry);

    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(12, meshGeometryDimensions.numnode);
    ASSERT_EQ(17, meshGeometryDimensions.numedge);
}

TEST_F(ApiTests, MakeGridThroughApi)
{
    // Prepare
    int meshKernelId;
    AllocateNewMesh(meshKernelId);

    meshkernelapi::MakeMeshParameters makeMeshParameters{};
    meshkernelapi::GeometryList geometryList{};

    makeMeshParameters.GridType = 0;
    makeMeshParameters.NumberOfColumns = 4;
    makeMeshParameters.NumberOfRows = 3;
    makeMeshParameters.GridAngle = 0.0;
    makeMeshParameters.GridBlockSize = 0.0;
    makeMeshParameters.OriginXCoordinate = 0.0;
    makeMeshParameters.OriginYCoordinate = 0.0;
    makeMeshParameters.OriginZCoordinate = 0.0;
    makeMeshParameters.XGridBlockSize = 0.0;
    makeMeshParameters.YGridBlockSize = 0.0;

    // Execute
    auto errorCode = mkernel_make_mesh(0, makeMeshParameters, geometryList);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::MeshGeometryDimensions meshGeometryDimensions{};
    meshkernelapi::MeshGeometry meshGeometry{};
    errorCode = mkernel_get_mesh2d(0, meshGeometryDimensions, meshGeometry);

    // Assert (nothing is done, just check that the api communication works)
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(20, meshGeometryDimensions.numnode);
    ASSERT_EQ(31, meshGeometryDimensions.numedge);
}
TEST_F(ApiTests, GenerateTransfiniteCurvilinearGridThroughApi)
{
    // Prepare
    int meshKernelId;
    AllocateNewMesh(meshKernelId);

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometrySeparator = meshkernel::doubleMissingValue;

    geometryListIn.xCoordinates = new double[]{1.340015E+02, 3.642529E+02, 6.927549E+02, meshkernel::doubleMissingValue,
                                               2.585022E+02, 4.550035E+02, 8.337558E+02, meshkernel::doubleMissingValue,
                                               1.002513E+02, 4.610035E+02, meshkernel::doubleMissingValue,
                                               6.522547E+02, 7.197551E+02};

    geometryListIn.yCoordinates = new double[]{
        2.546282E+02, 4.586302E+02, 5.441311E+02, meshkernel::doubleMissingValue,
        6.862631E+01, 2.726284E+02, 3.753794E+02, meshkernel::doubleMissingValue,
        4.068797E+02, 7.912642E+01, meshkernel::doubleMissingValue,
        6.026317E+02, 2.681283E+02};

    geometryListIn.zCoordinates = new double[]{
        0.0, 0.0, 0.0, meshkernel::doubleMissingValue,
        0.0, 0.0, 0.0, meshkernel::doubleMissingValue,
        0.0, 0.0, meshkernel::doubleMissingValue,
        0.0, 0.0};

    geometryListIn.numberOfCoordinates = 13;
    meshkernelapi::CurvilinearParameters curvilinearParameters;
    curvilinearParameters.MRefinement = 10;
    curvilinearParameters.NRefinement = 10;
    curvilinearParameters.SmoothingIterations = 10;
    curvilinearParameters.SmoothingParameter = 0.5;
    curvilinearParameters.AttractionParameter = 0.0;

    // Execute
    auto errorCode = mkernel_curvilinear_mesh_from_splines(0, geometryListIn, curvilinearParameters);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::MeshGeometryDimensions meshGeometryDimensions{};
    meshkernelapi::MeshGeometry meshGeometry{};
    errorCode = mkernel_get_mesh2d(0, meshGeometryDimensions, meshGeometry);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(121, meshGeometryDimensions.numnode);
    ASSERT_EQ(220, meshGeometryDimensions.numedge);

    // Delete dynamically allocated memory with operator new
    delete[] geometryListIn.xCoordinates;
    delete[] geometryListIn.yCoordinates;
    delete[] geometryListIn.zCoordinates;
}

TEST_F(ApiTests, GenerateOrthogonalCurvilinearGridThroughApi)
{
    // Prepare
    int meshKernelId;
    AllocateNewMesh(meshKernelId);

    meshkernelapi::GeometryList geometryListIn;

    geometryListIn.geometrySeparator = meshkernel::doubleMissingValue;

    geometryListIn.xCoordinates = new double[]{1.175014E+02, 3.755030E+02, 7.730054E+02, meshkernel::doubleMissingValue,
                                               4.100089E+01, 3.410027E+02};

    geometryListIn.yCoordinates = new double[]{2.437587E+01, 3.266289E+02, 4.563802E+02, meshkernel::doubleMissingValue,
                                               2.388780E+02, 2.137584E+01};

    geometryListIn.zCoordinates = new double[]{0.0, 0.0, 0.0, meshkernel::doubleMissingValue,
                                               0.0, 0.0, meshkernel::doubleMissingValue};
    geometryListIn.numberOfCoordinates = 6;

    meshkernelapi::CurvilinearParameters curvilinearParameters;
    curvilinearParameters.MRefinement = 40;
    curvilinearParameters.NRefinement = 10;
    meshkernelapi::SplinesToCurvilinearParameters splinesToCurvilinearParameters;
    splinesToCurvilinearParameters.AspectRatio = 0.1;
    splinesToCurvilinearParameters.AspectRatioGrowFactor = 1.1;
    splinesToCurvilinearParameters.AverageWidth = 500.0;
    splinesToCurvilinearParameters.CurvatureAdaptedGridSpacing = 1;
    splinesToCurvilinearParameters.GrowGridOutside = 1;
    splinesToCurvilinearParameters.MaximumNumberOfGridCellsInTheUniformPart = 5;
    splinesToCurvilinearParameters.GridsOnTopOfEachOtherTolerance = 0.0001;
    splinesToCurvilinearParameters.MinimumCosineOfCrossingAngles = 0.95;
    splinesToCurvilinearParameters.CheckFrontCollisions = 0;
    splinesToCurvilinearParameters.UniformGridSize = 0.0;
    splinesToCurvilinearParameters.DeleteSkinnyTriangles = 1;
    splinesToCurvilinearParameters.GrowGridOutside = 0;

    // Execute
    auto errorCode = mkernel_curvilinear_mesh_from_splines_ortho_initialize(0, geometryListIn, curvilinearParameters, splinesToCurvilinearParameters);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Grow grid, from the second layer
    for (auto layer = 1; layer <= curvilinearParameters.NRefinement; ++layer)
    {
        errorCode = meshkernelapi::mkernel_curvilinear_mesh_from_splines_ortho_iteration(0, layer);
        ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    }

    // Puts the computed curvilinear mesh into the mesh state (unstructured mesh)
    errorCode = meshkernelapi::mkernel_curvilinear_mesh_from_splines_ortho_refresh_mesh(0);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Delete the mesh curvilinearGridFromSplinesInstances vector entry
    errorCode = meshkernelapi::mkernel_curvilinear_mesh_from_splines_ortho_delete(0);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state
    meshkernelapi::MeshGeometryDimensions meshGeometryDimensions{};
    meshkernelapi::MeshGeometry meshGeometry{};
    errorCode = mkernel_get_mesh2d(0, meshGeometryDimensions, meshGeometry);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(21, meshGeometryDimensions.numnode);
    ASSERT_EQ(32, meshGeometryDimensions.numedge);

    // Delete dynamically allocated memory with operator new
    delete[] geometryListIn.xCoordinates;
    delete[] geometryListIn.yCoordinates;
    delete[] geometryListIn.zCoordinates;
}

TEST_F(ApiTests, GenerateTriangularGridThroughApi)
{
    // Prepare
    int meshKernelId;
    AllocateNewMesh(meshKernelId);

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometrySeparator = meshkernel::doubleMissingValue;
    geometryListIn.xCoordinates = new double[]{
        415.319672,
        390.271973,
        382.330048,
        392.715668,
        418.374268,
        453.807556,
        495.960968,
        532.005188,
        565.605774,
        590.653442,
        598.595398,
        593.708008,
        564.994812,
        514.899475,
        461.138611,
        422.039764,
        415.319672};

    geometryListIn.yCoordinates = new double[]{
        490.293762,
        464.024139,
        438.365448,
        411.484894,
        386.437103,
        366.276703,
        363.222107,
        370.553162,
        386.437103,
        412.095825,
        445.085571,
        481.129944,
        497.624817,
        504.955872,
        501.290344,
        493.348358,
        490.293762};

    geometryListIn.zCoordinates = new double[]{
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0};

    geometryListIn.numberOfCoordinates = 17;

    // Execute
    auto errorCode = mkernel_make_mesh_from_polygon(0, geometryListIn);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    // Get the new state
    meshkernelapi::MeshGeometryDimensions meshGeometryDimensions{};
    meshkernelapi::MeshGeometry meshGeometry{};
    errorCode = mkernel_get_mesh2d(0, meshGeometryDimensions, meshGeometry);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(38, meshGeometryDimensions.numnode);
    ASSERT_EQ(95, meshGeometryDimensions.numedge);

    // Delete dynamically allocated memory with operator new
    delete[] geometryListIn.xCoordinates;
    delete[] geometryListIn.yCoordinates;
    delete[] geometryListIn.zCoordinates;
}

TEST_F(ApiTests, GenerateTriangularGridFromSamplesThroughApi)
{
    // Prepare
    int meshKernelId;
    AllocateNewMesh(meshKernelId);

    meshkernelapi::GeometryList geometryListIn;

    geometryListIn.geometrySeparator = meshkernel::doubleMissingValue;

    geometryListIn.xCoordinates = new double[]{
        0.0,
        10.0,
        10.0,
        0.0,
        0.0};

    geometryListIn.yCoordinates = new double[]{
        0.0,
        0.0,
        10.0,
        10.0,
        0.0};

    geometryListIn.zCoordinates = new double[]{
        0.0,
        0.0,
        0.0,
        0.0,
        0.0};

    geometryListIn.numberOfCoordinates = 5;

    // Execute
    auto errorCode = mkernel_make_mesh_from_samples(0, geometryListIn);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state
    meshkernelapi::MeshGeometryDimensions meshGeometryDimensions{};
    meshkernelapi::MeshGeometry meshGeometry{};
    errorCode = mkernel_get_mesh2d(0, meshGeometryDimensions, meshGeometry);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(4, meshGeometryDimensions.numnode);
    ASSERT_EQ(5, meshGeometryDimensions.numedge);

    // Delete dynamically allocated memory with operator new
    delete[] geometryListIn.xCoordinates;
    delete[] geometryListIn.yCoordinates;
    delete[] geometryListIn.zCoordinates;
}

TEST_F(ApiTests, GetMeshBoundariesThroughApi)
{
    // Prepare
    MakeMesh();
    int numberOfpolygonNodes;
    auto errorCode = meshkernelapi::mkernel_copy_mesh_boundaries_to_polygon_count_nodes(0, numberOfpolygonNodes);
    ASSERT_EQ(11, numberOfpolygonNodes);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometrySeparator = meshkernel::doubleMissingValue;
    geometryListIn.numberOfCoordinates = numberOfpolygonNodes;
    geometryListIn.xCoordinates = new double[numberOfpolygonNodes];
    geometryListIn.yCoordinates = new double[numberOfpolygonNodes];
    geometryListIn.zCoordinates = new double[numberOfpolygonNodes];

    // Execute
    errorCode = mkernel_copy_mesh_boundaries_to_polygon(0, geometryListIn);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    const double tolerance = 1e-6;
    ASSERT_NEAR(0.0, geometryListIn.xCoordinates[0], tolerance);
    ASSERT_NEAR(0.0, geometryListIn.yCoordinates[0], tolerance);

    // Delete dynamically allocated memory with operator new
    delete[] geometryListIn.xCoordinates;
    delete[] geometryListIn.yCoordinates;
    delete[] geometryListIn.zCoordinates;
}

TEST_F(ApiTests, OffsetAPolygonThroughApi)
{
    // Prepare
    MakeMesh();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometrySeparator = meshkernel::doubleMissingValue;
    geometryListIn.numberOfCoordinates = 4;
    geometryListIn.xCoordinates = new double[]{
        0.0,
        1.0,
        1.0,
        0.0};

    geometryListIn.yCoordinates = new double[]{
        0.0,
        0.0,
        1.0,
        1.0};

    geometryListIn.zCoordinates = new double[]{
        0.0,
        0.0,
        0.0,
        0.0};

    // Execute
    int numberOfpolygonNodes;
    auto errorCode = mkernel_offsetted_polygon_count(0, geometryListIn, false, 0.5, numberOfpolygonNodes);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(4, numberOfpolygonNodes);

    meshkernelapi::GeometryList geometryListOut;

    geometryListOut.numberOfCoordinates = numberOfpolygonNodes;
    geometryListOut.geometrySeparator = meshkernel::doubleMissingValue;
    geometryListOut.xCoordinates = new double[numberOfpolygonNodes];
    geometryListOut.yCoordinates = new double[numberOfpolygonNodes];
    geometryListOut.zCoordinates = new double[numberOfpolygonNodes];
    errorCode = mkernel_offsetted_polygon(0, geometryListIn, false, 10.0, geometryListOut);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    const double tolerance = 1e-6;
    ASSERT_NEAR(0.0, geometryListOut.xCoordinates[0], tolerance);
    ASSERT_NEAR(-10.0, geometryListOut.yCoordinates[0], tolerance);

    // Delete dynamically allocated memory with operator new
    delete[] geometryListIn.xCoordinates;
    delete[] geometryListIn.yCoordinates;
    delete[] geometryListIn.zCoordinates;

    delete[] geometryListOut.xCoordinates;
    delete[] geometryListOut.yCoordinates;
    delete[] geometryListOut.zCoordinates;
}

TEST_F(ApiTests, RefineAPolygonThroughApi)
{
    // Prepare
    MakeMesh();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometrySeparator = meshkernel::doubleMissingValue;
    geometryListIn.numberOfCoordinates = 3;
    geometryListIn.xCoordinates = new double[]{
        76.251099,
        498.503723,
        505.253784};

    geometryListIn.yCoordinates = new double[]{
        92.626556,
        91.126541,
        490.130554};

    geometryListIn.zCoordinates = new double[]{
        0.0,
        0.0,
        0.0};

    // Execute
    int numberOfpolygonNodes;
    auto errorCode = mkernel_refine_polygon_count(0, geometryListIn, 0, 2, 40, numberOfpolygonNodes);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(22, numberOfpolygonNodes);

    meshkernelapi::GeometryList geometryListOut;
    geometryListOut.numberOfCoordinates = numberOfpolygonNodes;
    geometryListOut.geometrySeparator = meshkernel::doubleMissingValue;
    geometryListOut.xCoordinates = new double[numberOfpolygonNodes];
    geometryListOut.yCoordinates = new double[numberOfpolygonNodes];
    geometryListOut.zCoordinates = new double[numberOfpolygonNodes];
    errorCode = mkernel_refine_polygon(0, geometryListIn, false, 0, 2, geometryListOut);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    const double tolerance = 1e-6;
    ASSERT_NEAR(76.251099, geometryListOut.xCoordinates[0], tolerance);
    ASSERT_NEAR(92.626556, geometryListOut.yCoordinates[0], tolerance);

    // Delete dynamically allocated memory with operator new
    delete[] geometryListIn.xCoordinates;
    delete[] geometryListIn.yCoordinates;
    delete[] geometryListIn.zCoordinates;

    delete[] geometryListOut.xCoordinates;
    delete[] geometryListOut.yCoordinates;
    delete[] geometryListOut.zCoordinates;
}

TEST_F(ApiTests, RefineAGridBasedOnSamplesThroughApi)
{
    // Prepare
    MakeMesh();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometrySeparator = meshkernel::doubleMissingValue;

    geometryListIn.xCoordinates = new double[]{
        50.0,
        150.0,
        250.0,
        50.0,
        150.0,
        250.0,
        50.0,
        150.0,
        250.0};

    geometryListIn.yCoordinates = new double[]{
        50.0,
        50.0,
        50.0,
        150.0,
        150.0,
        150.0,
        250.0,
        250.0,
        250.0};

    geometryListIn.zCoordinates = new double[]{
        2.0,
        2.0,
        2.0,
        3.0,
        3.0,
        3.0,
        4.0,
        4.0,
        4.0};

    geometryListIn.numberOfCoordinates = 9;

    meshkernelapi::InterpolationParameters interpolationParameters;
    interpolationParameters.InterpolationType = 1;
    interpolationParameters.DisplayInterpolationProcess = 0;
    interpolationParameters.MaxNumberOfRefinementIterations = 2;
    interpolationParameters.AveragingMethod = 1;
    interpolationParameters.MinimumNumberOfPoints = 1;
    interpolationParameters.RelativeSearchRadius = 1.01;
    interpolationParameters.InterpolateTo = 3;
    interpolationParameters.RefineIntersected = 0;

    meshkernelapi::SampleRefineParameters samplesRefineParameters;
    samplesRefineParameters.SampleVectorDimension = 1;
    samplesRefineParameters.MinimumCellSize = 0.5;
    samplesRefineParameters.DirectionalRefinement = 0;
    samplesRefineParameters.RefinementType = 3;
    samplesRefineParameters.ConnectHangingNodes = 1;
    samplesRefineParameters.MaximumTimeStepInCourantGrid = 0.0;
    samplesRefineParameters.AccountForSamplesOutside = 0;

    // Execute
    auto errorCode = mkernel_refine_mesh_based_on_samples(0, geometryListIn, interpolationParameters, samplesRefineParameters);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state
    meshkernelapi::MeshGeometryDimensions meshGeometryDimensions{};
    meshkernelapi::MeshGeometry meshGeometry{};
    errorCode = mkernel_get_mesh2d(0, meshGeometryDimensions, meshGeometry);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(12, meshGeometryDimensions.numnode);
    ASSERT_EQ(17, meshGeometryDimensions.numedge);

    // Delete dynamically allocated memory with operator new
    delete[] geometryListIn.xCoordinates;
    delete[] geometryListIn.yCoordinates;
    delete[] geometryListIn.zCoordinates;
}

TEST_F(ApiTests, RefineAGridBasedOnPolygonThroughApi)
{
    // Prepare
    MakeMesh();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometrySeparator = meshkernel::doubleMissingValue;
    geometryListIn.xCoordinates = new double[]{
        50.0,
        150.0,
        250.0,
        50.0,
        150.0,
        250.0,
        50.0,
        150.0,
        250.0};

    geometryListIn.yCoordinates = new double[]{
        50.0,
        50.0,
        50.0,
        150.0,
        150.0,
        150.0,
        250.0,
        250.0,
        250.0};

    geometryListIn.zCoordinates = new double[]{
        2.0,
        2.0,
        2.0,
        3.0,
        3.0,
        3.0,
        4.0,
        4.0,
        4.0};

    geometryListIn.numberOfCoordinates = 9;

    meshkernelapi::InterpolationParameters interpolationParameters;
    interpolationParameters.InterpolationType = 1;
    interpolationParameters.DisplayInterpolationProcess = 0;
    interpolationParameters.MaxNumberOfRefinementIterations = 2;
    interpolationParameters.AveragingMethod = 1;
    interpolationParameters.MinimumNumberOfPoints = 1;
    interpolationParameters.RelativeSearchRadius = 1.01;
    interpolationParameters.InterpolateTo = 3;
    interpolationParameters.RefineIntersected = 0;

    // Execute
    auto errorCode = mkernel_refine_mesh_based_on_polygon(0, geometryListIn, interpolationParameters);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state
    meshkernelapi::MeshGeometryDimensions meshGeometryDimensions{};
    meshkernelapi::MeshGeometry meshGeometry{};
    errorCode = mkernel_get_mesh2d(0, meshGeometryDimensions, meshGeometry);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(12, meshGeometryDimensions.numnode);
    ASSERT_EQ(17, meshGeometryDimensions.numedge);

    // Delete dynamically allocated memory with operator new
    delete[] geometryListIn.xCoordinates;
    delete[] geometryListIn.yCoordinates;
    delete[] geometryListIn.zCoordinates;
}

TEST_F(ApiTests, MakeCurvilinearGridFromPolygonThroughApi)
{
    // Prepare
    MakeMesh();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometrySeparator = meshkernel::doubleMissingValue;
    geometryListIn.xCoordinates = new double[]{
        273.502319,
        274.252319,
        275.002350,
        458.003479,
        719.005127,
        741.505249,
        710.755066,
        507.503784,
        305.002533};

    geometryListIn.yCoordinates = new double[]{
        478.880432,
        325.128906,
        172.127350,
        157.127213,
        157.127213,
        328.128937,
        490.880554,
        494.630615,
        493.130615};

    geometryListIn.zCoordinates = new double[]{
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0};
    geometryListIn.numberOfCoordinates = 9;

    // Execute
    auto errorCode = meshkernelapi::mkernel_curvilinear_from_polygon(0, geometryListIn, 0, 2, 4, true);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state
    meshkernelapi::MeshGeometryDimensions meshGeometryDimensions{};
    meshkernelapi::MeshGeometry meshGeometry{};
    errorCode = mkernel_get_mesh2d(0, meshGeometryDimensions, meshGeometry);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(21, meshGeometryDimensions.numnode);
    ASSERT_EQ(29, meshGeometryDimensions.numedge);

    // Delete dynamically allocated memory with operator new
    delete[] geometryListIn.xCoordinates;
    delete[] geometryListIn.yCoordinates;
    delete[] geometryListIn.zCoordinates;
}

TEST_F(ApiTests, GetClosestMeshCoordinateThroughApi)
{
    // Prepare
    MakeMesh();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometrySeparator = meshkernel::doubleMissingValue;
    geometryListIn.xCoordinates = new double[]{-5.0};
    geometryListIn.yCoordinates = new double[]{5.0};
    geometryListIn.zCoordinates = new double[]{0.0};
    geometryListIn.numberOfCoordinates = 1;

    meshkernelapi::GeometryList geometryListOut;
    geometryListOut.geometrySeparator = meshkernel::doubleMissingValue;
    geometryListOut.xCoordinates = new double[]{meshkernel::doubleMissingValue};
    geometryListOut.yCoordinates = new double[]{meshkernel::doubleMissingValue};
    geometryListOut.zCoordinates = new double[]{meshkernel::doubleMissingValue};
    geometryListOut.numberOfCoordinates = 1;

    // Execute
    auto errorCode = mkernel_get_node_coordinate(0, geometryListIn, 10.0, geometryListOut);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(0.0, geometryListOut.xCoordinates[0]);

    // Delete dynamically allocated memory with operator new
    delete[] geometryListIn.xCoordinates;
    delete[] geometryListIn.yCoordinates;
    delete[] geometryListIn.zCoordinates;

    delete[] geometryListOut.xCoordinates;
    delete[] geometryListOut.yCoordinates;
    delete[] geometryListOut.zCoordinates;
}

TEST_F(ApiTests, MakeCurvilinearGridFromTriangleThroughApi)
{
    // Prepare
    MakeMesh();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometrySeparator = meshkernel::doubleMissingValue;
    geometryListIn.xCoordinates = new double[]{
        444.504791,
        427.731781,
        405.640503,
        381.094666,
        451.050354,
        528.778931,
        593.416260,
        558.643005,
        526.733398,
        444.095703};

    geometryListIn.yCoordinates = new double[]{
        437.155945,
        382.745758,
        317.699005,
        262.470612,
        262.879700,
        263.288788,
        266.561584,
        324.653687,
        377.836578,
        436.746857};

    geometryListIn.zCoordinates = new double[]{
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0};

    geometryListIn.numberOfCoordinates = 10;

    // Execute
    auto errorCode = mkernel_curvilinear_from_triangle(0, geometryListIn, 0, 3, 6);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state
    meshkernelapi::MeshGeometryDimensions meshGeometryDimensions{};
    meshkernelapi::MeshGeometry meshGeometry{};
    errorCode = mkernel_get_mesh2d(0, meshGeometryDimensions, meshGeometry);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(28, meshGeometryDimensions.numnode);
    ASSERT_EQ(40, meshGeometryDimensions.numedge);

    // Delete dynamically allocated memory with operator new
    delete[] geometryListIn.xCoordinates;
    delete[] geometryListIn.yCoordinates;
    delete[] geometryListIn.zCoordinates;
}

TEST_F(ApiTests, ComputeSingleContactsThroughApi)
{
    // Prepare
    MakeMesh(4, 4, 10);

    // Init 1d mesh
    meshkernelapi::Mesh1D mesh1d;
    std::unique_ptr<double> nodex(new double[]{
        1.73493900000000,
        2.35659313023165,
        5.38347452702839,
        14.2980910429074,
        22.9324017677239,
        25.3723169493137,
        25.8072280000000});
    std::unique_ptr<double> nodey(new double[]{
        -7.6626510000000,
        1.67281447902331,
        10.3513746546384,
        12.4797224193970,
        15.3007317677239,
        24.1623588554512,
        33.5111870000000});
    mesh1d.nodex = nodex.get();
    mesh1d.nodey = nodey.get();
    mesh1d.num_nodes = 7;

    std::unique_ptr<int> edge_nodes(new int[]{
        0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6});
    mesh1d.edge_nodes = edge_nodes.get();
    mesh1d.num_edges = 6;

    auto errorCode = mkernel_set_mesh1d(0, mesh1d, false);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Init 1d mask
    std::unique_ptr<int> onedNodeMask(new int[]{
        1, 1, 1, 1, 1, 1, 1});

    // Init polygon
    meshkernelapi::GeometryList polygon;
    polygon.geometrySeparator = meshkernel::doubleMissingValue;

    std::unique_ptr<double> xCoordinates(new double[]{
        -30, 40, 40, -40, -30});
    std::unique_ptr<double> yCoordinates(new double[]{
        -20, -20, 50, 50, -20});
    std::unique_ptr<double> zCoordinates(new double[]{
        0, 0, 0, 0, 0});
    polygon.xCoordinates = xCoordinates.get();
    polygon.yCoordinates = yCoordinates.get();
    polygon.zCoordinates = zCoordinates.get();

    polygon.numberOfCoordinates = 5;

    // Execute
    errorCode = mkernel_compute_single_contacts(0, onedNodeMask.get(), polygon);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state
    meshkernelapi::Contacts contacts{};
    errorCode = mkernel_get_contacts_size(0, contacts);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    std::unique_ptr<int> mesh1d_indices(new int[contacts.num_contacts]);
    std::unique_ptr<int> mesh2d_indices(new int[contacts.num_contacts]);
    contacts.mesh1d_indices = mesh1d_indices.get();
    contacts.mesh2d_indices = mesh2d_indices.get();

    errorCode = mkernel_get_contacts_data(0, contacts);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(5, contacts.num_contacts);

    ASSERT_EQ(1, contacts.mesh1d_indices[0]);
    ASSERT_EQ(2, contacts.mesh1d_indices[1]);
    ASSERT_EQ(3, contacts.mesh1d_indices[2]);
    ASSERT_EQ(4, contacts.mesh1d_indices[3]);
    ASSERT_EQ(5, contacts.mesh1d_indices[4]);

    ASSERT_EQ(0, contacts.mesh2d_indices[0]);
    ASSERT_EQ(1, contacts.mesh2d_indices[1]);
    ASSERT_EQ(4, contacts.mesh2d_indices[2]);
    ASSERT_EQ(7, contacts.mesh2d_indices[3]);
    ASSERT_EQ(8, contacts.mesh2d_indices[4]);
}

TEST_F(ApiTests, ComputeMultipleContactsThroughApi)
{
    // Prepare
    MakeMesh(4, 4, 10);

    // Init 1d mesh
    meshkernelapi::Mesh1D mesh1d;
    std::unique_ptr<double> nodex(new double[]{
        1.73493900000000,
        2.35659313023165,
        5.38347452702839,
        14.2980910429074,
        22.9324017677239,
        25.3723169493137,
        25.8072280000000});
    std::unique_ptr<double> nodey(new double[]{
        -7.6626510000000,
        1.67281447902331,
        10.3513746546384,
        12.4797224193970,
        15.3007317677239,
        24.1623588554512,
        33.5111870000000});
    mesh1d.nodex = nodex.get();
    mesh1d.nodey = nodey.get();
    mesh1d.num_nodes = 7;

    std::unique_ptr<int> edge_nodes(new int[]{
        0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6});
    mesh1d.edge_nodes = edge_nodes.get();
    mesh1d.num_edges = 6;

    auto errorCode = mkernel_set_mesh1d(0, mesh1d, false);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Init 1d mask
    std::unique_ptr<int> onedNodeMask(new int[]{
        1, 1, 1, 1, 1, 1, 1});

    // Execute
    // Why do we need to specify the namespace "meshkernelapi"?
    errorCode = meshkernelapi::mkernel_compute_multiple_contacts(0, onedNodeMask.get());
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state
    meshkernelapi::Contacts contacts{};
    errorCode = mkernel_get_contacts_size(0, contacts);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    std::unique_ptr<int> mesh1d_indices(new int[contacts.num_contacts]);
    std::unique_ptr<int> mesh2d_indices(new int[contacts.num_contacts]);
    contacts.mesh1d_indices = mesh1d_indices.get();
    contacts.mesh2d_indices = mesh2d_indices.get();

    errorCode = mkernel_get_contacts_data(0, contacts);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(5, contacts.num_contacts);

    ASSERT_EQ(1, contacts.mesh1d_indices[0]);
    ASSERT_EQ(2, contacts.mesh1d_indices[1]);
    ASSERT_EQ(3, contacts.mesh1d_indices[2]);
    ASSERT_EQ(4, contacts.mesh1d_indices[3]);
    ASSERT_EQ(5, contacts.mesh1d_indices[4]);

    ASSERT_EQ(0, contacts.mesh2d_indices[0]);
    ASSERT_EQ(1, contacts.mesh2d_indices[1]);
    ASSERT_EQ(4, contacts.mesh2d_indices[2]);
    ASSERT_EQ(7, contacts.mesh2d_indices[3]);
    ASSERT_EQ(8, contacts.mesh2d_indices[4]);
}

TEST_F(ApiTests, ComputeContactsWithPolygonsThroughApi)
{
    // Prepare
    MakeMesh(4, 4, 10);

    // Init 1d mesh
    meshkernelapi::Mesh1D mesh1d;
    std::unique_ptr<double> nodex(new double[]{
        1.73493900000000,
        2.35659313023165,
        5.38347452702839,
        14.2980910429074,
        22.9324017677239,
        25.3723169493137,
        25.8072280000000});
    std::unique_ptr<double> nodey(new double[]{
        -7.6626510000000,
        1.67281447902331,
        10.3513746546384,
        12.4797224193970,
        15.3007317677239,
        24.1623588554512,
        33.5111870000000});
    mesh1d.nodex = nodex.get();
    mesh1d.nodey = nodey.get();
    mesh1d.num_nodes = 7;

    std::unique_ptr<int> edge_nodes(new int[]{
        0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6});
    mesh1d.edge_nodes = edge_nodes.get();
    mesh1d.num_edges = 6;

    auto errorCode = mkernel_set_mesh1d(0, mesh1d, false);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Init 1d mask
    std::unique_ptr<int> onedNodeMask(new int[]{
        1, 1, 1, 1, 1, 1, 1});

    // Init polygon
    meshkernelapi::GeometryList polygon;
    polygon.geometrySeparator = meshkernel::doubleMissingValue;

    std::unique_ptr<double> xCoordinates(new double[]{
        25,
        50,
        50,
        25,
        25});
    std::unique_ptr<double> yCoordinates(new double[]{
        25,
        25,
        50,
        50,
        25});
    std::unique_ptr<double> zCoordinates(new double[]{
        0, 0, 0, 0, 0});
    polygon.xCoordinates = xCoordinates.get();
    polygon.yCoordinates = yCoordinates.get();
    polygon.zCoordinates = zCoordinates.get();

    polygon.numberOfCoordinates = 5;

    // Execute
    errorCode = mkernel_compute_contacts_with_polygons(0, onedNodeMask.get(), polygon);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state
    meshkernelapi::Contacts contacts{};
    errorCode = mkernel_get_contacts_size(0, contacts);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    std::unique_ptr<int> mesh1d_indices(new int[contacts.num_contacts]);
    std::unique_ptr<int> mesh2d_indices(new int[contacts.num_contacts]);
    contacts.mesh1d_indices = mesh1d_indices.get();
    contacts.mesh2d_indices = mesh2d_indices.get();

    errorCode = mkernel_get_contacts_data(0, contacts);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(1, contacts.num_contacts);
    ASSERT_EQ(5, contacts.mesh1d_indices[0]);
    ASSERT_EQ(8, contacts.mesh2d_indices[0]);
}

TEST_F(ApiTests, ComputeContactsWithPointsThroughApi)
{
    // Prepare
    MakeMesh(4, 4, 10);

    // Init 1d mesh
    meshkernelapi::Mesh1D mesh1d;
    std::unique_ptr<double> nodex(new double[]{
        1.73493900000000,
        2.35659313023165,
        5.38347452702839,
        14.2980910429074,
        22.9324017677239,
        25.3723169493137,
        25.8072280000000});
    std::unique_ptr<double> nodey(new double[]{
        -7.6626510000000,
        1.67281447902331,
        10.3513746546384,
        12.4797224193970,
        15.3007317677239,
        24.1623588554512,
        33.5111870000000});
    mesh1d.nodex = nodex.get();
    mesh1d.nodey = nodey.get();
    mesh1d.num_nodes = 7;

    std::unique_ptr<int> edge_nodes(new int[]{
        0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6});
    mesh1d.edge_nodes = edge_nodes.get();
    mesh1d.num_edges = 6;

    auto errorCode = mkernel_set_mesh1d(0, mesh1d, false);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Init 1d mask
    std::unique_ptr<int> onedNodeMask(new int[]{
        1, 1, 1, 1, 1, 1, 1});

    // Init polygon
    meshkernelapi::GeometryList points;
    points.geometrySeparator = meshkernel::doubleMissingValue;

    std::unique_ptr<double> xCoordinates(new double[]{
        5,
        15,
        25,
        35});
    std::unique_ptr<double> yCoordinates(new double[]{
        5,
        15,
        25,
        35});
    std::unique_ptr<double> zCoordinates(new double[]{
        0, 0, 0, 0});
    points.xCoordinates = xCoordinates.get();
    points.yCoordinates = yCoordinates.get();
    points.zCoordinates = zCoordinates.get();

    points.numberOfCoordinates = 4;

    // Execute
    errorCode = mkernel_compute_contacts_with_points(0, onedNodeMask.get(), points);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state
    meshkernelapi::Contacts contacts{};
    errorCode = mkernel_get_contacts_size(0, contacts);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    std::unique_ptr<int> mesh1d_indices(new int[contacts.num_contacts]);
    std::unique_ptr<int> mesh2d_indices(new int[contacts.num_contacts]);
    contacts.mesh1d_indices = mesh1d_indices.get();
    contacts.mesh2d_indices = mesh2d_indices.get();

    errorCode = mkernel_get_contacts_data(0, contacts);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(3, contacts.num_contacts);

    ASSERT_EQ(1, contacts.mesh1d_indices[0]);
    ASSERT_EQ(3, contacts.mesh1d_indices[1]);
    ASSERT_EQ(5, contacts.mesh1d_indices[2]);

    ASSERT_EQ(0, contacts.mesh2d_indices[0]);
    ASSERT_EQ(4, contacts.mesh2d_indices[1]);
    ASSERT_EQ(8, contacts.mesh2d_indices[2]);
}

TEST_F(ApiTests, ComputeBoundaryContactsThroughApi)
{
    // Prepare
    MakeMesh(4, 4, 10);

    // Init 1d mesh
    meshkernelapi::Mesh1D mesh1d;
    std::unique_ptr<double> nodex(new double[]{
        -16.1886410000000,
        -16.1464995876014,
        -16.1043581752028,
        -16.0622167628042,
        -15.7539488236928,
        -6.86476658679268,
        2.02441565010741,
        10.9135970000000,
    });
    std::unique_ptr<double> nodey(new double[]{
        0.89018900000000,
        9.78201442138723,
        18.6738398427745,
        27.5656652641617,
        36.1966603330179,
        36.4175095626911,
        36.6383587923643,
        36.8592080000000});
    mesh1d.nodex = nodex.get();
    mesh1d.nodey = nodey.get();
    mesh1d.num_nodes = 8;

    std::unique_ptr<int> edge_nodes(new int[]{
        0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7});
    mesh1d.edge_nodes = edge_nodes.get();
    mesh1d.num_edges = 7;

    auto errorCode = mkernel_set_mesh1d(0, mesh1d, false);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Init 1d mask
    std::unique_ptr<int> onedNodeMask(new int[]{
        1, 1, 1, 1, 1, 1, 1, 1});

    // Init polygon
    meshkernelapi::GeometryList polygon;
    polygon.geometrySeparator = meshkernel::doubleMissingValue;

    std::unique_ptr<double> xCoordinates(new double[]{
        -30,
        40,
        40,
        -40,
        -30});
    std::unique_ptr<double> yCoordinates(new double[]{
        -20,
        -20,
        50,
        50,
        -20});
    std::unique_ptr<double> zCoordinates(new double[]{
        0, 0, 0, 0, 0});
    polygon.xCoordinates = xCoordinates.get();
    polygon.yCoordinates = yCoordinates.get();
    polygon.zCoordinates = zCoordinates.get();

    polygon.numberOfCoordinates = 5;

    // Execute
    errorCode = mkernel_compute_boundary_contacts(0, onedNodeMask.get(), polygon, 200.0);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state
    meshkernelapi::Contacts contacts{};
    errorCode = mkernel_get_contacts_size(0, contacts);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    std::unique_ptr<int> mesh1d_indices(new int[contacts.num_contacts]);
    std::unique_ptr<int> mesh2d_indices(new int[contacts.num_contacts]);
    contacts.mesh1d_indices = mesh1d_indices.get();
    contacts.mesh2d_indices = mesh2d_indices.get();

    errorCode = mkernel_get_contacts_data(0, contacts);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(8, contacts.num_contacts);

    ASSERT_EQ(0, contacts.mesh1d_indices[0]);
    ASSERT_EQ(2, contacts.mesh1d_indices[1]);
    ASSERT_EQ(6, contacts.mesh1d_indices[2]);
    ASSERT_EQ(0, contacts.mesh1d_indices[3]);
    ASSERT_EQ(7, contacts.mesh1d_indices[4]);
    ASSERT_EQ(7, contacts.mesh1d_indices[5]);
    ASSERT_EQ(7, contacts.mesh1d_indices[6]);
    ASSERT_EQ(7, contacts.mesh1d_indices[7]);

    ASSERT_EQ(0, contacts.mesh2d_indices[0]);
    ASSERT_EQ(1, contacts.mesh2d_indices[1]);
    ASSERT_EQ(2, contacts.mesh2d_indices[2]);
    ASSERT_EQ(3, contacts.mesh2d_indices[3]);
    ASSERT_EQ(5, contacts.mesh2d_indices[4]);
    ASSERT_EQ(6, contacts.mesh2d_indices[5]);
    ASSERT_EQ(7, contacts.mesh2d_indices[6]);
    ASSERT_EQ(8, contacts.mesh2d_indices[7]);
}

TEST(ApiStatelessTests, GetSplinesThroughApi)
{
    // Prepare
    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.xCoordinates = new double[]{10.0, 20.0, 30.0};
    geometryListIn.yCoordinates = new double[]{-5.0, 5.0, -5.0};
    geometryListIn.zCoordinates = new double[]{0.0, 0.0, 0.0};
    geometryListIn.numberOfCoordinates = 3;
    geometryListIn.geometrySeparator = meshkernel::doubleMissingValue;

    meshkernelapi::GeometryList geometryListOut;
    int numberOfPointsBetweenNodes = 20;
    geometryListOut.xCoordinates = new double[(numberOfPointsBetweenNodes + 1) * 2 + 1];
    geometryListOut.yCoordinates = new double[(numberOfPointsBetweenNodes + 1) * 2 + 1];
    geometryListOut.zCoordinates = new double[(numberOfPointsBetweenNodes + 1) * 2 + 1];

    // Execute
    auto errorCode = mkernel_get_splines(geometryListIn, geometryListOut, numberOfPointsBetweenNodes);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ((numberOfPointsBetweenNodes + 1) * 2, geometryListOut.numberOfCoordinates);

    // Delete dynamically allocated memory with operator new
    delete[] geometryListIn.xCoordinates;
    delete[] geometryListIn.yCoordinates;
    delete[] geometryListIn.zCoordinates;

    delete[] geometryListOut.xCoordinates;
    delete[] geometryListOut.yCoordinates;
    delete[] geometryListOut.zCoordinates;
}

TEST(ApiStatelessTests, OrthogonalizingAnInvaliMeshShouldThrowAMeshGeometryError)
{
    // Prepare
    int meshKernelId;
    meshkernelapi::mkernel_new_mesh(meshKernelId);

    auto meshData = ReadLegacyMeshFromFileForApiTesting("../../../../tests/data/InvalidMeshes/invalid_orthogonalization_net.nc");
    meshkernelapi::Mesh1D mesh1d{};
    auto errorCode = mkernel_set_mesh2d(meshKernelId, std::get<1>(meshData), std::get<0>(meshData), false);
    DeleteRectangularMeshForApiTesting(std::get<0>(meshData));
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::OrthogonalizationParameters orthogonalizationParameters{};
    orthogonalizationParameters.OuterIterations = 1;
    orthogonalizationParameters.BoundaryIterations = 25;
    orthogonalizationParameters.InnerIterations = 25;
    orthogonalizationParameters.OrthogonalizationToSmoothingFactor = 0.975;

    meshkernelapi::GeometryList geometryList{};
    meshkernelapi::GeometryList landBoundaries{};

    // Execute
    errorCode = mkernel_orthogonalize_initialize(0,
                                                 1,
                                                 orthogonalizationParameters,
                                                 landBoundaries,
                                                 geometryList);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert there is a geometry error
    errorCode = meshkernelapi::mkernel_orthogonalize_prepare_outer_iteration(0);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::InvalidGeometry, errorCode);

    //Delete orthogonalization instance
    errorCode = meshkernelapi::mkernel_orthogonalize_delete(0);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the message
    const char* exceptionMessage;
    errorCode = meshkernelapi::mkernel_get_error(exceptionMessage);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the index of the invalid location
    int invalidIndex;
    int type;
    errorCode = meshkernelapi::mkernel_get_geometry_error(invalidIndex, type);
    ASSERT_EQ(static_cast<int>(meshkernel::MeshLocations::Nodes), type);
    ASSERT_EQ(478, invalidIndex);
}
