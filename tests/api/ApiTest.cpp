#include <gtest/gtest.h>

#include <MeshKernel/Constants.hpp>
#include <MeshKernel/GeometryList.hpp>
#include <MeshKernel/MakeMeshParameters.hpp>
#include <MeshKernel/MeshGeometry.hpp>
#include <MeshKernel/MeshGeometryDimensions.hpp>
#include <MeshKernel/MeshKernel.hpp>
#include <TestUtils/MakeMeshes.hpp>

class ApiTests : public ::testing::Test
{
public:
    void AllocateNewMesh(int& meshKernelId) const
    {
        // Allocate new mesh
        meshkernelapi::mkernel_new_mesh(meshKernelId);
        ASSERT_EQ(0, meshKernelId);
    }

    void MakeMesh(int n = 3, int m = 3, double delta = 1.0) const
    {
        // Allocate new mesh
        int meshKernelId;
        AllocateNewMesh(meshKernelId);

        // Set-up new mesh
        auto mesh = MakeRectangularMeshForApiTesting(n, m, delta);
        auto errorCode = mkernel_set_state(meshKernelId, std::get<1>(mesh), std::get<0>(mesh), false);
        ASSERT_EQ(0, errorCode);
        // Clean up client-allocated memory
        DeleteRectangularMeshForApiTesting(std::get<0>(mesh));
    }

protected:
    void TearDown() override
    {
        //De-allocate meshkernel
        const auto errorCode = meshkernelapi::mkernel_deallocate_state(0);
        ASSERT_EQ(0, errorCode);
    }
};

TEST_F(ApiTests, DeleteNodeThroughApi)
{
    // Set a new mesh in Mesh Kernel
    MakeMesh();

    // Execute
    auto errorCode = meshkernelapi::mkernel_delete_node(0, 0);
    ASSERT_EQ(0, errorCode);

    // Get the new state
    meshkernelapi::MeshGeometryDimensions meshGeometryDimensions{};
    meshkernelapi::MeshGeometry meshGeometry{};
    errorCode = mkernel_get_mesh(0, meshGeometryDimensions, meshGeometry);

    // Assert
    ASSERT_EQ(0, errorCode);
    ASSERT_EQ(8, meshGeometryDimensions.numnode);
    ASSERT_EQ(10, meshGeometryDimensions.numedge);
}

TEST_F(ApiTests, FlipEdgesThroughApi)
{
    // Set a new mesh in Mesh Kernel
    MakeMesh();

    // Execute
    auto errorCode = meshkernelapi::mkernel_flip_edges(0, 1, 1);
    ASSERT_EQ(0, errorCode);

    // Get the new state
    meshkernelapi::MeshGeometryDimensions meshGeometryDimensions{};
    meshkernelapi::MeshGeometry meshGeometry{};
    errorCode = mkernel_get_mesh(0, meshGeometryDimensions, meshGeometry);

    // Assert
    ASSERT_EQ(0, errorCode);
    ASSERT_EQ(9, meshGeometryDimensions.numnode);
    ASSERT_EQ(16, meshGeometryDimensions.numedge);
}

TEST_F(ApiTests, InsertEdgeThroughApi)
{
    // Set a new mesh in Mesh Kernel
    MakeMesh();

    // Execute
    int newEdgeIndex;
    auto errorCode = meshkernelapi::mkernel_insert_edge(0, 0, 4, newEdgeIndex);
    ASSERT_EQ(0, errorCode);
    ASSERT_EQ(12, newEdgeIndex);

    // Get the new state
    meshkernelapi::MeshGeometryDimensions meshGeometryDimensions{};
    meshkernelapi::MeshGeometry meshGeometry{};
    errorCode = mkernel_get_mesh(0, meshGeometryDimensions, meshGeometry);

    // Assert
    ASSERT_EQ(0, errorCode);
    ASSERT_EQ(9, meshGeometryDimensions.numnode);
    ASSERT_EQ(13, meshGeometryDimensions.numedge);
}

TEST_F(ApiTests, MergeTwoNodesThroughApi)
{
    // Set a new mesh in Mesh Kernel
    MakeMesh();

    // Execute
    auto errorCode = meshkernelapi::mkernel_merge_two_nodes(0, 0, 4);
    ASSERT_EQ(0, errorCode);

    // Get the new state
    meshkernelapi::MeshGeometryDimensions meshGeometryDimensions{};
    meshkernelapi::MeshGeometry meshGeometry{};
    errorCode = mkernel_get_mesh(0, meshGeometryDimensions, meshGeometry);

    // Assert
    ASSERT_EQ(0, errorCode);
    ASSERT_EQ(8, meshGeometryDimensions.numnode);
    ASSERT_EQ(10, meshGeometryDimensions.numedge);
}

TEST_F(ApiTests, MergeNodesThroughApi)
{
    // Prepare
    MakeMesh();
    meshkernelapi::GeometryList geometry_list{};

    // Execute
    auto errorCode = mkernel_merge_nodes(0, geometry_list);
    ASSERT_EQ(0, errorCode);

    // Get the new state
    meshkernelapi::MeshGeometryDimensions meshGeometryDimensions{};
    meshkernelapi::MeshGeometry meshGeometry{};
    errorCode = mkernel_get_mesh(0, meshGeometryDimensions, meshGeometry);

    // Assert (nothing is done, we just check that the api communication works)
    ASSERT_EQ(0, errorCode);
    ASSERT_EQ(9, meshGeometryDimensions.numnode);
    ASSERT_EQ(12, meshGeometryDimensions.numedge);
}

TEST_F(ApiTests, OrthogonalizationThroughApi)
{
    // Set a new mesh in Mesh Kernel
    MakeMesh();

    // Prepare
    meshkernelapi::OrthogonalizationParameters orthogonalizationParameters{};
    meshkernelapi::GeometryList geometryList{};
    meshkernelapi::GeometryList landBoundaries{};

    // Execute
    auto errorCode = mkernel_orthogonalize_initialize(0,
                                                      1,
                                                      orthogonalizationParameters,
                                                      landBoundaries,
                                                      geometryList);
    ASSERT_EQ(0, errorCode);

    errorCode = meshkernelapi::mkernel_orthogonalize_prepare_outer_iteration(0);
    ASSERT_EQ(0, errorCode);

    errorCode = meshkernelapi::mkernel_orthogonalize_inner_iteration(0);
    ASSERT_EQ(0, errorCode);

    errorCode = meshkernelapi::mkernel_orthogonalize_finalize_outer_iteration(0);
    ASSERT_EQ(0, errorCode);

    // Get the new state
    meshkernelapi::MeshGeometryDimensions meshGeometryDimensions{};
    meshkernelapi::MeshGeometry meshGeometry{};
    errorCode = mkernel_get_mesh(0, meshGeometryDimensions, meshGeometry);

    // Assert (nothing is done, we just check that the api communication works)
    ASSERT_EQ(0, errorCode);
    ASSERT_EQ(9, meshGeometryDimensions.numnode);
    ASSERT_EQ(12, meshGeometryDimensions.numedge);
}

TEST_F(ApiTests, MakeGridThroughApi)
{
    // Allocate a new mesh entry
    int meshKernelId;
    AllocateNewMesh(meshKernelId);

    // Prepare
    meshkernelapi::MakeMeshParameters makeMeshParameters{};
    meshkernelapi::GeometryList geometryList{};

    makeMeshParameters.GridType = 0;
    makeMeshParameters.NumberOfColumns = 3;
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
    ASSERT_EQ(0, errorCode);

    // Get the new state
    meshkernelapi::MeshGeometryDimensions meshGeometryDimensions{};
    meshkernelapi::MeshGeometry meshGeometry{};
    errorCode = mkernel_get_mesh(0, meshGeometryDimensions, meshGeometry);

    // Assert (nothing is done, we just check that the api communication works)
    ASSERT_EQ(0, errorCode);
    ASSERT_EQ(16, meshGeometryDimensions.numnode);
    ASSERT_EQ(24, meshGeometryDimensions.numedge);
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
    ASSERT_EQ(0, errorCode);
    ASSERT_EQ((numberOfPointsBetweenNodes + 1) * 2, geometryListOut.numberOfCoordinates);

    // Delete dynamically allocated memory with operator new
    delete[] geometryListIn.xCoordinates;
    delete[] geometryListIn.yCoordinates;
    delete[] geometryListIn.zCoordinates;

    delete[] geometryListOut.xCoordinates;
    delete[] geometryListOut.yCoordinates;
    delete[] geometryListOut.zCoordinates;
}

TEST_F(ApiTests, GenerateCurvilinearGridThroughApi)
{
    // Allocate a new mesh entry
    int meshKernelId;
    AllocateNewMesh(meshKernelId);

    // Execute
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

    auto errorCode = mkernel_curvilinear_mesh_from_splines(0, geometryListIn, curvilinearParameters);
    ASSERT_EQ(0, errorCode);

    // Get the new state
    meshkernelapi::MeshGeometryDimensions meshGeometryDimensions{};
    meshkernelapi::MeshGeometry meshGeometry{};
    errorCode = mkernel_get_mesh(0, meshGeometryDimensions, meshGeometry);

    // Assert (nothing is done, we just check that the api communication works)
    ASSERT_EQ(0, errorCode);
    ASSERT_EQ(121, meshGeometryDimensions.numnode);
    ASSERT_EQ(220, meshGeometryDimensions.numedge);
}