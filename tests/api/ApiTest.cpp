#include <gtest/gtest.h>

#include <MeshKernel/GeometryList.hpp>
#include <MeshKernel/MakeMeshParameters.hpp>
#include <MeshKernel/MeshGeometry.hpp>
#include <MeshKernel/MeshGeometryDimensions.hpp>
#include <MeshKernel/MeshKernel.hpp>
#include <TestUtils/MakeMeshes.hpp>

class ApiTest : public ::testing::Test
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

TEST_F(ApiTest, DeleteNodeThroughApi)
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

TEST_F(ApiTest, FlipEdgesThroughApi)
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

TEST_F(ApiTest, InsertEdgeThroughApi)
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

TEST_F(ApiTest, MergeTwoNodesThroughApi)
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

TEST_F(ApiTest, MergeNodesThroughApi)
{
    // Set a new mesh in Mesh Kernel
    MakeMesh();

    // Execute
    meshkernelapi::GeometryList geometry_list{};
    auto errorCode = meshkernelapi::mkernel_merge_nodes(0, geometry_list);
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

TEST_F(ApiTest, OrthogonalizationThroughApi)
{
    // Set a new mesh in Mesh Kernel
    MakeMesh();

    // Execute
    meshkernelapi::OrthogonalizationParameters orthogonalizationParameters{};
    meshkernelapi::GeometryList geometryList{};
    meshkernelapi::GeometryList landBoundaries{};

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

TEST_F(ApiTest, MakeGridThroughApi)
{
    // Allocate a new mesh entry
    int meshKernelId;
    AllocateNewMesh(meshKernelId);

    // Execute
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