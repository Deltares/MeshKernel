#include <gtest/gtest.h>

#include <MeshKernel/MeshGeometry.hpp>
#include <MeshKernel/MeshGeometryDimensions.hpp>
#include <MeshKernel/MeshKernel.hpp>
#include <TestUtils/MakeMeshes.hpp>

class ApiTest : public ::testing::Test
{
public:
    void MakeMesh(int n = 3, int m = 3, double delta = 1.0) const
    {
        // Allocate new mesh
        int meshKernelId;
        meshkernelapi::mkernel_new_mesh(meshKernelId);
        ASSERT_EQ(0, meshKernelId);

        // Set-up new mesh
        auto [meshGeometry, meshGeometryDimensions] = MakeRectangularMeshForApiTesting(n, m, delta);
        auto errorCode = mkernel_set_state(meshKernelId, meshGeometryDimensions, meshGeometry, false);
        ASSERT_EQ(0, errorCode);
        // Clean up client-allocated memory
        DeleteRectangularMeshForApiTesting(meshGeometry);
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