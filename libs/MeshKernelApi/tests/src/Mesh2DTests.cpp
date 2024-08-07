#include <algorithm>
#include <gtest/gtest.h>
#include <random>

#include "MeshKernel/Parameters.hpp"

#include "MeshKernelApi/BoundingBox.hpp"
#include "MeshKernelApi/Mesh2D.hpp"
#include "MeshKernelApi/MeshKernel.hpp"

#include "TestUtils/MakeMeshes.hpp"

// namespace aliases
namespace mk = meshkernel;
namespace mkapi = meshkernelapi;

TEST(Mesh2DTests, Mesh2dApiNodeEdgeDataTest)
{
    const int clgSize = 2;

    int meshKernelId = mkapi::mkernel_get_null_identifier();
    int errorCode;

    // Clear the meshkernel state and undo stack before starting the test.
    errorCode = mkapi::mkernel_clear_state();
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    errorCode = mkapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    meshkernel::MakeGridParameters makeGridParameters;

    // num_columns and num_rows indicate number of elements in each direction, so value = clgSize - 1
    makeGridParameters.num_columns = clgSize - 1;
    makeGridParameters.num_rows = clgSize - 1;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = 1.0;
    makeGridParameters.block_size_y = 1.0;

    // Generate curvilinear grid.
    errorCode = mkapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    // Convert the curvilinear grid to an unstructured grid.
    errorCode = mkapi::mkernel_curvilinear_convert_to_mesh2d(meshKernelId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    int node1;
    int node2;
    int node3;
    int node4;

    mkapi::BoundingBox boundingBox{-0.1, -0.1, 1.1, 1.1};

    errorCode = mkapi::mkernel_mesh2d_get_location_index(meshKernelId, 0.0, 0.0,
                                                         1 /*Location::Node*/,
                                                         boundingBox,
                                                         node1);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    errorCode = mkapi::mkernel_mesh2d_get_location_index(meshKernelId, 1.0, 0.0,
                                                         1 /*Location::Node*/,
                                                         boundingBox,
                                                         node2);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    errorCode = mkapi::mkernel_mesh2d_get_location_index(meshKernelId, 1.0, 1.0,
                                                         1 /*Location::Node*/,
                                                         boundingBox,
                                                         node3);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    errorCode = mkapi::mkernel_mesh2d_get_location_index(meshKernelId, 0.0, 1.0,
                                                         1 /*Location::Node*/,
                                                         boundingBox,
                                                         node4);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    int newNodeId;
    errorCode = mkapi::mkernel_mesh2d_insert_node(meshKernelId, 0.5, 0.5, newNodeId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    int edge1;
    int edge2;
    int edge3;
    int edge4;
    errorCode = mkapi::mkernel_mesh2d_insert_edge(meshKernelId, node1, newNodeId, edge1);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    errorCode = mkapi::mkernel_mesh2d_insert_edge(meshKernelId, node2, newNodeId, edge2);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    errorCode = mkapi::mkernel_mesh2d_insert_edge(meshKernelId, node3, newNodeId, edge3);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    errorCode = mkapi::mkernel_mesh2d_insert_edge(meshKernelId, node4, newNodeId, edge4);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    // All of the newly added edges should be removed.
    // But, there should still be space for them in the mesh node and edge lists.
    errorCode = mkapi::mkernel_mesh2d_delete_node(meshKernelId, newNodeId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    mkapi::Mesh2D mesh2d{};
    errorCode = mkapi::mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    std::vector<double> nodeX(mesh2d.num_nodes);
    std::vector<double> nodeY(mesh2d.num_nodes);
    std::vector<int> edges(2 * mesh2d.num_edges);

    mesh2d.node_x = nodeX.data();
    mesh2d.node_y = nodeY.data();
    mesh2d.edge_nodes = edges.data();

    errorCode = mkapi::mkernel_mesh2d_get_node_edge_data(meshKernelId, mesh2d);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    std::vector<double> expectedX{0.0, 1.0, 0.0, 1.0, meshkernel::constants::missing::doubleValue};
    std::vector<double> expectedY{0.0, 0.0, 1.0, 1.0, meshkernel::constants::missing::doubleValue};
    std::vector<int> expectedEdges{0, 2, 1, 3, 0, 1, 2, 3, -1, -1, -1, -1, -1, -1, -1, -1};

    ASSERT_EQ(mesh2d.num_nodes, clgSize * clgSize + 1);
    ASSERT_EQ(mesh2d.num_edges, 2 * (clgSize - 1) * clgSize + 4);

    for (size_t i = 0; i < nodeX.size(); ++i)
    {
        EXPECT_EQ(nodeX[i], expectedX[i]);
        EXPECT_EQ(nodeY[i], expectedY[i]);
    }

    for (size_t i = 0; i < edges.size(); ++i)
    {
        EXPECT_EQ(edges[i], expectedEdges[i]);
    }
}
