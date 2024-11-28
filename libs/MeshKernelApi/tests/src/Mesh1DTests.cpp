#include <algorithm>
#include <gtest/gtest.h>
#include <random>

#include "CartesianApiTestFixture.hpp"
#include "MeshKernel/Parameters.hpp"
#include "MeshKernelApi/BoundingBox.hpp"
#include "MeshKernelApi/Mesh1D.hpp"
#include "MeshKernelApi/MeshKernel.hpp"
#include "TestUtils/MakeMeshes.hpp"

// namespace aliases
namespace mk = meshkernel;
namespace mkapi = meshkernelapi;

TEST(Mesh1D, Mesh1DSetAndAdd)
{
    using namespace meshkernelapi;

    // Init 1d mesh
    meshkernelapi::Mesh1D mesh1d_1;
    meshkernelapi::Mesh1D mesh1d_2;

    std::vector<double> node_x{
        1.73493900000000,
        2.35659313023165,
        5.38347452702839,
        14.2980910429074,
        22.9324017677239,
        25.3723169493137,
        25.8072280000000};

    std::vector<double> node_y{
        -7.6626510000000,
        1.67281447902331,
        10.3513746546384,
        12.4797224193970,
        15.3007317677239,
        24.1623588554512,
        33.5111870000000};

    std::vector<int> edge_nodes{0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6};

    mesh1d_1.node_x = node_x.data();
    mesh1d_1.node_y = node_y.data();
    mesh1d_1.num_nodes = static_cast<int>(node_x.size());
    mesh1d_1.edge_nodes = edge_nodes.data();
    mesh1d_1.num_edges = static_cast<int>(edge_nodes.size()) / 2;

    // do not overwrite node_x
    std::vector<double> node_x_cp(node_x);
    double const offset = node_x_cp.back() + 1.0;
    std::transform(node_x_cp.begin(),
                   node_x_cp.end(),
                   node_x_cp.begin(),
                   [offset](double const val)
                   { return val + offset; });
    mesh1d_2.node_x = node_x_cp.data();
    mesh1d_2.node_y = node_y.data();
    mesh1d_2.num_nodes = static_cast<int>(node_x_cp.size());
    mesh1d_2.edge_nodes = edge_nodes.data();
    mesh1d_2.num_edges = static_cast<int>(edge_nodes.size()) / 2;

    // allocate state
    int mk_id = 0;
    int errorCode = mkernel_allocate_state(0, mk_id);

    // first initialise using the first mesh, mesh1d_1
    errorCode = mkernel_mesh1d_set(mk_id, mesh1d_1);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // then add the second mesh, mesh1d_2
    errorCode = mkernel_mesh1d_add(mk_id, mesh1d_2);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // get the dimensions and data of the resulting mesh
    Mesh1D mesh1d;
    errorCode = mkernel_mesh1d_get_dimensions(mk_id, mesh1d);
    std::vector<double> meshNodesX(mesh1d.num_nodes);
    std::vector<double> meshNodesY(mesh1d.num_nodes);
    std::vector<int> meshEdges(mesh1d.num_edges * 2);

    mesh1d.node_x = meshNodesX.data();
    mesh1d.node_y = meshNodesY.data();
    mesh1d.edge_nodes = meshEdges.data();
    errorCode = mkernel_mesh1d_get_data(mk_id, mesh1d);

    ASSERT_EQ(mesh1d.num_nodes, mesh1d_1.num_nodes + mesh1d_2.num_nodes);
    ASSERT_EQ(mesh1d.num_edges, mesh1d_1.num_edges + mesh1d_1.num_edges);

    for (int i = 0; i < mesh1d_1.num_nodes; ++i)
    {
        EXPECT_EQ(mesh1d.node_x[i], mesh1d_1.node_x[i]);
        EXPECT_EQ(mesh1d.node_y[i], mesh1d_1.node_y[i]);
    }

    for (int i = mesh1d_1.num_nodes; i < mesh1d_1.num_nodes + mesh1d_2.num_nodes; ++i)
    {
        EXPECT_EQ(mesh1d.node_x[i], mesh1d_2.node_x[i - mesh1d_1.num_nodes]);
        EXPECT_EQ(mesh1d.node_y[i], mesh1d_2.node_y[i - mesh1d_1.num_nodes]);
    }

    errorCode = mkernel_deallocate_state(mk_id);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}

TEST_F(CartesianApiTestFixture, GetDimensionsMesh1D_WithMesh1D_ShouldGetDimensionsMesh1D)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    std::vector nodesX{-16.1886410000000,
                       -16.1464995876014,
                       -16.1043581752028,
                       -16.0622167628042,
                       -15.7539488236928,
                       -6.86476658679268,
                       2.02441565010741,
                       10.9135970000000};

    std::vector nodesY{0.89018900000000,
                       9.78201442138723,
                       18.6738398427745,
                       27.5656652641617,
                       36.1966603330179,
                       36.4175095626911,
                       36.6383587923643,
                       36.8592080000000};

    std::vector edges{0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7};

    meshkernelapi::Mesh1D mesh1d;
    mesh1d.edge_nodes = edges.data();
    mesh1d.node_x = nodesX.data();
    mesh1d.node_y = nodesY.data();
    mesh1d.num_nodes = 8;
    mesh1d.num_edges = 7;

    auto errorCode = mkernel_mesh1d_set(GetMeshKernelId(), mesh1d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    meshkernelapi::Mesh1D mesh1dResults;
    errorCode = mkernel_mesh1d_get_dimensions(meshKernelId, mesh1dResults);

    // Assert
    ASSERT_EQ(8, mesh1dResults.num_nodes);
    ASSERT_EQ(7, mesh1dResults.num_edges);
}

TEST_F(CartesianApiTestFixture, GetDataMesh1D_WithMesh1D_ShouldGetDataMesh1D)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    std::vector nodesX{-16.1886410000000,
                       -16.1464995876014,
                       -16.1043581752028,
                       -16.0622167628042,
                       -15.7539488236928,
                       -6.86476658679268,
                       2.02441565010741,
                       10.9135970000000};

    std::vector nodesY{0.89018900000000,
                       9.78201442138723,
                       18.6738398427745,
                       27.5656652641617,
                       36.1966603330179,
                       36.4175095626911,
                       36.6383587923643,
                       36.8592080000000};

    std::vector edges{0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7};

    meshkernelapi::Mesh1D mesh1d;
    mesh1d.edge_nodes = edges.data();
    mesh1d.node_x = nodesX.data();
    mesh1d.node_y = nodesY.data();
    mesh1d.num_nodes = 8;
    mesh1d.num_edges = 7;

    auto errorCode = mkernel_mesh1d_set(GetMeshKernelId(), mesh1d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    meshkernelapi::Mesh1D mesh1dResults;
    errorCode = mkernel_mesh1d_get_dimensions(meshKernelId, mesh1dResults);
    std::vector<double> meshNodesX(mesh1dResults.num_nodes);
    std::vector<double> meshNodesY(mesh1dResults.num_nodes);
    std::vector<int> meshEdges(mesh1dResults.num_edges * 2);

    mesh1dResults.node_x = meshNodesX.data();
    mesh1dResults.node_y = meshNodesY.data();
    mesh1dResults.edge_nodes = meshEdges.data();
    errorCode = mkernel_mesh1d_get_data(meshKernelId, mesh1dResults);

    // Assert
    std::vector validMeshNodesX(nodesX.data(), nodesX.data() + mesh1d.num_nodes);
    std::vector computedMeshNodesX(meshNodesX.data(), meshNodesX.data() + mesh1dResults.num_nodes);
    ASSERT_THAT(computedMeshNodesX, ::testing::ContainerEq(validMeshNodesX));

    std::vector validMeshNodesY(nodesY.data(), nodesY.data() + mesh1d.num_nodes);
    std::vector computedMeshNodesY(meshNodesY.data(), meshNodesY.data() + mesh1dResults.num_nodes);
    ASSERT_THAT(computedMeshNodesY, ::testing::ContainerEq(validMeshNodesY));

    std::vector<double> validEdges(edges.data(), edges.data() + mesh1d.num_edges);
    std::vector<double> computedEdges(meshEdges.data(), meshEdges.data() + mesh1dResults.num_edges);
    ASSERT_THAT(computedEdges, ::testing::ContainerEq(validEdges));
}
