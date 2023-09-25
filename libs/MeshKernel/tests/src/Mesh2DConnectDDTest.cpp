#include <chrono>
#include <gtest/gtest.h>
#include <string>

#include <MeshKernel/ConnectMeshes.hpp>
#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <TestUtils/Definitions.hpp>
#include <TestUtils/MakeMeshes.hpp>

void CheckConnectGrids(const std::string& unconnectedGridName, const std::string& connectedGridName)
{
    static const std::string testDataDir = TEST_FOLDER + "/data/ConnectCurvilinearQuadsDDType/";

    meshkernel::ConnectMeshes connectCurviliearMeshes;

    // Grid to connect hanging node across the DD boundary.
    auto unconnectedGrid = ReadLegacyMesh2DFromFile(testDataDir + unconnectedGridName);

    // Expected grid after connecting hanging nodes.
    auto connectedGrid = ReadLegacyMesh2DFromFile(testDataDir + connectedGridName);

    // Connect hanging nodes
    connectCurviliearMeshes.Compute(*unconnectedGrid);

    // Check mesh entity counts are the same
    ASSERT_EQ(unconnectedGrid->GetNumNodes(), connectedGrid->GetNumNodes());
    ASSERT_EQ(unconnectedGrid->GetNumEdges(), connectedGrid->GetNumEdges());
    ASSERT_EQ(unconnectedGrid->GetNumFaces(), connectedGrid->GetNumFaces());

    constexpr double tolerance = 1.0e-10;

    // Check nodes
    for (meshkernel::UInt i = 0; i < unconnectedGrid->GetNumNodes(); ++i)
    {
        EXPECT_NEAR(unconnectedGrid->m_nodes[i].x, connectedGrid->m_nodes[i].x, tolerance);
        EXPECT_NEAR(unconnectedGrid->m_nodes[i].y, connectedGrid->m_nodes[i].y, tolerance);
    }

    // Check edges
    for (meshkernel::UInt i = 0; i < unconnectedGrid->GetNumEdges(); ++i)
    {
        EXPECT_EQ(unconnectedGrid->m_edges[i].first, connectedGrid->m_edges[i].first);
        EXPECT_EQ(unconnectedGrid->m_edges[i].second, connectedGrid->m_edges[i].second);
    }

    // Check edge-faces
    for (meshkernel::UInt i = 0; i < unconnectedGrid->GetNumEdges(); ++i)
    {
        EXPECT_EQ(unconnectedGrid->m_edgesFaces[i][0], connectedGrid->m_edgesFaces[i][0]);
        EXPECT_EQ(unconnectedGrid->m_edgesFaces[i][1], connectedGrid->m_edgesFaces[i][1]);
    }

    // Check number of nodes for each face
    for (meshkernel::UInt i = 0; i < unconnectedGrid->GetNumFaces(); ++i)
    {
        EXPECT_EQ(unconnectedGrid->m_numFacesNodes[i], connectedGrid->m_numFacesNodes[i]);
    }
}

// Tests are separated on number of hanging nodes (1-, 2-, 3- or 4-irregular mesh) or the complexity.

TEST(Mesh2DConnectDD, ConnectGridSimple1Level)
{
    // Test connecting edges with 1 hanging node per irregular edge
    CheckConnectGrids("unmerged_simple_net_east.nc", "merged_simple_net_east.nc");
    CheckConnectGrids("unmerged_simple_net_west.nc", "merged_simple_net_west.nc");
    CheckConnectGrids("unmerged_simple_net_north.nc", "merged_simple_net_north.nc");
}

TEST(Mesh2DConnectDD, ConnectGridSimple2Levels)
{
    // Test connecting edges with 2 hanging nodes per irregular edge
    CheckConnectGrids("unmatched_simple_2_levels_north.nc", "matched_simple_2_levels_north.nc");
}

TEST(Mesh2DConnectDD, ConnectGridSimple3Levels)
{
    // Test connecting edges with 3 hanging nodes per irregular edge
    CheckConnectGrids("unmatched_simple_3_levels_west.nc", "matched_simple_3_levels_west.nc");
}

TEST(Mesh2DConnectDD, ConnectGridSimple4Levels)
{
    // Test connecting edges with 4 hanging nodes per irregular edge
    CheckConnectGrids("unmatched_simple_4_levels_south.nc", "matched_simple_4_levels_south.nc");
}

TEST(Mesh2DConnectDD, ConnectGridComplexTest)
{
    // Test connecting edges with 1, 2, 3 and 4 hanging nodes per irregular edge
    // The main domain is a 10x10 element grid (11x11 nodes)
    // On each edge of this is a small rectangle each progressively more refined
    // than the previous.
    CheckConnectGrids("unmatched_all_sides.nc", "matched_all_sides.nc");
}

TEST(Mesh2DConnectDD, ConnectGridSimple4HangingNodes2Element)
{
    // Test connecting edges with 4 hanging nodes along irregular edge, with two master elements deep
    // When there are 4 hanging nodes, upto 3 elements on the coarse side will be affected
    CheckConnectGrids("unmatched_4_hanging_nodes_2_elements.nc", "matched_4_hanging_nodes_2_elements.nc");
}

TEST(Mesh2DConnectDD, ConnectGridSimple4HangingNodes1Element)
{
    // Test connecting edges with 4 hanging nodes along irregular edge, with single master element
    // When there are 4 hanging nodes, upto 3 elements on the coarse side will be affected
    CheckConnectGrids("unmatched_4_hanging_nodes_1_element.nc", "matched_4_hanging_nodes_1_element.nc");
}

TEST(Mesh2DConnectDD, ConnectGridSimple3HangingNodes)
{
    // Test connecting edges with 3 hanging nodes along irregular edge
    CheckConnectGrids("unmatched_3_hanging_nodes.nc", "matched_3_hanging_nodes.nc");
}

TEST(Mesh2DConnectDD, ConnectGridSimple2HangingNodes)
{
    // Test connecting edges with 2 hanging nodes along irregular edge
    CheckConnectGrids("unmatched_2_hanging_nodes.nc", "matched_2_hanging_nodes.nc");
}
