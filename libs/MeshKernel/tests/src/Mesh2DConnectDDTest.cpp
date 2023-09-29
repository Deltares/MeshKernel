#include <chrono>
#include <gtest/gtest.h>
#include <string>

#include <MeshKernel/ConnectMeshes.hpp>
#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Operations.hpp>
#include <TestUtils/Definitions.hpp>
#include <TestUtils/MakeMeshes.hpp>

void CheckConnectGrids(const std::string& unconnectedGridName, const std::string& connectedGridName)
{
    static const std::string testDataDir = TEST_FOLDER + "/data/ConnectCurvilinearQuadsDDType/";

    // Grid to connect hanging node across the DD boundary.
    auto unconnectedGrid = ReadLegacyMesh2DFromFile(testDataDir + unconnectedGridName);

    // Expected grid after connecting hanging nodes.
    auto connectedGrid = ReadLegacyMesh2DFromFile(testDataDir + connectedGridName);

    // Connect hanging nodes
    meshkernel::ConnectMeshes connectCurviliearMeshes;
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

meshkernel::Mesh2D generateMesh(const meshkernel::Point& origin, const meshkernel::Vector& delta, const int n, const int m)
{

    std::vector<std::vector<int>> indicesValues(n, std::vector<int>(m));
    std::vector<meshkernel::Point> nodes(n * m);
    std::size_t nodeIndex = 0;

    for (int j = 0; j < m; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            indicesValues[i][j] = i + j * n;
            nodes[nodeIndex] = {origin.x + i * delta.x(), origin.y + j * delta.y()};
            nodeIndex++;
        }
    }

    std::vector<meshkernel::Edge> edges((n - 1) * m + (m - 1) * n);
    std::size_t edgeIndex = 0;
    for (auto j = 0; j < m; ++j)
    {
        for (auto i = 0; i < n - 1; ++i)
        {
            edges[edgeIndex] = {indicesValues[i][j], indicesValues[i + 1][j]};
            edgeIndex++;
        }
    }

    for (auto j = 0; j < m - 1; ++j)
    {
        for (auto i = 0; i < n; ++i)
        {
            edges[edgeIndex] = {indicesValues[i][j + 1], indicesValues[i][j]};
            edgeIndex++;
        }
    }

    return meshkernel::Mesh2D(edges, nodes, meshkernel::Projection::cartesian);
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

TEST(Mesh2DConnectDD, MergeMeshes)
{

    // 3 nodes were removed because they are shared between both meshes
    // 2 node were added in order to free the hanging nodes
    const meshkernel::UInt NodeDifference = 3 - 2;

    // 18 faces were added when freeing the hanging nodes
    // 4 faces were removed and replaced by the new 18 faces.
    const meshkernel::UInt ExtraFaces = 18 - 4;

    meshkernel::Point origin{0.0, 0.0};
    meshkernel::Vector delta{10.0, 10.0};

    meshkernel::Mesh2D mesh1 = generateMesh(origin, delta, 11, 11);

    origin.x += 10.0 * delta.x();
    origin.y += delta.y();

    delta = meshkernel::Vector{2.5, 2.5};
    meshkernel::Mesh2D mesh2 = generateMesh(origin, delta, 9, 9);

    meshkernel::Mesh2D mergedMesh = meshkernel::Mesh2D::Merge(mesh1, mesh2);
    meshkernel::ConnectMeshes connectCurviliearMeshes;
    connectCurviliearMeshes.Compute(mergedMesh);

    EXPECT_EQ(mergedMesh.GetNumFaces(), mesh1.GetNumFaces() + mesh2.GetNumFaces() + ExtraFaces);
    EXPECT_EQ(mergedMesh.GetNumNodes(), mesh1.GetNumNodes() + mesh2.GetNumNodes() - NodeDifference);
}
