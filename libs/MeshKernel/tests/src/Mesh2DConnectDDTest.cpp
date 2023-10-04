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

std::shared_ptr<meshkernel::Mesh2D> generateMesh(const meshkernel::Point& origin, const meshkernel::Vector& delta, const int n, const int m)
{
    double dimX = static_cast<double>(n - 1) * delta.x();
    double dimY = static_cast<double>(n - 1) * delta.y();
    return MakeRectangularMeshForTesting(n, m, dimX, dimY, meshkernel::Projection::cartesian, origin);
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

    std::shared_ptr<meshkernel::Mesh2D> mesh1 = generateMesh(origin, delta, 11, 11);

    origin.x += 10.0 * delta.x();
    origin.y += delta.y();

    delta = meshkernel::Vector{2.5, 2.5};
    std::shared_ptr<meshkernel::Mesh2D> mesh2 = generateMesh(origin, delta, 9, 9);

    meshkernel::Mesh2D mergedMesh = meshkernel::Mesh2D::Merge(*mesh1, *mesh2);
    meshkernel::ConnectMeshes connectCurviliearMeshes;
    connectCurviliearMeshes.Compute(mergedMesh);

    EXPECT_EQ(mergedMesh.GetNumFaces(), mesh1->GetNumFaces() + mesh2->GetNumFaces() + ExtraFaces);
    EXPECT_EQ(mergedMesh.GetNumNodes(), mesh1->GetNumNodes() + mesh2->GetNumNodes() - NodeDifference);
}

TEST(Mesh2DConnectDD, MergeOneEmptyMesh)
{

    meshkernel::Point origin{0.0, 0.0};
    meshkernel::Vector delta{10.0, 10.0};

    std::shared_ptr<meshkernel::Mesh2D> mesh1 = generateMesh(origin, delta, 11, 11);
    meshkernel::Mesh2D mesh2;

    // The projection needs to be set, as there is no default
    mesh2.m_projection = meshkernel::Projection::cartesian;

    meshkernel::Mesh2D mergedMesh = meshkernel::Mesh2D::Merge(*mesh1, mesh2);

    EXPECT_EQ(mergedMesh.GetNumFaces(), mesh1->GetNumFaces());
    EXPECT_EQ(mergedMesh.GetNumNodes(), mesh1->GetNumNodes());

    // This time with the parameters reversed
    meshkernel::Mesh2D anotherMergedMesh = meshkernel::Mesh2D::Merge(mesh2, *mesh1);

    EXPECT_EQ(anotherMergedMesh.GetNumFaces(), mesh1->GetNumFaces());
    EXPECT_EQ(anotherMergedMesh.GetNumNodes(), mesh1->GetNumNodes());
}

TEST(Mesh2DConnectDD, MergeTwoEmptyMeshes)
{

    meshkernel::Mesh2D mesh1;
    meshkernel::Mesh2D mesh2;

    mesh1.m_projection = meshkernel::Projection::cartesian;
    mesh2.m_projection = meshkernel::Projection::cartesian;

    meshkernel::Mesh2D mergedMesh;

    EXPECT_THROW(mergedMesh = meshkernel::Mesh2D::Merge(mesh1, mesh2), meshkernel::MeshKernelError);
}

TEST(Mesh2DConnectDD, MergeTwoIncompatibleMeshes)
{

    meshkernel::Point origin{0.0, 0.0};
    meshkernel::Vector delta{10.0, 10.0};

    std::shared_ptr<meshkernel::Mesh2D> mesh1 = generateMesh(origin, delta, 11, 11);
    std::shared_ptr<meshkernel::Mesh2D> mesh2 = generateMesh(origin, delta, 11, 11);

    // change projection on mesh2.
    mesh2->m_projection = meshkernel::Projection::spherical;

    EXPECT_THROW([[maybe_unused]] meshkernel::Mesh2D mergedMesh = meshkernel::Mesh2D::Merge(*mesh1, *mesh2), meshkernel::MeshKernelError);
}
