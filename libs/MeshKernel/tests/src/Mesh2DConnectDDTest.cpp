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
    meshkernel::ConnectMeshes::Compute(*unconnectedGrid);

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
    meshkernel::ConnectMeshes::Compute(mergedMesh);

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

TEST(Mesh2DConnectDD, MergeTwoSameMeshesSmallNegativeOffset)
{
    // Merge two meshes that have the same resolution

    const int nx = 3;
    const int ny = 3;

    meshkernel::Point origin{0.0, 0.0};
    meshkernel::Vector delta{10.0, 10.0};

    std::shared_ptr<meshkernel::Mesh2D> mesh1 = generateMesh(origin, delta, nx, ny);

    origin.x += delta.x() * static_cast<double>(nx - 1) - 1.0;

    std::shared_ptr<meshkernel::Mesh2D> mesh2 = generateMesh(origin, delta, nx, ny);

    meshkernel::Mesh2D mergedMesh = meshkernel::Mesh2D::Merge(*mesh1, *mesh2);

    meshkernel::ConnectMeshes::Compute(mergedMesh);

    EXPECT_EQ(mergedMesh.GetNumNodes(), 15);
    EXPECT_EQ(mergedMesh.GetNumFaces(), 8);
    EXPECT_EQ(mergedMesh.GetNumEdges(), 22);

    const double tolerance = 1.0e-8;

    // Only comparing the nodes that were along the overlapping edge
    EXPECT_NEAR(mergedMesh.m_nodes[6].x, 19.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[7].x, 19.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[8].x, 19.0, tolerance);

    EXPECT_NEAR(mergedMesh.m_nodes[6].y, 0.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[7].y, 10.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[8].y, 20.0, tolerance);

    EXPECT_EQ(mergedMesh.m_nodesNumEdges[6], 3);
    EXPECT_EQ(mergedMesh.m_nodesNumEdges[7], 4);
    EXPECT_EQ(mergedMesh.m_nodesNumEdges[8], 3);
}

TEST(Mesh2DConnectDD, MergeTwoSameMeshesNoOffset)
{
    // Merge two meshes that have the same resolution

    const int nx = 3;
    const int ny = 3;

    meshkernel::Point origin{0.0, 0.0};
    meshkernel::Vector delta{10.0, 10.0};

    std::shared_ptr<meshkernel::Mesh2D> mesh1 = generateMesh(origin, delta, nx, ny);

    origin.x += delta.x() * static_cast<double>(nx - 1);

    std::shared_ptr<meshkernel::Mesh2D> mesh2 = generateMesh(origin, delta, nx, ny);

    meshkernel::Mesh2D mergedMesh = meshkernel::Mesh2D::Merge(*mesh1, *mesh2);

    meshkernel::ConnectMeshes::Compute(mergedMesh);

    EXPECT_EQ(mergedMesh.GetNumNodes(), 15);
    EXPECT_EQ(mergedMesh.GetNumFaces(), 8);
    EXPECT_EQ(mergedMesh.GetNumEdges(), 22);

    const double tolerance = 1.0e-8;

    // Only comparing the nodes that were along the overlapping edge
    EXPECT_NEAR(mergedMesh.m_nodes[6].x, 20.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[7].x, 20.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[8].x, 20.0, tolerance);

    EXPECT_NEAR(mergedMesh.m_nodes[6].y, 0.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[7].y, 10.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[8].y, 20.0, tolerance);

    EXPECT_EQ(mergedMesh.m_nodesNumEdges[6], 3);
    EXPECT_EQ(mergedMesh.m_nodesNumEdges[7], 4);
    EXPECT_EQ(mergedMesh.m_nodesNumEdges[8], 3);
}

TEST(Mesh2DConnectDD, MergeTwoSameMeshesSmallPositiveOffset)
{
    // Merge two meshes that have the same resolution

    const int nx = 3;
    const int ny = 3;

    meshkernel::Point origin{0.0, 0.0};
    meshkernel::Vector delta{10.0, 10.0};

    std::shared_ptr<meshkernel::Mesh2D> mesh1 = generateMesh(origin, delta, nx, ny);

    origin.x += delta.x() * static_cast<double>(nx - 1) + 1.0;

    std::shared_ptr<meshkernel::Mesh2D> mesh2 = generateMesh(origin, delta, nx, ny);

    meshkernel::Mesh2D mergedMesh = meshkernel::Mesh2D::Merge(*mesh1, *mesh2);

    meshkernel::ConnectMeshes::Compute(mergedMesh);

    EXPECT_EQ(mergedMesh.GetNumNodes(), 15);
    EXPECT_EQ(mergedMesh.GetNumFaces(), 8);
    EXPECT_EQ(mergedMesh.GetNumEdges(), 22);

    const double tolerance = 1.0e-8;

    // Only comparing the nodes that were along the overlapping edge
    EXPECT_NEAR(mergedMesh.m_nodes[6].x, 21.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[7].x, 21.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[8].x, 21.0, tolerance);

    EXPECT_NEAR(mergedMesh.m_nodes[6].y, 0.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[7].y, 10.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[8].y, 20.0, tolerance);

    EXPECT_EQ(mergedMesh.m_nodesNumEdges[6], 3);
    EXPECT_EQ(mergedMesh.m_nodesNumEdges[7], 4);
    EXPECT_EQ(mergedMesh.m_nodesNumEdges[8], 3);
}

TEST(Mesh2DConnectDD, MergeTwoMeshesWithSmallNegativeOffset)
{
    const int nx = 4;
    const int ny = 4;

    meshkernel::Point origin{0.0, 0.0};
    meshkernel::Vector delta{10.0, 10.0};

    std::shared_ptr<meshkernel::Mesh2D> mesh1 = generateMesh(origin, delta, nx, ny);

    origin.x += delta.x() * static_cast<double>(nx - 1) - 1.0;
    delta.y() *= 0.31;

    std::shared_ptr<meshkernel::Mesh2D> mesh2 = generateMesh(origin, delta, nx, ny);

    meshkernel::Mesh2D mergedMesh = meshkernel::Mesh2D::Merge(*mesh1, *mesh2);

    // Need to increase the search distance fraction
    meshkernel::ConnectMeshes::Compute(mergedMesh);

    // 8 triangles and 16 quadrilaterals
    EXPECT_EQ(mergedMesh.GetNumFaces(), 24);
    EXPECT_EQ(mergedMesh.GetNumNodes(), 31);
    EXPECT_EQ(mergedMesh.GetNumEdges(), 54);

    const double tolerance = 1.0e-8;

    // Only comparing the nodes that were along the irregular edge and the
    // edge where the node was created in order to free the hanging nodes
    EXPECT_NEAR(mergedMesh.m_nodes[14].x, 29.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[15].x, 29.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[16].x, 29.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[17].x, 29.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[8].x, 20.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[9].x, 20.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[30].x, 20.0, tolerance);

    EXPECT_NEAR(mergedMesh.m_nodes[14].y, 0.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[15].y, 3.1, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[16].y, 6.2, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[17].y, 9.3, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[8].y, 0.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[9].y, 10.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[30].y, 5.0, tolerance);

    const meshkernel::UInt nullValue = meshkernel::constants::missing::uintValue;

    // Allocate enough space for all edge, but will only check the edges around what was the irregular edge
    std::vector<meshkernel::Edge> expectedEdges(54, {nullValue, nullValue});

    // Only checking the edges connected to the (initial) irregular edge and the
    // edge where the node was created in order to free the hanging nodes
    expectedEdges[4].first = 4;
    expectedEdges[4].second = 8;
    expectedEdges[5].first = 5;
    expectedEdges[5].second = 9;
    expectedEdges[8].first = 8;
    expectedEdges[8].second = 14;
    expectedEdges[9].first = 9;
    expectedEdges[9].second = 17;
    expectedEdges[18].first = 10;
    expectedEdges[18].second = 9;
    expectedEdges[23].first = 15;
    expectedEdges[23].second = 19;
    expectedEdges[24].first = 16;
    expectedEdges[24].second = 20;
    expectedEdges[34].first = 15;
    expectedEdges[34].second = 14;
    expectedEdges[35].first = 16;
    expectedEdges[35].second = 15;
    expectedEdges[36].first = 17;
    expectedEdges[36].second = 16;
    expectedEdges[46].first = 16;
    expectedEdges[46].second = 30;
    expectedEdges[47].first = 16;
    expectedEdges[47].second = 9;
    expectedEdges[48].first = 15;
    expectedEdges[48].second = 30;
    expectedEdges[49].first = 15;
    expectedEdges[49].second = 8;
    expectedEdges[50].first = 30;
    expectedEdges[50].second = 9;
    expectedEdges[51].first = 30;
    expectedEdges[51].second = 8;
    expectedEdges[52].first = 30;
    expectedEdges[52].second = 5;
    expectedEdges[53].first = 30;
    expectedEdges[53].second = 4;

    for (meshkernel::UInt i = 0; i < mergedMesh.GetNumEdges(); ++i)
    {

        if (expectedEdges[i].first != nullValue)
        {
            EXPECT_EQ(expectedEdges[i].first, mergedMesh.m_edges[i].first);
            EXPECT_EQ(expectedEdges[i].second, mergedMesh.m_edges[i].second);
        }
    }
}

TEST(Mesh2DConnectDD, MergeTwoMeshesWithSmallPositiveOffset)
{
    const int nx = 4;
    const int ny = 4;

    meshkernel::Point origin{0.0, 0.0};
    meshkernel::Vector delta{10.0, 10.0};

    std::shared_ptr<meshkernel::Mesh2D> mesh1 = generateMesh(origin, delta, nx, ny);

    origin.x += delta.x() * static_cast<double>(nx - 1) + 1.0;
    delta.y() *= 0.31;

    std::shared_ptr<meshkernel::Mesh2D> mesh2 = generateMesh(origin, delta, nx, ny);

    meshkernel::Mesh2D mergedMesh = meshkernel::Mesh2D::Merge(*mesh1, *mesh2);

    // Need to increase the search distance fraction
    meshkernel::ConnectMeshes::Compute(mergedMesh);

    // 8 triangles and 16 quadrilaterals
    EXPECT_EQ(mergedMesh.GetNumFaces(), 24);
    EXPECT_EQ(mergedMesh.GetNumNodes(), 31);
    EXPECT_EQ(mergedMesh.GetNumEdges(), 54);

    const double tolerance = 1.0e-8;

    // Only comparing the nodes that were along the irregular edge and the
    // edge where the node was created in order to free the hanging nodes
    EXPECT_NEAR(mergedMesh.m_nodes[14].x, 31.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[15].x, 31.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[16].x, 31.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[17].x, 31.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[8].x, 20.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[9].x, 20.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[30].x, 20.0, tolerance);

    EXPECT_NEAR(mergedMesh.m_nodes[14].y, 0.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[15].y, 3.1, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[16].y, 6.2, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[17].y, 9.3, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[8].y, 0.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[9].y, 10.0, tolerance);
    EXPECT_NEAR(mergedMesh.m_nodes[30].y, 5.0, tolerance);

    const meshkernel::UInt nullValue = meshkernel::constants::missing::uintValue;

    // Allocate enough space for all edge, but will only check the edges around what was the irregular edge
    std::vector<meshkernel::Edge> expectedEdges(54, {nullValue, nullValue});

    // Only checking the edges connected to the (initial) irregular edge and the
    // edge where the node was created in order to free the hanging nodes
    expectedEdges[4].first = 4;
    expectedEdges[4].second = 8;
    expectedEdges[5].first = 5;
    expectedEdges[5].second = 9;
    expectedEdges[8].first = 8;
    expectedEdges[8].second = 14;
    expectedEdges[9].first = 9;
    expectedEdges[9].second = 17;
    expectedEdges[18].first = 10;
    expectedEdges[18].second = 9;
    expectedEdges[23].first = 15;
    expectedEdges[23].second = 19;
    expectedEdges[24].first = 16;
    expectedEdges[24].second = 20;
    expectedEdges[34].first = 15;
    expectedEdges[34].second = 14;
    expectedEdges[35].first = 16;
    expectedEdges[35].second = 15;
    expectedEdges[36].first = 17;
    expectedEdges[36].second = 16;
    expectedEdges[46].first = 16;
    expectedEdges[46].second = 30;
    expectedEdges[47].first = 16;
    expectedEdges[47].second = 9;
    expectedEdges[48].first = 15;
    expectedEdges[48].second = 30;
    expectedEdges[49].first = 15;
    expectedEdges[49].second = 8;
    expectedEdges[50].first = 30;
    expectedEdges[50].second = 9;
    expectedEdges[51].first = 30;
    expectedEdges[51].second = 8;
    expectedEdges[52].first = 30;
    expectedEdges[52].second = 5;
    expectedEdges[53].first = 30;
    expectedEdges[53].second = 4;

    for (meshkernel::UInt i = 0; i < mergedMesh.GetNumEdges(); ++i)
    {

        if (expectedEdges[i].first != nullValue)
        {
            EXPECT_EQ(expectedEdges[i].first, mergedMesh.m_edges[i].first);
            EXPECT_EQ(expectedEdges[i].second, mergedMesh.m_edges[i].second);
        }
    }
}

TEST(Mesh2DConnectDD, MergeTwoMeshesErrorInSeparationFraction)
{
    // The test checks that the meshkernel::ConnectMeshes::Compute function
    // fails when passed a separation fraction that is bout of a valid range

    const int nx = 3;
    const int ny = 3;

    meshkernel::Point origin{0.0, 0.0};
    meshkernel::Vector delta{10.0, 10.0};

    std::shared_ptr<meshkernel::Mesh2D> mesh1 = generateMesh(origin, delta, nx, ny);

    origin.x += delta.x() * static_cast<double>(nx - 1);

    std::shared_ptr<meshkernel::Mesh2D> mesh2 = generateMesh(origin, delta, nx, ny);

    meshkernel::Mesh2D mergedMesh = meshkernel::Mesh2D::Merge(*mesh1, *mesh2);

    EXPECT_THROW(meshkernel::ConnectMeshes::Compute(mergedMesh, 0.5), meshkernel::RangeError);
    EXPECT_THROW(meshkernel::ConnectMeshes::Compute(mergedMesh, 1.5), meshkernel::RangeError);
    EXPECT_THROW(meshkernel::ConnectMeshes::Compute(mergedMesh, -0.5), meshkernel::RangeError);
    EXPECT_THROW(meshkernel::ConnectMeshes::Compute(mergedMesh, 0.0), meshkernel::RangeError);
}
