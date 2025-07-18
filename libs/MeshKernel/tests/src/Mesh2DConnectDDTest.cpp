//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2025.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation version 3.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// contact: delft3d.support@deltares.nl
// Stichting Deltares
// P.O. Box 177
// 2600 MH Delft, The Netherlands
//
// All indications and logos of, and references to, "Delft3D" and "Deltares"
// are registered trademarks of Stichting Deltares, and remain the property of
// Stichting Deltares. All rights reserved.
//
//------------------------------------------------------------------------------

#include <chrono>
#include <gtest/gtest.h>
#include <string>

#include <MeshKernel/ConnectMeshes.hpp>
#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Utilities/Utilities.hpp>
#include <TestUtils/Definitions.hpp>
#include <TestUtils/MakeMeshes.hpp>

void CheckGridsConnectedCorrectly(const meshkernel::Mesh2D& connectedGrid,
                                  const meshkernel::Mesh2D& unconnectedGrid)
{

    constexpr double tolerance = 1.0e-10;

    meshkernel::UInt count = 0;

    // Check nodes
    for (meshkernel::UInt i = 0; i < unconnectedGrid.GetNumNodes(); ++i)
    {
        if (unconnectedGrid.Node(i).IsValid())
        {
            EXPECT_NEAR(unconnectedGrid.Node(i).x, connectedGrid.Node(count).x, tolerance);
            EXPECT_NEAR(unconnectedGrid.Node(i).y, connectedGrid.Node(count).y, tolerance);
            ++count;
        }
    }

    count = 0;

    for (meshkernel::UInt i = 0; i < unconnectedGrid.GetNumEdges(); ++i)
    {
        if (unconnectedGrid.GetEdge(i).first != meshkernel::constants::missing::uintValue)
        {
            EXPECT_TRUE(meshkernel::IsEqual(connectedGrid.Node(connectedGrid.GetEdge(count).first),
                                            unconnectedGrid.Node(unconnectedGrid.GetEdge(i).first),
                                            tolerance))
                << "edge.first indexes incorrect node";

            EXPECT_TRUE(meshkernel::IsEqual(connectedGrid.Node(connectedGrid.GetEdge(count).second),
                                            unconnectedGrid.Node(unconnectedGrid.GetEdge(i).second),
                                            tolerance))
                << "edge.second indexes incorrect node";

            ++count;
        }
    }

    std::vector<meshkernel::UInt> edgeMap(unconnectedGrid.GetNumEdges(), meshkernel::constants::missing::uintValue);
    std::vector<meshkernel::UInt> edgeMapInv(unconnectedGrid.GetNumValidEdges());

    count = 0;

    for (meshkernel::UInt i = 0; i < unconnectedGrid.GetNumEdges(); ++i)
    {
        if (unconnectedGrid.GetEdge(i).first != meshkernel::constants::missing::uintValue && unconnectedGrid.GetEdge(i).second != meshkernel::constants::missing::uintValue)
        {
            edgeMap[i] = count;
            edgeMapInv[count] = i;
            ++count;
        }
    }

    count = 0;

    // Check edge-faces
    for (meshkernel::UInt i = 0; i < unconnectedGrid.GetNumEdges(); ++i)
    {
        if (edgeMap[i] != meshkernel::constants::missing::uintValue)
        {
            EXPECT_EQ(unconnectedGrid.m_edgesFaces[i][0], connectedGrid.m_edgesFaces[count][0]);
            EXPECT_EQ(unconnectedGrid.m_edgesFaces[i][1], connectedGrid.m_edgesFaces[count][1]);
            ++count;
        }
    }

    // Check number of nodes for each face
    for (meshkernel::UInt i = 0; i < unconnectedGrid.GetNumFaces(); ++i)
    {
        EXPECT_EQ(unconnectedGrid.m_numFacesNodes[i], connectedGrid.m_numFacesNodes[i]);
    }
}

void CheckGridsDisconnectedCorrectly(const std::vector<meshkernel::Point>& originalNodes,
                                     const std::vector<meshkernel::Edge>& originalEdges,
                                     const meshkernel::Mesh2D& unconnectedGrid)
{
    meshkernel::UInt count = 0;

    ASSERT_EQ(originalNodes.size(), unconnectedGrid.GetNumValidNodes());
    ASSERT_EQ(originalEdges.size(), unconnectedGrid.GetNumValidEdges());

    count = 0;

    for (meshkernel::UInt i = 0; i < unconnectedGrid.Nodes().size(); ++i)
    {
        if (unconnectedGrid.Node(i).IsValid())
        {
            // Check valid nodes
            EXPECT_EQ(originalNodes[count].x, unconnectedGrid.Node(i).x);
            EXPECT_EQ(originalNodes[count].y, unconnectedGrid.Node(i).y);
            ++count;
        }
    }

    count = 0;

    for (meshkernel::UInt i = 0; i < unconnectedGrid.Edges().size(); ++i)
    {
        if (unconnectedGrid.IsValidEdge(i))
        {
            // Check valid edges
            EXPECT_EQ(originalEdges[count].first, unconnectedGrid.GetEdge(i).first);
            EXPECT_EQ(originalEdges[count].second, unconnectedGrid.GetEdge(i).second);
            ++count;
        }
    }
}

void CheckConnectGrids(const std::string& unconnectedGridName, const std::string& connectedGridName)
{
    static const std::string testDataDir = TEST_FOLDER + "/data/ConnectCurvilinearQuadsDDType/";

    // Grid to connect hanging node across the DD boundary.
    auto unconnectedGrid = ReadLegacyMesh2DFromFile(testDataDir + unconnectedGridName);

    const std::vector<meshkernel::Point> originalNodes(unconnectedGrid->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(unconnectedGrid->Edges());

    // Expected grid after connecting hanging nodes.
    auto connectedGrid = ReadLegacyMesh2DFromFile(testDataDir + connectedGridName);

    // Connect hanging nodes
    auto undoAction = meshkernel::ConnectMeshes::Compute(*unconnectedGrid);

    // Check mesh (active) entity counts are the same
    ASSERT_EQ(unconnectedGrid->GetNumValidNodes(), connectedGrid->GetNumNodes());
    ASSERT_EQ(unconnectedGrid->GetNumValidEdges(), connectedGrid->GetNumEdges());
    ASSERT_EQ(unconnectedGrid->GetNumFaces(), connectedGrid->GetNumFaces());

    CheckGridsConnectedCorrectly(*connectedGrid, *unconnectedGrid);

    ////////////////////////////////
    // Test the undo action has been performed correctly
    undoAction->Restore();
    // Recompute faces
    unconnectedGrid->Administrate();

    CheckGridsDisconnectedCorrectly(originalNodes, originalEdges, *unconnectedGrid);

    ////////////////////////////////
    // Test the commit (redo) action has been performed correctly
    // The merged grid should be in the state it was directly after the connect call
    undoAction->Commit();
    // Recompute faces
    unconnectedGrid->Administrate();

    CheckGridsConnectedCorrectly(*connectedGrid, *unconnectedGrid);
}

std::shared_ptr<meshkernel::Mesh2D> generateMesh(const meshkernel::Point& origin, const meshkernel::Vector& delta, const int n, const int m)
{
    double dimX = static_cast<double>(n - 1) * delta.x();
    double dimY = static_cast<double>(m - 1) * delta.y();
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

    const auto mergedMesh = meshkernel::Mesh2D::Merge(mesh1->Nodes(), mesh1->Edges(),
                                                      mesh2->Nodes(), mesh2->Edges(),
                                                      mesh1->m_projection);

    const std::vector<meshkernel::Point> originalNodes(mergedMesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mergedMesh->Edges());

    auto undoAction = meshkernel::ConnectMeshes::Compute(*mergedMesh);

    EXPECT_EQ(mergedMesh->GetNumFaces(), mesh1->GetNumFaces() + mesh2->GetNumFaces() + ExtraFaces);
    EXPECT_EQ(mergedMesh->GetNumValidNodes(), mesh1->GetNumNodes() + mesh2->GetNumNodes() - NodeDifference);

    // Test the undo action has been computed correctly
    undoAction->Restore();
    // Recompute faces
    mergedMesh->Administrate();

    ASSERT_EQ(originalNodes.size(), mergedMesh->GetNumValidNodes());
    ASSERT_EQ(originalEdges.size(), mergedMesh->GetNumValidEdges());

    meshkernel::UInt count = 0;

    for (meshkernel::UInt i = 0; i < mergedMesh->Nodes().size(); ++i)
    {
        if (mergedMesh->Node(i).IsValid())
        {
            // Check valid nodes
            EXPECT_EQ(originalNodes[count].x, mergedMesh->Node(i).x);
            EXPECT_EQ(originalNodes[count].y, mergedMesh->Node(i).y);
            ++count;
        }
    }

    count = 0;

    for (meshkernel::UInt i = 0; i < mergedMesh->Edges().size(); ++i)
    {
        if (mergedMesh->IsValidEdge(i))
        {
            EXPECT_EQ(originalEdges[count].first, mergedMesh->GetEdge(i).first);
            EXPECT_EQ(originalEdges[count].second, mergedMesh->GetEdge(i).second);
            ++count;
        }
    }
}

TEST(Mesh2DConnectDD, MergeOneEmptyMesh)
{

    meshkernel::Point origin{0.0, 0.0};
    meshkernel::Vector delta{10.0, 10.0};

    std::shared_ptr<meshkernel::Mesh2D> mesh1 = generateMesh(origin, delta, 11, 11);
    meshkernel::Mesh2D mesh2;

    // The projection needs to be set, as there is no default
    mesh2.m_projection = meshkernel::Projection::cartesian;

    const auto mergedMesh = meshkernel::Mesh2D::Merge(mesh1->Nodes(), mesh1->Edges(),
                                                      mesh2.Nodes(), mesh2.Edges(),
                                                      mesh1->m_projection);

    EXPECT_EQ(mergedMesh->GetNumFaces(), mesh1->GetNumFaces());
    EXPECT_EQ(mergedMesh->GetNumNodes(), mesh1->GetNumNodes());

    // This time with the parameters reversed
    const auto anotherMergedMesh = meshkernel::Mesh2D::Merge(mesh2.Nodes(), mesh2.Edges(),
                                                             mesh1->Nodes(), mesh1->Edges(),
                                                             mesh1->m_projection);

    EXPECT_EQ(anotherMergedMesh->GetNumFaces(), mesh1->GetNumFaces());
    EXPECT_EQ(anotherMergedMesh->GetNumValidNodes(), mesh1->GetNumNodes());
}

TEST(Mesh2DConnectDD, MergeTwoEmptyMeshes)
{

    meshkernel::Mesh2D mesh1;
    meshkernel::Mesh2D mesh2;

    mesh1.m_projection = meshkernel::Projection::cartesian;
    mesh2.m_projection = meshkernel::Projection::cartesian;

    EXPECT_THROW(const auto mergedMesh = meshkernel::Mesh2D::Merge(mesh1, mesh2), meshkernel::MeshKernelError);
}

TEST(Mesh2DConnectDD, MergeTwoIncompatibleMeshes)
{

    meshkernel::Point origin{0.0, 0.0};
    meshkernel::Vector delta{10.0, 10.0};

    std::shared_ptr<meshkernel::Mesh2D> mesh1 = generateMesh(origin, delta, 11, 11);
    std::shared_ptr<meshkernel::Mesh2D> mesh2 = generateMesh(origin, delta, 11, 11);

    // change projection on mesh2.
    mesh2->m_projection = meshkernel::Projection::spherical;

    EXPECT_THROW([[maybe_unused]] const auto mergedMesh = meshkernel::Mesh2D::Merge(*mesh1, *mesh2), meshkernel::MeshKernelError);
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

    const auto mergedMesh = meshkernel::Mesh2D::Merge(mesh1->Nodes(), mesh1->Edges(),
                                                      mesh2->Nodes(), mesh2->Edges(),
                                                      mesh1->m_projection);

    const std::vector<meshkernel::Point> originalNodes(mergedMesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mergedMesh->Edges());

    auto undoAction = meshkernel::ConnectMeshes::Compute(*mergedMesh);

    EXPECT_EQ(mergedMesh->GetNumValidNodes(), 15);
    EXPECT_EQ(mergedMesh->GetNumFaces(), 8);
    EXPECT_EQ(mergedMesh->GetNumValidEdges(), 22);

    const double tolerance = 1.0e-8;

    // Only comparing the nodes that were along the overlapping edge
    EXPECT_NEAR(mergedMesh->Node(9).x, 19.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(10).x, 19.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(11).x, 19.0, tolerance);

    EXPECT_NEAR(mergedMesh->Node(9).y, 0.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(10).y, 10.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(11).y, 20.0, tolerance);

    EXPECT_EQ(mergedMesh->GetNumNodesEdges(9), 3);
    EXPECT_EQ(mergedMesh->GetNumNodesEdges(10), 4);
    EXPECT_EQ(mergedMesh->GetNumNodesEdges(11), 3);

    // Test the undo action has been computed correctly
    undoAction->Restore();
    // Recompute faces
    mergedMesh->Administrate();

    ASSERT_EQ(originalNodes.size(), mergedMesh->GetNumValidNodes());
    ASSERT_EQ(originalEdges.size(), mergedMesh->GetNumValidEdges());

    meshkernel::UInt count = 0;

    for (meshkernel::UInt i = 0; i < mergedMesh->Nodes().size(); ++i)
    {
        if (mergedMesh->Node(i).IsValid())
        {
            // Check valid nodes
            EXPECT_EQ(originalNodes[count].x, mergedMesh->Node(i).x);
            EXPECT_EQ(originalNodes[count].y, mergedMesh->Node(i).y);
            ++count;
        }
    }

    count = 0;

    for (meshkernel::UInt i = 0; i < mergedMesh->Edges().size(); ++i)
    {
        if (mergedMesh->IsValidEdge(i))
        {
            EXPECT_EQ(originalEdges[count].first, mergedMesh->GetEdge(i).first);
            EXPECT_EQ(originalEdges[count].second, mergedMesh->GetEdge(i).second);
            ++count;
        }
    }
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

    const auto mergedMesh = meshkernel::Mesh2D::Merge(mesh1->Nodes(), mesh1->Edges(),
                                                      mesh2->Nodes(), mesh2->Edges(),
                                                      mesh1->m_projection);

    const std::vector<meshkernel::Point> originalNodes(mergedMesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mergedMesh->Edges());

    [[maybe_unused]] auto undoAction = meshkernel::ConnectMeshes::Compute(*mergedMesh);

    EXPECT_EQ(mergedMesh->GetNumValidNodes(), 15);
    EXPECT_EQ(mergedMesh->GetNumValidEdges(), 22);
    EXPECT_EQ(mergedMesh->GetNumFaces(), 8);

    const double tolerance = 1.0e-10;

    // Only comparing the nodes that were along the overlapping edge
    EXPECT_NEAR(mergedMesh->Node(9).x, 20.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(10).x, 20.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(11).x, 20.0, tolerance);

    EXPECT_NEAR(mergedMesh->Node(9).y, 0.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(10).y, 10.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(11).y, 20.0, tolerance);

    EXPECT_EQ(mergedMesh->GetNumNodesEdges(9), 3);
    EXPECT_EQ(mergedMesh->GetNumNodesEdges(10), 4);
    EXPECT_EQ(mergedMesh->GetNumNodesEdges(11), 3);

    // Test the undo action has been computed correctly
    undoAction->Restore();
    // Recompute faces
    mergedMesh->Administrate();

    ASSERT_EQ(originalNodes.size(), mergedMesh->GetNumValidNodes());
    ASSERT_EQ(originalEdges.size(), mergedMesh->GetNumValidEdges());

    meshkernel::UInt count = 0;

    for (meshkernel::UInt i = 0; i < mergedMesh->Nodes().size(); ++i)
    {
        if (mergedMesh->Node(i).IsValid())
        {
            // Check valid nodes
            EXPECT_EQ(originalNodes[count].x, mergedMesh->Node(i).x);
            EXPECT_EQ(originalNodes[count].y, mergedMesh->Node(i).y);
            ++count;
        }
    }

    count = 0;

    for (meshkernel::UInt i = 0; i < mergedMesh->Edges().size(); ++i)
    {
        if (mergedMesh->IsValidEdge(i))
        {
            EXPECT_EQ(originalEdges[count].first, mergedMesh->GetEdge(i).first);
            EXPECT_EQ(originalEdges[count].second, mergedMesh->GetEdge(i).second);
            ++count;
        }
    }
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

    const auto mergedMesh = meshkernel::Mesh2D::Merge(mesh1->Nodes(), mesh1->Edges(),
                                                      mesh2->Nodes(), mesh2->Edges(),
                                                      mesh1->m_projection);

    const std::vector<meshkernel::Point> originalNodes(mergedMesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mergedMesh->Edges());

    auto undoAction = meshkernel::ConnectMeshes::Compute(*mergedMesh);

    EXPECT_EQ(mergedMesh->GetNumValidNodes(), 15);
    EXPECT_EQ(mergedMesh->GetNumValidEdges(), 22);
    EXPECT_EQ(mergedMesh->GetNumFaces(), 8);

    const double tolerance = 1.0e-8;

    // Only comparing the nodes that were along the overlapping edge
    EXPECT_NEAR(mergedMesh->Node(9).x, 21.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(10).x, 21.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(11).x, 21.0, tolerance);

    EXPECT_NEAR(mergedMesh->Node(9).y, 0.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(10).y, 10.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(11).y, 20.0, tolerance);

    EXPECT_EQ(mergedMesh->GetNumNodesEdges(9), 3);
    EXPECT_EQ(mergedMesh->GetNumNodesEdges(10), 4);
    EXPECT_EQ(mergedMesh->GetNumNodesEdges(11), 3);

    // Test the undo action has been computed correctly
    undoAction->Restore();
    // Recompute faces
    mergedMesh->Administrate();

    ASSERT_EQ(originalNodes.size(), mergedMesh->GetNumValidNodes());
    ASSERT_EQ(originalEdges.size(), mergedMesh->GetNumValidEdges());

    meshkernel::UInt count = 0;

    for (meshkernel::UInt i = 0; i < mergedMesh->Nodes().size(); ++i)
    {
        if (mergedMesh->Node(i).IsValid())
        {
            // Check valid nodes
            EXPECT_EQ(originalNodes[count].x, mergedMesh->Node(i).x);
            EXPECT_EQ(originalNodes[count].y, mergedMesh->Node(i).y);
            ++count;
        }
    }

    count = 0;

    for (meshkernel::UInt i = 0; i < mergedMesh->Edges().size(); ++i)
    {
        if (mergedMesh->IsValidEdge(i))
        {
            EXPECT_EQ(originalEdges[count].first, mergedMesh->GetEdge(i).first);
            EXPECT_EQ(originalEdges[count].second, mergedMesh->GetEdge(i).second);
            ++count;
        }
    }
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

    const auto mergedMesh = meshkernel::Mesh2D::Merge(mesh1->Nodes(), mesh1->Edges(),
                                                      mesh2->Nodes(), mesh2->Edges(),
                                                      mesh1->m_projection);

    const std::vector<meshkernel::Point> originalNodes(mergedMesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mergedMesh->Edges());

    // Need to increase the search distance fraction
    auto undoAction = meshkernel::ConnectMeshes::Compute(*mergedMesh);

    // 8 triangles and 16 quadrilaterals
    EXPECT_EQ(mergedMesh->GetNumFaces(), 24);
    EXPECT_EQ(mergedMesh->GetNumValidNodes(), 31);
    EXPECT_EQ(mergedMesh->GetNumValidEdges(), 54);

    const double tolerance = 1.0e-10;

    // Only comparing the nodes that were along the irregular edge and the
    // edge where the node was created in order to free the hanging nodes
    EXPECT_NEAR(mergedMesh->Node(8).x, 20.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(9).x, 20.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(16).x, 29.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(17).x, 29.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(18).x, 29.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(19).x, 29.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(32).x, 20.0, tolerance);

    EXPECT_NEAR(mergedMesh->Node(8).y, 0.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(9).y, 10.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(16).y, 0.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(17).y, 3.1, tolerance);
    EXPECT_NEAR(mergedMesh->Node(18).y, 6.2, tolerance);
    EXPECT_NEAR(mergedMesh->Node(19).y, 9.3, tolerance);
    EXPECT_NEAR(mergedMesh->Node(32).y, 5.0, tolerance);

    const meshkernel::UInt nullValue = meshkernel::constants::missing::uintValue;

    // Allocate enough space for all edge, but will only check the edges around what was the irregular edge
    std::vector<meshkernel::Edge> expectedEdges(mergedMesh->GetNumEdges(), {nullValue, nullValue});

    expectedEdges[0].first = 0;
    expectedEdges[0].second = 4;
    expectedEdges[1].first = 1;
    expectedEdges[1].second = 5;
    expectedEdges[2].first = 2;
    expectedEdges[2].second = 6;
    expectedEdges[3].first = 3;
    expectedEdges[3].second = 7;
    expectedEdges[4].first = 4;
    expectedEdges[4].second = 8;
    expectedEdges[5].first = 5;
    expectedEdges[5].second = 9;
    expectedEdges[6].first = 6;
    expectedEdges[6].second = 10;
    expectedEdges[7].first = 7;
    expectedEdges[7].second = 11;
    expectedEdges[8].first = 8;
    expectedEdges[8].second = 16;
    expectedEdges[9].first = 9;
    expectedEdges[9].second = 19;
    expectedEdges[10].first = 10;
    expectedEdges[10].second = 14;
    expectedEdges[11].first = 11;
    expectedEdges[11].second = 15;
    expectedEdges[12].first = 1;
    expectedEdges[12].second = 0;
    expectedEdges[13].first = 2;
    expectedEdges[13].second = 1;
    expectedEdges[14].first = 3;
    expectedEdges[14].second = 2;
    expectedEdges[15].first = 5;
    expectedEdges[15].second = 4;
    expectedEdges[16].first = 6;
    expectedEdges[16].second = 5;
    expectedEdges[17].first = 7;
    expectedEdges[17].second = 6;
    expectedEdges[19].first = 10;
    expectedEdges[19].second = 9;
    expectedEdges[20].first = 11;
    expectedEdges[20].second = 10;
    expectedEdges[22].first = 14;
    expectedEdges[22].second = 19;
    expectedEdges[23].first = 15;
    expectedEdges[23].second = 14;
    expectedEdges[24].first = 16;
    expectedEdges[24].second = 20;
    expectedEdges[25].first = 17;
    expectedEdges[25].second = 21;
    expectedEdges[26].first = 18;
    expectedEdges[26].second = 22;
    expectedEdges[27].first = 19;
    expectedEdges[27].second = 23;
    expectedEdges[28].first = 20;
    expectedEdges[28].second = 24;
    expectedEdges[29].first = 21;
    expectedEdges[29].second = 25;
    expectedEdges[30].first = 22;
    expectedEdges[30].second = 26;
    expectedEdges[31].first = 23;
    expectedEdges[31].second = 27;
    expectedEdges[32].first = 24;
    expectedEdges[32].second = 28;
    expectedEdges[33].first = 25;
    expectedEdges[33].second = 29;
    expectedEdges[34].first = 26;
    expectedEdges[34].second = 30;
    expectedEdges[35].first = 27;
    expectedEdges[35].second = 31;
    expectedEdges[36].first = 17;
    expectedEdges[36].second = 16;
    expectedEdges[37].first = 18;
    expectedEdges[37].second = 17;
    expectedEdges[38].first = 19;
    expectedEdges[38].second = 18;
    expectedEdges[39].first = 21;
    expectedEdges[39].second = 20;
    expectedEdges[40].first = 22;
    expectedEdges[40].second = 21;
    expectedEdges[41].first = 23;
    expectedEdges[41].second = 22;
    expectedEdges[42].first = 25;
    expectedEdges[42].second = 24;
    expectedEdges[43].first = 26;
    expectedEdges[43].second = 25;
    expectedEdges[44].first = 27;
    expectedEdges[44].second = 26;
    expectedEdges[45].first = 29;
    expectedEdges[45].second = 28;
    expectedEdges[46].first = 30;
    expectedEdges[46].second = 29;
    expectedEdges[47].first = 31;
    expectedEdges[47].second = 30;
    expectedEdges[48].first = 18;
    expectedEdges[48].second = 32;
    expectedEdges[49].first = 18;
    expectedEdges[49].second = 9;
    expectedEdges[50].first = 17;
    expectedEdges[50].second = 32;
    expectedEdges[51].first = 17;
    expectedEdges[51].second = 8;
    expectedEdges[52].first = 32;
    expectedEdges[52].second = 9;
    expectedEdges[53].first = 32;
    expectedEdges[53].second = 8;
    expectedEdges[54].first = 32;
    expectedEdges[54].second = 5;
    expectedEdges[55].first = 32;
    expectedEdges[55].second = 4;

    for (meshkernel::UInt i = 0; i < mergedMesh->GetNumEdges(); ++i)
    {
        if (expectedEdges[i].first != nullValue)
        {
            EXPECT_EQ(expectedEdges[i].first, mergedMesh->GetEdge(i).first);
            EXPECT_EQ(expectedEdges[i].second, mergedMesh->GetEdge(i).second);
        }
    }

    // Test the undo action has been computed correctly
    undoAction->Restore();
    // Recompute faces
    mergedMesh->Administrate();

    ASSERT_EQ(originalNodes.size(), mergedMesh->GetNumValidNodes());
    ASSERT_EQ(originalEdges.size(), mergedMesh->GetNumValidEdges());

    meshkernel::UInt count = 0;

    for (meshkernel::UInt i = 0; i < mergedMesh->Nodes().size(); ++i)
    {
        if (mergedMesh->Node(i).IsValid())
        {
            // Check valid nodes
            EXPECT_EQ(originalNodes[count].x, mergedMesh->Node(i).x);
            EXPECT_EQ(originalNodes[count].y, mergedMesh->Node(i).y);
            ++count;
        }
    }

    count = 0;

    for (meshkernel::UInt i = 0; i < mergedMesh->Edges().size(); ++i)
    {
        if (mergedMesh->IsValidEdge(i))
        {
            EXPECT_EQ(originalEdges[count].first, mergedMesh->GetEdge(i).first);
            EXPECT_EQ(originalEdges[count].second, mergedMesh->GetEdge(i).second);
            ++count;
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

    const auto mergedMesh = meshkernel::Mesh2D::Merge(mesh1->Nodes(), mesh1->Edges(),
                                                      mesh2->Nodes(), mesh2->Edges(),
                                                      mesh1->m_projection);

    const std::vector<meshkernel::Point> originalNodes(mergedMesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mergedMesh->Edges());

    // Need to increase the search distance fraction
    auto undoAction = meshkernel::ConnectMeshes::Compute(*mergedMesh);

    // 8 triangles and 16 quadrilaterals
    EXPECT_EQ(mergedMesh->GetNumFaces(), 24);
    EXPECT_EQ(mergedMesh->GetNumValidNodes(), 31);
    EXPECT_EQ(mergedMesh->GetNumValidEdges(), 54);

    const double tolerance = 1.0e-8;

    // Only comparing the nodes that were along the irregular edge and the
    // edge where the node was created in order to free the hanging nodes
    EXPECT_NEAR(mergedMesh->Node(8).x, 20.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(9).x, 20.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(16).x, 31.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(17).x, 31.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(18).x, 31.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(19).x, 31.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(32).x, 20.0, tolerance);

    EXPECT_NEAR(mergedMesh->Node(8).y, 0.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(9).y, 10.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(16).y, 0.0, tolerance);
    EXPECT_NEAR(mergedMesh->Node(17).y, 3.1, tolerance);
    EXPECT_NEAR(mergedMesh->Node(18).y, 6.2, tolerance);
    EXPECT_NEAR(mergedMesh->Node(19).y, 9.3, tolerance);
    EXPECT_NEAR(mergedMesh->Node(32).y, 5.0, tolerance);

    const meshkernel::UInt nullValue = meshkernel::constants::missing::uintValue;

    // Allocate enough space for all edge, but will only check the edges around what was the irregular edge
    std::vector<meshkernel::Edge> expectedEdges(mergedMesh->GetNumEdges(), {nullValue, nullValue});

    expectedEdges[0].first = 0;
    expectedEdges[0].second = 4;
    expectedEdges[1].first = 1;
    expectedEdges[1].second = 5;
    expectedEdges[2].first = 2;
    expectedEdges[2].second = 6;
    expectedEdges[3].first = 3;
    expectedEdges[3].second = 7;
    expectedEdges[4].first = 4;
    expectedEdges[4].second = 8;
    expectedEdges[5].first = 5;
    expectedEdges[5].second = 9;
    expectedEdges[6].first = 6;
    expectedEdges[6].second = 10;
    expectedEdges[7].first = 7;
    expectedEdges[7].second = 11;
    expectedEdges[8].first = 8;
    expectedEdges[8].second = 16;
    expectedEdges[9].first = 9;
    expectedEdges[9].second = 19;
    expectedEdges[10].first = 10;
    expectedEdges[10].second = 14;
    expectedEdges[11].first = 11;
    expectedEdges[11].second = 15;
    expectedEdges[12].first = 1;
    expectedEdges[12].second = 0;
    expectedEdges[13].first = 2;
    expectedEdges[13].second = 1;
    expectedEdges[14].first = 3;
    expectedEdges[14].second = 2;
    expectedEdges[15].first = 5;
    expectedEdges[15].second = 4;
    expectedEdges[16].first = 6;
    expectedEdges[16].second = 5;
    expectedEdges[17].first = 7;
    expectedEdges[17].second = 6;
    expectedEdges[19].first = 10;
    expectedEdges[19].second = 9;
    expectedEdges[20].first = 11;
    expectedEdges[20].second = 10;
    expectedEdges[22].first = 14;
    expectedEdges[22].second = 19;
    expectedEdges[23].first = 15;
    expectedEdges[23].second = 14;
    expectedEdges[24].first = 16;
    expectedEdges[24].second = 20;
    expectedEdges[25].first = 17;
    expectedEdges[25].second = 21;
    expectedEdges[26].first = 18;
    expectedEdges[26].second = 22;
    expectedEdges[27].first = 19;
    expectedEdges[27].second = 23;
    expectedEdges[28].first = 20;
    expectedEdges[28].second = 24;
    expectedEdges[29].first = 21;
    expectedEdges[29].second = 25;
    expectedEdges[30].first = 22;
    expectedEdges[30].second = 26;
    expectedEdges[31].first = 23;
    expectedEdges[31].second = 27;
    expectedEdges[32].first = 24;
    expectedEdges[32].second = 28;
    expectedEdges[33].first = 25;
    expectedEdges[33].second = 29;
    expectedEdges[34].first = 26;
    expectedEdges[34].second = 30;
    expectedEdges[35].first = 27;
    expectedEdges[35].second = 31;
    expectedEdges[36].first = 17;
    expectedEdges[36].second = 16;
    expectedEdges[37].first = 18;
    expectedEdges[37].second = 17;
    expectedEdges[38].first = 19;
    expectedEdges[38].second = 18;
    expectedEdges[39].first = 21;
    expectedEdges[39].second = 20;
    expectedEdges[40].first = 22;
    expectedEdges[40].second = 21;
    expectedEdges[41].first = 23;
    expectedEdges[41].second = 22;
    expectedEdges[42].first = 25;
    expectedEdges[42].second = 24;
    expectedEdges[43].first = 26;
    expectedEdges[43].second = 25;
    expectedEdges[44].first = 27;
    expectedEdges[44].second = 26;
    expectedEdges[45].first = 29;
    expectedEdges[45].second = 28;
    expectedEdges[46].first = 30;
    expectedEdges[46].second = 29;
    expectedEdges[47].first = 31;
    expectedEdges[47].second = 30;
    expectedEdges[48].first = 18;
    expectedEdges[48].second = 32;
    expectedEdges[49].first = 18;
    expectedEdges[49].second = 9;
    expectedEdges[50].first = 17;
    expectedEdges[50].second = 32;
    expectedEdges[51].first = 17;
    expectedEdges[51].second = 8;
    expectedEdges[52].first = 32;
    expectedEdges[52].second = 9;
    expectedEdges[53].first = 32;
    expectedEdges[53].second = 8;
    expectedEdges[54].first = 32;
    expectedEdges[54].second = 5;
    expectedEdges[55].first = 32;
    expectedEdges[55].second = 4;

    for (meshkernel::UInt i = 0; i < mergedMesh->GetNumEdges(); ++i)
    {
        if (expectedEdges[i].first != nullValue)
        {
            EXPECT_EQ(expectedEdges[i].first, mergedMesh->GetEdge(i).first);
            EXPECT_EQ(expectedEdges[i].second, mergedMesh->GetEdge(i).second);
        }
    }

    // Test the undo action has been computed correctly
    undoAction->Restore();
    // Recompute faces
    mergedMesh->Administrate();

    ASSERT_EQ(originalNodes.size(), mergedMesh->GetNumValidNodes());
    ASSERT_EQ(originalEdges.size(), mergedMesh->GetNumValidEdges());

    meshkernel::UInt count = 0;

    for (meshkernel::UInt i = 0; i < mergedMesh->Nodes().size(); ++i)
    {
        if (mergedMesh->Node(i).IsValid())
        {
            // Check valid nodes
            EXPECT_EQ(originalNodes[count].x, mergedMesh->Node(i).x);
            EXPECT_EQ(originalNodes[count].y, mergedMesh->Node(i).y);
            ++count;
        }
    }

    count = 0;

    for (meshkernel::UInt i = 0; i < mergedMesh->Edges().size(); ++i)
    {
        if (mergedMesh->IsValidEdge(i))
        {
            EXPECT_EQ(originalEdges[count].first, mergedMesh->GetEdge(i).first);
            EXPECT_EQ(originalEdges[count].second, mergedMesh->GetEdge(i).second);
            ++count;
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

    const auto mergedMesh = meshkernel::Mesh2D::Merge(mesh1->Nodes(), mesh1->Edges(),
                                                      mesh2->Nodes(), mesh2->Edges(),
                                                      mesh1->m_projection);

    EXPECT_THROW([[maybe_unused]] auto undoAction = meshkernel::ConnectMeshes::Compute(*mergedMesh, 0.5), meshkernel::RangeError);
    EXPECT_THROW([[maybe_unused]] auto undoAction = meshkernel::ConnectMeshes::Compute(*mergedMesh, 1.5), meshkernel::RangeError);
    EXPECT_THROW([[maybe_unused]] auto undoAction = meshkernel::ConnectMeshes::Compute(*mergedMesh, -0.5), meshkernel::RangeError);
    EXPECT_THROW([[maybe_unused]] auto undoAction = meshkernel::ConnectMeshes::Compute(*mergedMesh, 0.0), meshkernel::RangeError);
}

TEST(Mesh2DConnectDD, SimpleMergeMeshes)
{

    // Should only concatinate the node and edge arrays.

    meshkernel::Point origin{0.0, 0.0};
    meshkernel::Vector delta{10.0, 10.0};

    std::shared_ptr<meshkernel::Mesh2D> mesh1 = generateMesh(origin, delta, 3, 3);

    origin.x += 2.0 * delta.x();

    delta = meshkernel::Vector{4.0, 4.0};
    std::shared_ptr<meshkernel::Mesh2D> mesh2 = generateMesh(origin, delta, 3, 5);

    const auto mergedMesh = meshkernel::Mesh2D::Merge(mesh1->Nodes(), mesh1->Edges(),
                                                      mesh2->Nodes(), mesh2->Edges(),
                                                      mesh1->m_projection);

    const std::vector<meshkernel::Point> expectedNodes{{0, 0}, {0, 10}, {0, 20}, {10, 0}, {10, 10}, {10, 20}, {20, 0}, {20, 10}, {20, 20}, {20, 0}, {20, 4}, {20, 8}, {20, 12}, {20, 16}, {24, 0}, {24, 4}, {24, 8}, {24, 12}, {24, 16}, {28, 0}, {28, 4}, {28, 8}, {28, 12}, {28, 16}};

    const std::vector<meshkernel::Edge> expectedEdges{{0, 3}, {1, 4}, {2, 5}, {3, 6}, {4, 7}, {5, 8}, {1, 0}, {2, 1}, {4, 3}, {5, 4}, {7, 6}, {8, 7}, {9, 14}, {10, 15}, {11, 16}, {12, 17}, {13, 18}, {14, 19}, {15, 20}, {16, 21}, {17, 22}, {18, 23}, {10, 9}, {11, 10}, {12, 11}, {13, 12}, {15, 14}, {16, 15}, {17, 16}, {18, 17}, {20, 19}, {21, 20}, {22, 21}, {23, 22}};

    const auto& nodes = mergedMesh->Nodes();
    const auto& edges = mergedMesh->Edges();

    ASSERT_EQ(expectedNodes.size(), nodes.size());
    ASSERT_EQ(expectedEdges.size(), edges.size());

    const double tolerance = 1.0e-10;

    for (size_t i = 0; i < nodes.size(); ++i)
    {
        EXPECT_NEAR(expectedNodes[i].x, nodes[i].x, tolerance);
        EXPECT_NEAR(expectedNodes[i].y, nodes[i].y, tolerance);
    }

    for (size_t i = 0; i < edges.size(); ++i)
    {
        EXPECT_EQ(expectedEdges[i].first, edges[i].first);
        EXPECT_EQ(expectedEdges[i].second, edges[i].second);
    }
}

TEST(Mesh2DConnectDD, ConnectMesheInsidePolygon)
{

    // Should only concatinate the node and edge arrays.

    meshkernel::Point origin{0.0, 0.0};
    meshkernel::Vector delta{10.0, 10.0};

    std::shared_ptr<meshkernel::Mesh2D> mesh1 = generateMesh(origin, delta, 5, 5);

    origin.x += 4.0 * delta.x();

    delta = meshkernel::Vector{2.5, 2.5};
    std::shared_ptr<meshkernel::Mesh2D> mesh2 = generateMesh(origin, delta, 3, 17);

    [[maybe_unused]] meshkernel::Polygons polygon({{35.0, -5.0}, {50.0, -5.0}, {50.0, 21.25}, {35.0, 21.25}, {35.0, -5.0}}, mesh1->m_projection);

    const auto mergedMesh = meshkernel::Mesh2D::Merge(mesh1->Nodes(), mesh1->Edges(),
                                                      mesh2->Nodes(), mesh2->Edges(),
                                                      mesh1->m_projection);

    auto undoAction = meshkernel::ConnectMeshes::Compute(*mergedMesh, polygon);

    const std::vector<meshkernel::Point> expectedNodes{{0.0, 0.0}, {0.0, 10.0}, {0.0, 20.0}, {0.0, 30.0}, {0.0, 40.0}, {10.0, 0.0}, {10.0, 10.0}, {10.0, 20.0}, {10.0, 30.0}, {10.0, 40.0}, {20.0, 0.0}, {20.0, 10.0}, {20.0, 20.0}, {20.0, 30.0}, {20.0, 40.0}, {30.0, 0.0}, {30.0, 10.0}, {30.0, 20.0}, {30.0, 30.0}, {30.0, 40.0}, {meshkernel::constants::missing::doubleValue, meshkernel::constants::missing::doubleValue}, {meshkernel::constants::missing::doubleValue, meshkernel::constants::missing::doubleValue}, {meshkernel::constants::missing::doubleValue, meshkernel::constants::missing::doubleValue}, {40.0, 30.0}, {40.0, 40.0}, {40.0, 0.0}, {40.0, 2.5}, {40.0, 5.0}, {40.0, 7.5}, {40.0, 10.0}, {40.0, 12.5}, {40.0, 15.0}, {40.0, 17.5}, {40.0, 20.0}, {40.0, 22.5}, {40.0, 25.0}, {40.0, 27.5}, {40.0, 30.0}, {40.0, 32.5}, {40.0, 35.0}, {40.0, 37.5}, {40.0, 40.0}, {42.5, 0.0}, {42.5, 2.5}, {42.5, 5.0}, {42.5, 7.5}, {42.5, 10.0}, {42.5, 12.5}, {42.5, 15.0}, {42.5, 17.5}, {42.5, 20.0}, {42.5, 22.5}, {42.5, 25.0}, {42.5, 27.5}, {42.5, 30.0}, {42.5, 32.5}, {42.5, 35.0}, {42.5, 37.5}, {42.5, 40.0}, {45.0, 0.0}, {45.0, 2.5}, {45.0, 5.0}, {45.0, 7.5}, {45.0, 10.0}, {45.0, 12.5}, {45.0, 15.0}, {45.0, 17.5}, {45.0, 20.0}, {45.0, 22.5}, {45.0, 25.0}, {45.0, 27.5}, {45.0, 30.0}, {45.0, 32.5}, {45.0, 35.0}, {45.0, 37.5}, {45.0, 40.0}, {30.0, 5.0}, {30.0, 15.0}};

    const std::vector<meshkernel::Edge> expectedEdges{{0, 5}, {1, 6}, {2, 7}, {3, 8}, {4, 9}, {5, 10}, {6, 11}, {7, 12}, {8, 13}, {9, 14}, {10, 15}, {11, 16}, {12, 17}, {13, 18}, {14, 19}, {15, 25}, {16, 29}, {17, 33}, {18, 23}, {19, 24}, {1, 0}, {2, 1}, {3, 2}, {4, 3}, {6, 5}, {7, 6}, {8, 7}, {9, 8}, {11, 10}, {12, 11}, {13, 12}, {14, 13}, {meshkernel::constants::missing::uintValue, meshkernel::constants::missing::uintValue}, {meshkernel::constants::missing::uintValue, meshkernel::constants::missing::uintValue}, {18, 17}, {19, 18}, {meshkernel::constants::missing::uintValue, meshkernel::constants::missing::uintValue}, {meshkernel::constants::missing::uintValue, meshkernel::constants::missing::uintValue}, {23, 33}, {24, 23}, {25, 42}, {26, 43}, {27, 44}, {28, 45}, {29, 46}, {30, 47}, {31, 48}, {32, 49}, {33, 50}, {34, 51}, {35, 52}, {36, 53}, {37, 54}, {38, 55}, {39, 56}, {40, 57}, {41, 58}, {42, 59}, {43, 60}, {44, 61}, {45, 62}, {46, 63}, {47, 64}, {48, 65}, {49, 66}, {50, 67}, {51, 68}, {52, 69}, {53, 70}, {54, 71}, {55, 72}, {56, 73}, {57, 74}, {58, 75}, {26, 25}, {27, 26}, {28, 27}, {29, 28}, {30, 29}, {31, 30}, {32, 31}, {33, 32}, {34, 33}, {35, 34}, {36, 35}, {37, 36}, {38, 37}, {39, 38}, {40, 39}, {41, 40}, {43, 42}, {44, 43}, {45, 44}, {46, 45}, {47, 46}, {48, 47}, {49, 48}, {50, 49}, {51, 50}, {52, 51}, {53, 52}, {54, 53}, {55, 54}, {56, 55}, {57, 56}, {58, 57}, {60, 59}, {61, 60}, {62, 61}, {63, 62}, {64, 63}, {65, 64}, {66, 65}, {67, 66}, {68, 67}, {69, 68}, {70, 69}, {71, 70}, {72, 71}, {73, 72}, {74, 73}, {75, 74}, {27, 76}, {76, 15}, {76, 16}, {28, 76}, {28, 16}, {26, 76}, {26, 15}, {76, 11}, {76, 10}, {31, 77}, {77, 16}, {77, 17}, {32, 77}, {32, 17}, {30, 77}, {30, 16}, {77, 12}, {77, 11}};

    const auto& nodes = mergedMesh->Nodes();
    const auto& edges = mergedMesh->Edges();

    ASSERT_EQ(expectedNodes.size(), nodes.size());
    ASSERT_EQ(expectedEdges.size(), edges.size());

    const double tolerance = 1.0e-10;

    for (size_t i = 0; i < nodes.size(); ++i)
    {
        EXPECT_NEAR(expectedNodes[i].x, nodes[i].x, tolerance);
        EXPECT_NEAR(expectedNodes[i].y, nodes[i].y, tolerance);
    }

    for (size_t i = 0; i < edges.size(); ++i)
    {
        EXPECT_EQ(expectedEdges[i].first, edges[i].first);
        EXPECT_EQ(expectedEdges[i].second, edges[i].second);
    }
}
