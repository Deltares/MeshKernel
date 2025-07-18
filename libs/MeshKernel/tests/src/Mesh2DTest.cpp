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
#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>
#include <random>

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Mesh.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/Mesh2DIntersections.hpp"
#include "MeshKernel/Mesh2DToCurvilinear.hpp"
#include "MeshKernel/MeshEdgeLength.hpp"
#include "MeshKernel/MeshOrthogonality.hpp"
#include "MeshKernel/MeshSmoothness.hpp"
#include "MeshKernel/Operations.hpp"
#include "MeshKernel/Polygons.hpp"
#include "MeshKernel/RemoveDisconnectedRegions.hpp"
#include "MeshKernel/Utilities/Utilities.hpp"

#include "TestUtils/Definitions.hpp"
#include "TestUtils/MakeCurvilinearGrids.hpp"
#include "TestUtils/MakeMeshes.hpp"

TEST(Mesh2D, OneQuadTestConstructor)
{
    // 1 Setup
    std::vector<meshkernel::Point> nodes;
    nodes.push_back({0.0, 0.0});
    nodes.push_back({0.0, 10.0});
    nodes.push_back({10.0, 0.0});
    nodes.push_back({10.0, 10.0});
    std::vector<meshkernel::Edge> edges;
    edges.push_back({0, 2});
    edges.push_back({1, 3});
    edges.push_back({0, 1});
    edges.push_back({2, 3});

    // 2 Execution
    const auto mesh = meshkernel::Mesh2D(edges, nodes, meshkernel::Projection::cartesian);

    // 3 Validation
    // expect nodesEdges to be sorted ccw
    ASSERT_EQ(0, mesh.m_nodesEdges[0][0]);
    ASSERT_EQ(2, mesh.m_nodesEdges[0][1]);

    ASSERT_EQ(1, mesh.m_nodesEdges[1][0]);
    ASSERT_EQ(2, mesh.m_nodesEdges[1][1]);

    ASSERT_EQ(0, mesh.m_nodesEdges[2][0]);
    ASSERT_EQ(3, mesh.m_nodesEdges[2][1]);

    ASSERT_EQ(1, mesh.m_nodesEdges[3][0]);
    ASSERT_EQ(3, mesh.m_nodesEdges[3][1]);

    // each node has two edges int this case
    ASSERT_EQ(2, mesh.GetNumNodesEdges(0));
    ASSERT_EQ(2, mesh.GetNumNodesEdges(1));
    ASSERT_EQ(2, mesh.GetNumNodesEdges(2));
    ASSERT_EQ(2, mesh.GetNumNodesEdges(3));

    // the nodes composing the face, in ccw order
    ASSERT_EQ(0, mesh.m_facesNodes[0][0]);
    ASSERT_EQ(2, mesh.m_facesNodes[0][1]);
    ASSERT_EQ(3, mesh.m_facesNodes[0][2]);
    ASSERT_EQ(1, mesh.m_facesNodes[0][3]);

    // the edges composing the face, in ccw order
    ASSERT_EQ(0, mesh.m_facesEdges[0][0]);
    ASSERT_EQ(3, mesh.m_facesEdges[0][1]);
    ASSERT_EQ(1, mesh.m_facesEdges[0][2]);
    ASSERT_EQ(2, mesh.m_facesEdges[0][3]);

    // // the found circumcenter for the face
    // ASSERT_DOUBLE_EQ(5.0, mesh.m_facesCircumcenters[0].x);
    // ASSERT_DOUBLE_EQ(5.0, mesh.m_facesCircumcenters[0].y);

    // each edge has only one face in this case
    ASSERT_EQ(1, mesh.GetNumEdgesFaces(0));
    ASSERT_EQ(1, mesh.GetNumEdgesFaces(1));
    ASSERT_EQ(1, mesh.GetNumEdgesFaces(2));
    ASSERT_EQ(1, mesh.GetNumEdgesFaces(3));

    // each edge is a boundary edge, so the second entry of edgesFaces is an invalid index (meshkernel::constants::missing::sizetValue)
    ASSERT_EQ(meshkernel::constants::missing::uintValue, mesh.m_edgesFaces[0][1]);
    ASSERT_EQ(meshkernel::constants::missing::uintValue, mesh.m_edgesFaces[1][1]);
    ASSERT_EQ(meshkernel::constants::missing::uintValue, mesh.m_edgesFaces[2][1]);
    ASSERT_EQ(meshkernel::constants::missing::uintValue, mesh.m_edgesFaces[3][1]);
}

TEST(Mesh2D, TriangulateSamples)
{
    // Prepare
    std::vector<meshkernel::Point> nodes;

    nodes.push_back({498.503152894023, 1645.82297461613});
    nodes.push_back({-5.90937355559299, 814.854361678898});
    nodes.push_back({851.30035347439, 150.079471329115});
    nodes.push_back({1411.11078745316, 1182.22995897746});
    nodes.push_back({501.418832237663, 1642.90729527249});
    nodes.push_back({498.503152894023, 1645.82297461613});

    meshkernel::Polygons polygons(nodes, meshkernel::Projection::cartesian);

    // Execute
    const auto generatedPoints = polygons.ComputePointsInPolygons();

    meshkernel::Mesh2D mesh(generatedPoints[0], polygons, meshkernel::Projection::cartesian);
}

TEST(Mesh2D, ComputePointsInPolygonGeneratesIncludedPoints)
{
    // Generate one polygon with a defined spacing
    std::vector<meshkernel::Point> polygonNodes{{-657.485056492302, 2363.36631085451},
                                                {-755.085460985301, 2079.6866305057},
                                                {-852.6858654783, 1796.00695015689},
                                                {-950.286269971299, 1512.32726980808},
                                                {-1047.8866744643, 1228.64758945927},
                                                {-1145.4870789573, 944.967909110461},
                                                {-1195.26299221778, 800.292123839341},
                                                {-896.081193037718, 778.150510542622},
                                                {-596.899393857656, 756.008897245902},
                                                {-297.717594677595, 733.867283949183},
                                                {1.46420450246637, 711.725670652464},
                                                {300.646003682528, 689.584057355744},
                                                {599.827802862589, 667.442444059025},
                                                {899.00960204265, 645.300830762306},
                                                {1198.19140122271, 623.159217465586},
                                                {1497.37320040277, 601.017604168867},
                                                {1589.119964716, 594.227681178177},
                                                {1352.01769916658, 778.027887030438},
                                                {1114.91543361716, 961.8280928827},
                                                {877.813168067745, 1145.62829873496},
                                                {640.710902518328, 1329.42850458722},
                                                {403.608636968911, 1513.22871043948},
                                                {292.421764555502, 1599.42008440337},
                                                {555.079172252551, 1744.37209777975},
                                                {817.7365799496, 1889.32411115613},
                                                {1080.39398764665, 2034.27612453251},
                                                {1343.0513953437, 2179.22813790889},
                                                {1503.67861044186, 2267.87303254812},
                                                {1203.97104330649, 2281.11592504945},
                                                {904.263476171129, 2294.35881755078},
                                                {604.555909035765, 2307.60171005211},
                                                {304.848341900402, 2320.84460255344},
                                                {5.14077476503894, 2334.08749505477},
                                                {-294.566792370324, 2347.3303875561},
                                                {-594.274359505687, 2360.57328005743},
                                                {-657.485056492302, 2363.36631085451}};

    meshkernel::Polygons polygons(polygonNodes, meshkernel::Projection::cartesian);

    // Execute
    const auto generatedPoints = polygons.ComputePointsInPolygons();

    // Check all points are included in the first polygon
    for (const auto& p : generatedPoints[0])
    {
        ASSERT_TRUE(polygons.IsPointInPolygon(p, 0));
    }
}

TEST(Mesh2D, TwoTrianglesDuplicatedEdges)
{
    // 1 Setup
    std::vector<meshkernel::Point> nodes;
    nodes.push_back({0.0, 0.0});
    nodes.push_back({5.0, -5.0});
    nodes.push_back({10.0, 0.0});
    nodes.push_back({5.0, 5.0});
    std::vector<meshkernel::Edge> edges;
    edges.push_back({0, 3});
    edges.push_back({0, 2});
    edges.push_back({2, 3});
    edges.push_back({0, 1});
    edges.push_back({2, 1});

    // 2 Execution
    const auto mesh = meshkernel::Mesh2D(edges, nodes, meshkernel::Projection::cartesian);

    // 3 Validation
    ASSERT_EQ(2, mesh.GetNumFaces());
}

TEST(Mesh2D, MeshBoundaryToPolygon)
{
    // 1 Setup
    std::vector<meshkernel::Point> nodes;
    nodes.push_back({0.0, 0.0});
    nodes.push_back({5.0, -5.0});
    nodes.push_back({10.0, 0.0});
    nodes.push_back({5.0, 5.0});
    std::vector<meshkernel::Edge> edges;
    edges.push_back({0, 3});
    edges.push_back({0, 2});
    edges.push_back({2, 3});
    edges.push_back({0, 1});
    edges.push_back({2, 1});

    auto mesh = meshkernel::Mesh2D(edges, nodes, meshkernel::Projection::cartesian);

    std::vector<meshkernel::Point> polygonNodes;
    const auto meshBoundaryPolygon = mesh.ComputeBoundaryPolygons(polygonNodes);

    const double tolerance = 1e-5;
    ASSERT_NEAR(0.0, meshBoundaryPolygon[0].x, tolerance);
    ASSERT_NEAR(5.0, meshBoundaryPolygon[1].x, tolerance);
    ASSERT_NEAR(10.0, meshBoundaryPolygon[2].x, tolerance);
    ASSERT_NEAR(5.0, meshBoundaryPolygon[3].x, tolerance);
    ASSERT_NEAR(0.0, meshBoundaryPolygon[4].x, tolerance);

    ASSERT_NEAR(0.0, meshBoundaryPolygon[0].y, tolerance);
    ASSERT_NEAR(5.0, meshBoundaryPolygon[1].y, tolerance);
    ASSERT_NEAR(0.0, meshBoundaryPolygon[2].y, tolerance);
    ASSERT_NEAR(-5.0, meshBoundaryPolygon[3].y, tolerance);
    ASSERT_NEAR(0.0, meshBoundaryPolygon[4].y, tolerance);
}

TEST(Mesh2D, HangingEdge)
{
    // 1 Setup
    std::vector<meshkernel::Point> nodes;
    nodes.push_back({0.0, 0.0});
    nodes.push_back({5.0, 0.0});
    nodes.push_back({3.0, 2.0});
    nodes.push_back({3.0, 4.0});

    std::vector<meshkernel::Edge> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 3});
    edges.push_back({3, 0});
    edges.push_back({2, 1});

    auto const mesh = meshkernel::Mesh2D(edges, nodes, meshkernel::Projection::cartesian);

    ASSERT_EQ(1, mesh.GetNumFaces());
}

TEST(Mesh2D, NodeMerging)
{
    // 1. Setup
    const int n = 10; // x
    const int m = 10; // y

    std::vector<std::vector<int>> indicesValues(n, std::vector<int>(m));
    std::vector<meshkernel::Point> nodes(n * m);
    meshkernel::UInt nodeIndex = 0;
    for (auto j = 0; j < m; ++j)
    {
        for (auto i = 0; i < n; ++i)
        {
            indicesValues[i][j] = i + j * n;
            nodes[nodeIndex] = {static_cast<double>(i), static_cast<double>(j)};
            nodeIndex++;
        }
    }

    std::vector<meshkernel::Edge> edges((n - 1) * m + (m - 1) * n);
    meshkernel::UInt edgeIndex = 0;
    for (meshkernel::UInt j = 0; j < m; ++j)
    {
        for (meshkernel::UInt i = 0; i < n - 1; ++i)
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

    auto mesh = std::make_unique<meshkernel::Mesh2D>(edges, nodes, meshkernel::Projection::cartesian);

    // Add overlapping nodes
    double generatingDistance = std::sqrt(std::pow(0.001 * 0.9, 2) / 2.0);
    std::uniform_real_distribution<double> x_distribution(0.0, generatingDistance);
    std::uniform_real_distribution<double> y_distribution(0.0, generatingDistance);
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());

    nodes.resize(mesh->GetNumNodes() * 2);
    edges.resize(mesh->GetNumEdges() + mesh->GetNumNodes() * 2);
    meshkernel::UInt originalNodeIndex = 0;
    for (meshkernel::UInt j = 0; j < m; ++j)
    {
        for (meshkernel::UInt i = 0; i < n; ++i)
        {
            nodes[nodeIndex] = {i + x_distribution(generator), j + y_distribution(generator)};

            edges[edgeIndex] = {nodeIndex, originalNodeIndex};
            edgeIndex++;

            nodeIndex++;
            originalNodeIndex++;
        }
    }

    nodes.resize(nodeIndex);
    edges.resize(edgeIndex);

    // re set with augmented nodes
    mesh = std::make_unique<meshkernel::Mesh2D>(edges, nodes, meshkernel::Projection::cartesian);

    // 2. Act
    meshkernel::Polygons polygon;
    [[maybe_unused]] auto action = mesh->MergeNodesInPolygon(polygon, 0.001);

    // 3. Assert
    ASSERT_EQ(mesh->GetNumValidNodes(), n * m);
    ASSERT_EQ(mesh->GetNumValidEdges(), (n - 1) * m + (m - 1) * n);
}

TEST(Mesh2D, MillionQuads)
{
    const int n = 4; // x
    const int m = 4; // y

    std::vector<std::vector<int>> indicesValues(n, std::vector<int>(m));
    std::vector<meshkernel::Point> nodes(n * m);
    std::size_t nodeIndex = 0;
    for (auto j = 0; j < m; ++j)
    {
        for (auto i = 0; i < n; ++i)
        {
            indicesValues[i][j] = i + j * n;
            nodes[nodeIndex] = {(double)i, (double)j};
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

    // now build node-edge mapping
    auto start(std::chrono::steady_clock::now());
    const auto mesh = meshkernel::Mesh2D(edges, nodes, meshkernel::Projection::cartesian);
    auto end(std::chrono::steady_clock::now());

    double elapsedTime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    std::cout << "Elapsed time " << elapsedTime << " s " << std::endl;
    std::cout << "Number of found cells " << mesh.GetNumFaces() << std::endl;

    EXPECT_LE(elapsedTime, 5.0);
}

TEST(Mesh2D, GetObtuseTriangles)
{
    // Setup a mesh with two triangles, one obtuse
    std::vector<meshkernel::Point> nodes{
        {0.0, 0.0},
        {3.0, 0.0},
        {-1.0, 2.0},
        {1.5, -2.0}};

    std::vector<meshkernel::Edge> edges{
        {0, 1},
        {1, 2},
        {2, 0},
        {0, 3},
        {3, 1}};

    const auto mesh = std::make_unique<meshkernel::Mesh2D>(edges, nodes, meshkernel::Projection::cartesian);

    // execute, only one obtuse triangle should be found
    const auto obtuseTrianglesCount = mesh->GetObtuseTrianglesCenters().size();

    // assert a small flow edge is found
    ASSERT_EQ(1, obtuseTrianglesCount);
}

TEST(Mesh2D, GetSmallFlowEdgeCenters)
{
    // Setup a mesh with two triangles
    std::vector<meshkernel::Point> nodes{
        {0.0, 0.0},
        {1.0, 0.0},
        {1.0, 0.3},
        {1.0, -0.3}};

    std::vector<meshkernel::Edge> edges{
        {0, 3},
        {3, 1},
        {1, 0},
        {1, 2},
        {2, 0},
    };

    meshkernel::Mesh2D mesh = meshkernel::Mesh2D(edges, nodes, meshkernel::Projection::cartesian);

    // execute, by setting the smallFlowEdgesThreshold high, a small flow edge will be found
    const auto numSmallFlowEdgeFirstQuery = mesh.GetEdgesCrossingSmallFlowEdges(100).size();

    // execute, by setting the smallFlowEdgesThreshold low, no small flow edge will be found
    const auto numSmallFlowEdgeSecondQuery = mesh.GetEdgesCrossingSmallFlowEdges(0.0).size();

    // assert a small flow edge is found
    ASSERT_EQ(1, numSmallFlowEdgeFirstQuery);
    ASSERT_EQ(0, numSmallFlowEdgeSecondQuery);
}

TEST(Mesh2D, DeleteSmallFlowEdge)
{
    // Setup a mesh with eight triangles
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/RemoveSmallFlowEdgesTests/remove_small_flow_edges_net.nc");

    std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    ASSERT_EQ(8, mesh->GetNumFaces());

    // After merging the number of faces is reduced
    auto undoAction = mesh->DeleteSmallFlowEdges(1.0);

    ASSERT_EQ(3, mesh->GetNumFaces());

    // Restore original mesh
    undoAction->Restore();
    mesh->Administrate();

    ASSERT_EQ(8, mesh->GetNumFaces());
    ASSERT_EQ(originalEdges.size(), mesh->GetNumEdges());
}

TEST(Mesh2D, DeleteSmallTrianglesAtBoundaries)
{
    // Setup a mesh with two triangles
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/RemoveSmallFlowEdgesTests/remove_small_flow_edges_quad_net.nc");

    ASSERT_EQ(2, mesh->GetNumFaces());

    // After merging
    auto undoAction = mesh->DeleteSmallTrianglesAtBoundaries(0.6);

    ASSERT_EQ(1, mesh->GetNumFaces());

    const double tolerance = 1e-8;
    ASSERT_NEAR(364.17013549804688, mesh->Node(0).x, tolerance);
    ASSERT_NEAR(meshkernel::constants::missing::doubleValue, mesh->Node(1).x, tolerance);
    ASSERT_NEAR(295.21142578125000, mesh->Node(2).x, tolerance);
    ASSERT_NEAR(421.46209716796875, mesh->Node(3).x, tolerance);
    ASSERT_NEAR(359.79510498046875, mesh->Node(4).x, tolerance);

    ASSERT_NEAR(374.00662231445313, mesh->Node(0).y, tolerance);
    ASSERT_NEAR(meshkernel::constants::missing::doubleValue, mesh->Node(1).y, tolerance);
    ASSERT_NEAR(300.48181152343750, mesh->Node(2).y, tolerance);
    ASSERT_NEAR(295.33038330078125, mesh->Node(3).y, tolerance);
    ASSERT_NEAR(398.59295654296875, mesh->Node(4).y, tolerance);

    // Restore original mesh
    undoAction->Restore();
    mesh->Administrate();
    ASSERT_EQ(2, mesh->GetNumFaces());
}

TEST(Mesh2D, DeleteHangingEdge)
{
    // 1 Setup
    std::vector<meshkernel::Point> nodes;
    nodes.push_back({0.0, 0.0});
    nodes.push_back({5.0, 0.0});
    nodes.push_back({3.0, 4.0});

    std::vector<meshkernel::Edge> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    edges.push_back({2, 0});

    // Execute
    const auto mesh = std::make_unique<meshkernel::Mesh2D>(edges, nodes, meshkernel::Projection::cartesian);

    // Add new node and connect with existing node, creating hanging edge
    nodes.push_back({3.0, 2.0});
    edges.push_back({3, 1});

    [[maybe_unused]] auto undoInsertNode = mesh->InsertNode(nodes[3]);
    [[maybe_unused]] auto undoConnectNodes = mesh->ConnectNodes(3, 1);

    // Assert
    ASSERT_EQ(1, mesh->GetNumFaces());
    ASSERT_EQ(4, mesh->GetNumEdges());

    // Execute
    auto hangingEdges = mesh->GetHangingEdges();

    // Assert
    ASSERT_EQ(1, hangingEdges.size());

    // Execute
    auto undoAction = mesh->DeleteHangingEdges();
    hangingEdges = mesh->GetHangingEdges();

    // Assert
    ASSERT_EQ(0, hangingEdges.size());

    // Restore original mesh
    undoAction->Restore();
    // Called to reconstruct the faces
    mesh->Administrate();

    ASSERT_EQ(1, mesh->GetNumFaces());
    ASSERT_EQ(4, mesh->GetNumValidEdges());

    for (meshkernel::UInt i = 0; i < mesh->Nodes().size(); ++i)
    {
        EXPECT_EQ(nodes[i].x, mesh->Node(i).x);
        EXPECT_EQ(nodes[i].y, mesh->Node(i).y);
    }

    for (meshkernel::UInt i = 0; i < mesh->Edges().size(); ++i)
    {
        EXPECT_EQ(edges[i].first, mesh->GetEdge(i).first);
        EXPECT_EQ(edges[i].second, mesh->GetEdge(i).second);
    }
}

TEST(Mesh2D, GetPolylineIntersectionsFromSimplePolylineShouldReturnCorrectIntersections)
{
    // 1. Setup
    auto mesh = MakeRectangularMeshForTesting(4, 4, 1.0, meshkernel::Projection::cartesian);

    std::vector<meshkernel::Point> boundaryPolygonNodes;
    boundaryPolygonNodes.emplace_back(0.5, 0.5);
    boundaryPolygonNodes.emplace_back(2.5, 0.5);
    boundaryPolygonNodes.emplace_back(2.5, 2.5);
    boundaryPolygonNodes.emplace_back(0.5, 2.5);
    boundaryPolygonNodes.emplace_back(0.5, 0.5);

    // 2. Execute
    const meshkernel::Polygons boundaryPolygon(boundaryPolygonNodes, mesh->m_projection);
    meshkernel::Mesh2DIntersections mesh2DIntersections(*mesh);
    mesh2DIntersections.Compute(boundaryPolygon);
    auto edgeIntersections = mesh2DIntersections.EdgeIntersections();
    auto faceIntersections = mesh2DIntersections.FaceIntersections();

    meshkernel::Mesh2DIntersections::sortAndEraseIntersections(edgeIntersections);
    meshkernel::Mesh2DIntersections::sortAndEraseIntersections(faceIntersections);

    // 3. Assert

    // edge intersections
    ASSERT_EQ(edgeIntersections[0].polylineSegmentIndex, 0);
    ASSERT_NEAR(edgeIntersections[0].adimensionalPolylineSegmentDistance, 0.25000000000000000, 1e-8);
    ASSERT_EQ(edgeIntersections[0].edgeIndex, 15);
    ASSERT_EQ(edgeIntersections[0].edgeFirstNode, 4);
    ASSERT_EQ(edgeIntersections[0].edgeSecondNode, 5);
    ASSERT_NEAR(edgeIntersections[0].edgeDistance, 0.50000000000000000, 1e-8);

    ASSERT_EQ(edgeIntersections[1].polylineSegmentIndex, 0);
    ASSERT_NEAR(edgeIntersections[1].adimensionalPolylineSegmentDistance, 0.75000000000000000, 1e-8);
    ASSERT_EQ(edgeIntersections[1].edgeIndex, 18);
    ASSERT_EQ(edgeIntersections[1].edgeFirstNode, 8);
    ASSERT_EQ(edgeIntersections[1].edgeSecondNode, 9);
    ASSERT_NEAR(edgeIntersections[1].edgeDistance, 0.50000000000000000, 1e-8);

    ASSERT_EQ(edgeIntersections[2].polylineSegmentIndex, 1);
    ASSERT_NEAR(edgeIntersections[2].adimensionalPolylineSegmentDistance, 0.25000000000000000, 1e-8);
    ASSERT_EQ(edgeIntersections[2].edgeIndex, 9);
    ASSERT_EQ(edgeIntersections[2].edgeFirstNode, 13);
    ASSERT_EQ(edgeIntersections[2].edgeSecondNode, 9);
    ASSERT_NEAR(edgeIntersections[2].edgeDistance, 0.50000000000000000, 1e-8);

    ASSERT_EQ(edgeIntersections[3].polylineSegmentIndex, 1);
    ASSERT_NEAR(edgeIntersections[3].adimensionalPolylineSegmentDistance, 0.75000000000000000, 1e-8);
    ASSERT_EQ(edgeIntersections[3].edgeIndex, 10);
    ASSERT_EQ(edgeIntersections[3].edgeFirstNode, 14);
    ASSERT_EQ(edgeIntersections[3].edgeSecondNode, 10);
    ASSERT_NEAR(edgeIntersections[3].edgeDistance, 0.50000000000000000, 1e-8);

    ASSERT_EQ(edgeIntersections[4].polylineSegmentIndex, 2);
    ASSERT_NEAR(edgeIntersections[4].adimensionalPolylineSegmentDistance, 0.25000000000000000, 1e-8);
    ASSERT_EQ(edgeIntersections[4].edgeIndex, 20);
    ASSERT_EQ(edgeIntersections[4].edgeFirstNode, 11);
    ASSERT_EQ(edgeIntersections[4].edgeSecondNode, 10);
    ASSERT_NEAR(edgeIntersections[4].edgeDistance, 0.50000000000000000, 1e-8);

    ASSERT_EQ(edgeIntersections[5].polylineSegmentIndex, 2);
    ASSERT_NEAR(edgeIntersections[5].adimensionalPolylineSegmentDistance, 0.75000000000000000, 1e-8);
    ASSERT_EQ(edgeIntersections[5].edgeIndex, 17);
    ASSERT_EQ(edgeIntersections[5].edgeFirstNode, 7);
    ASSERT_EQ(edgeIntersections[5].edgeSecondNode, 6);
    ASSERT_NEAR(edgeIntersections[5].edgeDistance, 0.50000000000000000, 1e-8);

    ASSERT_EQ(edgeIntersections[6].polylineSegmentIndex, 3);
    ASSERT_NEAR(edgeIntersections[6].adimensionalPolylineSegmentDistance, 0.25000000000000000, 1e-8);
    ASSERT_EQ(edgeIntersections[6].edgeIndex, 2);
    ASSERT_EQ(edgeIntersections[6].edgeFirstNode, 2);
    ASSERT_EQ(edgeIntersections[6].edgeSecondNode, 6);
    ASSERT_NEAR(edgeIntersections[6].edgeDistance, 0.50000000000000000, 1e-8);

    ASSERT_EQ(edgeIntersections[7].polylineSegmentIndex, 3);
    ASSERT_NEAR(edgeIntersections[7].adimensionalPolylineSegmentDistance, 0.75000000000000000, 1e-8);
    ASSERT_EQ(edgeIntersections[7].edgeIndex, 1);
    ASSERT_EQ(edgeIntersections[7].edgeFirstNode, 1);
    ASSERT_EQ(edgeIntersections[7].edgeSecondNode, 5);
    ASSERT_NEAR(edgeIntersections[7].edgeDistance, 0.50000000000000000, 1e-8);

    // face intersections
    ASSERT_EQ(faceIntersections[0].faceIndex, 3);
    ASSERT_NEAR(faceIntersections[0].polylineDistance, 1.0, 1e-8);
    ASSERT_EQ(faceIntersections[0].edgeNodes[0], 4);
    ASSERT_EQ(faceIntersections[0].edgeNodes[1], 5);
    ASSERT_EQ(faceIntersections[0].edgeNodes[2], 8);
    ASSERT_EQ(faceIntersections[0].edgeNodes[3], 9);
    ASSERT_EQ(faceIntersections[0].edgeIndices.size(), 2);

    ASSERT_EQ(faceIntersections[1].faceIndex, 6);
    ASSERT_NEAR(faceIntersections[1].polylineDistance, 2.0, 1e-8);
    ASSERT_EQ(faceIntersections[1].edgeNodes[0], 13);
    ASSERT_EQ(faceIntersections[1].edgeNodes[1], 9);
    ASSERT_EQ(faceIntersections[1].edgeNodes[2], 8);
    ASSERT_EQ(faceIntersections[1].edgeNodes[3], 9);
    ASSERT_EQ(faceIntersections[1].edgeIndices.size(), 2);

    ASSERT_EQ(faceIntersections[2].faceIndex, 7);
    ASSERT_NEAR(faceIntersections[2].polylineDistance, 3.0, 1e-8);
    ASSERT_EQ(faceIntersections[2].edgeNodes[0], 13);
    ASSERT_EQ(faceIntersections[2].edgeNodes[1], 9);
    ASSERT_EQ(faceIntersections[2].edgeNodes[2], 14);
    ASSERT_EQ(faceIntersections[2].edgeNodes[3], 10);
    ASSERT_EQ(faceIntersections[2].edgeIndices.size(), 2);

    ASSERT_EQ(faceIntersections[3].faceIndex, 0);
    ASSERT_NEAR(faceIntersections[3].polylineDistance, 4.0, 1e-8);
    ASSERT_EQ(faceIntersections[3].edgeNodes[0], 4);
    ASSERT_EQ(faceIntersections[3].edgeNodes[1], 5);
    ASSERT_EQ(faceIntersections[3].edgeNodes[2], 1);
    ASSERT_EQ(faceIntersections[3].edgeNodes[3], 5);
    ASSERT_EQ(faceIntersections[3].edgeIndices.size(), 2);

    ASSERT_EQ(faceIntersections[4].faceIndex, 8);
    ASSERT_NEAR(faceIntersections[4].polylineDistance, 4.0, 1e-8);
    ASSERT_EQ(faceIntersections[4].edgeNodes[0], 11);
    ASSERT_EQ(faceIntersections[4].edgeNodes[1], 10);
    ASSERT_EQ(faceIntersections[4].edgeNodes[2], 14);
    ASSERT_EQ(faceIntersections[4].edgeNodes[3], 10);
    ASSERT_EQ(faceIntersections[4].edgeIndices.size(), 2);

    ASSERT_EQ(faceIntersections[5].faceIndex, 5);
    ASSERT_NEAR(faceIntersections[5].polylineDistance, 5.0, 1e-8);
    ASSERT_EQ(faceIntersections[5].edgeNodes[0], 11);
    ASSERT_EQ(faceIntersections[5].edgeNodes[1], 10);
    ASSERT_EQ(faceIntersections[5].edgeNodes[2], 7);
    ASSERT_EQ(faceIntersections[5].edgeNodes[3], 6);
    ASSERT_EQ(faceIntersections[5].edgeIndices.size(), 2);

    ASSERT_EQ(faceIntersections[6].faceIndex, 2);
    ASSERT_NEAR(faceIntersections[6].polylineDistance, 6.0, 1e-8);
    ASSERT_EQ(faceIntersections[6].edgeNodes[0], 2);
    ASSERT_EQ(faceIntersections[6].edgeNodes[1], 6);
    ASSERT_EQ(faceIntersections[6].edgeNodes[2], 7);
    ASSERT_EQ(faceIntersections[6].edgeNodes[3], 6);
    ASSERT_EQ(faceIntersections[6].edgeIndices.size(), 2);

    ASSERT_EQ(faceIntersections[7].faceIndex, 1);
    ASSERT_NEAR(faceIntersections[7].polylineDistance, 7.0, 1e-8);
    ASSERT_EQ(faceIntersections[7].edgeNodes[0], 2);
    ASSERT_EQ(faceIntersections[7].edgeNodes[1], 6);
    ASSERT_EQ(faceIntersections[7].edgeNodes[2], 1);
    ASSERT_EQ(faceIntersections[7].edgeNodes[3], 5);
    ASSERT_EQ(faceIntersections[7].edgeIndices.size(), 2);
}

TEST(Mesh2D, GetPolylineIntersectionsFromObliqueLineShouldReturnCorrectIntersections)
{
    // 1. Setup
    auto mesh = MakeRectangularMeshForTesting(6, 6, 1.0, meshkernel::Projection::cartesian);

    std::vector<meshkernel::Point> polyLine;
    polyLine.emplace_back(3.9, 0.0);
    polyLine.emplace_back(0.0, 3.9);

    // 2. Execute
    meshkernel::Mesh2DIntersections mesh2DIntersections(*mesh);
    mesh2DIntersections.Compute(polyLine);
    auto edgeIntersections = mesh2DIntersections.EdgeIntersections();
    auto faceIntersections = mesh2DIntersections.FaceIntersections();

    meshkernel::Mesh2DIntersections::sortAndEraseIntersections(edgeIntersections);
    meshkernel::Mesh2DIntersections::sortAndEraseIntersections(faceIntersections);

    // 3. Assert

    // edge intersection
    ASSERT_EQ(edgeIntersections[0].polylineSegmentIndex, 0);
    ASSERT_NEAR(edgeIntersections[0].adimensionalPolylineSegmentDistance, 0.0, 1e-8);
    ASSERT_EQ(edgeIntersections[0].edgeIndex, 18);
    ASSERT_EQ(edgeIntersections[0].edgeFirstNode, 24);
    ASSERT_EQ(edgeIntersections[0].edgeSecondNode, 18);
    ASSERT_NEAR(edgeIntersections[0].edgeDistance, 0.89999999999999991, 1e-8);

    ASSERT_EQ(edgeIntersections[1].polylineSegmentIndex, 0);
    ASSERT_NEAR(edgeIntersections[1].adimensionalPolylineSegmentDistance, 0.23076923076923075, 1e-8);
    ASSERT_EQ(edgeIntersections[1].edgeIndex, 45);
    ASSERT_EQ(edgeIntersections[1].edgeFirstNode, 19);
    ASSERT_EQ(edgeIntersections[1].edgeSecondNode, 18);
    ASSERT_NEAR(edgeIntersections[1].edgeDistance, 0.10000000000000003, 1e-8);

    ASSERT_EQ(edgeIntersections[2].polylineSegmentIndex, 0);
    ASSERT_NEAR(edgeIntersections[2].adimensionalPolylineSegmentDistance, 0.25641025641025644, 1e-8);
    ASSERT_EQ(edgeIntersections[2].edgeIndex, 13);
    ASSERT_EQ(edgeIntersections[2].edgeFirstNode, 19);
    ASSERT_EQ(edgeIntersections[2].edgeSecondNode, 13);
    ASSERT_NEAR(edgeIntersections[2].edgeDistance, 0.89999999999999980, 1e-8);

    ASSERT_EQ(edgeIntersections[3].polylineSegmentIndex, 0);
    ASSERT_NEAR(edgeIntersections[3].adimensionalPolylineSegmentDistance, 0.48717948717948717, 1e-8);
    ASSERT_EQ(edgeIntersections[3].edgeIndex, 41);
    ASSERT_EQ(edgeIntersections[3].edgeFirstNode, 14);
    ASSERT_EQ(edgeIntersections[3].edgeSecondNode, 13);
    ASSERT_NEAR(edgeIntersections[3].edgeDistance, 0.10000000000000014, 1e-8);

    ASSERT_EQ(edgeIntersections[4].polylineSegmentIndex, 0);
    ASSERT_NEAR(edgeIntersections[4].adimensionalPolylineSegmentDistance, 0.51282051282051289, 1e-8);
    ASSERT_EQ(edgeIntersections[4].edgeIndex, 8);
    ASSERT_EQ(edgeIntersections[4].edgeFirstNode, 14);
    ASSERT_EQ(edgeIntersections[4].edgeSecondNode, 8);
    ASSERT_NEAR(edgeIntersections[4].edgeDistance, 0.89999999999999969, 1e-8);

    ASSERT_EQ(edgeIntersections[5].polylineSegmentIndex, 0);
    ASSERT_NEAR(edgeIntersections[5].adimensionalPolylineSegmentDistance, 0.74358974358974361, 1e-8);
    ASSERT_EQ(edgeIntersections[5].edgeIndex, 37);
    ASSERT_EQ(edgeIntersections[5].edgeFirstNode, 9);
    ASSERT_EQ(edgeIntersections[5].edgeSecondNode, 8);
    ASSERT_NEAR(edgeIntersections[5].edgeDistance, 0.10000000000000014, 1e-8);

    ASSERT_EQ(edgeIntersections[6].polylineSegmentIndex, 0);
    ASSERT_NEAR(edgeIntersections[6].adimensionalPolylineSegmentDistance, 0.76923076923076927, 1e-8);
    ASSERT_EQ(edgeIntersections[6].edgeIndex, 3);
    ASSERT_EQ(edgeIntersections[6].edgeFirstNode, 9);
    ASSERT_EQ(edgeIntersections[6].edgeSecondNode, 3);
    ASSERT_NEAR(edgeIntersections[6].edgeDistance, 0.89999999999999991, 1e-8);

    ASSERT_EQ(edgeIntersections[7].polylineSegmentIndex, 0);
    ASSERT_NEAR(edgeIntersections[7].adimensionalPolylineSegmentDistance, 1.0, 1e-8);
    ASSERT_EQ(edgeIntersections[7].edgeIndex, 33);
    ASSERT_EQ(edgeIntersections[7].edgeFirstNode, 4);
    ASSERT_EQ(edgeIntersections[7].edgeSecondNode, 3);
    ASSERT_NEAR(edgeIntersections[7].edgeDistance, 0.10000000000000014, 1e-8);

    // face intersections
    ASSERT_EQ(faceIntersections[0].faceIndex, 15);
    ASSERT_NEAR(faceIntersections[0].polylineDistance, 0.63639610306789274, 1e-8);
    ASSERT_EQ(faceIntersections[0].edgeNodes[0], 24);
    ASSERT_EQ(faceIntersections[0].edgeNodes[1], 18);
    ASSERT_EQ(faceIntersections[0].edgeNodes[2], 19);
    ASSERT_EQ(faceIntersections[0].edgeNodes[3], 18);
    ASSERT_EQ(faceIntersections[0].edgeIndices.size(), 2);

    ASSERT_EQ(faceIntersections[1].faceIndex, 10);
    ASSERT_NEAR(faceIntersections[1].polylineDistance, 1.3435028842544403, 1e-8);
    ASSERT_EQ(faceIntersections[1].edgeNodes[0], 19);
    ASSERT_EQ(faceIntersections[1].edgeNodes[1], 18);
    ASSERT_EQ(faceIntersections[1].edgeNodes[2], 19);
    ASSERT_EQ(faceIntersections[1].edgeNodes[3], 13);
    ASSERT_EQ(faceIntersections[1].edgeIndices.size(), 2);

    ASSERT_EQ(faceIntersections[2].faceIndex, 11);
    ASSERT_NEAR(faceIntersections[2].polylineDistance, 2.0506096654409878, 1e-8);
    ASSERT_EQ(faceIntersections[2].edgeNodes[0], 19);
    ASSERT_EQ(faceIntersections[2].edgeNodes[1], 13);
    ASSERT_EQ(faceIntersections[2].edgeNodes[2], 14);
    ASSERT_EQ(faceIntersections[2].edgeNodes[3], 13);
    ASSERT_EQ(faceIntersections[2].edgeIndices.size(), 2);

    ASSERT_EQ(faceIntersections[3].faceIndex, 6);
    ASSERT_NEAR(faceIntersections[3].polylineDistance, 2.7577164466275352, 1e-8);
    ASSERT_EQ(faceIntersections[3].edgeNodes[0], 14);
    ASSERT_EQ(faceIntersections[3].edgeNodes[1], 13);
    ASSERT_EQ(faceIntersections[3].edgeNodes[2], 14);
    ASSERT_EQ(faceIntersections[3].edgeNodes[3], 8);
    ASSERT_EQ(faceIntersections[3].edgeIndices.size(), 2);

    ASSERT_EQ(faceIntersections[4].faceIndex, 7);
    ASSERT_NEAR(faceIntersections[4].polylineDistance, 3.4648232278140831, 1e-8);
    ASSERT_EQ(faceIntersections[4].edgeNodes[0], 14);
    ASSERT_EQ(faceIntersections[4].edgeNodes[1], 8);
    ASSERT_EQ(faceIntersections[4].edgeNodes[2], 9);
    ASSERT_EQ(faceIntersections[4].edgeNodes[3], 8);
    ASSERT_EQ(faceIntersections[4].edgeIndices.size(), 2);

    ASSERT_EQ(faceIntersections[5].faceIndex, 2);
    ASSERT_NEAR(faceIntersections[5].polylineDistance, 4.1719300090006302, 1e-8);
    ASSERT_EQ(faceIntersections[5].edgeNodes[0], 9);
    ASSERT_EQ(faceIntersections[5].edgeNodes[1], 8);
    ASSERT_EQ(faceIntersections[5].edgeNodes[2], 9);
    ASSERT_EQ(faceIntersections[5].edgeNodes[3], 3);
    ASSERT_EQ(faceIntersections[5].edgeIndices.size(), 2);

    ASSERT_EQ(faceIntersections[6].faceIndex, 3);
    ASSERT_NEAR(faceIntersections[6].polylineDistance, 4.8790367901871772, 1e-8);
    ASSERT_EQ(faceIntersections[6].edgeNodes[0], 9);
    ASSERT_EQ(faceIntersections[6].edgeNodes[1], 3);
    ASSERT_EQ(faceIntersections[6].edgeNodes[2], 4);
    ASSERT_EQ(faceIntersections[6].edgeNodes[3], 3);
    ASSERT_EQ(faceIntersections[6].edgeIndices.size(), 2);
}

TEST(Mesh2D, GetPolylineIntersectionsFromComplexPolylineShouldReturnCorrectIntersections)
{
    // 1. Setup
    const meshkernel::Point origin{78.0, 45.0};
    auto mesh = MakeRectangularMeshForTesting(7, 7, 1.0, meshkernel::Projection::cartesian, origin);

    std::vector<meshkernel::Point> boundaryPolyline;
    boundaryPolyline.emplace_back(80.6623, 50.0074);
    boundaryPolyline.emplace_back(81.4075, 49.3843);
    boundaryPolyline.emplace_back(81.845, 48.885);
    boundaryPolyline.emplace_back(82.1464, 48.3577);
    boundaryPolyline.emplace_back(82.3599, 47.7658);
    boundaryPolyline.emplace_back(82.4847, 47.1451);
    boundaryPolyline.emplace_back(82.5261, 46.556);
    boundaryPolyline.emplace_back(82.5038, 46.0853);
    boundaryPolyline.emplace_back(82.0738, 45.8102);
    boundaryPolyline.emplace_back(81.0887, 45.2473);

    // 2. Execute
    meshkernel::Mesh2DIntersections mesh2DIntersections(*mesh);
    mesh2DIntersections.Compute(boundaryPolyline);
    auto edgeIntersections = mesh2DIntersections.EdgeIntersections();
    auto faceIntersections = mesh2DIntersections.FaceIntersections();

    meshkernel::Mesh2DIntersections::sortAndEraseIntersections(edgeIntersections);
    meshkernel::Mesh2DIntersections::sortAndEraseIntersections(faceIntersections);

    // 3. Assert

    // edge intersection
    ASSERT_EQ(edgeIntersections[0].polylineSegmentIndex, 0);
    ASSERT_NEAR(edgeIntersections[0].adimensionalPolylineSegmentDistance, 0.011876103354192005, 1e-8);
    ASSERT_EQ(edgeIntersections[0].edgeIndex, 19);
    ASSERT_EQ(edgeIntersections[0].edgeFirstNode, 19);
    ASSERT_EQ(edgeIntersections[0].edgeSecondNode, 26);
    ASSERT_NEAR(edgeIntersections[0].edgeDistance, 0.67115007221954570, 1e-8);

    ASSERT_EQ(edgeIntersections[1].polylineSegmentIndex, 0);
    ASSERT_NEAR(edgeIntersections[1].adimensionalPolylineSegmentDistance, 0.45316693505099231, 1e-8);
    ASSERT_EQ(edgeIntersections[1].edgeIndex, 64);
    ASSERT_EQ(edgeIntersections[1].edgeFirstNode, 25);
    ASSERT_EQ(edgeIntersections[1].edgeSecondNode, 26);
    ASSERT_NEAR(edgeIntersections[1].edgeDistance, 0.27496831723027354, 1e-8);

    ASSERT_EQ(edgeIntersections[2].polylineSegmentIndex, 1);
    ASSERT_NEAR(edgeIntersections[2].adimensionalPolylineSegmentDistance, 0.76967754856799364, 1e-8);
    ASSERT_EQ(edgeIntersections[2].edgeIndex, 25);
    ASSERT_EQ(edgeIntersections[2].edgeFirstNode, 25);
    ASSERT_EQ(edgeIntersections[2].edgeSecondNode, 32);
    ASSERT_NEAR(edgeIntersections[2].edgeDistance, 0.74423392749849604, 1e-8);

    ASSERT_EQ(edgeIntersections[3].polylineSegmentIndex, 2);
    ASSERT_NEAR(edgeIntersections[3].adimensionalPolylineSegmentDistance, 0.51426675514266962, 1e-8);
    ASSERT_EQ(edgeIntersections[3].edgeIndex, 69);
    ASSERT_EQ(edgeIntersections[3].edgeFirstNode, 31);
    ASSERT_EQ(edgeIntersections[3].edgeSecondNode, 32);
    ASSERT_NEAR(edgeIntersections[3].edgeDistance, 0.38617285998673001, 1e-8);

    ASSERT_EQ(edgeIntersections[4].polylineSegmentIndex, 3);
    ASSERT_NEAR(edgeIntersections[4].adimensionalPolylineSegmentDistance, 0.60432505490792310, 1e-8);
    ASSERT_EQ(edgeIntersections[4].edgeIndex, 31);
    ASSERT_EQ(edgeIntersections[4].edgeFirstNode, 31);
    ASSERT_EQ(edgeIntersections[4].edgeSecondNode, 38);
    ASSERT_NEAR(edgeIntersections[4].edgeDistance, 0.27542339922283915, 1e-8);

    ASSERT_EQ(edgeIntersections[5].polylineSegmentIndex, 5);
    ASSERT_NEAR(edgeIntersections[5].adimensionalPolylineSegmentDistance, 0.24630792734679827, 1e-8);
    ASSERT_EQ(edgeIntersections[5].edgeIndex, 30);
    ASSERT_EQ(edgeIntersections[5].edgeFirstNode, 30);
    ASSERT_EQ(edgeIntersections[5].edgeSecondNode, 37);
    ASSERT_NEAR(edgeIntersections[5].edgeDistance, 0.49489714819216007, 1e-8);

    ASSERT_EQ(edgeIntersections[6].polylineSegmentIndex, 7);
    ASSERT_NEAR(edgeIntersections[6].adimensionalPolylineSegmentDistance, 0.31006906579425014, 1e-8);
    ASSERT_EQ(edgeIntersections[6].edgeIndex, 29);
    ASSERT_EQ(edgeIntersections[6].edgeFirstNode, 29);
    ASSERT_EQ(edgeIntersections[6].edgeSecondNode, 36);
    ASSERT_NEAR(edgeIntersections[6].edgeDistance, 0.37047030170847295, 1e-8);

    ASSERT_EQ(edgeIntersections[7].polylineSegmentIndex, 8);
    ASSERT_NEAR(edgeIntersections[7].adimensionalPolylineSegmentDistance, 0.074916252157146923, 1e-8);
    ASSERT_EQ(edgeIntersections[7].edgeIndex, 66);
    ASSERT_EQ(edgeIntersections[7].edgeFirstNode, 29);
    ASSERT_EQ(edgeIntersections[7].edgeSecondNode, 28);
    ASSERT_NEAR(edgeIntersections[7].edgeDistance, 0.23197035833925611, 1e-8);

    // face intersections
    ASSERT_EQ(faceIntersections[0].faceIndex, 17);
    ASSERT_NEAR(faceIntersections[0].polylineDistance, 0.011536194272423500, 1e-8);
    ASSERT_EQ(faceIntersections[0].edgeNodes[0], 19);
    ASSERT_EQ(faceIntersections[0].edgeNodes[1], 26);
    ASSERT_EQ(faceIntersections[0].edgeIndices.size(), 1);

    ASSERT_EQ(faceIntersections[1].faceIndex, 16);
    ASSERT_NEAR(faceIntersections[1].polylineDistance, 0.22586645956506612, 1e-8);
    ASSERT_EQ(faceIntersections[1].edgeNodes[0], 19);
    ASSERT_EQ(faceIntersections[1].edgeNodes[1], 26);
    ASSERT_EQ(faceIntersections[1].edgeNodes[2], 25);
    ASSERT_EQ(faceIntersections[1].edgeNodes[3], 26);
    ASSERT_EQ(faceIntersections[1].edgeIndices.size(), 2);

    ASSERT_EQ(faceIntersections[2].faceIndex, 22);
    ASSERT_NEAR(faceIntersections[2].polylineDistance, 0.96126582566625318, 1e-8);
    ASSERT_EQ(faceIntersections[2].edgeNodes[0], 25);
    ASSERT_EQ(faceIntersections[2].edgeNodes[1], 26);
    ASSERT_EQ(faceIntersections[2].edgeNodes[2], 25);
    ASSERT_EQ(faceIntersections[2].edgeNodes[3], 32);
    ASSERT_EQ(faceIntersections[2].edgeIndices.size(), 2);

    ASSERT_EQ(faceIntersections[3].faceIndex, 21);
    ASSERT_NEAR(faceIntersections[3].polylineDistance, 1.7149583232814580, 1e-8);
    ASSERT_EQ(faceIntersections[3].edgeNodes[0], 31);
    ASSERT_EQ(faceIntersections[3].edgeNodes[1], 32);
    ASSERT_EQ(faceIntersections[3].edgeNodes[2], 25);
    ASSERT_EQ(faceIntersections[3].edgeNodes[3], 32);
    ASSERT_EQ(faceIntersections[3].edgeIndices.size(), 2);

    ASSERT_EQ(faceIntersections[4].faceIndex, 27);
    ASSERT_NEAR(faceIntersections[4].polylineDistance, 2.2852185268843637, 1e-8);
    ASSERT_EQ(faceIntersections[4].edgeNodes[0], 31);
    ASSERT_EQ(faceIntersections[4].edgeNodes[1], 32);
    ASSERT_EQ(faceIntersections[4].edgeNodes[2], 31);
    ASSERT_EQ(faceIntersections[4].edgeNodes[3], 38);
    ASSERT_EQ(faceIntersections[4].edgeIndices.size(), 2);

    ASSERT_EQ(faceIntersections[5].faceIndex, 26);
    ASSERT_NEAR(faceIntersections[5].polylineDistance, 3.1366301680701545, 1e-8);
    ASSERT_EQ(faceIntersections[5].edgeNodes[0], 30);
    ASSERT_EQ(faceIntersections[5].edgeNodes[1], 37);
    ASSERT_EQ(faceIntersections[5].edgeNodes[2], 31);
    ASSERT_EQ(faceIntersections[5].edgeNodes[3], 38);
    ASSERT_EQ(faceIntersections[5].edgeIndices.size(), 2);

    ASSERT_EQ(faceIntersections[6].faceIndex, 25);
    ASSERT_NEAR(faceIntersections[6].polylineDistance, 4.1877070472004538, 1e-8);
    ASSERT_EQ(faceIntersections[6].edgeNodes[0], 30);
    ASSERT_EQ(faceIntersections[6].edgeNodes[1], 37);
    ASSERT_EQ(faceIntersections[6].edgeNodes[2], 29);
    ASSERT_EQ(faceIntersections[6].edgeNodes[3], 36);
    ASSERT_EQ(faceIntersections[6].edgeIndices.size(), 2);

    ASSERT_EQ(faceIntersections[7].faceIndex, 24);
    ASSERT_NEAR(faceIntersections[7].polylineDistance, 4.9436030957613966, 1e-8);
    ASSERT_EQ(faceIntersections[7].edgeNodes[0], 29);
    ASSERT_EQ(faceIntersections[7].edgeNodes[1], 28);
    ASSERT_EQ(faceIntersections[7].edgeNodes[2], 29);
    ASSERT_EQ(faceIntersections[7].edgeNodes[3], 36);
    ASSERT_EQ(faceIntersections[7].edgeIndices.size(), 2);

    ASSERT_EQ(faceIntersections[8].faceIndex, 18);
    ASSERT_NEAR(faceIntersections[8].polylineDistance, 5.1621970995815847, 1e-8);
    ASSERT_EQ(faceIntersections[8].edgeNodes[0], 29);
    ASSERT_EQ(faceIntersections[8].edgeNodes[1], 28);
    ASSERT_EQ(faceIntersections[8].edgeIndices.size(), 1);
}

TEST(Mesh2D, RemoveSingleIsland)
{
    // Load mesh with 2 disconnected regions, first a 10x10 and the second is a smaller 2x2 mesh
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/RemoveDomainIslands/single_disconnected_region.nc");
    meshkernel::RemoveDisconnectedRegions removeDisconnectedRegions;

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    // Remove all smaller disconnected "island" regions.
    auto undoAction = removeDisconnectedRegions.Compute(*mesh);
    EXPECT_EQ(mesh->GetNumFaces(), 100);

    // Restore original mesh
    undoAction->Restore();
    mesh->Administrate();

    EXPECT_EQ(mesh->GetNumFaces(), 104);
    ASSERT_EQ(originalNodes.size(), mesh->Nodes().size());
    ASSERT_EQ(originalEdges.size(), mesh->Edges().size());

    for (meshkernel::UInt i = 0; i < mesh->Nodes().size(); ++i)
    {
        EXPECT_EQ(originalNodes[i].x, mesh->Node(i).x);
        EXPECT_EQ(originalNodes[i].y, mesh->Node(i).y);
    }

    for (meshkernel::UInt i = 0; i < mesh->Edges().size(); ++i)
    {
        EXPECT_EQ(originalEdges[i].first, mesh->GetEdge(i).first);
        EXPECT_EQ(originalEdges[i].second, mesh->GetEdge(i).second);
    }
}

TEST(Mesh2D, RemoveMultipleIslands)
{
    // Load mesh with 4 disconnected regions, the main domain is a 10x10, there are 3 other much small island regions,
    // each with a different shape and number of elements.
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/RemoveDomainIslands/multiple_disconnected_regions.nc");

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    meshkernel::RemoveDisconnectedRegions removeDisconnectedRegions;

    // Remove all smaller disconnected "island" regions.
    auto undoAction = removeDisconnectedRegions.Compute(*mesh);
    EXPECT_EQ(mesh->GetNumFaces(), 100);

    // Restore original mesh
    undoAction->Restore();
    mesh->Administrate();

    EXPECT_EQ(mesh->GetNumFaces(), 113);
    ASSERT_EQ(originalNodes.size(), mesh->Nodes().size());
    ASSERT_EQ(originalEdges.size(), mesh->Edges().size());

    for (meshkernel::UInt i = 0; i < mesh->Nodes().size(); ++i)
    {
        EXPECT_EQ(originalNodes[i].x, mesh->Node(i).x);
        EXPECT_EQ(originalNodes[i].y, mesh->Node(i).y);
    }

    for (meshkernel::UInt i = 0; i < mesh->Edges().size(); ++i)
    {
        EXPECT_EQ(originalEdges[i].first, mesh->GetEdge(i).first);
        EXPECT_EQ(originalEdges[i].second, mesh->GetEdge(i).second);
    }
}

TEST(Mesh2D, DeleteMesh_WhenFacesAreIntersected_ShouldNotDeleteFaces)
{
    // Prepare
    const auto mesh = MakeRectangularMeshForTesting(4, 4, 3, 3, meshkernel::Projection::cartesian, meshkernel::Point{0, 0});

    // a polygon including all nodes of a face, but also intersecting
    std::vector<meshkernel::Point> polygonNodes{
        {1.87622950819672, -0.299180327868853},
        {1.86885245901639, 0.187704918032786},
        {3.27049180327869, 0.195081967213114},
        {3.27049180327869, 0.320491803278688},
        {1.87622950819672, 0.320491803278688},
        {1.86147540983607, 1.16147540983606},
        {3.54344262295082, 1.18360655737705},
        {3.55081967213115, -0.358196721311476},
        {1.87622950819672, -0.299180327868853}};

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    auto polygon = meshkernel::Polygons(polygonNodes, meshkernel::Projection::cartesian);
    const auto deletion_option = meshkernel::Mesh2D::DeleteMeshOptions::InsideNotIntersected;

    // Execute
    auto undoAction = mesh->DeleteMesh(polygon, deletion_option, false);

    // Assert
    EXPECT_EQ(mesh->GetNumFaces(), 9);

    // Restore original mesh
    undoAction->Restore();

    ASSERT_EQ(originalNodes.size(), mesh->Nodes().size());
    ASSERT_EQ(originalEdges.size(), mesh->Edges().size());

    for (meshkernel::UInt i = 0; i < mesh->Nodes().size(); ++i)
    {
        EXPECT_EQ(originalNodes[i].x, mesh->Node(i).x);
        EXPECT_EQ(originalNodes[i].y, mesh->Node(i).y);
    }

    for (meshkernel::UInt i = 0; i < mesh->Edges().size(); ++i)
    {
        EXPECT_EQ(originalEdges[i].first, mesh->GetEdge(i).first);
        EXPECT_EQ(originalEdges[i].second, mesh->GetEdge(i).second);
    }
}

TEST(Mesh2D, DeleteMesh_WhenFacesAreIntersectedSpherical_ShouldNotDeleteFaces)
{
    // Prepare
    const auto mesh = MakeRectangularMeshForTesting(4, 4, 3, 3, meshkernel::Projection::spherical, meshkernel::Point{0, 0});

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    // a polygon including all nodes of a face, but also intersecting one
    std::vector<meshkernel::Point> polygonNodes{
        {1.87622950819672, -0.299180327868853},
        {1.86885245901639, 0.187704918032786},
        {3.27049180327869, 0.195081967213114},
        {3.27049180327869, 0.320491803278688},
        {1.87622950819672, 0.320491803278688},
        {1.86147540983607, 1.16147540983606},
        {3.54344262295082, 1.18360655737705},
        {3.55081967213115, -0.358196721311476},
        {1.87622950819672, -0.299180327868853}};

    auto polygon = meshkernel::Polygons(polygonNodes, meshkernel::Projection::spherical);
    const auto deletion_option = meshkernel::Mesh2D::DeleteMeshOptions::InsideNotIntersected;

    // Execute
    auto undoAction = mesh->DeleteMesh(polygon, deletion_option, false);

    // Assert
    EXPECT_EQ(mesh->GetNumFaces(), 9);

    // Restore original mesh
    undoAction->Restore();

    ASSERT_EQ(originalNodes.size(), mesh->Nodes().size());
    ASSERT_EQ(originalEdges.size(), mesh->Edges().size());

    for (meshkernel::UInt i = 0; i < mesh->Nodes().size(); ++i)
    {
        EXPECT_EQ(originalNodes[i].x, mesh->Node(i).x);
        EXPECT_EQ(originalNodes[i].y, mesh->Node(i).y);
    }

    for (meshkernel::UInt i = 0; i < mesh->Edges().size(); ++i)
    {
        EXPECT_EQ(originalEdges[i].first, mesh->GetEdge(i).first);
        EXPECT_EQ(originalEdges[i].second, mesh->GetEdge(i).second);
    }
}

TEST(Mesh2D, DeleteMesh_WithLargeSphericalPolygon_ShouldDeleteInnerMeshFaces)
{
    // Prepare
    const auto mesh = MakeRectangularMeshForTesting(4,
                                                    4,
                                                    2.0,
                                                    2.0,
                                                    meshkernel::Projection::spherical,
                                                    meshkernel::Point{-3.0, 48.5});

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    // a large polygon
    std::vector<meshkernel::Point> polygonNodes{
        {-2.29490103397341, 50.0126381093058},
        {179.33620776839, 50.3853885542098},
        {180.05965832319, -3.87340305583453},
        {-2.24988148655834, -3.14995250103394},
        {-2.29490103397341, 50.0126381093058}};

    auto polygon = meshkernel::Polygons(polygonNodes, meshkernel::Projection::spherical);
    const auto deletion_option = meshkernel::Mesh2D::DeleteMeshOptions::InsideNotIntersected;

    // Execute
    auto undoAction = mesh->DeleteMesh(polygon, deletion_option, false);

    // Assert
    EXPECT_EQ(mesh->GetNumFaces(), 7);

    // Restore original mesh
    undoAction->Restore();

    ASSERT_EQ(originalNodes.size(), mesh->Nodes().size());
    ASSERT_EQ(originalEdges.size(), mesh->Edges().size());

    for (meshkernel::UInt i = 0; i < mesh->Nodes().size(); ++i)
    {
        EXPECT_EQ(originalNodes[i].x, mesh->Node(i).x);
        EXPECT_EQ(originalNodes[i].y, mesh->Node(i).y);
    }

    for (meshkernel::UInt i = 0; i < mesh->Edges().size(); ++i)
    {
        EXPECT_EQ(originalEdges[i].first, mesh->GetEdge(i).first);
        EXPECT_EQ(originalEdges[i].second, mesh->GetEdge(i).second);
    }
}

TEST(Mesh2D, DeleteMesh_WithPolygonAndIncludedCircumcenters_ShouldDeleteInnerFaces)
{
    // Prepare
    const auto mesh = MakeRectangularMeshForTesting(5,
                                                    5,
                                                    8.0,
                                                    8.0,
                                                    meshkernel::Projection::cartesian);

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    // a large polygon
    std::vector<meshkernel::Point> polygonNodes{
        {2, 2},
        {6, 2},
        {6, 6},
        {2, 6},
        {2, 2}};

    auto polygon = meshkernel::Polygons(polygonNodes, meshkernel::Projection::cartesian);
    const auto deletion_option = meshkernel::Mesh2D::DeleteMeshOptions::FacesWithIncludedCircumcenters;

    // Execute
    auto undoAction = mesh->DeleteMesh(polygon, deletion_option, false);

    // Assert
    EXPECT_EQ(mesh->GetNumFaces(), 12);

    // Restore original mesh
    undoAction->Restore();

    ASSERT_EQ(originalNodes.size(), mesh->Nodes().size());
    ASSERT_EQ(originalEdges.size(), mesh->Edges().size());

    for (meshkernel::UInt i = 0; i < mesh->Nodes().size(); ++i)
    {
        EXPECT_EQ(originalNodes[i].x, mesh->Node(i).x);
        EXPECT_EQ(originalNodes[i].y, mesh->Node(i).y);
    }

    for (meshkernel::UInt i = 0; i < mesh->Edges().size(); ++i)
    {
        EXPECT_EQ(originalEdges[i].first, mesh->GetEdge(i).first);
        EXPECT_EQ(originalEdges[i].second, mesh->GetEdge(i).second);
    }
}

TEST(Mesh2D, Mesh2DToCurvilinear_WithRectangularMesh_ShouldCreateFullCurvilinearMesh)
{
    // Prepare
    const auto mesh = MakeRectangularMeshForTesting(3,
                                                    3,
                                                    10.0,
                                                    10.0,
                                                    meshkernel::Projection::cartesian);
    meshkernel::Mesh2DToCurvilinear mesh2DToCurvilinear(*mesh);

    // Execute
    const meshkernel::Point point(5.0, 5.0);
    const auto curvilinearGrid = mesh2DToCurvilinear.Compute(point);

    // Assert
    ASSERT_EQ(3, curvilinearGrid->NumM());
    ASSERT_EQ(3, curvilinearGrid->NumN());
    ASSERT_EQ(9, curvilinearGrid->GetNumNodes());

    ASSERT_EQ(10.0, curvilinearGrid->GetNode(0, 2).x);
    ASSERT_EQ(10.0, curvilinearGrid->GetNode(0, 2).y);

    ASSERT_EQ(10.0, curvilinearGrid->GetNode(0, 1).x);
    ASSERT_EQ(5.0, curvilinearGrid->GetNode(0, 1).y);

    ASSERT_EQ(10.0, curvilinearGrid->GetNode(0, 0).x);
    ASSERT_EQ(0.0, curvilinearGrid->GetNode(0, 0).y);
}

TEST(Mesh2D, Mesh2DToCurvilinear_WithStartingPointOutsideMesh_ShouldThrowAnException)
{
    // Prepare
    const auto mesh = MakeRectangularMeshForTesting(3,
                                                    3,
                                                    10.0,
                                                    10.0,
                                                    meshkernel::Projection::cartesian);
    meshkernel::Mesh2DToCurvilinear mesh2DToCurvilinear(*mesh);

    // Execute
    const meshkernel::Point point(-20.0, -20.0);

    // Assert
    EXPECT_THROW(mesh2DToCurvilinear.Compute(point), meshkernel::AlgorithmError);
}

TEST(Mesh2D, Mesh2DToCurvilinear_WithMixedMesh_ShouldCreatePartialCurvilinearMesh)
{
    // Prepare a mixed mesh with two triangles at the boundary
    std::vector<meshkernel::Point> nodes;
    nodes.push_back({-10.0, 0.0});
    nodes.push_back({0.0, 0.0});
    nodes.push_back({10.0, 0.0});
    nodes.push_back({20.0, 0.0});
    nodes.push_back({0.0, 10.0});
    nodes.push_back({10.0, 10.0});
    nodes.push_back({20.0, 10.0});
    nodes.push_back({0.0, 20.0});
    nodes.push_back({10.0, 20.0});

    std::vector<meshkernel::Edge> edges;
    edges.push_back({0, 4});
    edges.push_back({0, 1});
    edges.push_back({1, 4});
    edges.push_back({1, 2});
    edges.push_back({2, 5});
    edges.push_back({2, 3});
    edges.push_back({3, 6});
    edges.push_back({4, 7});
    edges.push_back({4, 5});
    edges.push_back({5, 8});
    edges.push_back({5, 6});
    edges.push_back({7, 8});
    edges.push_back({6, 8});

    // 2 Execute
    auto mesh = meshkernel::Mesh2D(edges, nodes, meshkernel::Projection::cartesian);
    meshkernel::Mesh2DToCurvilinear mesh2DToCurvilinear(mesh);
    const meshkernel::Point point(5.0, 5.0);
    const auto curvilinearGrid = mesh2DToCurvilinear.Compute(point);

    // Assert
    ASSERT_EQ(3, curvilinearGrid->NumM());
    ASSERT_EQ(3, curvilinearGrid->NumN());
    ASSERT_EQ(9, curvilinearGrid->GetNumNodes());

    ASSERT_EQ(20.0, curvilinearGrid->GetNode(0, 0).x);
    ASSERT_EQ(0.0, curvilinearGrid->GetNode(0, 0).y);
    ASSERT_EQ(10.0, curvilinearGrid->GetNode(1, 0).x);
    ASSERT_EQ(0.0, curvilinearGrid->GetNode(1, 0).y);
    ASSERT_EQ(0.0, curvilinearGrid->GetNode(2, 0).x);
    ASSERT_EQ(0.0, curvilinearGrid->GetNode(2, 0).y);

    ASSERT_EQ(20.0, curvilinearGrid->GetNode(0, 1).x);
    ASSERT_EQ(10.0, curvilinearGrid->GetNode(0, 1).y);
    ASSERT_EQ(10.0, curvilinearGrid->GetNode(1, 1).x);
    ASSERT_EQ(10.0, curvilinearGrid->GetNode(1, 1).y);
    ASSERT_EQ(0.0, curvilinearGrid->GetNode(2, 1).x);
    ASSERT_EQ(10.0, curvilinearGrid->GetNode(2, 1).y);

    ASSERT_EQ(-999.0, curvilinearGrid->GetNode(0, 2).x);
    ASSERT_EQ(-999.0, curvilinearGrid->GetNode(0, 2).y);
    ASSERT_EQ(10.0, curvilinearGrid->GetNode(1, 2).x);
    ASSERT_EQ(20.0, curvilinearGrid->GetNode(1, 2).y);
    ASSERT_EQ(0.0, curvilinearGrid->GetNode(2, 2).x);
    ASSERT_EQ(20.0, curvilinearGrid->GetNode(2, 2).y);
}

TEST(Mesh2D, GetBoundingBox_WithANonEmptyMesh_ShouldGetAValidBoundingBox)
{
    // Prepare
    const auto mesh = MakeRectangularMeshForTesting(10,
                                                    10,
                                                    10.0,
                                                    10.0,
                                                    meshkernel::Projection::cartesian);
    // Execute
    const auto boundingBox = mesh->GetBoundingBox();

    // Assert
    ASSERT_EQ(boundingBox.lowerLeft().x, 0.0);
    ASSERT_EQ(boundingBox.lowerLeft().y, 0.0);
    ASSERT_EQ(boundingBox.upperRight().x, 10.0);
    ASSERT_EQ(boundingBox.upperRight().y, 10.0);
}

TEST(Mesh2D, GetBoundingBox_WithAnEmptyMesh_ShouldGetAnEmptyBoundingBox)
{
    // Prepare
    const auto mesh = MakeRectangularMeshForTesting(10,
                                                    10,
                                                    10.0,
                                                    10.0,
                                                    meshkernel::Projection::cartesian);

    const auto polygon = meshkernel::Polygons({}, meshkernel::Projection::cartesian);
    const auto deletionOption = meshkernel::Mesh2D::DeleteMeshOptions::InsideNotIntersected;

    meshkernel::Point lowerLeft(std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
    meshkernel::Point upperRight(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest());

    // Execute
    auto undoAction = mesh->DeleteMesh(polygon, deletionOption, false);
    auto boundingBox = mesh->GetBoundingBox();

    // Assert
    ASSERT_EQ(boundingBox.lowerLeft().x, lowerLeft.x);
    ASSERT_EQ(boundingBox.lowerLeft().y, lowerLeft.y);
    ASSERT_EQ(boundingBox.upperRight().x, upperRight.x);
    ASSERT_EQ(boundingBox.upperRight().y, upperRight.y);
}

TEST(Mesh2D, GetEdgesBoundingBox_WithAnInvalidEdge_ShouldGetOneInvalidEdgeBoundingBox)
{
    // Prepare
    const auto mesh = MakeRectangularMeshForTesting(10,
                                                    10,
                                                    10.0,
                                                    10.0,
                                                    meshkernel::Projection::cartesian);
    // Execute
    const auto undoAction = mesh->DeleteEdge(0);
    const auto edgesBoundingBoxes = mesh->GetEdgesBoundingBoxes();

    // Assert
    meshkernel::Point lowerLeft(std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
    meshkernel::Point upperRight(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest());

    const double tolerance = 1e-6;
    ASSERT_NEAR(edgesBoundingBoxes[0].lowerLeft().x, lowerLeft.x, tolerance);
    ASSERT_NEAR(edgesBoundingBoxes[0].lowerLeft().y, lowerLeft.y, tolerance);
    ASSERT_NEAR(edgesBoundingBoxes[0].upperRight().x, upperRight.x, tolerance);
    ASSERT_NEAR(edgesBoundingBoxes[0].upperRight().y, upperRight.y, tolerance);

    ASSERT_NEAR(edgesBoundingBoxes[1].lowerLeft().x, 0.0, tolerance);
    ASSERT_NEAR(edgesBoundingBoxes[1].lowerLeft().y, 1.1111111111111112, tolerance);
    ASSERT_NEAR(edgesBoundingBoxes[1].upperRight().x, 1.1111111111111112, tolerance);
    ASSERT_NEAR(edgesBoundingBoxes[1].upperRight().y, 1.1111111111111112, tolerance);
}

TEST(Mesh2D, GetSmoothness_OnTriangularMesh_ShouldgetSmoothnessValues)
{
    // Setup
    const auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/TestOrthogonalizationMediumTriangularGrid_net.nc");

    // Execute
    const auto smoothness = meshkernel::MeshSmoothness::Compute(*mesh);

    // Assert
    const double tolerance = 1e-6;
    ASSERT_NEAR(1.0000000000000047, smoothness[0], tolerance);
    ASSERT_NEAR(1.5393847629344886, smoothness[10], tolerance);
    ASSERT_NEAR(1.1609660187036754, smoothness[20], tolerance);
    ASSERT_NEAR(1.4420158602682915, smoothness[30], tolerance);
}

TEST(Mesh2D, GetOrthogonality_OnTriangularMesh_ShouldGetOrthogonalityValues)
{
    // Setup
    const auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/TestOrthogonalizationMediumTriangularGrid_net.nc");

    // Execute
    const auto orthogonality = meshkernel::MeshOrthogonality::Compute(*mesh);

    // Assert
    const double tolerance = 1e-6;
    ASSERT_NEAR(1.0566340037701503e-15, orthogonality[0], tolerance);
    ASSERT_NEAR(0.052159566591519289, orthogonality[10], tolerance);
    ASSERT_NEAR(1.0342915752434056e-15, orthogonality[20], tolerance);
    ASSERT_NEAR(0.045878303256790140, orthogonality[30], tolerance);
}

TEST(Mesh2D, MeshToCurvilinear_OnRealMesh_ShouldConvertCurvilinearPart)
{
    // Prepare
    const auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/MeshToCurvilinear.nc");

    meshkernel::Mesh2DToCurvilinear mesh2DToCurvilinear(*mesh);

    // Execute
    const auto result = mesh2DToCurvilinear.Compute({10155.18, 391781.98});

    // Assert
    EXPECT_EQ(result->NumM(), 38);
    EXPECT_EQ(result->NumN(), 45);
}

TEST(Mesh2D, Mesh2DComputeAspectRatio)
{
    const double tolerance = 1.0e-12;
    const meshkernel::UInt size = 3;
    const double dimension = 10.0;

    const auto mesh = MakeRectangularMeshForTesting(size,
                                                    size,
                                                    dimension,
                                                    dimension,
                                                    meshkernel::Projection::cartesian);

    std::vector<meshkernel::Point> displacement{{0.573312665029025, 0.164422234673451},
                                                {0.848580895925857, 0.273698982765599},
                                                {0.649270465028434, 1.16836612028342},
                                                {0.662125241426321, 0.0432151380810593},
                                                {0.0835527965782094, 0.00962273257699897},
                                                {1.16304561871212, 0.858465890511307},
                                                {0.817398702725653, 0.658660972199883},
                                                {0.952747549998293, 0.876488243123971},
                                                {0.41029278270001, 0.0593306417328899}};

    for (meshkernel::UInt i = 0; i < mesh->GetNumNodes(); ++i)
    {
        mesh->SetNode(i, mesh->Node(i) + displacement[i]);
    }

    // The grid nodes have been displaced by some amount, the face circumcentres will need to be recomputed.
    mesh->Administrate();

    std::vector<double> aspectRatios;
    // Values calculated by the algorithm, not derived analytically
    std::vector<double> expectedAspectRatios{0.909799624513058, 1.36963998294878, 1.08331064761803,
                                             0.896631398796997, 0.839834200634414, 1.15475768474815,
                                             0.832219560848176, 0.819704195029218, 1.11720404458871,
                                             0.807269974963293, 0.972997967873087, 1.17917808082661};

    mesh->ComputeAspectRatios(aspectRatios);

    ASSERT_EQ(aspectRatios.size(), expectedAspectRatios.size());

    for (size_t i = 0; i < expectedAspectRatios.size(); ++i)
    {
        EXPECT_NEAR(aspectRatios[i], expectedAspectRatios[i], tolerance);
    }
}

TEST(Mesh2D, MeshToCurvilinear_SingleElement)
{
    // Test steps
    // - generate 3 x3 element mesh
    // - triangulate all elements except the centre element
    // - compute curvilinear grid with point in centre element, generated should be a single element
    // - add some tests checking for expected failure cases

    const int n = 4; // x
    const int m = 4; // y

    std::vector<std::vector<int>> indicesValues(n, std::vector<int>(m));
    std::vector<meshkernel::Point> nodes(n * m);
    std::size_t nodeIndex = 0;
    for (auto j = 0; j < m; ++j)
    {
        for (auto i = 0; i < n; ++i)
        {
            indicesValues[i][j] = i + j * n;
            nodes[nodeIndex] = {(double)i, (double)j};
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

    // leave centre quadrilateral element untouched, all others will be split into two triangles
    std::vector<meshkernel::Edge> toConnect{{0, 5}, {1, 6}, {2, 7}, {4, 9}, {6, 11}, {8, 13}, {9, 14}, {10, 15}};

    meshkernel::Mesh2D mesh(edges, nodes, meshkernel::Projection::cartesian);

    for (size_t i = 0; i < toConnect.size(); ++i)
    {
        [[maybe_unused]] auto undo = mesh.ConnectNodes(toConnect[i].first, toConnect[i].second, false);
    }

    meshkernel::Mesh2DToCurvilinear mesh2DToCurvilinear(mesh);

    // Execute
    const auto result = mesh2DToCurvilinear.Compute({1.5, 1.5});

    // Assert
    ASSERT_EQ(result->NumM(), 2);
    ASSERT_EQ(result->NumN(), 2);

    EXPECT_EQ(result->GetNode(0, 0).x, 2.0);
    EXPECT_EQ(result->GetNode(0, 0).y, 1.0);

    EXPECT_EQ(result->GetNode(0, 1).x, 2.0);
    EXPECT_EQ(result->GetNode(0, 1).y, 2.0);

    EXPECT_EQ(result->GetNode(1, 1).x, 1.0);
    EXPECT_EQ(result->GetNode(1, 1).y, 2.0);

    EXPECT_EQ(result->GetNode(1, 0).x, 1.0);
    EXPECT_EQ(result->GetNode(1, 0).y, 1.0);

    // Execute several expected failure tests
    // Element is not quadrilateral
    EXPECT_THROW([[maybe_unused]] auto result = mesh2DToCurvilinear.Compute({0.5, 0.45}), meshkernel::AlgorithmError);
    EXPECT_THROW([[maybe_unused]] auto result = mesh2DToCurvilinear.Compute({2.5, 0.45}), meshkernel::AlgorithmError);
    // Point is outside domain
    EXPECT_THROW([[maybe_unused]] auto result = mesh2DToCurvilinear.Compute({-1.0, -1.0}), meshkernel::AlgorithmError);
    EXPECT_THROW([[maybe_unused]] auto result = mesh2DToCurvilinear.Compute({4.0, 4.0}), meshkernel::AlgorithmError);
}
