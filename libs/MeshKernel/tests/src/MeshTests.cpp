#include <chrono>
#include <gtest/gtest.h>
#include <random>

#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Polygons.hpp>
#include <TestUtils/Definitions.hpp>
#include <TestUtils/MakeMeshes.hpp>

TEST(Mesh, OneQuadTestConstructor)
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
    ASSERT_EQ(2, mesh.m_nodesNumEdges[0]);
    ASSERT_EQ(2, mesh.m_nodesNumEdges[1]);
    ASSERT_EQ(2, mesh.m_nodesNumEdges[2]);
    ASSERT_EQ(2, mesh.m_nodesNumEdges[3]);

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

    // the found circumcenter for the face
    ASSERT_DOUBLE_EQ(5.0, mesh.m_facesCircumcenters[0].x);
    ASSERT_DOUBLE_EQ(5.0, mesh.m_facesCircumcenters[0].y);

    // each edge has only one face in this case
    ASSERT_EQ(1, mesh.m_edgesNumFaces[0]);
    ASSERT_EQ(1, mesh.m_edgesNumFaces[1]);
    ASSERT_EQ(1, mesh.m_edgesNumFaces[2]);
    ASSERT_EQ(1, mesh.m_edgesNumFaces[3]);

    // each edge is a boundary edge, so the second entry of edgesFaces is an invalid index (meshkernel::constants::missing::sizetValue)
    ASSERT_EQ(meshkernel::constants::missing::uintValue, mesh.m_edgesFaces[0][1]);
    ASSERT_EQ(meshkernel::constants::missing::uintValue, mesh.m_edgesFaces[1][1]);
    ASSERT_EQ(meshkernel::constants::missing::uintValue, mesh.m_edgesFaces[2][1]);
    ASSERT_EQ(meshkernel::constants::missing::uintValue, mesh.m_edgesFaces[3][1]);
}

TEST(Mesh2D, TriangulateSamplesWithSkinnyTriangle)
{
    // Prepare
    std::vector<meshkernel::Point> nodes;

    nodes.push_back({302.002502, 472.130371});
    nodes.push_back({144.501526, 253.128174});
    nodes.push_back({368.752930, 112.876755});
    nodes.push_back({707.755005, 358.879242});
    nodes.push_back({301.252502, 471.380371});
    nodes.push_back({302.002502, 472.130371});

    meshkernel::Polygons polygons(nodes, meshkernel::Projection::cartesian);

    // Execute
    const auto generatedPoints = polygons.ComputePointsInPolygons();

    meshkernel::Mesh2D mesh(generatedPoints[0], polygons, meshkernel::Projection::cartesian);

    // Assert
    ASSERT_EQ(5, mesh.GetNumNodes());
    ASSERT_EQ(6, mesh.GetNumEdges());

    ASSERT_EQ(3, mesh.GetEdge(0).first);
    ASSERT_EQ(0, mesh.GetEdge(0).second);

    ASSERT_EQ(0, mesh.GetEdge(1).first);
    ASSERT_EQ(1, mesh.GetEdge(1).second);

    ASSERT_EQ(1, mesh.GetEdge(2).first);
    ASSERT_EQ(3, mesh.GetEdge(2).second);

    ASSERT_EQ(4, mesh.GetEdge(3).first);
    ASSERT_EQ(1, mesh.GetEdge(3).second);

    ASSERT_EQ(1, mesh.GetEdge(4).first);
    ASSERT_EQ(2, mesh.GetEdge(4).second);

    ASSERT_EQ(2, mesh.GetEdge(5).first);
    ASSERT_EQ(4, mesh.GetEdge(5).second);
}

TEST(Mesh, TriangulateSamples)
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

TEST(Mesh, TwoTrianglesDuplicatedEdges)
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

TEST(Mesh, MeshBoundaryToPolygon)
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

TEST(Mesh, HangingEdge)
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

    auto mesh = meshkernel::Mesh2D(edges, nodes, meshkernel::Projection::cartesian);

    ASSERT_EQ(1, mesh.GetNumFaces());
}

TEST(Mesh, NodeMerging)
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

            // add artificial edges
            auto edge = mesh->GetEdge(mesh->m_nodesEdges[originalNodeIndex][0]);
            auto otherNode = edge.first + edge.second - originalNodeIndex;

            edges[edgeIndex] = {nodeIndex, otherNode};
            edgeIndex++;
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

TEST(Mesh, MillionQuads)
{
    const int n = 4; // x
    const int m = 4; // y

    std::vector<std::vector<meshkernel::UInt>> indicesValues(n, std::vector<meshkernel::UInt>(m));
    std::vector<meshkernel::Point> nodes(n * m);
    std::size_t nodeIndex = 0;
    for (meshkernel::UInt j = 0; j < m; ++j)
    {
        for (meshkernel::UInt i = 0; i < n; ++i)
        {
            indicesValues[i][j] = i + j * n;
            nodes[nodeIndex] = {static_cast<double>(i), static_cast<double>(j)};
            nodeIndex++;
        }
    }

    std::vector<meshkernel::Edge> edges((n - 1) * m + (m - 1) * n);
    std::size_t edgeIndex = 0;
    for (meshkernel::UInt j = 0; j < m; ++j)
    {
        for (meshkernel::UInt i = 0; i < n - 1; ++i)
        {
            edges[edgeIndex] = {indicesValues[i][j], indicesValues[i + 1][j]};
            edgeIndex++;
        }
    }

    for (meshkernel::UInt j = 0; j < m - 1; ++j)
    {
        for (meshkernel::UInt i = 0; i < n; ++i)
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
    // std::cout << "Elapsed time " << elapsedTime << " s " << std::endl;
    // std::cout << "Number of found cells " << mesh.GetNumFaces() << std::endl;

    EXPECT_LE(elapsedTime, 5.0);
}

TEST(Mesh, InsertNodeInMeshWithExistingNodesRtreeTriggersRTreeReBuild)
{
    // Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);
    mesh->BuildTree(meshkernel::Location::Nodes);

    // insert nodes modifies the number of nodes, m_nodesRTreeRequiresUpdate is set to true
    meshkernel::Point newPoint{10.0, 10.0};

    auto [newNodeIndex, insertAction] = mesh->InsertNode(newPoint);

    [[maybe_unused]] auto connectAction = mesh->ConnectNodes(0, newNodeIndex);

    // when m_nodesRTreeRequiresUpdate = true m_nodesRTree is not empty the mesh.m_nodesRTree is re-build
    mesh->Administrate();

    // builds edges RTree
    mesh->BuildTree(meshkernel::Location::Edges);
    const auto& rtreeEdges = mesh->GetRTree(meshkernel::Location::Edges);

    // even if m_edgesRTreeRequiresUpdate = true, m_edgesRTree is initially empty, so it is assumed that is not needed for searches
    ASSERT_EQ(5, rtreeEdges.Size());
}

TEST(Mesh, DeleteNodeInMeshWithExistingNodesRtreeTriggersRTreeReBuild)
{
    // Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    meshkernel::Point newPoint{10.0, 10.0};
    mesh->BuildTree(meshkernel::Location::Nodes);
    auto& rtree = mesh->GetRTree(meshkernel::Location::Nodes);

    [[maybe_unused]] auto [nodeId, indertAction] = mesh->InsertNode(newPoint);

    // delete nodes modifies the number of nodes, m_nodesRTreeRequiresUpdate is set to true
    [[maybe_unused]] auto deleteAction = mesh->DeleteNode(0);

    // when m_nodesRTreeRequiresUpdate
    mesh->DeleteInvalidNodesAndEdges();
    mesh->Administrate();

    // building a tree based on nodes
    rtree.BuildTree(mesh->Nodes());

    // After deleting a node, the nodes RTree is reduced
    ASSERT_EQ(3, rtree.Size());
}

TEST(Mesh, ConnectNodesInMeshWithExistingEdgesRtreeTriggersRTreeReBuild)
{
    // 1 Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    meshkernel::Point newPoint{10.0, 10.0};

    auto [newNodeIndex, insertAction] = mesh->InsertNode(newPoint);

    // connect nodes modifies the number of edges, m_nodesRTreeRequiresUpdate is set to true
    [[maybe_unused]] auto connectAction = mesh->ConnectNodes(0, newNodeIndex);

    // re-do mesh adminstration
    mesh->Administrate();

    // re-build tree
    mesh->BuildTree(meshkernel::Location::Edges);
    const auto& rtree = mesh->GetRTree(meshkernel::Location::Edges);

    // even if m_nodesRTreeRequiresUpdate = true, m_nodesRTree is initially empty, so it is assumed that is not needed for searches
    ASSERT_EQ(5, rtree.Size());
}

TEST(Mesh, DeleteEdgeInMeshWithExistingEdgesRtreeTriggersRTreeReBuild)
{
    // 1 Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);
    mesh->BuildTree(meshkernel::Location::Edges);
    const auto& rtree = mesh->GetRTree(meshkernel::Location::Edges);

    // DeleteEdge modifies the number of edges, m_edgesRTreeRequiresUpdate is set to true
    [[maybe_unused]] auto action = mesh->DeleteEdge(0);

    // re-do mesh administration
    mesh->Administrate();

    // re-build tree
    mesh->BuildTree(meshkernel::Location::Edges);

    // deleting an edge produces an edges rtree of size 3
    ASSERT_EQ(3, rtree.Size());
}

TEST(Mesh, InsertUnconnectedNodeInMeshIsSetToInvalid)
{
    // Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);
    mesh->BuildTree(meshkernel::Location::Nodes);
    const auto& rtreesNodes = mesh->GetRTree(meshkernel::Location::Nodes);

    // insert nodes modifies the number of nodes, m_nodesRTreeRequiresUpdate is set to true
    meshkernel::Point newPoint{10.0, 10.0};

    [[maybe_unused]] auto [nodeId, action] = mesh->InsertNode(newPoint);

    // when m_nodesRTreeRequiresUpdate = true m_nodesRTree is not empty the mesh.m_nodesRTree is re-build
    mesh->Administrate();

    // building a tree based on nodes
    mesh->BuildTree(meshkernel::Location::Nodes);

    // building the edges rtree
    mesh->BuildTree(meshkernel::Location::Edges);
    const auto& rtreeEdges = mesh->GetRTree(meshkernel::Location::Edges);

    // building a tree based on edges
    EXPECT_EQ(5, mesh->GetNumNodes());
    EXPECT_EQ(4, mesh->GetNumValidNodes());
    EXPECT_EQ(4, rtreesNodes.Size());
    EXPECT_EQ(4, rtreeEdges.Size());
    // Administrate should set the unconnected node to be invalid.
    EXPECT_FALSE(mesh->Node(4).IsValid());
}

TEST(Mesh, EdgeConnectedToInvalidNodeInMeshIsSetToInvalid)
{
    // Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);
    mesh->BuildTree(meshkernel::Location::Nodes);
    const auto& nodesRtree = mesh->GetRTree(meshkernel::Location::Nodes);
    const auto& edgesRtree = mesh->GetRTree(meshkernel::Location::Edges);

    meshkernel::Point newPoint{meshkernel::constants::missing::doubleValue,
                               meshkernel::constants::missing::doubleValue};
    auto [nodeIndex, nodeAction] = mesh->InsertNode(newPoint);
    auto [edgeIndex, edgeAction] = mesh->ConnectNodes(0, nodeIndex);

    EXPECT_EQ(mesh->GetEdge(edgeIndex).first, 0);
    EXPECT_EQ(mesh->GetEdge(edgeIndex).second, nodeIndex);

    // when m_nodesRTreeRequiresUpdate = true m_nodesRTree is not empty the mesh.m_nodesRTree is re-build
    mesh->Administrate();

    // building a tree based on nodes
    mesh->BuildTree(meshkernel::Location::Nodes);

    // building a tree based on edges
    mesh->BuildTree(meshkernel::Location::Edges);

    EXPECT_EQ(5, mesh->GetNumNodes());
    EXPECT_EQ(4, mesh->GetNumValidNodes());

    EXPECT_EQ(5, mesh->GetNumEdges());
    EXPECT_EQ(4, mesh->GetNumValidEdges());

    EXPECT_EQ(4, nodesRtree.Size());
    EXPECT_EQ(4, edgesRtree.Size());

    // Administrate should set the unconnected node to be invalid.
    EXPECT_FALSE(mesh->Node(4).IsValid());

    // Administrate should set the edge connecting an invalid node to be invalid.
    EXPECT_EQ(mesh->GetEdge(4).first, meshkernel::constants::missing::uintValue);
    EXPECT_EQ(mesh->GetEdge(4).second, meshkernel::constants::missing::uintValue);
}

TEST(Mesh, GetNodeIndexShouldTriggerNodesRTreeBuild)
{
    // 1 Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    // By default, no nodesRTree is build
    const auto& edgesRTree = mesh->GetRTree(meshkernel::Location::Edges);
    const auto& nodesRTree = mesh->GetRTree(meshkernel::Location::Nodes);

    ASSERT_EQ(0, nodesRTree.Size());
    ASSERT_EQ(0, edgesRTree.Size());

    // FindNodeCloseToAPoint builds m_nodesRTree for searching the nodes
    const auto index = mesh->FindNodeCloseToAPoint({1.5, 1.5}, 10.0);
    ASSERT_EQ(3, index);

    // m_nodesRTree is build
    ASSERT_EQ(4, nodesRTree.Size());

    // m_edgesRTree is not build when searching for nodes
    ASSERT_EQ(0, edgesRTree.Size());
}

TEST(Mesh, FindEdgeCloseToAPointShouldTriggerEdgesRTreeBuild)
{
    // 1 Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    // FindEdgeCloseToAPoint builds m_edgesRTree for searching the edges

    mesh->BuildTree(meshkernel::Location::Edges);
    const auto& edgesRTree = mesh->GetRTree(meshkernel::Location::Edges);
    const auto& nodesRTree = mesh->GetRTree(meshkernel::Location::Nodes);

    const auto index = mesh->FindLocationIndex({1.5, 1.5}, meshkernel::Location::Edges);
    ASSERT_EQ(1, index);

    // m_nodesRTree is not build when searching for edges
    ASSERT_EQ(0, nodesRTree.Size());

    // m_edgesRTree is build
    ASSERT_EQ(4, edgesRTree.Size());
}

TEST(Mesh, GetObtuseTriangles)
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

    auto mesh = std::make_unique<meshkernel::Mesh2D>(edges, nodes, meshkernel::Projection::cartesian);

    // execute, only one obtuse triangle should be found
    const auto obtuseTrianglesCount = mesh->GetObtuseTrianglesCenters().size();

    // assert a small flow edge is found
    ASSERT_EQ(1, obtuseTrianglesCount);
}

TEST(Mesh, GetSmallFlowEdgeCenters)
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

TEST(Mesh, DeleteSmallFlowEdge)
{
    // Setup a mesh with eight triangles
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/RemoveSmallFlowEdgesTests/remove_small_flow_edges_net.nc");

    ASSERT_EQ(8, mesh->GetNumFaces());

    // After merging the number of faces is reduced
    auto undoAction = mesh->DeleteSmallFlowEdges(1.0);

    ASSERT_EQ(3, mesh->GetNumFaces());

    // Test the undo action has been computed correctly
    undoAction->Restore();
    // Recompute faces
    mesh->Administrate();
    ASSERT_EQ(8, mesh->GetNumFaces());
}

TEST(Mesh, DeleteSmallTrianglesAtBoundaries)
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

    // Test the undo action has been computed correctly
    undoAction->Restore();
    // Recompute faces
    mesh->Administrate();
    ASSERT_EQ(2, mesh->GetNumFaces());
}

TEST(Mesh, DeleteHangingEdge)
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

    // Execute
    auto mesh = meshkernel::Mesh2D(edges, nodes, meshkernel::Projection::cartesian);

    // Assert
    ASSERT_EQ(1, mesh.GetNumFaces());
    ASSERT_EQ(4, mesh.GetNumEdges());

    // Execute
    auto hangingEdges = mesh.GetHangingEdges();

    // Assert
    ASSERT_EQ(1, hangingEdges.size());

    // Execute
    auto undoAction = mesh.DeleteHangingEdges();
    hangingEdges = mesh.GetHangingEdges();

    // Assert
    ASSERT_EQ(0, hangingEdges.size());

    // Test the undo action has been computed correctly
    undoAction->Restore();
    // Recompute faces
    mesh.Administrate();
    // Assert
    ASSERT_EQ(1, mesh.GetNumFaces());
    ASSERT_EQ(4, mesh.GetNumEdges());
}

class MeshDeletion : public ::testing::TestWithParam<std::tuple<meshkernel::Mesh2D::DeleteMeshOptions, bool, int>>
{
public:
    [[nodiscard]] static std::vector<std::tuple<meshkernel::Mesh2D::DeleteMeshOptions, bool, int>> GetData()
    {
        return {
            {meshkernel::Mesh2D::DeleteMeshOptions::InsideNotIntersected, false, 16},
            {meshkernel::Mesh2D::DeleteMeshOptions::InsideAndIntersected, false, 14},
            {meshkernel::Mesh2D::DeleteMeshOptions::InsideNotIntersected, true, 6},
            {meshkernel::Mesh2D::DeleteMeshOptions::InsideAndIntersected, true, 0}

        };
    }
};

TEST_P(MeshDeletion, expected_results)
{
    // Get the test parameters
    auto const [deleteOption, invertSelection, numNodes] = GetParam();

    // Setup
    auto mesh = MakeRectangularMeshForTesting(4, 4, 1.0, meshkernel::Projection::cartesian);

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    std::vector<meshkernel::Point> polygonNodes{
        {-0.5, -1.0},
        {0.8, -1.0},
        {0.8, 1.8},
        {-0.5, 1.8},
        {-0.5, -1.0}};

    const meshkernel::Polygons polygon(polygonNodes, meshkernel::Projection::cartesian);

    // Execute
    auto undoAction = mesh->DeleteMesh(polygon, deleteOption, invertSelection);

    // Assert
    const auto numValidNodes = mesh->GetNumValidNodes();
    ASSERT_EQ(numNodes, numValidNodes);

    // Test the undo action has been computed correctly
    undoAction->Restore();
    // Recompute faces
    mesh->Administrate();

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

INSTANTIATE_TEST_SUITE_P(Mesh, MeshDeletion, ::testing::ValuesIn(MeshDeletion::GetData()));

class MeshDeletionWithInnerPolygons : public ::testing::TestWithParam<std::tuple<meshkernel::Mesh2D::DeleteMeshOptions, bool, std::vector<meshkernel::Point>, int>>
{

    static inline std::vector<meshkernel::Point> single_polygon_{
        {-0.722886114680926, 2.22765832371444},
        {1.50244688883463, 2.65710855246306},
        {3.04456361934103, 1.91533088462454},
        {3.47401384808964, 5.56565782898778},
        {-1.42562285263321, 5.70230108358961},
        {-0.722886114680926, 2.22765832371444},
    };

    static inline std::vector<meshkernel::Point> double_polygon_{
        {-0.722886114680926, 2.22765832371444},
        {1.50244688883463, 2.65710855246306},
        {3.04456361934103, 1.91533088462454},
        {3.47401384808964, 5.56565782898778},
        {-1.42562285263321, 5.70230108358961},
        {-0.722886114680926, 2.22765832371444},
        {meshkernel::constants::missing::doubleValue, meshkernel::constants::missing::doubleValue},
        {3.59021337251456, -1.90949318892618},
        {6.02802927483664, -1.86121960670198},
        {5.88320852816404, 2.09721413568238},
        {2.55233135469427, 0.600733086732198},
        {3.59021337251456, -1.90949318892618},
    };

public:
    [[nodiscard]] static std::vector<std::tuple<meshkernel::Mesh2D::DeleteMeshOptions, bool, std::vector<meshkernel::Point>, int>> GetData()
    {
        return {
            {meshkernel::Mesh2D::DeleteMeshOptions::InsideAndIntersected, false, single_polygon_, 23},
            {meshkernel::Mesh2D::DeleteMeshOptions::InsideNotIntersected, false, single_polygon_, 30},
            {meshkernel::Mesh2D::DeleteMeshOptions::InsideAndIntersected, true, single_polygon_, 12},
            {meshkernel::Mesh2D::DeleteMeshOptions::InsideNotIntersected, true, single_polygon_, 23},
            {meshkernel::Mesh2D::DeleteMeshOptions::InsideAndIntersected, false, double_polygon_, 17},
            {meshkernel::Mesh2D::DeleteMeshOptions::InsideNotIntersected, false, double_polygon_, 29},
            {meshkernel::Mesh2D::DeleteMeshOptions::InsideAndIntersected, true, double_polygon_, 16},
            {meshkernel::Mesh2D::DeleteMeshOptions::InsideNotIntersected, true, double_polygon_, 29}

        };
    }
};

TEST_P(MeshDeletionWithInnerPolygons, expected_results)
{
    // Get the test parameters
    auto const& [deleteOption, invertSelection, polygonNodes, numNodes] = GetParam();

    // Setup
    auto mesh = MakeRectangularMeshForTesting(6, 6, 1.0, meshkernel::Projection::cartesian);

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    const meshkernel::Polygons polygon(polygonNodes, meshkernel::Projection::cartesian);

    // Execute
    auto undoAction = mesh->DeleteMesh(polygon, deleteOption, invertSelection);

    // Assert
    const auto numValidNodes = mesh->GetNumValidNodes();
    ASSERT_EQ(numNodes, numValidNodes);

    // Test the undo action has been computed correctly
    undoAction->Restore();
    // Recompute faces
    mesh->Administrate();

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

INSTANTIATE_TEST_SUITE_P(Mesh, MeshDeletionWithInnerPolygons, ::testing::ValuesIn(MeshDeletionWithInnerPolygons::GetData()));
