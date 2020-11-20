#include <MeshKernel/Mesh.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/Constants.hpp>
#include <TestUtils/MakeMeshes.hpp>
#include <gtest/gtest.h>
#include <chrono>
#include <random>

TEST(Mesh, OneQuadTestConstructor)
{
    //1 Setup
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

    meshkernel::Mesh mesh;

    // 2 Execution
    mesh.Set(edges, nodes, meshkernel::Projections::cartesian);

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

    //each edge is a boundary edge, so the second entry of edgesFaces is an invalid index (in C++, -1)
    ASSERT_EQ(-1, mesh.m_edgesFaces[0][1]);
    ASSERT_EQ(-1, mesh.m_edgesFaces[1][1]);
    ASSERT_EQ(-1, mesh.m_edgesFaces[2][1]);
    ASSERT_EQ(-1, mesh.m_edgesFaces[3][1]);
}

TEST(Mesh, MakeMeshInPolygon)
{
    //1 Setup
    std::vector<meshkernel::Point> nodes;

    nodes.push_back({302.002502, 472.130371});
    nodes.push_back({144.501526, 253.128174});
    nodes.push_back({368.752930, 112.876755});
    nodes.push_back({707.755005, 358.879242});
    nodes.push_back({301.252502, 471.380371});
    nodes.push_back({302.002502, 472.130371});

    meshkernel::Polygons polygons(nodes, meshkernel::Projections::cartesian);

    meshkernel::Mesh mesh;
    meshkernelapi::MakeGridParametersNative makeGridParametersNative;
    makeGridParametersNative.GridType = 0;
    makeGridParametersNative.GridAngle = 0.0;
    makeGridParametersNative.OriginXCoordinate = 0.0;
    makeGridParametersNative.OriginYCoordinate = 0.0;
    makeGridParametersNative.OriginZCoordinate = 0.0;
    makeGridParametersNative.NumberOfColumns = 3;
    makeGridParametersNative.NumberOfRows = 3;
    makeGridParametersNative.XGridBlockSize = 100.0;
    makeGridParametersNative.YGridBlockSize = 100.0;

    // 2 Execution
    mesh.MakeMesh(makeGridParametersNative, polygons);
    ASSERT_EQ(43, mesh.GetNumEdges());
    ASSERT_EQ(27, mesh.GetNumNodes());
}

TEST(Mesh, MakeMeshInPolygonSpherical)
{
    //1 Setup
    std::vector<meshkernel::Point> nodes;

    nodes.push_back({302.002502, 472.130371});
    nodes.push_back({144.501526, 253.128174});
    nodes.push_back({368.752930, 112.876755});
    nodes.push_back({707.755005, 358.879242});
    nodes.push_back({301.252502, 471.380371});
    nodes.push_back({302.002502, 472.130371});

    meshkernel::Polygons polygons(nodes, meshkernel::Projections::spherical);

    meshkernel::Mesh mesh;
    meshkernelapi::MakeGridParametersNative makeGridParametersNative;
    makeGridParametersNative.GridType = 0;
    makeGridParametersNative.GridAngle = 0.0;
    makeGridParametersNative.OriginXCoordinate = 0.0;
    makeGridParametersNative.OriginYCoordinate = 0.0;
    makeGridParametersNative.OriginZCoordinate = 0.0;
    makeGridParametersNative.NumberOfColumns = 3;
    makeGridParametersNative.NumberOfRows = 3;
    makeGridParametersNative.XGridBlockSize = 5000000.0; //resolution in meters (when using spherical coordinates distances are usually much larger)
    makeGridParametersNative.YGridBlockSize = 5000000.0;

    // 2 Execution: function not producing grid points (points gets transformed in meters, therfore everything is outside)
    mesh.MakeMesh(makeGridParametersNative, polygons);
    ASSERT_EQ(0, mesh.GetNumEdges());
    ASSERT_EQ(0, mesh.GetNumNodes());
}

TEST(Mesh, MakeMeshInEmptyPolygonSpherical)
{
    //1 Setup
    std::vector<meshkernel::Point> nodes;
    meshkernel::Polygons polygons(nodes, meshkernel::Projections::spherical);

    meshkernel::Mesh mesh;
    meshkernelapi::MakeGridParametersNative makeGridParametersNative;
    makeGridParametersNative.GridType = 0;
    makeGridParametersNative.GridAngle = 0.0;
    makeGridParametersNative.OriginXCoordinate = 0.0;
    makeGridParametersNative.OriginYCoordinate = 0.0;
    makeGridParametersNative.OriginZCoordinate = 0.0;
    makeGridParametersNative.NumberOfColumns = 3;
    makeGridParametersNative.NumberOfRows = 3;
    makeGridParametersNative.XGridBlockSize = 5000000.0; //resolution in meters (when using spherical coordinates distances are usually much larger)
    makeGridParametersNative.YGridBlockSize = 5000000.0;

    // 2 Execution
    mesh.MakeMesh(makeGridParametersNative, polygons);
    ASSERT_EQ(24, mesh.GetNumEdges());
    ASSERT_EQ(16, mesh.GetNumNodes());

    // x coordinate
    ASSERT_EQ(0.0, mesh.m_nodes[0].x);
    ASSERT_EQ(makeGridParametersNative.XGridBlockSize, mesh.m_nodes[1].x);
    ASSERT_EQ(makeGridParametersNative.XGridBlockSize * 2, mesh.m_nodes[2].x);
    ASSERT_EQ(makeGridParametersNative.XGridBlockSize * 3, mesh.m_nodes[3].x);

    ASSERT_EQ(0.0, mesh.m_nodes[4].x);
    ASSERT_EQ(makeGridParametersNative.XGridBlockSize, mesh.m_nodes[5].x);
    ASSERT_EQ(makeGridParametersNative.XGridBlockSize * 2, mesh.m_nodes[6].x);
    ASSERT_EQ(makeGridParametersNative.XGridBlockSize * 3, mesh.m_nodes[7].x);

    ASSERT_EQ(0.0, mesh.m_nodes[8].x);
    ASSERT_EQ(makeGridParametersNative.XGridBlockSize, mesh.m_nodes[9].x);
    ASSERT_EQ(makeGridParametersNative.XGridBlockSize * 2, mesh.m_nodes[10].x);
    ASSERT_EQ(makeGridParametersNative.XGridBlockSize * 3, mesh.m_nodes[11].x);

    // y coordinate
    ASSERT_EQ(0.0, mesh.m_nodes[0].y);
    ASSERT_EQ(0.0, mesh.m_nodes[1].y);
    ASSERT_EQ(0.0, mesh.m_nodes[2].y);
    ASSERT_EQ(0.0, mesh.m_nodes[3].y);

    ASSERT_EQ(makeGridParametersNative.XGridBlockSize, mesh.m_nodes[4].y);
    ASSERT_EQ(makeGridParametersNative.XGridBlockSize, mesh.m_nodes[5].y);
    ASSERT_EQ(makeGridParametersNative.XGridBlockSize, mesh.m_nodes[6].y);
    ASSERT_EQ(makeGridParametersNative.XGridBlockSize, mesh.m_nodes[7].y);

    ASSERT_EQ(3830222.2156113400, mesh.m_nodes[8].y);
    ASSERT_EQ(3830222.2156113400, mesh.m_nodes[9].y);
    ASSERT_EQ(3830222.2156113400, mesh.m_nodes[10].y);
    ASSERT_EQ(3830222.2156113400, mesh.m_nodes[11].y);
}

TEST(Mesh, TriangulateSamplesWithSkinnyTriangle)
{
    // Prepare
    std::vector<meshkernel::Point> nodes;

    nodes.push_back({302.002502, 472.130371});
    nodes.push_back({144.501526, 253.128174});
    nodes.push_back({368.752930, 112.876755});
    nodes.push_back({707.755005, 358.879242});
    nodes.push_back({301.252502, 471.380371});
    nodes.push_back({302.002502, 472.130371});

    meshkernel::Polygons polygons(nodes, meshkernel::Projections::cartesian);

    // Execute
    std::vector<std::vector<meshkernel::Point>> generatedPoints;
    polygons.CreatePointsInPolygons(generatedPoints);

    meshkernel::Mesh mesh(generatedPoints[0], polygons, meshkernel::Projections::cartesian);

    //// Assert
    constexpr double tolerance = 1e-5;

    ASSERT_EQ(6, mesh.GetNumEdges());

    ASSERT_EQ(4, mesh.m_edges[0].first);
    ASSERT_EQ(1, mesh.m_edges[0].second);

    ASSERT_EQ(1, mesh.m_edges[1].first);
    ASSERT_EQ(2, mesh.m_edges[1].second);

    ASSERT_EQ(2, mesh.m_edges[2].first);
    ASSERT_EQ(4, mesh.m_edges[2].second);

    ASSERT_EQ(0, mesh.m_edges[3].first);
    ASSERT_EQ(2, mesh.m_edges[3].second);

    ASSERT_EQ(2, mesh.m_edges[4].first);
    ASSERT_EQ(3, mesh.m_edges[4].second);

    ASSERT_EQ(3, mesh.m_edges[5].first);
    ASSERT_EQ(0, mesh.m_edges[5].second);
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

    meshkernel::Polygons polygons(nodes, meshkernel::Projections::cartesian);

    // Execute
    std::vector<std::vector<meshkernel::Point>> generatedPoints;
    polygons.CreatePointsInPolygons(generatedPoints);

    meshkernel::Mesh mesh(generatedPoints[0], polygons, meshkernel::Projections::cartesian);
}

TEST(Mesh, TwoTrianglesDuplicatedEdges)
{
    //1 Setup
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

    meshkernel::Mesh mesh;
    // 2 Execution
    mesh.Set(edges, nodes, meshkernel::Projections::cartesian);

    // 3 Validation
    ASSERT_EQ(2, mesh.GetNumFaces());
}

TEST(Mesh, MeshBoundaryToPolygon)
{
    //1 Setup
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

    meshkernel::Mesh mesh;
    mesh.Set(edges, nodes, meshkernel::Projections::cartesian);

    meshkernel::Polygons polygons;

    std::vector<meshkernel::Point> meshBoundaryPolygon;
    int numNodesBoundaryPolygons = 0;
    polygons.MeshBoundaryToPolygon(mesh, meshBoundaryPolygon, numNodesBoundaryPolygons);

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
    //1 Setup
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

    meshkernel::Mesh mesh;
    mesh.Set(edges, nodes, meshkernel::Projections::cartesian);

    ASSERT_EQ(1, mesh.GetNumFaces());
}

TEST(Mesh, NodeMerging)
{
    // 1. Setup
    const int n = 10; // x
    const int m = 10; // y

    std::vector<std::vector<int>> indexesValues(n, std::vector<int>(m));
    std::vector<meshkernel::Point> nodes(n * m);
    std::size_t nodeIndex = 0;
    for (int j = 0; j < m; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            indexesValues[i][j] = i + j * n;
            nodes[nodeIndex] = {(double)i, (double)j};
            nodeIndex++;
        }
    }

    std::vector<meshkernel::Edge> edges((n - 1) * m + (m - 1) * n);
    std::size_t edgeIndex = 0;
    for (int j = 0; j < m; ++j)
    {
        for (int i = 0; i < n - 1; ++i)
        {
            edges[edgeIndex] = {indexesValues[i][j], indexesValues[i + 1][j]};
            edgeIndex++;
        }
    }

    for (int j = 0; j < m - 1; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            edges[edgeIndex] = {indexesValues[i][j + 1], indexesValues[i][j]};
            edgeIndex++;
        }
    }

    meshkernel::Mesh mesh;
    mesh.Set(edges, nodes, meshkernel::Projections::cartesian);

    // Add overlapping nodes
    double generatingDistance = std::sqrt(std::pow(meshkernel::mergingDistance * 0.9, 2) / 2.0);
    std::uniform_real_distribution<double> xDistrution(0.0, generatingDistance);
    std::uniform_real_distribution<double> yDistrution(0.0, generatingDistance);
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());

    nodes.resize(mesh.GetNumNodes() * 2);
    edges.resize(mesh.GetNumEdges() + mesh.GetNumNodes() * 2);
    int originalNodeIndex = 0;
    for (int j = 0; j < m; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            nodes[nodeIndex] = {i + xDistrution(generator), j + yDistrution(generator)};

            // add artificial edges
            auto edge = mesh.m_edges[mesh.m_nodesEdges[originalNodeIndex][0]];
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
    mesh.Set(edges, nodes, meshkernel::Projections::cartesian);

    // 2. Act
    meshkernel::Polygons polygon;
    mesh.MergeNodesInPolygon(polygon);

    // 3. Assert
    ASSERT_EQ(mesh.GetNumNodes(), n * m);
    ASSERT_EQ(mesh.GetNumEdges(), (n - 1) * m + (m - 1) * n);
}

TEST(Mesh, MillionQuads)
{
    const int n = 11; // x
    const int m = 11; // y

    std::vector<std::vector<int>> indexesValues(n, std::vector<int>(m));
    std::vector<meshkernel::Point> nodes(n * m);
    std::size_t nodeIndex = 0;
    for (int j = 0; j < m; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            indexesValues[i][j] = i + j * n;
            nodes[nodeIndex] = {(double)i, (double)j};
            nodeIndex++;
        }
    }

    std::vector<meshkernel::Edge> edges((n - 1) * m + (m - 1) * n);
    std::size_t edgeIndex = 0;
    for (int j = 0; j < m; ++j)
    {
        for (int i = 0; i < n - 1; ++i)
        {
            edges[edgeIndex] = {indexesValues[i][j], indexesValues[i + 1][j]};
            edgeIndex++;
        }
    }

    for (int j = 0; j < m - 1; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            edges[edgeIndex] = {indexesValues[i][j + 1], indexesValues[i][j]};
            edgeIndex++;
        }
    }

    meshkernel::Mesh mesh;
    // now build node-edge mapping
    auto start(std::chrono::steady_clock::now());
    mesh.Set(edges, nodes, meshkernel::Projections::cartesian);
    auto end(std::chrono::steady_clock::now());

    double elapsedTime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    // std::cout << "Elapsed time " << elapsedTime << " s " << std::endl;
    // std::cout << "Number of found cells " << mesh.GetNumFaces() << std::endl;

    EXPECT_LE(elapsedTime, 5.0);
}

TEST(Mesh, InsertNodeInMeshWithExistingNodesRtreeTriggersRTreeReBuild)
{
    // Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projections::cartesian);
    mesh->m_nodesRTree.BuildTree(mesh->m_nodes);

    // insert nodes modifies the number of nodes, m_nodesRTreeRequiresUpdate is set to true
    meshkernel::Point newPoint{10.0, 10.0};
    int newNodeIndex;
    mesh->InsertNode(newPoint, newNodeIndex);

    int newEdgeIndex;
    mesh->ConnectNodes(0, newNodeIndex, newEdgeIndex);

    // when m_nodesRTreeRequiresUpdate = true m_nodesRTree is not empty the mesh.m_nodesRTree is re-build
    mesh->Administrate(meshkernel::Mesh::AdministrationOptions::AdministrateMeshEdgesAndFaces);

    ASSERT_EQ(5, mesh->m_nodesRTree.Size());

    // even if m_edgesRTreeRequiresUpdate = true, m_edgesRTree is initially empty, so it is assumed that is not needed for searches
    ASSERT_EQ(0, mesh->m_edgesRTree.Size());
}

TEST(Mesh, DeleteNodeInMeshWithExistingNodesRtreeTriggersRTreeReBuild)
{
    // Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projections::cartesian);

    meshkernel::Point newPoint{10.0, 10.0};
    mesh->m_nodesRTree.BuildTree(mesh->m_nodes);
    int newNodeIndex;
    mesh->InsertNode(newPoint, newNodeIndex);

    // delete nodes modifies the number of nodes, m_nodesRTreeRequiresUpdate is set to true
    mesh->DeleteNode(0);

    // when m_nodesRTreeRequiresUpdate = true and m_nodesRTree is not empty the mesh.m_nodesRTree is re-build
    mesh->Administrate(meshkernel::Mesh::AdministrationOptions::AdministrateMeshEdgesAndFaces);

    ASSERT_EQ(3, mesh->m_nodesRTree.Size());
}

TEST(Mesh, ConnectNodesInMeshWithExistingEdgesRtreeTriggersRTreeReBuild)
{
    //1 Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projections::cartesian);
    mesh->ComputeEdgesCenters();
    mesh->m_edgesRTree.BuildTree(mesh->m_edgesCenters);

    meshkernel::Point newPoint{10.0, 10.0};
    int newNodeIndex;
    mesh->InsertNode(newPoint, newNodeIndex);

    // connect nodes modifies the number of edges, m_nodesRTreeRequiresUpdate is set to true
    int newEdgeIndex;
    mesh->ConnectNodes(0, newNodeIndex, newEdgeIndex);

    // when m_nodesRTreeRequiresUpdate = true m_nodesRTree is not empty the mesh.m_nodesRTree is re-build
    mesh->Administrate(meshkernel::Mesh::AdministrationOptions::AdministrateMeshEdgesAndFaces);

    // even if m_nodesRTreeRequiresUpdate = true, m_nodesRTree is initially empty, so it is assumed that is not needed for searches
    ASSERT_EQ(0, mesh->m_nodesRTree.Size());

    ASSERT_EQ(5, mesh->m_edgesRTree.Size());
}

TEST(Mesh, DeleteEdgeeInMeshWithExistingEdgesRtreeTriggersRTreeReBuild)
{
    //1 Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projections::cartesian);
    mesh->ComputeEdgesCenters();
    mesh->m_edgesRTree.BuildTree(mesh->m_edgesCenters);

    // DeleteEdge modifies the number of edges, m_edgesRTreeRequiresUpdate is set to true
    mesh->DeleteEdge(0);

    // when m_edgesRTreeRequiresUpdate = true the mesh.m_edgesRTree is re-build with one less edge
    mesh->Administrate(meshkernel::Mesh::AdministrationOptions::AdministrateMeshEdgesAndFaces);

    ASSERT_EQ(3, mesh->m_edgesRTree.Size());
}

TEST(Mesh, GetNodeIndexShouldTriggerNodesRTreeBuild)
{
    //1 Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projections::cartesian);

    // By default, no nodesRTree is build
    ASSERT_EQ(0, mesh->m_nodesRTree.Size());

    // GetNodeIndex builds m_nodesRTree for searching the nodes
    int nodeIndex = mesh->GetNodeIndex({1.5, 1.5}, 10);

    // m_nodesRTree is build
    ASSERT_EQ(4, mesh->m_nodesRTree.Size());

    // m_edgesRTree is not build when searching for nodes
    ASSERT_EQ(0, mesh->m_edgesRTree.Size());
}

TEST(Mesh, FindEdgeCloseToAPointShouldTriggerEdgesRTreeBuild)
{
    //1 Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projections::cartesian);

    // FindEdgeCloseToAPoint builds m_edgesRTree for searching the edges
    const auto edgeIndex = mesh->FindEdgeCloseToAPoint({1.5, 1.5});

    // m_nodesRTree is not build when searching for edges
    ASSERT_EQ(0, mesh->m_nodesRTree.Size());

    // m_edgesRTree is build
    ASSERT_EQ(4, mesh->m_edgesRTree.Size());
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

    meshkernel::Mesh mesh;
    mesh.Set(edges, nodes, meshkernel::Projections::cartesian);

    // execute, only one obtuse triangle should be found
    const auto obtuseTrianglesCount = mesh.GetObtuseTrianglesCenters().size();

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

    meshkernel::Mesh mesh;
    mesh.Set(edges, nodes, meshkernel::Projections::cartesian);

    // execute, by setting the smallFlowEdgesThreshold high, a small flow edge will be found
    const auto numSmallFlowEdgeFirstQuery = mesh.GetEdgesCrossingSmallFlowEdges(100).size();

    // execute, by setting the smallFlowEdgesThreshold low, no small flow edge will be found
    const auto numSmallFlowEdgeSecondQuery = mesh.GetEdgesCrossingSmallFlowEdges(0.0).size();

    // assert a small flow edge is found
    ASSERT_EQ(1, numSmallFlowEdgeFirstQuery);
    ASSERT_EQ(0, numSmallFlowEdgeSecondQuery);
}

TEST(Mesh, RemoveSmallFlowEdge)
{
    // Setup a mesh with eight triangles
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/RemoveSmallFlowEdgesTests/remove_small_flow_edges_net.nc");

    ASSERT_EQ(8, mesh->GetNumFaces());

    // After merging the number of faces is reduced
    mesh->RemoveSmallFlowEdges(1.0);

    ASSERT_EQ(3, mesh->GetNumFaces());
}

TEST(Mesh, RemoveSmallTrianglesAtBoundaries)
{
    // Setup a mesh with two triangles
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/RemoveSmallFlowEdgesTests/remove_small_flow_edges_quad_net.nc");

    ASSERT_EQ(2, mesh->GetNumFaces());

    // After merging
    mesh->RemoveSmallTrianglesAtBoundaries(0.6);

    ASSERT_EQ(1, mesh->GetNumFaces());

    const double tolerance = 1e-8;
    ASSERT_NEAR(364.17013549804688, mesh->m_nodes[0].x, tolerance);
    ASSERT_NEAR(295.21142578125000, mesh->m_nodes[1].x, tolerance);
    ASSERT_NEAR(421.46209716796875, mesh->m_nodes[2].x, tolerance);
    ASSERT_NEAR(359.79510498046875, mesh->m_nodes[3].x, tolerance);

    ASSERT_NEAR(374.00662231445313, mesh->m_nodes[0].y, tolerance);
    ASSERT_NEAR(300.48181152343750, mesh->m_nodes[1].y, tolerance);
    ASSERT_NEAR(295.33038330078125, mesh->m_nodes[2].y, tolerance);
    ASSERT_NEAR(398.59295654296875, mesh->m_nodes[3].y, tolerance);
}
