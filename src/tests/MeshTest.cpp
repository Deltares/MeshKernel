#include "../Mesh.cpp"
#include "../Entities.hpp"
#include "../Polygons.cpp"
#include "../Constants.cpp"
#include <gtest/gtest.h>
#include <chrono>
#include <random>

TEST(Mesh, OneQuadTestConstructor) 
{
    //1 Setup
    std::vector<GridGeom::Point> nodes;
    nodes.push_back({ 0.0,0.0 });
    nodes.push_back({ 0.0,10.0 });
    nodes.push_back({ 10.0,0.0 });
    nodes.push_back({ 10.0,10.0 });
    std::vector<GridGeom::Edge> edges;
    edges.push_back({ 0, 2 });
    edges.push_back({ 1, 3 });
    edges.push_back({ 0, 1 });
    edges.push_back({ 2, 3 });
    
    GridGeom::Mesh mesh;

    // 2 Execution
    mesh.Set(edges, nodes, GridGeom::Projections::cartesian);

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
    GridGeom::Polygons polygons;
    std::vector<GridGeom::Point> nodes;

    nodes.push_back({ 302.002502,472.130371 });
    nodes.push_back({ 144.501526, 253.128174 });
    nodes.push_back({ 368.752930, 112.876755 });
    nodes.push_back({ 707.755005, 358.879242 });
    nodes.push_back({ 301.252502, 471.380371 });
    nodes.push_back({ 302.002502, 472.130371 });

    polygons.Set(nodes, GridGeom::Projections::cartesian);
    
    GridGeom::Mesh mesh;
    GridGeomApi::MakeGridParametersNative makeGridParametersNative;
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
    ASSERT_EQ(17, mesh.GetNumFaces());
}

TEST(Mesh, TriangulateSamplesWithSkinnyTriangle)
{
    // Prepare
    GridGeom::Polygons polygons;
    std::vector<GridGeom::Point> nodes;

    nodes.push_back({ 302.002502,472.130371 });
    nodes.push_back({ 144.501526, 253.128174 });
    nodes.push_back({ 368.752930, 112.876755 });
    nodes.push_back({ 707.755005, 358.879242 });
    nodes.push_back({ 301.252502, 471.380371 });
    nodes.push_back({ 302.002502, 472.130371 });

    polygons.Set(nodes, GridGeom::Projections::cartesian);

    // Execute
    std::vector<std::vector<GridGeom::Point>> generatedPoints;
    bool success = polygons.CreatePointsInPolygons(generatedPoints);
    ASSERT_TRUE(success);

    GridGeom::Mesh mesh(generatedPoints[0], polygons, GridGeom::Projections::cartesian);

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
    GridGeom::Polygons polygons;
    std::vector<GridGeom::Point> nodes;

    nodes.push_back({ 498.503152894023, 1645.82297461613 });
    nodes.push_back({ -5.90937355559299, 814.854361678898 });
    nodes.push_back({ 851.30035347439, 150.079471329115 });
    nodes.push_back({ 1411.11078745316, 1182.22995897746 });
    nodes.push_back({ 501.418832237663, 1642.90729527249 });
    nodes.push_back({ 498.503152894023, 1645.82297461613 });

    polygons.Set(nodes, GridGeom::Projections::cartesian);

    // Execute
    std::vector<std::vector<GridGeom::Point>> generatedPoints;
    bool success = polygons.CreatePointsInPolygons(generatedPoints);
    ASSERT_TRUE(success);

    GridGeom::Mesh mesh(generatedPoints[0], polygons, GridGeom::Projections::cartesian);
}


TEST(Mesh, TwoTrianglesDuplicatedEdges)
{
    //1 Setup
    std::vector<GridGeom::Point> nodes;
    nodes.push_back({ 0.0, 0.0 });
    nodes.push_back({ 5.0, -5.0 });
    nodes.push_back({ 10.0, 0.0 });
    nodes.push_back({ 5.0, 5.0 });
    std::vector<GridGeom::Edge> edges;
    edges.push_back({ 0, 3 });
    edges.push_back({ 0, 2 });
    edges.push_back({ 2, 3 });
    edges.push_back({ 0, 1 });
    edges.push_back({ 2, 1 });

    GridGeom::Mesh mesh;
    // 2 Execution
    mesh.Set(edges, nodes, GridGeom::Projections::cartesian);

    // 3 Validation
    ASSERT_EQ(2, mesh.GetNumFaces());
}

TEST(Mesh, MeshBoundaryToPolygon)
{
    //1 Setup
    std::vector<GridGeom::Point> nodes;
    nodes.push_back({ 0.0, 0.0 });
    nodes.push_back({ 5.0, -5.0 });
    nodes.push_back({ 10.0, 0.0 });
    nodes.push_back({ 5.0, 5.0 });
    std::vector<GridGeom::Edge> edges;
    edges.push_back({ 0, 3 });
    edges.push_back({ 0, 2 });
    edges.push_back({ 2, 3 });
    edges.push_back({ 0, 1 });
    edges.push_back({ 2, 1 });

    GridGeom::Mesh mesh;
    mesh.Set(edges, nodes, GridGeom::Projections::cartesian);

    GridGeom::Polygons polygons;

    int counterClockWise = 0;
    int setMeshState = 0;
    std::vector<GridGeom::Point> meshBoundaryPolygon;
    int numNodesBoundaryPolygons = 0;
    polygons.MeshBoundaryToPolygon(mesh, counterClockWise, setMeshState, meshBoundaryPolygon, numNodesBoundaryPolygons);

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
    std::vector<GridGeom::Point> nodes;
    nodes.push_back({ 0.0, 0.0 });
    nodes.push_back({ 5.0, 0.0 });
    nodes.push_back({ 3.0, 2.0 });
    nodes.push_back({ 3.0, 4.0 });

    std::vector<GridGeom::Edge> edges;
    edges.push_back({ 0, 1 });
    edges.push_back({ 1, 3 });
    edges.push_back({ 3, 0 });
    edges.push_back({ 2, 1 });


    GridGeom::Mesh mesh;
    mesh.Set(edges, nodes, GridGeom::Projections::cartesian);

    ASSERT_EQ(1, mesh.GetNumFaces());
}

TEST(Mesh, NodeMerging)
{
    // 1. Setup
    const int n = 10; // x
    const int m = 10; // y

    std::vector<std::vector<int>> indexesValues(n, std::vector<int>(m));
    std::vector<GridGeom::Point> nodes(n * m);
    std::size_t nodeIndex = 0;
    for (int j = 0; j < m; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            indexesValues[i][j] = i + j * n;
            nodes[nodeIndex] = { (double)i, (double)j };
            nodeIndex++;
        }
    }

    std::vector<GridGeom::Edge> edges((n - 1) * m + (m - 1) * n);
    std::size_t edgeIndex = 0;
    for (int j = 0; j < m; ++j)
    {
        for (int i = 0; i < n - 1; ++i)
        {
            edges[edgeIndex] = { indexesValues[i][j], indexesValues[i + 1][j] };
            edgeIndex++;
        }
    }

    for (int j = 0; j < m - 1; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            edges[edgeIndex] = { indexesValues[i][j + 1], indexesValues[i][j] };
            edgeIndex++;
        }
    }

    GridGeom::Mesh mesh;
    mesh.Set(edges, nodes, GridGeom::Projections::cartesian);

    // Add overlapping nodes
    double generatingDistance = std::sqrt(std::pow(GridGeom::mergingDistance*0.9, 2) / 2.0);
    std::uniform_real_distribution<double>  xDistrution(0.0, generatingDistance);
    std::uniform_real_distribution<double>  yDistrution(0.0, generatingDistance);
    std::random_device                      rand_dev;
    std::mt19937                            generator(rand_dev());
    
    nodes.resize(nodes.size()*2);
    edges.resize(edges.size()+ nodes.size() * 2);
    int originalNodeIndex = 0;
    for (int j = 0; j < m; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            nodes[nodeIndex] = { i + xDistrution(generator), j + yDistrution(generator)};
            
            // add artificial edges
            auto edge = mesh.m_edges[mesh.m_nodesEdges[originalNodeIndex][0]];
            auto otherNode = edge.first + edge.second - originalNodeIndex;

            edges[edgeIndex] = { nodeIndex, otherNode };
            edgeIndex++;
            edges[edgeIndex] = { nodeIndex, originalNodeIndex };
            edgeIndex++;

            nodeIndex++;
            originalNodeIndex++;
        }
    }

    nodes.resize(nodeIndex);
    edges.resize(edgeIndex);

    // re set with augmented nodes
    mesh.Set(edges, nodes, GridGeom::Projections::cartesian);

    // 2. Act
    GridGeom::Polygons polygon;
    //auto start(std::chrono::steady_clock::now());
    mesh.MergeNodesInPolygon(polygon);
    //auto end(std::chrono::steady_clock::now());

    //double elapsedTime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    //std::cout << "Elapsed time NodeMerging " << elapsedTime << " s " << std::endl;

    // 3. Assert
    ASSERT_EQ(mesh.GetNumFaces(), (n-1)*(m-1));
}

TEST(Mesh, MillionQuads)
{
    const int n = 3; // x
    const int m = 3; // y

    std::vector<std::vector<int>> indexesValues(n, std::vector<int>(m));
    std::vector<GridGeom::Point> nodes(n * m);
    std::size_t nodeIndex = 0;
    for (int j = 0; j < m; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            indexesValues[i][j] = i + j * n;
            nodes[nodeIndex] = { (double)i, (double)j };
            nodeIndex++;
        }
    }

    std::vector<GridGeom::Edge> edges((n - 1) * m + (m - 1) * n);
    std::size_t edgeIndex = 0;
    for (int j = 0; j < m; ++j)
    {
        for (int i = 0; i < n - 1; ++i)
        {
            edges[edgeIndex] = { indexesValues[i][j], indexesValues[i + 1][j] };
            edgeIndex++;
        }
    }

    for (int j = 0; j < m - 1; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            edges[edgeIndex] = { indexesValues[i][j + 1], indexesValues[i][j] };
            edgeIndex++;
        }
    }

    //std::cout << "Elapsed time " << std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count() << " s " << std::endl;

    //std::cout << "start finding cells " << std::endl;
    auto start(std::chrono::steady_clock::now());
    // now build node-edge mapping
    GridGeom::Mesh mesh;
    mesh.Set(edges, nodes, GridGeom::Projections::cartesian);

    auto end(std::chrono::steady_clock::now());

    double elapsedTime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    std::cout << "Elapsed time " << elapsedTime << " s " << std::endl;

    // the number of found faces is
    //auto faces = mesh.m_facesNodes;
    std::cout << "Number of found cells " << mesh.GetNumFaces() << std::endl;
    //std::cout << "First face " << faces[0][0] << " " << faces[0][1] << " " << faces[0][2] << " " << faces[0][3] << std::endl;
    //std::cout << "Second face " << faces[1][0] << " " << faces[1][1] << " " << faces[1][2] << " " << faces[1][3] << std::endl;

    // to beat fortran interactor, we need to perform the entire administration in less than 1.5 seconds
    EXPECT_LE(elapsedTime, 2.0);
}
