#include "../Mesh.cpp"
#include "../Entities.hpp"
#include "../Polygons.hpp"
#include <gtest/gtest.h>
#include <chrono>

TEST(TestMesh, OneQuadTestConstructor) 
{
    //1 Setup
    std::vector<GridGeom::Point> nodes;
    nodes.push_back(GridGeom::Point{ 0.0,0.0 });
    nodes.push_back(GridGeom::Point{ 0.0,10.0 });
    nodes.push_back(GridGeom::Point{ 10.0,0.0 });
    nodes.push_back(GridGeom::Point{ 10.0,10.0 });
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
    EXPECT_EQ(0, mesh.m_nodesEdges[0][0]);
    EXPECT_EQ(2, mesh.m_nodesEdges[0][1]);

    EXPECT_EQ(1, mesh.m_nodesEdges[1][0]);
    EXPECT_EQ(2, mesh.m_nodesEdges[1][1]);

    EXPECT_EQ(0, mesh.m_nodesEdges[2][0]);
    EXPECT_EQ(3, mesh.m_nodesEdges[2][1]);

    EXPECT_EQ(1, mesh.m_nodesEdges[3][0]);
    EXPECT_EQ(3, mesh.m_nodesEdges[3][1]);

    // each node has two edges int this case
    EXPECT_EQ(2, mesh.m_nodesNumEdges[0]);
    EXPECT_EQ(2, mesh.m_nodesNumEdges[1]);
    EXPECT_EQ(2, mesh.m_nodesNumEdges[2]);
    EXPECT_EQ(2, mesh.m_nodesNumEdges[3]);

    // the nodes composing the face, in ccw order
    EXPECT_EQ(0, mesh.m_facesNodes[0][0]);
    EXPECT_EQ(2, mesh.m_facesNodes[0][1]);
    EXPECT_EQ(3, mesh.m_facesNodes[0][2]);
    EXPECT_EQ(1, mesh.m_facesNodes[0][3]);

    // the edges composing the face, in ccw order
    EXPECT_EQ(0, mesh.m_facesEdges[0][0]);
    EXPECT_EQ(3, mesh.m_facesEdges[0][1]);
    EXPECT_EQ(1, mesh.m_facesEdges[0][2]);
    EXPECT_EQ(2, mesh.m_facesEdges[0][3]);

    // the found circumcenter for the face
    EXPECT_DOUBLE_EQ(5.0, mesh.m_facesCircumcenters[0].x);
    EXPECT_DOUBLE_EQ(5.0, mesh.m_facesCircumcenters[0].y);

    // each edge has only one face in this case
    EXPECT_EQ(1, mesh.m_edgesNumFaces[0]);
    EXPECT_EQ(1, mesh.m_edgesNumFaces[1]);
    EXPECT_EQ(1, mesh.m_edgesNumFaces[2]);
    EXPECT_EQ(1, mesh.m_edgesNumFaces[3]);

    //each edge is a boundary edge, so the second entry of edgesFaces is an invalid index (in C++, -1)
    EXPECT_EQ(-1, mesh.m_edgesFaces[0][1]);
    EXPECT_EQ(-1, mesh.m_edgesFaces[1][1]);
    EXPECT_EQ(-1, mesh.m_edgesFaces[2][1]);
    EXPECT_EQ(-1, mesh.m_edgesFaces[3][1]);
}

TEST(TestMesh, MakeMeshInPolygon)
{
    //1 Setup
    GridGeom::Polygons polygons;
    std::vector<GridGeom::Point> nodes;

    nodes.push_back(GridGeom::Point{ 302.002502,472.130371 });
    nodes.push_back(GridGeom::Point{ 144.501526, 253.128174 });
    nodes.push_back(GridGeom::Point{ 368.752930, 112.876755 });
    nodes.push_back(GridGeom::Point{ 707.755005, 358.879242 });
    nodes.push_back(GridGeom::Point{ 301.252502, 471.380371 });
    nodes.push_back(GridGeom::Point{ 302.002502, 472.130371 });

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
    EXPECT_EQ(17, mesh.m_numFaces);

}

TEST(TestMesh, TriangulateSamplesWithSkinnyTriangle)
{
    // Prepare
    GridGeom::Polygons polygons;
    std::vector<GridGeom::Point> nodes;

    nodes.push_back(GridGeom::Point{ 302.002502,472.130371 });
    nodes.push_back(GridGeom::Point{ 144.501526, 253.128174 });
    nodes.push_back(GridGeom::Point{ 368.752930, 112.876755 });
    nodes.push_back(GridGeom::Point{ 707.755005, 358.879242 });
    nodes.push_back(GridGeom::Point{ 301.252502, 471.380371 });
    nodes.push_back(GridGeom::Point{ 302.002502, 472.130371 });

    polygons.Set(nodes, GridGeom::Projections::cartesian);

    // Execute
    std::vector<std::vector<GridGeom::Point>> generatedPoints;
    bool success = polygons.CreatePointsInPolygons(generatedPoints);
    ASSERT_TRUE(success);

    GridGeom::Mesh mesh(generatedPoints[0], polygons, GridGeom::Projections::cartesian);

    //// Assert
    constexpr double tolerance = 1e-5;

    EXPECT_EQ(6, mesh.m_edges.size());

    EXPECT_EQ(4, mesh.m_edges[0].first);
    EXPECT_EQ(1, mesh.m_edges[0].second);

    EXPECT_EQ(1, mesh.m_edges[1].first);
    EXPECT_EQ(2, mesh.m_edges[1].second);

    EXPECT_EQ(2, mesh.m_edges[2].first);
    EXPECT_EQ(4, mesh.m_edges[2].second);

    EXPECT_EQ(0, mesh.m_edges[3].first);
    EXPECT_EQ(2, mesh.m_edges[3].second);

    EXPECT_EQ(2, mesh.m_edges[4].first);
    EXPECT_EQ(3, mesh.m_edges[4].second);

    EXPECT_EQ(3, mesh.m_edges[5].first);
    EXPECT_EQ(0, mesh.m_edges[5].second);

}

TEST(TestMesh, TriangulateSamples)
{
    // Prepare
    GridGeom::Polygons polygons;
    std::vector<GridGeom::Point> nodes;

    nodes.push_back(GridGeom::Point{ 498.503152894023, 1645.82297461613 });
    nodes.push_back(GridGeom::Point{ -5.90937355559299, 814.854361678898 });
    nodes.push_back(GridGeom::Point{ 851.30035347439, 150.079471329115 });
    nodes.push_back(GridGeom::Point{ 1411.11078745316, 1182.22995897746 });
    nodes.push_back(GridGeom::Point{ 501.418832237663, 1642.90729527249 });
    nodes.push_back(GridGeom::Point{ 498.503152894023, 1645.82297461613 });

    polygons.Set(nodes, GridGeom::Projections::cartesian);

    // Execute
    std::vector<std::vector<GridGeom::Point>> generatedPoints;
    bool success = polygons.CreatePointsInPolygons(generatedPoints);
    ASSERT_TRUE(success);

    GridGeom::Mesh mesh(generatedPoints[0], polygons, GridGeom::Projections::cartesian);
}


TEST(TestMesh, TwoTrianglesDuplicatedEdges)
{
    //1 Setup
    std::vector<GridGeom::Point> nodes;
    nodes.push_back(GridGeom::Point{ -5.90937355559299, 814.854361678898 });
    nodes.push_back(GridGeom::Point{ 851.30035347439, 150.079471329115 });
    nodes.push_back(GridGeom::Point{ 1411.11078745316, 1182.22995897746 });
    nodes.push_back(GridGeom::Point{ 501.418832237663, 1642.90729527249 });
    nodes.push_back(GridGeom::Point{ 498.503152894023, 1645.82297461613 });
    std::vector<GridGeom::Edge> edges;
    edges.push_back({ 0, 3 });
    edges.push_back({ 0, 2 });
    edges.push_back({ 2, 3 });
    edges.push_back({ 0, 1 });
    edges.push_back({ 1, 2 });
    edges.push_back({ 3, 0 });


    GridGeom::Mesh mesh;

    // 2 Execution
    mesh.Set(edges, nodes, GridGeom::Projections::cartesian);

    // 3 Validation
}


TEST(PerformanceTest, MillionQuads)
{
    const int n = 11; //x
    const int m = 11; //y

    //std::cout << "start adding edges " << std::endl;
    auto start(std::chrono::steady_clock::now());

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
    auto end(std::chrono::steady_clock::now());
    //std::cout << "Elapsed time " << std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count() << " s " << std::endl;

    //std::cout << "start finding cells " << std::endl;
    start = std::chrono::steady_clock::now();
    // now build node-edge mapping
    GridGeom::Mesh mesh;
    mesh.Set(edges, nodes, GridGeom::Projections::cartesian);

    end = std::chrono::steady_clock::now();

    double elapsedTime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    std::cout << "Elapsed time " << elapsedTime << " s " << std::endl;

    // the number of found faces is
    //auto faces = mesh.m_facesNodes;
    std::cout << "Number of found cells " << mesh.m_facesNodes.size() << std::endl;
    //std::cout << "First face " << faces[0][0] << " " << faces[0][1] << " " << faces[0][2] << " " << faces[0][3] << std::endl;
    //std::cout << "Second face " << faces[1][0] << " " << faces[1][1] << " " << faces[1][2] << " " << faces[1][3] << std::endl;

    // to beat fortran interactor, we need to perform the entire administration in less than 1.5 seconds
    EXPECT_LE(elapsedTime, 1.5);
}
