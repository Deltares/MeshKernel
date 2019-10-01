#include "GridGeomTest.hpp"
#include "../Orthogonalization.cpp"


TEST(OrthogonalizationTest, TestOrthogonalizationOneQuadOneTriangle)
{
    using Mesh = Mesh<cartesianPoint>;

    //One gets the edges
    std::vector<cartesianPoint> nodes;
    nodes.push_back(cartesianPoint{ 0.0,0.0 });
    nodes.push_back(cartesianPoint{ 0.0,10.0 });
    nodes.push_back(cartesianPoint{ 10.0,0.0 });
    nodes.push_back(cartesianPoint{ 10.0,10.0 });
    nodes.push_back(cartesianPoint{ 20.0,0.0 });

    std::vector<Edge> edges;
    // Local edges
    edges.push_back({ 1, 0 });
    edges.push_back({ 0, 2 });
    edges.push_back({ 2, 4 });
    edges.push_back({ 4, 3 });
    edges.push_back({ 3, 2 });
    edges.push_back({ 3, 1 });

    // now build node-edge mapping
    Mesh mesh;
    mesh.setState(edges, nodes);
    Orthogonalization<Mesh> orthogonalization;

    std::vector<double> aspectRatio;
    orthogonalization.initialize(mesh);

    constexpr double tolerance = 1e-8;
    constexpr double largeTolerance = 10.0;

    //centers of mass
    ASSERT_NEAR(13.333333333, mesh.m_facesMassCenters[0].x, tolerance);
    ASSERT_NEAR(3.333333333, mesh.m_facesMassCenters[0].y, tolerance);

    ASSERT_NEAR(5.0, mesh.m_facesCircumcenters[1].x, tolerance);
    ASSERT_NEAR(5.0, mesh.m_facesCircumcenters[1].y, tolerance);

    //aspect ratios (not correct yet)
    ASSERT_NEAR(1.0, orthogonalization.m_aspectRatios[0], tolerance);
    ASSERT_NEAR(1.0, orthogonalization.m_aspectRatios[1], tolerance);
    ASSERT_NEAR(1.0, orthogonalization.m_aspectRatios[2], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_aspectRatios[3], tolerance);
    ASSERT_NEAR(1.0, orthogonalization.m_aspectRatios[4], tolerance);
    ASSERT_NEAR(1.0, orthogonalization.m_aspectRatios[5], tolerance);

    //node 0
    EXPECT_EQ(2, orthogonalization.m_numTopologyFaces[0]);
    EXPECT_EQ(4, orthogonalization.m_numTopologyNodes[0]);

    ASSERT_NEAR(0.0, orthogonalization.m_topologyXi[0][0], tolerance);
    ASSERT_NEAR(1.0, orthogonalization.m_topologyXi[0][1], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_topologyXi[0][2], tolerance);
    ASSERT_NEAR(1.0, orthogonalization.m_topologyXi[0][3], tolerance);

    ASSERT_NEAR(0.0, orthogonalization.m_topologyEta[0][0], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_topologyEta[0][1], tolerance);
    ASSERT_NEAR(1.0, orthogonalization.m_topologyEta[0][2], tolerance);
    ASSERT_NEAR(1.0, orthogonalization.m_topologyEta[0][3], tolerance);

    //node 1
    EXPECT_EQ(2, orthogonalization.m_numTopologyFaces[1]);
    EXPECT_EQ(4, orthogonalization.m_numTopologyNodes[1]);

    ASSERT_NEAR(0.0, orthogonalization.m_topologyXi[1][0], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_topologyXi[1][1], tolerance);
    ASSERT_NEAR(1.0, orthogonalization.m_topologyXi[1][2], tolerance);
    ASSERT_NEAR(1.0, orthogonalization.m_topologyXi[1][3], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_topologyXi[1][4], tolerance);

    ASSERT_NEAR(0.0, orthogonalization.m_topologyEta[1][0], tolerance);
    ASSERT_NEAR(1.0, orthogonalization.m_topologyEta[1][1], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_topologyEta[1][2], tolerance);
    ASSERT_NEAR(1.0, orthogonalization.m_topologyEta[1][3], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_topologyEta[1][4], tolerance);

    //node 2
    EXPECT_EQ(3, orthogonalization.m_numTopologyFaces[2]);
    EXPECT_EQ(5, orthogonalization.m_numTopologyNodes[2]);

    ASSERT_NEAR(0.0, orthogonalization.m_topologyXi[2][0], tolerance);
    ASSERT_NEAR(1.0, orthogonalization.m_topologyXi[2][1], tolerance);
    ASSERT_NEAR(-1.0, orthogonalization.m_topologyXi[2][2], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_topologyXi[2][3], tolerance);
    ASSERT_NEAR(1.0, orthogonalization.m_topologyXi[2][4], tolerance);

    ASSERT_NEAR(0.0, orthogonalization.m_topologyEta[2][0], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_topologyEta[2][1], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_topologyEta[2][2], tolerance);
    ASSERT_NEAR(-1.0, orthogonalization.m_topologyEta[2][3], tolerance);
    ASSERT_NEAR(-1.0, orthogonalization.m_topologyEta[2][4], tolerance);


    //node 3
    EXPECT_EQ(3, orthogonalization.m_numTopologyFaces[3]);
    EXPECT_EQ(5, orthogonalization.m_numTopologyNodes[3]);

    ASSERT_NEAR(0.0, orthogonalization.m_topologyXi[3][0], tolerance);
    ASSERT_NEAR(1.0, orthogonalization.m_topologyXi[3][1], tolerance);
    ASSERT_NEAR(-1.0, orthogonalization.m_topologyXi[3][2], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_topologyXi[3][3], tolerance);
    ASSERT_NEAR(-1.0, orthogonalization.m_topologyXi[3][4], tolerance);

    ASSERT_NEAR(0.0, orthogonalization.m_topologyEta[3][0], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_topologyEta[3][1], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_topologyEta[3][2], tolerance);
    ASSERT_NEAR(-1.0, orthogonalization.m_topologyEta[3][3], tolerance);
    ASSERT_NEAR(-1.0, orthogonalization.m_topologyEta[3][4], tolerance);

    //node 4
    EXPECT_EQ(2, orthogonalization.m_numTopologyFaces[4]);
    EXPECT_EQ(3, orthogonalization.m_numTopologyNodes[4]);

    ASSERT_NEAR(0.0, orthogonalization.m_topologyXi[4][0], tolerance);
    ASSERT_NEAR(1.0, orthogonalization.m_topologyXi[4][1], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_topologyXi[4][2], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_topologyXi[4][3], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_topologyXi[4][4], tolerance);

    ASSERT_NEAR(0.0, orthogonalization.m_topologyEta[4][0], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_topologyEta[4][1], tolerance);
    ASSERT_NEAR(1.0, orthogonalization.m_topologyEta[4][2], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_topologyEta[4][3], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_topologyEta[4][4], tolerance);
    
    // first topology
    ASSERT_NEAR(0.0, orthogonalization.m_Az[0][0][0], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Az[0][0][1], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Az[0][0][2], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Az[0][0][3], tolerance);

    ASSERT_NEAR(0.0, orthogonalization.m_Az[0][1][0], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Az[0][1][1], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Az[0][1][2], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Az[0][1][3], tolerance);

    ASSERT_NEAR(0.0, orthogonalization.m_Gxi[0][0][0], largeTolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Gxi[0][0][1], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Gxi[0][0][2], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Gxi[0][0][3], tolerance);

    ASSERT_NEAR(-2e16, orthogonalization.m_Gxi[0][1][0], largeTolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Gxi[0][1][1], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Gxi[0][1][2], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Gxi[0][1][3], tolerance);

    ASSERT_NEAR(-2e16, orthogonalization.m_Geta[0][0][0], largeTolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Geta[0][0][1], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Geta[0][0][2], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Geta[0][0][3], tolerance);

    ASSERT_NEAR(0.0, orthogonalization.m_Geta[0][1][0], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Geta[0][1][1], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Geta[0][1][2], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Geta[0][1][3], tolerance);

    ASSERT_NEAR(0.0, orthogonalization.m_Divxi[0][0], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Divxi[0][1], tolerance);

    ASSERT_NEAR(0.0, orthogonalization.m_Diveta[0][0], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Diveta[0][1], tolerance);

    ASSERT_NEAR(0.0, orthogonalization.m_Jxi[0][0], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Jxi[0][1], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Jxi[0][2], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Jxi[0][3], tolerance);

    ASSERT_NEAR(0.0, orthogonalization.m_Jeta[0][0], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Jeta[0][1], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Jeta[0][2], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Jeta[0][3], tolerance);

    ASSERT_NEAR(0.0, orthogonalization.m_ww2[0][0], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_ww2[0][1], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_ww2[0][2], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_ww2[0][3], tolerance);

    // second topology
    ASSERT_NEAR(0.0, orthogonalization.m_Az[1][0][0], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Az[1][0][1], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Az[1][0][2], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Az[1][0][3], tolerance);

    ASSERT_NEAR(0.0, orthogonalization.m_Az[1][1][0], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Az[1][1][1], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Az[1][1][2], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Az[1][1][3], tolerance);

    ASSERT_NEAR(2e16, orthogonalization.m_Gxi[1][0][0], largeTolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Gxi[1][0][1], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Gxi[1][0][2], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Gxi[1][0][3], tolerance);

    ASSERT_NEAR(0.0, orthogonalization.m_Gxi[1][1][0], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Gxi[1][1][1], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Gxi[1][1][2], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Gxi[1][1][3], tolerance);

    ASSERT_NEAR(0.0, orthogonalization.m_Geta[1][0][0], largeTolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Geta[1][0][1], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Geta[1][0][2], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Geta[1][0][3], tolerance);

    ASSERT_NEAR(2e16, orthogonalization.m_Geta[1][1][0], largeTolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Geta[1][1][1], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Geta[1][1][2], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Geta[1][1][3], tolerance);

    ASSERT_NEAR(0.0, orthogonalization.m_Divxi[1][0], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Divxi[1][1], tolerance);

    ASSERT_NEAR(0.0, orthogonalization.m_Diveta[1][0], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Diveta[1][1], tolerance);

    ASSERT_NEAR(0.0, orthogonalization.m_Jxi[1][0], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Jxi[1][1], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Jxi[1][2], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Jxi[1][3], tolerance);

    ASSERT_NEAR(0.0, orthogonalization.m_Jeta[1][0], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Jeta[1][1], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Jeta[1][2], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_Jeta[1][3], tolerance);

    ASSERT_NEAR(0.0, orthogonalization.m_ww2[1][0], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_ww2[1][1], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_ww2[1][2], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_ww2[1][3], tolerance);

    orthogonalization.iterate(mesh);

    //check nodes are not moved after orthogonalization for thi specific case
    ASSERT_NEAR(0.0, mesh.m_nodes[0].x, tolerance);
    ASSERT_NEAR(0.0, mesh.m_nodes[1].x, tolerance);
    ASSERT_NEAR(10.0, mesh.m_nodes[2].x, tolerance);
    ASSERT_NEAR(10.0, mesh.m_nodes[3].x, tolerance);
    ASSERT_NEAR(20.0, mesh.m_nodes[4].x, tolerance);

    ASSERT_NEAR(0.0, mesh.m_nodes[0].y, tolerance);
    ASSERT_NEAR(10.0, mesh.m_nodes[1].y, tolerance);
    ASSERT_NEAR(0.0, mesh.m_nodes[2].y, tolerance);
    ASSERT_NEAR(10.0, mesh.m_nodes[3].y, tolerance);
    ASSERT_NEAR(0.0, mesh.m_nodes[4].y, tolerance);
}


TEST(OrthogonalizationTest, TestOrthogonalizationFunctions)
{
    using Mesh = Mesh<cartesianPoint>;

    //One gets the edges
    std::vector<cartesianPoint> nodes;
    nodes.push_back(cartesianPoint{ 0.0,0.0 });
    nodes.push_back(cartesianPoint{ 0.0,10.0 });
    nodes.push_back(cartesianPoint{ 10.0,0.0 });
    nodes.push_back(cartesianPoint{ 10.0,10.0 });

    std::vector<Edge> edges;
    // Local edges
    edges.push_back({ 0, 2 });
    edges.push_back({ 1, 3 });
    edges.push_back({ 0, 1 });
    edges.push_back({ 2, 3 });

    // now build node-edge mapping
    Mesh mesh;
    mesh.setState(edges, nodes);
    Orthogonalization<Mesh> orthogonalization;

    std::vector<double> aspectRatio;
    orthogonalization.initialize(mesh);

    EXPECT_EQ(1, orthogonalization.m_aspectRatios[0]);
    EXPECT_EQ(1, orthogonalization.m_aspectRatios[1]);
    EXPECT_EQ(1, orthogonalization.m_aspectRatios[2]);
    EXPECT_EQ(1, orthogonalization.m_aspectRatios[3]);
}

TEST(OrthogonalizationTest, TestOrthogonalizationTriangularGrid)
{
    using Mesh = Mesh<cartesianPoint>;

    //One gets the edges
    std::vector<cartesianPoint> nodes;

    nodes.push_back(cartesianPoint{ 322.252624511719,454.880187988281 });
    nodes.push_back(cartesianPoint{ 227.002044677734,360.379241943359 });
    nodes.push_back(cartesianPoint{ 259.252227783203,241.878051757813 });
    nodes.push_back(cartesianPoint{ 428.003295898438,210.377746582031 });
    nodes.push_back(cartesianPoint{ 536.003967285156,310.878753662109 });
    nodes.push_back(cartesianPoint{ 503.753784179688,432.379974365234 });
    nodes.push_back(cartesianPoint{ 350.752807617188,458.630249023438 });
    nodes.push_back(cartesianPoint{ 343.15053976393,406.232256102912 });
    nodes.push_back(cartesianPoint{ 310.300984548069,319.41005739802 });
    nodes.push_back(cartesianPoint{ 423.569603308318,326.17986967523 });

    std::vector<Edge> edges;
    // Local edges
    edges.push_back({ 3, 9 });
    edges.push_back({ 9, 2 });
    edges.push_back({ 2, 3 });
    edges.push_back({ 3, 4 });
    edges.push_back({ 4, 9 });
    edges.push_back({ 2, 8 });
    edges.push_back({ 8, 1 });
    edges.push_back({ 1, 2 });
    edges.push_back({ 9, 8 });
    edges.push_back({ 8, 7 });
    edges.push_back({ 7, 1 });
    edges.push_back({ 9, 10 });
    edges.push_back({ 10, 8 });
    edges.push_back({ 4, 5 });
    edges.push_back({ 5, 10 });
    edges.push_back({ 10, 4 });
    edges.push_back({ 8, 6 });
    edges.push_back({ 6, 7 });
    edges.push_back({ 10, 6 });
    edges.push_back({ 5, 6 });

    for (int i = 0; i < edges.size(); i++) 
    {
        edges[i].first -= 1;
        edges[i].second -= 1;
    }

    // now build node-edge mapping
    Mesh mesh;
    mesh.setState(edges, nodes);
    Orthogonalization<Mesh> orthogonalization;

    orthogonalization.initialize(mesh);

    orthogonalization.iterate(mesh);

    constexpr double tolerance = 1e-2;

    ASSERT_NEAR(325.590101919525, mesh.m_nodes[0].x, tolerance);
    ASSERT_NEAR(229.213730481198, mesh.m_nodes[1].x, tolerance);
    ASSERT_NEAR(263.439319753147, mesh.m_nodes[2].x, tolerance);
    ASSERT_NEAR(429.191105834504, mesh.m_nodes[3].x, tolerance);
    ASSERT_NEAR(535.865215426468, mesh.m_nodes[4].x, tolerance);
    ASSERT_NEAR(503.753784179688, mesh.m_nodes[5].x, tolerance);
    ASSERT_NEAR(354.048340705929, mesh.m_nodes[6].x, tolerance);
    ASSERT_NEAR(346.790050854504, mesh.m_nodes[7].x, tolerance);
    ASSERT_NEAR(315.030130405285, mesh.m_nodes[8].x, tolerance);
    ASSERT_NEAR(424.314957449766, mesh.m_nodes[9].x, tolerance);

    ASSERT_NEAR(455.319334078551, mesh.m_nodes[0].y, tolerance);
    ASSERT_NEAR(362.573521507281, mesh.m_nodes[1].y, tolerance);
    ASSERT_NEAR(241.096458631763, mesh.m_nodes[2].y, tolerance);
    ASSERT_NEAR(211.483073921775, mesh.m_nodes[3].y, tolerance);
    ASSERT_NEAR(311.401495506714, mesh.m_nodes[4].y, tolerance);
    ASSERT_NEAR(432.379974365234, mesh.m_nodes[5].y, tolerance);
    ASSERT_NEAR(458.064836627594, mesh.m_nodes[6].y, tolerance);
    ASSERT_NEAR(405.311585650679, mesh.m_nodes[7].y, tolerance);
    ASSERT_NEAR(319.612138503550, mesh.m_nodes[8].y, tolerance);
    ASSERT_NEAR(327.102805172725, mesh.m_nodes[9].y, tolerance);
}

TEST(OrthogonalizationTest, TestOrthogonalizationFourQuads)
{
    using Mesh = Mesh<cartesianPoint>;

    const int n = 3; //x
    const int m = 3; //y

    std::vector<std::vector<int>> indexesValues(n, std::vector<int>(m));
    std::vector<cartesianPoint> nodes(n * m);
    size_t nodeIndex = 0;
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
    size_t edgeIndex = 0;
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

    // now build node-edge mapping
    Mesh mesh;
    mesh.setState(edges, nodes);
    Orthogonalization<Mesh> orthogonalization;
    orthogonalization.initialize(mesh);
}

