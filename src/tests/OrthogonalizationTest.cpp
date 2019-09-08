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
    edges.push_back({ 0, 1 });
    edges.push_back({ 0, 2 });
    edges.push_back({ 2, 4 });
    edges.push_back({ 3, 4 });
    edges.push_back({ 2, 3 });
    edges.push_back({ 1, 3 });

    // now build node-edge mapping
    Mesh mesh(edges, nodes);
    Orthogonalization<cartesianPoint> orthogonalization;

    std::vector<double> aspectRatio;
    orthogonalization.initialize(mesh);
    orthogonalization.aspectRatio(mesh);
    orthogonalization.solveWeights(mesh);

    constexpr  double tolerance = 1e-8;
    //node 0
    EXPECT_EQ(2, orthogonalization.m_topologyFaces[0]);
    EXPECT_EQ(4, orthogonalization.m_topologyNodes[0]);

    ASSERT_NEAR(0.0, orthogonalization.m_topologyXi[0][0], tolerance);
    ASSERT_NEAR(1.0, orthogonalization.m_topologyXi[0][1], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_topologyXi[0][2], tolerance);
    ASSERT_NEAR(1.0, orthogonalization.m_topologyXi[0][3], tolerance);

    ASSERT_NEAR(0.0, orthogonalization.m_topologyEta[0][0], tolerance);
    ASSERT_NEAR(0.0, orthogonalization.m_topologyEta[0][1], tolerance);
    ASSERT_NEAR(1.0, orthogonalization.m_topologyEta[0][2], tolerance);
    ASSERT_NEAR(1.0, orthogonalization.m_topologyEta[0][3], tolerance);

    //node 1
    EXPECT_EQ(2, orthogonalization.m_topologyFaces[1]);
    EXPECT_EQ(4, orthogonalization.m_topologyNodes[1]);

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
    EXPECT_EQ(3, orthogonalization.m_topologyFaces[2]);
    EXPECT_EQ(4, orthogonalization.m_topologyNodes[2]);

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
    EXPECT_EQ(3, orthogonalization.m_topologyFaces[3]);
    EXPECT_EQ(4, orthogonalization.m_topologyNodes[3]);

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
    EXPECT_EQ(2, orthogonalization.m_topologyFaces[4]);
    EXPECT_EQ(3, orthogonalization.m_topologyNodes[4]);

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
    Mesh mesh(edges, nodes);
    Orthogonalization<cartesianPoint> orthogonalization;

    std::vector<double> aspectRatio;
    orthogonalization.initialize(mesh);
    orthogonalization.aspectRatio(mesh);
    
    EXPECT_EQ(1, orthogonalization.m_aspectRatios[0]);
    EXPECT_EQ(1, orthogonalization.m_aspectRatios[1]);
    EXPECT_EQ(1, orthogonalization.m_aspectRatios[2]);
    EXPECT_EQ(1, orthogonalization.m_aspectRatios[3]);

    orthogonalization.solveWeights(mesh);

}

TEST(OrthogonalizationTest, TestOrthogonalizationFourQuads)
{
    using Mesh = Mesh<cartesianPoint>;

    const int n = 3; //x
    const int m = 3; //y

    std::vector<std::vector<int>> indexesValues(n, std::vector<int>(m));
    std::vector<GridGeom::cartesianPoint> nodes(n * m);
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
    Mesh mesh(edges, nodes);
    Orthogonalization<cartesianPoint> orthogonalization;

    std::vector<double> aspectRatio;
    orthogonalization.initialize(mesh);
    orthogonalization.aspectRatio(mesh);
    orthogonalization.solveWeights(mesh);
}

