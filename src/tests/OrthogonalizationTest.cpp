#include "GridGeomTest.hpp"
#include "../Orthogonaization.cpp"

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
    
    //EXPECT_EQ(1, orthogonalization.m_aspectRatios[0]);
    //EXPECT_EQ(1, orthogonalization.m_aspectRatios[1]);
    //EXPECT_EQ(1, orthogonalization.m_aspectRatios[2]);
    //EXPECT_EQ(1, orthogonalization.m_aspectRatios[3]);

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

    //EXPECT_EQ(1, orthogonalization.m_aspectRatios[0]);
    //EXPECT_EQ(1, orthogonalization.m_aspectRatios[1]);
    //EXPECT_EQ(1, orthogonalization.m_aspectRatios[2]);
    //EXPECT_EQ(1, orthogonalization.m_aspectRatios[3]);

    orthogonalization.solveWeights(mesh);
}

