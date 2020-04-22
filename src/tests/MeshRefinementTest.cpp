#include "../MeshRefinement.cpp"
#include "../Mesh.hpp"
#include "../Polygons.hpp"
#include "../SampleRefineParametersNative.hpp"
#include "../InterpolationParametersNative.hpp"
#include <gtest/gtest.h>

TEST(MeshRefinement, FourByFourWithFourSamples) 
{
    //1 Setup
    const int n = 5; // x
    const int m = 5; // y
    double delta = 10.0;

    std::vector<std::vector<int>> indexesValues(n, std::vector<int>(m));
    std::vector<GridGeom::Point> nodes(n * m);
    std::size_t nodeIndex = 0;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            indexesValues[i][j] = i * m + j;
            nodes[nodeIndex] = { i *delta, j*delta };
            nodeIndex++;
        }
    }

    std::vector<GridGeom::Edge> edges((n - 1) * m + (m - 1) * n);
    std::size_t edgeIndex = 0;

    for (int i = 0; i < n - 1; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            edges[edgeIndex] = { indexesValues[i][j], indexesValues[i + 1][j] };
            edgeIndex++;
        }
    }

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m - 1; ++j)
        {
            edges[edgeIndex] = { indexesValues[i][j + 1], indexesValues[i][j] };
            edgeIndex++;
        }
    }

    GridGeom::Mesh mesh;
    mesh.Set(edges, nodes, GridGeom::Projections::cartesian);

    //sample points
    std::vector<GridGeom::Sample> samples;
    samples.push_back({ 14.7153645, 14.5698833, 1.0 });
    samples.push_back({ 24.7033062, 14.4729137, 1.0 });
    samples.push_back({ 15.5396099, 24.2669525, 1.0 });
    samples.push_back({ 23.8305721, 23.9275551, 1.0 });


    GridGeom::MeshRefinement  meshRefinement(mesh);
    GridGeom::Polygons polygon;
    GridGeomApi::SampleRefineParametersNative sampleRefineParametersNative;
    sampleRefineParametersNative.MaximumTimeStepInCourantGrid = 0.32;
    sampleRefineParametersNative.MinimumCellSize = 1.0;
    sampleRefineParametersNative.AccountForSamplesOutside = false;
    sampleRefineParametersNative.MaxNumberOfRefinementIterations = 1;
    sampleRefineParametersNative.ConnectHangingNodes = 1;

    GridGeomApi::InterpolationParametersNative interpolationParametersNative;
    
    meshRefinement.RefineMeshBasedOnSamples(samples, polygon, sampleRefineParametersNative, interpolationParametersNative);
       
    // 3 Validation edges connecting hanging nodes

    //bottom side
    ASSERT_EQ(5, mesh.m_edges[73].first);
    ASSERT_EQ(25, mesh.m_edges[73].second);

    ASSERT_EQ(10, mesh.m_edges[72].first);
    ASSERT_EQ(25, mesh.m_edges[72].second);

    ASSERT_EQ(10, mesh.m_edges[77].first);
    ASSERT_EQ(28, mesh.m_edges[77].second);

    ASSERT_EQ(15, mesh.m_edges[76].first);
    ASSERT_EQ(28, mesh.m_edges[76].second);

    //right side
    ASSERT_EQ(21, mesh.m_edges[81].first);
    ASSERT_EQ(35, mesh.m_edges[81].second);

    ASSERT_EQ(22, mesh.m_edges[80].first);
    ASSERT_EQ(35, mesh.m_edges[80].second);

    ASSERT_EQ(22, mesh.m_edges[83].first);
    ASSERT_EQ(36, mesh.m_edges[83].second);

    ASSERT_EQ(23, mesh.m_edges[82].first);
    ASSERT_EQ(36, mesh.m_edges[82].second);

    //upper side
    ASSERT_EQ(19, mesh.m_edges[79].first);
    ASSERT_EQ(30, mesh.m_edges[79].second);

    ASSERT_EQ(14, mesh.m_edges[78].first);
    ASSERT_EQ(30, mesh.m_edges[78].second);

    ASSERT_EQ(14, mesh.m_edges[75].first);
    ASSERT_EQ(27, mesh.m_edges[75].second);

    ASSERT_EQ(9, mesh.m_edges[74].first);
    ASSERT_EQ(27, mesh.m_edges[74].second);


    //left side
    ASSERT_EQ(3, mesh.m_edges[71].first);
    ASSERT_EQ(32, mesh.m_edges[71].second);

    ASSERT_EQ(2, mesh.m_edges[70].first);
    ASSERT_EQ(32, mesh.m_edges[70].second);

    ASSERT_EQ(2, mesh.m_edges[69].first);
    ASSERT_EQ(31, mesh.m_edges[69].second);

    ASSERT_EQ(1, mesh.m_edges[68].first);
    ASSERT_EQ(31, mesh.m_edges[68].second);

}
