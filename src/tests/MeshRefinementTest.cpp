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

    GridGeomApi::InterpolationParametersNative interpolationParametersNative;
    
    meshRefinement.RefineMeshBasedOnSamples(samples, polygon, sampleRefineParametersNative, interpolationParametersNative);
       
    // 3 Validation
    // expect nodesEdges to be sorted ccw
    ASSERT_EQ(0, mesh.m_nodesEdges[0][0]);
  
}
