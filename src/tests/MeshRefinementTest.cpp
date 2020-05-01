#include "../MeshRefinement.cpp"
#include "../Mesh.hpp"
#include "../Polygons.hpp"
#include "../SampleRefineParametersNative.hpp"
#include "../InterpolationParametersNative.hpp"
#include <gtest/gtest.h>
#include <fstream>

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
    sampleRefineParametersNative.ConnectHangingNodes = 1;

    GridGeomApi::InterpolationParametersNative interpolationParametersNative;
    interpolationParametersNative.MaxNumberOfRefinementIterations = 1;
    
    meshRefinement.Refine(samples, polygon, sampleRefineParametersNative, interpolationParametersNative);
       
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

    // total number of edges
    ASSERT_EQ(84, mesh.GetNumEdges()); 
}


TEST(MeshRefinement, SmallTriangualMeshTwoSamples)
{
    // Prepare
    std::vector<GridGeom::Point> nodes;

    nodes.push_back({ 322.252624511719,454.880187988281 });
    nodes.push_back({ 227.002044677734,360.379241943359 });
    nodes.push_back({ 259.252227783203,241.878051757813 });
    nodes.push_back({ 428.003295898438,210.377746582031 });
    nodes.push_back({ 536.003967285156,310.878753662109 });
    nodes.push_back({ 503.753784179688,432.379974365234 });
    nodes.push_back({ 350.752807617188,458.630249023438 });
    nodes.push_back({ 343.15053976393,406.232256102912 });
    nodes.push_back({ 310.300984548069,319.41005739802 });
    nodes.push_back({ 423.569603308318,326.17986967523 });

    std::vector<GridGeom::Edge> edges;
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

    GridGeom::Mesh mesh;
    mesh.Set(edges, nodes, GridGeom::Projections::cartesian);

    //sample points
    std::vector<GridGeom::Sample> samples;
    samples.push_back({ 359.8657532,350.3144836, 1.0000000 });
    samples.push_back({ 387.5152588 ,299.2614746, 1.0000000 });


    GridGeom::MeshRefinement  meshRefinement(mesh);
    GridGeom::Polygons polygon;
    GridGeomApi::SampleRefineParametersNative sampleRefineParametersNative;
    sampleRefineParametersNative.MaximumTimeStepInCourantGrid = 15.97;
    sampleRefineParametersNative.MinimumCellSize = 50.0;
    sampleRefineParametersNative.AccountForSamplesOutside = false;
    sampleRefineParametersNative.ConnectHangingNodes = 1;

    GridGeomApi::InterpolationParametersNative interpolationParametersNative;
    interpolationParametersNative.MaxNumberOfRefinementIterations = 1;

    meshRefinement.Refine(samples, polygon, sampleRefineParametersNative, interpolationParametersNative);

    // edges connecting hanging nodes
    ASSERT_EQ(10, mesh.m_edges[32].first);
    ASSERT_EQ(2, mesh.m_edges[32].second);

    ASSERT_EQ(14, mesh.m_edges[33].first);
    ASSERT_EQ(4, mesh.m_edges[33].second);

    ASSERT_EQ(13, mesh.m_edges[34].first);
    ASSERT_EQ(5, mesh.m_edges[34].second);

    ASSERT_EQ(11, mesh.m_edges[31].first);
    ASSERT_EQ(1, mesh.m_edges[31].second);


    // total number of edges
    ASSERT_EQ(35, mesh.GetNumEdges());

}

TEST(MeshRefinement, ThreeBythreeWithThreeSamplesPerface)
{
    // Prepare
    //1 Setup
    const int n = 4; // x
    const int m = 4; // y
    double delta = 10.0;

    std::vector<std::vector<int>> indexesValues(n, std::vector<int>(m));
    std::vector<GridGeom::Point> nodes(n * m);
    std::size_t nodeIndex = 0;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            indexesValues[i][j] = i * m + j;
            nodes[nodeIndex] = { i * delta, j * delta };
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

    samples.push_back({ 2.7091951,       5.4000854,       0.0000000 });
    samples.push_back({ 6.4910383,       2.4182367,       0.0000000 });
    samples.push_back({ 8.0910482,      6.7091894,       0.0000000 });
    samples.push_back({ 13.2910795,       5.2182646,       0.0000000 });
    samples.push_back({ 16.1274605,      2.0909605,       0.0000000 });
    samples.push_back({ 18.7820244,       7.5091972,       0.0000000 });
    samples.push_back({ 23.5456886,       8.1637497,       0.0000000 });
    samples.push_back({ 24.6366081,       1.5818644,       0.0000000 });
    samples.push_back({ 27.8729897,       6.5273695,       0.0000000 });
    samples.push_back({ 28.0184441,      14.7092705,       0.0000000 });
    samples.push_back({ 23.8366013,      12.4910660,       0.0000000 });
    samples.push_back({ 22.6002312,      17.2183857,       0.0000000 });
    samples.push_back({ 24.0184212,      23.7639065,       0.0000000 });
    samples.push_back({ 27.9457169,      25.9093838,       0.0000000 });
    samples.push_back({ 24.6366081,      27.3275795,       0.0000000 });
    samples.push_back({ 17.4365616,      27.5821285,       0.0000000 });
    samples.push_back({ 16.6001930,      22.8184433,       0.0000000 });
    samples.push_back({ 12.7092581,      27.2548523,       0.0000000 });
    samples.push_back({ 4.7455721,      27.6912193,       0.0000000 });
    samples.push_back({ 2.6728315,      24.7457352,       0.0000000 });
    samples.push_back({ 7.5819540,      22.9638996,       0.0000000 });
    samples.push_back({ 3.5455656,      15.1820030,       0.0000000 });
    samples.push_back({ 4.8546629,      11.8365135,       0.0000000 });
    samples.push_back({ 8.7455969,      17.2183857,       0.0000000 });
    samples.push_back({ 11.8404741,      17.6817989,       3.0000000 });
    samples.push_back({ 13.5837603,     12.1783361,       3.0000000 });
    samples.push_back({ 17.2156067,      16.9106121,       3.0000000 });

    GridGeom::MeshRefinement  meshRefinement(mesh);
    GridGeom::Polygons polygon;
    GridGeomApi::SampleRefineParametersNative sampleRefineParametersNative;
    sampleRefineParametersNative.MaximumTimeStepInCourantGrid = 0.96;
    sampleRefineParametersNative.MinimumCellSize = 3.0;
    sampleRefineParametersNative.AccountForSamplesOutside = false;
    sampleRefineParametersNative.ConnectHangingNodes = 1;

    GridGeomApi::InterpolationParametersNative interpolationParametersNative;
    interpolationParametersNative.MaxNumberOfRefinementIterations = 2;

    meshRefinement.Refine(samples, polygon, sampleRefineParametersNative, interpolationParametersNative);

    // total number of edges
    ASSERT_EQ(150, mesh.GetNumNodes());
    ASSERT_EQ(293, mesh.GetNumEdges());

    // assert on newly generated edges
    ASSERT_EQ(109, mesh.m_edges[282].first);
    ASSERT_EQ(44, mesh.m_edges[282].second);

    ASSERT_EQ(44, mesh.m_edges[283].first);
    ASSERT_EQ(67, mesh.m_edges[283].second);

    ASSERT_EQ(92, mesh.m_edges[284].first);
    ASSERT_EQ(91, mesh.m_edges[284].second);

    ASSERT_EQ(91, mesh.m_edges[285].first);
    ASSERT_EQ(12, mesh.m_edges[285].second);

    ASSERT_EQ(12, mesh.m_edges[286].first);
    ASSERT_EQ(92, mesh.m_edges[286].second);

    ASSERT_EQ(100, mesh.m_edges[287].first);
    ASSERT_EQ(99, mesh.m_edges[287].second);

    ASSERT_EQ(99, mesh.m_edges[288].first);
    ASSERT_EQ(14, mesh.m_edges[288].second);

    ASSERT_EQ(14, mesh.m_edges[289].first);
    ASSERT_EQ(100, mesh.m_edges[289].second);

    ASSERT_EQ(97, mesh.m_edges[290].first);
    ASSERT_EQ(96, mesh.m_edges[290].second);

    ASSERT_EQ(96, mesh.m_edges[291].first);
    ASSERT_EQ(14, mesh.m_edges[291].second);

    ASSERT_EQ(14, mesh.m_edges[292].first);
    ASSERT_EQ(97, mesh.m_edges[292].second);

}

TEST(MeshRefinement, WindowOfRefinementFile)
{
    // Prepare
    //1 Setup
    const int n = 4; // x
    const int m = 4; // y
    double delta = 40.0;
    double originX = 197253;
    double originY = 442281;

    std::vector<std::vector<int>> indexesValues(n, std::vector<int>(m));
    std::vector<GridGeom::Point> nodes(n * m);
    std::size_t nodeIndex = 0;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            indexesValues[i][j] = i * m + j;
            nodes[nodeIndex] = { originX + i * delta, originY + j * delta };
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

    // read sample file
    std::string line;
    std::ifstream infile("..\\..\\tests\\WindowedSamplesForMeshRefinement.xyz");
    while (std::getline(infile, line))  
    {
        std::istringstream iss(line);
        double sampleX;
        double sampleY;
        double sampleValue;

        iss >> sampleX;
        iss >> sampleY;
        iss >> sampleValue;

        samples.push_back({ sampleX , sampleY , sampleValue });
    }


    GridGeom::MeshRefinement  meshRefinement(mesh);
    GridGeom::Polygons polygon;
    GridGeomApi::SampleRefineParametersNative sampleRefineParametersNative;
    sampleRefineParametersNative.MaximumTimeStepInCourantGrid = 0.96;
    sampleRefineParametersNative.MinimumCellSize = 3.0;
    sampleRefineParametersNative.AccountForSamplesOutside = false;
    sampleRefineParametersNative.ConnectHangingNodes = 1;

    GridGeomApi::InterpolationParametersNative interpolationParametersNative;
    interpolationParametersNative.MaxNumberOfRefinementIterations = 4;

    meshRefinement.Refine(samples, polygon, sampleRefineParametersNative, interpolationParametersNative);

    // total number of edges
    ASSERT_EQ(1614, mesh.GetNumNodes());
    ASSERT_EQ(3216, mesh.GetNumEdges());

    ASSERT_EQ(170, mesh.m_edges[3113].first);
    ASSERT_EQ(804, mesh.m_edges[3113].second);

    ASSERT_EQ(1, mesh.m_edges[3114].first);
    ASSERT_EQ(804, mesh.m_edges[3114].second);

    ASSERT_EQ(462, mesh.m_edges[3115].first);
    ASSERT_EQ(856, mesh.m_edges[3115].second);

    ASSERT_EQ(9, mesh.m_edges[3116].first);
    ASSERT_EQ(856, mesh.m_edges[3116].second);

    ASSERT_EQ(211, mesh.m_edges[3117].first);
    ASSERT_EQ(1256, mesh.m_edges[3117].second);

    ASSERT_EQ(538, mesh.m_edges[3118].first);
    ASSERT_EQ(1256, mesh.m_edges[3118].second);

    ASSERT_EQ(77, mesh.m_edges[3119].first);
    ASSERT_EQ(1052, mesh.m_edges[3119].second);

    ASSERT_EQ(258, mesh.m_edges[3120].first);
    ASSERT_EQ(1052, mesh.m_edges[3120].second);

    ASSERT_EQ(290, mesh.m_edges[3121].first);
    ASSERT_EQ(769, mesh.m_edges[3121].second);

    ASSERT_EQ(593, mesh.m_edges[3122].first);
    ASSERT_EQ(769, mesh.m_edges[3122].second);
}

TEST(MeshRefinement, RefineBasedOnPolygon)
{
    // Prepare
    //1 Setup
    const int n = 5; 
    const int m = 5; 
    double delta = 10.0;
    double originX = 0.0;
    double originY = 0.0;

    std::vector<std::vector<int>> indexesValues(n, std::vector<int>(m));
    std::vector<GridGeom::Point> nodes(n * m);
    std::size_t nodeIndex = 0;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            indexesValues[i][j] = i * m + j;
            nodes[nodeIndex] = { originX + i * delta, originY + j * delta };
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

    GridGeom::MeshRefinement  meshRefinement(mesh);

    std::vector<GridGeom::Point> point;
    point.push_back({ 25.0 ,-10.0 });
    point.push_back({ 25.0  ,15.0 });
    point.push_back({ 45.0  ,15.0 });
    point.push_back({ 45.0 ,-10.0 });
    point.push_back({ 25.0 ,-10.0 });

    GridGeom::Polygons polygon;
    polygon.Set(point,mesh.m_projection);

    GridGeomApi::SampleRefineParametersNative sampleRefineParametersNative;
    sampleRefineParametersNative.MaximumTimeStepInCourantGrid = 0.96;
    sampleRefineParametersNative.MinimumCellSize = 3.0;
    sampleRefineParametersNative.AccountForSamplesOutside = false;
    sampleRefineParametersNative.ConnectHangingNodes = 1;

    GridGeomApi::InterpolationParametersNative interpolationParametersNative;
    interpolationParametersNative.MaxNumberOfRefinementIterations = 1;

    meshRefinement.Refine(samples, polygon, sampleRefineParametersNative, interpolationParametersNative);

    // total number of edges
    ASSERT_EQ(30, mesh.GetNumNodes());
    ASSERT_EQ(52, mesh.GetNumEdges());

    ASSERT_EQ(25, mesh.m_edges[40].first);
    ASSERT_EQ(29, mesh.m_edges[40].second);

    ASSERT_EQ(28, mesh.m_edges[41].first);
    ASSERT_EQ(29, mesh.m_edges[41].second);

    ASSERT_EQ(26, mesh.m_edges[42].first);
    ASSERT_EQ(29, mesh.m_edges[42].second);

    ASSERT_EQ(27, mesh.m_edges[43].first);
    ASSERT_EQ(29, mesh.m_edges[43].second);

    ASSERT_EQ(25, mesh.m_edges[44].first);
    ASSERT_EQ(20, mesh.m_edges[44].second);

    ASSERT_EQ(26, mesh.m_edges[45].first);
    ASSERT_EQ(21, mesh.m_edges[45].second);

    ASSERT_EQ(27, mesh.m_edges[46].first);
    ASSERT_EQ(15, mesh.m_edges[46].second);

    ASSERT_EQ(28, mesh.m_edges[47].first);
    ASSERT_EQ(20, mesh.m_edges[47].second);

    ASSERT_EQ(10, mesh.m_edges[48].first);
    ASSERT_EQ(27, mesh.m_edges[48].second);

}