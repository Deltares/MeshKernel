#include "../Mesh.hpp"
#include "../Entities.hpp"
#include "../Polygons.hpp"
#include "../Constants.cpp"
#include "../Orthogonalization.cpp"
#include "MakeMeshes.cpp"
#include <gtest/gtest.h>
#include <chrono>

#if defined(_WIN32)
#include <Windows.h>
#endif

TEST(Orthogonalization, TestOrthogonalizationOneQuadOneTriangle)
{
    // Preparation
    std::vector<GridGeom::Point> nodes;
    nodes.push_back(GridGeom::Point{ 0.0,0.0 });
    nodes.push_back(GridGeom::Point{ 0.0,10.0 });
    nodes.push_back(GridGeom::Point{ 10.0,0.0 });
    nodes.push_back(GridGeom::Point{ 10.0,10.0 });
    nodes.push_back(GridGeom::Point{ 20.0,0.0 });

    std::vector<GridGeom::Edge> edges;
    edges.push_back({ 1, 0 });
    edges.push_back({ 0, 2 });
    edges.push_back({ 2, 4 });
    edges.push_back({ 4, 3 });
    edges.push_back({ 3, 2 });
    edges.push_back({ 3, 1 });

    int isTriangulationRequired = 0;
    int isAccountingForLandBoundariesRequired = 0;
    int projectToLandBoundaryOption = 0;
    GridGeomApi::OrthogonalizationParametersNative orthogonalizationParametersNative;
    orthogonalizationParametersNative.InnerIterations = 2;
    orthogonalizationParametersNative.BoundaryIterations = 25;
    orthogonalizationParametersNative.OuterIterations = 25;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactor = 0.975;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactorBoundary = 0.975;

    // Execute
    GridGeom::Mesh mesh;
    mesh.Set(edges, nodes, GridGeom::Projections::cartesian);

    GridGeom::Orthogonalization orthogonalization;
    std::vector<GridGeom::Point> polygon;
    std::vector<GridGeom::Point> landBoundary;

    orthogonalization.Set(mesh, 
        isTriangulationRequired, 
        isAccountingForLandBoundariesRequired, 
        projectToLandBoundaryOption,
        orthogonalizationParametersNative,
        polygon,
        landBoundary);

    orthogonalization.Iterate(mesh);

    // Assert
    constexpr double tolerance = 1e-8;
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

TEST(Orthogonalization, TestOrthogonalizationSmallTriangularGrid)
{
   
    // now build node-edge mapping
    auto mesh = MakeSmallSizeTriangularMesh();

    int isTriangulationRequired = 0;
    int isAccountingForLandBoundariesRequired = 0;
    int projectToLandBoundaryOption = 0;
    GridGeomApi::OrthogonalizationParametersNative orthogonalizationParametersNative;
    orthogonalizationParametersNative.InnerIterations = 2;
    orthogonalizationParametersNative.BoundaryIterations = 25;
    orthogonalizationParametersNative.OuterIterations = 25;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactor = 0.975;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactorBoundary = 0.975;


    GridGeom::Orthogonalization orthogonalization;
    std::vector<GridGeom::Point> polygon;
    std::vector<GridGeom::Point> landBoundary;

    orthogonalization.Set(mesh,
        isTriangulationRequired,
        isAccountingForLandBoundariesRequired,
        projectToLandBoundaryOption,
        orthogonalizationParametersNative,
        polygon,
        landBoundary);

    orthogonalization.Iterate(mesh);

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

TEST(Orthogonalization, TestOrthogonalizationMediumTriangularGrid)
{
    // now build node-edge mapping
    auto mesh = MakeMediumSizeTriangularMesh();

    int isTriangulationRequired = 0;
    int isAccountingForLandBoundariesRequired = 0;
    int projectToLandBoundaryOption = 0;
    GridGeomApi::OrthogonalizationParametersNative orthogonalizationParametersNative;
    orthogonalizationParametersNative.InnerIterations = 2;
    orthogonalizationParametersNative.BoundaryIterations = 25;
    orthogonalizationParametersNative.OuterIterations = 25;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactor = 0.975;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactorBoundary = 0.5;

    GridGeom::Orthogonalization orthogonalization;
    std::vector<GridGeom::Point> polygon;
    std::vector<GridGeom::Point> landBoundary;

    orthogonalization.Set(mesh,
        isTriangulationRequired,
        isAccountingForLandBoundariesRequired,
        projectToLandBoundaryOption,
        orthogonalizationParametersNative,
        polygon,
        landBoundary);

    orthogonalization.Iterate(mesh);

    constexpr double tolerance = 1.5;

    // check the first 10 points
    ASSERT_NEAR(68.771705432835475, mesh.m_nodes[0].x, tolerance);
    ASSERT_NEAR(169.49338272334273, mesh.m_nodes[1].x, tolerance);
    ASSERT_NEAR(262.80128484924921, mesh.m_nodes[2].x, tolerance);
    ASSERT_NEAR(361.60010033352023, mesh.m_nodes[3].x, tolerance);
    ASSERT_NEAR(468.13991812406925, mesh.m_nodes[4].x, tolerance);
    ASSERT_NEAR(549.89461192844624, mesh.m_nodes[5].x, tolerance);
    ASSERT_NEAR(653.02704974527421, mesh.m_nodes[6].x, tolerance);
    ASSERT_NEAR(747.81537706979441, mesh.m_nodes[7].x, tolerance);
    ASSERT_NEAR(853.40641427112951, mesh.m_nodes[8].x, tolerance);
    ASSERT_NEAR(938.69752431820143, mesh.m_nodes[9].x, tolerance);

    ASSERT_NEAR(1399.7751472360221, mesh.m_nodes[0].y, tolerance);
    ASSERT_NEAR(1426.5945287630802, mesh.m_nodes[1].y, tolerance);
    ASSERT_NEAR(1451.4398281457179, mesh.m_nodes[2].y, tolerance);
    ASSERT_NEAR(1477.7472050498141, mesh.m_nodes[3].y, tolerance);
    ASSERT_NEAR(1506.1157955857589, mesh.m_nodes[4].y, tolerance);
    ASSERT_NEAR(1527.8847968946166, mesh.m_nodes[5].y, tolerance);
    ASSERT_NEAR(1555.3460969050145, mesh.m_nodes[6].y, tolerance);
    ASSERT_NEAR(1580.5855923464549, mesh.m_nodes[7].y, tolerance);
    ASSERT_NEAR(1608.7015489976982, mesh.m_nodes[8].y, tolerance);
    ASSERT_NEAR(1631.412199601948, mesh.m_nodes[9].y, tolerance);

}

TEST(Orthogonalization, TestOrthogonalizationFourQuads)
{

    const int n = 3; //x
    const int m = 3; //y

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

    int isTriangulationRequired = 0;
    int isAccountingForLandBoundariesRequired = 0;
    int projectToLandBoundaryOption = 0;
    GridGeomApi::OrthogonalizationParametersNative orthogonalizationParametersNative;
    orthogonalizationParametersNative.InnerIterations = 2;
    orthogonalizationParametersNative.BoundaryIterations = 25;
    orthogonalizationParametersNative.OuterIterations = 25;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactor = 0.975;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactorBoundary = 0.975;

    // now build node-edge mapping
    GridGeom::Mesh mesh;
    mesh.Set(edges, nodes, GridGeom::Projections::cartesian);

    std::vector<GridGeom::Point> polygon;
    std::vector<GridGeom::Point> landBoundary;

    GridGeom::Orthogonalization orthogonalization;
    orthogonalization.Set(mesh,
        isTriangulationRequired,
        isAccountingForLandBoundariesRequired,
        projectToLandBoundaryOption,
        orthogonalizationParametersNative,
        polygon,
        landBoundary);
}

TEST(Orthogonalization, OrthogonalizeAndSnapToLandBoundaries)
{
    // Prepare
    auto mesh = MakeSmallSizeTriangularMesh();

    // the land boundary to use
    std::vector<GridGeom::Point> landBoundary
    {
        { 235.561218, 290.571899 },
        { 265.953522, 436.515747 },
        { 429.349854, 450.959656 },
        { 535.271545, 386.262909 },
        { GridGeom::doubleMissingValue, GridGeom::doubleMissingValue },
        { 246.995941, 262.285858 },
        { 351.112183, 237.309906 },
        { 443.191895, 262.285858 },
        { 553.627319, 327.283539 },
    };

    // no enclosing polygon
    std::vector<GridGeom::Point> polygon;

    int isTriangulationRequired = 0;
    int isAccountingForLandBoundariesRequired = 1;

    // snap to land boundaries
    int projectToLandBoundaryOption = 2;
    GridGeomApi::OrthogonalizationParametersNative orthogonalizationParametersNative;
    orthogonalizationParametersNative.InnerIterations = 2;
    orthogonalizationParametersNative.BoundaryIterations = 25;
    orthogonalizationParametersNative.OuterIterations = 25;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactor = 0.975;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactorBoundary = 0.975;

    GridGeom::Orthogonalization orthogonalization;
    orthogonalization.Set(mesh,
        isTriangulationRequired,
        isAccountingForLandBoundariesRequired,
        projectToLandBoundaryOption,
        orthogonalizationParametersNative,
        polygon,
        landBoundary);

    orthogonalization.Iterate(mesh);

    // check the values
    constexpr double tolerance = 0.15;
    ASSERT_NEAR(313.081472564480, mesh.m_nodes[0].x, tolerance);
    ASSERT_NEAR(253.641466857330, mesh.m_nodes[1].x, tolerance);
    ASSERT_NEAR(254.777224294204, mesh.m_nodes[2].x, tolerance);
    ASSERT_NEAR(443.191895000000, mesh.m_nodes[3].x, tolerance);
    ASSERT_NEAR(535.240231516760, mesh.m_nodes[4].x, tolerance);
    ASSERT_NEAR(480.436129612752, mesh.m_nodes[5].x, tolerance);
    ASSERT_NEAR(345.948240805397, mesh.m_nodes[6].x, tolerance);
    ASSERT_NEAR(342.668434889472, mesh.m_nodes[7].x, tolerance);
    ASSERT_NEAR(318.414413615199, mesh.m_nodes[8].x, tolerance);
    ASSERT_NEAR(424.616311031376, mesh.m_nodes[9].x, tolerance);

    ASSERT_NEAR(440.681763586650, mesh.m_nodes[0].y, tolerance);
    ASSERT_NEAR(377.393256506700, mesh.m_nodes[1].y, tolerance);
    ASSERT_NEAR(260.419242817573, mesh.m_nodes[2].y, tolerance);
    ASSERT_NEAR(262.285858000000, mesh.m_nodes[3].y, tolerance);
    ASSERT_NEAR(316.461666783032, mesh.m_nodes[4].y, tolerance);
    ASSERT_NEAR(419.756265860671, mesh.m_nodes[5].y, tolerance);
    ASSERT_NEAR(443.587120174434, mesh.m_nodes[6].y, tolerance);
    ASSERT_NEAR(402.913858250569, mesh.m_nodes[7].y, tolerance);
    ASSERT_NEAR(336.831643075189, mesh.m_nodes[8].y, tolerance);
    ASSERT_NEAR(340.875100904741, mesh.m_nodes[9].y, tolerance);
}