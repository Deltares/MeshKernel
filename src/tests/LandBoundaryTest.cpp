#include "GridGeomTest.hpp"
#include "../Mesh.cpp"
#include "../LandBoundaries.hpp"
#include "../Constants.cpp"
#include <gtest/gtest.h>

TEST(LandBoundariesTests, OneLandBoundary)
{
    // Prepare
    std::vector<GridGeom::Point> nodes;

    nodes.push_back(GridGeom::Point{ 322.252624511719,454.880187988281 });
    nodes.push_back(GridGeom::Point{ 227.002044677734,360.379241943359 });
    nodes.push_back(GridGeom::Point{ 259.252227783203,241.878051757813 });
    nodes.push_back(GridGeom::Point{ 428.003295898438,210.377746582031 });
    nodes.push_back(GridGeom::Point{ 536.003967285156,310.878753662109 });
    nodes.push_back(GridGeom::Point{ 503.753784179688,432.379974365234 });
    nodes.push_back(GridGeom::Point{ 350.752807617188,458.630249023438 });
    nodes.push_back(GridGeom::Point{ 343.15053976393,406.232256102912 });
    nodes.push_back(GridGeom::Point{ 310.300984548069,319.41005739802 });
    nodes.push_back(GridGeom::Point{ 423.569603308318,326.17986967523 });

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

    GridGeom::LandBoundaries landboundaries(mesh.m_projection);
    std::vector<GridGeom::Point> landBoundaryPolygon
    {
        { 222.621918, 382.651917 },
        { 316.206177, 461.190796 },
        { 350.811279, 465.102692 },
        { 510.295715, 438.923065 }
    };

    GridGeom::Polygons polygons;

    // Execute
    landboundaries.Set(landBoundaryPolygon);
    landboundaries.Administrate(mesh, polygons);
    landboundaries.FindNearestMeshBoundary(mesh, polygons, 2);

    // Checks
    EXPECT_EQ(1, landboundaries.m_meshNodesLandBoundarySegments[0]);
    EXPECT_EQ(0, landboundaries.m_meshNodesLandBoundarySegments[1]);
    EXPECT_EQ(2, landboundaries.m_meshNodesLandBoundarySegments[2]);
    EXPECT_EQ(2, landboundaries.m_meshNodesLandBoundarySegments[3]);
    EXPECT_EQ(3, landboundaries.m_meshNodesLandBoundarySegments[4]);
    EXPECT_EQ(1, landboundaries.m_meshNodesLandBoundarySegments[5]);
    EXPECT_EQ(1, landboundaries.m_meshNodesLandBoundarySegments[6]);
    EXPECT_EQ(-1, landboundaries.m_meshNodesLandBoundarySegments[7]);
    EXPECT_EQ(-1, landboundaries.m_meshNodesLandBoundarySegments[8]);
    EXPECT_EQ(-1, landboundaries.m_meshNodesLandBoundarySegments[9]);
}

TEST(LandBoundariesTests, TwoLandBoundaries)
{
    // Prepare
    std::vector<GridGeom::Point> nodes;

    nodes.push_back(GridGeom::Point{ 322.252624511719,454.880187988281 });
    nodes.push_back(GridGeom::Point{ 227.002044677734,360.379241943359 });
    nodes.push_back(GridGeom::Point{ 259.252227783203,241.878051757813 });
    nodes.push_back(GridGeom::Point{ 428.003295898438,210.377746582031 });
    nodes.push_back(GridGeom::Point{ 536.003967285156,310.878753662109 });
    nodes.push_back(GridGeom::Point{ 503.753784179688,432.379974365234 });
    nodes.push_back(GridGeom::Point{ 350.752807617188,458.630249023438 });
    nodes.push_back(GridGeom::Point{ 343.15053976393,406.232256102912 });
    nodes.push_back(GridGeom::Point{ 310.300984548069,319.41005739802 });
    nodes.push_back(GridGeom::Point{ 423.569603308318,326.17986967523 });

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


    GridGeom::LandBoundaries landboundaries(mesh.m_projection);
    std::vector<GridGeom::Point> landBoundaryPolygon
    {
        { 222.621918, 382.651917 },
        { 316.206177, 461.190796 },
        { 350.811279, 465.102692 },
        { 510.295715, 438.923065 },
        { GridGeom::doubleMissingValue, GridGeom::doubleMissingValue },
        { 215.980743, 363.986420 },
        { 250.253036, 235.233246 },
        { 423.158325, 200.652054 },
        { 559.938782, 312.732147 },
        { 518.873718, 421.415894 }
    };

    GridGeom::Polygons polygons;

    // Execute
    landboundaries.Set(landBoundaryPolygon);
    landboundaries.Administrate(mesh, polygons);
    landboundaries.FindNearestMeshBoundary(mesh, polygons, 2);

    // Checks
    EXPECT_EQ(2, landboundaries.m_meshNodesLandBoundarySegments[0]);
    EXPECT_EQ(1, landboundaries.m_meshNodesLandBoundarySegments[1]);
    EXPECT_EQ(1, landboundaries.m_meshNodesLandBoundarySegments[2]);
    EXPECT_EQ(3, landboundaries.m_meshNodesLandBoundarySegments[3]);
    EXPECT_EQ(3, landboundaries.m_meshNodesLandBoundarySegments[4]);
    EXPECT_EQ(2, landboundaries.m_meshNodesLandBoundarySegments[5]);
    EXPECT_EQ(2, landboundaries.m_meshNodesLandBoundarySegments[6]);
    EXPECT_EQ(-1, landboundaries.m_meshNodesLandBoundarySegments[7]);
    EXPECT_EQ(-1, landboundaries.m_meshNodesLandBoundarySegments[8]);
    EXPECT_EQ(-1, landboundaries.m_meshNodesLandBoundarySegments[9]);
}

TEST(LandBoundariesTests, OneCrossingLandBoundary)
{
    // Prepare
    std::vector<GridGeom::Point> nodes;

    nodes.push_back(GridGeom::Point{ 322.252624511719,454.880187988281 });
    nodes.push_back(GridGeom::Point{ 227.002044677734,360.379241943359 });
    nodes.push_back(GridGeom::Point{ 259.252227783203,241.878051757813 });
    nodes.push_back(GridGeom::Point{ 428.003295898438,210.377746582031 });
    nodes.push_back(GridGeom::Point{ 536.003967285156,310.878753662109 });
    nodes.push_back(GridGeom::Point{ 503.753784179688,432.379974365234 });
    nodes.push_back(GridGeom::Point{ 350.752807617188,458.630249023438 });
    nodes.push_back(GridGeom::Point{ 343.15053976393,406.232256102912 });
    nodes.push_back(GridGeom::Point{ 310.300984548069,319.41005739802 });
    nodes.push_back(GridGeom::Point{ 423.569603308318,326.17986967523 });

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


    GridGeom::LandBoundaries landboundaries(mesh.m_projection);
    std::vector<GridGeom::Point> landBoundaryPolygon
    {
        { 221.418243, 315.848755 },
        { 248.801422, 375.129028 },
        { 337.571045, 459.686218 },
        { 516.313965, 419.062683 },
        { 528.651428, 292.377380 }
    };

    GridGeom::Polygons polygons;

    // Execute
    landboundaries.Set(landBoundaryPolygon);
    landboundaries.Administrate(mesh, polygons);
    landboundaries.FindNearestMeshBoundary(mesh, polygons, 2);

    // Checks
    EXPECT_EQ(0, landboundaries.m_meshNodesLandBoundarySegments[0]);
    EXPECT_EQ(0, landboundaries.m_meshNodesLandBoundarySegments[1]);
    EXPECT_EQ(2, landboundaries.m_meshNodesLandBoundarySegments[2]);
    EXPECT_EQ(2, landboundaries.m_meshNodesLandBoundarySegments[3]);
    EXPECT_EQ(1, landboundaries.m_meshNodesLandBoundarySegments[4]);
    EXPECT_EQ(1, landboundaries.m_meshNodesLandBoundarySegments[5]);
    EXPECT_EQ(1, landboundaries.m_meshNodesLandBoundarySegments[6]);
    EXPECT_EQ(-1, landboundaries.m_meshNodesLandBoundarySegments[7]);
    EXPECT_EQ(-1, landboundaries.m_meshNodesLandBoundarySegments[8]);
    EXPECT_EQ(-1, landboundaries.m_meshNodesLandBoundarySegments[9]);
}

TEST(LandBoundariesTests, TwoCrossingLandBoundary)
{
    // Prepare
    std::vector<GridGeom::Point> nodes;

    nodes.push_back(GridGeom::Point{ 322.252624511719,454.880187988281 });
    nodes.push_back(GridGeom::Point{ 227.002044677734,360.379241943359 });
    nodes.push_back(GridGeom::Point{ 259.252227783203,241.878051757813 });
    nodes.push_back(GridGeom::Point{ 428.003295898438,210.377746582031 });
    nodes.push_back(GridGeom::Point{ 536.003967285156,310.878753662109 });
    nodes.push_back(GridGeom::Point{ 503.753784179688,432.379974365234 });
    nodes.push_back(GridGeom::Point{ 350.752807617188,458.630249023438 });
    nodes.push_back(GridGeom::Point{ 343.15053976393,406.232256102912 });
    nodes.push_back(GridGeom::Point{ 310.300984548069,319.41005739802 });
    nodes.push_back(GridGeom::Point{ 423.569603308318,326.17986967523 });

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


    GridGeom::LandBoundaries landboundaries(mesh.m_projection);
    std::vector<GridGeom::Point> landBoundaryPolygon
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

    GridGeom::Polygons polygons;

    // Execute
    landboundaries.Set(landBoundaryPolygon);
    landboundaries.Administrate(mesh, polygons);
    landboundaries.FindNearestMeshBoundary(mesh, polygons, 2);

    // Checks
    EXPECT_EQ(2, landboundaries.m_meshNodesLandBoundarySegments[0]);
    EXPECT_EQ(0, landboundaries.m_meshNodesLandBoundarySegments[1]);
    EXPECT_EQ(1, landboundaries.m_meshNodesLandBoundarySegments[2]);
    EXPECT_EQ(3, landboundaries.m_meshNodesLandBoundarySegments[3]);
    EXPECT_EQ(3, landboundaries.m_meshNodesLandBoundarySegments[4]);
    EXPECT_EQ(2, landboundaries.m_meshNodesLandBoundarySegments[5]);
    EXPECT_EQ(2, landboundaries.m_meshNodesLandBoundarySegments[6]);
    EXPECT_EQ(-1, landboundaries.m_meshNodesLandBoundarySegments[7]);
    EXPECT_EQ(-1, landboundaries.m_meshNodesLandBoundarySegments[8]);
    EXPECT_EQ(-1, landboundaries.m_meshNodesLandBoundarySegments[9]);
}