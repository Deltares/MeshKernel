#include <MeshKernel/Mesh.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/Constants.hpp>
#include <MeshKernel/LandBoundaries.hpp>
#include <TestUtils/MakeMeshes.hpp>
#include <gtest/gtest.h>

TEST(LandBoundaries, OneLandBoundary)
{
    // Prepare
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/SmallTriangularGrid_net.nc");

    std::vector<meshkernel::Point> landBoundaryPolygon{
        {222.621918, 382.651917},
        {316.206177, 461.190796},
        {350.811279, 465.102692},
        {510.295715, 438.923065},
        {meshkernel::doubleMissingValue, meshkernel::doubleMissingValue}};

    auto polygons = std::make_shared<meshkernel::Polygons>();

    // Execute
    auto landboundaries = std::make_shared<meshkernel::LandBoundaries>(landBoundaryPolygon, mesh, polygons);
    landboundaries->FindNearestMeshBoundary(2);

    // Checks
    EXPECT_EQ(1, landboundaries->m_meshNodesLandBoundarySegments[0]);
    EXPECT_EQ(0, landboundaries->m_meshNodesLandBoundarySegments[1]);
    EXPECT_EQ(2, landboundaries->m_meshNodesLandBoundarySegments[2]);
    EXPECT_EQ(2, landboundaries->m_meshNodesLandBoundarySegments[3]);
    EXPECT_EQ(3, landboundaries->m_meshNodesLandBoundarySegments[4]);
    EXPECT_EQ(1, landboundaries->m_meshNodesLandBoundarySegments[5]);
    EXPECT_EQ(1, landboundaries->m_meshNodesLandBoundarySegments[6]);
    EXPECT_EQ(-1, landboundaries->m_meshNodesLandBoundarySegments[7]);
    EXPECT_EQ(-1, landboundaries->m_meshNodesLandBoundarySegments[8]);
    EXPECT_EQ(-1, landboundaries->m_meshNodesLandBoundarySegments[9]);
}

TEST(LandBoundaries, TwoLandBoundaries)
{
    // Prepare
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/SmallTriangularGrid_net.nc");

    std::vector<meshkernel::Point> landBoundaryPolygon{
        {222.621918, 382.651917},
        {316.206177, 461.190796},
        {350.811279, 465.102692},
        {510.295715, 438.923065},
        {meshkernel::doubleMissingValue, meshkernel::doubleMissingValue},
        {215.980743, 363.986420},
        {250.253036, 235.233246},
        {423.158325, 200.652054},
        {559.938782, 312.732147},
        {518.873718, 421.415894},
        {meshkernel::doubleMissingValue, meshkernel::doubleMissingValue}};

    auto polygons = std::make_shared<meshkernel::Polygons>();

    // Execute
    auto landboundaries = std::make_shared<meshkernel::LandBoundaries>(landBoundaryPolygon, mesh, polygons);
    landboundaries->FindNearestMeshBoundary(2);

    // Checks
    EXPECT_EQ(2, landboundaries->m_meshNodesLandBoundarySegments[0]);
    EXPECT_EQ(1, landboundaries->m_meshNodesLandBoundarySegments[1]);
    EXPECT_EQ(1, landboundaries->m_meshNodesLandBoundarySegments[2]);
    EXPECT_EQ(3, landboundaries->m_meshNodesLandBoundarySegments[3]);
    EXPECT_EQ(3, landboundaries->m_meshNodesLandBoundarySegments[4]);
    EXPECT_EQ(2, landboundaries->m_meshNodesLandBoundarySegments[5]);
    EXPECT_EQ(2, landboundaries->m_meshNodesLandBoundarySegments[6]);
    EXPECT_EQ(-1, landboundaries->m_meshNodesLandBoundarySegments[7]);
    EXPECT_EQ(-1, landboundaries->m_meshNodesLandBoundarySegments[8]);
    EXPECT_EQ(-1, landboundaries->m_meshNodesLandBoundarySegments[9]);
}

TEST(LandBoundaries, OneCrossingLandBoundary)
{
    // Prepare
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/SmallTriangularGrid_net.nc");

    std::vector<meshkernel::Point> landBoundaryPolygon{
        {221.418243, 315.848755},
        {248.801422, 375.129028},
        {337.571045, 459.686218},
        {516.313965, 419.062683},
        {528.651428, 292.377380},
        {meshkernel::doubleMissingValue, meshkernel::doubleMissingValue}};

    auto polygons = std::make_shared<meshkernel::Polygons>();

    // Execute
    auto landboundaries = std::make_shared<meshkernel::LandBoundaries>(landBoundaryPolygon, mesh, polygons);
    landboundaries->FindNearestMeshBoundary(2);

    // Checks
    EXPECT_EQ(0, landboundaries->m_meshNodesLandBoundarySegments[0]);
    EXPECT_EQ(0, landboundaries->m_meshNodesLandBoundarySegments[1]);
    EXPECT_EQ(2, landboundaries->m_meshNodesLandBoundarySegments[2]);
    EXPECT_EQ(2, landboundaries->m_meshNodesLandBoundarySegments[3]);
    EXPECT_EQ(1, landboundaries->m_meshNodesLandBoundarySegments[4]);
    EXPECT_EQ(1, landboundaries->m_meshNodesLandBoundarySegments[5]);
    EXPECT_EQ(1, landboundaries->m_meshNodesLandBoundarySegments[6]);
    EXPECT_EQ(-1, landboundaries->m_meshNodesLandBoundarySegments[7]);
    EXPECT_EQ(-1, landboundaries->m_meshNodesLandBoundarySegments[8]);
    EXPECT_EQ(-1, landboundaries->m_meshNodesLandBoundarySegments[9]);
}

TEST(LandBoundaries, TwoCrossingLandBoundary)
{
    // Prepare
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/SmallTriangularGrid_net.nc");

    std::vector<meshkernel::Point> landBoundaryPolygon{
        {235.561218, 290.571899},
        {265.953522, 436.515747},
        {429.349854, 450.959656},
        {535.271545, 386.262909},
        {meshkernel::doubleMissingValue, meshkernel::doubleMissingValue},
        {246.995941, 262.285858},
        {351.112183, 237.309906},
        {443.191895, 262.285858},
        {553.627319, 327.283539},
        {meshkernel::doubleMissingValue, meshkernel::doubleMissingValue}};

    auto polygons = std::make_shared<meshkernel::Polygons>();

    // Execute
    auto landboundaries = std::make_shared<meshkernel::LandBoundaries>(landBoundaryPolygon, mesh, polygons);
    landboundaries->FindNearestMeshBoundary(2);

    // Checks
    EXPECT_EQ(2, landboundaries->m_meshNodesLandBoundarySegments[0]);
    EXPECT_EQ(0, landboundaries->m_meshNodesLandBoundarySegments[1]);
    EXPECT_EQ(1, landboundaries->m_meshNodesLandBoundarySegments[2]);
    EXPECT_EQ(3, landboundaries->m_meshNodesLandBoundarySegments[3]);
    EXPECT_EQ(3, landboundaries->m_meshNodesLandBoundarySegments[4]);
    EXPECT_EQ(2, landboundaries->m_meshNodesLandBoundarySegments[5]);
    EXPECT_EQ(2, landboundaries->m_meshNodesLandBoundarySegments[6]);
    EXPECT_EQ(-1, landboundaries->m_meshNodesLandBoundarySegments[7]);
    EXPECT_EQ(-1, landboundaries->m_meshNodesLandBoundarySegments[8]);
    EXPECT_EQ(-1, landboundaries->m_meshNodesLandBoundarySegments[9]);
}
