#include <gtest/gtest.h>

#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/LandBoundaries.hpp>
#include <MeshKernel/LandBoundary.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Polygons.hpp>
#include <TestUtils/Definitions.hpp>
#include <TestUtils/MakeMeshes.hpp>

TEST(LandBoundaries, OneLandBoundary)
{
    // Prepare
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/SmallTriangularGrid_net.nc");

    std::vector<meshkernel::Point> landBoundaryPolygon{
        {222.621918, 382.651917},
        {316.206177, 461.190796},
        {350.811279, 465.102692},
        {510.295715, 438.923065},
        {meshkernel::constants::missing::doubleValue, meshkernel::constants::missing::doubleValue}};

    auto polygons = std::make_shared<meshkernel::Polygons>();

    // Execute
    auto landboundaries = std::make_shared<meshkernel::LandBoundaries>(landBoundaryPolygon, mesh, polygons);
    landboundaries->FindNearestMeshBoundary(meshkernel::LandBoundaries::ProjectToLandBoundaryOption ::OuterMeshBoundaryToLandBoundary);

    // Checks
    EXPECT_EQ(1, landboundaries->m_meshNodesLandBoundarySegments[0]);
    EXPECT_EQ(0, landboundaries->m_meshNodesLandBoundarySegments[1]);
    EXPECT_EQ(2, landboundaries->m_meshNodesLandBoundarySegments[2]);
    EXPECT_EQ(2, landboundaries->m_meshNodesLandBoundarySegments[3]);
    EXPECT_EQ(3, landboundaries->m_meshNodesLandBoundarySegments[4]);
    EXPECT_EQ(1, landboundaries->m_meshNodesLandBoundarySegments[5]);
    EXPECT_EQ(1, landboundaries->m_meshNodesLandBoundarySegments[6]);
    EXPECT_EQ(meshkernel::constants::missing::uintValue, landboundaries->m_meshNodesLandBoundarySegments[7]);
    EXPECT_EQ(meshkernel::constants::missing::uintValue, landboundaries->m_meshNodesLandBoundarySegments[8]);
    EXPECT_EQ(meshkernel::constants::missing::uintValue, landboundaries->m_meshNodesLandBoundarySegments[9]);
}

TEST(LandBoundaries, TwoLandBoundaries)
{
    // Prepare
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/SmallTriangularGrid_net.nc");

    std::vector<meshkernel::Point> landBoundaryPolygon{
        {222.621918, 382.651917},
        {316.206177, 461.190796},
        {350.811279, 465.102692},
        {510.295715, 438.923065},
        {meshkernel::constants::missing::doubleValue, meshkernel::constants::missing::doubleValue},
        {215.980743, 363.986420},
        {250.253036, 235.233246},
        {423.158325, 200.652054},
        {559.938782, 312.732147},
        {518.873718, 421.415894},
        {meshkernel::constants::missing::doubleValue, meshkernel::constants::missing::doubleValue}};

    auto polygons = std::make_shared<meshkernel::Polygons>();

    // Execute
    auto landboundaries = std::make_shared<meshkernel::LandBoundaries>(landBoundaryPolygon, mesh, polygons);
    landboundaries->FindNearestMeshBoundary(meshkernel::LandBoundaries::ProjectToLandBoundaryOption::OuterMeshBoundaryToLandBoundary);

    // Checks
    EXPECT_EQ(2, landboundaries->m_meshNodesLandBoundarySegments[0]);
    EXPECT_EQ(1, landboundaries->m_meshNodesLandBoundarySegments[1]);
    EXPECT_EQ(1, landboundaries->m_meshNodesLandBoundarySegments[2]);
    EXPECT_EQ(3, landboundaries->m_meshNodesLandBoundarySegments[3]);
    EXPECT_EQ(3, landboundaries->m_meshNodesLandBoundarySegments[4]);
    EXPECT_EQ(2, landboundaries->m_meshNodesLandBoundarySegments[5]);
    EXPECT_EQ(2, landboundaries->m_meshNodesLandBoundarySegments[6]);
    EXPECT_EQ(meshkernel::constants::missing::uintValue, landboundaries->m_meshNodesLandBoundarySegments[7]);
    EXPECT_EQ(meshkernel::constants::missing::uintValue, landboundaries->m_meshNodesLandBoundarySegments[8]);
    EXPECT_EQ(meshkernel::constants::missing::uintValue, landboundaries->m_meshNodesLandBoundarySegments[9]);
}

TEST(LandBoundaries, OneCrossingLandBoundary)
{
    // Prepare
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/SmallTriangularGrid_net.nc");

    std::vector<meshkernel::Point> landBoundaryPolygon{
        {221.418243, 315.848755},
        {248.801422, 375.129028},
        {337.571045, 459.686218},
        {516.313965, 419.062683},
        {528.651428, 292.377380},
        {meshkernel::constants::missing::doubleValue, meshkernel::constants::missing::doubleValue}};

    auto polygons = std::make_shared<meshkernel::Polygons>();

    // Execute
    auto landboundaries = std::make_shared<meshkernel::LandBoundaries>(landBoundaryPolygon, mesh, polygons);
    landboundaries->FindNearestMeshBoundary(meshkernel::LandBoundaries::ProjectToLandBoundaryOption ::OuterMeshBoundaryToLandBoundary);

    // Checks
    EXPECT_EQ(0, landboundaries->m_meshNodesLandBoundarySegments[0]);
    EXPECT_EQ(0, landboundaries->m_meshNodesLandBoundarySegments[1]);
    EXPECT_EQ(2, landboundaries->m_meshNodesLandBoundarySegments[2]);
    EXPECT_EQ(2, landboundaries->m_meshNodesLandBoundarySegments[3]);
    EXPECT_EQ(1, landboundaries->m_meshNodesLandBoundarySegments[4]);
    EXPECT_EQ(1, landboundaries->m_meshNodesLandBoundarySegments[5]);
    EXPECT_EQ(1, landboundaries->m_meshNodesLandBoundarySegments[6]);
    EXPECT_EQ(meshkernel::constants::missing::uintValue, landboundaries->m_meshNodesLandBoundarySegments[7]);
    EXPECT_EQ(meshkernel::constants::missing::uintValue, landboundaries->m_meshNodesLandBoundarySegments[8]);
    EXPECT_EQ(meshkernel::constants::missing::uintValue, landboundaries->m_meshNodesLandBoundarySegments[9]);
}

TEST(LandBoundaries, TwoCrossingLandBoundary)
{
    // Prepare
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/SmallTriangularGrid_net.nc");

    std::vector<meshkernel::Point> landBoundaryPolygon{
        {235.561218, 290.571899},
        {265.953522, 436.515747},
        {429.349854, 450.959656},
        {535.271545, 386.262909},
        {meshkernel::constants::missing::doubleValue, meshkernel::constants::missing::doubleValue},
        {246.995941, 262.285858},
        {351.112183, 237.309906},
        {443.191895, 262.285858},
        {553.627319, 327.283539},
        {meshkernel::constants::missing::doubleValue, meshkernel::constants::missing::doubleValue}};

    auto polygons = std::make_shared<meshkernel::Polygons>();

    // Execute
    auto landboundaries = std::make_shared<meshkernel::LandBoundaries>(landBoundaryPolygon, mesh, polygons);
    landboundaries->FindNearestMeshBoundary(meshkernel::LandBoundaries::ProjectToLandBoundaryOption ::OuterMeshBoundaryToLandBoundary);

    // Checks
    EXPECT_EQ(2, landboundaries->m_meshNodesLandBoundarySegments[0]);
    EXPECT_EQ(0, landboundaries->m_meshNodesLandBoundarySegments[1]);
    EXPECT_EQ(1, landboundaries->m_meshNodesLandBoundarySegments[2]);
    EXPECT_EQ(3, landboundaries->m_meshNodesLandBoundarySegments[3]);
    EXPECT_EQ(3, landboundaries->m_meshNodesLandBoundarySegments[4]);
    EXPECT_EQ(2, landboundaries->m_meshNodesLandBoundarySegments[5]);
    EXPECT_EQ(2, landboundaries->m_meshNodesLandBoundarySegments[6]);
    EXPECT_EQ(meshkernel::constants::missing::uintValue, landboundaries->m_meshNodesLandBoundarySegments[7]);
    EXPECT_EQ(meshkernel::constants::missing::uintValue, landboundaries->m_meshNodesLandBoundarySegments[8]);
    EXPECT_EQ(meshkernel::constants::missing::uintValue, landboundaries->m_meshNodesLandBoundarySegments[9]);
}

TEST(LandBoundaries, LandBoundaryConstructorTestSinglePolyline)
{
    std::vector<meshkernel::Point> controlPoints{{235.561218, 290.571899},
                                                 {265.953522, 436.515747},
                                                 {429.349854, 450.959656},
                                                 {535.271545, 386.262909}};

    meshkernel::LandBoundary landBoundary(controlPoints);

    EXPECT_EQ(landBoundary.GetNumNodes(), controlPoints.size());
    EXPECT_FALSE(landBoundary.IsEmpty());

    // Check the points are the same as those used during construction
    for (size_t i = 0; i < landBoundary.GetNumNodes(); ++i)
    {
        EXPECT_EQ(landBoundary.Node(i).x, controlPoints[i].x);
        EXPECT_EQ(landBoundary.Node(i).y, controlPoints[i].y);
    }

    // Check the point array values are the same as those used during construction
    const std::vector<meshkernel::Point>& points(landBoundary.GetNodes());

    for (size_t i = 0; i < landBoundary.GetNumNodes(); ++i)
    {
        EXPECT_EQ(points[i].x, controlPoints[i].x);
        EXPECT_EQ(points[i].y, controlPoints[i].y);
    }
}

TEST(LandBoundaries, LandBoundaryConstructorTestMultiPolyline)
{
    std::vector<meshkernel::Point> controlPoints{{235.561218, 290.571899},
                                                 {265.953522, 436.515747},
                                                 {429.349854, 450.959656},
                                                 {535.271545, 386.262909},
                                                 {meshkernel::constants::missing::doubleValue, meshkernel::constants::missing::doubleValue},
                                                 {246.995941, 262.285858},
                                                 {351.112183, 237.309906},
                                                 {443.191895, 262.285858},
                                                 {553.627319, 327.283539},
                                                 {meshkernel::constants::missing::doubleValue, meshkernel::constants::missing::doubleValue}};

    meshkernel::LandBoundary landBoundary(controlPoints);

    std::vector<std::pair<meshkernel::UInt, meshkernel::UInt>> indices = landBoundary.FindPolylineIndices();

    EXPECT_EQ(indices.size(), 2);

    // Indices of the first polyline
    EXPECT_EQ(indices[0].first, 0);
    EXPECT_EQ(indices[0].second, 3);

    // Indices of the second polyline
    EXPECT_EQ(indices[1].first, 5);
    EXPECT_EQ(indices[1].second, 8);
}

TEST(LandBoundaries, LandBoundaryConstructorTestFindClosestPoint)
{
    std::vector<meshkernel::Point> controlPoints{{0.0, 0.0},
                                                 {10.0, 0.0},
                                                 {10.0, 10.0}};

    meshkernel::LandBoundary landBoundary(controlPoints);
    meshkernel::Projection::Type projection = meshkernel::Projection::Cartesian;

    // Expect  this point to be closest to the first point (0, 0)
    meshkernel::Point samplePoint({1.0, 1.0});
    meshkernel::Point closestPoint = landBoundary.ClosestPoint(samplePoint, 0, 1, projection);
    EXPECT_EQ(closestPoint.x, controlPoints[0].x);
    EXPECT_EQ(closestPoint.y, controlPoints[0].y);

    // Expect  this point to be closest to the first point (10, 0)
    samplePoint = meshkernel::Point({11.0, 1.0});
    closestPoint = landBoundary.ClosestPoint(samplePoint, 0, 1, projection);
    EXPECT_EQ(closestPoint.x, controlPoints[1].x);
    EXPECT_EQ(closestPoint.y, controlPoints[1].y);

    // Expect  this point to be closest to the first point (10, 0)
    samplePoint = meshkernel::Point({11.0, 10.0});
    closestPoint = landBoundary.ClosestPoint(samplePoint, 0, 1, projection);
    EXPECT_EQ(closestPoint.x, controlPoints[1].x);
    EXPECT_EQ(closestPoint.y, controlPoints[1].y);

    // Expect  this point to be closest to the first point (10, 0)
    samplePoint = meshkernel::Point({11.0, 1.0});
    closestPoint = landBoundary.ClosestPoint(samplePoint, 1, 2, projection);
    EXPECT_EQ(closestPoint.x, controlPoints[1].x);
    EXPECT_EQ(closestPoint.y, controlPoints[1].y);

    // Expect  this point to be closest to the first point (10, 10)
    samplePoint = meshkernel::Point({9.0, 5.1});
    closestPoint = landBoundary.ClosestPoint(samplePoint, 1, 2, projection);
    EXPECT_EQ(closestPoint.x, controlPoints[2].x);
    EXPECT_EQ(closestPoint.y, controlPoints[2].y);

    // Expect  this point to be closest to the first point (10, 10)
    samplePoint = meshkernel::Point({11.0, 5.1});
    closestPoint = landBoundary.ClosestPoint(samplePoint, 1, 2, projection);
    EXPECT_EQ(closestPoint.x, controlPoints[2].x);
    EXPECT_EQ(closestPoint.y, controlPoints[2].y);
}

TEST(LandBoundaries, LandBoundaryConstructorTestAddSegment)
{
    std::vector<meshkernel::Point> controlPoints{{235.561218, 290.571899},
                                                 {265.953522, 436.515747},
                                                 {429.349854, 450.959656},
                                                 {535.271545, 386.262909}};

    std::vector<meshkernel::Point> controlPointsAfterAddition{{235.561218, 290.571899},
                                                              {265.953522, 436.515747},
                                                              {429.349854, 450.959656},
                                                              {535.271545, 386.262909},
                                                              {meshkernel::constants::missing::doubleValue, meshkernel::constants::missing::doubleValue},
                                                              {246.995941, 262.285858},
                                                              {351.112183, 237.309906},
                                                              {meshkernel::constants::missing::doubleValue, meshkernel::constants::missing::doubleValue}};

    meshkernel::LandBoundary landBoundary(controlPoints);

    EXPECT_EQ(landBoundary.GetNumNodes(), controlPoints.size());
    EXPECT_FALSE(landBoundary.IsEmpty());

    // Check the points are the same as those used during construction
    for (size_t i = 0; i < landBoundary.GetNumNodes(); ++i)
    {
        EXPECT_EQ(landBoundary.Node(i).x, controlPoints[i].x);
        EXPECT_EQ(landBoundary.Node(i).y, controlPoints[i].y);
    }

    meshkernel::Point p1{246.995941, 262.285858};
    meshkernel::Point p2{351.112183, 237.309906};

    // Now add the new segment.
    landBoundary.AddSegment(p1, p2);

    EXPECT_EQ(landBoundary.GetNumNodes(), controlPointsAfterAddition.size());
    EXPECT_FALSE(landBoundary.IsEmpty());

    // Check the points are the same as those used during construction
    for (size_t i = 0; i < landBoundary.GetNumNodes(); ++i)
    {
        EXPECT_EQ(landBoundary.Node(i).x, controlPointsAfterAddition[i].x);
        EXPECT_EQ(landBoundary.Node(i).y, controlPointsAfterAddition[i].y);
    }
}
