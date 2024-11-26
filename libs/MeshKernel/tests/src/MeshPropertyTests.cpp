#include <chrono>
#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>
#include <random>
#include <span>

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/MeshTriangulation.hpp"
#include "MeshKernel/SampleInterpolator.hpp"

namespace mk = meshkernel;

TEST(MeshPropertyTests, TriangulationTest)
{
    // const double tolerance = 1.0e-13;
    std::vector<double> xValues{0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0};
    std::vector<double> yValues{0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0};

    std::vector<mk::Point> nodes{{0.25, 0.25}, {1.25, 0.25}, {2.25, 0.25}, {0.25, 1.25}, {1.25, 1.25}, {2.25, 1.25}, {0.25, 2.25}, {1.25, 2.25}, {2.25, 2.25}};
    std::vector<double> interpolated(nodes.size(), 0.0);

    mk::MeshTriangulation triangulation(xValues, yValues, mk::Projection::cartesian);

    EXPECT_EQ (triangulation.GetProjection (), mk::Projection::cartesian);
    ASSERT_EQ (triangulation.NumberOfNodes(), 16);
    ASSERT_EQ (triangulation.NumberOfEdges(), 33);
    ASSERT_EQ (triangulation.NumberOfFaces(), 18);

    std::vector<mk::Edge> expectedEdgeIds{{4, 0}, {0, 1}, {1, 4}, {1, 5}, {5, 4}, {1, 2}, {2, 5}, {5, 9},
                                          {9, 4}, {12, 8}, {8, 9}, {9, 12}, {9, 13}, {13, 12}, {9, 10},
                                          {10, 13}, {8, 4}, {5, 6}, {6, 9}, {2, 6}, {2, 3}, {3, 6}, {3, 7},
                                          {7, 6}, {7, 11}, {11, 6}, {10, 14}, {14, 13}, {10, 11}, {11, 14},
                                          {11, 15}, {15, 14}, {10, 6}};
    std::vector<std::array<mk::UInt, 3>> expectedTriangleNodeIds{{4, 0, 1}, {1, 5, 4}, {5, 1, 2},
                                                                 {4, 5, 9}, {12, 8, 9}, {9, 13, 12},
                                                                 {13, 9, 10}, {8, 4, 9}, {5, 6, 9},
                                                                 {2, 6, 5}, {6, 2, 3}, {3, 7, 6},
                                                                 {6, 7, 11}, {10, 14, 13}, {14, 10, 11},
                                                                 {11, 15, 14}, {10, 6, 11}, {6, 10, 9}};


    for (mk::UInt i = 0; i < triangulation.NumberOfEdges(); ++i) {
        auto edge = triangulation.GetEdge (i);
        EXPECT_EQ (edge.first, expectedEdgeIds [i].first);
        EXPECT_EQ (edge.second, expectedEdgeIds [i].second);
    }


    for (mk::UInt i = 0; i < triangulation.NumberOfFaces(); ++i) {
        auto face = triangulation.GetNodeIds (i);
        EXPECT_EQ (face [0], expectedTriangleNodeIds [i][0]);
        EXPECT_EQ (face [1], expectedTriangleNodeIds [i][1]);
        EXPECT_EQ (face [2], expectedTriangleNodeIds [i][2]);
    }

    for (mk::UInt i = 0; i < triangulation.NumberOfFaces(); ++i) {
        auto faceNodes = triangulation.GetNodes (i);
        EXPECT_EQ (faceNodes [0].x, triangulation.GetNode (expectedTriangleNodeIds [i][0]).x);
        EXPECT_EQ (faceNodes [0].y, triangulation.GetNode (expectedTriangleNodeIds [i][0]).y);

        EXPECT_EQ (faceNodes [1].x, triangulation.GetNode (expectedTriangleNodeIds [i][1]).x);
        EXPECT_EQ (faceNodes [1].y, triangulation.GetNode (expectedTriangleNodeIds [i][1]).y);

        EXPECT_EQ (faceNodes [2].x, triangulation.GetNode (expectedTriangleNodeIds [i][2]).x);
        EXPECT_EQ (faceNodes [2].y, triangulation.GetNode (expectedTriangleNodeIds [i][2]).y);
    }

    // Origin of mesh is {0, 0}, with delta {1,1} and there are 3x3 nodes

    std::vector<mk::Point> pointsToCheck{{-1.0, -1.0}, {0.25, 0.25}, {0.25, 2.25}, {3.5, 0.25}, {4.0, 4.0}, {2.5, 2.5}};
    std::vector<mk::UInt> expectedElementId{0, 0, 4, 11, 15, 15};
    std::vector<bool> expectedIsInElement{false, true, true, false, false, true};

    for (size_t i = 0; i < pointsToCheck.size (); ++i)
    {
        mk::UInt faceId = triangulation.FindNearestFace (pointsToCheck [i]);
        EXPECT_EQ (faceId, expectedElementId[i]);
        bool isInFace = triangulation.PointIsInElement (pointsToCheck [i], faceId);
        EXPECT_EQ (isInFace, expectedIsInElement[i]);
    }

}

TEST(MeshPropertyTests, TriangulationFailureTest)
{
    // Not enough points
    EXPECT_THROW (mk::MeshTriangulation triangulation(std::vector<double>{0.0, 1.0}, std::vector<double>{0.0, 1.0}, mk::Projection::cartesian), mk::ConstraintError);
    // x- and y-points not same size
    EXPECT_THROW (mk::MeshTriangulation triangulation(std::vector<double>{0.0, 1.0, 2.0}, std::vector<double>{0.0, 1.0}, mk::Projection::cartesian), mk::ConstraintError);

    std::vector<double> xValues{0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0};
    std::vector<double> yValues{0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0};
    mk::MeshTriangulation triangulation(xValues, yValues, mk::Projection::cartesian);

    // Invalid point passed to PointIsInElement
    EXPECT_THROW (triangulation.PointIsInElement({mk::constants::missing::doubleValue, mk::constants::missing::doubleValue}, 1), mk::ConstraintError);
    // Face id is null
    EXPECT_THROW (triangulation.PointIsInElement({0.25, 0.25}, mk::constants::missing::uintValue), mk::ConstraintError);
    // Face id is out of range
    EXPECT_THROW (triangulation.PointIsInElement({0.25, 0.25}, 100), mk::ConstraintError);

}

TEST(MeshPropertyTests, SimpleSampleInterpolationTest)
{
    const double tolerance = 1.0e-13;

    std::vector<double> xValues{0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0};
    std::vector<double> yValues{0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0};

    std::vector<mk::Point> nodes{{0.25, 0.25}, {1.25, 0.25}, {2.25, 0.25}, {0.25, 1.25}, {1.25, 1.25}, {2.25, 1.25}, {0.25, 2.25}, {1.25, 2.25}, {2.25, 2.25}};
    std::vector<double> interpolated(nodes.size(), 0.0);

    mk::SampleInterpolator properties(xValues, yValues, mk::Projection::cartesian);

    std::vector<double> expectedXDepths{0.25, 1.25, 2.25, 0.25, 1.25, 2.25, 0.25, 1.25, 2.25};
    std::vector<double> expectedYDepths{0.25, 0.25, 0.25, 1.25, 1.25, 1.25, 2.25, 2.25, 2.25};

    properties.SetData("xdepth", std::vector{0.0, 1.0, 2.0, 3.0,
                                             0.0, 1.0, 2.0, 3.0,
                                             0.0, 1.0, 2.0, 3.0,
                                             0.0, 1.0, 2.0, 3.0});

    properties.SetData("ydepth", std::vector{0.0, 0.0, 0.0, 0.0,
                                             1.0, 1.0, 1.0, 1.0,
                                             2.0, 2.0, 2.0, 2.0,
                                             3.0, 3.0, 3.0, 3.0});

    properties.Interpolate("xdepth", nodes, interpolated);

    for (size_t i = 0; i < interpolated.size(); ++i)
    {
        EXPECT_NEAR(interpolated[i], expectedXDepths[i], tolerance);
    }

    std::span interpolatedData(interpolated.data(), interpolated.size());

    properties.Interpolate("ydepth", nodes, interpolatedData);

    for (size_t i = 0; i < interpolated.size(); ++i)
    {
        EXPECT_NEAR(interpolatedData[i], expectedYDepths[i], tolerance);
    }
}

TEST(MeshPropertyTests, AveragePointTest)
{
    const double tolerance = 1.0e-13;

    std::vector<mk::Point> pnts{{1.0, 2.0}, {-999.0, -999.0}, {1.0, 2.0}, {1.0, 2.0}};
    mk::Point average = mk::ComputeAverageCoordinate(pnts, mk::Projection::sphericalAccurate);

    EXPECT_NEAR(average.x, 1.0, tolerance);
    EXPECT_NEAR(average.y, 2.0, tolerance);
}
