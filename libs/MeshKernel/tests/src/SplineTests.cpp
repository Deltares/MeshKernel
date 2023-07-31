#include <gtest/gtest.h>

#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/LandBoundary.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/SplineAlgorithms.hpp>
#include <MeshKernel/Splines.hpp>

TEST(Splines, SetSpline)
{
    // a valid spline
    std::vector<meshkernel::Point> splineNodes({{212.001953125000, 155.627197265625},
                                                {529.253906250000, 432.379974365234},
                                                {930.506469726562, 453.380187988281},
                                                {1030.506469726562, 653.380187988281}});

    meshkernel::Splines splines(meshkernel::Projection::cartesian);
    splines.AddSpline(splineNodes, 0, static_cast<meshkernel::UInt>(splineNodes.size()));

    ASSERT_EQ(1, splines.GetNumSplines());
    ASSERT_EQ(4, splines.m_splineNodes[0].size());
}

TEST(Splines, CubicSplineInterpolation)
{
    // One gets the edges
    std::vector<meshkernel::Point> splineNodes;
    splineNodes.push_back(meshkernel::Point{212.001953125000, 155.627197265625});
    splineNodes.push_back(meshkernel::Point{529.253906250000, 432.379974365234});
    splineNodes.push_back(meshkernel::Point{930.506469726562, 453.380187988281});

    int pointsBetweenNodes = 20;
    auto coordinatesDerivatives = meshkernel::SplineAlgorithms::SecondOrderDerivative(splineNodes, 0, static_cast<meshkernel::UInt>(splineNodes.size()) - 1);
    std::vector<meshkernel::Point> splineCoordinates;

    for (size_t n = 0; n < splineNodes.size() - 1; n++)
    {
        for (auto p = 0; p <= pointsBetweenNodes; p++)
        {
            const double pointAdimensionalCoordinate = n + double(p) / double(pointsBetweenNodes);
            auto pointCoordinate = ComputePointOnSplineAtAdimensionalDistance(splineNodes, coordinatesDerivatives, pointAdimensionalCoordinate);
            ASSERT_TRUE(pointCoordinate.IsValid());

            splineCoordinates.push_back({pointCoordinate.x, pointCoordinate.y});
        }
    }

    const double tolerance = 1e-6;
    ASSERT_NEAR(226.817168170929, splineCoordinates[1].x, tolerance);
    ASSERT_NEAR(241.648133331299, splineCoordinates[2].x, tolerance);
    ASSERT_NEAR(256.510598720551, splineCoordinates[3].x, tolerance);
    ASSERT_NEAR(271.420314453125, splineCoordinates[4].x, tolerance);
    ASSERT_NEAR(286.393030643463, splineCoordinates[5].x, tolerance);
    ASSERT_NEAR(930.506469726562, splineCoordinates.back().x, tolerance);

    ASSERT_NEAR(172.653750896454, splineCoordinates[1].y, tolerance);
    ASSERT_NEAR(189.632350921631, splineCoordinates[2].y, tolerance);
    ASSERT_NEAR(206.515043735504, splineCoordinates[3].y, tolerance);
    ASSERT_NEAR(223.253875732422, splineCoordinates[4].y, tolerance);
    ASSERT_NEAR(453.380187988281, splineCoordinates.back().y, tolerance);
}

TEST(Splines, SplineIntersection)
{
    std::vector<meshkernel::Point> firstSpline;
    firstSpline.push_back(meshkernel::Point{152.001571655273, 86.6264953613281});
    firstSpline.push_back(meshkernel::Point{374.752960205078, 336.378997802734});
    firstSpline.push_back(meshkernel::Point{850.255920410156, 499.130676269531});

    meshkernel::Splines splines(meshkernel::Projection::cartesian);

    splines.AddSpline(firstSpline, 0, static_cast<meshkernel::UInt>(firstSpline.size()));

    std::vector<meshkernel::Point> secondSpline;
    secondSpline.push_back(meshkernel::Point{72.5010681152344, 391.129577636719});
    secondSpline.push_back(meshkernel::Point{462.503479003906, 90.3765411376953});
    splines.AddSpline(secondSpline, 0, static_cast<meshkernel::UInt>(secondSpline.size()));

    double crossProductIntersection;
    meshkernel::Point dimensionalIntersection;
    double firstSplineRatio;
    double secondSplineRatio;

    splines.GetSplinesIntersection(0, 1, crossProductIntersection, dimensionalIntersection, firstSplineRatio, secondSplineRatio);

    const double tolerance = 1e-6;
    ASSERT_NEAR(261.736770097059, dimensionalIntersection.x, tolerance);
    ASSERT_NEAR(245.199166962145, dimensionalIntersection.y, tolerance);
    ASSERT_NEAR(0.601498208554790, firstSplineRatio, tolerance);
    ASSERT_NEAR(0.485216749175026, secondSplineRatio, tolerance);
    ASSERT_NEAR(-0.996215079635043, crossProductIntersection, tolerance);
}

TEST(Spline, SnapToLandBoundaryTest)
{
    // Test the algorithm for snapping splines to land boundaries.

    constexpr double tolerance = 1.0e-8;

    constexpr int iterationCountForTest = 5;

    // The land boundary to which the spline is to be snapped.
    std::vector<meshkernel::Point> landBoundaryPoints{{257.002197, 442.130066},
                                                      {518.753845, 301.128662},
                                                      {938.006470, 416.629822}};

    meshkernel::LandBoundary landBoundary(landBoundaryPoints);

    // The original spline points.
    std::vector<meshkernel::Point> splinePoints{{281.0023, 447.3801},
                                                {367.2529, 401.6296},
                                                {461.7534, 354.3792},
                                                {517.2538, 318.3788},
                                                {614.0045, 338.629},
                                                {720.5051, 377.6294},
                                                {827.7558, 417.3798},
                                                {923.7563, 424.1299}};

    // The expected spline values after snapping to land boundary.
    std::vector<meshkernel::Point> expectedSplinePoints{{273.5868719643935, 434.2730022174478},
                                                        {359.5998304717778, 386.1712239047134},
                                                        {451.5303458337523, 338.3551703843473},
                                                        {517.7962262926076, 306.3259738916997},
                                                        {616.7325138813335, 327.9627689164845},
                                                        {725.7358644094627, 358.0902879743862},
                                                        {836.262785315633, 388.6415116416172},
                                                        {923.500177844106, 412.5818685325169}};

    meshkernel::Splines splines(meshkernel::Projection::cartesian);
    splines.AddSpline(splinePoints, 0, static_cast<meshkernel::UInt>(splinePoints.size()));

    // Snap the spline to the land boundary
    splines.SnapSpline(0, landBoundary, iterationCountForTest);

    // The number of splines should be unchanged.
    ASSERT_EQ(splines.GetNumSplines(), 1) << ", expected number of splines to be unchenged, there should be only 1 spline";

    // The number of points in the spline should be unchanged
    ASSERT_EQ(splines.m_splineNodes[0].size(), expectedSplinePoints.size()) << ", expected the number of points to be unchanged.";

    // Now test the snapped spline points are close to the expected values
    for (size_t i = 0; i < splinePoints.size(); ++i)
    {
        EXPECT_TRUE(meshkernel::IsEqual(expectedSplinePoints[i].x, splines.m_splineNodes[0][i].x, tolerance))
            << "Expected x-value: " << expectedSplinePoints[i].x << ", actual: " << splines.m_splineNodes[0][i].x << ", relative tolerance: " << tolerance;
        EXPECT_TRUE(meshkernel::IsEqual(expectedSplinePoints[i].y, splines.m_splineNodes[0][i].y, tolerance))
            << "Expected y-value: " << expectedSplinePoints[i].y << ", actual: " << splines.m_splineNodes[0][i].y << ", relative tolerance: " << tolerance;
    }
}

TEST(Spline, SplineSnappingIndexOutOfRangeTest)
{
    // Test the algorithm for snapping splines to land boundaries.

    // The land boundary to which the spline is to be snapped.
    std::vector<meshkernel::Point> landBoundaryPoints{{257.002197, 442.130066},
                                                      {518.753845, 301.128662},
                                                      {938.006470, 416.629822}};

    meshkernel::LandBoundary landBoundary(landBoundaryPoints);

    // The original spline points.
    std::vector<meshkernel::Point> firstSpline{{281.0023, 447.3801},
                                               {367.2529, 401.6296},
                                               {461.7534, 354.3792},
                                               {517.2538, 318.3788},
                                               {614.0045, 338.629},
                                               {720.5051, 377.6294},
                                               {827.7558, 417.3798},
                                               {923.7563, 424.1299}};

    std::vector<meshkernel::Point> secondSpline({{212.001953125000, 155.627197265625},
                                                 {529.253906250000, 432.379974365234},
                                                 {930.506469726562, 453.380187988281},
                                                 {1030.506469726562, 653.380187988281}});

    meshkernel::Splines splines(meshkernel::Projection::cartesian);
    splines.AddSpline(firstSpline, 0, static_cast<meshkernel::UInt>(firstSpline.size()));
    splines.AddSpline(secondSpline, 0, static_cast<meshkernel::UInt>(secondSpline.size()));

    // Spline 2 is out of range, so should throw an exception.
    EXPECT_THROW(splines.SnapSpline(2, landBoundary), meshkernel::ConstraintError);
}
