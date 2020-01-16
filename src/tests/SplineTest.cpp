#include "GridGeomTest.hpp"
#include "../Splines.hpp"
#include "../Entities.hpp"

TEST(SplineTests, SetSpline)
{
    //One gets the edges
    std::vector<GridGeom::Point> splineNodes;
    splineNodes.push_back(GridGeom::Point{GridGeom::doubleMissingValue, GridGeom::doubleMissingValue });
    splineNodes.push_back(GridGeom::Point{ GridGeom::doubleMissingValue, GridGeom::doubleMissingValue });
    splineNodes.push_back(GridGeom::Point{ GridGeom::doubleMissingValue, GridGeom::doubleMissingValue });
    splineNodes.push_back(GridGeom::Point{ GridGeom::doubleMissingValue, GridGeom::doubleMissingValue });

    GridGeom::Splines splines(GridGeom::Projections::cartesian);
    bool success = splines.Set(splineNodes);

    ASSERT_TRUE(success);
    EXPECT_EQ(splines.m_numSplines, 1);
    EXPECT_EQ(splines.m_splines[0].size(), splineNodes.size());
    EXPECT_EQ(splines.m_numAllocatedSplines, 5);
    EXPECT_EQ(splines.m_numAllocatedSplineNodes[0], 10);
}

TEST(SplineTests, CubicSplineInterpolation)
{
    //One gets the edges
    std::vector<GridGeom::Point> splineNodes;
    
    splineNodes.push_back(GridGeom::Point{ 212.001953125000, 155.627197265625 });
    splineNodes.push_back(GridGeom::Point{ 529.253906250000, 432.379974365234 });
    splineNodes.push_back(GridGeom::Point{ 930.506469726562, 453.380187988281 });

    int pointsBetweenVertices = 20;
    std::vector<GridGeom::Point> coordinatesDerivatives(splineNodes.size());
    GridGeom::Splines::SecondOrderDerivative(splineNodes, splineNodes.size(), coordinatesDerivatives);
    std::vector<GridGeom::Point> splineCoordinates;

    for (int n = 0; n < splineNodes.size() - 1; n++)
    {
        for (int p = 0; p <= pointsBetweenVertices; p++)
        {
            const double pointAdimensionalCoordinate = n + double(p) / double(pointsBetweenVertices);
            GridGeom::Point pointCoordinate;
            GridGeom::Splines::Interpolate(splineNodes, coordinatesDerivatives, pointAdimensionalCoordinate, pointCoordinate);
            splineCoordinates.push_back({ pointCoordinate.x, pointCoordinate.y });
        }
    }

    const double tolerance = 1e-3;
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

TEST(SplineTests, SplineIntersection)
{
    std::vector<GridGeom::Point> firstSpline;
    firstSpline.push_back(GridGeom::Point{ 152.001571655273, 86.6264953613281 });
    firstSpline.push_back(GridGeom::Point{ 374.752960205078, 336.378997802734 });
    firstSpline.push_back(GridGeom::Point{ 850.255920410156, 499.130676269531 });
    GridGeom::Splines splines(GridGeom::Projections::cartesian);
    bool success = splines.Set(firstSpline);

    std::vector<GridGeom::Point> secondSpline;
    secondSpline.push_back(GridGeom::Point{ 72.5010681152344,391.129577636719 });
    secondSpline.push_back(GridGeom::Point{ 462.503479003906, 90.3765411376953 });
    success = splines.Set(secondSpline);

    double crossProductIntersection;
    GridGeom::Point dimensionalIntersection;
    double firstSplineRatio;
    double secondSplineRatio;

    splines.GetSplinesIntersection(0, 1, GridGeom::Projections::cartesian,
        crossProductIntersection, dimensionalIntersection, firstSplineRatio, secondSplineRatio);

    const double tolerance = 1e-5;
    ASSERT_NEAR(261.736770097059, dimensionalIntersection.x, tolerance);
    ASSERT_NEAR(245.199166962145, dimensionalIntersection.y, tolerance);
    ASSERT_NEAR(0.601498208554790, firstSplineRatio, tolerance);
    ASSERT_NEAR(0.485216749175026, secondSplineRatio, tolerance);
    ASSERT_NEAR(-0.996215079635043, crossProductIntersection, tolerance);
}

TEST(SplineTests, ComputeSplinesProperties)
{
    std::vector<GridGeom::Point> firstSpline;
    firstSpline.push_back(GridGeom::Point{ 152.001571655273, 86.6264953613281 });
    firstSpline.push_back(GridGeom::Point{ 374.752960205078, 336.378997802734 });
    firstSpline.push_back(GridGeom::Point{ 850.255920410156, 499.130676269531 });
    GridGeom::Splines splines(GridGeom::Projections::cartesian);
    bool success = splines.Set(firstSpline);
    ASSERT_EQ(true, success);

    std::vector<GridGeom::Point> secondSpline;
    secondSpline.push_back(GridGeom::Point{ 72.5010681152344,391.129577636719 });
    secondSpline.push_back(GridGeom::Point{ 462.503479003906, 90.3765411376953 });
    success = splines.Set(secondSpline);

    success = splines.ComputeSplineProperties();
    ASSERT_EQ(true, success);
    
    const double tolerance = 1e-3;
    ASSERT_NEAR(253.52971595547601, splines.m_maximumGridHeight[0], tolerance);
    ASSERT_NEAR(0.0, splines.m_maximumGridHeight[1], tolerance);
}

TEST(SplineTests, ComputeBoundingBox)
{
    std::vector<GridGeom::Point> firstSpline;
    firstSpline.push_back(GridGeom::Point{ 91.2511901855469, 299.628631591797 });
    firstSpline.push_back(GridGeom::Point{ 354.502838134766, 518.630859375000 });
    firstSpline.push_back(GridGeom::Point{ 770.755432128906, 607.881774902344 });
    GridGeom::Splines splines(GridGeom::Projections::cartesian);
    bool success = splines.Set(firstSpline);
    ASSERT_EQ(true, success);

    std::vector<GridGeom::Point> secondSpline;
    secondSpline.push_back(GridGeom::Point{ 273.502319335938, 86.6264953613281 });
    secondSpline.push_back(GridGeom::Point{ 557.004089355469, 316.128814697266 });
    secondSpline.push_back(GridGeom::Point{ 847.255920410156, 409.129730224609 });
    success = splines.Set(secondSpline);
    ASSERT_EQ(true, success);

    std::vector<GridGeom::Point> thirdpline;
    thirdpline.push_back(GridGeom::Point{ 62.7510070800781, 396.379608154297 });
    thirdpline.push_back(GridGeom::Point{ 350.752807617188, 73.8763732910156 });
    success = splines.Set(thirdpline);
    ASSERT_EQ(true, success);
    
    std::vector<GridGeom::Point> fourthSpline;
    fourthSpline.push_back(GridGeom::Point{ 704.755004882812, 636.382019042969 });
    fourthSpline.push_back(GridGeom::Point{ 845.005859375000, 285.378509521484 });
    success = splines.Set(fourthSpline);
    ASSERT_EQ(true, success);

    success = splines.ComputeSplineProperties();
    ASSERT_EQ(true, success);

    const double tolerance = 1e-2;
    ASSERT_NEAR(345.967070088532, splines.m_maximumGridHeight[0], tolerance);
    ASSERT_NEAR(370.339417298715, splines.m_maximumGridHeight[1], tolerance);
    ASSERT_NEAR(0.0, splines.m_maximumGridHeight[2], tolerance);
    ASSERT_NEAR(0.0, splines.m_maximumGridHeight[3], tolerance);
}