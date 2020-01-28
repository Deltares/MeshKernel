#include "GridGeomTest.hpp"
#include "../Splines.hpp"
#include "../Polygons.hpp"
#include "../Entities.hpp"
#include "../Constants.cpp"

TEST(SplineTests, SetSpline)
{
    //One gets the edges
    std::vector<GridGeom::Point> splineNodes;
    splineNodes.push_back(GridGeom::Point{ GridGeom::doubleMissingValue, GridGeom::doubleMissingValue });
    splineNodes.push_back(GridGeom::Point{ GridGeom::doubleMissingValue, GridGeom::doubleMissingValue });
    splineNodes.push_back(GridGeom::Point{ GridGeom::doubleMissingValue, GridGeom::doubleMissingValue });
    splineNodes.push_back(GridGeom::Point{ GridGeom::doubleMissingValue, GridGeom::doubleMissingValue });

    GridGeom::Polygons polygon;
    GridGeom::Splines splines(GridGeom::Projections::cartesian, polygon);
    bool success = splines.AddSpline(splineNodes);
    ASSERT_TRUE(success);

    ASSERT_EQ(1, splines.m_numSplines);
    ASSERT_EQ(10, splines.m_splineCornerPoints[0].size());
    ASSERT_EQ(5, splines.m_numAllocatedSplines);
    ASSERT_EQ(10, splines.m_numAllocatedSplineNodes[0]);
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

    GridGeom::Polygons polygon;
    GridGeom::Splines splines(GridGeom::Projections::cartesian, polygon);

    bool success = splines.AddSpline(firstSpline);
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> secondSpline;
    secondSpline.push_back(GridGeom::Point{ 72.5010681152344,391.129577636719 });
    secondSpline.push_back(GridGeom::Point{ 462.503479003906, 90.3765411376953 });
    success = splines.AddSpline(secondSpline);

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

    GridGeom::Polygons polygon;
    GridGeom::Splines splines(GridGeom::Projections::cartesian, polygon);
    bool success = splines.AddSpline(firstSpline);
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> secondSpline;
    secondSpline.push_back(GridGeom::Point{ 72.5010681152344,391.129577636719 });
    secondSpline.push_back(GridGeom::Point{ 462.503479003906, 90.3765411376953 });
    success = splines.AddSpline(secondSpline);

    success = splines.ComputeSplineProperties(false);
    ASSERT_TRUE(success);
    success = splines.MakeAllGridLines(true);
    ASSERT_TRUE(success);
    ASSERT_EQ(7, splines.m_numGridLines);

    const double tolerance = 1e-3;
    ASSERT_NEAR(253.52971595547601, splines.m_maximumGridHeights[0], tolerance);
    ASSERT_NEAR(0.0, splines.m_maximumGridHeights[1], tolerance);

    ASSERT_NEAR(152.001571655273, splines.m_gridLine[0].x, tolerance);
    ASSERT_NEAR(407.924702423872, splines.m_gridLine[1].x, tolerance);
    ASSERT_NEAR(850.250753009544, splines.m_gridLine[2].x, tolerance);
    ASSERT_NEAR(-999.000000000000, splines.m_gridLine[3].x, tolerance);
    ASSERT_NEAR(850.250753009544, splines.m_gridLine[4].x, tolerance);
    ASSERT_NEAR(407.924702423872, splines.m_gridLine[5].x, tolerance);
    ASSERT_NEAR(152.001571655273, splines.m_gridLine[6].x, tolerance);

    ASSERT_NEAR(86.6264953613281, splines.m_gridLine[0].y, tolerance);
    ASSERT_NEAR(354.562246561336, splines.m_gridLine[1].y, tolerance);
    ASSERT_NEAR(499.129323710654, splines.m_gridLine[2].y, tolerance);
    ASSERT_NEAR(-999.000000000000, splines.m_gridLine[3].y, tolerance);
    ASSERT_NEAR(499.129323710654, splines.m_gridLine[4].y, tolerance);
    ASSERT_NEAR(354.562246561336, splines.m_gridLine[5].y, tolerance);
    ASSERT_NEAR(86.6264953613281, splines.m_gridLine[6].y, tolerance);

    ASSERT_NEAR(4.214936859441928e-14, splines.m_gridLineDimensionalCoordinates[0], tolerance);
    ASSERT_NEAR(1.09068327269294, splines.m_gridLineDimensionalCoordinates[1], tolerance);
    ASSERT_NEAR(1.99999040748403, splines.m_gridLineDimensionalCoordinates[2], tolerance);
    ASSERT_NEAR(-999.000000000000, splines.m_gridLineDimensionalCoordinates[3], tolerance);
    ASSERT_NEAR(1.99999040748403, splines.m_gridLineDimensionalCoordinates[4], tolerance);
    ASSERT_NEAR(1.09068327269294, splines.m_gridLineDimensionalCoordinates[5], tolerance);
    ASSERT_NEAR(4.214936859441928e-14, splines.m_gridLineDimensionalCoordinates[6], tolerance);
}

TEST(SplineTests, ComputeBoundingBox)
{
    std::vector<GridGeom::Point> firstSpline;
    firstSpline.push_back(GridGeom::Point{ 91.2511901855469, 299.628631591797 });
    firstSpline.push_back(GridGeom::Point{ 354.502838134766, 518.630859375000 });
    firstSpline.push_back(GridGeom::Point{ 770.755432128906, 607.881774902344 });

    GridGeom::Polygons polygon;
    GridGeom::Splines splines(GridGeom::Projections::cartesian, polygon);
    bool success = splines.AddSpline(firstSpline);
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> secondSpline;
    secondSpline.push_back(GridGeom::Point{ 273.502319335938, 86.6264953613281 });
    secondSpline.push_back(GridGeom::Point{ 557.004089355469, 316.128814697266 });
    secondSpline.push_back(GridGeom::Point{ 847.255920410156, 409.129730224609 });
    success = splines.AddSpline(secondSpline);
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> thirdpline;
    thirdpline.push_back(GridGeom::Point{ 62.7510070800781, 396.379608154297 });
    thirdpline.push_back(GridGeom::Point{ 350.752807617188, 73.8763732910156 });
    success = splines.AddSpline(thirdpline);
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> fourthSpline;
    fourthSpline.push_back(GridGeom::Point{ 704.755004882812, 636.382019042969 });
    fourthSpline.push_back(GridGeom::Point{ 845.005859375000, 285.378509521484 });
    success = splines.AddSpline(fourthSpline);
    ASSERT_TRUE(success);

    success = splines.ComputeSplineProperties(false);
    ASSERT_TRUE(success);

    const double tolerance = 1e-2;
    ASSERT_NEAR(345.967070088532, splines.m_maximumGridHeights[0], tolerance);
    ASSERT_NEAR(370.339417298715, splines.m_maximumGridHeights[1], tolerance);
    ASSERT_NEAR(0.0, splines.m_maximumGridHeights[2], tolerance);
    ASSERT_NEAR(0.0, splines.m_maximumGridHeights[3], tolerance);
}

TEST(SplineTests, OrthogonalCurvilinearMeshFromSplines)
{
    std::vector<GridGeom::Point> firstSpline;
    firstSpline.push_back(GridGeom::Point{ 152.001571655273, 86.6264953613281 });
    firstSpline.push_back(GridGeom::Point{ 374.752960205078, 336.378997802734 });
    firstSpline.push_back(GridGeom::Point{ 850.255920410156, 499.130676269531 });

    GridGeom::Polygons polygon;
    GridGeom::Splines splines(GridGeom::Projections::cartesian, polygon);
    bool success = splines.AddSpline(firstSpline);
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> secondSpline;
    secondSpline.push_back(GridGeom::Point{ 72.5010681152344,391.129577636719 });
    secondSpline.push_back(GridGeom::Point{ 462.503479003906, 90.3765411376953 });
    success = splines.AddSpline(secondSpline);

    success = splines.OrthogonalCurvilinearMeshFromSplines();
    ASSERT_TRUE(success);

    ASSERT_EQ(0, splines.m_numCrossSplineLeftHeights[0][0]);
    ASSERT_EQ(1, splines.m_numCrossSplineLeftHeights[0][1]);
    ASSERT_EQ(0, splines.m_numCrossSplineLeftHeights[0][2]);
    ASSERT_EQ(GridGeom::intMissingValue, splines.m_numCrossSplineLeftHeights[0][3]);

    for (int s = 1; s < splines.m_numSplines; ++s)
    {
        for (int ss = 0; ss < splines.m_numSplines; ++ss)
        {
            ASSERT_EQ(GridGeom::intMissingValue, splines.m_numCrossSplineLeftHeights[s][ss]);
        }
    }

    const double tolerance = 1e-3;
    ASSERT_NEAR(238.968273342307, splines.m_gridHeights[0][0], tolerance);
    ASSERT_NEAR(238.968273342307, splines.m_gridHeights[0][1], tolerance);
    ASSERT_NEAR(GridGeom::doubleMissingValue, splines.m_maximumGridHeights[2], tolerance);
    ASSERT_NEAR(GridGeom::doubleMissingValue, splines.m_maximumGridHeights[3], tolerance);
    ASSERT_NEAR(253.529715955476, splines.m_maximumGridHeights[4], tolerance);
    ASSERT_NEAR(253.529715955476, splines.m_maximumGridHeights[5], tolerance);
}