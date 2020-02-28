#include "../Mesh.hpp"
#include "../Entities.hpp"
#include "../Polygons.hpp"
#include "../Constants.cpp"
#include "../Splines.hpp"
#include <gtest/gtest.h>

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
    bool success = splines.AddSpline(splineNodes, 0, splineNodes.size());
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

    bool success = splines.AddSpline(firstSpline, 0, firstSpline.size());
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> secondSpline;
    secondSpline.push_back(GridGeom::Point{ 72.5010681152344,391.129577636719 });
    secondSpline.push_back(GridGeom::Point{ 462.503479003906, 90.3765411376953 });
    success = splines.AddSpline(secondSpline, 0, secondSpline.size());

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
    bool success = splines.AddSpline(firstSpline,0, firstSpline.size());
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> secondSpline;
    secondSpline.push_back(GridGeom::Point{ 72.5010681152344,391.129577636719 });
    secondSpline.push_back(GridGeom::Point{ 462.503479003906, 90.3765411376953 });
    success = splines.AddSpline(secondSpline,0, secondSpline.size());

    success = splines.ComputeSplineProperties(false);
    ASSERT_TRUE(success);
    success = splines.MakeAllGridLines(true);
    ASSERT_TRUE(success);
    ASSERT_EQ(7, splines.m_numM);

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
    bool success = splines.AddSpline(firstSpline,0, firstSpline.size());
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> secondSpline;
    secondSpline.push_back(GridGeom::Point{ 273.502319335938, 86.6264953613281 });
    secondSpline.push_back(GridGeom::Point{ 557.004089355469, 316.128814697266 });
    secondSpline.push_back(GridGeom::Point{ 847.255920410156, 409.129730224609 });
    success = splines.AddSpline(secondSpline,0, secondSpline.size());
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> thirdSpline;
    thirdSpline.push_back(GridGeom::Point{ 62.7510070800781, 396.379608154297 });
    thirdSpline.push_back(GridGeom::Point{ 350.752807617188, 73.8763732910156 });
    success = splines.AddSpline(thirdSpline, 0, thirdSpline.size());
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> fourthSpline;
    fourthSpline.push_back(GridGeom::Point{ 704.755004882812, 636.382019042969 });
    fourthSpline.push_back(GridGeom::Point{ 845.005859375000, 285.378509521484 });
    success = splines.AddSpline(fourthSpline,0, fourthSpline.size());
    ASSERT_TRUE(success);

    success = splines.ComputeSplineProperties(false);
    ASSERT_TRUE(success);

    const double tolerance = 1e-2;
    ASSERT_NEAR(345.967070088532, splines.m_maximumGridHeights[0], tolerance);
    ASSERT_NEAR(370.339417298715, splines.m_maximumGridHeights[1], tolerance);
    ASSERT_NEAR(0.0, splines.m_maximumGridHeights[2], tolerance);
    ASSERT_NEAR(0.0, splines.m_maximumGridHeights[3], tolerance);
}

TEST(SplineTests, OrthogonalCurvilinearMeshTwoCrossingCurvatureAdapted)
{
    std::vector<GridGeom::Point> firstSpline;
    firstSpline.push_back(GridGeom::Point{ 26.9847544156713,  327.556992634968 });
    firstSpline.push_back(GridGeom::Point{ 340.139086327869, 819.656657068422 });
    firstSpline.push_back(GridGeom::Point{ 2048.50780774173, 1644.48279915859 });

    GridGeom::Polygons polygon;
    GridGeom::Splines splines(GridGeom::Projections::cartesian, polygon);
    bool success = splines.AddSpline(firstSpline, 0, firstSpline.size());
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> secondSpline;
    secondSpline.push_back(GridGeom::Point{ -179.920786312031,1068.50251010579 });
    secondSpline.push_back(GridGeom::Point{ 600.169022647819, 321.964950993679 });
    success = splines.AddSpline(secondSpline, 0, secondSpline.size());
    ASSERT_TRUE(success);

    GridGeom::CurvilinearGrid curvilinearGrid;
    success = splines.OrthogonalCurvilinearGridFromSplines(curvilinearGrid);
    ASSERT_TRUE(success);

    GridGeom::Mesh mesh(curvilinearGrid, GridGeom::Projections::cartesian);

    const double tolerance = 1e-5;
    ASSERT_NEAR(588.142743143705, mesh.m_nodes[0].x, tolerance);
    ASSERT_NEAR(469.414124846559, mesh.m_nodes[1].x, tolerance);
    ASSERT_NEAR(366.687368582373, mesh.m_nodes[2].x, tolerance);
    ASSERT_NEAR(277.805794755322, mesh.m_nodes[3].x, tolerance);
    ASSERT_NEAR(200.903394122161, mesh.m_nodes[4].x, tolerance);
    ASSERT_NEAR(134.365652175209, mesh.m_nodes[5].x, tolerance);
    ASSERT_NEAR(76.7956534892387, mesh.m_nodes[6].x, tolerance);
    ASSERT_NEAR(26.9847544156713, mesh.m_nodes[7].x, tolerance);
    ASSERT_NEAR(-22.8261446578961, mesh.m_nodes[8].x, tolerance);
    ASSERT_NEAR(-78.9948453544669, mesh.m_nodes[9].x, tolerance);
    ASSERT_NEAR(-142.332849618982, mesh.m_nodes[10].x, tolerance);
    ASSERT_NEAR(-213.755238502133, mesh.m_nodes[11].x, tolerance);
    ASSERT_NEAR(-294.293892869459, mesh.m_nodes[12].x, tolerance);
    ASSERT_NEAR(-385.112401585425, mesh.m_nodes[13].x, tolerance);
    ASSERT_NEAR(-487.522872559699, mesh.m_nodes[14].x, tolerance);
    ASSERT_NEAR(589.751674409702, mesh.m_nodes[15].x, tolerance);
    ASSERT_NEAR(475.485400123041, mesh.m_nodes[16].x, tolerance);
    ASSERT_NEAR(375.293682537274, mesh.m_nodes[17].x, tolerance);
    ASSERT_NEAR(288.605483013514, mesh.m_nodes[18].x, tolerance);
    ASSERT_NEAR(213.600840736189, mesh.m_nodes[19].x, tolerance);

    ASSERT_NEAR(278.613317010633, mesh.m_nodes[0].y, tolerance);
    ASSERT_NEAR(288.968716199430, mesh.m_nodes[1].y, tolerance);
    ASSERT_NEAR(297.928447928333, mesh.m_nodes[2].y, tolerance);
    ASSERT_NEAR(305.680615777235, mesh.m_nodes[3].y, tolerance);
    ASSERT_NEAR(312.387971329374, mesh.m_nodes[4].y, tolerance);
    ASSERT_NEAR(318.191331032065, mesh.m_nodes[5].y, tolerance);
    ASSERT_NEAR(323.212532543900, mesh.m_nodes[6].y, tolerance);
    ASSERT_NEAR(327.556992634968, mesh.m_nodes[7].y, tolerance);
    ASSERT_NEAR(331.901452726036, mesh.m_nodes[8].y, tolerance);
    ASSERT_NEAR(336.800434336126, mesh.m_nodes[9].y, tolerance);
    ASSERT_NEAR(342.324715906886, mesh.m_nodes[10].y, tolerance);
    ASSERT_NEAR(348.554109952822, mesh.m_nodes[11].y, tolerance);
    ASSERT_NEAR(355.578616159194, mesh.m_nodes[12].y, tolerance);
    ASSERT_NEAR(363.499721659902, mesh.m_nodes[13].y, tolerance);
    ASSERT_NEAR(372.431867281248, mesh.m_nodes[14].y, tolerance);
    ASSERT_NEAR(297.060330270498, mesh.m_nodes[15].y, tolerance);
    ASSERT_NEAR(358.578212823972, mesh.m_nodes[16].y, tolerance);
    ASSERT_NEAR(396.603134038195, mesh.m_nodes[17].y, tolerance);
    ASSERT_NEAR(429.503178438293, mesh.m_nodes[18].y, tolerance);
    ASSERT_NEAR(457.969060469458, mesh.m_nodes[19].y, tolerance);
}


TEST(SplineTests, OrthogonalCurvilinearMeshTwoCrossingCurvatureNotAdapted)
{
    std::vector<GridGeom::Point> firstSpline;
    firstSpline.push_back(GridGeom::Point{ 26.9847544156713,  327.556992634968 });
    firstSpline.push_back(GridGeom::Point{ 340.139086327869, 819.656657068422 });
    firstSpline.push_back(GridGeom::Point{ 2048.50780774173, 1644.48279915859 });

    GridGeom::Polygons polygon;
    GridGeom::Splines splines(GridGeom::Projections::cartesian, polygon);
    bool success = splines.AddSpline(firstSpline, 0, firstSpline.size());
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> secondSpline;
    secondSpline.push_back(GridGeom::Point{ -179.920786312031,1068.50251010579 });
    secondSpline.push_back(GridGeom::Point{ 600.169022647819, 321.964950993679 });
    success = splines.AddSpline(secondSpline, 0, secondSpline.size());
    ASSERT_TRUE(success);


    GridGeomApi::CurvilinearParametersNative curvilinearParametersNative;
    GridGeomApi::SplinesToCurvilinearParametersNative splinesToCurvilinearParametersNative;

    splinesToCurvilinearParametersNative.AspectRatio = 0.1;
    splinesToCurvilinearParametersNative.AspectRatioGrowFactor = 1.1;
    splinesToCurvilinearParametersNative.AverageWidth = 500.0;
    splinesToCurvilinearParametersNative.GridsOnTopOfEachOtherTolerance = 1e-4;
    splinesToCurvilinearParametersNative.MinimumCosineOfCrossingAngles = 0.95;
    splinesToCurvilinearParametersNative.CheckFrontCollisions = false;
    splinesToCurvilinearParametersNative.CurvatureAdapetedGridSpacing = false;
    curvilinearParametersNative.MRefinement = 20;
    curvilinearParametersNative.NRefinement = 40;

    splines.SetParameters(curvilinearParametersNative, splinesToCurvilinearParametersNative);
    GridGeom::CurvilinearGrid curvilinearGrid;
    success = splines.OrthogonalCurvilinearGridFromSplines(curvilinearGrid);
    ASSERT_TRUE(success);

    GridGeom::Mesh mesh(curvilinearGrid, GridGeom::Projections::cartesian);

    const double tolerance = 1e-5;
    ASSERT_NEAR(548.641052099198, mesh.m_nodes[0].x, tolerance);
    ASSERT_NEAR(438.270117898728, mesh.m_nodes[1].x, tolerance);
    ASSERT_NEAR(342.774623898401, mesh.m_nodes[2].x, tolerance);
    ASSERT_NEAR(260.149706020995, mesh.m_nodes[3].x, tolerance);
    ASSERT_NEAR(188.660709333624, mesh.m_nodes[4].x, tolerance);
    ASSERT_NEAR(126.806770126691, mesh.m_nodes[5].x, tolerance);
    ASSERT_NEAR(73.2893062833241, mesh.m_nodes[6].x, tolerance);
    ASSERT_NEAR(26.9847544156713, mesh.m_nodes[7].x, tolerance);
    ASSERT_NEAR(-19.3197974519815, mesh.m_nodes[8].x, tolerance);
    ASSERT_NEAR(-71.5346051169009, mesh.m_nodes[9].x, tolerance);
    ASSERT_NEAR(-130.414046327714, mesh.m_nodes[10].x, tolerance);
    ASSERT_NEAR(-196.808786676951, mesh.m_nodes[11].x, tolerance);
    ASSERT_NEAR(-271.678069662485, mesh.m_nodes[12].x, tolerance);
    ASSERT_NEAR(-356.103575437243, mesh.m_nodes[13].x, tolerance);
    ASSERT_NEAR(-451.305048472597, mesh.m_nodes[14].x, tolerance);
    ASSERT_NEAR(634.363749815992, mesh.m_nodes[15].x, tolerance);
    ASSERT_NEAR(537.793931651381, mesh.m_nodes[16].x, tolerance);
    ASSERT_NEAR(454.239483571770, mesh.m_nodes[17].x, tolerance);
    ASSERT_NEAR(381.946235778882, mesh.m_nodes[18].x, tolerance);
    ASSERT_NEAR(319.396439848637, mesh.m_nodes[19].x, tolerance);


    ASSERT_NEAR(115.028220973275, mesh.m_nodes[0].y, tolerance);
    ASSERT_NEAR(159.994608184438, mesh.m_nodes[1].y, tolerance);
    ASSERT_NEAR(198.900570077816, mesh.m_nodes[2].y, tolerance);
    ASSERT_NEAR(232.562911322639, mesh.m_nodes[3].y, tolerance);
    ASSERT_NEAR(261.688350276912, mesh.m_nodes[4].y, tolerance);
    ASSERT_NEAR(286.888356067050, mesh.m_nodes[5].y, tolerance);
    ASSERT_NEAR(308.691985973537, mesh.m_nodes[6].y, tolerance);
    ASSERT_NEAR(327.556992634968, mesh.m_nodes[7].y, tolerance);
    ASSERT_NEAR(346.421999296400, mesh.m_nodes[8].y, tolerance);
    ASSERT_NEAR(367.694912102887, mesh.m_nodes[9].y, tolerance);
    ASSERT_NEAR(391.683073220138, mesh.m_nodes[10].y, tolerance);
    ASSERT_NEAR(418.733053588210, mesh.m_nodes[11].y, tolerance);
    ASSERT_NEAR(449.235660033806, mesh.m_nodes[12].y, tolerance);
    ASSERT_NEAR(483.631581484201, mesh.m_nodes[13].y, tolerance);
    ASSERT_NEAR(522.417755856934, mesh.m_nodes[14].y, tolerance);
    ASSERT_NEAR(325.436368391402, mesh.m_nodes[15].y, tolerance);
    ASSERT_NEAR(404.277882433049, mesh.m_nodes[16].y, tolerance);
    ASSERT_NEAR(472.493390308208, mesh.m_nodes[17].y, tolerance);
    ASSERT_NEAR(531.515031016526, mesh.m_nodes[18].y, tolerance);
    ASSERT_NEAR(582.581924460712, mesh.m_nodes[19].y, tolerance);
}

TEST(SplineTests, OrthogonalCurvilinearMeshTwoCrossingCurvatureAdaptedLargeMRefinement)
{
    std::vector<GridGeom::Point> firstSpline;
    firstSpline.push_back(GridGeom::Point{ 26.9847544156713,  327.556992634968 });
    firstSpline.push_back(GridGeom::Point{ 340.139086327869, 819.656657068422 });
    firstSpline.push_back(GridGeom::Point{ 2048.50780774173, 1644.48279915859 });

    GridGeom::Polygons polygon;
    GridGeom::Splines splines(GridGeom::Projections::cartesian, polygon);
    bool success = splines.AddSpline(firstSpline, 0, firstSpline.size());
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> secondSpline;
    secondSpline.push_back(GridGeom::Point{ -179.920786312031,1068.50251010579 });
    secondSpline.push_back(GridGeom::Point{ 600.169022647819, 321.964950993679 });
    success = splines.AddSpline(secondSpline, 0, secondSpline.size());
    ASSERT_TRUE(success);

    GridGeomApi::CurvilinearParametersNative curvilinearParametersNative;
    GridGeomApi::SplinesToCurvilinearParametersNative splinesToCurvilinearParametersNative;

    splinesToCurvilinearParametersNative.AspectRatio = 0.1;
    splinesToCurvilinearParametersNative.AspectRatioGrowFactor = 1.1;
    splinesToCurvilinearParametersNative.AverageWidth = 100.0;
    splinesToCurvilinearParametersNative.GridsOnTopOfEachOtherTolerance = 1e-4;
    splinesToCurvilinearParametersNative.MinimumCosineOfCrossingAngles = 0.95;
    splinesToCurvilinearParametersNative.CheckFrontCollisions = false;
    splinesToCurvilinearParametersNative.CurvatureAdapetedGridSpacing = true;
    curvilinearParametersNative.MRefinement = 10;
    curvilinearParametersNative.NRefinement = 20;

    splines.SetParameters(curvilinearParametersNative, splinesToCurvilinearParametersNative);

    GridGeom::CurvilinearGrid curvilinearGrid;
    success = splines.OrthogonalCurvilinearGridFromSplines(curvilinearGrid);
    ASSERT_TRUE(success);

    GridGeom::Mesh mesh(curvilinearGrid, GridGeom::Projections::cartesian);

    const double tolerance = 1e-5;
}