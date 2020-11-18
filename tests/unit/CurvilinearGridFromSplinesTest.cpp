#include <MeshKernel/Mesh.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/Constants.hpp>
#include <MeshKernel/CurvilinearGrid.hpp>
#include <MeshKernel/Splines.hpp>
#include <MeshKernel/CurvilinearGridFromSplines.hpp>
#include <gtest/gtest.h>

TEST(CurvilinearGridFromSplines, ComputeSplinesProperties)
{
    std::vector<meshkernel::Point> firstSpline{{152.001571655273, 86.6264953613281},
                                               {374.752960205078, 336.378997802734},
                                               {850.255920410156, 499.130676269531}};

    auto splines = std::make_shared<meshkernel::Splines>(meshkernel::Projections::cartesian);
    splines->AddSpline(firstSpline, 0, firstSpline.size());

    std::vector<meshkernel::Point> secondSpline{{72.5010681152344, 391.129577636719},
                                                {462.503479003906, 90.3765411376953}};

    splines->AddSpline(secondSpline, 0, secondSpline.size());

    // construct CurvilinearGridFromSplines
    meshkernelapi::CurvilinearParametersNative curvilinearParametersNative;
    meshkernelapi::SplinesToCurvilinearParametersNative splinesToCurvilinearParametersNative;

    splinesToCurvilinearParametersNative.AspectRatio = 0.1;
    splinesToCurvilinearParametersNative.AspectRatioGrowFactor = 1.1;
    splinesToCurvilinearParametersNative.AverageWidth = 500.0;
    splinesToCurvilinearParametersNative.GridsOnTopOfEachOtherTolerance = 1e-4;
    splinesToCurvilinearParametersNative.MinimumCosineOfCrossingAngles = 0.95;
    splinesToCurvilinearParametersNative.CheckFrontCollisions = false;
    splinesToCurvilinearParametersNative.CurvatureAdapetedGridSpacing = true;
    splinesToCurvilinearParametersNative.RemoveSkinnyTriangles = 0;
    curvilinearParametersNative.MRefinement = 20;
    curvilinearParametersNative.NRefinement = 40;

    meshkernel::CurvilinearGridFromSplines curvilinearGridFromSplines(splines, curvilinearParametersNative, splinesToCurvilinearParametersNative);

    curvilinearGridFromSplines.ComputeSplineProperties(false);
    curvilinearGridFromSplines.MakeAllGridLines();
    ASSERT_EQ(7, curvilinearGridFromSplines.m_numM);

    const double tolerance = 1e-6;
    ASSERT_NEAR(253.52971595547601, curvilinearGridFromSplines.m_maximumGridHeights[0], tolerance);
    ASSERT_NEAR(0.0, curvilinearGridFromSplines.m_maximumGridHeights[1], tolerance);

    ASSERT_NEAR(152.001571655273, curvilinearGridFromSplines.m_gridLine[0].x, tolerance);
    ASSERT_NEAR(407.924702423872, curvilinearGridFromSplines.m_gridLine[1].x, tolerance);
    ASSERT_NEAR(850.250753009544, curvilinearGridFromSplines.m_gridLine[2].x, tolerance);
    ASSERT_NEAR(-999.000000000000, curvilinearGridFromSplines.m_gridLine[3].x, tolerance);
    ASSERT_NEAR(850.250753009544, curvilinearGridFromSplines.m_gridLine[4].x, tolerance);
    ASSERT_NEAR(407.924702423872, curvilinearGridFromSplines.m_gridLine[5].x, tolerance);
    ASSERT_NEAR(152.001571655273, curvilinearGridFromSplines.m_gridLine[6].x, tolerance);

    ASSERT_NEAR(86.6264953613281, curvilinearGridFromSplines.m_gridLine[0].y, tolerance);
    ASSERT_NEAR(354.562246561336, curvilinearGridFromSplines.m_gridLine[1].y, tolerance);
    ASSERT_NEAR(499.129323710654, curvilinearGridFromSplines.m_gridLine[2].y, tolerance);
    ASSERT_NEAR(-999.000000000000, curvilinearGridFromSplines.m_gridLine[3].y, tolerance);
    ASSERT_NEAR(499.129323710654, curvilinearGridFromSplines.m_gridLine[4].y, tolerance);
    ASSERT_NEAR(354.562246561336, curvilinearGridFromSplines.m_gridLine[5].y, tolerance);
    ASSERT_NEAR(86.6264953613281, curvilinearGridFromSplines.m_gridLine[6].y, tolerance);

    ASSERT_NEAR(4.214936859441928e-14, curvilinearGridFromSplines.m_gridLineDimensionalCoordinates[0], tolerance);
    ASSERT_NEAR(1.09068327269294, curvilinearGridFromSplines.m_gridLineDimensionalCoordinates[1], tolerance);
    ASSERT_NEAR(1.99999040748403, curvilinearGridFromSplines.m_gridLineDimensionalCoordinates[2], tolerance);
    ASSERT_NEAR(-999.000000000000, curvilinearGridFromSplines.m_gridLineDimensionalCoordinates[3], tolerance);
    ASSERT_NEAR(1.99999040748403, curvilinearGridFromSplines.m_gridLineDimensionalCoordinates[4], tolerance);
    ASSERT_NEAR(1.09068327269294, curvilinearGridFromSplines.m_gridLineDimensionalCoordinates[5], tolerance);
    ASSERT_NEAR(4.214936859441928e-14, curvilinearGridFromSplines.m_gridLineDimensionalCoordinates[6], tolerance);
}

TEST(CurvilinearGridFromSplines, ComputeBoundingBox)
{
    std::vector<meshkernel::Point> firstSpline{{91.2511901855469, 299.628631591797},
                                               {354.502838134766, 518.630859375000},
                                               {770.755432128906, 607.881774902344}};

    auto splines = std::make_shared<meshkernel::Splines>(meshkernel::Projections::cartesian);
    splines->AddSpline(firstSpline, 0, firstSpline.size());

    std::vector<meshkernel::Point> secondSpline;
    secondSpline.emplace_back(meshkernel::Point{273.502319335938, 86.6264953613281});
    secondSpline.emplace_back(meshkernel::Point{557.004089355469, 316.128814697266});
    secondSpline.emplace_back(meshkernel::Point{847.255920410156, 409.129730224609});
    splines->AddSpline(secondSpline, 0, secondSpline.size());

    std::vector<meshkernel::Point> thirdSpline;
    thirdSpline.emplace_back(meshkernel::Point{62.7510070800781, 396.379608154297});
    thirdSpline.emplace_back(meshkernel::Point{350.752807617188, 73.8763732910156});
    splines->AddSpline(thirdSpline, 0, thirdSpline.size());

    std::vector<meshkernel::Point> fourthSpline;
    fourthSpline.emplace_back(meshkernel::Point{704.755004882812, 636.382019042969});
    fourthSpline.emplace_back(meshkernel::Point{845.005859375000, 285.378509521484});
    splines->AddSpline(fourthSpline, 0, fourthSpline.size());

    // construct CurvilinearGridFromSplines
    meshkernelapi::CurvilinearParametersNative curvilinearParametersNative;
    meshkernelapi::SplinesToCurvilinearParametersNative splinesToCurvilinearParametersNative;

    splinesToCurvilinearParametersNative.AspectRatio = 0.1;
    splinesToCurvilinearParametersNative.AspectRatioGrowFactor = 1.1;
    splinesToCurvilinearParametersNative.AverageWidth = 500.0;
    splinesToCurvilinearParametersNative.GridsOnTopOfEachOtherTolerance = 1e-4;
    splinesToCurvilinearParametersNative.MinimumCosineOfCrossingAngles = 0.95;
    splinesToCurvilinearParametersNative.CheckFrontCollisions = false;
    splinesToCurvilinearParametersNative.CurvatureAdapetedGridSpacing = false;
    splinesToCurvilinearParametersNative.RemoveSkinnyTriangles = 0;
    curvilinearParametersNative.MRefinement = 20;
    curvilinearParametersNative.NRefinement = 40;

    meshkernel::CurvilinearGridFromSplines curvilinearGridFromSplines(splines, curvilinearParametersNative, splinesToCurvilinearParametersNative);
    curvilinearGridFromSplines.ComputeSplineProperties(false);

    const double tolerance = 1e-6;
    ASSERT_NEAR(345.967070088532, curvilinearGridFromSplines.m_maximumGridHeights[0], tolerance);
    ASSERT_NEAR(370.339417298715, curvilinearGridFromSplines.m_maximumGridHeights[1], tolerance);
    ASSERT_NEAR(0.0, curvilinearGridFromSplines.m_maximumGridHeights[2], tolerance);
    ASSERT_NEAR(0.0, curvilinearGridFromSplines.m_maximumGridHeights[3], tolerance);
}

TEST(CurvilinearGridFromSplines, OrthogonalCurvilinearMeshTwoCrossingCurvatureAdapted)
{
    std::vector<meshkernel::Point> firstSpline{{26.9847544156713, 327.556992634968},
                                               {340.139086327869, 819.656657068422},
                                               {2048.50780774173, 1644.48279915859}};

    auto splines = std::make_shared<meshkernel::Splines>(meshkernel::Projections::cartesian);
    splines->AddSpline(firstSpline, 0, firstSpline.size());

    std::vector<meshkernel::Point> secondSpline;
    secondSpline.emplace_back(meshkernel::Point{-179.920786312031, 1068.50251010579});
    secondSpline.emplace_back(meshkernel::Point{600.169022647819, 321.964950993679});
    splines->AddSpline(secondSpline, 0, secondSpline.size());

    // construct CurvilinearGridFromSplines
    meshkernelapi::CurvilinearParametersNative curvilinearParametersNative;
    meshkernelapi::SplinesToCurvilinearParametersNative splinesToCurvilinearParametersNative;

    splinesToCurvilinearParametersNative.AspectRatio = 0.1;
    splinesToCurvilinearParametersNative.AspectRatioGrowFactor = 1.1;
    splinesToCurvilinearParametersNative.AverageWidth = 500.0;
    splinesToCurvilinearParametersNative.GridsOnTopOfEachOtherTolerance = 1e-4;
    splinesToCurvilinearParametersNative.MinimumCosineOfCrossingAngles = 0.95;
    splinesToCurvilinearParametersNative.CheckFrontCollisions = false;
    splinesToCurvilinearParametersNative.CurvatureAdapetedGridSpacing = true;
    splinesToCurvilinearParametersNative.RemoveSkinnyTriangles = 0;
    curvilinearParametersNative.MRefinement = 20;
    curvilinearParametersNative.NRefinement = 40;

    meshkernel::CurvilinearGridFromSplines curvilinearGridFromSplines(splines, curvilinearParametersNative, splinesToCurvilinearParametersNative);
    meshkernel::CurvilinearGrid curvilinearGrid;
    curvilinearGridFromSplines.Compute(curvilinearGrid);

    meshkernel::Mesh mesh(curvilinearGrid, meshkernel::Projections::cartesian);

    const double tolerance = 1e-6;
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

TEST(CurvilinearGridFromSplines, OrthogonalCurvilinearMeshTwoCrossingCurvatureNotAdapted)
{
    std::vector<meshkernel::Point> firstSpline{{26.9847544156713, 327.556992634968},
                                               {340.139086327869, 819.656657068422},
                                               {2048.50780774173, 1644.48279915859}};

    auto splines = std::make_shared<meshkernel::Splines>(meshkernel::Projections::cartesian);
    splines->AddSpline(firstSpline, 0, firstSpline.size());

    std::vector<meshkernel::Point> secondSpline;
    secondSpline.emplace_back(meshkernel::Point{-179.920786312031, 1068.50251010579});
    secondSpline.emplace_back(meshkernel::Point{600.169022647819, 321.964950993679});
    splines->AddSpline(secondSpline, 0, secondSpline.size());

    meshkernelapi::CurvilinearParametersNative curvilinearParametersNative;
    meshkernelapi::SplinesToCurvilinearParametersNative splinesToCurvilinearParametersNative;

    splinesToCurvilinearParametersNative.AspectRatio = 0.1;
    splinesToCurvilinearParametersNative.AspectRatioGrowFactor = 1.1;
    splinesToCurvilinearParametersNative.AverageWidth = 500.0;
    splinesToCurvilinearParametersNative.GridsOnTopOfEachOtherTolerance = 1e-4;
    splinesToCurvilinearParametersNative.MinimumCosineOfCrossingAngles = 0.95;
    splinesToCurvilinearParametersNative.CheckFrontCollisions = false;
    splinesToCurvilinearParametersNative.CurvatureAdapetedGridSpacing = false;
    splinesToCurvilinearParametersNative.RemoveSkinnyTriangles = 0;
    curvilinearParametersNative.MRefinement = 20;
    curvilinearParametersNative.NRefinement = 40;

    meshkernel::CurvilinearGridFromSplines curvilinearGridFromSplines(splines, curvilinearParametersNative, splinesToCurvilinearParametersNative);

    meshkernel::CurvilinearGrid curvilinearGrid;
    curvilinearGridFromSplines.Compute(curvilinearGrid);

    meshkernel::Mesh mesh(curvilinearGrid, meshkernel::Projections::cartesian);

    const double tolerance = 1e-6;
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

TEST(CurvilinearGridFromSplines, OrthogonalCurvilinearMeshTwoCrossingCurvatureAdaptedLargeMRefinement)
{
    std::vector<meshkernel::Point> firstSpline{{26.9847544156713, 327.556992634968},
                                               {340.139086327869, 819.656657068422},
                                               {2048.50780774173, 1644.48279915859}};

    auto splines = std::make_shared<meshkernel::Splines>(meshkernel::Projections::cartesian);
    splines->AddSpline(firstSpline, 0, firstSpline.size());

    std::vector<meshkernel::Point> secondSpline;
    secondSpline.emplace_back(meshkernel::Point{-179.920786312031, 1068.50251010579});
    secondSpline.emplace_back(meshkernel::Point{600.169022647819, 321.964950993679});
    splines->AddSpline(secondSpline, 0, secondSpline.size());

    meshkernelapi::CurvilinearParametersNative curvilinearParametersNative;
    meshkernelapi::SplinesToCurvilinearParametersNative splinesToCurvilinearParametersNative;

    splinesToCurvilinearParametersNative.AspectRatio = 0.1;
    splinesToCurvilinearParametersNative.AspectRatioGrowFactor = 1.1;
    splinesToCurvilinearParametersNative.AverageWidth = 500.0;
    splinesToCurvilinearParametersNative.GridsOnTopOfEachOtherTolerance = 1e-4;
    splinesToCurvilinearParametersNative.MinimumCosineOfCrossingAngles = 0.95;
    splinesToCurvilinearParametersNative.CheckFrontCollisions = false;
    splinesToCurvilinearParametersNative.CurvatureAdapetedGridSpacing = true;
    splinesToCurvilinearParametersNative.RemoveSkinnyTriangles = 0;
    curvilinearParametersNative.MRefinement = 10;
    curvilinearParametersNative.NRefinement = 20;

    // create the algorithm
    meshkernel::CurvilinearGridFromSplines curvilinearGridFromSplines(splines, curvilinearParametersNative, splinesToCurvilinearParametersNative);

    // compute
    meshkernel::CurvilinearGrid curvilinearGrid;
    curvilinearGridFromSplines.Compute(curvilinearGrid);

    meshkernel::Mesh mesh(curvilinearGrid, meshkernel::Projections::cartesian);

    const double tolerance = 1e-6;

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
    ASSERT_NEAR(148.705083366035, mesh.m_nodes[20].x, tolerance);
    ASSERT_NEAR(92.5557678466786, mesh.m_nodes[21].x, tolerance);
    ASSERT_NEAR(43.9740768231890, mesh.m_nodes[22].x, tolerance);
    ASSERT_NEAR(-4.60761420030062, mesh.m_nodes[23].x, tolerance);
    ASSERT_NEAR(-59.3902122498050, mesh.m_nodes[24].x, tolerance);
    ASSERT_NEAR(-121.165193437728, mesh.m_nodes[25].x, tolerance);
    ASSERT_NEAR(-190.825056909698, mesh.m_nodes[26].x, tolerance);
    ASSERT_NEAR(-269.376219299725, mesh.m_nodes[27].x, tolerance);
    ASSERT_NEAR(-357.953555017692, mesh.m_nodes[28].x, tolerance);
    ASSERT_NEAR(-457.836792441174, mesh.m_nodes[29].x, tolerance);
    ASSERT_NEAR(589.751674409702, mesh.m_nodes[30].x, tolerance);
    ASSERT_NEAR(481.252352883142, mesh.m_nodes[31].x, tolerance);
    ASSERT_NEAR(408.289995573200, mesh.m_nodes[32].x, tolerance);
    ASSERT_NEAR(345.161270558643, mesh.m_nodes[33].x, tolerance);
    ASSERT_NEAR(290.540832446502, mesh.m_nodes[34].x, tolerance);
    ASSERT_NEAR(243.281961641267, mesh.m_nodes[35].x, tolerance);
    ASSERT_NEAR(202.392489733323, mesh.m_nodes[36].x, tolerance);
    ASSERT_NEAR(167.013969586259, mesh.m_nodes[37].x, tolerance);
    ASSERT_NEAR(131.635449439195, mesh.m_nodes[38].x, tolerance);
    ASSERT_NEAR(91.7412586860200, mesh.m_nodes[39].x, tolerance);

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
    ASSERT_NEAR(482.598402301283, mesh.m_nodes[20].y, tolerance);
    ASSERT_NEAR(503.908280505820, mesh.m_nodes[21].y, tolerance);
    ASSERT_NEAR(522.346081734654, mesh.m_nodes[22].y, tolerance);
    ASSERT_NEAR(540.783882963489, mesh.m_nodes[23].y, tolerance);
    ASSERT_NEAR(561.575062363719, mesh.m_nodes[24].y, tolerance);
    ASSERT_NEAR(585.020002217855, mesh.m_nodes[25].y, tolerance);
    ASSERT_NEAR(611.457425231865, mesh.m_nodes[26].y, tolerance);
    ASSERT_NEAR(641.269288259460, mesh.m_nodes[27].y, tolerance);
    ASSERT_NEAR(674.886300655303, mesh.m_nodes[28].y, tolerance);
    ASSERT_NEAR(712.794146983999, mesh.m_nodes[29].y, tolerance);
    ASSERT_NEAR(297.060330270498, mesh.m_nodes[30].y, tolerance);
    ASSERT_NEAR(366.348803073513, mesh.m_nodes[31].y, tolerance);
    ASSERT_NEAR(441.063499679079, mesh.m_nodes[32].y, tolerance);
    ASSERT_NEAR(505.708389324335, mesh.m_nodes[33].y, tolerance);
    ASSERT_NEAR(561.640648265421, mesh.m_nodes[34].y, tolerance);
    ASSERT_NEAR(610.034536899468, mesh.m_nodes[35].y, tolerance);
    ASSERT_NEAR(651.906052576876, mesh.m_nodes[36].y, tolerance);
    ASSERT_NEAR(688.134259786700, mesh.m_nodes[37].y, tolerance);
    ASSERT_NEAR(724.362466996525, mesh.m_nodes[38].y, tolerance);
    ASSERT_NEAR(765.214797819475, mesh.m_nodes[39].y, tolerance);
}

TEST(CurvilinearGridFromSplines, OrthogonalCurvilinearMeshFourSplineCrossingFront)
{
    std::vector<meshkernel::Point> firstSpline{{-93.379753864775, 231.718492951018},
                                               {72.8139687226708, 724.468302026077},
                                               {746.335897103372, 234.634172294657},
                                               {1498.58116776234, 776.950530211586}};

    auto splines = std::make_shared<meshkernel::Splines>(meshkernel::Projections::cartesian);
    splines->AddSpline(firstSpline, 0, firstSpline.size());

    std::vector<meshkernel::Point> secondSpline{{-250.826438421303, 394.996536194825},
                                                {101.970762159065, 292.947759167446}};

    splines->AddSpline(secondSpline, 0, secondSpline.size());

    std::vector<meshkernel::Point> thirdSpline{{647.202799419633, 482.466916504007},
                                               {486.840435519465, 144.248112641836}};

    splines->AddSpline(thirdSpline, 0, thirdSpline.size());

    std::vector<meshkernel::Point> fourthSpline{{1224.50730946023, 747.793736775192},
                                                {1568.55747200968, 464.97284044217}};

    splines->AddSpline(fourthSpline, 0, fourthSpline.size());

    meshkernelapi::CurvilinearParametersNative curvilinearParametersNative;
    meshkernelapi::SplinesToCurvilinearParametersNative splinesToCurvilinearParametersNative;

    splinesToCurvilinearParametersNative.AspectRatio = 0.1;
    splinesToCurvilinearParametersNative.AspectRatioGrowFactor = 1.1;
    splinesToCurvilinearParametersNative.AverageWidth = 500.0;
    splinesToCurvilinearParametersNative.GridsOnTopOfEachOtherTolerance = 1e-4;
    splinesToCurvilinearParametersNative.MinimumCosineOfCrossingAngles = 0.95;
    splinesToCurvilinearParametersNative.CheckFrontCollisions = false;
    splinesToCurvilinearParametersNative.CurvatureAdapetedGridSpacing = true;
    splinesToCurvilinearParametersNative.RemoveSkinnyTriangles = 0;
    curvilinearParametersNative.MRefinement = 10;
    curvilinearParametersNative.NRefinement = 20;

    // create the algorithm
    meshkernel::CurvilinearGridFromSplines curvilinearGridFromSplines(splines, curvilinearParametersNative, splinesToCurvilinearParametersNative);

    // compute
    meshkernel::CurvilinearGrid curvilinearGrid;
    curvilinearGridFromSplines.Compute(curvilinearGrid);

    meshkernel::Mesh mesh(curvilinearGrid, meshkernel::Projections::cartesian);

    const double tolerance = 1e-6;
    ASSERT_NEAR(100.529546838174, mesh.m_nodes[0].x, tolerance);
    ASSERT_NEAR(53.2575912534712, mesh.m_nodes[1].x, tolerance);
    ASSERT_NEAR(5.19099854865588, mesh.m_nodes[2].x, tolerance);
    ASSERT_NEAR(-43.6835890514023, mesh.m_nodes[3].x, tolerance);
    ASSERT_NEAR(-93.3797538647750, mesh.m_nodes[4].x, tolerance);
    ASSERT_NEAR(-143.075918678148, mesh.m_nodes[5].x, tolerance);
    ASSERT_NEAR(-199.566533061236, mesh.m_nodes[6].x, tolerance);
    ASSERT_NEAR(-263.780532780588, mesh.m_nodes[7].x, tolerance);
    ASSERT_NEAR(125.956639277939, mesh.m_nodes[8].x, tolerance);
    ASSERT_NEAR(87.1910799819445, mesh.m_nodes[9].x, tolerance);
    ASSERT_NEAR(46.9123124226315, mesh.m_nodes[10].x, tolerance);
    ASSERT_NEAR(-0.211223295124292, mesh.m_nodes[11].x, tolerance);
    ASSERT_NEAR(-47.9144545514771, mesh.m_nodes[12].x, tolerance);
    ASSERT_NEAR(-95.6176858078300, mesh.m_nodes[13].x, tolerance);
    ASSERT_NEAR(-149.564378093483, mesh.m_nodes[14].x, tolerance);
    ASSERT_NEAR(-210.561033694981, mesh.m_nodes[15].x, tolerance);
    ASSERT_NEAR(125.956639277939, mesh.m_nodes[16].x, tolerance);
    ASSERT_NEAR(87.1910799819445, mesh.m_nodes[17].x, tolerance);
    ASSERT_NEAR(53.1366668535683, mesh.m_nodes[18].x, tolerance);
    ASSERT_NEAR(48.8694421214994, mesh.m_nodes[19].x, tolerance);

    ASSERT_NEAR(210.243450776517, mesh.m_nodes[0].y, tolerance);
    ASSERT_NEAR(215.478719204450, mesh.m_nodes[1].y, tolerance);
    ASSERT_NEAR(220.801992001466, mesh.m_nodes[2].y, tolerance);
    ASSERT_NEAR(226.214748512720, mesh.m_nodes[3].y, tolerance);
    ASSERT_NEAR(231.718492951018, mesh.m_nodes[4].y, tolerance);
    ASSERT_NEAR(237.222237389316, mesh.m_nodes[5].y, tolerance);
    ASSERT_NEAR(243.478452650484, mesh.m_nodes[6].y, tolerance);
    ASSERT_NEAR(250.590016392746, mesh.m_nodes[7].y, tolerance);
    ASSERT_NEAR(439.837864321305, mesh.m_nodes[8].y, tolerance);
    ASSERT_NEAR(521.881799407137, mesh.m_nodes[9].y, tolerance);
    ASSERT_NEAR(597.525387092956, mesh.m_nodes[10].y, tolerance);
    ASSERT_NEAR(618.749299793464, mesh.m_nodes[11].y, tolerance);
    ASSERT_NEAR(642.248275591448, mesh.m_nodes[12].y, tolerance);
    ASSERT_NEAR(665.747251389432, mesh.m_nodes[13].y, tolerance);
    ASSERT_NEAR(694.973859488393, mesh.m_nodes[14].y, tolerance);
    ASSERT_NEAR(731.136492976853, mesh.m_nodes[15].y, tolerance);
    ASSERT_NEAR(439.837864321305, mesh.m_nodes[16].y, tolerance);
    ASSERT_NEAR(521.881799407137, mesh.m_nodes[17].y, tolerance);
    ASSERT_NEAR(603.508778997112, mesh.m_nodes[18].y, tolerance);
    ASSERT_NEAR(665.929912556253, mesh.m_nodes[19].y, tolerance);
}
TEST(CurvilinearGridFromSplines, OrthogonalCurvilinearGridFromSplineWithSevenSplies)
{
    auto splines = std::make_shared<meshkernel::Splines>(meshkernel::Projections::cartesian);

    std::vector<meshkernel::Point> firstCentralSpline{{1.542516E+02, 7.687640E+01},
                                                      {3.192526E+02, 2.733784E+02},
                                                      {6.350046E+02, 5.231309E+02}};
    splines->AddSpline(firstCentralSpline, 0, firstCentralSpline.size());

    std::vector<meshkernel::Point> secondCentralSpline{{1.310014E+02, 9.637659E+01},
                                                       {2.960025E+02, 2.891285E+02},
                                                       {6.222545E+02, 5.418811E+02}};
    splines->AddSpline(secondCentralSpline, 0, secondCentralSpline.size());

    std::vector<meshkernel::Point> thirdCentralSpline{{1.782517E+02, 5.662619E+01},
                                                      {3.335027E+02, 2.448781E+02},
                                                      {6.455046E+02, 4.983806E+02}};

    splines->AddSpline(thirdCentralSpline, 0, thirdCentralSpline.size());

    std::vector<meshkernel::Point> fourthBankSpline{{3.500084E+01, 1.901275E+02},
                                                    {2.195020E+02, 3.813795E+02},
                                                    {5.727542E+02, 6.221319E+02}};

    splines->AddSpline(fourthBankSpline, 0, fourthBankSpline.size());

    std::vector<meshkernel::Point> fifthBankSpline{{3.177526E+02, -3.712475E+01},
                                                   {4.377534E+02, 1.218768E+02},
                                                   {7.445052E+02, 4.136298E+02}};

    splines->AddSpline(fifthBankSpline, 0, fifthBankSpline.size());

    std::vector<meshkernel::Point> sixthBankSpline{{1.250633E+00, 2.748784E+02},
                                                   {3.620029E+02, -3.337471E+01}};
    splines->AddSpline(sixthBankSpline, 0, sixthBankSpline.size());

    std::vector<meshkernel::Point> seventhBankSpline{{5.030038E+02, 6.476321E+02},
                                                     {7.542553E+02, 3.378790E+02}};

    splines->AddSpline(seventhBankSpline, 0, seventhBankSpline.size());

    meshkernelapi::CurvilinearParametersNative curvilinearParametersNative;
    meshkernelapi::SplinesToCurvilinearParametersNative splinesToCurvilinearParametersNative;

    splinesToCurvilinearParametersNative.AspectRatio = 0.1;
    splinesToCurvilinearParametersNative.AspectRatioGrowFactor = 1.1;
    splinesToCurvilinearParametersNative.AverageWidth = 50.0;
    splinesToCurvilinearParametersNative.GridsOnTopOfEachOtherTolerance = 1e-4;
    splinesToCurvilinearParametersNative.MinimumCosineOfCrossingAngles = 0.95;
    splinesToCurvilinearParametersNative.CheckFrontCollisions = false;
    splinesToCurvilinearParametersNative.CurvatureAdapetedGridSpacing = true;
    splinesToCurvilinearParametersNative.RemoveSkinnyTriangles = 0;
    curvilinearParametersNative.MRefinement = 20;
    curvilinearParametersNative.NRefinement = 40;

    // create the algorithm
    meshkernel::CurvilinearGridFromSplines curvilinearGridFromSplines(splines, curvilinearParametersNative, splinesToCurvilinearParametersNative);

    // compute
    meshkernel::CurvilinearGrid curvilinearGrid;
    curvilinearGridFromSplines.Compute(curvilinearGrid);

    meshkernel::Mesh mesh(curvilinearGrid, meshkernel::Projections::cartesian);

    const double tolerance = 1e-6;

    ASSERT_NEAR(318.512311735701, mesh.m_nodes[0].x, tolerance);
    ASSERT_NEAR(299.483532213548, mesh.m_nodes[1].x, tolerance);
    ASSERT_NEAR(282.447339060415, mesh.m_nodes[2].x, tolerance);
    ASSERT_NEAR(267.192527548330, mesh.m_nodes[3].x, tolerance);
    ASSERT_NEAR(253.530637237784, mesh.m_nodes[4].x, tolerance);
    ASSERT_NEAR(241.293455383702, mesh.m_nodes[5].x, tolerance);
    ASSERT_NEAR(230.330801968344, mesh.m_nodes[6].x, tolerance);
    ASSERT_NEAR(220.508562206121, mesh.m_nodes[7].x, tolerance);
    ASSERT_NEAR(211.706937639478, mesh.m_nodes[8].x, tolerance);
    ASSERT_NEAR(203.818890903662, mesh.m_nodes[9].x, tolerance);
    ASSERT_NEAR(196.748762361429, mesh.m_nodes[10].x, tolerance);
    ASSERT_NEAR(190.411039392839, mesh.m_nodes[11].x, tolerance);
    ASSERT_NEAR(184.729261287761, mesh.m_nodes[12].x, tolerance);
    ASSERT_NEAR(179.635044611790, mesh.m_nodes[13].x, tolerance);
    ASSERT_NEAR(174.546633360283, mesh.m_nodes[14].x, tolerance);
    ASSERT_NEAR(169.464055776080, mesh.m_nodes[15].x, tolerance);
    ASSERT_NEAR(164.387339610100, mesh.m_nodes[16].x, tolerance);
    ASSERT_NEAR(159.316512108475, mesh.m_nodes[17].x, tolerance);
    ASSERT_NEAR(154.251600000000, mesh.m_nodes[18].x, tolerance);
    ASSERT_NEAR(149.317118005829, mesh.m_nodes[19].x, tolerance);

    ASSERT_NEAR(-32.9072706211971, mesh.m_nodes[0].y, tolerance);
    ASSERT_NEAR(-20.6818375059078, mesh.m_nodes[1].y, tolerance);
    ASSERT_NEAR(-9.63982764890134, mesh.m_nodes[2].y, tolerance);
    ASSERT_NEAR(0.328334959383014, mesh.m_nodes[3].y, tolerance);
    ASSERT_NEAR(9.32293085146806, mesh.m_nodes[4].y, tolerance);
    ASSERT_NEAR(17.4355807455482, mesh.m_nodes[5].y, tolerance);
    ASSERT_NEAR(24.7498708103273, mesh.m_nodes[6].y, tolerance);
    ASSERT_NEAR(31.3419568890594, mesh.m_nodes[7].y, tolerance);
    ASSERT_NEAR(37.2811417394207, mesh.m_nodes[8].y, tolerance);
    ASSERT_NEAR(42.6304221098499, mesh.m_nodes[9].y, tolerance);
    ASSERT_NEAR(47.4470041338451, mesh.m_nodes[10].y, tolerance);
    ASSERT_NEAR(51.7827865718035, mesh.m_nodes[11].y, tolerance);
    ASSERT_NEAR(55.6848120567432, mesh.m_nodes[12].y, tolerance);
    ASSERT_NEAR(59.1956869754190, mesh.m_nodes[13].y, tolerance);
    ASSERT_NEAR(62.7149705917065, mesh.m_nodes[14].y, tolerance);
    ASSERT_NEAR(66.2426739999994, mesh.m_nodes[15].y, tolerance);
    ASSERT_NEAR(69.7788073694307, mesh.m_nodes[16].y, tolerance);
    ASSERT_NEAR(73.3233799367651, mesh.m_nodes[17].y, tolerance);
    ASSERT_NEAR(76.8764000000000, mesh.m_nodes[18].y, tolerance);
    ASSERT_NEAR(80.3379237444341, mesh.m_nodes[19].y, tolerance);
}

TEST(CurvilinearGridFromSplines, OrthogonalCurvilinearMeshTwoCrossingCurvatureAdaptedSpherical)
{
    std::vector<meshkernel::Point> firstSpline{{4.109727E+01, 4.110174E+01},
                                               {4.109865E+01, 4.110418E+01},
                                               {4.110644E+01, 4.110904E+01}};

    auto splines = std::make_shared<meshkernel::Splines>(meshkernel::Projections::spherical);
    splines->AddSpline(firstSpline, 0, firstSpline.size());

    std::vector<meshkernel::Point> secondSpline{{4.109612E+01, 4.110473E+01},
                                                {4.109923E+01, 4.110212E+01}};
    splines->AddSpline(secondSpline, 0, secondSpline.size());

    meshkernelapi::CurvilinearParametersNative curvilinearParametersNative;
    meshkernelapi::SplinesToCurvilinearParametersNative splinesToCurvilinearParametersNative;

    splinesToCurvilinearParametersNative.AspectRatio = 0.1;
    splinesToCurvilinearParametersNative.AspectRatioGrowFactor = 1.1;
    splinesToCurvilinearParametersNative.AverageWidth = 100.0;
    splinesToCurvilinearParametersNative.GridsOnTopOfEachOtherTolerance = 1e-4;
    splinesToCurvilinearParametersNative.MinimumCosineOfCrossingAngles = 0.95;
    splinesToCurvilinearParametersNative.CheckFrontCollisions = false;
    splinesToCurvilinearParametersNative.CurvatureAdapetedGridSpacing = true;
    splinesToCurvilinearParametersNative.RemoveSkinnyTriangles = 0;
    curvilinearParametersNative.MRefinement = 10;
    curvilinearParametersNative.NRefinement = 20;

    // create the algorithm
    meshkernel::CurvilinearGridFromSplines curvilinearGridFromSplines(splines, curvilinearParametersNative, splinesToCurvilinearParametersNative);

    // compute
    meshkernel::CurvilinearGrid curvilinearGrid;
    curvilinearGridFromSplines.Compute(curvilinearGrid);

    meshkernel::Mesh mesh(curvilinearGrid, meshkernel::Projections::spherical);

    const double tolerance = 1e-6;

    ASSERT_NEAR(41.0994668963040, mesh.m_nodes[0].x, tolerance);
    ASSERT_NEAR(41.0991084975934, mesh.m_nodes[1].x, tolerance);
    ASSERT_NEAR(41.0987913585479, mesh.m_nodes[2].x, tolerance);
    ASSERT_NEAR(41.0985107292634, mesh.m_nodes[3].x, tolerance);
    ASSERT_NEAR(41.0982624066548, mesh.m_nodes[4].x, tolerance);
    ASSERT_NEAR(41.0980426715056, mesh.m_nodes[5].x, tolerance);
    ASSERT_NEAR(41.0978482327633, mesh.m_nodes[6].x, tolerance);
    ASSERT_NEAR(41.0976761782490, mesh.m_nodes[7].x, tolerance);
    ASSERT_NEAR(41.0975239310399, mesh.m_nodes[8].x, tolerance);
    ASSERT_NEAR(41.0973892108745, mesh.m_nodes[9].x, tolerance);
    ASSERT_NEAR(41.0972700000000, mesh.m_nodes[10].x, tolerance);
    ASSERT_NEAR(41.0971507891255, mesh.m_nodes[11].x, tolerance);
    ASSERT_NEAR(41.0970173108791, mesh.m_nodes[12].x, tolerance);
    ASSERT_NEAR(41.0968678577157, mesh.m_nodes[13].x, tolerance);
    ASSERT_NEAR(41.0967005177284, mesh.m_nodes[14].x, tolerance);
    ASSERT_NEAR(41.0965131501894, mesh.m_nodes[15].x, tolerance);
    ASSERT_NEAR(41.0963033581646, mesh.m_nodes[16].x, tolerance);
    ASSERT_NEAR(41.0960684578508, mesh.m_nodes[17].x, tolerance);
    ASSERT_NEAR(41.0958054442416, mesh.m_nodes[18].x, tolerance);
    ASSERT_NEAR(41.0955109526861, mesh.m_nodes[19].x, tolerance);

    ASSERT_NEAR(41.1017323378167, mesh.m_nodes[0].y, tolerance);
    ASSERT_NEAR(41.1017335906480, mesh.m_nodes[1].y, tolerance);
    ASSERT_NEAR(41.1017346983120, mesh.m_nodes[2].y, tolerance);
    ASSERT_NEAR(41.1017356777341, mesh.m_nodes[3].y, tolerance);
    ASSERT_NEAR(41.1017365438410, mesh.m_nodes[4].y, tolerance);
    ASSERT_NEAR(41.1017373098033, mesh.m_nodes[5].y, tolerance);
    ASSERT_NEAR(41.1017379872469, mesh.m_nodes[6].y, tolerance);
    ASSERT_NEAR(41.1017385864373, mesh.m_nodes[7].y, tolerance);
    ASSERT_NEAR(41.1017391164415, mesh.m_nodes[8].y, tolerance);
    ASSERT_NEAR(41.1017395852698, mesh.m_nodes[9].y, tolerance);
    ASSERT_NEAR(41.1017400000000, mesh.m_nodes[10].y, tolerance);
    ASSERT_NEAR(41.1017404147302, mesh.m_nodes[11].y, tolerance);
    ASSERT_NEAR(41.1017408789613, mesh.m_nodes[12].y, tolerance);
    ASSERT_NEAR(41.1017413985836, mesh.m_nodes[13].y, tolerance);
    ASSERT_NEAR(41.1017419801843, mesh.m_nodes[14].y, tolerance);
    ASSERT_NEAR(41.1017426311279, mesh.m_nodes[15].y, tolerance);
    ASSERT_NEAR(41.1017433596473, mesh.m_nodes[16].y, tolerance);
    ASSERT_NEAR(41.1017441749444, mesh.m_nodes[17].y, tolerance);
    ASSERT_NEAR(41.1017450873013, mesh.m_nodes[18].y, tolerance);
    ASSERT_NEAR(41.1017461082051, mesh.m_nodes[19].y, tolerance);
}

TEST(CurvilinearGridFromSplines, OrthogonalCurvilinearMeshTwoCrossingHighCurvature)
{
    std::vector<meshkernel::Point> firstSpline{{-103.468336664918, 420.606335102062},
                                               {111.984583950396, 845.689124424167},
                                               {1465.84415268176, 1608.50892444055}};

    auto splines = std::make_shared<meshkernel::Splines>(meshkernel::Projections::cartesian);
    splines->AddSpline(firstSpline, 0, firstSpline.size());

    std::vector<meshkernel::Point> secondSpline{{-333.478887051536, 921.388799234953},
                                                {167.303577081355, 490.482958004326}};

    splines->AddSpline(secondSpline, 0, secondSpline.size());

    meshkernelapi::CurvilinearParametersNative curvilinearParametersNative;
    meshkernelapi::SplinesToCurvilinearParametersNative splinesToCurvilinearParametersNative;

    splinesToCurvilinearParametersNative.AspectRatio = 0.5;
    splinesToCurvilinearParametersNative.AspectRatioGrowFactor = 1.1;
    splinesToCurvilinearParametersNative.AverageWidth = 120.0;
    splinesToCurvilinearParametersNative.GridsOnTopOfEachOtherTolerance = 1e-4;
    splinesToCurvilinearParametersNative.MinimumCosineOfCrossingAngles = 0.95;
    splinesToCurvilinearParametersNative.CheckFrontCollisions = false;
    splinesToCurvilinearParametersNative.CurvatureAdapetedGridSpacing = true;
    splinesToCurvilinearParametersNative.RemoveSkinnyTriangles = 0;
    curvilinearParametersNative.MRefinement = 20;
    curvilinearParametersNative.NRefinement = 3;

    // create the algorithm
    meshkernel::CurvilinearGridFromSplines curvilinearGridFromSplines(splines, curvilinearParametersNative, splinesToCurvilinearParametersNative);

    // compute
    meshkernel::CurvilinearGrid curvilinearGrid;
    curvilinearGridFromSplines.Compute(curvilinearGrid);

    meshkernel::Mesh mesh(curvilinearGrid, meshkernel::Projections::cartesian);

    const double tolerance = 1e-6;

    ASSERT_NEAR(189.317899173845, mesh.m_nodes[0].x, tolerance);
    ASSERT_NEAR(47.4845436454654, mesh.m_nodes[1].x, tolerance);
    ASSERT_NEAR(-44.2005274315104, mesh.m_nodes[2].x, tolerance);
    ASSERT_NEAR(-103.468336664918, mesh.m_nodes[3].x, tolerance);
    ASSERT_NEAR(-162.736145898326, mesh.m_nodes[4].x, tolerance);
    ASSERT_NEAR(-269.814914053481, mesh.m_nodes[5].x, tolerance);
    ASSERT_NEAR(-463.273432900883, mesh.m_nodes[6].x, tolerance);
    ASSERT_NEAR(182.838365135208, mesh.m_nodes[7].x, tolerance);
    ASSERT_NEAR(37.2764147141659, mesh.m_nodes[8].x, tolerance);
    ASSERT_NEAR(-55.7300943458967, mesh.m_nodes[9].x, tolerance);
    ASSERT_NEAR(-115.852118266335, mesh.m_nodes[10].x, tolerance);
    ASSERT_NEAR(-175.974142186772, mesh.m_nodes[11].x, tolerance);
    ASSERT_NEAR(-284.596214499957, mesh.m_nodes[12].x, tolerance);
    ASSERT_NEAR(-480.843011002095, mesh.m_nodes[13].x, tolerance);
    ASSERT_NEAR(182.838365135208, mesh.m_nodes[14].x, tolerance);
    ASSERT_NEAR(37.8679827389320, mesh.m_nodes[15].x, tolerance);
    ASSERT_NEAR(-54.6628549419648, mesh.m_nodes[16].x, tolerance);
    ASSERT_NEAR(-114.477391501297, mesh.m_nodes[17].x, tolerance);
    ASSERT_NEAR(-174.291928060630, mesh.m_nodes[18].x, tolerance);
    ASSERT_NEAR(-282.358464944443, mesh.m_nodes[19].x, tolerance);
    ASSERT_NEAR(-477.601579174313, mesh.m_nodes[20].x, tolerance);
    ASSERT_NEAR(182.838365135208, mesh.m_nodes[21].x, tolerance);
    ASSERT_NEAR(42.1724119247128, mesh.m_nodes[22].x, tolerance);
    ASSERT_NEAR(-45.7933216922166, mesh.m_nodes[23].x, tolerance);
    ASSERT_NEAR(-102.656846558853, mesh.m_nodes[24].x, tolerance);

    ASSERT_NEAR(466.770557740332, mesh.m_nodes[0].y, tolerance);
    ASSERT_NEAR(444.407393763762, mesh.m_nodes[1].y, tolerance);
    ASSERT_NEAR(429.951215447599, mesh.m_nodes[2].y, tolerance);
    ASSERT_NEAR(420.606335102062, mesh.m_nodes[3].y, tolerance);
    ASSERT_NEAR(411.261454756525, mesh.m_nodes[4].y, tolerance);
    ASSERT_NEAR(394.378119791050, mesh.m_nodes[5].y, tolerance);
    ASSERT_NEAR(363.875107551570, mesh.m_nodes[6].y, tolerance);
    ASSERT_NEAR(507.865550198607, mesh.m_nodes[7].y, tolerance);
    ASSERT_NEAR(509.150163682609, mesh.m_nodes[8].y, tolerance);
    ASSERT_NEAR(503.074909615816, mesh.m_nodes[9].y, tolerance);
    ASSERT_NEAR(499.147694477674, mesh.m_nodes[10].y, tolerance);
    ASSERT_NEAR(495.220479339532, mesh.m_nodes[11].y, tolerance);
    ASSERT_NEAR(488.125205113747, mesh.m_nodes[12].y, tolerance);
    ASSERT_NEAR(475.306218995532, mesh.m_nodes[13].y, tolerance);
    ASSERT_NEAR(507.865550198607, mesh.m_nodes[14].y, tolerance);
    ASSERT_NEAR(531.960403892052, mesh.m_nodes[15].y, tolerance);
    ASSERT_NEAR(544.226537400807, mesh.m_nodes[16].y, tolerance);
    ASSERT_NEAR(552.155711172380, mesh.m_nodes[17].y, tolerance);
    ASSERT_NEAR(560.084884943953, mesh.m_nodes[18].y, tolerance);
    ASSERT_NEAR(574.410471985723, mesh.m_nodes[19].y, tolerance);
    ASSERT_NEAR(600.292417570909, mesh.m_nodes[20].y, tolerance);
    ASSERT_NEAR(507.865550198607, mesh.m_nodes[21].y, tolerance);
    ASSERT_NEAR(549.729828974096, mesh.m_nodes[22].y, tolerance);
    ASSERT_NEAR(580.841498361190, mesh.m_nodes[23].y, tolerance);
    ASSERT_NEAR(600.952956688129, mesh.m_nodes[24].y, tolerance);
}

TEST(CurvilinearGridFromSplines, OrthogonalCurvilinearMeshTwoCrossingHighCurvatureRemoveSkinnyTriangles)
{
    std::vector<meshkernel::Point> firstSpline{{-103.468336664918, 420.606335102062},
                                               {111.984583950396, 845.689124424167},
                                               {1465.84415268176, 1608.50892444055}};

    auto splines = std::make_shared<meshkernel::Splines>(meshkernel::Projections::cartesian);
    splines->AddSpline(firstSpline, 0, firstSpline.size());

    std::vector<meshkernel::Point> secondSpline{{-333.478887051536, 921.388799234953},
                                                {167.303577081355, 490.482958004326}};

    splines->AddSpline(secondSpline, 0, secondSpline.size());

    meshkernelapi::CurvilinearParametersNative curvilinearParametersNative;
    meshkernelapi::SplinesToCurvilinearParametersNative splinesToCurvilinearParametersNative;

    splinesToCurvilinearParametersNative.AspectRatio = 0.5;
    splinesToCurvilinearParametersNative.AspectRatioGrowFactor = 1.1;
    splinesToCurvilinearParametersNative.AverageWidth = 120.0;
    splinesToCurvilinearParametersNative.GridsOnTopOfEachOtherTolerance = 1e-4;
    splinesToCurvilinearParametersNative.MinimumCosineOfCrossingAngles = 0.95;
    splinesToCurvilinearParametersNative.CheckFrontCollisions = false;
    splinesToCurvilinearParametersNative.CurvatureAdapetedGridSpacing = true;
    splinesToCurvilinearParametersNative.RemoveSkinnyTriangles = true;
    curvilinearParametersNative.MRefinement = 20;
    curvilinearParametersNative.NRefinement = 3;

    // create the algorithm
    meshkernel::CurvilinearGridFromSplines curvilinearGridFromSplines(splines, curvilinearParametersNative, splinesToCurvilinearParametersNative);

    // compute
    meshkernel::CurvilinearGrid curvilinearGrid;
    curvilinearGridFromSplines.Compute(curvilinearGrid);

    meshkernel::Mesh mesh(curvilinearGrid, meshkernel::Projections::cartesian);

    const double tolerance = 1e-6;

    ASSERT_NEAR(189.317899173845, mesh.m_nodes[0].x, tolerance);
    ASSERT_NEAR(47.4845436454654, mesh.m_nodes[1].x, tolerance);
    ASSERT_NEAR(-44.2005274315104, mesh.m_nodes[2].x, tolerance);
    ASSERT_NEAR(-103.468336664918, mesh.m_nodes[3].x, tolerance);
    ASSERT_NEAR(-162.736145898326, mesh.m_nodes[4].x, tolerance);
    ASSERT_NEAR(-269.814914053481, mesh.m_nodes[5].x, tolerance);
    ASSERT_NEAR(-463.273432900883, mesh.m_nodes[6].x, tolerance);
    ASSERT_NEAR(182.838365135208, mesh.m_nodes[7].x, tolerance);
    ASSERT_NEAR(42.7645089001448, mesh.m_nodes[8].x, tolerance);
    ASSERT_NEAR(-55.7300943458967, mesh.m_nodes[9].x, tolerance);
    ASSERT_NEAR(-115.852118266335, mesh.m_nodes[10].x, tolerance);
    ASSERT_NEAR(-175.974142186772, mesh.m_nodes[11].x, tolerance);
    ASSERT_NEAR(-284.596214499957, mesh.m_nodes[12].x, tolerance);
    ASSERT_NEAR(-480.843011002095, mesh.m_nodes[13].x, tolerance);
    ASSERT_NEAR(182.838365135208, mesh.m_nodes[14].x, tolerance);
    ASSERT_NEAR(42.7645089001448, mesh.m_nodes[15].x, tolerance);
    ASSERT_NEAR(-54.6628549419648, mesh.m_nodes[16].x, tolerance);
    ASSERT_NEAR(-114.477391501297, mesh.m_nodes[17].x, tolerance);
    ASSERT_NEAR(-174.291928060630, mesh.m_nodes[18].x, tolerance);
    ASSERT_NEAR(-282.358464944443, mesh.m_nodes[19].x, tolerance);
    ASSERT_NEAR(-477.601579174313, mesh.m_nodes[20].x, tolerance);
    ASSERT_NEAR(182.838365135208, mesh.m_nodes[21].x, tolerance);
    ASSERT_NEAR(42.7645089001448, mesh.m_nodes[22].x, tolerance);
    ASSERT_NEAR(-45.7933216922166, mesh.m_nodes[23].x, tolerance);
    ASSERT_NEAR(-102.656846558853, mesh.m_nodes[24].x, tolerance);

    ASSERT_NEAR(466.770557740332, mesh.m_nodes[0].y, tolerance);
    ASSERT_NEAR(444.407393763762, mesh.m_nodes[1].y, tolerance);
    ASSERT_NEAR(429.951215447599, mesh.m_nodes[2].y, tolerance);
    ASSERT_NEAR(420.606335102062, mesh.m_nodes[3].y, tolerance);
    ASSERT_NEAR(411.261454756525, mesh.m_nodes[4].y, tolerance);
    ASSERT_NEAR(394.378119791050, mesh.m_nodes[5].y, tolerance);
    ASSERT_NEAR(363.875107551570, mesh.m_nodes[6].y, tolerance);
    ASSERT_NEAR(507.865550198607, mesh.m_nodes[7].y, tolerance);
    ASSERT_NEAR(541.250471813772, mesh.m_nodes[8].y, tolerance);
    ASSERT_NEAR(503.074909615816, mesh.m_nodes[9].y, tolerance);
    ASSERT_NEAR(499.147694477674, mesh.m_nodes[10].y, tolerance);
    ASSERT_NEAR(495.220479339532, mesh.m_nodes[11].y, tolerance);
    ASSERT_NEAR(488.125205113747, mesh.m_nodes[12].y, tolerance);
    ASSERT_NEAR(475.306218995532, mesh.m_nodes[13].y, tolerance);
    ASSERT_NEAR(507.865550198607, mesh.m_nodes[14].y, tolerance);
    ASSERT_NEAR(541.250471813772, mesh.m_nodes[15].y, tolerance);
    ASSERT_NEAR(544.226537400807, mesh.m_nodes[16].y, tolerance);
    ASSERT_NEAR(552.155711172380, mesh.m_nodes[17].y, tolerance);
    ASSERT_NEAR(560.084884943953, mesh.m_nodes[18].y, tolerance);
    ASSERT_NEAR(574.410471985723, mesh.m_nodes[19].y, tolerance);
    ASSERT_NEAR(600.292417570909, mesh.m_nodes[20].y, tolerance);
    ASSERT_NEAR(507.865550198607, mesh.m_nodes[21].y, tolerance);
    ASSERT_NEAR(541.250471813772, mesh.m_nodes[22].y, tolerance);
    ASSERT_NEAR(580.841498361190, mesh.m_nodes[23].y, tolerance);
    ASSERT_NEAR(600.952956688129, mesh.m_nodes[24].y, tolerance);
}
