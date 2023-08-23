#include <gtest/gtest.h>

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridFromSplines.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Parameters.hpp>
#include <MeshKernel/Splines.hpp>

using namespace meshkernel;

TEST(CurvilinearGridFromSplines, ComputeSplinesProperties)
{
    std::vector<Point> firstSpline{{152.001571655273, 86.6264953613281},
                                   {374.752960205078, 336.378997802734},
                                   {850.255920410156, 499.130676269531}};

    auto splines = std::make_shared<Splines>(Projection::cartesian);
    splines->AddSpline(firstSpline, 0, static_cast<UInt>(firstSpline.size()));

    std::vector<Point> secondSpline{{72.5010681152344, 391.129577636719},
                                    {462.503479003906, 90.3765411376953}};

    splines->AddSpline(secondSpline, 0, static_cast<UInt>(secondSpline.size()));

    // construct CurvilinearGridFromSplines
    CurvilinearParameters curvilinearParameters;
    SplinesToCurvilinearParameters splinesToCurvilinearParameters;

    splinesToCurvilinearParameters.aspect_ratio = 0.1;
    splinesToCurvilinearParameters.aspect_ratio_grow_factor = 1.1;
    splinesToCurvilinearParameters.average_width = 500.0;
    splinesToCurvilinearParameters.nodes_on_top_of_each_other_tolerance = 1e-4;
    splinesToCurvilinearParameters.min_cosine_crossing_angles = 0.95;
    splinesToCurvilinearParameters.check_front_collisions = false;
    splinesToCurvilinearParameters.curvature_adapted_grid_spacing = true;
    splinesToCurvilinearParameters.remove_skinny_triangles = 0;
    curvilinearParameters.m_refinement = 20;
    curvilinearParameters.n_refinement = 40;

    CurvilinearGridFromSplines curvilinearGridFromSplines(splines, curvilinearParameters, splinesToCurvilinearParameters);

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
    std::vector<Point> firstSpline{{91.2511901855469, 299.628631591797},
                                   {354.502838134766, 518.630859375000},
                                   {770.755432128906, 607.881774902344}};

    auto splines = std::make_shared<Splines>(Projection::cartesian);
    splines->AddSpline(firstSpline, 0, static_cast<UInt>(firstSpline.size()));

    std::vector<Point> secondSpline;
    secondSpline.emplace_back(Point{273.502319335938, 86.6264953613281});
    secondSpline.emplace_back(Point{557.004089355469, 316.128814697266});
    secondSpline.emplace_back(Point{847.255920410156, 409.129730224609});
    splines->AddSpline(secondSpline, 0, static_cast<UInt>(secondSpline.size()));

    std::vector<Point> thirdSpline;
    thirdSpline.emplace_back(Point{62.7510070800781, 396.379608154297});
    thirdSpline.emplace_back(Point{350.752807617188, 73.8763732910156});
    splines->AddSpline(thirdSpline, 0, static_cast<UInt>(thirdSpline.size()));

    std::vector<Point> fourthSpline;
    fourthSpline.emplace_back(Point{704.755004882812, 636.382019042969});
    fourthSpline.emplace_back(Point{845.005859375000, 285.378509521484});
    splines->AddSpline(fourthSpline, 0, static_cast<UInt>(fourthSpline.size()));

    // construct CurvilinearGridFromSplines
    CurvilinearParameters curvilinearParameters;
    SplinesToCurvilinearParameters splinesToCurvilinearParameters;

    splinesToCurvilinearParameters.aspect_ratio = 0.1;
    splinesToCurvilinearParameters.aspect_ratio_grow_factor = 1.1;
    splinesToCurvilinearParameters.average_width = 500.0;
    splinesToCurvilinearParameters.nodes_on_top_of_each_other_tolerance = 1e-4;
    splinesToCurvilinearParameters.min_cosine_crossing_angles = 0.95;
    splinesToCurvilinearParameters.check_front_collisions = false;
    splinesToCurvilinearParameters.curvature_adapted_grid_spacing = false;
    splinesToCurvilinearParameters.remove_skinny_triangles = 0;
    curvilinearParameters.m_refinement = 20;
    curvilinearParameters.n_refinement = 40;

    CurvilinearGridFromSplines curvilinearGridFromSplines(splines, curvilinearParameters, splinesToCurvilinearParameters);
    curvilinearGridFromSplines.ComputeSplineProperties(false);

    const double tolerance = 1e-6;
    ASSERT_NEAR(345.967070088532, curvilinearGridFromSplines.m_maximumGridHeights[0], tolerance);
    ASSERT_NEAR(370.339417298715, curvilinearGridFromSplines.m_maximumGridHeights[1], tolerance);
    ASSERT_NEAR(0.0, curvilinearGridFromSplines.m_maximumGridHeights[2], tolerance);
    ASSERT_NEAR(0.0, curvilinearGridFromSplines.m_maximumGridHeights[3], tolerance);
}

TEST(CurvilinearGridFromSplines, OrthogonalCurvilinearMeshTwoCrossingCurvatureAdapted)
{
    std::vector<Point> firstSpline{{26.9847544156713, 327.556992634968},
                                   {340.139086327869, 819.656657068422},
                                   {2048.50780774173, 1644.48279915859}};

    auto splines = std::make_shared<Splines>(Projection::cartesian);
    splines->AddSpline(firstSpline, 0, static_cast<UInt>(firstSpline.size()));

    std::vector<Point> secondSpline;
    secondSpline.emplace_back(Point{-179.920786312031, 1068.50251010579});
    secondSpline.emplace_back(Point{600.169022647819, 321.964950993679});
    splines->AddSpline(secondSpline, 0, static_cast<UInt>(secondSpline.size()));

    // construct CurvilinearGridFromSplines
    CurvilinearParameters curvilinearParameters;
    SplinesToCurvilinearParameters splinesToCurvilinearParameters;

    splinesToCurvilinearParameters.aspect_ratio = 0.1;
    splinesToCurvilinearParameters.aspect_ratio_grow_factor = 1.1;
    splinesToCurvilinearParameters.average_width = 500.0;
    splinesToCurvilinearParameters.nodes_on_top_of_each_other_tolerance = 1e-4;
    splinesToCurvilinearParameters.min_cosine_crossing_angles = 0.95;
    splinesToCurvilinearParameters.check_front_collisions = false;
    splinesToCurvilinearParameters.curvature_adapted_grid_spacing = true;
    splinesToCurvilinearParameters.remove_skinny_triangles = 0;
    curvilinearParameters.m_refinement = 20;
    curvilinearParameters.n_refinement = 40;

    CurvilinearGridFromSplines curvilinearGridFromSplines(splines, curvilinearParameters, splinesToCurvilinearParameters);
    auto curviGrid = curvilinearGridFromSplines.Compute();

    const double tolerance = 1e-6;
    ASSERT_NEAR(588.14274314862507, curviGrid.m_gridNodes(0, 0).x, tolerance);
    ASSERT_NEAR(469.41412484972938, curviGrid.m_gridNodes(0, 1).x, tolerance);
    ASSERT_NEAR(366.68736858428133, curviGrid.m_gridNodes(0, 2).x, tolerance);
    ASSERT_NEAR(277.80579475635676, curviGrid.m_gridNodes(0, 3).x, tolerance);
    ASSERT_NEAR(200.90339412262884, curviGrid.m_gridNodes(0, 4).x, tolerance);
    ASSERT_NEAR(134.36565217535070, curviGrid.m_gridNodes(0, 5).x, tolerance);
    ASSERT_NEAR(76.795653489238660, curviGrid.m_gridNodes(0, 6).x, tolerance);
    ASSERT_NEAR(26.984754415671301, curviGrid.m_gridNodes(0, 7).x, tolerance);
    ASSERT_NEAR(-22.826144657896062, curviGrid.m_gridNodes(0, 8).x, tolerance);
    ASSERT_NEAR(-78.994845354601750, curviGrid.m_gridNodes(0, 9).x, tolerance);
    ASSERT_NEAR(-142.33284961942132, curviGrid.m_gridNodes(0, 10).x, tolerance);
    ASSERT_NEAR(-213.75523850308639, curviGrid.m_gridNodes(0, 11).x, tolerance);
    ASSERT_NEAR(-294.29389287118539, curviGrid.m_gridNodes(0, 12).x, tolerance);
    ASSERT_NEAR(-385.11240158824154, curviGrid.m_gridNodes(0, 13).x, tolerance);
    ASSERT_NEAR(-487.52287256398955, curviGrid.m_gridNodes(0, 14).x, tolerance);

    ASSERT_NEAR(278.61331701020396, curviGrid.m_gridNodes(0, 0).y, tolerance);
    ASSERT_NEAR(288.96871619915311, curviGrid.m_gridNodes(0, 1).y, tolerance);
    ASSERT_NEAR(297.92844792816703, curviGrid.m_gridNodes(0, 2).y, tolerance);
    ASSERT_NEAR(305.68061577714451, curviGrid.m_gridNodes(0, 3).y, tolerance);
    ASSERT_NEAR(312.38797132933286, curviGrid.m_gridNodes(0, 4).y, tolerance);
    ASSERT_NEAR(318.19133103205218, curviGrid.m_gridNodes(0, 5).y, tolerance);
    ASSERT_NEAR(323.21253254389961, curviGrid.m_gridNodes(0, 6).y, tolerance);
    ASSERT_NEAR(327.55699263496803, curviGrid.m_gridNodes(0, 7).y, tolerance);
    ASSERT_NEAR(331.90145272603644, curviGrid.m_gridNodes(0, 8).y, tolerance);
    ASSERT_NEAR(336.80043433613747, curviGrid.m_gridNodes(0, 9).y, tolerance);
    ASSERT_NEAR(342.32471590692433, curviGrid.m_gridNodes(0, 10).y, tolerance);
    ASSERT_NEAR(348.55410995290555, curviGrid.m_gridNodes(0, 11).y, tolerance);
    ASSERT_NEAR(355.57861615934416, curviGrid.m_gridNodes(0, 12).y, tolerance);
    ASSERT_NEAR(363.49972166014766, curviGrid.m_gridNodes(0, 13).y, tolerance);
    ASSERT_NEAR(372.43186728162192, curviGrid.m_gridNodes(0, 14).y, tolerance);
}

TEST(CurvilinearGridFromSplines, OrthogonalCurvilinearMeshTwoCrossingCurvatureNotAdapted)
{
    std::vector<Point> firstSpline{{26.9847544156713, 327.556992634968},
                                   {340.139086327869, 819.656657068422},
                                   {2048.50780774173, 1644.48279915859}};

    auto splines = std::make_shared<Splines>(Projection::cartesian);
    splines->AddSpline(firstSpline, 0, static_cast<UInt>(firstSpline.size()));

    std::vector<Point> secondSpline;
    secondSpline.emplace_back(Point{-179.920786312031, 1068.50251010579});
    secondSpline.emplace_back(Point{600.169022647819, 321.964950993679});
    splines->AddSpline(secondSpline, 0, static_cast<UInt>(secondSpline.size()));

    CurvilinearParameters curvilinearParameters;
    SplinesToCurvilinearParameters splinesToCurvilinearParameters;

    splinesToCurvilinearParameters.aspect_ratio = 0.1;
    splinesToCurvilinearParameters.aspect_ratio_grow_factor = 1.1;
    splinesToCurvilinearParameters.average_width = 500.0;
    splinesToCurvilinearParameters.nodes_on_top_of_each_other_tolerance = 1e-4;
    splinesToCurvilinearParameters.min_cosine_crossing_angles = 0.95;
    splinesToCurvilinearParameters.check_front_collisions = false;
    splinesToCurvilinearParameters.curvature_adapted_grid_spacing = false;
    splinesToCurvilinearParameters.remove_skinny_triangles = 0;
    curvilinearParameters.m_refinement = 20;
    curvilinearParameters.n_refinement = 40;

    CurvilinearGridFromSplines curvilinearGridFromSplines(splines, curvilinearParameters, splinesToCurvilinearParameters);

    const auto curviGrid = curvilinearGridFromSplines.Compute();

    const double tolerance = 1e-6;
    ASSERT_NEAR(548.64105210377159, curviGrid.m_gridNodes(0, 0).x, tolerance);
    ASSERT_NEAR(438.27011790167523, curviGrid.m_gridNodes(0, 1).x, tolerance);
    ASSERT_NEAR(342.77462390017541, curviGrid.m_gridNodes(0, 2).x, tolerance);
    ASSERT_NEAR(260.14970602195672, curviGrid.m_gridNodes(0, 3).x, tolerance);
    ASSERT_NEAR(188.66070933405891, curviGrid.m_gridNodes(0, 4).x, tolerance);
    ASSERT_NEAR(126.80677012682224, curviGrid.m_gridNodes(0, 5).x, tolerance);
    ASSERT_NEAR(73.289306283324080, curviGrid.m_gridNodes(0, 6).x, tolerance);
    ASSERT_NEAR(26.984754415671301, curviGrid.m_gridNodes(0, 7).x, tolerance);
    ASSERT_NEAR(-19.319797451981476, curviGrid.m_gridNodes(0, 8).x, tolerance);
    ASSERT_NEAR(-71.534605117026189, curviGrid.m_gridNodes(0, 9).x, tolerance);
    ASSERT_NEAR(-130.41404632812231, curviGrid.m_gridNodes(0, 10).x, tolerance);
    ASSERT_NEAR(-196.80878667783674, curviGrid.m_gridNodes(0, 11).x, tolerance);
    ASSERT_NEAR(-271.67806966408966, curviGrid.m_gridNodes(0, 12).x, tolerance);
    ASSERT_NEAR(-356.10357543986055, curviGrid.m_gridNodes(0, 13).x, tolerance);
    ASSERT_NEAR(-451.30504847658597, curviGrid.m_gridNodes(0, 14).x, tolerance);

    ASSERT_NEAR(115.02822097141149, curviGrid.m_gridNodes(0, 0).y, tolerance);
    ASSERT_NEAR(159.99460818323763, curviGrid.m_gridNodes(0, 1).y, tolerance);
    ASSERT_NEAR(198.90057007709373, curviGrid.m_gridNodes(0, 2).y, tolerance);
    ASSERT_NEAR(232.56291132224692, curviGrid.m_gridNodes(0, 3).y, tolerance);
    ASSERT_NEAR(261.68835027673492, curviGrid.m_gridNodes(0, 4).y, tolerance);
    ASSERT_NEAR(286.88835606699604, curviGrid.m_gridNodes(0, 5).y, tolerance);
    ASSERT_NEAR(308.69198597353653, curviGrid.m_gridNodes(0, 6).y, tolerance);
    ASSERT_NEAR(327.55699263496803, curviGrid.m_gridNodes(0, 7).y, tolerance);
    ASSERT_NEAR(346.42199929639952, curviGrid.m_gridNodes(0, 8).y, tolerance);
    ASSERT_NEAR(367.69491210293836, curviGrid.m_gridNodes(0, 9).y, tolerance);
    ASSERT_NEAR(391.68307322030449, curviGrid.m_gridNodes(0, 10).y, tolerance);
    ASSERT_NEAR(418.73305358857078, curviGrid.m_gridNodes(0, 11).y, tolerance);
    ASSERT_NEAR(449.23566003445961, curviGrid.m_gridNodes(0, 12).y, tolerance);
    ASSERT_NEAR(483.63158148526736, curviGrid.m_gridNodes(0, 13).y, tolerance);
    ASSERT_NEAR(522.41775585855908, curviGrid.m_gridNodes(0, 14).y, tolerance);
}

TEST(CurvilinearGridFromSplines, OrthogonalCurvilinearMeshTwoCrossingCurvatureAdaptedLargeMRefinement)
{
    std::vector<Point> firstSpline{{26.9847544156713, 327.556992634968},
                                   {340.139086327869, 819.656657068422},
                                   {2048.50780774173, 1644.48279915859}};

    auto splines = std::make_shared<Splines>(Projection::cartesian);
    splines->AddSpline(firstSpline, 0, static_cast<UInt>(firstSpline.size()));

    std::vector<Point> secondSpline;
    secondSpline.emplace_back(Point{-179.920786312031, 1068.50251010579});
    secondSpline.emplace_back(Point{600.169022647819, 321.964950993679});
    splines->AddSpline(secondSpline, 0, static_cast<UInt>(secondSpline.size()));

    CurvilinearParameters curvilinearParameters;
    SplinesToCurvilinearParameters splinesToCurvilinearParameters;

    splinesToCurvilinearParameters.aspect_ratio = 0.1;
    splinesToCurvilinearParameters.aspect_ratio_grow_factor = 1.1;
    splinesToCurvilinearParameters.average_width = 500.0;
    splinesToCurvilinearParameters.nodes_on_top_of_each_other_tolerance = 1e-4;
    splinesToCurvilinearParameters.min_cosine_crossing_angles = 0.95;
    splinesToCurvilinearParameters.check_front_collisions = false;
    splinesToCurvilinearParameters.curvature_adapted_grid_spacing = true;
    splinesToCurvilinearParameters.remove_skinny_triangles = 0;
    curvilinearParameters.m_refinement = 10;
    curvilinearParameters.n_refinement = 20;

    // create the algorithm
    CurvilinearGridFromSplines curvilinearGridFromSplines(splines, curvilinearParameters, splinesToCurvilinearParameters);

    // compute
    const auto curviGrid = curvilinearGridFromSplines.Compute();

    const double tolerance = 1e-6;
    ASSERT_NEAR(588.14274314862507, curviGrid.m_gridNodes(0, 0).x, tolerance);
    ASSERT_NEAR(469.41412484972938, curviGrid.m_gridNodes(0, 1).x, tolerance);
    ASSERT_NEAR(366.68736858428133, curviGrid.m_gridNodes(0, 2).x, tolerance);
    ASSERT_NEAR(277.80579475635676, curviGrid.m_gridNodes(0, 3).x, tolerance);
    ASSERT_NEAR(200.90339412262884, curviGrid.m_gridNodes(0, 4).x, tolerance);
    ASSERT_NEAR(134.36565217535070, curviGrid.m_gridNodes(0, 5).x, tolerance);
    ASSERT_NEAR(76.795653489238660, curviGrid.m_gridNodes(0, 6).x, tolerance);
    ASSERT_NEAR(26.984754415671301, curviGrid.m_gridNodes(0, 7).x, tolerance);
    ASSERT_NEAR(-22.826144657896062, curviGrid.m_gridNodes(0, 8).x, tolerance);
    ASSERT_NEAR(-78.994845354601750, curviGrid.m_gridNodes(0, 9).x, tolerance);
    ASSERT_NEAR(-142.33284961942132, curviGrid.m_gridNodes(0, 10).x, tolerance);
    ASSERT_NEAR(-213.75523850308639, curviGrid.m_gridNodes(0, 11).x, tolerance);
    ASSERT_NEAR(-294.29389287118539, curviGrid.m_gridNodes(0, 12).x, tolerance);
    ASSERT_NEAR(-385.11240158824154, curviGrid.m_gridNodes(0, 13).x, tolerance);
    ASSERT_NEAR(-487.52287256398955, curviGrid.m_gridNodes(0, 14).x, tolerance);

    ASSERT_NEAR(278.61331701020396, curviGrid.m_gridNodes(0, 0).y, tolerance);
    ASSERT_NEAR(288.96871619915311, curviGrid.m_gridNodes(0, 1).y, tolerance);
    ASSERT_NEAR(297.92844792816703, curviGrid.m_gridNodes(0, 2).y, tolerance);
    ASSERT_NEAR(305.68061577714451, curviGrid.m_gridNodes(0, 3).y, tolerance);
    ASSERT_NEAR(312.38797132933286, curviGrid.m_gridNodes(0, 4).y, tolerance);
    ASSERT_NEAR(318.19133103205218, curviGrid.m_gridNodes(0, 5).y, tolerance);
    ASSERT_NEAR(323.21253254389961, curviGrid.m_gridNodes(0, 6).y, tolerance);
    ASSERT_NEAR(327.55699263496803, curviGrid.m_gridNodes(0, 7).y, tolerance);
    ASSERT_NEAR(331.90145272603644, curviGrid.m_gridNodes(0, 8).y, tolerance);
    ASSERT_NEAR(336.80043433613747, curviGrid.m_gridNodes(0, 9).y, tolerance);
    ASSERT_NEAR(342.32471590692433, curviGrid.m_gridNodes(0, 10).y, tolerance);
    ASSERT_NEAR(348.55410995290555, curviGrid.m_gridNodes(0, 11).y, tolerance);
    ASSERT_NEAR(355.57861615934416, curviGrid.m_gridNodes(0, 12).y, tolerance);
    ASSERT_NEAR(363.49972166014766, curviGrid.m_gridNodes(0, 13).y, tolerance);
    ASSERT_NEAR(372.43186728162192, curviGrid.m_gridNodes(0, 14).y, tolerance);

    ASSERT_NEAR(589.75167441442306, curviGrid.m_gridNodes(1, 0).x, tolerance);
    ASSERT_NEAR(475.48540012613313, curviGrid.m_gridNodes(1, 1).x, tolerance);
    ASSERT_NEAR(375.29368253913503, curviGrid.m_gridNodes(1, 2).x, tolerance);
    ASSERT_NEAR(288.60548301452314, curviGrid.m_gridNodes(1, 3).x, tolerance);
    ASSERT_NEAR(213.60084073664558, curviGrid.m_gridNodes(1, 4).x, tolerance);
    ASSERT_NEAR(148.70508336617249, curviGrid.m_gridNodes(1, 5).x, tolerance);
    ASSERT_NEAR(92.555767846678634, curviGrid.m_gridNodes(1, 6).x, tolerance);
    ASSERT_NEAR(43.974076823189009, curviGrid.m_gridNodes(1, 7).x, tolerance);
    ASSERT_NEAR(-4.6076142003006169, curviGrid.m_gridNodes(1, 8).x, tolerance);
    ASSERT_NEAR(-59.390212249936532, curviGrid.m_gridNodes(1, 9).x, tolerance);
    ASSERT_NEAR(-121.16519343815622, curviGrid.m_gridNodes(1, 10).x, tolerance);
    ASSERT_NEAR(-190.82505691062727, curviGrid.m_gridNodes(1, 11).x, tolerance);
    ASSERT_NEAR(-269.37621930140858, curviGrid.m_gridNodes(1, 12).x, tolerance);
    ASSERT_NEAR(-357.95355502043907, curviGrid.m_gridNodes(1, 13).x, tolerance);
    ASSERT_NEAR(-457.83679244535881, curviGrid.m_gridNodes(1, 14).x, tolerance);

    ASSERT_NEAR(297.06033026778510, curviGrid.m_gridNodes(1, 0).y, tolerance);
    ASSERT_NEAR(358.57821282279917, curviGrid.m_gridNodes(1, 1).y, tolerance);
    ASSERT_NEAR(396.60313403748916, curviGrid.m_gridNodes(1, 2).y, tolerance);
    ASSERT_NEAR(429.50317843791038, curviGrid.m_gridNodes(1, 3).y, tolerance);
    ASSERT_NEAR(457.96906046928518, curviGrid.m_gridNodes(1, 4).y, tolerance);
    ASSERT_NEAR(482.59840230123029, curviGrid.m_gridNodes(1, 5).y, tolerance);
    ASSERT_NEAR(503.90828050581950, curviGrid.m_gridNodes(1, 6).y, tolerance);
    ASSERT_NEAR(522.34608173465415, curviGrid.m_gridNodes(1, 7).y, tolerance);
    ASSERT_NEAR(540.78388296348885, curviGrid.m_gridNodes(1, 8).y, tolerance);
    ASSERT_NEAR(561.57506236376855, curviGrid.m_gridNodes(1, 9).y, tolerance);
    ASSERT_NEAR(585.02000221801700, curviGrid.m_gridNodes(1, 10).y, tolerance);
    ASSERT_NEAR(611.45742523221747, curviGrid.m_gridNodes(1, 11).y, tolerance);
    ASSERT_NEAR(641.26928826009907, curviGrid.m_gridNodes(1, 12).y, tolerance);
    ASSERT_NEAR(674.88630065634504, curviGrid.m_gridNodes(1, 13).y, tolerance);
    ASSERT_NEAR(712.79414698558730, curviGrid.m_gridNodes(1, 14).y, tolerance);
}

TEST(CurvilinearGridFromSplines, OrthogonalCurvilinearMeshFourSplineCrossingFront)
{
    std::vector<Point> firstSpline{{-93.379753864775, 231.718492951018},
                                   {72.8139687226708, 724.468302026077},
                                   {746.335897103372, 234.634172294657},
                                   {1498.58116776234, 776.950530211586}};

    auto splines = std::make_shared<Splines>(Projection::cartesian);
    splines->AddSpline(firstSpline, 0, static_cast<UInt>(firstSpline.size()));

    std::vector<Point> secondSpline{{-250.826438421303, 394.996536194825},
                                    {101.970762159065, 292.947759167446}};

    splines->AddSpline(secondSpline, 0, static_cast<UInt>(secondSpline.size()));

    std::vector<Point> thirdSpline{{647.202799419633, 482.466916504007},
                                   {486.840435519465, 144.248112641836}};

    splines->AddSpline(thirdSpline, 0, static_cast<UInt>(thirdSpline.size()));

    std::vector<Point> fourthSpline{{1224.50730946023, 747.793736775192},
                                    {1568.55747200968, 464.97284044217}};

    splines->AddSpline(fourthSpline, 0, static_cast<UInt>(fourthSpline.size()));

    CurvilinearParameters curvilinearParameters;
    SplinesToCurvilinearParameters splinesToCurvilinearParameters;

    splinesToCurvilinearParameters.aspect_ratio = 0.1;
    splinesToCurvilinearParameters.aspect_ratio_grow_factor = 1.1;
    splinesToCurvilinearParameters.average_width = 500.0;
    splinesToCurvilinearParameters.nodes_on_top_of_each_other_tolerance = 1e-4;
    splinesToCurvilinearParameters.min_cosine_crossing_angles = 0.95;
    splinesToCurvilinearParameters.check_front_collisions = false;
    splinesToCurvilinearParameters.curvature_adapted_grid_spacing = true;
    splinesToCurvilinearParameters.remove_skinny_triangles = 0;
    curvilinearParameters.m_refinement = 10;
    curvilinearParameters.n_refinement = 20;

    // create the algorithm
    CurvilinearGridFromSplines curvilinearGridFromSplines(splines, curvilinearParameters, splinesToCurvilinearParameters);

    // compute
    const auto curviGrid = curvilinearGridFromSplines.Compute();

    const double tolerance = 1e-6;
    ASSERT_NEAR(100.52954683251164, curviGrid.m_gridNodes(0, 0).x, tolerance);
    ASSERT_NEAR(53.257591250608634, curviGrid.m_gridNodes(0, 1).x, tolerance);
    ASSERT_NEAR(5.1909985476910663, curviGrid.m_gridNodes(0, 2).x, tolerance);
    ASSERT_NEAR(-43.683589051402279, curviGrid.m_gridNodes(0, 3).x, tolerance);
    ASSERT_NEAR(-93.379753864774997, curviGrid.m_gridNodes(0, 4).x, tolerance);
    ASSERT_NEAR(-143.07591867814773, curviGrid.m_gridNodes(0, 5).x, tolerance);
    ASSERT_NEAR(-199.56653306317651, curviGrid.m_gridNodes(0, 6).x, tolerance);
    ASSERT_NEAR(-263.78053278694000, curviGrid.m_gridNodes(0, 7).x, tolerance);

    ASSERT_NEAR(210.24345077714435, curviGrid.m_gridNodes(0, 0).y, tolerance);
    ASSERT_NEAR(215.47871920476678, curviGrid.m_gridNodes(0, 1).y, tolerance);
    ASSERT_NEAR(220.80199200157290, curviGrid.m_gridNodes(0, 2).y, tolerance);
    ASSERT_NEAR(226.21474851271995, curviGrid.m_gridNodes(0, 3).y, tolerance);
    ASSERT_NEAR(231.71849295101799, curviGrid.m_gridNodes(0, 4).y, tolerance);
    ASSERT_NEAR(237.22223738931604, curviGrid.m_gridNodes(0, 5).y, tolerance);
    ASSERT_NEAR(243.47845265069932, curviGrid.m_gridNodes(0, 6).y, tolerance);
    ASSERT_NEAR(250.59001639344984, curviGrid.m_gridNodes(0, 7).y, tolerance);

    ASSERT_NEAR(125.95663927357324, curviGrid.m_gridNodes(1, 0).x, tolerance);
    ASSERT_NEAR(87.191079979748508, curviGrid.m_gridNodes(1, 1).x, tolerance);
    ASSERT_NEAR(46.912312421706382, curviGrid.m_gridNodes(1, 2).x, tolerance);
    ASSERT_NEAR(-0.21122329512429161, curviGrid.m_gridNodes(1, 3).x, tolerance);
    ASSERT_NEAR(-47.914454551477142, curviGrid.m_gridNodes(1, 4).x, tolerance);
    ASSERT_NEAR(-95.617685807829986, curviGrid.m_gridNodes(1, 5).x, tolerance);
    ASSERT_NEAR(-149.56437809547685, curviGrid.m_gridNodes(1, 6).x, tolerance);
    ASSERT_NEAR(-210.56103370149361, curviGrid.m_gridNodes(1, 7).x, tolerance);

    ASSERT_NEAR(439.83786433363935, curviGrid.m_gridNodes(1, 0).y, tolerance);
    ASSERT_NEAR(521.88179941347357, curviGrid.m_gridNodes(1, 1).y, tolerance);
    ASSERT_NEAR(597.52538709342150, curviGrid.m_gridNodes(1, 2).y, tolerance);
    ASSERT_NEAR(618.74929979346371, curviGrid.m_gridNodes(1, 3).y, tolerance);
    ASSERT_NEAR(642.24827559144808, curviGrid.m_gridNodes(1, 4).y, tolerance);
    ASSERT_NEAR(665.74725138943245, curviGrid.m_gridNodes(1, 5).y, tolerance);
    ASSERT_NEAR(694.97385948812678, curviGrid.m_gridNodes(1, 6).y, tolerance);
    ASSERT_NEAR(731.13649297609936, curviGrid.m_gridNodes(1, 7).y, tolerance);

    ASSERT_NEAR(125.95663927357324, curviGrid.m_gridNodes(2, 0).x, tolerance);
    ASSERT_NEAR(87.191079979748508, curviGrid.m_gridNodes(2, 1).x, tolerance);
    ASSERT_NEAR(53.136666853719120, curviGrid.m_gridNodes(2, 2).x, tolerance);
    ASSERT_NEAR(48.869442121499411, curviGrid.m_gridNodes(2, 3).x, tolerance);
    ASSERT_NEAR(44.963973598369989, curviGrid.m_gridNodes(2, 4).x, tolerance);
    ASSERT_NEAR(41.058505075240568, curviGrid.m_gridNodes(2, 5).x, tolerance);
    ASSERT_NEAR(37.375052143743041, curviGrid.m_gridNodes(2, 6).x, tolerance);
    ASSERT_NEAR(34.121747416893278, curviGrid.m_gridNodes(2, 7).x, tolerance);

    ASSERT_NEAR(439.83786433363935, curviGrid.m_gridNodes(2, 0).y, tolerance);
    ASSERT_NEAR(521.88179941347357, curviGrid.m_gridNodes(2, 1).y, tolerance);
    ASSERT_NEAR(603.50877899861121, curviGrid.m_gridNodes(2, 2).y, tolerance);
    ASSERT_NEAR(665.92991255625338, curviGrid.m_gridNodes(2, 3).y, tolerance);
    ASSERT_NEAR(731.53111467598717, curviGrid.m_gridNodes(2, 4).y, tolerance);
    ASSERT_NEAR(797.13231679572095, curviGrid.m_gridNodes(2, 5).y, tolerance);
    ASSERT_NEAR(874.67633078399956, curviGrid.m_gridNodes(2, 6).y, tolerance);
    ASSERT_NEAR(966.34690521633274, curviGrid.m_gridNodes(2, 7).y, tolerance);
}
TEST(CurvilinearGridFromSplines, OrthogonalCurvilinearGridFromSplineWithSevenSplies)
{
    auto splines = std::make_shared<Splines>(Projection::cartesian);

    std::vector<Point> firstCentralSpline{{1.542516E+02, 7.687640E+01},
                                          {3.192526E+02, 2.733784E+02},
                                          {6.350046E+02, 5.231309E+02}};
    splines->AddSpline(firstCentralSpline, 0, static_cast<UInt>(firstCentralSpline.size()));

    std::vector<Point> secondCentralSpline{{1.310014E+02, 9.637659E+01},
                                           {2.960025E+02, 2.891285E+02},
                                           {6.222545E+02, 5.418811E+02}};
    splines->AddSpline(secondCentralSpline, 0, static_cast<UInt>(secondCentralSpline.size()));

    std::vector<Point> thirdCentralSpline{{1.782517E+02, 5.662619E+01},
                                          {3.335027E+02, 2.448781E+02},
                                          {6.455046E+02, 4.983806E+02}};

    splines->AddSpline(thirdCentralSpline, 0, static_cast<UInt>(thirdCentralSpline.size()));

    std::vector<Point> fourthBankSpline{{3.500084E+01, 1.901275E+02},
                                        {2.195020E+02, 3.813795E+02},
                                        {5.727542E+02, 6.221319E+02}};

    splines->AddSpline(fourthBankSpline, 0, static_cast<UInt>(fourthBankSpline.size()));

    std::vector<Point> fifthBankSpline{{3.177526E+02, -3.712475E+01},
                                       {4.377534E+02, 1.218768E+02},
                                       {7.445052E+02, 4.136298E+02}};

    splines->AddSpline(fifthBankSpline, 0, static_cast<UInt>(fifthBankSpline.size()));

    std::vector<Point> sixthBankSpline{{1.250633E+00, 2.748784E+02},
                                       {3.620029E+02, -3.337471E+01}};
    splines->AddSpline(sixthBankSpline, 0, static_cast<UInt>(sixthBankSpline.size()));

    std::vector<Point> seventhBankSpline{{5.030038E+02, 6.476321E+02},
                                         {7.542553E+02, 3.378790E+02}};

    splines->AddSpline(seventhBankSpline, 0, static_cast<UInt>(seventhBankSpline.size()));

    CurvilinearParameters curvilinearParameters;
    SplinesToCurvilinearParameters splinesToCurvilinearParameters;

    splinesToCurvilinearParameters.aspect_ratio = 0.1;
    splinesToCurvilinearParameters.aspect_ratio_grow_factor = 1.1;
    splinesToCurvilinearParameters.average_width = 50.0;
    splinesToCurvilinearParameters.nodes_on_top_of_each_other_tolerance = 1e-4;
    splinesToCurvilinearParameters.min_cosine_crossing_angles = 0.95;
    splinesToCurvilinearParameters.check_front_collisions = false;
    splinesToCurvilinearParameters.curvature_adapted_grid_spacing = true;
    splinesToCurvilinearParameters.remove_skinny_triangles = 0;
    splinesToCurvilinearParameters.grow_grid_outside = 1;
    curvilinearParameters.m_refinement = 20;
    curvilinearParameters.n_refinement = 40;

    // create the algorithm
    CurvilinearGridFromSplines curvilinearGridFromSplines(splines, curvilinearParameters, splinesToCurvilinearParameters);

    // compute
    const auto curviGrid = curvilinearGridFromSplines.Compute();

    const double tolerance = 1e-6;
    ASSERT_NEAR(318.51231173935355, curviGrid.m_gridNodes(0, 0).x, tolerance);
    ASSERT_NEAR(299.48353221639496, curviGrid.m_gridNodes(0, 1).x, tolerance);
    ASSERT_NEAR(282.44733906260205, curviGrid.m_gridNodes(0, 2).x, tolerance);
    ASSERT_NEAR(267.19252754998109, curviGrid.m_gridNodes(0, 3).x, tolerance);
    ASSERT_NEAR(253.53063723900343, curviGrid.m_gridNodes(0, 4).x, tolerance);
    ASSERT_NEAR(241.29345538457824, curviGrid.m_gridNodes(0, 5).x, tolerance);
    ASSERT_NEAR(230.33080196895105, curviGrid.m_gridNodes(0, 6).x, tolerance);
    ASSERT_NEAR(220.50856220652270, curviGrid.m_gridNodes(0, 7).x, tolerance);
    ASSERT_NEAR(211.70693763972625, curviGrid.m_gridNodes(0, 8).x, tolerance);
    ASSERT_NEAR(203.81889090380005, curviGrid.m_gridNodes(0, 9).x, tolerance);
    ASSERT_NEAR(196.74876236149268, curviGrid.m_gridNodes(0, 10).x, tolerance);
    ASSERT_NEAR(190.41103939285921, curviGrid.m_gridNodes(0, 11).x, tolerance);
    ASSERT_NEAR(184.72926128776140, curviGrid.m_gridNodes(0, 12).x, tolerance);
    ASSERT_NEAR(179.63504461179019, curviGrid.m_gridNodes(0, 13).x, tolerance);
    ASSERT_NEAR(174.54663336028344, curviGrid.m_gridNodes(0, 14).x, tolerance);
    ASSERT_NEAR(169.46405577607990, curviGrid.m_gridNodes(0, 15).x, tolerance);
    ASSERT_NEAR(164.38733961009973, curviGrid.m_gridNodes(0, 16).x, tolerance);
    ASSERT_NEAR(159.31651210847537, curviGrid.m_gridNodes(0, 17).x, tolerance);
    ASSERT_NEAR(154.25160000000000, curviGrid.m_gridNodes(0, 18).x, tolerance);
    ASSERT_NEAR(149.31711800582937, curviGrid.m_gridNodes(0, 19).x, tolerance);
    ASSERT_NEAR(144.38562161142110, curviGrid.m_gridNodes(0, 20).x, tolerance);
    ASSERT_NEAR(139.45710784084523, curviGrid.m_gridNodes(0, 21).x, tolerance);
    ASSERT_NEAR(134.53157358568254, curviGrid.m_gridNodes(0, 22).x, tolerance);
    ASSERT_NEAR(129.60901561173370, curviGrid.m_gridNodes(0, 23).x, tolerance);
    ASSERT_NEAR(124.68943056574528, curviGrid.m_gridNodes(0, 24).x, tolerance);
    ASSERT_NEAR(119.26411278039352, curviGrid.m_gridNodes(0, 25).x, tolerance);
    ASSERT_NEAR(113.33983596520923, curviGrid.m_gridNodes(0, 26).x, tolerance);
    ASSERT_NEAR(106.94724906818072, curviGrid.m_gridNodes(0, 27).x, tolerance);
    ASSERT_NEAR(100.15448186391468, curviGrid.m_gridNodes(0, 28).x, tolerance);
    ASSERT_NEAR(93.085794486715571, curviGrid.m_gridNodes(0, 29).x, tolerance);
    ASSERT_NEAR(85.945168726353714, curviGrid.m_gridNodes(0, 30).x, tolerance);
    ASSERT_NEAR(79.042711413094352, curviGrid.m_gridNodes(0, 31).x, tolerance);
    ASSERT_NEAR(72.818165313086013, curviGrid.m_gridNodes(0, 32).x, tolerance);
    ASSERT_NEAR(67.850573321365587, curviGrid.m_gridNodes(0, 33).x, tolerance);
    ASSERT_NEAR(64.836764457842776, curviGrid.m_gridNodes(0, 34).x, tolerance);
    ASSERT_NEAR(64.515162424795946, curviGrid.m_gridNodes(0, 35).x, tolerance);
    ASSERT_NEAR(67.506775188377858, curviGrid.m_gridNodes(0, 36).x, tolerance);
    ASSERT_NEAR(74.607619951405468, curviGrid.m_gridNodes(0, 37).x, tolerance);
    ASSERT_NEAR(84.055547192477718, curviGrid.m_gridNodes(0, 38).x, tolerance);
    ASSERT_NEAR(96.661189343119005, curviGrid.m_gridNodes(0, 39).x, tolerance);
    ASSERT_NEAR(111.19388303065054, curviGrid.m_gridNodes(0, 40).x, tolerance);
    ASSERT_NEAR(130.07126718065518, curviGrid.m_gridNodes(0, 41).x, tolerance);
    ASSERT_NEAR(151.46031591233805, curviGrid.m_gridNodes(0, 42).x, tolerance);
    ASSERT_NEAR(176.73644628514995, curviGrid.m_gridNodes(0, 43).x, tolerance);
    ASSERT_NEAR(207.42707538337629, curviGrid.m_gridNodes(0, 44).x, tolerance);
    ASSERT_NEAR(241.61953469189427, curviGrid.m_gridNodes(0, 45).x, tolerance);

    ASSERT_NEAR(-32.907270623515025, curviGrid.m_gridNodes(0, 0).y, tolerance);
    ASSERT_NEAR(-20.681837507733704, curviGrid.m_gridNodes(0, 1).y, tolerance);
    ASSERT_NEAR(-9.6398276503176721, curviGrid.m_gridNodes(0, 2).y, tolerance);
    ASSERT_NEAR(0.32833495830444726, curviGrid.m_gridNodes(0, 3).y, tolerance);
    ASSERT_NEAR(9.3229308506649744, curviGrid.m_gridNodes(0, 4).y, tolerance);
    ASSERT_NEAR(17.435580744966774, curviGrid.m_gridNodes(0, 5).y, tolerance);
    ASSERT_NEAR(24.749870809921365, curviGrid.m_gridNodes(0, 6).y, tolerance);
    ASSERT_NEAR(31.341956888789504, curviGrid.m_gridNodes(0, 7).y, tolerance);
    ASSERT_NEAR(37.281141739253201, curviGrid.m_gridNodes(0, 8).y, tolerance);
    ASSERT_NEAR(42.630422109756289, curviGrid.m_gridNodes(0, 9).y, tolerance);
    ASSERT_NEAR(47.447004133801570, curviGrid.m_gridNodes(0, 10).y, tolerance);
    ASSERT_NEAR(51.782786571790027, curviGrid.m_gridNodes(0, 11).y, tolerance);
    ASSERT_NEAR(55.684812056743255, curviGrid.m_gridNodes(0, 12).y, tolerance);
    ASSERT_NEAR(59.195686975419036, curviGrid.m_gridNodes(0, 13).y, tolerance);
    ASSERT_NEAR(62.714970591706496, curviGrid.m_gridNodes(0, 14).y, tolerance);
    ASSERT_NEAR(66.242673999999440, curviGrid.m_gridNodes(0, 15).y, tolerance);
    ASSERT_NEAR(69.778807369430695, curviGrid.m_gridNodes(0, 16).y, tolerance);
    ASSERT_NEAR(73.323379936765136, curviGrid.m_gridNodes(0, 17).y, tolerance);
    ASSERT_NEAR(76.876400000000004, curviGrid.m_gridNodes(0, 18).y, tolerance);
    ASSERT_NEAR(80.337923744434121, curviGrid.m_gridNodes(0, 19).y, tolerance);
    ASSERT_NEAR(83.803699631016158, curviGrid.m_gridNodes(0, 20).y, tolerance);
    ASSERT_NEAR(87.273715653647628, curviGrid.m_gridNodes(0, 21).y, tolerance);
    ASSERT_NEAR(90.747959670019270, curviGrid.m_gridNodes(0, 22).y, tolerance);
    ASSERT_NEAR(94.226419412087282, curviGrid.m_gridNodes(0, 23).y, tolerance);
    ASSERT_NEAR(97.709082496475787, curviGrid.m_gridNodes(0, 24).y, tolerance);
    ASSERT_NEAR(101.55670701883272, curviGrid.m_gridNodes(0, 25).y, tolerance);
    ASSERT_NEAR(105.88904155040480, curviGrid.m_gridNodes(0, 26).y, tolerance);
    ASSERT_NEAR(110.86139387022313, curviGrid.m_gridNodes(0, 27).y, tolerance);
    ASSERT_NEAR(116.66849852789267, curviGrid.m_gridNodes(0, 28).y, tolerance);
    ASSERT_NEAR(123.54448633755437, curviGrid.m_gridNodes(0, 29).y, tolerance);
    ASSERT_NEAR(131.75553064619231, curviGrid.m_gridNodes(0, 30).y, tolerance);
    ASSERT_NEAR(141.58081563533725, curviGrid.m_gridNodes(0, 31).y, tolerance);
    ASSERT_NEAR(153.27757201630283, curviGrid.m_gridNodes(0, 32).y, tolerance);
    ASSERT_NEAR(167.02858202611705, curviGrid.m_gridNodes(0, 33).y, tolerance);
    ASSERT_NEAR(182.87812073196093, curviGrid.m_gridNodes(0, 34).y, tolerance);
    ASSERT_NEAR(200.67802082725191, curviGrid.m_gridNodes(0, 35).y, tolerance);
    ASSERT_NEAR(220.09368423101697, curviGrid.m_gridNodes(0, 36).y, tolerance);
    ASSERT_NEAR(240.56147175543990, curviGrid.m_gridNodes(0, 37).y, tolerance);
    ASSERT_NEAR(262.53677014280578, curviGrid.m_gridNodes(0, 38).y, tolerance);
    ASSERT_NEAR(285.72723326832892, curviGrid.m_gridNodes(0, 39).y, tolerance);
    ASSERT_NEAR(310.96873531044707, curviGrid.m_gridNodes(0, 40).y, tolerance);
    ASSERT_NEAR(336.98028385317122, curviGrid.m_gridNodes(0, 41).y, tolerance);
    ASSERT_NEAR(365.26582661002578, curviGrid.m_gridNodes(0, 42).y, tolerance);
    ASSERT_NEAR(395.14268835984223, curviGrid.m_gridNodes(0, 43).y, tolerance);
    ASSERT_NEAR(425.52204184241054, curviGrid.m_gridNodes(0, 44).y, tolerance);
    ASSERT_NEAR(458.71059422143117, curviGrid.m_gridNodes(0, 45).y, tolerance);
}

TEST(CurvilinearGridFromSplines, OrthogonalCurvilinearMeshTwoCrossingCurvatureAdaptedSpherical)
{
    std::vector<Point> firstSpline{{4.109727E+01, 4.110174E+01},
                                   {4.109865E+01, 4.110418E+01},
                                   {4.110644E+01, 4.110904E+01}};

    auto splines = std::make_shared<Splines>(Projection::spherical);
    splines->AddSpline(firstSpline, 0, static_cast<UInt>(firstSpline.size()));

    std::vector<Point> secondSpline{{4.109612E+01, 4.110473E+01},
                                    {4.109923E+01, 4.110212E+01}};
    splines->AddSpline(secondSpline, 0, static_cast<UInt>(secondSpline.size()));

    CurvilinearParameters curvilinearParameters;
    SplinesToCurvilinearParameters splinesToCurvilinearParameters;

    splinesToCurvilinearParameters.aspect_ratio = 0.1;
    splinesToCurvilinearParameters.aspect_ratio_grow_factor = 1.1;
    splinesToCurvilinearParameters.average_width = 100.0;
    splinesToCurvilinearParameters.nodes_on_top_of_each_other_tolerance = 1e-4;
    splinesToCurvilinearParameters.min_cosine_crossing_angles = 0.95;
    splinesToCurvilinearParameters.check_front_collisions = false;
    splinesToCurvilinearParameters.curvature_adapted_grid_spacing = true;
    splinesToCurvilinearParameters.remove_skinny_triangles = 0;
    curvilinearParameters.m_refinement = 10;
    curvilinearParameters.n_refinement = 20;

    // create the algorithm
    CurvilinearGridFromSplines curvilinearGridFromSplines(splines, curvilinearParameters, splinesToCurvilinearParameters);

    // compute
    const auto curviGrid = curvilinearGridFromSplines.Compute();

    // Mesh2D mesh(edges, nodes, Projection::cartesian);

    const double tolerance = 1e-6;

    ASSERT_NEAR(41.099466896304122, curviGrid.m_gridNodes(0, 0).x, tolerance);
    ASSERT_NEAR(41.099108497593477, curviGrid.m_gridNodes(0, 1).x, tolerance);
    ASSERT_NEAR(41.098791358547970, curviGrid.m_gridNodes(0, 2).x, tolerance);
    ASSERT_NEAR(41.098510729263410, curviGrid.m_gridNodes(0, 3).x, tolerance);
    ASSERT_NEAR(41.098262406654840, curviGrid.m_gridNodes(0, 4).x, tolerance);
    ASSERT_NEAR(41.098042671505588, curviGrid.m_gridNodes(0, 5).x, tolerance);
    ASSERT_NEAR(41.097848232763340, curviGrid.m_gridNodes(0, 6).x, tolerance);
    ASSERT_NEAR(41.097676178248982, curviGrid.m_gridNodes(0, 7).x, tolerance);
    ASSERT_NEAR(41.097523931039888, curviGrid.m_gridNodes(0, 8).x, tolerance);
    ASSERT_NEAR(41.097389210874489, curviGrid.m_gridNodes(0, 9).x, tolerance);
    ASSERT_NEAR(41.097270000000002, curviGrid.m_gridNodes(0, 10).x, tolerance);
    ASSERT_NEAR(41.097150789125514, curviGrid.m_gridNodes(0, 11).x, tolerance);
    ASSERT_NEAR(41.097017310879082, curviGrid.m_gridNodes(0, 12).x, tolerance);
    ASSERT_NEAR(41.096867857715743, curviGrid.m_gridNodes(0, 13).x, tolerance);
    ASSERT_NEAR(41.096700517728394, curviGrid.m_gridNodes(0, 14).x, tolerance);
    ASSERT_NEAR(41.096513150189345, curviGrid.m_gridNodes(0, 15).x, tolerance);
    ASSERT_NEAR(41.096303358164640, curviGrid.m_gridNodes(0, 16).x, tolerance);
    ASSERT_NEAR(41.096068457850762, curviGrid.m_gridNodes(0, 17).x, tolerance);
    ASSERT_NEAR(41.095805444241556, curviGrid.m_gridNodes(0, 18).x, tolerance);
    ASSERT_NEAR(41.095510952686034, curviGrid.m_gridNodes(0, 19).x, tolerance);
    ASSERT_NEAR(41.095181215845365, curviGrid.m_gridNodes(0, 20).x, tolerance);
    ASSERT_NEAR(41.094812015498384, curviGrid.m_gridNodes(0, 21).x, tolerance);

    ASSERT_NEAR(41.101732337816657, curviGrid.m_gridNodes(0, 0).y, tolerance);
    ASSERT_NEAR(41.101733590647953, curviGrid.m_gridNodes(0, 1).y, tolerance);
    ASSERT_NEAR(41.101734698311994, curviGrid.m_gridNodes(0, 2).y, tolerance);
    ASSERT_NEAR(41.101735677734105, curviGrid.m_gridNodes(0, 3).y, tolerance);
    ASSERT_NEAR(41.101736543840993, curviGrid.m_gridNodes(0, 4).y, tolerance);
    ASSERT_NEAR(41.101737309803312, curviGrid.m_gridNodes(0, 5).y, tolerance);
    ASSERT_NEAR(41.101737987246892, curviGrid.m_gridNodes(0, 6).y, tolerance);
    ASSERT_NEAR(41.101738586437314, curviGrid.m_gridNodes(0, 7).y, tolerance);
    ASSERT_NEAR(41.101739116441543, curviGrid.m_gridNodes(0, 8).y, tolerance);
    ASSERT_NEAR(41.101739585269776, curviGrid.m_gridNodes(0, 9).y, tolerance);
    ASSERT_NEAR(41.101739999999999, curviGrid.m_gridNodes(0, 10).y, tolerance);
    ASSERT_NEAR(41.101740414730223, curviGrid.m_gridNodes(0, 11).y, tolerance);
    ASSERT_NEAR(41.101740878961273, curviGrid.m_gridNodes(0, 12).y, tolerance);
    ASSERT_NEAR(41.101741398583627, curviGrid.m_gridNodes(0, 13).y, tolerance);
    ASSERT_NEAR(41.101741980184251, curviGrid.m_gridNodes(0, 14).y, tolerance);
    ASSERT_NEAR(41.101742631127856, curviGrid.m_gridNodes(0, 15).y, tolerance);
    ASSERT_NEAR(41.101743359647344, curviGrid.m_gridNodes(0, 16).y, tolerance);
    ASSERT_NEAR(41.101744174944386, curviGrid.m_gridNodes(0, 17).y, tolerance);
    ASSERT_NEAR(41.101745087301268, curviGrid.m_gridNodes(0, 18).y, tolerance);
    ASSERT_NEAR(41.101746108205063, curviGrid.m_gridNodes(0, 19).y, tolerance);
    ASSERT_NEAR(41.101747250485388, curviGrid.m_gridNodes(0, 20).y, tolerance);
    ASSERT_NEAR(41.101748528467041, curviGrid.m_gridNodes(0, 21).y, tolerance);
}

TEST(CurvilinearGridFromSplines, OrthogonalCurvilinearMeshTwoCrossingHighCurvature)
{
    std::vector<Point> firstSpline{{-103.468336664918, 420.606335102062},
                                   {111.984583950396, 845.689124424167},
                                   {1465.84415268176, 1608.50892444055}};

    auto splines = std::make_shared<Splines>(Projection::cartesian);
    splines->AddSpline(firstSpline, 0, static_cast<UInt>(firstSpline.size()));

    std::vector<Point> secondSpline{{-333.478887051536, 921.388799234953},
                                    {167.303577081355, 490.482958004326}};

    splines->AddSpline(secondSpline, 0, static_cast<UInt>(secondSpline.size()));

    CurvilinearParameters curvilinearParameters;
    SplinesToCurvilinearParameters splinesToCurvilinearParameters;

    splinesToCurvilinearParameters.aspect_ratio = 0.5;
    splinesToCurvilinearParameters.aspect_ratio_grow_factor = 1.1;
    splinesToCurvilinearParameters.average_width = 120.0;
    splinesToCurvilinearParameters.nodes_on_top_of_each_other_tolerance = 1e-4;
    splinesToCurvilinearParameters.min_cosine_crossing_angles = 0.95;
    splinesToCurvilinearParameters.check_front_collisions = false;
    splinesToCurvilinearParameters.curvature_adapted_grid_spacing = true;
    splinesToCurvilinearParameters.remove_skinny_triangles = 0;
    curvilinearParameters.m_refinement = 20;
    curvilinearParameters.n_refinement = 3;

    // create the algorithm
    CurvilinearGridFromSplines curvilinearGridFromSplines(splines, curvilinearParameters, splinesToCurvilinearParameters);

    // compute
    const auto curviGrid = curvilinearGridFromSplines.Compute();

    const double tolerance = 1e-6;
    ASSERT_NEAR(189.31789918088032, curviGrid.m_gridNodes(0, 0).x, tolerance);
    ASSERT_NEAR(47.484543647183727, curviGrid.m_gridNodes(0, 1).x, tolerance);
    ASSERT_NEAR(-44.200527431510437, curviGrid.m_gridNodes(0, 2).x, tolerance);
    ASSERT_NEAR(-103.46833666491800, curviGrid.m_gridNodes(0, 3).x, tolerance);
    ASSERT_NEAR(-162.73614589832556, curviGrid.m_gridNodes(0, 4).x, tolerance);
    ASSERT_NEAR(-269.81491405451612, curviGrid.m_gridNodes(0, 5).x, tolerance);
    ASSERT_NEAR(-463.27343290566068, curviGrid.m_gridNodes(0, 6).x, tolerance);

    ASSERT_NEAR(466.77055774144100, curviGrid.m_gridNodes(0, 0).y, tolerance);
    ASSERT_NEAR(444.40739376403246, curviGrid.m_gridNodes(0, 1).y, tolerance);
    ASSERT_NEAR(429.95121544759934, curviGrid.m_gridNodes(0, 2).y, tolerance);
    ASSERT_NEAR(420.60633510206202, curviGrid.m_gridNodes(0, 3).y, tolerance);
    ASSERT_NEAR(411.26145475652470, curviGrid.m_gridNodes(0, 4).y, tolerance);
    ASSERT_NEAR(394.37811979088713, curviGrid.m_gridNodes(0, 5).y, tolerance);
    ASSERT_NEAR(363.87510755081644, curviGrid.m_gridNodes(0, 6).y, tolerance);

    ASSERT_NEAR(182.83836514269635, curviGrid.m_gridNodes(1, 0).x, tolerance);
    ASSERT_NEAR(37.276414715909084, curviGrid.m_gridNodes(1, 1).x, tolerance);
    ASSERT_NEAR(-55.730094345896688, curviGrid.m_gridNodes(1, 2).x, tolerance);
    ASSERT_NEAR(-115.85211826633454, curviGrid.m_gridNodes(1, 3).x, tolerance);
    ASSERT_NEAR(-175.97414218677238, curviGrid.m_gridNodes(1, 4).x, tolerance);
    ASSERT_NEAR(-284.59621450100786, curviGrid.m_gridNodes(1, 5).x, tolerance);
    ASSERT_NEAR(-480.84301100694080, curviGrid.m_gridNodes(1, 6).x, tolerance);

    ASSERT_NEAR(507.86555019683885, curviGrid.m_gridNodes(1, 0).y, tolerance);
    ASSERT_NEAR(509.15016368272256, curviGrid.m_gridNodes(1, 1).y, tolerance);
    ASSERT_NEAR(503.07490961581618, curviGrid.m_gridNodes(1, 2).y, tolerance);
    ASSERT_NEAR(499.14769447767401, curviGrid.m_gridNodes(1, 3).y, tolerance);
    ASSERT_NEAR(495.22047933953183, curviGrid.m_gridNodes(1, 4).y, tolerance);
    ASSERT_NEAR(488.12520511367825, curviGrid.m_gridNodes(1, 5).y, tolerance);
    ASSERT_NEAR(475.30621899521509, curviGrid.m_gridNodes(1, 6).y, tolerance);
}

TEST(CurvilinearGridFromSplines, OrthogonalCurvilinearMeshTwoCrossingHighCurvatureRemoveSkinnyTriangles)
{
    std::vector<Point> firstSpline{{-103.468336664918, 420.606335102062},
                                   {111.984583950396, 845.689124424167},
                                   {1465.84415268176, 1608.50892444055}};

    auto splines = std::make_shared<Splines>(Projection::cartesian);
    splines->AddSpline(firstSpline, 0, static_cast<UInt>(firstSpline.size()));

    std::vector<Point> secondSpline{{-333.478887051536, 921.388799234953},
                                    {167.303577081355, 490.482958004326}};

    splines->AddSpline(secondSpline, 0, static_cast<UInt>(secondSpline.size()));

    CurvilinearParameters curvilinearParameters;
    SplinesToCurvilinearParameters splinesToCurvilinearParameters;

    splinesToCurvilinearParameters.aspect_ratio = 0.5;
    splinesToCurvilinearParameters.aspect_ratio_grow_factor = 1.1;
    splinesToCurvilinearParameters.average_width = 120.0;
    splinesToCurvilinearParameters.nodes_on_top_of_each_other_tolerance = 1e-4;
    splinesToCurvilinearParameters.min_cosine_crossing_angles = 0.95;
    splinesToCurvilinearParameters.check_front_collisions = false;
    splinesToCurvilinearParameters.curvature_adapted_grid_spacing = true;
    splinesToCurvilinearParameters.remove_skinny_triangles = true;
    curvilinearParameters.m_refinement = 20;
    curvilinearParameters.n_refinement = 3;

    // create the algorithm
    CurvilinearGridFromSplines curvilinearGridFromSplines(splines, curvilinearParameters, splinesToCurvilinearParameters);

    // compute
    const auto curviGrid = curvilinearGridFromSplines.Compute();

    const double tolerance = 1e-6;
    ASSERT_NEAR(189.31789918088032, curviGrid.m_gridNodes(0, 0).x, tolerance);
    ASSERT_NEAR(47.484543647183727, curviGrid.m_gridNodes(0, 1).x, tolerance);
    ASSERT_NEAR(-44.200527431510437, curviGrid.m_gridNodes(0, 2).x, tolerance);
    ASSERT_NEAR(-103.46833666491800, curviGrid.m_gridNodes(0, 3).x, tolerance);
    ASSERT_NEAR(-162.73614589832556, curviGrid.m_gridNodes(0, 4).x, tolerance);
    ASSERT_NEAR(-269.81491405451612, curviGrid.m_gridNodes(0, 5).x, tolerance);
    ASSERT_NEAR(-463.27343290566068, curviGrid.m_gridNodes(0, 6).x, tolerance);

    ASSERT_NEAR(466.77055774144100, curviGrid.m_gridNodes(0, 0).y, tolerance);
    ASSERT_NEAR(444.40739376403246, curviGrid.m_gridNodes(0, 1).y, tolerance);
    ASSERT_NEAR(429.95121544759934, curviGrid.m_gridNodes(0, 2).y, tolerance);
    ASSERT_NEAR(420.60633510206202, curviGrid.m_gridNodes(0, 3).y, tolerance);
    ASSERT_NEAR(411.26145475652470, curviGrid.m_gridNodes(0, 4).y, tolerance);
    ASSERT_NEAR(394.37811979088713, curviGrid.m_gridNodes(0, 5).y, tolerance);
    ASSERT_NEAR(363.87510755081644, curviGrid.m_gridNodes(0, 6).y, tolerance);

    ASSERT_NEAR(182.83836514269635, curviGrid.m_gridNodes(1, 0).x, tolerance);
    ASSERT_NEAR(42.764508901802778, curviGrid.m_gridNodes(1, 1).x, tolerance);
    ASSERT_NEAR(-55.730094345896688, curviGrid.m_gridNodes(1, 2).x, tolerance);
    ASSERT_NEAR(-115.85211826633454, curviGrid.m_gridNodes(1, 3).x, tolerance);
    ASSERT_NEAR(-175.97414218677238, curviGrid.m_gridNodes(1, 4).x, tolerance);
    ASSERT_NEAR(-284.59621450100786, curviGrid.m_gridNodes(1, 5).x, tolerance);
    ASSERT_NEAR(-480.84301100694080, curviGrid.m_gridNodes(1, 6).x, tolerance);

    ASSERT_NEAR(507.86555019683885, curviGrid.m_gridNodes(1, 0).y, tolerance);
    ASSERT_NEAR(541.25047181337595, curviGrid.m_gridNodes(1, 1).y, tolerance);
    ASSERT_NEAR(503.07490961581618, curviGrid.m_gridNodes(1, 2).y, tolerance);
    ASSERT_NEAR(499.14769447767401, curviGrid.m_gridNodes(1, 3).y, tolerance);
    ASSERT_NEAR(495.22047933953183, curviGrid.m_gridNodes(1, 4).y, tolerance);
    ASSERT_NEAR(488.12520511367825, curviGrid.m_gridNodes(1, 5).y, tolerance);
    ASSERT_NEAR(475.30621899521509, curviGrid.m_gridNodes(1, 6).y, tolerance);

    ASSERT_NEAR(182.83836514269635, curviGrid.m_gridNodes(2, 0).x, tolerance);
    ASSERT_NEAR(42.764508901802778, curviGrid.m_gridNodes(2, 1).x, tolerance);
    ASSERT_NEAR(-54.662854941964767, curviGrid.m_gridNodes(2, 2).x, tolerance);
    ASSERT_NEAR(-114.47739150129740, curviGrid.m_gridNodes(2, 3).x, tolerance);
    ASSERT_NEAR(-174.29192806063003, curviGrid.m_gridNodes(2, 4).x, tolerance);
    ASSERT_NEAR(-282.35846494548809, curviGrid.m_gridNodes(2, 5).x, tolerance);
    ASSERT_NEAR(-477.60157917913409, curviGrid.m_gridNodes(2, 6).x, tolerance);

    ASSERT_NEAR(507.86555019683885, curviGrid.m_gridNodes(2, 0).y, tolerance);
    ASSERT_NEAR(541.25047181337595, curviGrid.m_gridNodes(2, 1).y, tolerance);
    ASSERT_NEAR(544.22653740080693, curviGrid.m_gridNodes(2, 2).y, tolerance);
    ASSERT_NEAR(552.15571117238005, curviGrid.m_gridNodes(2, 3).y, tolerance);
    ASSERT_NEAR(560.08488494395317, curviGrid.m_gridNodes(2, 4).y, tolerance);
    ASSERT_NEAR(574.41047198586193, curviGrid.m_gridNodes(2, 5).y, tolerance);
    ASSERT_NEAR(600.29241757154784, curviGrid.m_gridNodes(2, 6).y, tolerance);
}
TEST(CurvilinearGridFromSplines, Compute_ThreeLongitudinalSplinesTwoCrossingSplines_ShouldComputeMesh)
{
    // Setup
    std::vector<Point> firstSpline{{7.7979524E+04, 3.7127829E+05},
                                   {7.7979524E+04, 3.7025723E+05},
                                   {7.8302860E+04, 3.6898090E+05},
                                   {7.9732343E+04, 3.6809598E+05},
                                   {8.0889543E+04, 3.6698984E+05},
                                   {8.0668314E+04, 3.6578158E+05},
                                   {7.9579184E+04, 3.6419894E+05}};

    auto splines = std::make_shared<Splines>(Projection::cartesian);
    splines->AddSpline(firstSpline, 0, static_cast<UInt>(firstSpline.size()));

    std::vector<Point> secondSpline{{7.6618112E+04, 3.7136337E+05},
                                    {7.6754253E+04, 3.7005301E+05},
                                    {7.7179694E+04, 3.6874265E+05},
                                    {7.8404966E+04, 3.6780668E+05},
                                    {7.9681290E+04, 3.6721107E+05},
                                    {8.0140766E+04, 3.6636018E+05},
                                    {7.9477078E+04, 3.6544123E+05},
                                    {7.8779354E+04, 3.6452228E+05}};

    splines->AddSpline(secondSpline, 0, static_cast<UInt>(secondSpline.size()));

    std::vector<Point> thirdSpline{{7.7281800E+04, 3.7144846E+05},
                                   {7.7366889E+04, 3.6984880E+05},
                                   {7.7928471E+04, 3.6874265E+05},
                                   {7.9153742E+04, 3.6792581E+05},
                                   {8.0242872E+04, 3.6722808E+05},
                                   {8.0481119E+04, 3.6641124E+05},
                                   {7.9970590E+04, 3.6542421E+05},
                                   {7.9579184E+04, 3.6484561E+05},
                                   {7.9170760E+04, 3.6431806E+05}};

    splines->AddSpline(thirdSpline, 0, static_cast<UInt>(thirdSpline.size()));

    std::vector<Point> fourthSpline{{7.613792E+04, 3.712157E+05},
                                    {7.831719E+04, 3.710751E+05}};

    splines->AddSpline(fourthSpline, 0, static_cast<UInt>(fourthSpline.size()));

    std::vector<Point> fifthSpline{{7.857202E+04, 3.649151E+05},
                                   {8.003072E+04, 3.641506E+05}};

    splines->AddSpline(fifthSpline, 0, static_cast<UInt>(fifthSpline.size()));

    SplinesToCurvilinearParameters splinesToCurvilinearParameters;
    CurvilinearParameters curvilinearParameters;

    curvilinearParameters.m_refinement = 200;
    curvilinearParameters.n_refinement = 40;
    splinesToCurvilinearParameters.aspect_ratio = 0.5;
    splinesToCurvilinearParameters.aspect_ratio_grow_factor = 1.0;
    splinesToCurvilinearParameters.average_width = 400.0;
    splinesToCurvilinearParameters.curvature_adapted_grid_spacing = true;
    splinesToCurvilinearParameters.grow_grid_outside = 0;
    splinesToCurvilinearParameters.maximum_num_faces_in_uniform_part = 8;

    splinesToCurvilinearParameters.nodes_on_top_of_each_other_tolerance = 1e-4;
    splinesToCurvilinearParameters.min_cosine_crossing_angles = 0.95;
    splinesToCurvilinearParameters.check_front_collisions = false;
    splinesToCurvilinearParameters.remove_skinny_triangles = true;

    CurvilinearGridFromSplines curvilinearGridFromSplines(splines, curvilinearParameters, splinesToCurvilinearParameters);

    // Compute
    const auto curviGrid = curvilinearGridFromSplines.Compute();

    ASSERT_EQ(curviGrid.GetNumNodes(), 279);
    ASSERT_EQ(curviGrid.GetNumEdges(), 518);

    const double tolerance = 1e-6;
    ASSERT_NEAR(76628.277886551819, curviGrid.m_gridNodes(0, 0).x, tolerance);
    ASSERT_NEAR(76791.455599638575, curviGrid.m_gridNodes(0, 1).x, tolerance);
    ASSERT_NEAR(76954.811639895925, curviGrid.m_gridNodes(0, 2).x, tolerance);
    ASSERT_NEAR(77118.281322492985, curviGrid.m_gridNodes(0, 3).x, tolerance);
    ASSERT_NEAR(77281.800000000003, curviGrid.m_gridNodes(0, 4).x, tolerance);
    ASSERT_NEAR(77456.012786563486, curviGrid.m_gridNodes(0, 5).x, tolerance);
    ASSERT_NEAR(77630.227785337833, curviGrid.m_gridNodes(0, 6).x, tolerance);
    ASSERT_NEAR(77804.442019475464, curviGrid.m_gridNodes(0, 7).x, tolerance);
    ASSERT_NEAR(77978.652474172515, curviGrid.m_gridNodes(0, 8).x, tolerance);

    ASSERT_NEAR(371425.60232370096, curviGrid.m_gridNodes(0, 0).y, tolerance);
    ASSERT_NEAR(371436.18997859408, curviGrid.m_gridNodes(0, 1).y, tolerance);
    ASSERT_NEAR(371443.52950145339, curviGrid.m_gridNodes(0, 2).y, tolerance);
    ASSERT_NEAR(371447.61941049609, curviGrid.m_gridNodes(0, 3).y, tolerance);
    ASSERT_NEAR(371448.46000000002, curviGrid.m_gridNodes(0, 4).y, tolerance);
    ASSERT_NEAR(371449.35556399345, curviGrid.m_gridNodes(0, 5).y, tolerance);
    ASSERT_NEAR(371449.53230994282, curviGrid.m_gridNodes(0, 6).y, tolerance);
    ASSERT_NEAR(371448.98672684340, curviGrid.m_gridNodes(0, 7).y, tolerance);
    ASSERT_NEAR(371447.71609262202, curviGrid.m_gridNodes(0, 8).y, tolerance);

    ASSERT_NEAR(76662.032448408805, curviGrid.m_gridNodes(1, 0).x, tolerance);
    ASSERT_NEAR(76817.286847389085, curviGrid.m_gridNodes(1, 1).x, tolerance);
    ASSERT_NEAR(76972.714709055581, curviGrid.m_gridNodes(1, 2).x, tolerance);
    ASSERT_NEAR(77128.255737655345, curviGrid.m_gridNodes(1, 3).x, tolerance);
    ASSERT_NEAR(77283.849613552913, curviGrid.m_gridNodes(1, 4).x, tolerance);
    ASSERT_NEAR(77456.417295877734, curviGrid.m_gridNodes(1, 5).x, tolerance);
    ASSERT_NEAR(77628.979126778970, curviGrid.m_gridNodes(1, 6).x, tolerance);
    ASSERT_NEAR(77801.533934459163, curviGrid.m_gridNodes(1, 7).x, tolerance);
    ASSERT_NEAR(77974.080718107303, curviGrid.m_gridNodes(1, 8).x, tolerance);

    ASSERT_NEAR(371028.00333830336, curviGrid.m_gridNodes(1, 0).y, tolerance);
    ASSERT_NEAR(371038.07689532253, curviGrid.m_gridNodes(1, 1).y, tolerance);
    ASSERT_NEAR(371045.06020822527, curviGrid.m_gridNodes(1, 2).y, tolerance);
    ASSERT_NEAR(371048.95174731233, curviGrid.m_gridNodes(1, 3).y, tolerance);
    ASSERT_NEAR(371049.75159831985, curviGrid.m_gridNodes(1, 4).y, tolerance);
    ASSERT_NEAR(371050.63870543620, curviGrid.m_gridNodes(1, 5).y, tolerance);
    ASSERT_NEAR(371050.81377420085, curviGrid.m_gridNodes(1, 6).y, tolerance);
    ASSERT_NEAR(371050.27338789281, curviGrid.m_gridNodes(1, 7).y, tolerance);
    ASSERT_NEAR(371049.01488794532, curviGrid.m_gridNodes(1, 8).y, tolerance);

    ASSERT_NEAR(76696.427555678631, curviGrid.m_gridNodes(2, 0).x, tolerance);
    ASSERT_NEAR(76845.029240917880, curviGrid.m_gridNodes(2, 1).x, tolerance);
    ASSERT_NEAR(76993.806439139778, curviGrid.m_gridNodes(2, 2).x, tolerance);
    ASSERT_NEAR(77142.721409325459, curviGrid.m_gridNodes(2, 3).x, tolerance);
    ASSERT_NEAR(77291.736508965405, curviGrid.m_gridNodes(2, 4).x, tolerance);
    ASSERT_NEAR(77461.424214244718, curviGrid.m_gridNodes(2, 5).x, tolerance);
    ASSERT_NEAR(77631.102117621107, curviGrid.m_gridNodes(2, 6).x, tolerance);
    ASSERT_NEAR(77800.766971634963, curviGrid.m_gridNodes(2, 7).x, tolerance);
    ASSERT_NEAR(77970.415790794781, curviGrid.m_gridNodes(2, 8).x, tolerance);

    ASSERT_NEAR(370641.72934693948, curviGrid.m_gridNodes(2, 0).y, tolerance);
    ASSERT_NEAR(370652.41904621699, curviGrid.m_gridNodes(2, 1).y, tolerance);
    ASSERT_NEAR(370660.58064652071, curviGrid.m_gridNodes(2, 2).y, tolerance);
    ASSERT_NEAR(370666.20885042200, curviGrid.m_gridNodes(2, 3).y, tolerance);
    ASSERT_NEAR(370669.29796934803, curviGrid.m_gridNodes(2, 4).y, tolerance);
    ASSERT_NEAR(370672.81563637080, curviGrid.m_gridNodes(2, 5).y, tolerance);
    ASSERT_NEAR(370675.06421111501, curviGrid.m_gridNodes(2, 6).y, tolerance);
    ASSERT_NEAR(370676.02282016393, curviGrid.m_gridNodes(2, 7).y, tolerance);
    ASSERT_NEAR(370675.67515379097, curviGrid.m_gridNodes(2, 8).y, tolerance);
}
