//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2025.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation version 3.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// contact: delft3d.support@deltares.nl
// Stichting Deltares
// P.O. Box 177
// 2600 MH Delft, The Netherlands
//
// All indications and logos of, and references to, "Delft3D" and "Deltares"
// are registered trademarks of Stichting Deltares, and remain the property of
// Stichting Deltares. All rights reserved.
//
//------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridFromSplines.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridSplineToGrid.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Parameters.hpp>
#include <MeshKernel/Splines.hpp>
#include <MeshKernel/Utilities/Utilities.hpp>
#include <TestUtils/Definitions.hpp>

#include <fstream>
#include <iomanip>

using namespace meshkernel;

TEST(CurvilinearGridFromSplines, ComputeSplinesProperties)
{
    std::vector<Point> firstSpline{{152.001571655273, 86.6264953613281},
                                   {374.752960205078, 336.378997802734},
                                   {850.255920410156, 499.130676269531}};

    auto splines = std::make_shared<Splines>(Projection::cartesian);
    splines->AddSpline(firstSpline);

    std::vector<Point> secondSpline{{72.5010681152344, 391.129577636719},
                                    {462.503479003906, 90.3765411376953}};

    splines->AddSpline(secondSpline);

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
    EXPECT_NEAR(253.52971595547601, curvilinearGridFromSplines.m_maximumGridHeights[0], tolerance);
    EXPECT_NEAR(0.0, curvilinearGridFromSplines.m_maximumGridHeights[1], tolerance);

    ASSERT_NEAR(152.001571655273, curvilinearGridFromSplines.m_gridLine[0].x, tolerance);
    ASSERT_NEAR(448.3149638146893, curvilinearGridFromSplines.m_gridLine[1].x, tolerance);
    ASSERT_NEAR(850.2507530095439, curvilinearGridFromSplines.m_gridLine[2].x, tolerance);
    ASSERT_NEAR(-999.000000000000, curvilinearGridFromSplines.m_gridLine[3].x, tolerance);
    ASSERT_NEAR(850.2507530095439, curvilinearGridFromSplines.m_gridLine[4].x, tolerance);
    ASSERT_NEAR(448.3149638146893, curvilinearGridFromSplines.m_gridLine[5].x, tolerance);
    ASSERT_NEAR(152.001571655273, curvilinearGridFromSplines.m_gridLine[6].x, tolerance);

    ASSERT_NEAR(86.6264953613281, curvilinearGridFromSplines.m_gridLine[0].y, tolerance);
    ASSERT_NEAR(373.7229667546401, curvilinearGridFromSplines.m_gridLine[1].y, tolerance);
    ASSERT_NEAR(499.1293237106543, curvilinearGridFromSplines.m_gridLine[2].y, tolerance);
    ASSERT_NEAR(-999.000000000000, curvilinearGridFromSplines.m_gridLine[3].y, tolerance);
    ASSERT_NEAR(499.1293237106543, curvilinearGridFromSplines.m_gridLine[4].y, tolerance);
    ASSERT_NEAR(373.7229667546401, curvilinearGridFromSplines.m_gridLine[5].y, tolerance);
    ASSERT_NEAR(86.6264953613281, curvilinearGridFromSplines.m_gridLine[6].y, tolerance);

    ASSERT_NEAR(4.214936859441928e-14, curvilinearGridFromSplines.m_gridLineDimensionalCoordinates[0], tolerance);
    ASSERT_NEAR(1.191972880809887, curvilinearGridFromSplines.m_gridLineDimensionalCoordinates[1], tolerance);
    ASSERT_NEAR(1.999990407484028, curvilinearGridFromSplines.m_gridLineDimensionalCoordinates[2], tolerance);
    ASSERT_NEAR(-999.000000000000, curvilinearGridFromSplines.m_gridLineDimensionalCoordinates[3], tolerance);
    ASSERT_NEAR(1.999990407484028, curvilinearGridFromSplines.m_gridLineDimensionalCoordinates[4], tolerance);
    ASSERT_NEAR(1.191972880809887, curvilinearGridFromSplines.m_gridLineDimensionalCoordinates[5], tolerance);
    ASSERT_NEAR(4.214936859441928e-14, curvilinearGridFromSplines.m_gridLineDimensionalCoordinates[6], tolerance);
}

TEST(CurvilinearGridFromSplines, ComputeBoundingBox)
{
    std::vector<Point> firstSpline{{91.2511901855469, 299.628631591797},
                                   {354.502838134766, 518.630859375000},
                                   {770.755432128906, 607.881774902344}};

    const auto splines = std::make_shared<Splines>(Projection::cartesian);
    splines->AddSpline(firstSpline);

    std::vector<Point> secondSpline;
    secondSpline.emplace_back(Point{273.502319335938, 86.6264953613281});
    secondSpline.emplace_back(Point{557.004089355469, 316.128814697266});
    secondSpline.emplace_back(Point{847.255920410156, 409.129730224609});
    splines->AddSpline(secondSpline);

    std::vector<Point> thirdSpline;
    thirdSpline.emplace_back(Point{62.7510070800781, 396.379608154297});
    thirdSpline.emplace_back(Point{350.752807617188, 73.8763732910156});
    splines->AddSpline(thirdSpline);

    std::vector<Point> fourthSpline;
    fourthSpline.emplace_back(Point{704.755004882812, 636.382019042969});
    fourthSpline.emplace_back(Point{845.005859375000, 285.378509521484});
    splines->AddSpline(fourthSpline);

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
    splines->AddSpline(firstSpline);

    std::vector<Point> secondSpline;
    secondSpline.emplace_back(Point{-179.920786312031, 1068.50251010579});
    secondSpline.emplace_back(Point{600.169022647819, 321.964950993679});
    splines->AddSpline(secondSpline);

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

    ASSERT_NEAR(548.6410521037716, curviGrid->GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(438.2701179016752, curviGrid->GetNode(0, 1).x, tolerance);
    ASSERT_NEAR(342.7746239001754, curviGrid->GetNode(0, 2).x, tolerance);
    ASSERT_NEAR(260.1497060219567, curviGrid->GetNode(0, 3).x, tolerance);
    ASSERT_NEAR(188.6607093340589, curviGrid->GetNode(0, 4).x, tolerance);
    ASSERT_NEAR(126.8067701268222, curviGrid->GetNode(0, 5).x, tolerance);
    ASSERT_NEAR(73.28930628332408, curviGrid->GetNode(0, 6).x, tolerance);
    ASSERT_NEAR(26.9847544156713, curviGrid->GetNode(0, 7).x, tolerance);
    ASSERT_NEAR(-19.31979745198148, curviGrid->GetNode(0, 8).x, tolerance);
    ASSERT_NEAR(-71.53460511702619, curviGrid->GetNode(0, 9).x, tolerance);
    ASSERT_NEAR(-130.4140463281223, curviGrid->GetNode(0, 10).x, tolerance);
    ASSERT_NEAR(-196.8087866778367, curviGrid->GetNode(0, 11).x, tolerance);
    ASSERT_NEAR(-271.6780696640897, curviGrid->GetNode(0, 12).x, tolerance);
    ASSERT_NEAR(-356.1035754398605, curviGrid->GetNode(0, 13).x, tolerance);
    ASSERT_NEAR(-451.305048476586, curviGrid->GetNode(0, 14).x, tolerance);

    ASSERT_NEAR(115.0282209714115, curviGrid->GetNode(0, 0).y, tolerance);
    ASSERT_NEAR(159.9946081832376, curviGrid->GetNode(0, 1).y, tolerance);
    ASSERT_NEAR(198.9005700770937, curviGrid->GetNode(0, 2).y, tolerance);
    ASSERT_NEAR(232.5629113222469, curviGrid->GetNode(0, 3).y, tolerance);
    ASSERT_NEAR(261.6883502767349, curviGrid->GetNode(0, 4).y, tolerance);
    ASSERT_NEAR(286.888356066996, curviGrid->GetNode(0, 5).y, tolerance);
    ASSERT_NEAR(308.6919859735365, curviGrid->GetNode(0, 6).y, tolerance);
    ASSERT_NEAR(327.556992634968, curviGrid->GetNode(0, 7).y, tolerance);
    ASSERT_NEAR(346.4219992963995, curviGrid->GetNode(0, 8).y, tolerance);
    ASSERT_NEAR(367.6949121029384, curviGrid->GetNode(0, 9).y, tolerance);
    ASSERT_NEAR(391.6830732203045, curviGrid->GetNode(0, 10).y, tolerance);
    ASSERT_NEAR(418.7330535885708, curviGrid->GetNode(0, 11).y, tolerance);
    ASSERT_NEAR(449.2356600344596, curviGrid->GetNode(0, 12).y, tolerance);
    ASSERT_NEAR(483.6315814852674, curviGrid->GetNode(0, 13).y, tolerance);
    ASSERT_NEAR(522.4177558585591, curviGrid->GetNode(0, 14).y, tolerance);
}

TEST(CurvilinearGridFromSplines, OrthogonalCurvilinearMeshTwoCrossingCurvatureNotAdapted)
{
    std::vector<Point> firstSpline{{26.9847544156713, 327.556992634968},
                                   {340.139086327869, 819.656657068422},
                                   {2048.50780774173, 1644.48279915859}};

    auto splines = std::make_shared<Splines>(Projection::cartesian);
    splines->AddSpline(firstSpline);

    std::vector<Point> secondSpline;
    secondSpline.emplace_back(Point{-179.920786312031, 1068.50251010579});
    secondSpline.emplace_back(Point{600.169022647819, 321.964950993679});
    splines->AddSpline(secondSpline);

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
    ASSERT_NEAR(548.64105210377159, curviGrid->GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(438.27011790167523, curviGrid->GetNode(0, 1).x, tolerance);
    ASSERT_NEAR(342.77462390017541, curviGrid->GetNode(0, 2).x, tolerance);
    ASSERT_NEAR(260.14970602195672, curviGrid->GetNode(0, 3).x, tolerance);
    ASSERT_NEAR(188.66070933405891, curviGrid->GetNode(0, 4).x, tolerance);
    ASSERT_NEAR(126.80677012682224, curviGrid->GetNode(0, 5).x, tolerance);
    ASSERT_NEAR(73.289306283324080, curviGrid->GetNode(0, 6).x, tolerance);
    ASSERT_NEAR(26.984754415671301, curviGrid->GetNode(0, 7).x, tolerance);
    ASSERT_NEAR(-19.319797451981476, curviGrid->GetNode(0, 8).x, tolerance);
    ASSERT_NEAR(-71.534605117026189, curviGrid->GetNode(0, 9).x, tolerance);
    ASSERT_NEAR(-130.41404632812231, curviGrid->GetNode(0, 10).x, tolerance);
    ASSERT_NEAR(-196.80878667783674, curviGrid->GetNode(0, 11).x, tolerance);
    ASSERT_NEAR(-271.67806966408966, curviGrid->GetNode(0, 12).x, tolerance);
    ASSERT_NEAR(-356.10357543986055, curviGrid->GetNode(0, 13).x, tolerance);
    ASSERT_NEAR(-451.30504847658597, curviGrid->GetNode(0, 14).x, tolerance);

    ASSERT_NEAR(115.02822097141149, curviGrid->GetNode(0, 0).y, tolerance);
    ASSERT_NEAR(159.99460818323763, curviGrid->GetNode(0, 1).y, tolerance);
    ASSERT_NEAR(198.90057007709373, curviGrid->GetNode(0, 2).y, tolerance);
    ASSERT_NEAR(232.56291132224692, curviGrid->GetNode(0, 3).y, tolerance);
    ASSERT_NEAR(261.68835027673492, curviGrid->GetNode(0, 4).y, tolerance);
    ASSERT_NEAR(286.88835606699604, curviGrid->GetNode(0, 5).y, tolerance);
    ASSERT_NEAR(308.69198597353653, curviGrid->GetNode(0, 6).y, tolerance);
    ASSERT_NEAR(327.55699263496803, curviGrid->GetNode(0, 7).y, tolerance);
    ASSERT_NEAR(346.42199929639952, curviGrid->GetNode(0, 8).y, tolerance);
    ASSERT_NEAR(367.69491210293836, curviGrid->GetNode(0, 9).y, tolerance);
    ASSERT_NEAR(391.68307322030449, curviGrid->GetNode(0, 10).y, tolerance);
    ASSERT_NEAR(418.73305358857078, curviGrid->GetNode(0, 11).y, tolerance);
    ASSERT_NEAR(449.23566003445961, curviGrid->GetNode(0, 12).y, tolerance);
    ASSERT_NEAR(483.63158148526736, curviGrid->GetNode(0, 13).y, tolerance);
    ASSERT_NEAR(522.41775585855908, curviGrid->GetNode(0, 14).y, tolerance);
}

TEST(CurvilinearGridFromSplines, OrthogonalCurvilinearMeshTwoCrossingCurvatureAdaptedLargeMRefinement)
{
    std::vector<Point> firstSpline{{26.9847544156713, 327.556992634968},
                                   {340.139086327869, 819.656657068422},
                                   {2048.50780774173, 1644.48279915859}};

    auto splines = std::make_shared<Splines>(Projection::cartesian);
    splines->AddSpline(firstSpline);

    std::vector<Point> secondSpline;
    secondSpline.emplace_back(Point{-179.920786312031, 1068.50251010579});
    secondSpline.emplace_back(Point{600.169022647819, 321.964950993679});
    splines->AddSpline(secondSpline);

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

    ASSERT_NEAR(548.6410521037716, curviGrid->GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(438.2701179016752, curviGrid->GetNode(0, 1).x, tolerance);
    ASSERT_NEAR(342.7746239001754, curviGrid->GetNode(0, 2).x, tolerance);
    ASSERT_NEAR(260.1497060219567, curviGrid->GetNode(0, 3).x, tolerance);
    ASSERT_NEAR(188.6607093340589, curviGrid->GetNode(0, 4).x, tolerance);
    ASSERT_NEAR(126.8067701268222, curviGrid->GetNode(0, 5).x, tolerance);
    ASSERT_NEAR(73.28930628332408, curviGrid->GetNode(0, 6).x, tolerance);
    ASSERT_NEAR(26.9847544156713, curviGrid->GetNode(0, 7).x, tolerance);
    ASSERT_NEAR(-19.31979745198148, curviGrid->GetNode(0, 8).x, tolerance);
    ASSERT_NEAR(-71.53460511702619, curviGrid->GetNode(0, 9).x, tolerance);
    ASSERT_NEAR(-130.4140463281223, curviGrid->GetNode(0, 10).x, tolerance);
    ASSERT_NEAR(-196.8087866778367, curviGrid->GetNode(0, 11).x, tolerance);
    ASSERT_NEAR(-271.6780696640897, curviGrid->GetNode(0, 12).x, tolerance);
    ASSERT_NEAR(-356.1035754398605, curviGrid->GetNode(0, 13).x, tolerance);
    ASSERT_NEAR(-451.305048476586, curviGrid->GetNode(0, 14).x, tolerance);

    ASSERT_NEAR(115.0282209714115, curviGrid->GetNode(0, 0).y, tolerance);
    ASSERT_NEAR(159.9946081832376, curviGrid->GetNode(0, 1).y, tolerance);
    ASSERT_NEAR(198.9005700770937, curviGrid->GetNode(0, 2).y, tolerance);
    ASSERT_NEAR(232.5629113222469, curviGrid->GetNode(0, 3).y, tolerance);
    ASSERT_NEAR(261.6883502767349, curviGrid->GetNode(0, 4).y, tolerance);
    ASSERT_NEAR(286.888356066996, curviGrid->GetNode(0, 5).y, tolerance);
    ASSERT_NEAR(308.6919859735365, curviGrid->GetNode(0, 6).y, tolerance);
    ASSERT_NEAR(327.556992634968, curviGrid->GetNode(0, 7).y, tolerance);
    ASSERT_NEAR(346.4219992963995, curviGrid->GetNode(0, 8).y, tolerance);
    ASSERT_NEAR(367.6949121029384, curviGrid->GetNode(0, 9).y, tolerance);
    ASSERT_NEAR(391.6830732203045, curviGrid->GetNode(0, 10).y, tolerance);
    ASSERT_NEAR(418.7330535885708, curviGrid->GetNode(0, 11).y, tolerance);
    ASSERT_NEAR(449.2356600344596, curviGrid->GetNode(0, 12).y, tolerance);
    ASSERT_NEAR(483.6315814852674, curviGrid->GetNode(0, 13).y, tolerance);
    ASSERT_NEAR(522.4177558585591, curviGrid->GetNode(0, 14).y, tolerance);

    ASSERT_NEAR(634.3637498199938, curviGrid->GetNode(1, 0).x, tolerance);
    ASSERT_NEAR(537.7939316539596, curviGrid->GetNode(1, 1).x, tolerance);
    ASSERT_NEAR(454.2394835733226, curviGrid->GetNode(1, 2).x, tolerance);
    ASSERT_NEAR(381.9462357797242, curviGrid->GetNode(1, 3).x, tolerance);
    ASSERT_NEAR(319.3964398490182, curviGrid->GetNode(1, 4).x, tolerance);
    ASSERT_NEAR(265.2769046174108, curviGrid->GetNode(1, 5).x, tolerance);
    ASSERT_NEAR(218.4514266105773, curviGrid->GetNode(1, 6).x, tolerance);
    ASSERT_NEAR(177.9369362110388, curviGrid->GetNode(1, 7).x, tolerance);
    ASSERT_NEAR(137.4224458115003, curviGrid->GetNode(1, 8).x, tolerance);
    ASSERT_NEAR(91.73673590752415, curviGrid->GetNode(1, 9).x, tolerance);
    ASSERT_NEAR(40.21975842990438, curviGrid->GetNode(1, 10).x, tolerance);
    ASSERT_NEAR(-17.87278241071843, curviGrid->GetNode(1, 11).x, tolerance);
    ASSERT_NEAR(-83.38018339878106, curviGrid->GetNode(1, 12).x, tolerance);
    ASSERT_NEAR(-157.2488681234484, curviGrid->GetNode(1, 13).x, tolerance);
    ASSERT_NEAR(-240.5460605116284, curviGrid->GetNode(1, 14).x, tolerance);

    ASSERT_NEAR(325.4363683881346, curviGrid->GetNode(1, 0).y, tolerance);
    ASSERT_NEAR(404.2778824309436, curviGrid->GetNode(1, 1).y, tolerance);
    ASSERT_NEAR(472.4933903069406, curviGrid->GetNode(1, 2).y, tolerance);
    ASSERT_NEAR(531.5150310158385, curviGrid->GetNode(1, 3).y, tolerance);
    ASSERT_NEAR(582.5819244604007, curviGrid->GetNode(1, 4).y, tolerance);
    ASSERT_NEAR(626.7661859407009, curviGrid->GetNode(1, 5).y, tolerance);
    ASSERT_NEAR(664.9954344983574, curviGrid->GetNode(1, 6).y, tolerance);
    ASSERT_NEAR(698.0722676583888, curviGrid->GetNode(1, 7).y, tolerance);
    ASSERT_NEAR(731.1491008184202, curviGrid->GetNode(1, 8).y, tolerance);
    ASSERT_NEAR(768.4478201010555, curviGrid->GetNode(1, 9).y, tolerance);
    ASSERT_NEAR(810.5073018354163, curviGrid->GetNode(1, 10).y, tolerance);
    ASSERT_NEAR(857.9352038595642, curviGrid->GetNode(1, 11).y, tolerance);
    ASSERT_NEAR(911.4167447073323, curviGrid->GetNode(1, 12).y, tolerance);
    ASSERT_NEAR(971.7246033597175, curviGrid->GetNode(1, 13).y, tolerance);
    ASSERT_NEAR(1039.730082588311, curviGrid->GetNode(1, 14).y, tolerance);
}

TEST(CurvilinearGridFromSplines, OrthogonalCurvilinearMeshFourSplineCrossingFront)
{
    std::vector<Point> firstSpline{{-93.379753864775, 231.718492951018},
                                   {72.8139687226708, 724.468302026077},
                                   {746.335897103372, 234.634172294657},
                                   {1498.58116776234, 776.950530211586}};

    auto splines = std::make_shared<Splines>(Projection::cartesian);
    splines->AddSpline(firstSpline);

    std::vector<Point> secondSpline{{-250.826438421303, 394.996536194825},
                                    {101.970762159065, 292.947759167446}};

    splines->AddSpline(secondSpline);

    std::vector<Point> thirdSpline{{647.202799419633, 482.466916504007},
                                   {486.840435519465, 144.248112641836}};

    splines->AddSpline(thirdSpline);

    std::vector<Point> fourthSpline{{1224.50730946023, 747.793736775192},
                                    {1568.55747200968, 464.97284044217}};

    splines->AddSpline(fourthSpline);

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

    ASSERT_NEAR(98.50110148378664, curviGrid->GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(51.93985307872587, curviGrid->GetNode(0, 1).x, tolerance);
    ASSERT_NEAR(4.451658119478935, curviGrid->GetNode(0, 2).x, tolerance);
    ASSERT_NEAR(-43.98193714989393, curviGrid->GetNode(0, 3).x, tolerance);
    ASSERT_NEAR(-93.379753864775, curviGrid->GetNode(0, 4).x, tolerance);
    ASSERT_NEAR(-142.7775705796561, curviGrid->GetNode(0, 5).x, tolerance);
    ASSERT_NEAR(-199.1190763654492, curviGrid->GetNode(0, 6).x, tolerance);
    ASSERT_NEAR(-263.3803228195721, curviGrid->GetNode(0, 7).x, tolerance);

    ASSERT_NEAR(201.6662658402313, curviGrid->GetNode(0, 0).y, tolerance);
    ASSERT_NEAR(208.9586515714122, curviGrid->GetNode(0, 1).y, tolerance);
    ASSERT_NEAR(216.3962149354287, curviGrid->GetNode(0, 2).y, tolerance);
    ASSERT_NEAR(223.981846144874, curviGrid->GetNode(0, 3).y, tolerance);
    ASSERT_NEAR(231.718492951018, curviGrid->GetNode(0, 4).y, tolerance);
    ASSERT_NEAR(239.455139757162, curviGrid->GetNode(0, 5).y, tolerance);
    ASSERT_NEAR(248.279301628666, curviGrid->GetNode(0, 6).y, tolerance);
    ASSERT_NEAR(258.3438469912214, curviGrid->GetNode(0, 7).y, tolerance);

    ASSERT_NEAR(129.8414588783787, curviGrid->GetNode(1, 0).x, tolerance);
    ASSERT_NEAR(93.03884627050913, curviGrid->GetNode(1, 1).x, tolerance);
    ASSERT_NEAR(55.77349982752719, curviGrid->GetNode(1, 2).x, tolerance);
    ASSERT_NEAR(18.05202128720908, curviGrid->GetNode(1, 3).x, tolerance);
    ASSERT_NEAR(-20.11792974782675, curviGrid->GetNode(1, 4).x, tolerance);
    ASSERT_NEAR(-58.28788078286259, curviGrid->GetNode(1, 5).x, tolerance);
    ASSERT_NEAR(-101.3825221570333, curviGrid->GetNode(1, 6).x, tolerance);
    ASSERT_NEAR(-150.0086874565013, curviGrid->GetNode(1, 7).x, tolerance);

    ASSERT_NEAR(401.7717209745524, curviGrid->GetNode(1, 0).y, tolerance);
    ASSERT_NEAR(471.3721472012758, curviGrid->GetNode(1, 1).y, tolerance);
    ASSERT_NEAR(544.0817088823699, curviGrid->GetNode(1, 2).y, tolerance);
    ASSERT_NEAR(620.0632733778014, curviGrid->GetNode(1, 3).y, tolerance);
    ASSERT_NEAR(699.4888656469648, curviGrid->GetNode(1, 4).y, tolerance);
    ASSERT_NEAR(778.9144579161282, curviGrid->GetNode(1, 5).y, tolerance);
    ASSERT_NEAR(872.3186964860478, curviGrid->GetNode(1, 6).y, tolerance);
    ASSERT_NEAR(982.2118746597305, curviGrid->GetNode(1, 7).y, tolerance);

    ASSERT_NEAR(259.9214173562335, curviGrid->GetNode(2, 0).x, tolerance);
    ASSERT_NEAR(282.4167213795628, curviGrid->GetNode(2, 1).x, tolerance);
    ASSERT_NEAR(306.1914533288371, curviGrid->GetNode(2, 2).x, tolerance);
    ASSERT_NEAR(331.3376897997447, curviGrid->GetNode(2, 3).x, tolerance);
    ASSERT_NEAR(357.0573938989936, curviGrid->GetNode(2, 4).x, tolerance);
    ASSERT_NEAR(382.7770979982425, curviGrid->GetNode(2, 5).x, tolerance);
    ASSERT_NEAR(414.5553116685546, curviGrid->GetNode(2, 6).x, tolerance);
    ASSERT_NEAR(452.2132521545221, curviGrid->GetNode(2, 7).x, tolerance);

    ASSERT_NEAR(327.3853043312594, curviGrid->GetNode(2, 0).y, tolerance);
    ASSERT_NEAR(363.0761292310531, curviGrid->GetNode(2, 1).y, tolerance);
    ASSERT_NEAR(400.8798318055067, curviGrid->GetNode(2, 2).y, tolerance);
    ASSERT_NEAR(440.9104004906728, curviGrid->GetNode(2, 3).y, tolerance);
    ASSERT_NEAR(483.8005989787546, curviGrid->GetNode(2, 4).y, tolerance);
    ASSERT_NEAR(526.6907974668364, curviGrid->GetNode(2, 5).y, tolerance);
    ASSERT_NEAR(577.2788831429193, curviGrid->GetNode(2, 6).y, tolerance);
    ASSERT_NEAR(637.8303678093578, curviGrid->GetNode(2, 7).y, tolerance);
}

TEST(CurvilinearGridFromSplines, OrthogonalCurvilinearGridFromSplineWithSevenSplies)
{
    auto splines = std::make_shared<Splines>(Projection::cartesian);

    std::vector<Point> firstCentralSpline{{1.542516E+02, 7.687640E+01},
                                          {3.192526E+02, 2.733784E+02},
                                          {6.350046E+02, 5.231309E+02}};
    splines->AddSpline(firstCentralSpline);

    std::vector<Point> secondCentralSpline{{1.310014E+02, 9.637659E+01},
                                           {2.960025E+02, 2.891285E+02},
                                           {6.222545E+02, 5.418811E+02}};
    splines->AddSpline(secondCentralSpline);

    std::vector<Point> thirdCentralSpline{{1.782517E+02, 5.662619E+01},
                                          {3.335027E+02, 2.448781E+02},
                                          {6.455046E+02, 4.983806E+02}};

    splines->AddSpline(thirdCentralSpline);

    std::vector<Point> fourthBankSpline{{3.500084E+01, 1.901275E+02},
                                        {2.195020E+02, 3.813795E+02},
                                        {5.727542E+02, 6.221319E+02}};

    splines->AddSpline(fourthBankSpline);

    std::vector<Point> fifthBankSpline{{3.177526E+02, -3.712475E+01},
                                       {4.377534E+02, 1.218768E+02},
                                       {7.445052E+02, 4.136298E+02}};

    splines->AddSpline(fifthBankSpline);

    std::vector<Point> sixthBankSpline{{1.250633E+00, 2.748784E+02},
                                       {3.620029E+02, -3.337471E+01}};
    splines->AddSpline(sixthBankSpline);

    std::vector<Point> seventhBankSpline{{5.030038E+02, 6.476321E+02},
                                         {7.542553E+02, 3.378790E+02}};

    splines->AddSpline(seventhBankSpline);

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

    ASSERT_NEAR(318.468387062165, curviGrid->GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(299.4447771760493, curviGrid->GetNode(0, 1).x, tolerance);
    ASSERT_NEAR(282.412979269927, curviGrid->GetNode(0, 2).x, tolerance);
    ASSERT_NEAR(267.1619561462167, curviGrid->GetNode(0, 3).x, tolerance);
    ASSERT_NEAR(253.5033805924582, curviGrid->GetNode(0, 4).x, tolerance);
    ASSERT_NEAR(241.2691419540521, curviGrid->GetNode(0, 5).x, tolerance);
    ASSERT_NEAR(230.3091363229219, curviGrid->GetNode(0, 6).x, tolerance);
    ASSERT_NEAR(220.4893045144444, curviGrid->GetNode(0, 7).x, tolerance);
    ASSERT_NEAR(211.6898881368979, curviGrid->GetNode(0, 8).x, tolerance);
    ASSERT_NEAR(203.8038785127272, curviGrid->GetNode(0, 9).x, tolerance);
    ASSERT_NEAR(196.735636622716, curviGrid->GetNode(0, 10).x, tolerance);
    ASSERT_NEAR(190.3996648898901, curviGrid->GetNode(0, 11).x, tolerance);
    ASSERT_NEAR(184.7195138822011, curviGrid->GetNode(0, 12).x, tolerance);
    ASSERT_NEAR(179.6268089728142, curviGrid->GetNode(0, 13).x, tolerance);
    ASSERT_NEAR(174.539949309664, curviGrid->GetNode(0, 14).x, tolerance);
    ASSERT_NEAR(169.4589672147332, curviGrid->GetNode(0, 15).x, tolerance);
    ASSERT_NEAR(164.3838943799172, curviGrid->GetNode(0, 16).x, tolerance);
    ASSERT_NEAR(159.3147618489951, curviGrid->GetNode(0, 17).x, tolerance);
    ASSERT_NEAR(154.2516, curviGrid->GetNode(0, 18).x, tolerance);
    ASSERT_NEAR(149.3181855320145, curviGrid->GetNode(0, 19).x, tolerance);
    ASSERT_NEAR(144.3879716820646, curviGrid->GetNode(0, 20).x, tolerance);
    ASSERT_NEAR(139.4609542816169, curviGrid->GetNode(0, 21).x, tolerance);
    ASSERT_NEAR(134.5371289834641, curviGrid->GetNode(0, 22).x, tolerance);
    ASSERT_NEAR(129.6164912698897, curviGrid->GetNode(0, 23).x, tolerance);
    ASSERT_NEAR(124.6990364608645, curviGrid->GetNode(0, 24).x, tolerance);
    ASSERT_NEAR(119.2762775136987, curviGrid->GetNode(0, 25).x, tolerance);
    ASSERT_NEAR(113.3537691475495, curviGrid->GetNode(0, 26).x, tolerance);
    ASSERT_NEAR(106.960261242034, curviGrid->GetNode(0, 27).x, tolerance);
    ASSERT_NEAR(100.1608923364466, curviGrid->GetNode(0, 28).x, tolerance);
    ASSERT_NEAR(93.07531884850053, curviGrid->GetNode(0, 29).x, tolerance);
    ASSERT_NEAR(85.90078959857036, curviGrid->GetNode(0, 30).x, tolerance);
    ASSERT_NEAR(78.93831744370044, curviGrid->GetNode(0, 31).x, tolerance);
    ASSERT_NEAR(72.61677070646537, curviGrid->GetNode(0, 32).x, tolerance);
    ASSERT_NEAR(67.50472563521177, curviGrid->GetNode(0, 33).x, tolerance);
    ASSERT_NEAR(64.2937046022718, curviGrid->GetNode(0, 34).x, tolerance);
    ASSERT_NEAR(63.73013613390841, curviGrid->GetNode(0, 35).x, tolerance);
    ASSERT_NEAR(66.46811887085525, curviGrid->GetNode(0, 36).x, tolerance);
    ASSERT_NEAR(72.81119152384244, curviGrid->GetNode(0, 37).x, tolerance);
    ASSERT_NEAR(82.57976793462781, curviGrid->GetNode(0, 38).x, tolerance);
    ASSERT_NEAR(95.25946118871572, curviGrid->GetNode(0, 39).x, tolerance);
    ASSERT_NEAR(109.9269290027926, curviGrid->GetNode(0, 40).x, tolerance);
    ASSERT_NEAR(128.2528847992245, curviGrid->GetNode(0, 41).x, tolerance);
    ASSERT_NEAR(149.3171906217304, curviGrid->GetNode(0, 42).x, tolerance);
    ASSERT_NEAR(175.8597855185216, curviGrid->GetNode(0, 43).x, tolerance);
    ASSERT_NEAR(205.521942487327, curviGrid->GetNode(0, 44).x, tolerance);
    ASSERT_NEAR(240.4004412835495, curviGrid->GetNode(0, 45).x, tolerance);

    ASSERT_NEAR(-32.90806030860536, curviGrid->GetNode(0, 0).y, tolerance);
    ASSERT_NEAR(-20.68148651232686, curviGrid->GetNode(0, 1).y, tolerance);
    ASSERT_NEAR(-9.638970248771939, curviGrid->GetNode(0, 2).y, tolerance);
    ASSERT_NEAR(0.3292874672802313, curviGrid->GetNode(0, 3).y, tolerance);
    ASSERT_NEAR(9.323735061703081, curviGrid->GetNode(0, 4).y, tolerance);
    ASSERT_NEAR(17.43611295346604, curviGrid->GetNode(0, 5).y, tolerance);
    ASSERT_NEAR(24.75008814194286, curviGrid->GetNode(0, 6).y, tolerance);
    ASSERT_NEAR(31.34186768952031, curviGrid->GetNode(0, 7).y, tolerance);
    ASSERT_NEAR(37.280784037788, curviGrid->GetNode(0, 8).y, tolerance);
    ASSERT_NEAR(42.62984862033718, curviGrid->GetNode(0, 9).y, tolerance);
    ASSERT_NEAR(47.44627231302959, curviGrid->GetNode(0, 10).y, tolerance);
    ASSERT_NEAR(51.7819523803277, curviGrid->GetNode(0, 11).y, tolerance);
    ASSERT_NEAR(55.68392632885328, curviGrid->GetNode(0, 12).y, tolerance);
    ASSERT_NEAR(59.19479358839903, curviGrid->GetNode(0, 13).y, tolerance);
    ASSERT_NEAR(62.71412462932994, curviGrid->GetNode(0, 14).y, tolerance);
    ASSERT_NEAR(66.24193601306862, curviGrid->GetNode(0, 15).y, tolerance);
    ASSERT_NEAR(69.77824311646854, curviGrid->GetNode(0, 16).y, tolerance);
    ASSERT_NEAR(73.32306012051703, curviGrid->GetNode(0, 17).y, tolerance);
    ASSERT_NEAR(76.8764, curviGrid->GetNode(0, 18).y, tolerance);
    ASSERT_NEAR(80.3386828370755, curviGrid->GetNode(0, 19).y, tolerance);
    ASSERT_NEAR(83.805521764773, curviGrid->GetNode(0, 20).y, tolerance);
    ASSERT_NEAR(87.27690193590402, curviGrid->GetNode(0, 21).y, tolerance);
    ASSERT_NEAR(90.7528083198717, curviGrid->GetNode(0, 22).y, tolerance);
    ASSERT_NEAR(94.23322571561403, curviGrid->GetNode(0, 23).y, tolerance);
    ASSERT_NEAR(97.718138764462, curviGrid->GetNode(0, 24).y, tolerance);
    ASSERT_NEAR(101.5685808932841, curviGrid->GetNode(0, 25).y, tolerance);
    ASSERT_NEAR(105.9025306902949, curviGrid->GetNode(0, 26).y, tolerance);
    ASSERT_NEAR(110.8729007407594, curviGrid->GetNode(0, 27).y, tolerance);
    ASSERT_NEAR(116.6714972582315, curviGrid->GetNode(0, 28).y, tolerance);
    ASSERT_NEAR(123.5293404924591, curviGrid->GetNode(0, 29).y, tolerance);
    ASSERT_NEAR(131.7100785363568, curviGrid->GetNode(0, 30).y, tolerance);
    ASSERT_NEAR(141.4922747854787, curviGrid->GetNode(0, 31).y, tolerance);
    ASSERT_NEAR(153.1362868882886, curviGrid->GetNode(0, 32).y, tolerance);
    ASSERT_NEAR(166.8336877862578, curviGrid->GetNode(0, 33).y, tolerance);
    ASSERT_NEAR(182.6439238287426, curviGrid->GetNode(0, 34).y, tolerance);
    ASSERT_NEAR(200.4373073456774, curviGrid->GetNode(0, 35).y, tolerance);
    ASSERT_NEAR(219.8898825613308, curviGrid->GetNode(0, 36).y, tolerance);
    ASSERT_NEAR(240.617990216747, curviGrid->GetNode(0, 37).y, tolerance);
    ASSERT_NEAR(262.4521008631941, curviGrid->GetNode(0, 38).y, tolerance);
    ASSERT_NEAR(285.6018943537209, curviGrid->GetNode(0, 39).y, tolerance);
    ASSERT_NEAR(310.7616821554465, curviGrid->GetNode(0, 40).y, tolerance);
    ASSERT_NEAR(337.1645788035793, curviGrid->GetNode(0, 41).y, tolerance);
    ASSERT_NEAR(365.6938042847398, curviGrid->GetNode(0, 42).y, tolerance);
    ASSERT_NEAR(394.451353236255, curviGrid->GetNode(0, 43).y, tolerance);
    ASSERT_NEAR(425.8362375946482, curviGrid->GetNode(0, 44).y, tolerance);
    ASSERT_NEAR(458.3049361663383, curviGrid->GetNode(0, 45).y, tolerance);
}

TEST(CurvilinearGridFromSplines, OrthogonalCurvilinearMeshTwoCrossingCurvatureAdaptedSpherical)
{
    std::vector<Point> firstSpline{{4.109727E+01, 4.110174E+01},
                                   {4.109865E+01, 4.110418E+01},
                                   {4.110644E+01, 4.110904E+01}};

    auto splines = std::make_shared<Splines>(Projection::spherical);
    splines->AddSpline(firstSpline);

    std::vector<Point> secondSpline{{4.109612E+01, 4.110473E+01},
                                    {4.109923E+01, 4.110212E+01}};
    splines->AddSpline(secondSpline);

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

    ASSERT_NEAR(41.09946055879178, curviGrid->GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(41.09910319449384, curviGrid->GetNode(0, 1).x, tolerance);
    ASSERT_NEAR(41.09878697060671, curviGrid->GetNode(0, 2).x, tolerance);
    ASSERT_NEAR(41.09850715099284, curviGrid->GetNode(0, 3).x, tolerance);
    ASSERT_NEAR(41.09825954474053, curviGrid->GetNode(0, 4).x, tolerance);
    ASSERT_NEAR(41.09804044339841, curviGrid->GetNode(0, 5).x, tolerance);
    ASSERT_NEAR(41.09784656543511, curviGrid->GetNode(0, 6).x, tolerance);
    ASSERT_NEAR(41.09767500709243, curviGrid->GetNode(0, 7).x, tolerance);
    ASSERT_NEAR(41.09752319889625, curviGrid->GetNode(0, 8).x, tolerance);
    ASSERT_NEAR(41.0973888671736, curviGrid->GetNode(0, 9).x, tolerance);
    ASSERT_NEAR(41.09727, curviGrid->GetNode(0, 10).x, tolerance);
    ASSERT_NEAR(41.0971511328264, curviGrid->GetNode(0, 11).x, tolerance);
    ASSERT_NEAR(41.09701803938985, curviGrid->GetNode(0, 12).x, tolerance);
    ASSERT_NEAR(41.09686901705888, curviGrid->GetNode(0, 13).x, tolerance);
    ASSERT_NEAR(41.09670215942629, curviGrid->GetNode(0, 14).x, tolerance);
    ASSERT_NEAR(41.09651533192051, curviGrid->GetNode(0, 15).x, tolerance);
    ASSERT_NEAR(41.09630614449787, curviGrid->GetNode(0, 16).x, tolerance);
    ASSERT_NEAR(41.0960719210666, curviGrid->GetNode(0, 17).x, tolerance);
    ASSERT_NEAR(41.0958096652512, curviGrid->GetNode(0, 18).x, tolerance);
    ASSERT_NEAR(41.09551602205924, curviGrid->GetNode(0, 19).x, tolerance);
    ASSERT_NEAR(41.09518723496014, curviGrid->GetNode(0, 20).x, tolerance);
    ASSERT_NEAR(41.09481909782657, curviGrid->GetNode(0, 21).x, tolerance);

    ASSERT_NEAR(41.10161413055265, curviGrid->GetNode(0, 0).y, tolerance);
    ASSERT_NEAR(41.10163466753467, curviGrid->GetNode(0, 1).y, tolerance);
    ASSERT_NEAR(41.10165283932361, curviGrid->GetNode(0, 2).y, tolerance);
    ASSERT_NEAR(41.10166891841637, curviGrid->GetNode(0, 3).y, tolerance);
    ASSERT_NEAR(41.10168314589021, curviGrid->GetNode(0, 4).y, tolerance);
    ASSERT_NEAR(41.10169573503197, curviGrid->GetNode(0, 5).y, tolerance);
    ASSERT_NEAR(41.10170687454622, curviGrid->GetNode(0, 6).y, tolerance);
    ASSERT_NEAR(41.10171673139183, curviGrid->GetNode(0, 7).y, tolerance);
    ASSERT_NEAR(41.10172545329034, curviGrid->GetNode(0, 8).y, tolerance);
    ASSERT_NEAR(41.10173317094443, curviGrid->GetNode(0, 9).y, tolerance);
    ASSERT_NEAR(41.10174, curviGrid->GetNode(0, 10).y, tolerance);
    ASSERT_NEAR(41.10174682905557, curviGrid->GetNode(0, 11).y, tolerance);
    ASSERT_NEAR(41.10175447529053, curviGrid->GetNode(0, 12).y, tolerance);
    ASSERT_NEAR(41.10176303647263, curviGrid->GetNode(0, 13).y, tolerance);
    ASSERT_NEAR(41.10177262206223, curviGrid->GetNode(0, 14).y, tolerance);
    ASSERT_NEAR(41.10178335460949, curviGrid->GetNode(0, 15).y, tolerance);
    ASSERT_NEAR(41.10179537131844, curviGrid->GetNode(0, 16).y, tolerance);
    ASSERT_NEAR(41.10180882579744, curviGrid->GetNode(0, 17).y, tolerance);
    ASSERT_NEAR(41.10182389001839, curviGrid->GetNode(0, 18).y, tolerance);
    ASSERT_NEAR(41.10184075650935, curviGrid->GetNode(0, 19).y, tolerance);
    ASSERT_NEAR(41.10185964080839, curviGrid->GetNode(0, 20).y, tolerance);
    ASSERT_NEAR(41.10188078420936, curviGrid->GetNode(0, 21).y, tolerance);
}

TEST(CurvilinearGridFromSplines, OrthogonalCurvilinearMeshTwoCrossingHighCurvature)
{
    std::vector<Point> firstSpline{{-103.468336664918, 420.606335102062},
                                   {111.984583950396, 845.689124424167},
                                   {1465.84415268176, 1608.50892444055}};

    auto splines = std::make_shared<Splines>(Projection::cartesian);
    splines->AddSpline(firstSpline);

    std::vector<Point> secondSpline{{-333.478887051536, 921.388799234953},
                                    {167.303577081355, 490.482958004326}};

    splines->AddSpline(secondSpline);

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

    ASSERT_NEAR(191.1109837821975, curviGrid->GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(48.40901089350142, curviGrid->GetNode(0, 1).x, tolerance);
    ASSERT_NEAR(-43.83755887944302, curviGrid->GetNode(0, 2).x, tolerance);
    ASSERT_NEAR(-103.468336664918, curviGrid->GetNode(0, 3).x, tolerance);
    ASSERT_NEAR(-163.099114450393, curviGrid->GetNode(0, 4).x, tolerance);
    ASSERT_NEAR(-270.8336555453428, curviGrid->GetNode(0, 5).x, tolerance);
    ASSERT_NEAR(-465.4769551295826, curviGrid->GetNode(0, 6).x, tolerance);

    ASSERT_NEAR(453.4383212979242, curviGrid->GetNode(0, 0).y, tolerance);
    ASSERT_NEAR(437.5336434105126, curviGrid->GetNode(0, 1).y, tolerance);
    ASSERT_NEAR(427.2524121930461, curviGrid->GetNode(0, 2).y, tolerance);
    ASSERT_NEAR(420.606335102062, curviGrid->GetNode(0, 3).y, tolerance);
    ASSERT_NEAR(413.9602580110779, curviGrid->GetNode(0, 4).y, tolerance);
    ASSERT_NEAR(401.9528334550086, curviGrid->GetNode(0, 5).y, tolerance);
    ASSERT_NEAR(380.2590970604539, curviGrid->GetNode(0, 6).y, tolerance);

    ASSERT_NEAR(185.1393938946496, curviGrid->GetNode(1, 0).x, tolerance);
    ASSERT_NEAR(39.18313770594886, curviGrid->GetNode(1, 1).x, tolerance);
    ASSERT_NEAR(-55.16709234174845, curviGrid->GetNode(1, 2).x, tolerance);
    ASSERT_NEAR(-116.1577353090799, curviGrid->GetNode(1, 3).x, tolerance);
    ASSERT_NEAR(-177.1483782764113, curviGrid->GetNode(1, 4).x, tolerance);
    ASSERT_NEAR(-287.3397790126249, curviGrid->GetNode(1, 5).x, tolerance);
    ASSERT_NEAR(-486.4218710429149, curviGrid->GetNode(1, 6).x, tolerance);

    ASSERT_NEAR(507.0173792497823, curviGrid->GetNode(1, 0).y, tolerance);
    ASSERT_NEAR(520.311195700755, curviGrid->GetNode(1, 1).y, tolerance);
    ASSERT_NEAR(528.9046926368475, curviGrid->GetNode(1, 2).y, tolerance);
    ASSERT_NEAR(534.4597708089422, curviGrid->GetNode(1, 3).y, tolerance);
    ASSERT_NEAR(540.0148489810368, curviGrid->GetNode(1, 4).y, tolerance);
    ASSERT_NEAR(550.0511728366511, curviGrid->GetNode(1, 5).y, tolerance);
    ASSERT_NEAR(568.1837343262024, curviGrid->GetNode(1, 6).y, tolerance);
}

TEST(CurvilinearGridFromSplines, OrthogonalCurvilinearMeshTwoCrossingHighCurvatureRemoveSkinnyTriangles)
{
    std::vector<Point> firstSpline{{-103.468336664918, 420.606335102062},
                                   {111.984583950396, 845.689124424167},
                                   {1465.84415268176, 1608.50892444055}};

    auto splines = std::make_shared<Splines>(Projection::cartesian);
    splines->AddSpline(firstSpline);

    std::vector<Point> secondSpline{{-333.478887051536, 921.388799234953},
                                    {167.303577081355, 490.482958004326}};

    splines->AddSpline(secondSpline);

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

    ASSERT_NEAR(191.1109837821975, curviGrid->GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(48.40901089350142, curviGrid->GetNode(0, 1).x, tolerance);
    ASSERT_NEAR(-43.83755887944302, curviGrid->GetNode(0, 2).x, tolerance);
    ASSERT_NEAR(-103.468336664918, curviGrid->GetNode(0, 3).x, tolerance);
    ASSERT_NEAR(-163.099114450393, curviGrid->GetNode(0, 4).x, tolerance);
    ASSERT_NEAR(-270.8336555453428, curviGrid->GetNode(0, 5).x, tolerance);
    ASSERT_NEAR(-465.4769551295826, curviGrid->GetNode(0, 6).x, tolerance);

    ASSERT_NEAR(453.4383212979242, curviGrid->GetNode(0, 0).y, tolerance);
    ASSERT_NEAR(437.5336434105126, curviGrid->GetNode(0, 1).y, tolerance);
    ASSERT_NEAR(427.2524121930461, curviGrid->GetNode(0, 2).y, tolerance);
    ASSERT_NEAR(420.606335102062, curviGrid->GetNode(0, 3).y, tolerance);
    ASSERT_NEAR(413.9602580110779, curviGrid->GetNode(0, 4).y, tolerance);
    ASSERT_NEAR(401.9528334550086, curviGrid->GetNode(0, 5).y, tolerance);
    ASSERT_NEAR(380.2590970604539, curviGrid->GetNode(0, 6).y, tolerance);

    ASSERT_NEAR(185.1393938946496, curviGrid->GetNode(1, 0).x, tolerance);
    ASSERT_NEAR(39.18313770594886, curviGrid->GetNode(1, 1).x, tolerance);
    ASSERT_NEAR(-55.16709234174845, curviGrid->GetNode(1, 2).x, tolerance);
    ASSERT_NEAR(-116.1577353090799, curviGrid->GetNode(1, 3).x, tolerance);
    ASSERT_NEAR(-177.1483782764113, curviGrid->GetNode(1, 4).x, tolerance);
    ASSERT_NEAR(-287.3397790126249, curviGrid->GetNode(1, 5).x, tolerance);
    ASSERT_NEAR(-486.4218710429149, curviGrid->GetNode(1, 6).x, tolerance);

    ASSERT_NEAR(507.0173792497823, curviGrid->GetNode(1, 0).y, tolerance);
    ASSERT_NEAR(520.311195700755, curviGrid->GetNode(1, 1).y, tolerance);
    ASSERT_NEAR(528.9046926368475, curviGrid->GetNode(1, 2).y, tolerance);
    ASSERT_NEAR(534.4597708089422, curviGrid->GetNode(1, 3).y, tolerance);
    ASSERT_NEAR(540.0148489810368, curviGrid->GetNode(1, 4).y, tolerance);
    ASSERT_NEAR(550.0511728366511, curviGrid->GetNode(1, 5).y, tolerance);
    ASSERT_NEAR(568.1837343262024, curviGrid->GetNode(1, 6).y, tolerance);

    ASSERT_NEAR(185.2535774620607, curviGrid->GetNode(2, 0).x, tolerance);
    ASSERT_NEAR(55.17693856050623, curviGrid->GetNode(2, 1).x, tolerance);
    ASSERT_NEAR(-28.90826014674191, curviGrid->GetNode(2, 2).x, tolerance);
    ASSERT_NEAR(-83.26329807932888, curviGrid->GetNode(2, 3).x, tolerance);
    ASSERT_NEAR(-137.6183360119159, curviGrid->GetNode(2, 4).x, tolerance);
    ASSERT_NEAR(-235.8212319884337, curviGrid->GetNode(2, 5).x, tolerance);
    ASSERT_NEAR(-413.2437698980817, curviGrid->GetNode(2, 6).x, tolerance);

    ASSERT_NEAR(507.3963391703687, curviGrid->GetNode(2, 0).y, tolerance);
    ASSERT_NEAR(573.3924720975141, curviGrid->GetNode(2, 1).y, tolerance);
    ASSERT_NEAR(616.0542290275101, curviGrid->GetNode(2, 2).y, tolerance);
    ASSERT_NEAR(643.6319888693365, curviGrid->GetNode(2, 3).y, tolerance);
    ASSERT_NEAR(671.2097487111629, curviGrid->GetNode(2, 4).y, tolerance);
    ASSERT_NEAR(721.03430930426, curviGrid->GetNode(2, 5).y, tolerance);
    ASSERT_NEAR(811.0520211319655, curviGrid->GetNode(2, 6).y, tolerance);

    //--------------------------------
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
    splines->AddSpline(firstSpline);

    std::vector<Point> secondSpline{{7.6618112E+04, 3.7136337E+05},
                                    {7.6754253E+04, 3.7005301E+05},
                                    {7.7179694E+04, 3.6874265E+05},
                                    {7.8404966E+04, 3.6780668E+05},
                                    {7.9681290E+04, 3.6721107E+05},
                                    {8.0140766E+04, 3.6636018E+05},
                                    {7.9477078E+04, 3.6544123E+05},
                                    {7.8779354E+04, 3.6452228E+05}};

    splines->AddSpline(secondSpline);

    std::vector<Point> thirdSpline{{7.7281800E+04, 3.7144846E+05},
                                   {7.7366889E+04, 3.6984880E+05},
                                   {7.7928471E+04, 3.6874265E+05},
                                   {7.9153742E+04, 3.6792581E+05},
                                   {8.0242872E+04, 3.6722808E+05},
                                   {8.0481119E+04, 3.6641124E+05},
                                   {7.9970590E+04, 3.6542421E+05},
                                   {7.9579184E+04, 3.6484561E+05},
                                   {7.9170760E+04, 3.6431806E+05}};

    splines->AddSpline(thirdSpline);

    std::vector<Point> fourthSpline{{7.613792E+04, 3.712157E+05},
                                    {7.831719E+04, 3.710751E+05}};

    splines->AddSpline(fourthSpline);

    std::vector<Point> fifthSpline{{7.857202E+04, 3.649151E+05},
                                   {8.003072E+04, 3.641506E+05}};

    splines->AddSpline(fifthSpline);

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

    ASSERT_EQ(curviGrid->GetNumNodes(), 216);
    ASSERT_EQ(curviGrid->GetNumEdges(), 399);

    const double tolerance = 1e-6;

    ASSERT_NEAR(76628.0244427564, curviGrid->GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(76791.25671189229, curviGrid->GetNode(0, 1).x, tolerance);
    ASSERT_NEAR(76954.67526993192, curviGrid->GetNode(0, 2).x, tolerance);
    ASSERT_NEAR(77118.21229417545, curviGrid->GetNode(0, 3).x, tolerance);
    ASSERT_NEAR(77281.8, curviGrid->GetNode(0, 4).x, tolerance);
    ASSERT_NEAR(77456.02341662964, curviGrid->GetNode(0, 5).x, tolerance);
    ASSERT_NEAR(77630.24899749605, curviGrid->GetNode(0, 6).x, tolerance);
    ASSERT_NEAR(77804.47364128292, curviGrid->GetNode(0, 7).x, tolerance);
    ASSERT_NEAR(77978.6942050769, curviGrid->GetNode(0, 8).x, tolerance);

    ASSERT_NEAR(371425.1853797209, curviGrid->GetNode(0, 0).y, tolerance);
    ASSERT_NEAR(371435.9949955427, curviGrid->GetNode(0, 1).y, tolerance);
    ASSERT_NEAR(371443.4780677535, curviGrid->GetNode(0, 2).y, tolerance);
    ASSERT_NEAR(371447.6329664325, curviGrid->GetNode(0, 3).y, tolerance);
    ASSERT_NEAR(371448.46, curviGrid->GetNode(0, 4).y, tolerance);
    ASSERT_NEAR(371449.3408034386, curviGrid->GetNode(0, 5).y, tolerance);
    ASSERT_NEAR(371449.4880738311, curviGrid->GetNode(0, 6).y, tolerance);
    ASSERT_NEAR(371448.8979771081, curviGrid->GetNode(0, 7).y, tolerance);
    ASSERT_NEAR(371447.5676132051, curviGrid->GetNode(0, 8).y, tolerance);

    ASSERT_NEAR(76661.93274954233, curviGrid->GetNode(1, 0).x, tolerance);
    ASSERT_NEAR(76817.18684812896, curviGrid->GetNode(1, 1).x, tolerance);
    ASSERT_NEAR(76972.6219613088, curviGrid->GetNode(1, 2).x, tolerance);
    ASSERT_NEAR(77128.17494466474, curviGrid->GetNode(1, 3).x, tolerance);
    ASSERT_NEAR(77283.78265752396, curviGrid->GetNode(1, 4).x, tolerance);
    ASSERT_NEAR(77456.35492003536, curviGrid->GetNode(1, 5).x, tolerance);
    ASSERT_NEAR(77628.92068688275, curviGrid->GetNode(1, 6).x, tolerance);
    ASSERT_NEAR(77801.47895858958, curviGrid->GetNode(1, 7).x, tolerance);
    ASSERT_NEAR(77974.02895307761, curviGrid->GetNode(1, 8).x, tolerance);

    ASSERT_NEAR(371034.1516952458, curviGrid->GetNode(1, 0).y, tolerance);
    ASSERT_NEAR(371044.4329782871, curviGrid->GetNode(1, 1).y, tolerance);
    ASSERT_NEAR(371051.5504818924, curviGrid->GetNode(1, 2).y, tolerance);
    ASSERT_NEAR(371055.5025342674, curviGrid->GetNode(1, 3).y, tolerance);
    ASSERT_NEAR(371056.2892242038, curviGrid->GetNode(1, 4).y, tolerance);
    ASSERT_NEAR(371057.161680072, curviGrid->GetNode(1, 5).y, tolerance);
    ASSERT_NEAR(371057.3075474474, curviGrid->GetNode(1, 6).y, tolerance);
    ASSERT_NEAR(371056.7230947062, curviGrid->GetNode(1, 7).y, tolerance);
    ASSERT_NEAR(371055.4054874207, curviGrid->GetNode(1, 8).y, tolerance);

    ASSERT_NEAR(76697.93711080145, curviGrid->GetNode(2, 0).x, tolerance);
    ASSERT_NEAR(76846.19231052838, curviGrid->GetNode(2, 1).x, tolerance);
    ASSERT_NEAR(76994.63120583622, curviGrid->GetNode(2, 2).x, tolerance);
    ASSERT_NEAR(77143.21644135848, curviGrid->GetNode(2, 3).x, tolerance);
    ASSERT_NEAR(77291.91098375124, curviGrid->GetNode(2, 4).x, tolerance);
    ASSERT_NEAR(77461.38217435485, curviGrid->GetNode(2, 5).x, tolerance);
    ASSERT_NEAR(77630.83599384772, curviGrid->GetNode(2, 6).x, tolerance);
    ASSERT_NEAR(77800.27088568854, curviGrid->GetNode(2, 7).x, tolerance);
    ASSERT_NEAR(77969.68564989329, curviGrid->GetNode(2, 8).x, tolerance);

    ASSERT_NEAR(370636.3843794758, curviGrid->GetNode(2, 0).y, tolerance);
    ASSERT_NEAR(370647.2100657219, curviGrid->GetNode(2, 1).y, tolerance);
    ASSERT_NEAR(370655.4582144818, curviGrid->GetNode(2, 2).y, tolerance);
    ASSERT_NEAR(370661.1252352109, curviGrid->GetNode(2, 3).y, tolerance);
    ASSERT_NEAR(370664.2078549, curviGrid->GetNode(2, 4).y, tolerance);
    ASSERT_NEAR(370667.7211999126, curviGrid->GetNode(2, 5).y, tolerance);
    ASSERT_NEAR(370669.9086649796, curviGrid->GetNode(2, 6).y, tolerance);
    ASSERT_NEAR(370670.7463540983, curviGrid->GetNode(2, 7).y, tolerance);
    ASSERT_NEAR(370670.2161009402, curviGrid->GetNode(2, 8).y, tolerance);
}

meshkernel::Splines LoadSplines(const std::string& fileName)
{

    std::ifstream splineFile;
    splineFile.open(fileName.c_str());

    meshkernel::Splines splines(Projection::cartesian);
    std::string line;

    std::vector<meshkernel::Point> splinePoints;

    while (std::getline(splineFile, line))
    {
        size_t found = line.find("S00");

        if (found == std::string::npos)
        {
            found = line.find("L00");
        }

        if (found != std::string::npos)
        {
            std::getline(splineFile, line);
            std::istringstream sizes(line);

            meshkernel::UInt numPoints = 0;
            meshkernel::UInt numDim = 0;

            sizes >> numPoints;
            sizes >> numDim;

            splinePoints.clear();
            splinePoints.reserve(numPoints);

            for (meshkernel::UInt i = 0; i < numPoints; ++i)
            {
                std::getline(splineFile, line);
                std::istringstream values(line);
                double x;
                double y;
                values >> x;
                values >> y;
                splinePoints.emplace_back(meshkernel::Point(x, y));
            }

            splines.AddSpline(splinePoints);
        }
    }

    splineFile.close();
    return splines;
}

meshkernel::CurvilinearGrid LoadCurvilinearGrid(const std::string& fileName)
{
    std::ifstream splineFile;
    splineFile.open(fileName.c_str());

    meshkernel::Splines splines(Projection::cartesian);
    std::string line;

    std::vector<meshkernel::Point> splinePoints;

    // read header
    while (std::getline(splineFile, line))
    {
        if (line[0] != '*')
        {
            break;
        }
    }

    if (size_t found = line.find("Missing Value="); found != std::string::npos)
    {
        found = line.find("= ");
        std::istringstream missingValueLine(line.substr(found + 2));
        double missingValue;
        missingValueLine >> missingValue;
        std::cout << "found = " << found << "  " << missingValue << std::endl;
    }

    UInt rows = 0;
    UInt cols = 0;

    if (std::getline(splineFile, line))
    {
        std::istringstream sizeLine(line);
        sizeLine >> rows;
        sizeLine >> cols;
        std::cout << " sizes: " << rows << "  " << cols << std::endl;
    }

    if (std::getline(splineFile, line))
    {
        std::istringstream sizeLine(line);
        UInt z1, z2, z3;
        sizeLine >> z1;
        sizeLine >> z2;
        sizeLine >> z3;
        std::cout << " zeros: " << z1 << "  " << z2 << "  " << z3 << std::endl;
    }

    lin_alg::Matrix<Point> gridNodes(rows, cols);
    UInt currentRow = 0;
    UInt currentCol = 0;

    for (UInt c = 0; c < cols; ++c)
    {
        std::getline(splineFile, line);

        if (size_t found = line.find("ETA="); found != std::string::npos)
        {
            // Ignore ETA=
            found += 4;
            std::stringstream values(line.substr(found));
            values >> currentCol;
            --currentCol;
            currentRow = 0;
            double coord;

            while (values >> coord)
            {
                gridNodes(currentRow, currentCol).x = coord;
                ++currentRow;
            }

            while (currentRow < rows)
            {
                std::getline(splineFile, line);
                std::stringstream values(line);
                double coord;

                while (values >> coord)
                {
                    gridNodes(currentRow, currentCol).x = coord;
                    ++currentRow;
                }
            }
        }
    }

    for (UInt c = 0; c < cols; ++c)
    {
        std::getline(splineFile, line);

        if (size_t found = line.find("ETA="); found != std::string::npos)
        {

            // Ignore ETA=
            found += 4;
            std::stringstream values(line.substr(found));
            values >> currentCol;
            --currentCol;
            currentRow = 0;
            double coord;

            while (values >> coord)
            {
                gridNodes(currentRow, currentCol).y = coord;
                ++currentRow;
            }

            while (currentRow < rows)
            {
                std::getline(splineFile, line);
                std::stringstream values(line);
                double coord;

                while (values >> coord)
                {
                    gridNodes(currentRow, currentCol).y = coord;
                    ++currentRow;
                }
            }
        }
    }

    meshkernel::CurvilinearGrid grid(gridNodes, meshkernel::Projection::cartesian);

    return grid;
}

TEST(CurvilinearGridFromSplines, GenerateSimpleGridFromSplines)
{
    // Generate a grid from four splines forming a square
    namespace mk = meshkernel;

    // Bottom boundary
    std::vector<mk::Point> spline1{{-1.0, 0.0}, {11.0, 0.0}};

    // right boundary
    std::vector<mk::Point> spline2{{10.0, 11.0}, {10.0, -1.0}};

    // top boundary
    std::vector<mk::Point> spline3{{11.0, 10.0}, {-1.0, 10.0}};

    // left boundary
    std::vector<mk::Point> spline4{{0.0, 11.0}, {0.0, -1.0}};

    auto splines = std::make_shared<Splines>(Projection::cartesian);
    splines->AddSpline(spline1);
    splines->AddSpline(spline2);
    splines->AddSpline(spline3);
    splines->AddSpline(spline4);

    mk::CurvilinearParameters curvilinearParameters;

    curvilinearParameters.m_refinement = 5;
    curvilinearParameters.n_refinement = 20;

    mk::CurvilinearGridSplineToGrid splinesToGrid;

    mk::CurvilinearGrid grid(splinesToGrid.Compute(*splines, curvilinearParameters));

    lin_alg::Matrix<Point> gridNodes = grid.GetNodes();

    constexpr double tolerance = 1.0e-4;

    double meshDeltaX = 10.0 / static_cast<double>(curvilinearParameters.m_refinement);
    double meshDeltaY = 10.0 / static_cast<double>(curvilinearParameters.n_refinement);
    double yCoord = 0.0;

    for (UInt i = 0; i < gridNodes.rows(); ++i)
    {
        double xCoord = 0.0;

        for (UInt j = 0; j < gridNodes.cols(); ++j)
        {
            EXPECT_NEAR(xCoord, gridNodes(i, j).x, tolerance);
            EXPECT_NEAR(yCoord, gridNodes(i, j).y, tolerance);
            xCoord += meshDeltaX;
        }

        yCoord += meshDeltaY;
    }
}

TEST(CurvilinearGridFromSplines, GridFromSeventySplines)
{
    // Test generating a more complicated grid from a set of much more complcated splines
    // Only check the size of the grid is correct and the number of valid grid nodes.
    namespace mk = meshkernel;
    auto splines = std::make_shared<Splines>(LoadSplines(TEST_FOLDER + "/data/CurvilinearGrids/seventy_splines.spl"));

    mk::CurvilinearParameters curvilinearParameters;

    curvilinearParameters.m_refinement = 5;
    curvilinearParameters.n_refinement = 5;

    mk::CurvilinearGridSplineToGrid splinesToGrid;

    mk::CurvilinearGrid grid(splinesToGrid.Compute(*splines, curvilinearParameters));

    auto computedPoints = grid.ComputeNodes();
    auto computedEdges = grid.ComputeEdges();

    mk::UInt validPointCount = 0;

    for (size_t i = 0; i < computedPoints.size(); ++i)
    {
        if (computedPoints[i].IsValid())
        {
            ++validPointCount;
        }
    }

    EXPECT_EQ(grid.NumN(), 36);
    EXPECT_EQ(grid.NumM(), 301);
    EXPECT_EQ(validPointCount, 4941);
}

TEST(CurvilinearGridFromSplines, GenerateGridWithIllDefinedSplines)
{
    // Attempt to generate a grid with 5 splines, 2 of which form a single
    // boundary along the top of the domain.
    // Should raise a AlgorithmError exception

    namespace mk = meshkernel;

    // Bottom boundary
    std::vector<mk::Point> spline1{{-1.0, 0.0}, {11.0, 0.0}};

    // right boundary
    std::vector<mk::Point> spline2{{10.0, 11.0}, {10.0, -1.0}};

    // top boundary 1
    std::vector<mk::Point> spline3{{11.0, 9.0}, {4.0, 11.5}};

    // top boundary 2
    std::vector<mk::Point> spline4{{6.0, 11.5}, {-1.0, 9.0}};

    // left boundary
    std::vector<mk::Point> spline5{{0.0, 11.0}, {0.0, -1.0}};

    auto splines = std::make_shared<Splines>(Projection::cartesian);
    splines->AddSpline(spline1);
    splines->AddSpline(spline2);
    splines->AddSpline(spline3);
    splines->AddSpline(spline4);
    splines->AddSpline(spline5);

    mk::CurvilinearParameters curvilinearParameters;

    curvilinearParameters.m_refinement = 5;
    curvilinearParameters.n_refinement = 10;

    mk::CurvilinearGridSplineToGrid splinesToGrid;

    EXPECT_THROW([[maybe_unused]] auto grid = splinesToGrid.Compute(*splines, curvilinearParameters), mk::AlgorithmError);
}

TEST(CurvilinearGridFromSplines, GenerateGridWithDisconectedSpline)
{
    // Attempt to generate a grid with 4 splines forming a square and
    // a fifth disconnected spline.
    // Should raise a AlgorithmError exception

    namespace mk = meshkernel;

    // Bottom boundary
    std::vector<mk::Point> spline1{{-1.0, 0.0}, {11.0, 0.0}};

    // right boundary
    std::vector<mk::Point> spline2{{10.0, -1.0}, {10.0, 11.0}};

    // top boundary
    std::vector<mk::Point> spline3{{11.0, 10.0}, {-1.0, 10.0}};

    // left boundary
    std::vector<mk::Point> spline4{{0.0, 11.0}, {0.0, -1.0}};

    // unattached spline
    std::vector<mk::Point> spline5{{15.0, 15.0}, {15.0, -1.0}};

    auto splines = std::make_shared<Splines>(Projection::cartesian);
    splines->AddSpline(spline1);
    splines->AddSpline(spline2);
    splines->AddSpline(spline3);
    splines->AddSpline(spline4);
    splines->AddSpline(spline5);

    mk::CurvilinearParameters curvilinearParameters;

    curvilinearParameters.m_refinement = 5;
    curvilinearParameters.n_refinement = 10;

    mk::CurvilinearGridSplineToGrid splinesToGrid;
    mk::CurvilinearGrid grid(Projection::cartesian);

    EXPECT_THROW([[maybe_unused]] auto grid = splinesToGrid.Compute(*splines, curvilinearParameters), mk::AlgorithmError);
}

TEST(CurvilinearGridFromSplines, GenerateGridWithThreeSplines)
{
    // Attempt to generate a grid with 3 splines
    // Should raise a ConstraintError exception.
    // At least four splines are required to generate a curvilinear grid

    namespace mk = meshkernel;

    std::vector<mk::Point> spline1{{-1.0, 0.0}, {11.0, 0.0}};
    std::vector<mk::Point> spline2{{10.0, -1.0}, {10.0, 11.0}};
    std::vector<mk::Point> spline3{{11.0, 10.0}, {-1.0, 10.0}};

    auto splines = std::make_shared<Splines>(Projection::cartesian);
    splines->AddSpline(spline1);
    splines->AddSpline(spline2);
    splines->AddSpline(spline3);

    mk::CurvilinearParameters curvilinearParameters;

    curvilinearParameters.m_refinement = 5;
    curvilinearParameters.n_refinement = 5;

    mk::CurvilinearGridSplineToGrid splinesToGrid;

    EXPECT_THROW([[maybe_unused]] auto grid = splinesToGrid.Compute(*splines, curvilinearParameters), mk::ConstraintError);
}

TEST(CurvilinearGridFromSplines, GenerateGridWithSplinesTooShort)
{
    // Attempt to generate a grid with 4 splines, 1 of which has only a single point
    // Should raise a ConstraintError exception.
    // All splines are required to have at least 2 points in order to generate a curvilinear grid

    namespace mk = meshkernel;

    std::vector<mk::Point> spline1{{-1.0, 0.0}, {11.0, 0.0}};
    std::vector<mk::Point> spline2{{10.0, -1.0}, {10.0, 11.0}};
    std::vector<mk::Point> spline3{{11.0, 10.0}, {-1.0, 10.0}};
    std::vector<mk::Point> spline4{{0.0, 11.0}};

    auto splines = std::make_shared<Splines>(Projection::cartesian);
    splines->AddSpline(spline1);
    splines->AddSpline(spline2);
    splines->AddSpline(spline3);
    splines->AddSpline(spline4);

    mk::CurvilinearParameters curvilinearParameters;

    curvilinearParameters.m_refinement = 5;
    curvilinearParameters.n_refinement = 5;

    mk::CurvilinearGridSplineToGrid splinesToGrid;

    EXPECT_THROW([[maybe_unused]] auto grid = splinesToGrid.Compute(*splines, curvilinearParameters), mk::ConstraintError);
}

TEST(CurvilinearGridFromSplines, GenerateGridWithHighRefinementFactor)
{
    // Attempt to generate a grid with 4 splines forming a square
    // with the refinement factor being too high
    // Should raise a ConstraintError exception

    namespace mk = meshkernel;

    std::vector<mk::Point> spline1{{-1.0, 0.0}, {11.0, 0.0}};
    std::vector<mk::Point> spline2{{10.0, -1.0}, {10.0, 11.0}};
    std::vector<mk::Point> spline3{{11.0, 10.0}, {-1.0, 10.0}};
    std::vector<mk::Point> spline4{{0.0, 11.0}, {0.0, -1.0}};

    auto splines = std::make_shared<Splines>(Projection::cartesian);
    splines->AddSpline(spline1);
    splines->AddSpline(spline2);
    splines->AddSpline(spline3);
    splines->AddSpline(spline4);

    mk::CurvilinearParameters curvilinearParameters;

    mk::CurvilinearGridSplineToGrid splinesToGrid;

    // First check the m-refinement
    curvilinearParameters.m_refinement = 2000;
    curvilinearParameters.n_refinement = 10;

    EXPECT_THROW([[maybe_unused]] auto gridM = splinesToGrid.Compute(*splines, curvilinearParameters), mk::ConstraintError);

    // Then check the n-refinement
    curvilinearParameters.m_refinement = 10;
    curvilinearParameters.n_refinement = 2000;

    EXPECT_THROW([[maybe_unused]] auto gridN = splinesToGrid.Compute(*splines, curvilinearParameters), mk::ConstraintError);
}

TEST(CurvilinearGridFromSplines, WithoutOverlappingElements)
{

    namespace mk = meshkernel;
    mk::CurvilinearParameters curvilinearParameters;
    mk::SplinesToCurvilinearParameters splinesToCurvilinearParameters;

    splinesToCurvilinearParameters.aspect_ratio = 0.2;
    splinesToCurvilinearParameters.aspect_ratio_grow_factor = 1.1;
    splinesToCurvilinearParameters.average_width = 100.0;
    splinesToCurvilinearParameters.nodes_on_top_of_each_other_tolerance = 1e-4;
    splinesToCurvilinearParameters.min_cosine_crossing_angles = 0.95;

    splinesToCurvilinearParameters.check_front_collisions = 1;
    splinesToCurvilinearParameters.grow_grid_outside = 0;
    splinesToCurvilinearParameters.curvature_adapted_grid_spacing = 0;
    splinesToCurvilinearParameters.remove_skinny_triangles = 0;

    curvilinearParameters.m_refinement = 10;
    curvilinearParameters.n_refinement = 3;

    std::vector<mk::Point> spline1{{9.500679E+00, 6.187624E+01},
                                   {4.130032E+02, 5.362616E+01},
                                   {7.370052E+02, 5.212615E+01},
                                   {8.352558E+02, 1.248769E+02},
                                   {7.865055E+02, 1.788774E+02},
                                   {5.472540E+02, 2.073777E+02},
                                   {2.825024E+02, 2.171278E+02}};

    std::vector<mk::Point> spline2{{6.500102E+01, 1.848775E+02},
                                   {4.700092E+01, -2.137459E+01}};

    auto splines = std::make_shared<Splines>(Projection::cartesian);
    splines->AddSpline(spline1);
    splines->AddSpline(spline2);

    const double tolerance = 1.0e-8;

    mk::CurvilinearGridFromSplines splinesToGrid(splines, curvilinearParameters, splinesToCurvilinearParameters);

    [[maybe_unused]] auto grid = splinesToGrid.Compute();
    std::cout.precision(12);

    std::vector<double> xNodes({8.49759082988, 8.9343870837, 9.25917578417, 9.500679, 9.74218221583, 10.1851367696,
                                10.9975844882, 153.068185819, 153.635954468, 154.058130345, 154.372047768, 154.685965192,
                                155.261738693, 156.157323803, 296.983176781, 297.960035756, 298.686399031, 299.226501153,
                                299.766603276, 300.757234784, -999, 440.860416881, 442.236601148, 443.259890829,
                                444.020778596, 444.781666362, 446.177253165, -999, 588.547149155, 588.646587946,
                                588.720527673, 588.775507054, 588.830486436, 588.931327194, -999, 764.201804305,
                                750.779207691, 740.798564145, 733.377254509, 725.955944873, 712.344105847, -999,
                                994.126426047, 925.996638475, 875.337355168, 837.668619119, 799.999883069, 730.909547721,
                                -999, 737.988284591, 728.488463362, 721.424678897, 716.172258903, 710.91983891,
                                701.286081667, -999, 578.243881409, 575.530316198, 573.512589831, 572.012268528,
                                570.511947224, 567.760124019, -999, 430.604350702, 429.17794234, 428.117307514,
                                427.32865101, 426.539994506, 425.093475507, -999, 284.744766133, 283.771578515,
                                283.047945152, 282.50987291, 281.971800668, 280.984892278, -999});

    std::vector<double> yNodes({-21.1880814511, 14.9824023768, 41.8776981482, 61.87624, 81.8747818518, 118.555225846,
                                185.832879467, -22.9339219161, 13.2349802838, 40.1291000037, 60.1267673785, 80.1244347534,
                                116.803274821, 184.079924471, -25.7143270936, 10.4466715734, 37.3349144613, 57.3282119987,
                                77.3215095362, 113.99233464, -999, -30.7085566668, 5.43858072396, 32.3165167792,
                                52.3021504624, 72.2877841457, 108.944552546, -999, -36.827600308, -0.627564047093,
                                26.2897060344, 46.3045873294, 66.3194686244, 103.029881705, -999, -30.5162375246,
                                5.1979449971, 31.7539486913, 51.5002030477, 71.2464574042, 107.464166795, -999,
                                156.5921184, 147.78661234, 141.239099981, 136.370564566, 131.502029152, 122.572376084,
                                -999, 272.831103885, 237.418518661, 211.086773745, 191.507271545, 171.927769345,
                                136.01590941, -999, 288.562178287, 252.481306123, 225.652642919, 205.703647062,
                                185.754651206, 149.165082292, -999, 296.242370889, 260.094540963, 233.21608996,
                                213.230073376, 193.244056793, 156.586586094, -999, 300.167941502, 264.007913886,
                                237.120393043, 217.127632397, 197.134871751, 160.465031388, -999});

    std::vector<mk::UInt> edgeStart({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
                                     12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
                                     24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35,
                                     36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
                                     48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
                                     60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 0, 1,
                                     2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 14, 15,
                                     16, 17, 18, 19, 21, 22, 23, 24, 25, 26, 28, 29,
                                     30, 31, 32, 33, 35, 36, 37, 38, 39, 40, 42, 43,
                                     44, 45, 46, 47, 49, 50, 51, 52, 53, 54, 56, 57,
                                     58, 59, 60, 61, 63, 64, 65, 66, 67, 68, 70, 71,
                                     72, 73, 74, 75});

    std::vector<mk::UInt> edgeEnd({7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
                                   19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
                                   31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42,
                                   43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54,
                                   55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66,
                                   67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 1, 2,
                                   3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 15, 16,
                                   17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 29, 30,
                                   31, 32, 33, 34, 36, 37, 38, 39, 40, 41, 43, 44,
                                   45, 46, 47, 48, 50, 51, 52, 53, 54, 55, 57, 58,
                                   59, 60, 61, 62, 64, 65, 66, 67, 68, 69, 71, 72,
                                   73, 74, 75, 76});

    auto nodes = grid->ComputeNodes();
    auto edges = grid->ComputeEdges();

    ASSERT_EQ(nodes.size(), xNodes.size());

    for (size_t i = 0; i < nodes.size(); ++i)
    {
        EXPECT_NEAR(nodes[i].x, xNodes[i], tolerance);
        EXPECT_NEAR(nodes[i].y, yNodes[i], tolerance);
    }

    ASSERT_EQ(edges.size(), edgeStart.size());

    for (size_t i = 0; i < edges.size(); ++i)
    {
        EXPECT_EQ(edges[i].first, edgeStart[i]);
        EXPECT_EQ(edges[i].second, edgeEnd[i]);
    }
}

TEST(CurvilinearGridFromSplines, GridFromFFiveSplines)
{
    // Test generating a more complicated grid from a set of much more complcated splines
    // Only check the size of the grid is correct and the number of valid grid nodes.
    namespace mk = meshkernel;
    auto splines = std::make_shared<Splines>(LoadSplines(TEST_FOLDER + "/data/CurvilinearGrids/five_splines.spl"));

    mk::CurvilinearParameters curvilinearParameters;
    mk::SplinesToCurvilinearParameters splineParameters{.aspect_ratio = 1.0,
                                                        .aspect_ratio_grow_factor = 1.0,
                                                        .average_width = 25.0};

    curvilinearParameters.m_refinement = 40;
    curvilinearParameters.n_refinement = 40;

    mk::CurvilinearGridFromSplines gridFromSplines(splines, curvilinearParameters, splineParameters);

    auto grid = gridFromSplines.Compute();

    auto computedPoints = grid->ComputeNodes();
    auto computedEdges = grid->ComputeEdges();

    mk::UInt validPointCount = 0;

    for (size_t i = 0; i < computedPoints.size(); ++i)
    {
        if (computedPoints[i].IsValid())
        {
            ++validPointCount;
        }
    }

    mk::Print(computedPoints, computedEdges);

    // EXPECT_EQ(grid->NumN(), 36);
    // EXPECT_EQ(grid->NumM(), 301);
    // EXPECT_EQ(validPointCount, 4941);
}
