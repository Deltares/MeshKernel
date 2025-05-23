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
#include <MeshKernel/CurvilinearGrid/CurvilinearGridFromPolygon.hpp>
#include <MeshKernel/Polygon.hpp>
#include <MeshKernel/Polygons.hpp>

TEST(CurvilinearGridFromPolygon, ComputeGridInPolygonWithFourthSide)
{
    std::vector<meshkernel::Point> polygonPoints{{273.502319, 478.880432},
                                                 {274.252319, 325.128906},
                                                 {275.002350, 172.127350},
                                                 {458.003479, 157.127213},
                                                 {719.005127, 157.127213},
                                                 {741.505249, 328.128937},
                                                 {710.755066, 490.880554},
                                                 {507.503784, 494.630615},
                                                 {305.002533, 493.130615},
                                                 {273.502319, 478.880432}};

    meshkernel::Polygon polygon(polygonPoints, meshkernel::Projection::cartesian);

    meshkernel::CurvilinearGridFromPolygon curvilinearGridFromPolygon(polygon);

    const auto curvilinearGrid = curvilinearGridFromPolygon.Compute(0, 2, 4, true);

    // check the values
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(273.50231900000000, curvilinearGrid->GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(305.00253300000003, curvilinearGrid->GetNode(0, 1).x, tolerance);
    ASSERT_NEAR(507.50378400000000, curvilinearGrid->GetNode(0, 2).x, tolerance);

    ASSERT_NEAR(478.88043199999998, curvilinearGrid->GetNode(0, 0).y, tolerance);
    ASSERT_NEAR(493.13061499999998, curvilinearGrid->GetNode(0, 1).y, tolerance);
    ASSERT_NEAR(494.63061499999998, curvilinearGrid->GetNode(0, 2).y, tolerance);

    ASSERT_NEAR(274.25231900000000, curvilinearGrid->GetNode(1, 0).x, tolerance);
    ASSERT_NEAR(410.51616175207897, curvilinearGrid->GetNode(1, 1).x, tolerance);
    ASSERT_NEAR(741.50524900000005, curvilinearGrid->GetNode(1, 2).x, tolerance);

    ASSERT_NEAR(325.12890599999997, curvilinearGrid->GetNode(1, 0).y, tolerance);
    ASSERT_NEAR(314.33420290273324, curvilinearGrid->GetNode(1, 1).y, tolerance);
    ASSERT_NEAR(328.12893700000001, curvilinearGrid->GetNode(1, 2).y, tolerance);
}

TEST(CurvilinearGridFromPolygon, ComputeGridInPolygonWithoutFourthSide)
{
    std::vector<meshkernel::Point> polygonPoints{{273.502319, 478.880432},
                                                 {274.252319, 325.128906},
                                                 {275.002350, 172.127350},
                                                 {458.003479, 157.127213},
                                                 {719.005127, 157.127213},
                                                 {741.505249, 328.128937},
                                                 {710.755066, 490.880554},
                                                 {507.503784, 494.630615},
                                                 {305.002533, 493.130615},
                                                 {273.502319, 478.880432}};

    meshkernel::Polygon polygon(polygonPoints, meshkernel::Projection::cartesian);

    meshkernel::CurvilinearGridFromPolygon curvilinearGridFromPolygon(polygon);

    const auto curvilinearGrid = curvilinearGridFromPolygon.Compute(0, 2, 4, false);

    // check the values
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(273.50231900000000, curvilinearGrid->GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(492.12869250000006, curvilinearGrid->GetNode(0, 1).x, tolerance);
    ASSERT_NEAR(710.75506600000006, curvilinearGrid->GetNode(0, 2).x, tolerance);

    ASSERT_NEAR(478.88043199999998, curvilinearGrid->GetNode(0, 0).y, tolerance);
    ASSERT_NEAR(484.88049300000000, curvilinearGrid->GetNode(0, 1).y, tolerance);
    ASSERT_NEAR(490.88055400000002, curvilinearGrid->GetNode(0, 2).y, tolerance);

    ASSERT_NEAR(274.25231900000000, curvilinearGrid->GetNode(1, 0).x, tolerance);
    ASSERT_NEAR(481.37408996173241, curvilinearGrid->GetNode(1, 1).x, tolerance);
    ASSERT_NEAR(741.50524900000005, curvilinearGrid->GetNode(1, 2).x, tolerance);

    ASSERT_NEAR(325.12890599999997, curvilinearGrid->GetNode(1, 0).y, tolerance);
    ASSERT_NEAR(322.93773596204318, curvilinearGrid->GetNode(1, 1).y, tolerance);
    ASSERT_NEAR(328.12893700000001, curvilinearGrid->GetNode(1, 2).y, tolerance);
}

TEST(CurvilinearGridFromPolygon, ComputeGridTriangle)
{
    std::vector<meshkernel::Point> const polygonPoints{{444.504791, 437.155945},
                                                       {427.731781, 382.745758},
                                                       {405.640503, 317.699005},
                                                       {381.094666, 262.470612},
                                                       {451.050354, 262.879700},
                                                       {528.778931, 263.288788},
                                                       {593.416260, 266.561584},
                                                       {558.643005, 324.653687},
                                                       {526.733398, 377.836578},
                                                       {444.504791, 437.155945}};

    meshkernel::Polygon polygon(polygonPoints, meshkernel::Projection::cartesian);

    meshkernel::CurvilinearGridFromPolygon curvilinearGridFromPolygon(polygon);

    const auto curvilinearGrid = curvilinearGridFromPolygon.Compute(0, 3, 6);

    // check the values
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(444.50479100000001, curvilinearGrid->GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(558.64300500000002, curvilinearGrid->GetNode(0, 1).x, tolerance);
    ASSERT_NEAR(593.41625999999997, curvilinearGrid->GetNode(0, 2).x, tolerance);

    ASSERT_NEAR(437.15594499999997, curvilinearGrid->GetNode(0, 0).y, tolerance);
    ASSERT_NEAR(324.65368699999999, curvilinearGrid->GetNode(0, 1).y, tolerance);
    ASSERT_NEAR(266.56158399999998, curvilinearGrid->GetNode(0, 2).y, tolerance);

    ASSERT_NEAR(427.73178100000001, curvilinearGrid->GetNode(1, 0).x, tolerance);
    ASSERT_NEAR(473.00523900000007, curvilinearGrid->GetNode(1, 1).x, tolerance);
    ASSERT_NEAR(528.77893099999994, curvilinearGrid->GetNode(1, 2).x, tolerance);

    ASSERT_NEAR(382.74575800000002, curvilinearGrid->GetNode(1, 0).y, tolerance);
    ASSERT_NEAR(322.06271366666670, curvilinearGrid->GetNode(1, 1).y, tolerance);
    ASSERT_NEAR(263.28878800000001, curvilinearGrid->GetNode(1, 2).y, tolerance);
}
