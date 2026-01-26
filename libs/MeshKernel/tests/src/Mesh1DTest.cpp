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

#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh1D.hpp>
#include <MeshKernel/Network1D.hpp>
#include <TestUtils/MakeMeshes.hpp>

TEST(Mesh1D, GenerateMeshFromPolyLines_WithOverlappingNodes_ShouldRemoveOverlappingNodes)
{
    // 1 Setup
    std::vector<std::vector<meshkernel::Point>> polyLines{
        {{0.0, 0.0},
         {10.0, 0.0},
         {20.0, 0.0}},
        {{10.0, -10.0},
         {10.0, 0.0},
         {10.0, 10.0}}};

    std::vector<std::vector<double>> fixedChaninagesOnPolyline(polyLines.size());
    double const offset = 5.0;
    double const minFaceSize = 0.01;
    double const offsetFromFixedChainages = 1.0;

    meshkernel::Network1D Network1D(polyLines, meshkernel::Projection::cartesian);
    Network1D.ComputeFixedChainages(fixedChaninagesOnPolyline, minFaceSize, offsetFromFixedChainages);
    Network1D.ComputeOffsettedChainages(offset);

    // 2 Execution
    const auto mesh = meshkernel::Mesh1D(Network1D, minFaceSize);

    // 3 Assertion
    const auto tolerance = 1e-6;
    EXPECT_NEAR(0.00, mesh.Node(0).x, tolerance);
    EXPECT_NEAR(5.00, mesh.Node(1).x, tolerance);
    EXPECT_NEAR(10.0, mesh.Node(2).x, tolerance);
    EXPECT_NEAR(15.0, mesh.Node(3).x, tolerance);
    EXPECT_NEAR(20.0, mesh.Node(4).x, tolerance);
    EXPECT_NEAR(10.0, mesh.Node(5).x, tolerance);
    EXPECT_NEAR(10.0, mesh.Node(6).x, tolerance);
    EXPECT_NEAR(meshkernel::constants::missing::doubleValue, mesh.Node(7).x, tolerance);
    EXPECT_NEAR(10.0, mesh.Node(8).x, tolerance);
    EXPECT_NEAR(10.0, mesh.Node(9).x, tolerance);

    EXPECT_NEAR(0.0, mesh.Node(0).y, tolerance);
    EXPECT_NEAR(0.0, mesh.Node(1).y, tolerance);
    EXPECT_NEAR(0.0, mesh.Node(2).y, tolerance);
    EXPECT_NEAR(0.0, mesh.Node(3).y, tolerance);
    EXPECT_NEAR(0.0, mesh.Node(4).y, tolerance);
    EXPECT_NEAR(-10.0, mesh.Node(5).y, tolerance);
    EXPECT_NEAR(-5.0, mesh.Node(6).y, tolerance);
    EXPECT_NEAR(meshkernel::constants::missing::doubleValue, mesh.Node(7).y, tolerance);
    EXPECT_NEAR(5.0, mesh.Node(8).y, tolerance);
    EXPECT_NEAR(10.0, mesh.Node(9).y, tolerance);
}

TEST(Mesh1D, GenerateMeshFromPolyLines_WithInexactOffset_ShouldGenerateMesh)
{
    // 1 Setup
    std::vector<std::vector<meshkernel::Point>> polyLines{
        {{0.0, 0.0},
         {10.0, 0.0},
         {20.0, 0.0}},
        {{10.0, -10.0},
         {10.0, 0.0},
         {10.0, 10.0}}};

    std::vector<std::vector<double>> fixedChaninagesOnPolyline(polyLines.size());
    double const offset = 7.0;
    double const minFaceSize = 0.01;
    double const offsetFromFixedChainages = 1.0;

    meshkernel::Network1D Network1D(polyLines, meshkernel::Projection::cartesian);
    Network1D.ComputeFixedChainages(fixedChaninagesOnPolyline, minFaceSize, offsetFromFixedChainages);
    Network1D.ComputeOffsettedChainages(offset);

    // 2 Execution
    const auto mesh = meshkernel::Mesh1D(Network1D, minFaceSize);

    // 3 Assertion
    const auto tolerance = 1e-6;
    ASSERT_NEAR(0.0000000000000000, mesh.Node(0).x, tolerance);
    ASSERT_NEAR(6.6666666666666670, mesh.Node(1).x, tolerance);
    ASSERT_NEAR(13.333333333333334, mesh.Node(2).x, tolerance);
    ASSERT_NEAR(20.000000000000000, mesh.Node(3).x, tolerance);
    ASSERT_NEAR(10.000000000000000, mesh.Node(4).x, tolerance);
    ASSERT_NEAR(10.000000000000000, mesh.Node(5).x, tolerance);
    ASSERT_NEAR(10.000000000000000, mesh.Node(6).x, tolerance);
    ASSERT_NEAR(10.000000000000000, mesh.Node(7).x, tolerance);

    ASSERT_NEAR(0.00000000000000000, mesh.Node(0).y, tolerance);
    ASSERT_NEAR(0.00000000000000000, mesh.Node(1).y, tolerance);
    ASSERT_NEAR(0.00000000000000000, mesh.Node(2).y, tolerance);
    ASSERT_NEAR(0.00000000000000000, mesh.Node(3).y, tolerance);
    ASSERT_NEAR(-10.000000000000000, mesh.Node(4).y, tolerance);
    ASSERT_NEAR(-3.3333333333333330, mesh.Node(5).y, tolerance);
    ASSERT_NEAR(3.3333333333333344, mesh.Node(6).y, tolerance);
    ASSERT_NEAR(10.000000000000000, mesh.Node(7).y, tolerance);
}

TEST(Mesh1D, GenerateMeshFromPolyLines_WithFixedChainages_ShouldGenerateMesh)
{
    // 1 Setup
    std::vector<std::vector<meshkernel::Point>> polyLines{{{0.0, 0.0},
                                                           {10.0, 0.0},
                                                           {20.0, 0.0}}};

    std::vector<std::vector<double>> fixedChainages(polyLines.size());
    fixedChainages[0].emplace_back(8.0);
    fixedChainages[0].emplace_back(15.0);

    double const offset = 1.0;
    double const minFaceSize = 0.01;
    double const offsetFromFixedChainages = 0.2;

    meshkernel::Network1D Network1D(polyLines, meshkernel::Projection::cartesian);
    Network1D.ComputeFixedChainages(fixedChainages, minFaceSize, offsetFromFixedChainages);
    Network1D.ComputeOffsettedChainages(offset);

    // 2 Execution
    const auto mesh = meshkernel::Mesh1D(Network1D, minFaceSize);

    // 3 Assertion
    const auto tolerance = 1e-6;
    ASSERT_NEAR(0.0000000000000000, mesh.Node(0).x, tolerance);
    ASSERT_NEAR(0.9749999999999999, mesh.Node(1).x, tolerance);
    ASSERT_NEAR(1.9500000000000000, mesh.Node(2).x, tolerance);
    ASSERT_NEAR(2.9249999999999998, mesh.Node(3).x, tolerance);
    ASSERT_NEAR(3.8999999999999999, mesh.Node(4).x, tolerance);
    ASSERT_NEAR(4.8750000000000000, mesh.Node(5).x, tolerance);
    ASSERT_NEAR(5.8499999999999996, mesh.Node(6).x, tolerance);
    ASSERT_NEAR(6.8250000000000002, mesh.Node(7).x, tolerance);
    ASSERT_NEAR(7.7999999999999998, mesh.Node(8).x, tolerance);
    ASSERT_NEAR(8.1999999999999993, mesh.Node(9).x, tolerance);
    ASSERT_NEAR(9.1428571428571423, mesh.Node(10).x, tolerance);
    ASSERT_NEAR(10.085714285714285, mesh.Node(11).x, tolerance);
    ASSERT_NEAR(11.028571428571428, mesh.Node(12).x, tolerance);
    ASSERT_NEAR(11.971428571428572, mesh.Node(13).x, tolerance);
    ASSERT_NEAR(12.914285714285715, mesh.Node(14).x, tolerance);
    ASSERT_NEAR(13.857142857142858, mesh.Node(15).x, tolerance);
    ASSERT_NEAR(14.800000000000001, mesh.Node(16).x, tolerance);
    ASSERT_NEAR(15.199999999999999, mesh.Node(17).x, tolerance);
    ASSERT_NEAR(16.160000000000000, mesh.Node(18).x, tolerance);
    ASSERT_NEAR(17.120000000000001, mesh.Node(19).x, tolerance);
    ASSERT_NEAR(18.079999999999998, mesh.Node(20).x, tolerance);
    ASSERT_NEAR(19.039999999999999, mesh.Node(21).x, tolerance);
    ASSERT_NEAR(20.000000000000000, mesh.Node(22).x, tolerance);
}

TEST(Mesh1D, GenerateMeshFromPolyLines_WithChainagesWithinOffset_ShouldGenerateMesh)
{
    // 1 Setup
    std::vector<std::vector<meshkernel::Point>> polyLines{{{0.0, 0.0},
                                                           {10.0, 0.0},
                                                           {20.0, 0.0}}};

    std::vector<std::vector<double>> fixedChainages(polyLines.size());
    fixedChainages[0].emplace_back(8.0);
    fixedChainages[0].emplace_back(8.1);

    double const offset = 1.0;
    double const minFaceSize = 0.01;
    double const offsetFromFixedChainages = 0.2;

    meshkernel::Network1D Network1D(polyLines, meshkernel::Projection::cartesian);
    Network1D.ComputeFixedChainages(fixedChainages, minFaceSize, offsetFromFixedChainages);
    Network1D.ComputeOffsettedChainages(offset);

    // 2 Execution
    const auto mesh = meshkernel::Mesh1D(Network1D, minFaceSize);

    // 3 Assertion
    const auto tolerance = 1e-6;
    ASSERT_NEAR(0.00000000000000000, mesh.Node(0).x, tolerance);
    ASSERT_NEAR(0.97499999999999998, mesh.Node(1).x, tolerance);
    ASSERT_NEAR(1.9500000000000000, mesh.Node(2).x, tolerance);
    ASSERT_NEAR(2.9249999999999998, mesh.Node(3).x, tolerance);
    ASSERT_NEAR(3.8999999999999999, mesh.Node(4).x, tolerance);
    ASSERT_NEAR(4.8750000000000000, mesh.Node(5).x, tolerance);
    ASSERT_NEAR(5.8499999999999996, mesh.Node(6).x, tolerance);
    ASSERT_NEAR(6.8250000000000002, mesh.Node(7).x, tolerance);
    ASSERT_NEAR(7.7999999999999998, mesh.Node(8).x, tolerance);
    ASSERT_NEAR(8.0499999999999989, mesh.Node(9).x, tolerance);
    ASSERT_NEAR(8.2999999999999989, mesh.Node(10).x, tolerance);
    ASSERT_NEAR(9.2749999999999986, mesh.Node(11).x, tolerance);
    ASSERT_NEAR(10.250000000000000, mesh.Node(12).x, tolerance);
    ASSERT_NEAR(11.225000000000000, mesh.Node(13).x, tolerance);
    ASSERT_NEAR(12.199999999999999, mesh.Node(14).x, tolerance);
    ASSERT_NEAR(13.174999999999999, mesh.Node(15).x, tolerance);
    ASSERT_NEAR(14.149999999999999, mesh.Node(16).x, tolerance);
    ASSERT_NEAR(15.125000000000000, mesh.Node(17).x, tolerance);
    ASSERT_NEAR(16.100000000000001, mesh.Node(18).x, tolerance);
    ASSERT_NEAR(17.074999999999999, mesh.Node(19).x, tolerance);
    ASSERT_NEAR(18.049999999999997, mesh.Node(20).x, tolerance);
    ASSERT_NEAR(19.024999999999999, mesh.Node(21).x, tolerance);
    ASSERT_NEAR(20.000000000000000, mesh.Node(22).x, tolerance);
}
