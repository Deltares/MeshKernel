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
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLineAttractionRepulsion.hpp>
#include <TestUtils/MakeCurvilinearGrids.hpp>

TEST(CurvilinearLineAttraction, Compute_OnMLine_ShouldAttractMLines)
{
    // Set-up
    auto curvilinearGrid = MakeSmallCurvilinearGrid();
    meshkernel::CurvilinearGridLineAttractionRepulsion curvilinearLineAttractionRepulsion(*curvilinearGrid, 0.5);
    curvilinearLineAttractionRepulsion.SetLine({80266.8, 367104.0}, {80419.3, 366566.2});
    curvilinearLineAttractionRepulsion.SetBlock(meshkernel::Point{80198.2, 366750.6},
                                                meshkernel::Point{80583.1, 366889.8});

    // Execute
    [[maybe_unused]] auto dummyUndoAction = curvilinearLineAttractionRepulsion.Compute();

    // Asserts
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(80113.927180594226, curvilinearGrid->GetNode(0, 2).x, tolerance);
    ASSERT_NEAR(80202.937709587452, curvilinearGrid->GetNode(1, 2).x, tolerance);
    ASSERT_NEAR(80262.314053574693, curvilinearGrid->GetNode(2, 2).x, tolerance);
    ASSERT_NEAR(80290.291966405988, curvilinearGrid->GetNode(3, 2).x, tolerance);
    ASSERT_NEAR(80299.521456868417, curvilinearGrid->GetNode(4, 2).x, tolerance);

    ASSERT_NEAR(367025.12442372262, curvilinearGrid->GetNode(0, 2).y, tolerance);
    ASSERT_NEAR(366900.80385695840, curvilinearGrid->GetNode(1, 2).y, tolerance);
    ASSERT_NEAR(366781.77878604224, curvilinearGrid->GetNode(2, 2).y, tolerance);
    ASSERT_NEAR(366674.26710998092, curvilinearGrid->GetNode(3, 2).y, tolerance);
    ASSERT_NEAR(366549.70842920581, curvilinearGrid->GetNode(4, 2).y, tolerance);
}

TEST(CurvilinearLineAttraction, Compute_OnNLine_ShouldAttractNLines)
{
    // Set-up
    auto curvilinearGrid = MakeSmallCurvilinearGrid();
    meshkernel::CurvilinearGridLineAttractionRepulsion curvilinearLineAttractionRepulsion(*curvilinearGrid, 0.5);
    curvilinearLineAttractionRepulsion.SetLine({80198.2, 366750.6}, {80583.1, 366889.8});
    curvilinearLineAttractionRepulsion.SetBlock(meshkernel::Point{80266.8, 367104.0},
                                                meshkernel::Point{80419.3, 366566.2});

    // Execute
    [[maybe_unused]] auto dummyUndoAction = curvilinearLineAttractionRepulsion.Compute();

    // Asserts
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(80145.970831448722, curvilinearGrid->GetNode(0, 2).x, tolerance);
    ASSERT_NEAR(80247.575117740766, curvilinearGrid->GetNode(1, 2).x, tolerance);
    ASSERT_NEAR(80292.449108019704, curvilinearGrid->GetNode(2, 2).x, tolerance);
    ASSERT_NEAR(80316.537053694148, curvilinearGrid->GetNode(3, 2).x, tolerance);
    ASSERT_NEAR(80331.200564142913, curvilinearGrid->GetNode(4, 2).x, tolerance);

    ASSERT_NEAR(367047.36276461056, curvilinearGrid->GetNode(0, 2).y, tolerance);
    ASSERT_NEAR(366897.17707224732, curvilinearGrid->GetNode(1, 2).y, tolerance);
    ASSERT_NEAR(366792.50812354451, curvilinearGrid->GetNode(2, 2).y, tolerance);
    ASSERT_NEAR(366703.03776077798, curvilinearGrid->GetNode(3, 2).y, tolerance);
    ASSERT_NEAR(366552.40947499714, curvilinearGrid->GetNode(4, 2).y, tolerance);
}

TEST(CurvilinearLineRepulsion, Compute_OnMLine_ShouldRepulseMLines)
{
    // Set-up
    auto curvilinearGrid = MakeSmallCurvilinearGrid();
    meshkernel::CurvilinearGridLineAttractionRepulsion curvilinearLineAttractionRepulsion(*curvilinearGrid, -0.5);
    curvilinearLineAttractionRepulsion.SetLine({80266.8, 367104.0}, {80419.3, 366566.2});
    curvilinearLineAttractionRepulsion.SetBlock(meshkernel::Point{80198.2, 366750.6}, meshkernel::Point{80583.1, 366889.8});

    // Execute
    [[maybe_unused]] auto dummyUndoAction = curvilinearLineAttractionRepulsion.Compute();

    // Asserts
    const double tolerance = 1e-6;

    ASSERT_NEAR(80178.014482303217, curvilinearGrid->GetNode(0, 2).x, tolerance);
    ASSERT_NEAR(80266.910680413363, curvilinearGrid->GetNode(1, 2).x, tolerance);
    ASSERT_NEAR(80322.584162464715, curvilinearGrid->GetNode(2, 2).x, tolerance);
    ASSERT_NEAR(80350.500795549306, curvilinearGrid->GetNode(3, 2).x, tolerance);
    ASSERT_NEAR(80362.879671417410, curvilinearGrid->GetNode(4, 2).x, tolerance);

    ASSERT_NEAR(367069.60110549850, curvilinearGrid->GetNode(0, 2).y, tolerance);
    ASSERT_NEAR(366937.57246542675, curvilinearGrid->GetNode(1, 2).y, tolerance);
    ASSERT_NEAR(366803.23746104678, curvilinearGrid->GetNode(2, 2).y, tolerance);
    ASSERT_NEAR(366683.98469820933, curvilinearGrid->GetNode(3, 2).y, tolerance);
    ASSERT_NEAR(366555.11052078847, curvilinearGrid->GetNode(4, 2).y, tolerance);
}

TEST(CurvilinearLineRepulsion, Compute_OnNLine_ShouldRepulseNLines)
{
    // Set-up
    auto curvilinearGrid = MakeSmallCurvilinearGrid();
    meshkernel::CurvilinearGridLineAttractionRepulsion curvilinearLineAttractionRepulsion(*curvilinearGrid, -0.5);
    curvilinearLineAttractionRepulsion.SetLine({80198.2, 366750.6}, {80583.1, 366889.8});
    curvilinearLineAttractionRepulsion.SetBlock(meshkernel::Point{80266.8, 367104.0}, meshkernel::Point{80419.3, 366566.2});

    // Execute
    [[maybe_unused]] auto dummyUndoAction = curvilinearLineAttractionRepulsion.Compute();

    // Asserts
    const double tolerance = 1e-6;

    ASSERT_NEAR(80145.970831448722, curvilinearGrid->GetNode(0, 2).x, tolerance);
    ASSERT_NEAR(80222.273272260049, curvilinearGrid->GetNode(1, 2).x, tolerance);
    ASSERT_NEAR(80292.449108019704, curvilinearGrid->GetNode(2, 2).x, tolerance);
    ASSERT_NEAR(80324.255708261146, curvilinearGrid->GetNode(3, 2).x, tolerance);
    ASSERT_NEAR(80331.200564142913, curvilinearGrid->GetNode(4, 2).x, tolerance);

    ASSERT_NEAR(367047.36276461056, curvilinearGrid->GetNode(0, 2).y, tolerance);
    ASSERT_NEAR(366941.19925013784, curvilinearGrid->GetNode(1, 2).y, tolerance);
    ASSERT_NEAR(366792.50812354451, curvilinearGrid->GetNode(2, 2).y, tolerance);
    ASSERT_NEAR(366655.21404741227, curvilinearGrid->GetNode(3, 2).y, tolerance);
    ASSERT_NEAR(366552.40947499714, curvilinearGrid->GetNode(4, 2).y, tolerance);
}
