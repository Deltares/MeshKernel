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
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLineShift.hpp>
#include <TestUtils/MakeCurvilinearGrids.hpp>

TEST(CurvilinearLineShift, Compute_OnMGridlineShiftingOneNode_ShouldShiftLine)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();
    meshkernel::CurvilinearGridLineShift curvilinearLineShift(*curvilinearGrid);
    curvilinearLineShift.SetLine({79982.0, 366934.0}, {80155.0, 366530.0});
    curvilinearLineShift.SetBlock(meshkernel::Point{80108.0, 366707.0}, meshkernel::Point{80291.0, 366792.0});
    [[maybe_unused]] auto dummyUndoAction1 = curvilinearLineShift.MoveNode({79982.0, 366934.0}, {79872.0, 366876.0});

    // Execute
    [[maybe_unused]] auto dummyUndoAction2 = curvilinearLineShift.Compute();

    // Asserts
    const double tolerance = 1e-6;

    ASSERT_NEAR(79872.000000000000, curvilinearGrid->GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(80010.039799507853, curvilinearGrid->GetNode(0, 1).x, tolerance);
    ASSERT_NEAR(80145.970831448722, curvilinearGrid->GetNode(0, 2).x, tolerance);
    ASSERT_NEAR(80225.900042018140, curvilinearGrid->GetNode(0, 3).x, tolerance);
    ASSERT_NEAR(80305.243756829266, curvilinearGrid->GetNode(0, 4).x, tolerance);
    ASSERT_NEAR(80381.747982750283, curvilinearGrid->GetNode(0, 5).x, tolerance);
    ASSERT_NEAR(80458.252208671300, curvilinearGrid->GetNode(0, 6).x, tolerance);
    ASSERT_NEAR(80534.756434592317, curvilinearGrid->GetNode(0, 7).x, tolerance);
    ASSERT_NEAR(80611.260660513333, curvilinearGrid->GetNode(0, 8).x, tolerance);

    ASSERT_NEAR(79970.149644452977, curvilinearGrid->GetNode(1, 0).x, tolerance);
    ASSERT_NEAR(80103.062377666603, curvilinearGrid->GetNode(1, 1).x, tolerance);
    ASSERT_NEAR(80234.924195000407, curvilinearGrid->GetNode(1, 2).x, tolerance);
    ASSERT_NEAR(80324.671765221428, curvilinearGrid->GetNode(1, 3).x, tolerance);
    ASSERT_NEAR(80414.057391982613, curvilinearGrid->GetNode(1, 4).x, tolerance);
    ASSERT_NEAR(80505.096476712482, curvilinearGrid->GetNode(1, 5).x, tolerance);
    ASSERT_NEAR(80595.183339827883, curvilinearGrid->GetNode(1, 6).x, tolerance);
    ASSERT_NEAR(80684.333994102650, curvilinearGrid->GetNode(1, 7).x, tolerance);
    ASSERT_NEAR(80772.567299473958, curvilinearGrid->GetNode(1, 8).x, tolerance);

    ASSERT_NEAR(366876.00000000000, curvilinearGrid->GetNode(0, 0).y, tolerance);
    ASSERT_NEAR(366959.82623907487, curvilinearGrid->GetNode(0, 1).y, tolerance);
    ASSERT_NEAR(367047.36276461056, curvilinearGrid->GetNode(0, 2).y, tolerance);
    ASSERT_NEAR(367104.62934968271, curvilinearGrid->GetNode(0, 3).y, tolerance);
    ASSERT_NEAR(367163.01691965276, curvilinearGrid->GetNode(0, 4).y, tolerance);
    ASSERT_NEAR(367224.10904462705, curvilinearGrid->GetNode(0, 5).y, tolerance);
    ASSERT_NEAR(367285.20116960135, curvilinearGrid->GetNode(0, 6).y, tolerance);
    ASSERT_NEAR(367346.29329457565, curvilinearGrid->GetNode(0, 7).y, tolerance);
    ASSERT_NEAR(367407.38541954994, curvilinearGrid->GetNode(0, 8).y, tolerance);

    ASSERT_NEAR(366781.50715811126, curvilinearGrid->GetNode(1, 0).y, tolerance);
    ASSERT_NEAR(366849.28921837400, curvilinearGrid->GetNode(1, 1).y, tolerance);
    ASSERT_NEAR(366919.18816119258, curvilinearGrid->GetNode(1, 2).y, tolerance);
    ASSERT_NEAR(366966.76979594346, curvilinearGrid->GetNode(1, 3).y, tolerance);
    ASSERT_NEAR(367015.14849423966, curvilinearGrid->GetNode(1, 4).y, tolerance);
    ASSERT_NEAR(367056.48898898275, curvilinearGrid->GetNode(1, 5).y, tolerance);
    ASSERT_NEAR(367099.12347147451, curvilinearGrid->GetNode(1, 6).y, tolerance);
    ASSERT_NEAR(367143.03018172452, curvilinearGrid->GetNode(1, 7).y, tolerance);
    ASSERT_NEAR(367188.18349069095, curvilinearGrid->GetNode(1, 8).y, tolerance);
}

TEST(CurvilinearLineShift, Compute_OnMGridlineShiftingTwoNodes_ShouldShiftLine)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();
    meshkernel::CurvilinearGridLineShift curvilinearLineShift(*curvilinearGrid);
    curvilinearLineShift.SetLine({79982.0, 366934.0}, {80155.0, 366530.0});
    curvilinearLineShift.SetBlock(meshkernel::Point{80108.0, 366707.0}, meshkernel::Point{80291.0, 366792.0});

    // Move two nodes
    [[maybe_unused]] auto dummyUndoAction1 = curvilinearLineShift.MoveNode({79982.0, 366934.0}, {79872.0, 366876.0});
    [[maybe_unused]] auto dummyUndoAction2 = curvilinearLineShift.MoveNode({80053.0, 366823.0}, {79932.0, 366773.0});

    // Execute
    [[maybe_unused]] auto dummyUndoAction3 = curvilinearLineShift.Compute();

    // Asserts
    const double tolerance = 1e-6;

    ASSERT_NEAR(79872.000000000000, curvilinearGrid->GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(80010.039799507853, curvilinearGrid->GetNode(0, 1).x, tolerance);
    ASSERT_NEAR(80145.970831448722, curvilinearGrid->GetNode(0, 2).x, tolerance);
    ASSERT_NEAR(80225.900042018140, curvilinearGrid->GetNode(0, 3).x, tolerance);
    ASSERT_NEAR(80305.243756829266, curvilinearGrid->GetNode(0, 4).x, tolerance);
    ASSERT_NEAR(80381.747982750283, curvilinearGrid->GetNode(0, 5).x, tolerance);
    ASSERT_NEAR(80458.252208671300, curvilinearGrid->GetNode(0, 6).x, tolerance);
    ASSERT_NEAR(80534.756434592317, curvilinearGrid->GetNode(0, 7).x, tolerance);
    ASSERT_NEAR(80611.260660513333, curvilinearGrid->GetNode(0, 8).x, tolerance);

    ASSERT_NEAR(79932.000000000000, curvilinearGrid->GetNode(1, 0).x, tolerance);
    ASSERT_NEAR(80084.035373914361, curvilinearGrid->GetNode(1, 1).x, tolerance);
    ASSERT_NEAR(80234.924195000407, curvilinearGrid->GetNode(1, 2).x, tolerance);
    ASSERT_NEAR(80324.671765221428, curvilinearGrid->GetNode(1, 3).x, tolerance);
    ASSERT_NEAR(80414.057391982613, curvilinearGrid->GetNode(1, 4).x, tolerance);
    ASSERT_NEAR(80505.096476712482, curvilinearGrid->GetNode(1, 5).x, tolerance);
    ASSERT_NEAR(80595.183339827883, curvilinearGrid->GetNode(1, 6).x, tolerance);
    ASSERT_NEAR(80684.333994102650, curvilinearGrid->GetNode(1, 7).x, tolerance);
    ASSERT_NEAR(80772.567299473958, curvilinearGrid->GetNode(1, 8).x, tolerance);

    ASSERT_NEAR(366876.00000000000, curvilinearGrid->GetNode(0, 0).y, tolerance);
    ASSERT_NEAR(366959.82623907487, curvilinearGrid->GetNode(0, 1).y, tolerance);
    ASSERT_NEAR(367047.36276461056, curvilinearGrid->GetNode(0, 2).y, tolerance);
    ASSERT_NEAR(367104.62934968271, curvilinearGrid->GetNode(0, 3).y, tolerance);
    ASSERT_NEAR(367163.01691965276, curvilinearGrid->GetNode(0, 4).y, tolerance);
    ASSERT_NEAR(367224.10904462705, curvilinearGrid->GetNode(0, 5).y, tolerance);
    ASSERT_NEAR(367285.20116960135, curvilinearGrid->GetNode(0, 6).y, tolerance);
    ASSERT_NEAR(367346.29329457565, curvilinearGrid->GetNode(0, 7).y, tolerance);
    ASSERT_NEAR(367407.38541954994, curvilinearGrid->GetNode(0, 8).y, tolerance);

    ASSERT_NEAR(366773.00000000000, curvilinearGrid->GetNode(1, 0).y, tolerance);
    ASSERT_NEAR(366844.82660636236, curvilinearGrid->GetNode(1, 1).y, tolerance);
    ASSERT_NEAR(366919.18816119258, curvilinearGrid->GetNode(1, 2).y, tolerance);
    ASSERT_NEAR(366966.76979594346, curvilinearGrid->GetNode(1, 3).y, tolerance);
    ASSERT_NEAR(367015.14849423966, curvilinearGrid->GetNode(1, 4).y, tolerance);
    ASSERT_NEAR(367056.48898898275, curvilinearGrid->GetNode(1, 5).y, tolerance);
    ASSERT_NEAR(367099.12347147451, curvilinearGrid->GetNode(1, 6).y, tolerance);
    ASSERT_NEAR(367143.03018172452, curvilinearGrid->GetNode(1, 7).y, tolerance);
    ASSERT_NEAR(367188.18349069095, curvilinearGrid->GetNode(1, 8).y, tolerance);
}
