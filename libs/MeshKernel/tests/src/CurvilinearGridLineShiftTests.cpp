#include <gtest/gtest.h>

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLineShift.hpp>
#include <TestUtils/MakeCurvilinearGrids.hpp>

TEST(CurvilinearLineShift, Compute_OnMGridlineShiftingOneNode_ShouldShiftLine)
{
    // Set-up
    auto curvilinearGrid = MakeSmallCurvilinearGrid();
    meshkernel::CurvilinearGridLineShift curvilinearLineShift(curvilinearGrid);
    curvilinearLineShift.SetLine({79982.0, 366934.0}, {80155.0, 366530.0});
    curvilinearLineShift.SetBlock({80108.0, 366707.0}, {80291.0, 366792.0});
    curvilinearLineShift.MoveNode({79982.0, 366934.0}, {79872.0, 366876.0});

    // Execute
    curvilinearLineShift.Compute();

    // Asserts
    const double tolerance = 1e-6;

    ASSERT_NEAR(79872.000000000000, curvilinearGrid.m_gridNodes(0, 0).x, tolerance);
    ASSERT_NEAR(80010.039799507853, curvilinearGrid.m_gridNodes(0, 1).x, tolerance);
    ASSERT_NEAR(80145.970831448722, curvilinearGrid.m_gridNodes(0, 2).x, tolerance);
    ASSERT_NEAR(80225.900042018140, curvilinearGrid.m_gridNodes(0, 3).x, tolerance);
    ASSERT_NEAR(80305.243756829266, curvilinearGrid.m_gridNodes(0, 4).x, tolerance);
    ASSERT_NEAR(80381.747982750283, curvilinearGrid.m_gridNodes(0, 5).x, tolerance);
    ASSERT_NEAR(80458.252208671300, curvilinearGrid.m_gridNodes(0, 6).x, tolerance);
    ASSERT_NEAR(80534.756434592317, curvilinearGrid.m_gridNodes(0, 7).x, tolerance);
    ASSERT_NEAR(80611.260660513333, curvilinearGrid.m_gridNodes(0, 8).x, tolerance);

    ASSERT_NEAR(79970.149644452977, curvilinearGrid.m_gridNodes(1, 0).x, tolerance);
    ASSERT_NEAR(80103.062377666603, curvilinearGrid.m_gridNodes(1, 1).x, tolerance);
    ASSERT_NEAR(80234.924195000407, curvilinearGrid.m_gridNodes(1, 2).x, tolerance);
    ASSERT_NEAR(80324.671765221428, curvilinearGrid.m_gridNodes(1, 3).x, tolerance);
    ASSERT_NEAR(80414.057391982613, curvilinearGrid.m_gridNodes(1, 4).x, tolerance);
    ASSERT_NEAR(80505.096476712482, curvilinearGrid.m_gridNodes(1, 5).x, tolerance);
    ASSERT_NEAR(80595.183339827883, curvilinearGrid.m_gridNodes(1, 6).x, tolerance);
    ASSERT_NEAR(80684.333994102650, curvilinearGrid.m_gridNodes(1, 7).x, tolerance);
    ASSERT_NEAR(80772.567299473958, curvilinearGrid.m_gridNodes(1, 8).x, tolerance);

    ASSERT_NEAR(366876.00000000000, curvilinearGrid.m_gridNodes(0, 0).y, tolerance);
    ASSERT_NEAR(366959.82623907487, curvilinearGrid.m_gridNodes(0, 1).y, tolerance);
    ASSERT_NEAR(367047.36276461056, curvilinearGrid.m_gridNodes(0, 2).y, tolerance);
    ASSERT_NEAR(367104.62934968271, curvilinearGrid.m_gridNodes(0, 3).y, tolerance);
    ASSERT_NEAR(367163.01691965276, curvilinearGrid.m_gridNodes(0, 4).y, tolerance);
    ASSERT_NEAR(367224.10904462705, curvilinearGrid.m_gridNodes(0, 5).y, tolerance);
    ASSERT_NEAR(367285.20116960135, curvilinearGrid.m_gridNodes(0, 6).y, tolerance);
    ASSERT_NEAR(367346.29329457565, curvilinearGrid.m_gridNodes(0, 7).y, tolerance);
    ASSERT_NEAR(367407.38541954994, curvilinearGrid.m_gridNodes(0, 8).y, tolerance);

    ASSERT_NEAR(366781.50715811126, curvilinearGrid.m_gridNodes(1, 0).y, tolerance);
    ASSERT_NEAR(366849.28921837400, curvilinearGrid.m_gridNodes(1, 1).y, tolerance);
    ASSERT_NEAR(366919.18816119258, curvilinearGrid.m_gridNodes(1, 2).y, tolerance);
    ASSERT_NEAR(366966.76979594346, curvilinearGrid.m_gridNodes(1, 3).y, tolerance);
    ASSERT_NEAR(367015.14849423966, curvilinearGrid.m_gridNodes(1, 4).y, tolerance);
    ASSERT_NEAR(367056.48898898275, curvilinearGrid.m_gridNodes(1, 5).y, tolerance);
    ASSERT_NEAR(367099.12347147451, curvilinearGrid.m_gridNodes(1, 6).y, tolerance);
    ASSERT_NEAR(367143.03018172452, curvilinearGrid.m_gridNodes(1, 7).y, tolerance);
    ASSERT_NEAR(367188.18349069095, curvilinearGrid.m_gridNodes(1, 8).y, tolerance);
}

TEST(CurvilinearLineShift, Compute_OnMGridlineShiftingTwoNodes_ShouldShiftLine)
{
    // Set-up
    auto curvilinearGrid = MakeSmallCurvilinearGrid();
    meshkernel::CurvilinearGridLineShift curvilinearLineShift(curvilinearGrid);
    curvilinearLineShift.SetLine({79982.0, 366934.0}, {80155.0, 366530.0});
    curvilinearLineShift.SetBlock({80108.0, 366707.0}, {80291.0, 366792.0});

    // Move two nodes
    curvilinearLineShift.MoveNode({79982.0, 366934.0}, {79872.0, 366876.0});
    curvilinearLineShift.MoveNode({80053.0, 366823.0}, {79932.0, 366773.0});

    // Execute
    curvilinearLineShift.Compute();

    // Asserts
    const double tolerance = 1e-6;

    ASSERT_NEAR(79872.000000000000, curvilinearGrid.m_gridNodes(0, 0).x, tolerance);
    ASSERT_NEAR(80010.039799507853, curvilinearGrid.m_gridNodes(0, 1).x, tolerance);
    ASSERT_NEAR(80145.970831448722, curvilinearGrid.m_gridNodes(0, 2).x, tolerance);
    ASSERT_NEAR(80225.900042018140, curvilinearGrid.m_gridNodes(0, 3).x, tolerance);
    ASSERT_NEAR(80305.243756829266, curvilinearGrid.m_gridNodes(0, 4).x, tolerance);
    ASSERT_NEAR(80381.747982750283, curvilinearGrid.m_gridNodes(0, 5).x, tolerance);
    ASSERT_NEAR(80458.252208671300, curvilinearGrid.m_gridNodes(0, 6).x, tolerance);
    ASSERT_NEAR(80534.756434592317, curvilinearGrid.m_gridNodes(0, 7).x, tolerance);
    ASSERT_NEAR(80611.260660513333, curvilinearGrid.m_gridNodes(0, 8).x, tolerance);

    ASSERT_NEAR(79932.000000000000, curvilinearGrid.m_gridNodes(1, 0).x, tolerance);
    ASSERT_NEAR(80084.035373914361, curvilinearGrid.m_gridNodes(1, 1).x, tolerance);
    ASSERT_NEAR(80234.924195000407, curvilinearGrid.m_gridNodes(1, 2).x, tolerance);
    ASSERT_NEAR(80324.671765221428, curvilinearGrid.m_gridNodes(1, 3).x, tolerance);
    ASSERT_NEAR(80414.057391982613, curvilinearGrid.m_gridNodes(1, 4).x, tolerance);
    ASSERT_NEAR(80505.096476712482, curvilinearGrid.m_gridNodes(1, 5).x, tolerance);
    ASSERT_NEAR(80595.183339827883, curvilinearGrid.m_gridNodes(1, 6).x, tolerance);
    ASSERT_NEAR(80684.333994102650, curvilinearGrid.m_gridNodes(1, 7).x, tolerance);
    ASSERT_NEAR(80772.567299473958, curvilinearGrid.m_gridNodes(1, 8).x, tolerance);

    ASSERT_NEAR(366876.00000000000, curvilinearGrid.m_gridNodes(0, 0).y, tolerance);
    ASSERT_NEAR(366959.82623907487, curvilinearGrid.m_gridNodes(0, 1).y, tolerance);
    ASSERT_NEAR(367047.36276461056, curvilinearGrid.m_gridNodes(0, 2).y, tolerance);
    ASSERT_NEAR(367104.62934968271, curvilinearGrid.m_gridNodes(0, 3).y, tolerance);
    ASSERT_NEAR(367163.01691965276, curvilinearGrid.m_gridNodes(0, 4).y, tolerance);
    ASSERT_NEAR(367224.10904462705, curvilinearGrid.m_gridNodes(0, 5).y, tolerance);
    ASSERT_NEAR(367285.20116960135, curvilinearGrid.m_gridNodes(0, 6).y, tolerance);
    ASSERT_NEAR(367346.29329457565, curvilinearGrid.m_gridNodes(0, 7).y, tolerance);
    ASSERT_NEAR(367407.38541954994, curvilinearGrid.m_gridNodes(0, 8).y, tolerance);

    ASSERT_NEAR(366773.00000000000, curvilinearGrid.m_gridNodes(1, 0).y, tolerance);
    ASSERT_NEAR(366844.82660636236, curvilinearGrid.m_gridNodes(1, 1).y, tolerance);
    ASSERT_NEAR(366919.18816119258, curvilinearGrid.m_gridNodes(1, 2).y, tolerance);
    ASSERT_NEAR(366966.76979594346, curvilinearGrid.m_gridNodes(1, 3).y, tolerance);
    ASSERT_NEAR(367015.14849423966, curvilinearGrid.m_gridNodes(1, 4).y, tolerance);
    ASSERT_NEAR(367056.48898898275, curvilinearGrid.m_gridNodes(1, 5).y, tolerance);
    ASSERT_NEAR(367099.12347147451, curvilinearGrid.m_gridNodes(1, 6).y, tolerance);
    ASSERT_NEAR(367143.03018172452, curvilinearGrid.m_gridNodes(1, 7).y, tolerance);
    ASSERT_NEAR(367188.18349069095, curvilinearGrid.m_gridNodes(1, 8).y, tolerance);
}
