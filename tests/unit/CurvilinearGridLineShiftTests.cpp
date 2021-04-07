#include <gtest/gtest.h>

#include <MeshKernel/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGridLineShift.hpp>
#include <TestUtils/MakeCurvilinearGrids.hpp>

TEST(CurvilinearLineShift, Compute_OnMGridline_ShouldShiftLine)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();
    meshkernel::CurvilinearGridLineShift curvilinearLineShift(curvilinearGrid);
    curvilinearLineShift.SetLine({79982, 366934}, {80155, 366530});
    curvilinearLineShift.SetBlock({80108, 366707}, {80291, 366792});
    curvilinearLineShift.MoveNode({79982, 366934}, {79872, 366876});

    // Execute
    const auto smoothedGrid = curvilinearLineShift.Compute();
}