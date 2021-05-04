#include <gtest/gtest.h>

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLineAttraction.hpp>
#include <TestUtils/MakeCurvilinearGrids.hpp>

TEST(CurvilinearLineAttraction, Compute_OnNLinesAttraction_ShouldAtractLines)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();
    meshkernel::CurvilinearGridLineAttraction curvilinearLineAttraction(curvilinearGrid, 0.5);
    curvilinearLineAttraction.SetLine({79982.0, 366934.0}, {80155.0, 366530.0});
    curvilinearLineAttraction.SetBlock({80108.0, 366707.0}, {80291.0, 366792.0});

    // Execute
    const auto shiftedGrid = curvilinearLineAttraction.Compute();

    // Asserts
    const double tolerance = 1e-6;
}