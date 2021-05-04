#include <gtest/gtest.h>

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLineAttraction.hpp>
#include <TestUtils/MakeCurvilinearGrids.hpp>

TEST(CurvilinearLineAttraction, Compute_OnNLinesAttraction_ShouldAtractLines)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();
    meshkernel::CurvilinearGridLineAttraction curvilinearLineAttraction(curvilinearGrid, 0.5);
    curvilinearLineAttraction.SetLine({80266.8, 367104.0}, {80419.3, 366566.2});
    curvilinearLineAttraction.SetBlock({80198.2, 366750.6}, {80583.1, 366889.8});

    // Execute
    const auto shiftedGrid = curvilinearLineAttraction.Compute();

    // Asserts
    const double tolerance = 1e-6;
}