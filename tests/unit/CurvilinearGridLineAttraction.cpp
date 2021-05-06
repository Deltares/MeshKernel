#include <gtest/gtest.h>

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLineAttraction.hpp>
#include <TestUtils/MakeCurvilinearGrids.hpp>

TEST(CurvilinearLineAttraction, Compute_OnMLine_ShouldAttractMLines)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();
    meshkernel::CurvilinearGridLineAttraction curvilinearLineAttraction(curvilinearGrid, 0.5);
    curvilinearLineAttraction.SetLine({80266.8, 367104.0}, {80419.3, 366566.2});
    curvilinearLineAttraction.SetBlock({80198.2, 366750.6}, {80583.1, 366889.8});

    // Execute
    const auto modifiedGrid = curvilinearLineAttraction.Compute();

    // Asserts
    const double tolerance = 1e-6;

    ASSERT_NEAR(80178.014482303217, modifiedGrid.m_gridNodes[0][2].x, tolerance);
    ASSERT_NEAR(80266.910680413363, modifiedGrid.m_gridNodes[1][2].x, tolerance);
    ASSERT_NEAR(80322.584162464715, modifiedGrid.m_gridNodes[2][2].x, tolerance);
    ASSERT_NEAR(80350.500795549306, modifiedGrid.m_gridNodes[3][2].x, tolerance);
    ASSERT_NEAR(80362.879671417410, modifiedGrid.m_gridNodes[4][2].x, tolerance);

    ASSERT_NEAR(367069.60110549850, modifiedGrid.m_gridNodes[0][2].y, tolerance);
    ASSERT_NEAR(366937.57246542675, modifiedGrid.m_gridNodes[1][2].y, tolerance);
    ASSERT_NEAR(366803.23746104678, modifiedGrid.m_gridNodes[2][2].y, tolerance);
    ASSERT_NEAR(366683.98469820933, modifiedGrid.m_gridNodes[3][2].y, tolerance);
    ASSERT_NEAR(366555.11052078847, modifiedGrid.m_gridNodes[4][2].y, tolerance);

}