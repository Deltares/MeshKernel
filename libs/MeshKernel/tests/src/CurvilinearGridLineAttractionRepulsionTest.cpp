#include <gtest/gtest.h>

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLineAttractionRepulsion.hpp>
#include <TestUtils/MakeCurvilinearGrids.hpp>

TEST(CurvilinearLineAttraction, Compute_OnMLine_ShouldAttractMLines)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();
    meshkernel::CurvilinearGridLineAttractionRepulsion curvilinearLineAttractionRepulsion(curvilinearGrid, 0.5);
    curvilinearLineAttractionRepulsion.SetLine({80266.8, 367104.0}, {80419.3, 366566.2});
    curvilinearLineAttractionRepulsion.SetBlock({80198.2, 366750.6}, {80583.1, 366889.8});

    // Execute
    const auto modifiedGrid = curvilinearLineAttractionRepulsion.Compute();

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

TEST(CurvilinearLineAttraction, Compute_OnNLine_ShouldAttractNLines)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();
    meshkernel::CurvilinearGridLineAttractionRepulsion curvilinearLineAttractionRepulsion(curvilinearGrid, 0.5);
    curvilinearLineAttractionRepulsion.SetLine({80198.2, 366750.6}, {80583.1, 366889.8});
    curvilinearLineAttractionRepulsion.SetBlock({80266.8, 367104.0}, {80419.3, 366566.2});

    // Execute
    const auto modifiedGrid = curvilinearLineAttractionRepulsion.Compute();

    // Asserts
    const double tolerance = 1e-6;

    ASSERT_NEAR(80145.970831448722, modifiedGrid.m_gridNodes[0][2].x, tolerance);
    ASSERT_NEAR(80247.575117740766, modifiedGrid.m_gridNodes[1][2].x, tolerance);
    ASSERT_NEAR(80292.449108019704, modifiedGrid.m_gridNodes[2][2].x, tolerance);
    ASSERT_NEAR(80316.537053694148, modifiedGrid.m_gridNodes[3][2].x, tolerance);
    ASSERT_NEAR(80331.200564142913, modifiedGrid.m_gridNodes[4][2].x, tolerance);

    ASSERT_NEAR(367047.36276461056, modifiedGrid.m_gridNodes[0][2].y, tolerance);
    ASSERT_NEAR(366897.17707224732, modifiedGrid.m_gridNodes[1][2].y, tolerance);
    ASSERT_NEAR(366792.50812354451, modifiedGrid.m_gridNodes[2][2].y, tolerance);
    ASSERT_NEAR(366703.03776077798, modifiedGrid.m_gridNodes[3][2].y, tolerance);
    ASSERT_NEAR(366552.40947499714, modifiedGrid.m_gridNodes[4][2].y, tolerance);
}

TEST(CurvilinearLineRepulsion, Compute_OnMLine_ShouldRepulseMLines)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();
    meshkernel::CurvilinearGridLineAttractionRepulsion curvilinearLineAttractionRepulsion(curvilinearGrid, -0.5);
    curvilinearLineAttractionRepulsion.SetLine({80266.8, 367104.0}, {80419.3, 366566.2});
    curvilinearLineAttractionRepulsion.SetBlock({80198.2, 366750.6}, {80583.1, 366889.8});

    // Execute
    const auto modifiedGrid = curvilinearLineAttractionRepulsion.Compute();

    // Asserts
    const double tolerance = 1e-6;

    ASSERT_NEAR(80113.927180594226, modifiedGrid.m_gridNodes[0][2].x, tolerance);
    ASSERT_NEAR(80202.937709587452, modifiedGrid.m_gridNodes[1][2].x, tolerance);
    ASSERT_NEAR(80262.314053574693, modifiedGrid.m_gridNodes[2][2].x, tolerance);
    ASSERT_NEAR(80290.291966405988, modifiedGrid.m_gridNodes[3][2].x, tolerance);
    ASSERT_NEAR(80299.521456868417, modifiedGrid.m_gridNodes[4][2].x, tolerance);

    ASSERT_NEAR(367025.12442372262, modifiedGrid.m_gridNodes[0][2].y, tolerance);
    ASSERT_NEAR(366900.80385695840, modifiedGrid.m_gridNodes[1][2].y, tolerance);
    ASSERT_NEAR(366781.77878604224, modifiedGrid.m_gridNodes[2][2].y, tolerance);
    ASSERT_NEAR(366674.26710998092, modifiedGrid.m_gridNodes[3][2].y, tolerance);
    ASSERT_NEAR(366549.70842920581, modifiedGrid.m_gridNodes[4][2].y, tolerance);
}

TEST(CurvilinearLineRepulsion, Compute_OnNLine_ShouldRepulseNLines)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();
    meshkernel::CurvilinearGridLineAttractionRepulsion curvilinearLineAttractionRepulsion(curvilinearGrid, -0.5);
    curvilinearLineAttractionRepulsion.SetLine({80198.2, 366750.6}, {80583.1, 366889.8});
    curvilinearLineAttractionRepulsion.SetBlock({80266.8, 367104.0}, {80419.3, 366566.2});

    // Execute
    const auto modifiedGrid = curvilinearLineAttractionRepulsion.Compute();

    // Asserts
    const double tolerance = 1e-6;

    ASSERT_NEAR(80145.970831448722, modifiedGrid.m_gridNodes[0][2].x, tolerance);
    ASSERT_NEAR(80222.273272260049, modifiedGrid.m_gridNodes[1][2].x, tolerance);
    ASSERT_NEAR(80292.449108019704, modifiedGrid.m_gridNodes[2][2].x, tolerance);
    ASSERT_NEAR(80324.255708261146, modifiedGrid.m_gridNodes[3][2].x, tolerance);
    ASSERT_NEAR(80331.200564142913, modifiedGrid.m_gridNodes[4][2].x, tolerance);

    ASSERT_NEAR(367047.36276461056, modifiedGrid.m_gridNodes[0][2].y, tolerance);
    ASSERT_NEAR(366941.19925013784, modifiedGrid.m_gridNodes[1][2].y, tolerance);
    ASSERT_NEAR(366792.50812354451, modifiedGrid.m_gridNodes[2][2].y, tolerance);
    ASSERT_NEAR(366655.21404741227, modifiedGrid.m_gridNodes[3][2].y, tolerance);
    ASSERT_NEAR(366552.40947499714, modifiedGrid.m_gridNodes[4][2].y, tolerance);
}
