#include <gtest/gtest.h>

#include <MeshKernel/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGridOrthogonalization.hpp>
#include <MeshKernel/Entities.hpp>
#include <TestUtils/MakeCurvilinearGrids.hpp>

TEST(CurvilinearGridOrthogonalization, Compute_OnOrthogonalCurvilinearGrid_ShouldNotModifyGrid)
{
    // Set-up
    std::vector<std::vector<meshkernel::Point>> grid{
        {{0, 0}, {0, 10}, {0, 20}, {0, 30}},
        {{10, 0}, {10, 10}, {10, 20}, {10, 30}},
        {{20, 0}, {20, 10}, {20, 20}, {20, 30}},
        {{30, 0}, {30, 10}, {30, 20}, {30, 30}}};

    const auto curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(std::move(grid), meshkernel::Projection::cartesian);

    meshkernelapi::OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.OuterIterations = 1;
    orthogonalizationParameters.BoundaryIterations = 25;
    orthogonalizationParameters.InnerIterations = 25;
    orthogonalizationParameters.OrthogonalizationToSmoothingFactor = 0.975;
    meshkernel::CurvilinearGridOrthogonalization curvilinearGridOrthogonalization(curvilinearGrid, orthogonalizationParameters, {0, 0}, {30, 30});

    // Execute
    curvilinearGridOrthogonalization.Compute();

    // Assert nodes are on the same location because the grid is already orthogonal
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(0.0, curvilinearGrid->m_gridNodes[0][0].x, tolerance);
    ASSERT_NEAR(0.0, curvilinearGrid->m_gridNodes[0][1].x, tolerance);
    ASSERT_NEAR(0.0, curvilinearGrid->m_gridNodes[0][2].x, tolerance);
    ASSERT_NEAR(0.0, curvilinearGrid->m_gridNodes[0][3].x, tolerance);

    ASSERT_NEAR(10.0, curvilinearGrid->m_gridNodes[1][0].x, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid->m_gridNodes[1][1].x, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid->m_gridNodes[1][2].x, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid->m_gridNodes[1][3].x, tolerance);

    ASSERT_NEAR(20.0, curvilinearGrid->m_gridNodes[2][0].x, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid->m_gridNodes[2][1].x, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid->m_gridNodes[2][2].x, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid->m_gridNodes[2][3].x, tolerance);

    ASSERT_NEAR(30.0, curvilinearGrid->m_gridNodes[3][0].x, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid->m_gridNodes[3][1].x, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid->m_gridNodes[3][2].x, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid->m_gridNodes[3][3].x, tolerance);

    ASSERT_NEAR(0.0, curvilinearGrid->m_gridNodes[0][0].y, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid->m_gridNodes[0][1].y, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid->m_gridNodes[0][2].y, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid->m_gridNodes[0][3].y, tolerance);

    ASSERT_NEAR(0.0, curvilinearGrid->m_gridNodes[1][0].y, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid->m_gridNodes[1][1].y, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid->m_gridNodes[1][2].y, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid->m_gridNodes[1][3].y, tolerance);

    ASSERT_NEAR(0.0, curvilinearGrid->m_gridNodes[2][0].y, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid->m_gridNodes[2][1].y, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid->m_gridNodes[2][2].y, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid->m_gridNodes[2][3].y, tolerance);

    ASSERT_NEAR(0.0, curvilinearGrid->m_gridNodes[3][0].y, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid->m_gridNodes[3][1].y, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid->m_gridNodes[3][2].y, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid->m_gridNodes[3][3].y, tolerance);
}

TEST(CurvilinearGridOrthogonalization, Compute_OnONonOrthogonalCurvilinearGrid_ShouldOrthogonalizeGrid)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();

    meshkernelapi::OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.OuterIterations = 2;
    orthogonalizationParameters.BoundaryIterations = 25;
    orthogonalizationParameters.InnerIterations = 25;
    orthogonalizationParameters.OrthogonalizationToSmoothingFactor = 0.975;
    meshkernel::CurvilinearGridOrthogonalization curvilinearGridOrthogonalization(curvilinearGrid, orthogonalizationParameters, {80154, 366530}, {80610, 367407});

    // Execute
    curvilinearGridOrthogonalization.Compute();

    // Assert
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(79983.796374595549, curvilinearGrid->m_gridNodes[0][0].x, tolerance);
    ASSERT_NEAR(80062.458148124875, curvilinearGrid->m_gridNodes[0][1].x, tolerance);
    ASSERT_NEAR(80139.231311212163, curvilinearGrid->m_gridNodes[0][2].x, tolerance);
    ASSERT_NEAR(80215.764666262330, curvilinearGrid->m_gridNodes[0][3].x, tolerance);
    ASSERT_NEAR(80293.275800678923, curvilinearGrid->m_gridNodes[0][4].x, tolerance);
    ASSERT_NEAR(80369.452240930041, curvilinearGrid->m_gridNodes[0][5].x, tolerance);
    ASSERT_NEAR(80446.051457860507, curvilinearGrid->m_gridNodes[0][6].x, tolerance);
    ASSERT_NEAR(80525.600108469254, curvilinearGrid->m_gridNodes[0][7].x, tolerance);
    ASSERT_NEAR(80611.260660513333, curvilinearGrid->m_gridNodes[0][8].x, tolerance);

    ASSERT_NEAR(80046.074998067590, curvilinearGrid->m_gridNodes[1][0].x, tolerance);
    ASSERT_NEAR(80126.368056777312, curvilinearGrid->m_gridNodes[1][1].x, tolerance);
    ASSERT_NEAR(80206.601210925364, curvilinearGrid->m_gridNodes[1][2].x, tolerance);
    ASSERT_NEAR(80286.287381013084, curvilinearGrid->m_gridNodes[1][3].x, tolerance);
    ASSERT_NEAR(80365.556974053165, curvilinearGrid->m_gridNodes[1][4].x, tolerance);
    ASSERT_NEAR(80443.587651289563, curvilinearGrid->m_gridNodes[1][5].x, tolerance);
    ASSERT_NEAR(80518.402282663534, curvilinearGrid->m_gridNodes[1][6].x, tolerance);
    ASSERT_NEAR(80591.265795857922, curvilinearGrid->m_gridNodes[1][7].x, tolerance);
    ASSERT_NEAR(80662.822394333387, curvilinearGrid->m_gridNodes[1][8].x, tolerance);

    ASSERT_NEAR(366936.89538054139, curvilinearGrid->m_gridNodes[0][0].y, tolerance);
    ASSERT_NEAR(366989.44915708830, curvilinearGrid->m_gridNodes[0][1].y, tolerance);
    ASSERT_NEAR(367042.59840712498, curvilinearGrid->m_gridNodes[0][2].y, tolerance);
    ASSERT_NEAR(367097.32852206373, curvilinearGrid->m_gridNodes[0][3].y, tolerance);
    ASSERT_NEAR(367153.90433931275, curvilinearGrid->m_gridNodes[0][4].y, tolerance);
    ASSERT_NEAR(367214.16882124962, curvilinearGrid->m_gridNodes[0][5].y, tolerance);
    ASSERT_NEAR(367275.49042299920, curvilinearGrid->m_gridNodes[0][6].y, tolerance);
    ASSERT_NEAR(367338.97499635734, curvilinearGrid->m_gridNodes[0][7].y, tolerance);
    ASSERT_NEAR(367407.38541954994, curvilinearGrid->m_gridNodes[0][8].y, tolerance);

    ASSERT_NEAR(366841.35146471433, curvilinearGrid->m_gridNodes[1][0].y, tolerance);
    ASSERT_NEAR(366887.54534671927, curvilinearGrid->m_gridNodes[1][1].y, tolerance);
    ASSERT_NEAR(366936.88064741943, curvilinearGrid->m_gridNodes[1][2].y, tolerance);
    ASSERT_NEAR(366990.14406133827, curvilinearGrid->m_gridNodes[1][3].y, tolerance);
    ASSERT_NEAR(367047.87322551024, curvilinearGrid->m_gridNodes[1][4].y, tolerance);
    ASSERT_NEAR(367110.98782724276, curvilinearGrid->m_gridNodes[1][5].y, tolerance);
    ASSERT_NEAR(367179.02693904703, curvilinearGrid->m_gridNodes[1][6].y, tolerance);
    ASSERT_NEAR(367254.16935723333, curvilinearGrid->m_gridNodes[1][7].y, tolerance);
    ASSERT_NEAR(367340.26390535629, curvilinearGrid->m_gridNodes[1][8].y, tolerance);
}