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
    ASSERT_NEAR(80067.930920933941, curvilinearGrid->m_gridNodes[0][1].x, tolerance);
    ASSERT_NEAR(80150.164245250286, curvilinearGrid->m_gridNodes[0][2].x, tolerance);
    ASSERT_NEAR(80231.651392830841, curvilinearGrid->m_gridNodes[0][3].x, tolerance);
    ASSERT_NEAR(80312.801762204603, curvilinearGrid->m_gridNodes[0][4].x, tolerance);
    ASSERT_NEAR(80390.860703898958, curvilinearGrid->m_gridNodes[0][5].x, tolerance);
    ASSERT_NEAR(80466.956631008565, curvilinearGrid->m_gridNodes[0][6].x, tolerance);
    ASSERT_NEAR(80540.510527743099, curvilinearGrid->m_gridNodes[0][7].x, tolerance);
    ASSERT_NEAR(80611.260660513333, curvilinearGrid->m_gridNodes[0][8].x, tolerance);

    ASSERT_NEAR(80056.139417226193, curvilinearGrid->m_gridNodes[1][0].x, tolerance);
    ASSERT_NEAR(80143.392341684361, curvilinearGrid->m_gridNodes[1][1].x, tolerance);
    ASSERT_NEAR(80232.701973556381, curvilinearGrid->m_gridNodes[1][2].x, tolerance);
    ASSERT_NEAR(80323.273756955823, curvilinearGrid->m_gridNodes[1][3].x, tolerance);
    ASSERT_NEAR(80415.126141636414, curvilinearGrid->m_gridNodes[1][4].x, tolerance);
    ASSERT_NEAR(80506.879117047443, curvilinearGrid->m_gridNodes[1][5].x, tolerance);
    ASSERT_NEAR(80595.173874885382, curvilinearGrid->m_gridNodes[1][6].x, tolerance);
    ASSERT_NEAR(80680.093103208215, curvilinearGrid->m_gridNodes[1][7].x, tolerance);
    ASSERT_NEAR(80761.372120477128, curvilinearGrid->m_gridNodes[1][8].x, tolerance);

    ASSERT_NEAR(366936.89538054139, curvilinearGrid->m_gridNodes[0][0].y, tolerance);
    ASSERT_NEAR(366993.16520344437, curvilinearGrid->m_gridNodes[0][1].y, tolerance);
    ASSERT_NEAR(367050.33461904334, curvilinearGrid->m_gridNodes[0][2].y, tolerance);
    ASSERT_NEAR(367108.77222890750, curvilinearGrid->m_gridNodes[0][3].y, tolerance);
    ASSERT_NEAR(367168.86295006215, curvilinearGrid->m_gridNodes[0][4].y, tolerance);
    ASSERT_NEAR(367231.44573248079, curvilinearGrid->m_gridNodes[0][5].y, tolerance);
    ASSERT_NEAR(367292.13654755038, curvilinearGrid->m_gridNodes[0][6].y, tolerance);
    ASSERT_NEAR(367350.89140520856, curvilinearGrid->m_gridNodes[0][7].y, tolerance);
    ASSERT_NEAR(367407.38541954994, curvilinearGrid->m_gridNodes[0][8].y, tolerance);

    ASSERT_NEAR(366823.20595668437, curvilinearGrid->m_gridNodes[1][0].y, tolerance);
    ASSERT_NEAR(366869.41001402331, curvilinearGrid->m_gridNodes[1][1].y, tolerance);
    ASSERT_NEAR(366915.26601702173, curvilinearGrid->m_gridNodes[1][2].y, tolerance);
    ASSERT_NEAR(366961.36595671845, curvilinearGrid->m_gridNodes[1][3].y, tolerance);
    ASSERT_NEAR(367008.20690690872, curvilinearGrid->m_gridNodes[1][4].y, tolerance);
    ASSERT_NEAR(367056.34885035310, curvilinearGrid->m_gridNodes[1][5].y, tolerance);
    ASSERT_NEAR(367105.00611276925, curvilinearGrid->m_gridNodes[1][6].y, tolerance);
    ASSERT_NEAR(367154.47604494903, curvilinearGrid->m_gridNodes[1][7].y, tolerance);
    ASSERT_NEAR(367204.88319783867, curvilinearGrid->m_gridNodes[1][8].y, tolerance);
}

TEST(CurvilinearGridOrthogonalization, Compute_OnONonOrthogonalCurvilinearGridWithMissingElements_ShouldOrthogonalizeGrid)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGridWithMissingFaces();

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
    ASSERT_NEAR(80069.224559424794, curvilinearGrid->m_gridNodes[0][1].x, tolerance);
    ASSERT_NEAR(80152.818443899014, curvilinearGrid->m_gridNodes[0][2].x, tolerance);
    ASSERT_NEAR(80235.610419184784, curvilinearGrid->m_gridNodes[0][3].x, tolerance);
    ASSERT_NEAR(80317.467404251467, curvilinearGrid->m_gridNodes[0][4].x, tolerance);
    ASSERT_NEAR(80395.149908970227, curvilinearGrid->m_gridNodes[0][5].x, tolerance);
    ASSERT_NEAR(80470.021886679329, curvilinearGrid->m_gridNodes[0][6].x, tolerance);
    ASSERT_NEAR(80542.049673235932, curvilinearGrid->m_gridNodes[0][7].x, tolerance);
    ASSERT_NEAR(80611.260660513333, curvilinearGrid->m_gridNodes[0][8].x, tolerance);

    ASSERT_NEAR(80055.199856352847, curvilinearGrid->m_gridNodes[1][0].x, tolerance);
    ASSERT_NEAR(80143.721270425682, curvilinearGrid->m_gridNodes[1][1].x, tolerance);
    ASSERT_NEAR(80234.545805664267, curvilinearGrid->m_gridNodes[1][2].x, tolerance);
    ASSERT_NEAR(80326.294071799057, curvilinearGrid->m_gridNodes[1][3].x, tolerance);
    ASSERT_NEAR(80418.809908811469, curvilinearGrid->m_gridNodes[1][4].x, tolerance);
    ASSERT_NEAR(80510.300359235334, curvilinearGrid->m_gridNodes[1][5].x, tolerance);
    ASSERT_NEAR(80597.036737342059, curvilinearGrid->m_gridNodes[1][6].x, tolerance);
    ASSERT_NEAR(80680.228644184652, curvilinearGrid->m_gridNodes[1][7].x, tolerance);
    ASSERT_NEAR(80759.822363775238, curvilinearGrid->m_gridNodes[1][8].x, tolerance);

    ASSERT_NEAR(80104.514488846587, curvilinearGrid->m_gridNodes[2][0].x, tolerance);
    ASSERT_NEAR(80197.383457403659, curvilinearGrid->m_gridNodes[2][1].x, tolerance);
    ASSERT_NEAR(80292.449108019704, curvilinearGrid->m_gridNodes[2][2].x, tolerance);
    ASSERT_NEAR(80387.175475807715, curvilinearGrid->m_gridNodes[2][3].x, tolerance);
    ASSERT_NEAR(80480.608251576487, curvilinearGrid->m_gridNodes[2][4].x, tolerance);
    ASSERT_NEAR(80583.120607369667, curvilinearGrid->m_gridNodes[2][5].x, tolerance);
    ASSERT_NEAR(80682.512842750177, curvilinearGrid->m_gridNodes[2][6].x, tolerance);
    ASSERT_NEAR(80780.894793453088, curvilinearGrid->m_gridNodes[2][7].x, tolerance);
    ASSERT_NEAR(80877.580909293247, curvilinearGrid->m_gridNodes[2][8].x, tolerance);

    ASSERT_NEAR(366936.89538054139, curvilinearGrid->m_gridNodes[0][0].y, tolerance);
    ASSERT_NEAR(366994.04530248570, curvilinearGrid->m_gridNodes[0][1].y, tolerance);
    ASSERT_NEAR(367052.21844658040, curvilinearGrid->m_gridNodes[0][2].y, tolerance);
    ASSERT_NEAR(367111.62516794214, curvilinearGrid->m_gridNodes[0][3].y, tolerance);
    ASSERT_NEAR(367172.50539615576, curvilinearGrid->m_gridNodes[0][4].y, tolerance);
    ASSERT_NEAR(367234.88964184484, curvilinearGrid->m_gridNodes[0][5].y, tolerance);
    ASSERT_NEAR(367294.58040378935, curvilinearGrid->m_gridNodes[0][6].y, tolerance);
    ASSERT_NEAR(367352.12121038162, curvilinearGrid->m_gridNodes[0][7].y, tolerance);
    ASSERT_NEAR(367407.38541954994, curvilinearGrid->m_gridNodes[0][8].y, tolerance);

    ASSERT_NEAR(366824.95581170503, curvilinearGrid->m_gridNodes[1][0].y, tolerance);
    ASSERT_NEAR(366871.97101464303, curvilinearGrid->m_gridNodes[1][1].y, tolerance);
    ASSERT_NEAR(366918.90347823157, curvilinearGrid->m_gridNodes[1][2].y, tolerance);
    ASSERT_NEAR(366966.74293227377, curvilinearGrid->m_gridNodes[1][3].y, tolerance);
    ASSERT_NEAR(367014.95766231039, curvilinearGrid->m_gridNodes[1][4].y, tolerance);
    ASSERT_NEAR(367062.00799990579, curvilinearGrid->m_gridNodes[1][5].y, tolerance);
    ASSERT_NEAR(367109.50803491072, curvilinearGrid->m_gridNodes[1][6].y, tolerance);
    ASSERT_NEAR(367157.80634231429, curvilinearGrid->m_gridNodes[1][7].y, tolerance);
    ASSERT_NEAR(367207.16259613103, curvilinearGrid->m_gridNodes[1][8].y, tolerance);

    ASSERT_NEAR(366718.00768776325, curvilinearGrid->m_gridNodes[2][0].y, tolerance);
    ASSERT_NEAR(366755.34796923219, curvilinearGrid->m_gridNodes[2][1].y, tolerance);
    ASSERT_NEAR(366792.50812354451, curvilinearGrid->m_gridNodes[2][2].y, tolerance);
    ASSERT_NEAR(366832.52615748829, curvilinearGrid->m_gridNodes[2][3].y, tolerance);
    ASSERT_NEAR(366870.46128228144, curvilinearGrid->m_gridNodes[2][4].y, tolerance);
    ASSERT_NEAR(366891.62363194284, curvilinearGrid->m_gridNodes[2][5].y, tolerance);
    ASSERT_NEAR(366923.13989736629, curvilinearGrid->m_gridNodes[2][6].y, tolerance);
    ASSERT_NEAR(366957.61315351113, curvilinearGrid->m_gridNodes[2][7].y, tolerance);
    ASSERT_NEAR(366996.07892524434, curvilinearGrid->m_gridNodes[2][8].y, tolerance);
}
