#include <gtest/gtest.h>

#include <MeshKernel/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGridSmoothing.hpp>
#include <MeshKernel/Entities.hpp>
#include <TestUtils/MakeCurvilinearGrids.hpp>

TEST(CurvilinearGridSmoothing, Compute_OnSmoothCurvilinearGrid_ShouldNotSmoothGrid)
{
    // Set-up
    std::vector<std::vector<meshkernel::Point>> grid{
        {{0, 0}, {0, 10}, {0, 20}, {0, 30}},
        {{10, 0}, {10, 10}, {10, 20}, {10, 30}},
        {{20, 0}, {20, 10}, {20, 20}, {20, 30}},
        {{30, 0}, {30, 10}, {30, 20}, {30, 30}}};

    const auto curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(std::move(grid), meshkernel::Projection::cartesian);

    meshkernel::CurvilinearGridSmoothing curvilinearGridSmoothing(curvilinearGrid, 10);

    // Execute
    curvilinearGridSmoothing.SetBlock({0, 0}, {30, 30});
    const auto smoothedGrid = curvilinearGridSmoothing.Compute();

    // Assert nodes are on the same location because the grid is already smooth
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(0.0, smoothedGrid.m_gridNodes[0][0].x, tolerance);
    ASSERT_NEAR(0.0, smoothedGrid.m_gridNodes[0][1].x, tolerance);
    ASSERT_NEAR(0.0, smoothedGrid.m_gridNodes[0][2].x, tolerance);
    ASSERT_NEAR(0.0, smoothedGrid.m_gridNodes[0][3].x, tolerance);

    ASSERT_NEAR(10.0, smoothedGrid.m_gridNodes[1][0].x, tolerance);
    ASSERT_NEAR(10.0, smoothedGrid.m_gridNodes[1][1].x, tolerance);
    ASSERT_NEAR(10.0, smoothedGrid.m_gridNodes[1][2].x, tolerance);
    ASSERT_NEAR(10.0, smoothedGrid.m_gridNodes[1][3].x, tolerance);

    ASSERT_NEAR(20.0, smoothedGrid.m_gridNodes[2][0].x, tolerance);
    ASSERT_NEAR(20.0, smoothedGrid.m_gridNodes[2][1].x, tolerance);
    ASSERT_NEAR(20.0, smoothedGrid.m_gridNodes[2][2].x, tolerance);
    ASSERT_NEAR(20.0, smoothedGrid.m_gridNodes[2][3].x, tolerance);

    ASSERT_NEAR(30.0, smoothedGrid.m_gridNodes[3][0].x, tolerance);
    ASSERT_NEAR(30.0, smoothedGrid.m_gridNodes[3][1].x, tolerance);
    ASSERT_NEAR(30.0, smoothedGrid.m_gridNodes[3][2].x, tolerance);
    ASSERT_NEAR(30.0, smoothedGrid.m_gridNodes[3][3].x, tolerance);

    ASSERT_NEAR(0.0, smoothedGrid.m_gridNodes[0][0].y, tolerance);
    ASSERT_NEAR(10.0, smoothedGrid.m_gridNodes[0][1].y, tolerance);
    ASSERT_NEAR(20.0, smoothedGrid.m_gridNodes[0][2].y, tolerance);
    ASSERT_NEAR(30.0, smoothedGrid.m_gridNodes[0][3].y, tolerance);

    ASSERT_NEAR(0.0, smoothedGrid.m_gridNodes[1][0].y, tolerance);
    ASSERT_NEAR(10.0, smoothedGrid.m_gridNodes[1][1].y, tolerance);
    ASSERT_NEAR(20.0, smoothedGrid.m_gridNodes[1][2].y, tolerance);
    ASSERT_NEAR(30.0, smoothedGrid.m_gridNodes[1][3].y, tolerance);

    ASSERT_NEAR(0.0, smoothedGrid.m_gridNodes[2][0].y, tolerance);
    ASSERT_NEAR(10.0, smoothedGrid.m_gridNodes[2][1].y, tolerance);
    ASSERT_NEAR(20.0, smoothedGrid.m_gridNodes[2][2].y, tolerance);
    ASSERT_NEAR(30.0, smoothedGrid.m_gridNodes[2][3].y, tolerance);

    ASSERT_NEAR(0.0, smoothedGrid.m_gridNodes[3][0].y, tolerance);
    ASSERT_NEAR(10.0, smoothedGrid.m_gridNodes[3][1].y, tolerance);
    ASSERT_NEAR(20.0, smoothedGrid.m_gridNodes[3][2].y, tolerance);
    ASSERT_NEAR(30.0, smoothedGrid.m_gridNodes[3][3].y, tolerance);
}
TEST(CurvilinearGridSmoothing, Compute_OnONonSmoothCurvilinearGrid_ShouldSmoothGrid)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();

    meshkernel::CurvilinearGridSmoothing curvilinearGridSmoothing(curvilinearGrid, 10);

    // Execute
    curvilinearGridSmoothing.SetBlock({80154, 366530}, {80610, 367407});
    const auto smoothedGrid = curvilinearGridSmoothing.Compute();

    // Assert
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(79983.796374595549, smoothedGrid.m_gridNodes[0][0].x, tolerance);
    ASSERT_NEAR(80060.011105254613, smoothedGrid.m_gridNodes[0][1].x, tolerance);
    ASSERT_NEAR(80137.131323077789, smoothedGrid.m_gridNodes[0][2].x, tolerance);
    ASSERT_NEAR(80214.289133171449, smoothedGrid.m_gridNodes[0][3].x, tolerance);
    ASSERT_NEAR(80290.623313344797, smoothedGrid.m_gridNodes[0][4].x, tolerance);
    ASSERT_NEAR(80364.411459799783, smoothedGrid.m_gridNodes[0][5].x, tolerance);
    ASSERT_NEAR(80440.632944042736, smoothedGrid.m_gridNodes[0][6].x, tolerance);
    ASSERT_NEAR(80521.863482944755, smoothedGrid.m_gridNodes[0][7].x, tolerance);
    ASSERT_NEAR(80611.260660513333, smoothedGrid.m_gridNodes[0][8].x, tolerance);

    ASSERT_NEAR(80049.839878894229, smoothedGrid.m_gridNodes[1][0].x, tolerance);
    ASSERT_NEAR(80126.599628477867, smoothedGrid.m_gridNodes[1][1].x, tolerance);
    ASSERT_NEAR(80209.094215192861, smoothedGrid.m_gridNodes[1][2].x, tolerance);
    ASSERT_NEAR(80293.744880346814, smoothedGrid.m_gridNodes[1][3].x, tolerance);
    ASSERT_NEAR(80379.324140649813, smoothedGrid.m_gridNodes[1][4].x, tolerance);
    ASSERT_NEAR(80465.542584239694, smoothedGrid.m_gridNodes[1][5].x, tolerance);
    ASSERT_NEAR(80554.276324305538, smoothedGrid.m_gridNodes[1][6].x, tolerance);
    ASSERT_NEAR(80648.676638649602, smoothedGrid.m_gridNodes[1][7].x, tolerance);
    ASSERT_NEAR(80756.869964790196, smoothedGrid.m_gridNodes[1][8].x, tolerance);

    ASSERT_NEAR(366936.89538054139, smoothedGrid.m_gridNodes[0][0].y, tolerance);
    ASSERT_NEAR(366987.83821424743, smoothedGrid.m_gridNodes[0][1].y, tolerance);
    ASSERT_NEAR(367041.23176860309, smoothedGrid.m_gridNodes[0][2].y, tolerance);
    ASSERT_NEAR(367096.32094200718, smoothedGrid.m_gridNodes[0][3].y, tolerance);
    ASSERT_NEAR(367152.27534010477, smoothedGrid.m_gridNodes[0][4].y, tolerance);
    ASSERT_NEAR(367210.34624999913, smoothedGrid.m_gridNodes[0][5].y, tolerance);
    ASSERT_NEAR(367271.13533544831, smoothedGrid.m_gridNodes[0][6].y, tolerance);
    ASSERT_NEAR(367335.99776061886, smoothedGrid.m_gridNodes[0][7].y, tolerance);
    ASSERT_NEAR(367407.38541954994, smoothedGrid.m_gridNodes[0][8].y, tolerance);

    ASSERT_NEAR(366833.73477307899, smoothedGrid.m_gridNodes[1][0].y, tolerance);
    ASSERT_NEAR(366875.95830515033, smoothedGrid.m_gridNodes[1][1].y, tolerance);
    ASSERT_NEAR(366920.06335723272, smoothedGrid.m_gridNodes[1][2].y, tolerance);
    ASSERT_NEAR(366964.51932115341, smoothedGrid.m_gridNodes[1][3].y, tolerance);
    ASSERT_NEAR(367008.36165438319, smoothedGrid.m_gridNodes[1][4].y, tolerance);
    ASSERT_NEAR(367051.96991358045, smoothedGrid.m_gridNodes[1][5].y, tolerance);
    ASSERT_NEAR(367097.07280781888, smoothedGrid.m_gridNodes[1][6].y, tolerance);
    ASSERT_NEAR(367147.66697242024, smoothedGrid.m_gridNodes[1][7].y, tolerance);
    ASSERT_NEAR(367208.53384889866, smoothedGrid.m_gridNodes[1][8].y, tolerance);
}

TEST(CurvilinearGridSmoothing, Compute_OnONonSmoothCurvilinearGridWithMissingElements_ShouldSmoothGrid)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGridWithMissingFaces();

    meshkernel::CurvilinearGridSmoothing curvilinearGridSmoothing(curvilinearGrid, 10);

    // Execute
    curvilinearGridSmoothing.SetBlock({80154, 366530}, {80610, 367407});
    const auto smoothedGrid = curvilinearGridSmoothing.Compute();

    // Assert
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(79983.796374595549, smoothedGrid.m_gridNodes[0][0].x, tolerance);
    ASSERT_NEAR(80060.788799524904, smoothedGrid.m_gridNodes[0][1].x, tolerance);
    ASSERT_NEAR(80138.938885926123, smoothedGrid.m_gridNodes[0][2].x, tolerance);
    ASSERT_NEAR(80216.491022095070, smoothedGrid.m_gridNodes[0][3].x, tolerance);
    ASSERT_NEAR(80293.293375308422, smoothedGrid.m_gridNodes[0][4].x, tolerance);
    ASSERT_NEAR(80367.153256326288, smoothedGrid.m_gridNodes[0][5].x, tolerance);
    ASSERT_NEAR(80441.975289048860, smoothedGrid.m_gridNodes[0][6].x, tolerance);
    ASSERT_NEAR(80522.272230797375, smoothedGrid.m_gridNodes[0][7].x, tolerance);
    ASSERT_NEAR(80611.260660513333, smoothedGrid.m_gridNodes[0][8].x, tolerance);

    ASSERT_NEAR(80050.309286359756, smoothedGrid.m_gridNodes[1][0].x, tolerance);
    ASSERT_NEAR(80129.710469510028, smoothedGrid.m_gridNodes[1][1].x, tolerance);
    ASSERT_NEAR(80217.157476954410, smoothedGrid.m_gridNodes[1][2].x, tolerance);
    ASSERT_NEAR(80301.619921323407, smoothedGrid.m_gridNodes[1][3].x, tolerance);
    ASSERT_NEAR(80387.883340666624, smoothedGrid.m_gridNodes[1][4].x, tolerance);
    ASSERT_NEAR(80476.978522055375, smoothedGrid.m_gridNodes[1][5].x, tolerance);
    ASSERT_NEAR(80558.979514368868, smoothedGrid.m_gridNodes[1][6].x, tolerance);
    ASSERT_NEAR(80649.960798514221, smoothedGrid.m_gridNodes[1][7].x, tolerance);
    ASSERT_NEAR(80756.940464565720, smoothedGrid.m_gridNodes[1][8].x, tolerance);

    ASSERT_NEAR(80098.932488001548, smoothedGrid.m_gridNodes[2][0].x, tolerance);
    ASSERT_NEAR(80186.124532376241, smoothedGrid.m_gridNodes[2][1].x, tolerance);
    ASSERT_NEAR(80292.449108019704, smoothedGrid.m_gridNodes[2][2].x, tolerance);
    ASSERT_NEAR(80375.822861451976, smoothedGrid.m_gridNodes[2][3].x, tolerance);
    ASSERT_NEAR(80468.946809326764, smoothedGrid.m_gridNodes[2][4].x, tolerance);
    ASSERT_NEAR(80583.120607369667, smoothedGrid.m_gridNodes[2][5].x, tolerance);
    ASSERT_NEAR(80652.102309464681, smoothedGrid.m_gridNodes[2][6].x, tolerance);
    ASSERT_NEAR(80749.982910696970, smoothedGrid.m_gridNodes[2][7].x, tolerance);
    ASSERT_NEAR(80871.931419427943, smoothedGrid.m_gridNodes[2][8].x, tolerance);

    ASSERT_NEAR(366936.89538054139, smoothedGrid.m_gridNodes[0][0].y, tolerance);
    ASSERT_NEAR(366988.35803436005, smoothedGrid.m_gridNodes[0][1].y, tolerance);
    ASSERT_NEAR(367042.48402813828, smoothedGrid.m_gridNodes[0][2].y, tolerance);
    ASSERT_NEAR(367097.89403811248, smoothedGrid.m_gridNodes[0][3].y, tolerance);
    ASSERT_NEAR(367154.23426078091, smoothedGrid.m_gridNodes[0][4].y, tolerance);
    ASSERT_NEAR(367212.51011039014, smoothedGrid.m_gridNodes[0][5].y, tolerance);
    ASSERT_NEAR(367272.20589329768, smoothedGrid.m_gridNodes[0][6].y, tolerance);
    ASSERT_NEAR(367336.32413904829, smoothedGrid.m_gridNodes[0][7].y, tolerance);
    ASSERT_NEAR(367407.38541954994, smoothedGrid.m_gridNodes[0][8].y, tolerance);

    ASSERT_NEAR(366833.00155397120, smoothedGrid.m_gridNodes[1][0].y, tolerance);
    ASSERT_NEAR(366875.36856874428, smoothedGrid.m_gridNodes[1][1].y, tolerance);
    ASSERT_NEAR(366919.22743417800, smoothedGrid.m_gridNodes[1][2].y, tolerance);
    ASSERT_NEAR(366964.89754991967, smoothedGrid.m_gridNodes[1][3].y, tolerance);
    ASSERT_NEAR(367010.93790833943, smoothedGrid.m_gridNodes[1][4].y, tolerance);
    ASSERT_NEAR(367054.20807785576, smoothedGrid.m_gridNodes[1][5].y, tolerance);
    ASSERT_NEAR(367097.95891073684, smoothedGrid.m_gridNodes[1][6].y, tolerance);
    ASSERT_NEAR(367147.87982618762, smoothedGrid.m_gridNodes[1][7].y, tolerance);
    ASSERT_NEAR(367208.42653128889, smoothedGrid.m_gridNodes[1][8].y, tolerance);

    ASSERT_NEAR(366729.04871992749, smoothedGrid.m_gridNodes[2][0].y, tolerance);
    ASSERT_NEAR(366761.65380535304, smoothedGrid.m_gridNodes[2][1].y, tolerance);
    ASSERT_NEAR(366792.50812354451, smoothedGrid.m_gridNodes[2][2].y, tolerance);
    ASSERT_NEAR(366827.70632165857, smoothedGrid.m_gridNodes[2][3].y, tolerance);
    ASSERT_NEAR(366865.78532145533, smoothedGrid.m_gridNodes[2][4].y, tolerance);
    ASSERT_NEAR(366891.62363194284, smoothedGrid.m_gridNodes[2][5].y, tolerance);
    ASSERT_NEAR(366916.02815869858, smoothedGrid.m_gridNodes[2][6].y, tolerance);
    ASSERT_NEAR(366948.34436165588, smoothedGrid.m_gridNodes[2][7].y, tolerance);
    ASSERT_NEAR(366996.75152488949, smoothedGrid.m_gridNodes[2][8].y, tolerance);
}

TEST(CurvilinearGridSmoothing, ComputedDirectionalSmooth_OnMDrirection_ShouldSmoothGrid)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();

    meshkernel::CurvilinearGridSmoothing curvilinearGridSmoothing(curvilinearGrid, 10);
    curvilinearGridSmoothing.SetLine({80143, 367041}, {80333, 366553});
    curvilinearGridSmoothing.SetBlock({80199, 366749}, {80480, 366869});

    // Execute
    const auto smoothedGrid = curvilinearGridSmoothing.ComputeDirectional();

    // Assert
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(80053.996925399639, smoothedGrid.m_gridNodes[1][0].x, tolerance);
    ASSERT_NEAR(80144.732940478571, smoothedGrid.m_gridNodes[1][1].x, tolerance);
    ASSERT_NEAR(80222.990051062123, smoothedGrid.m_gridNodes[1][2].x, tolerance);
    ASSERT_NEAR(80314.427829118606, smoothedGrid.m_gridNodes[1][3].x, tolerance);
    ASSERT_NEAR(80414.057391982613, smoothedGrid.m_gridNodes[1][4].x, tolerance);
    ASSERT_NEAR(80505.096476712482, smoothedGrid.m_gridNodes[1][5].x, tolerance);
    ASSERT_NEAR(80595.183339827883, smoothedGrid.m_gridNodes[1][6].x, tolerance);
    ASSERT_NEAR(80684.333994102650, smoothedGrid.m_gridNodes[1][7].x, tolerance);
    ASSERT_NEAR(80772.567299473958, smoothedGrid.m_gridNodes[1][8].x, tolerance);

    ASSERT_NEAR(366827.17869351729, smoothedGrid.m_gridNodes[1][0].y, tolerance);
    ASSERT_NEAR(366872.58359778317, smoothedGrid.m_gridNodes[1][1].y, tolerance);
    ASSERT_NEAR(366936.38429752010, smoothedGrid.m_gridNodes[1][2].y, tolerance);
    ASSERT_NEAR(366981.06765786058, smoothedGrid.m_gridNodes[1][3].y, tolerance);
    ASSERT_NEAR(367015.14849423966, smoothedGrid.m_gridNodes[1][4].y, tolerance);
    ASSERT_NEAR(367056.48898898275, smoothedGrid.m_gridNodes[1][5].y, tolerance);
    ASSERT_NEAR(367099.12347147451, smoothedGrid.m_gridNodes[1][6].y, tolerance);
    ASSERT_NEAR(367143.03018172452, smoothedGrid.m_gridNodes[1][7].y, tolerance);
    ASSERT_NEAR(367188.18349069095, smoothedGrid.m_gridNodes[1][8].y, tolerance);
}

TEST(CurvilinearGridSmoothing, ComputedDirectionalSmooth_OnNDrirection_ShouldSmoothGrid)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();

    meshkernel::CurvilinearGridSmoothing curvilinearGridSmoothing(curvilinearGrid, 10);
    curvilinearGridSmoothing.SetLine({80199, 366749}, {80480, 366869});
    curvilinearGridSmoothing.SetBlock({80143, 367041}, {80333, 366553});

    // Execute
    const auto smoothedGrid = curvilinearGridSmoothing.ComputeDirectional();

    // Assert
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(80053.996925399639, smoothedGrid.m_gridNodes[1][0].x, tolerance);
    ASSERT_NEAR(80144.692538122108, smoothedGrid.m_gridNodes[1][1].x, tolerance);
    ASSERT_NEAR(80234.737077531070, smoothedGrid.m_gridNodes[1][2].x, tolerance);
    ASSERT_NEAR(80324.189331719666, smoothedGrid.m_gridNodes[1][3].x, tolerance);
    ASSERT_NEAR(80413.125110631052, smoothedGrid.m_gridNodes[1][4].x, tolerance);
    ASSERT_NEAR(80505.096476712482, smoothedGrid.m_gridNodes[1][5].x, tolerance);
    ASSERT_NEAR(80595.183339827883, smoothedGrid.m_gridNodes[1][6].x, tolerance);
    ASSERT_NEAR(80684.333994102650, smoothedGrid.m_gridNodes[1][7].x, tolerance);
    ASSERT_NEAR(80772.567299473958, smoothedGrid.m_gridNodes[1][8].x, tolerance);

    ASSERT_NEAR(366827.17869351729, smoothedGrid.m_gridNodes[1][0].y, tolerance);
    ASSERT_NEAR(366872.56390471710, smoothedGrid.m_gridNodes[1][1].y, tolerance);
    ASSERT_NEAR(366919.09182483586, smoothedGrid.m_gridNodes[1][2].y, tolerance);
    ASSERT_NEAR(366966.51416593988, smoothedGrid.m_gridNodes[1][3].y, tolerance);
    ASSERT_NEAR(367014.64392142039, smoothedGrid.m_gridNodes[1][4].y, tolerance);
    ASSERT_NEAR(367056.48898898275, smoothedGrid.m_gridNodes[1][5].y, tolerance);
    ASSERT_NEAR(367099.12347147451, smoothedGrid.m_gridNodes[1][6].y, tolerance);
    ASSERT_NEAR(367143.03018172452, smoothedGrid.m_gridNodes[1][7].y, tolerance);
    ASSERT_NEAR(367188.18349069095, smoothedGrid.m_gridNodes[1][8].y, tolerance);
}
