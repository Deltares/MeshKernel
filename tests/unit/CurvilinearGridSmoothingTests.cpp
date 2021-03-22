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

    meshkernel::CurvilinearGridSmoothing curvilinearGridSmoothing(curvilinearGrid, 10, {0, 0}, {30, 30});

    // Execute
    curvilinearGridSmoothing.Compute();

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

// test these two
TEST(CurvilinearGridSmoothing, Compute_OnONonSmoothCurvilinearGrid_ShouldSmoothGrid)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();

    meshkernel::CurvilinearGridSmoothing curvilinearGridSmoothing(curvilinearGrid, 10, {80154, 366530}, {80610, 367407});

    // Execute
    curvilinearGridSmoothing.Compute();

    // Assert
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(79983.796374595549, curvilinearGrid->m_gridNodes[0][0].x, tolerance);
    ASSERT_NEAR(80060.011105254613, curvilinearGrid->m_gridNodes[0][1].x, tolerance);
    ASSERT_NEAR(80137.131323077789, curvilinearGrid->m_gridNodes[0][2].x, tolerance);
    ASSERT_NEAR(80214.289133171449, curvilinearGrid->m_gridNodes[0][3].x, tolerance);
    ASSERT_NEAR(80290.623313344797, curvilinearGrid->m_gridNodes[0][4].x, tolerance);
    ASSERT_NEAR(80364.411459799783, curvilinearGrid->m_gridNodes[0][5].x, tolerance);
    ASSERT_NEAR(80440.632944042736, curvilinearGrid->m_gridNodes[0][6].x, tolerance);
    ASSERT_NEAR(80521.863482944755, curvilinearGrid->m_gridNodes[0][7].x, tolerance);
    ASSERT_NEAR(80611.260660513333, curvilinearGrid->m_gridNodes[0][8].x, tolerance);

    ASSERT_NEAR(80049.839878894229, curvilinearGrid->m_gridNodes[1][0].x, tolerance);
    ASSERT_NEAR(80126.599628477867, curvilinearGrid->m_gridNodes[1][1].x, tolerance);
    ASSERT_NEAR(80209.094215192861, curvilinearGrid->m_gridNodes[1][2].x, tolerance);
    ASSERT_NEAR(80293.744880346814, curvilinearGrid->m_gridNodes[1][3].x, tolerance);
    ASSERT_NEAR(80379.324140649813, curvilinearGrid->m_gridNodes[1][4].x, tolerance);
    ASSERT_NEAR(80465.542584239694, curvilinearGrid->m_gridNodes[1][5].x, tolerance);
    ASSERT_NEAR(80554.276324305538, curvilinearGrid->m_gridNodes[1][6].x, tolerance);
    ASSERT_NEAR(80648.676638649602, curvilinearGrid->m_gridNodes[1][7].x, tolerance);
    ASSERT_NEAR(80756.869964790196, curvilinearGrid->m_gridNodes[1][8].x, tolerance);

    ASSERT_NEAR(366936.89538054139, curvilinearGrid->m_gridNodes[0][0].y, tolerance);
    ASSERT_NEAR(366987.83821424743, curvilinearGrid->m_gridNodes[0][1].y, tolerance);
    ASSERT_NEAR(367041.23176860309, curvilinearGrid->m_gridNodes[0][2].y, tolerance);
    ASSERT_NEAR(367096.32094200718, curvilinearGrid->m_gridNodes[0][3].y, tolerance);
    ASSERT_NEAR(367152.27534010477, curvilinearGrid->m_gridNodes[0][4].y, tolerance);
    ASSERT_NEAR(367210.34624999913, curvilinearGrid->m_gridNodes[0][5].y, tolerance);
    ASSERT_NEAR(367271.13533544831, curvilinearGrid->m_gridNodes[0][6].y, tolerance);
    ASSERT_NEAR(367335.99776061886, curvilinearGrid->m_gridNodes[0][7].y, tolerance);
    ASSERT_NEAR(367407.38541954994, curvilinearGrid->m_gridNodes[0][8].y, tolerance);

    ASSERT_NEAR(366833.73477307899, curvilinearGrid->m_gridNodes[1][0].y, tolerance);
    ASSERT_NEAR(366875.95830515033, curvilinearGrid->m_gridNodes[1][1].y, tolerance);
    ASSERT_NEAR(366920.06335723272, curvilinearGrid->m_gridNodes[1][2].y, tolerance);
    ASSERT_NEAR(366964.51932115341, curvilinearGrid->m_gridNodes[1][3].y, tolerance);
    ASSERT_NEAR(367008.36165438319, curvilinearGrid->m_gridNodes[1][4].y, tolerance);
    ASSERT_NEAR(367051.96991358045, curvilinearGrid->m_gridNodes[1][5].y, tolerance);
    ASSERT_NEAR(367097.07280781888, curvilinearGrid->m_gridNodes[1][6].y, tolerance);
    ASSERT_NEAR(367147.66697242024, curvilinearGrid->m_gridNodes[1][7].y, tolerance);
    ASSERT_NEAR(367208.53384889866, curvilinearGrid->m_gridNodes[1][8].y, tolerance);
}

TEST(CurvilinearGridSmoothing, Compute_OnONonSmoothCurvilinearGridWithMissingElements_ShouldSmoothGrid)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGridWithMissingFaces();

    meshkernel::CurvilinearGridSmoothing curvilinearGridSmoothing(curvilinearGrid, 10, {80154, 366530}, {80610, 367407});

    // Execute
    curvilinearGridSmoothing.Compute();

    // Assert
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(79983.796374595549, curvilinearGrid->m_gridNodes[0][0].x, tolerance);
    ASSERT_NEAR(80060.788799524904, curvilinearGrid->m_gridNodes[0][1].x, tolerance);
    ASSERT_NEAR(80138.938885926123, curvilinearGrid->m_gridNodes[0][2].x, tolerance);
    ASSERT_NEAR(80216.491022095070, curvilinearGrid->m_gridNodes[0][3].x, tolerance);
    ASSERT_NEAR(80293.293375308422, curvilinearGrid->m_gridNodes[0][4].x, tolerance);
    ASSERT_NEAR(80367.153256326288, curvilinearGrid->m_gridNodes[0][5].x, tolerance);
    ASSERT_NEAR(80441.975289048860, curvilinearGrid->m_gridNodes[0][6].x, tolerance);
    ASSERT_NEAR(80522.272230797375, curvilinearGrid->m_gridNodes[0][7].x, tolerance);
    ASSERT_NEAR(80611.260660513333, curvilinearGrid->m_gridNodes[0][8].x, tolerance);

    ASSERT_NEAR(80050.309286359756, curvilinearGrid->m_gridNodes[1][0].x, tolerance);
    ASSERT_NEAR(80129.710469510028, curvilinearGrid->m_gridNodes[1][1].x, tolerance);
    ASSERT_NEAR(80217.157476954410, curvilinearGrid->m_gridNodes[1][2].x, tolerance);
    ASSERT_NEAR(80301.619921323407, curvilinearGrid->m_gridNodes[1][3].x, tolerance);
    ASSERT_NEAR(80387.883340666624, curvilinearGrid->m_gridNodes[1][4].x, tolerance);
    ASSERT_NEAR(80476.978522055375, curvilinearGrid->m_gridNodes[1][5].x, tolerance);
    ASSERT_NEAR(80558.979514368868, curvilinearGrid->m_gridNodes[1][6].x, tolerance);
    ASSERT_NEAR(80649.960798514221, curvilinearGrid->m_gridNodes[1][7].x, tolerance);
    ASSERT_NEAR(80756.940464565720, curvilinearGrid->m_gridNodes[1][8].x, tolerance);

    ASSERT_NEAR(80098.932488001548, curvilinearGrid->m_gridNodes[2][0].x, tolerance);
    ASSERT_NEAR(80186.124532376241, curvilinearGrid->m_gridNodes[2][1].x, tolerance);
    ASSERT_NEAR(80292.449108019704, curvilinearGrid->m_gridNodes[2][2].x, tolerance);
    ASSERT_NEAR(80375.822861451976, curvilinearGrid->m_gridNodes[2][3].x, tolerance);
    ASSERT_NEAR(80468.946809326764, curvilinearGrid->m_gridNodes[2][4].x, tolerance);
    ASSERT_NEAR(80583.120607369667, curvilinearGrid->m_gridNodes[2][5].x, tolerance);
    ASSERT_NEAR(80652.102309464681, curvilinearGrid->m_gridNodes[2][6].x, tolerance);
    ASSERT_NEAR(80749.982910696970, curvilinearGrid->m_gridNodes[2][7].x, tolerance);
    ASSERT_NEAR(80871.931419427943, curvilinearGrid->m_gridNodes[2][8].x, tolerance);

    ASSERT_NEAR(366936.89538054139, curvilinearGrid->m_gridNodes[0][0].y, tolerance);
    ASSERT_NEAR(366988.35803436005, curvilinearGrid->m_gridNodes[0][1].y, tolerance);
    ASSERT_NEAR(367042.48402813828, curvilinearGrid->m_gridNodes[0][2].y, tolerance);
    ASSERT_NEAR(367097.89403811248, curvilinearGrid->m_gridNodes[0][3].y, tolerance);
    ASSERT_NEAR(367154.23426078091, curvilinearGrid->m_gridNodes[0][4].y, tolerance);
    ASSERT_NEAR(367212.51011039014, curvilinearGrid->m_gridNodes[0][5].y, tolerance);
    ASSERT_NEAR(367272.20589329768, curvilinearGrid->m_gridNodes[0][6].y, tolerance);
    ASSERT_NEAR(367336.32413904829, curvilinearGrid->m_gridNodes[0][7].y, tolerance);
    ASSERT_NEAR(367407.38541954994, curvilinearGrid->m_gridNodes[0][8].y, tolerance);

    ASSERT_NEAR(366833.00155397120, curvilinearGrid->m_gridNodes[1][0].y, tolerance);
    ASSERT_NEAR(366875.36856874428, curvilinearGrid->m_gridNodes[1][1].y, tolerance);
    ASSERT_NEAR(366919.22743417800, curvilinearGrid->m_gridNodes[1][2].y, tolerance);
    ASSERT_NEAR(366964.89754991967, curvilinearGrid->m_gridNodes[1][3].y, tolerance);
    ASSERT_NEAR(367010.93790833943, curvilinearGrid->m_gridNodes[1][4].y, tolerance);
    ASSERT_NEAR(367054.20807785576, curvilinearGrid->m_gridNodes[1][5].y, tolerance);
    ASSERT_NEAR(367097.95891073684, curvilinearGrid->m_gridNodes[1][6].y, tolerance);
    ASSERT_NEAR(367147.87982618762, curvilinearGrid->m_gridNodes[1][7].y, tolerance);
    ASSERT_NEAR(367208.42653128889, curvilinearGrid->m_gridNodes[1][8].y, tolerance);

    ASSERT_NEAR(366729.04871992749, curvilinearGrid->m_gridNodes[2][0].y, tolerance);
    ASSERT_NEAR(366761.65380535304, curvilinearGrid->m_gridNodes[2][1].y, tolerance);
    ASSERT_NEAR(366792.50812354451, curvilinearGrid->m_gridNodes[2][2].y, tolerance);
    ASSERT_NEAR(366827.70632165857, curvilinearGrid->m_gridNodes[2][3].y, tolerance);
    ASSERT_NEAR(366865.78532145533, curvilinearGrid->m_gridNodes[2][4].y, tolerance);
    ASSERT_NEAR(366891.62363194284, curvilinearGrid->m_gridNodes[2][5].y, tolerance);
    ASSERT_NEAR(366916.02815869858, curvilinearGrid->m_gridNodes[2][6].y, tolerance);
    ASSERT_NEAR(366948.34436165588, curvilinearGrid->m_gridNodes[2][7].y, tolerance);
    ASSERT_NEAR(366996.75152488949, curvilinearGrid->m_gridNodes[2][8].y, tolerance);
}
