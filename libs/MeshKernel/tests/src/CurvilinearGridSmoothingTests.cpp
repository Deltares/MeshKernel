#include <gtest/gtest.h>

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridSmoothing.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridSnapping.hpp>
#include <MeshKernel/Entities.hpp>
#include <TestUtils/MakeCurvilinearGrids.hpp>

using namespace meshkernel;

TEST(CurvilinearGridSmoothing, Compute_OnSmoothCurvilinearGrid_ShouldNotSmoothGrid)
{
    // Set-up
    lin_alg::Matrix<Point> grid(4, 4);
    grid << Point{0, 0}, Point{0, 10}, Point{0, 20}, Point{0, 30},
        Point{10, 0}, Point{10, 10}, Point{10, 20}, Point{10, 30},
        Point{20, 0}, Point{20, 10}, Point{20, 20}, Point{20, 30},
        Point{30, 0}, Point{30, 10}, Point{30, 20}, Point{30, 30};

    CurvilinearGrid curvilinearGrid(grid, Projection::cartesian);
    CurvilinearGridSmoothing curvilinearGridSmoothing(curvilinearGrid, 10);

    // Execute
    curvilinearGridSmoothing.SetBlock({0, 0}, {30, 30});
    curvilinearGridSmoothing.Compute();

    // Assert nodes are on the same location because the grid is already smooth
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(0.0, curvilinearGrid.GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(0.0, curvilinearGrid.GetNode(0, 1).x, tolerance);
    ASSERT_NEAR(0.0, curvilinearGrid.GetNode(0, 2).x, tolerance);
    ASSERT_NEAR(0.0, curvilinearGrid.GetNode(0, 3).x, tolerance);

    ASSERT_NEAR(10.0, curvilinearGrid.GetNode(1, 0).x, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid.GetNode(1, 1).x, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid.GetNode(1, 2).x, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid.GetNode(1, 3).x, tolerance);

    ASSERT_NEAR(20.0, curvilinearGrid.GetNode(2, 0).x, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid.GetNode(2, 1).x, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid.GetNode(2, 2).x, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid.GetNode(2, 3).x, tolerance);

    ASSERT_NEAR(30.0, curvilinearGrid.GetNode(3, 0).x, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid.GetNode(3, 1).x, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid.GetNode(3, 2).x, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid.GetNode(3, 3).x, tolerance);

    ASSERT_NEAR(0.0, curvilinearGrid.GetNode(0, 0).y, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid.GetNode(0, 1).y, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid.GetNode(0, 2).y, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid.GetNode(0, 3).y, tolerance);

    ASSERT_NEAR(0.0, curvilinearGrid.GetNode(1, 0).y, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid.GetNode(1, 1).y, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid.GetNode(1, 2).y, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid.GetNode(1, 3).y, tolerance);

    ASSERT_NEAR(0.0, curvilinearGrid.GetNode(2, 0).y, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid.GetNode(2, 1).y, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid.GetNode(2, 2).y, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid.GetNode(2, 3).y, tolerance);

    ASSERT_NEAR(0.0, curvilinearGrid.GetNode(3, 0).y, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid.GetNode(3, 1).y, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid.GetNode(3, 2).y, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid.GetNode(3, 3).y, tolerance);
}
TEST(CurvilinearGridSmoothing, Compute_OnONonSmoothCurvilinearGrid_ShouldSmoothGrid)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();

    CurvilinearGridSmoothing curvilinearGridSmoothing(*curvilinearGrid, 10);

    // Execute
    curvilinearGridSmoothing.SetBlock({80154, 366530}, {80610, 367407});
    curvilinearGridSmoothing.Compute();

    // Assert
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(79983.796374595549, curvilinearGrid->GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(80060.011105254613, curvilinearGrid->GetNode(0, 1).x, tolerance);
    ASSERT_NEAR(80137.131323077789, curvilinearGrid->GetNode(0, 2).x, tolerance);
    ASSERT_NEAR(80214.289133171449, curvilinearGrid->GetNode(0, 3).x, tolerance);
    ASSERT_NEAR(80290.623313344797, curvilinearGrid->GetNode(0, 4).x, tolerance);
    ASSERT_NEAR(80364.411459799783, curvilinearGrid->GetNode(0, 5).x, tolerance);
    ASSERT_NEAR(80440.632944042736, curvilinearGrid->GetNode(0, 6).x, tolerance);
    ASSERT_NEAR(80521.863482944755, curvilinearGrid->GetNode(0, 7).x, tolerance);
    ASSERT_NEAR(80611.260660513333, curvilinearGrid->GetNode(0, 8).x, tolerance);

    ASSERT_NEAR(80049.839878894229, curvilinearGrid->GetNode(1, 0).x, tolerance);
    ASSERT_NEAR(80126.599628477867, curvilinearGrid->GetNode(1, 1).x, tolerance);
    ASSERT_NEAR(80209.094215192861, curvilinearGrid->GetNode(1, 2).x, tolerance);
    ASSERT_NEAR(80293.744880346814, curvilinearGrid->GetNode(1, 3).x, tolerance);
    ASSERT_NEAR(80379.324140649813, curvilinearGrid->GetNode(1, 4).x, tolerance);
    ASSERT_NEAR(80465.542584239694, curvilinearGrid->GetNode(1, 5).x, tolerance);
    ASSERT_NEAR(80554.276324305538, curvilinearGrid->GetNode(1, 6).x, tolerance);
    ASSERT_NEAR(80648.676638649602, curvilinearGrid->GetNode(1, 7).x, tolerance);
    ASSERT_NEAR(80756.869964790196, curvilinearGrid->GetNode(1, 8).x, tolerance);

    ASSERT_NEAR(366936.89538054139, curvilinearGrid->GetNode(0, 0).y, tolerance);
    ASSERT_NEAR(366987.83821424743, curvilinearGrid->GetNode(0, 1).y, tolerance);
    ASSERT_NEAR(367041.23176860309, curvilinearGrid->GetNode(0, 2).y, tolerance);
    ASSERT_NEAR(367096.32094200718, curvilinearGrid->GetNode(0, 3).y, tolerance);
    ASSERT_NEAR(367152.27534010477, curvilinearGrid->GetNode(0, 4).y, tolerance);
    ASSERT_NEAR(367210.34624999913, curvilinearGrid->GetNode(0, 5).y, tolerance);
    ASSERT_NEAR(367271.13533544831, curvilinearGrid->GetNode(0, 6).y, tolerance);
    ASSERT_NEAR(367335.99776061886, curvilinearGrid->GetNode(0, 7).y, tolerance);
    ASSERT_NEAR(367407.38541954994, curvilinearGrid->GetNode(0, 8).y, tolerance);

    ASSERT_NEAR(366833.73477307899, curvilinearGrid->GetNode(1, 0).y, tolerance);
    ASSERT_NEAR(366875.95830515033, curvilinearGrid->GetNode(1, 1).y, tolerance);
    ASSERT_NEAR(366920.06335723272, curvilinearGrid->GetNode(1, 2).y, tolerance);
    ASSERT_NEAR(366964.51932115341, curvilinearGrid->GetNode(1, 3).y, tolerance);
    ASSERT_NEAR(367008.36165438319, curvilinearGrid->GetNode(1, 4).y, tolerance);
    ASSERT_NEAR(367051.96991358045, curvilinearGrid->GetNode(1, 5).y, tolerance);
    ASSERT_NEAR(367097.07280781888, curvilinearGrid->GetNode(1, 6).y, tolerance);
    ASSERT_NEAR(367147.66697242024, curvilinearGrid->GetNode(1, 7).y, tolerance);
    ASSERT_NEAR(367208.53384889866, curvilinearGrid->GetNode(1, 8).y, tolerance);
}

TEST(CurvilinearGridSmoothing, Compute_OnONonSmoothCurvilinearGridWithMissingElements_ShouldSmoothGrid)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGridWithMissingFaces();

    CurvilinearGridSmoothing curvilinearGridSmoothing(*curvilinearGrid, 10);

    // Execute
    curvilinearGridSmoothing.SetBlock({80154, 366530}, {80610, 367407});
    curvilinearGridSmoothing.Compute();

    // Assert
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(79983.796374595549, curvilinearGrid->GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(80060.788799524904, curvilinearGrid->GetNode(0, 1).x, tolerance);
    ASSERT_NEAR(80138.938885926123, curvilinearGrid->GetNode(0, 2).x, tolerance);
    ASSERT_NEAR(80216.491022095070, curvilinearGrid->GetNode(0, 3).x, tolerance);
    ASSERT_NEAR(80293.293375308422, curvilinearGrid->GetNode(0, 4).x, tolerance);
    ASSERT_NEAR(80367.153256326288, curvilinearGrid->GetNode(0, 5).x, tolerance);
    ASSERT_NEAR(80441.975289048860, curvilinearGrid->GetNode(0, 6).x, tolerance);
    ASSERT_NEAR(80522.272230797375, curvilinearGrid->GetNode(0, 7).x, tolerance);
    ASSERT_NEAR(80611.260660513333, curvilinearGrid->GetNode(0, 8).x, tolerance);

    ASSERT_NEAR(80050.309286359756, curvilinearGrid->GetNode(1, 0).x, tolerance);
    ASSERT_NEAR(80129.710469510028, curvilinearGrid->GetNode(1, 1).x, tolerance);
    ASSERT_NEAR(80217.157476954410, curvilinearGrid->GetNode(1, 2).x, tolerance);
    ASSERT_NEAR(80301.619921323407, curvilinearGrid->GetNode(1, 3).x, tolerance);
    ASSERT_NEAR(80387.883340666624, curvilinearGrid->GetNode(1, 4).x, tolerance);
    ASSERT_NEAR(80476.978522055375, curvilinearGrid->GetNode(1, 5).x, tolerance);
    ASSERT_NEAR(80558.979514368868, curvilinearGrid->GetNode(1, 6).x, tolerance);
    ASSERT_NEAR(80649.960798514221, curvilinearGrid->GetNode(1, 7).x, tolerance);
    ASSERT_NEAR(80756.940464565720, curvilinearGrid->GetNode(1, 8).x, tolerance);

    ASSERT_NEAR(80098.932488001548, curvilinearGrid->GetNode(2, 0).x, tolerance);
    ASSERT_NEAR(80186.124532376241, curvilinearGrid->GetNode(2, 1).x, tolerance);
    ASSERT_NEAR(80292.449108019704, curvilinearGrid->GetNode(2, 2).x, tolerance);
    ASSERT_NEAR(80375.822861451976, curvilinearGrid->GetNode(2, 3).x, tolerance);
    ASSERT_NEAR(80468.946809326764, curvilinearGrid->GetNode(2, 4).x, tolerance);
    ASSERT_NEAR(80583.120607369667, curvilinearGrid->GetNode(2, 5).x, tolerance);
    ASSERT_NEAR(80652.102309464681, curvilinearGrid->GetNode(2, 6).x, tolerance);
    ASSERT_NEAR(80749.982910696970, curvilinearGrid->GetNode(2, 7).x, tolerance);
    ASSERT_NEAR(80871.931419427943, curvilinearGrid->GetNode(2, 8).x, tolerance);

    ASSERT_NEAR(366936.89538054139, curvilinearGrid->GetNode(0, 0).y, tolerance);
    ASSERT_NEAR(366988.35803436005, curvilinearGrid->GetNode(0, 1).y, tolerance);
    ASSERT_NEAR(367042.48402813828, curvilinearGrid->GetNode(0, 2).y, tolerance);
    ASSERT_NEAR(367097.89403811248, curvilinearGrid->GetNode(0, 3).y, tolerance);
    ASSERT_NEAR(367154.23426078091, curvilinearGrid->GetNode(0, 4).y, tolerance);
    ASSERT_NEAR(367212.51011039014, curvilinearGrid->GetNode(0, 5).y, tolerance);
    ASSERT_NEAR(367272.20589329768, curvilinearGrid->GetNode(0, 6).y, tolerance);
    ASSERT_NEAR(367336.32413904829, curvilinearGrid->GetNode(0, 7).y, tolerance);
    ASSERT_NEAR(367407.38541954994, curvilinearGrid->GetNode(0, 8).y, tolerance);

    ASSERT_NEAR(366833.00155397120, curvilinearGrid->GetNode(1, 0).y, tolerance);
    ASSERT_NEAR(366875.36856874428, curvilinearGrid->GetNode(1, 1).y, tolerance);
    ASSERT_NEAR(366919.22743417800, curvilinearGrid->GetNode(1, 2).y, tolerance);
    ASSERT_NEAR(366964.89754991967, curvilinearGrid->GetNode(1, 3).y, tolerance);
    ASSERT_NEAR(367010.93790833943, curvilinearGrid->GetNode(1, 4).y, tolerance);
    ASSERT_NEAR(367054.20807785576, curvilinearGrid->GetNode(1, 5).y, tolerance);
    ASSERT_NEAR(367097.95891073684, curvilinearGrid->GetNode(1, 6).y, tolerance);
    ASSERT_NEAR(367147.87982618762, curvilinearGrid->GetNode(1, 7).y, tolerance);
    ASSERT_NEAR(367208.42653128889, curvilinearGrid->GetNode(1, 8).y, tolerance);

    ASSERT_NEAR(366729.04871992749, curvilinearGrid->GetNode(2, 0).y, tolerance);
    ASSERT_NEAR(366761.65380535304, curvilinearGrid->GetNode(2, 1).y, tolerance);
    ASSERT_NEAR(366792.50812354451, curvilinearGrid->GetNode(2, 2).y, tolerance);
    ASSERT_NEAR(366827.70632165857, curvilinearGrid->GetNode(2, 3).y, tolerance);
    ASSERT_NEAR(366865.78532145533, curvilinearGrid->GetNode(2, 4).y, tolerance);
    ASSERT_NEAR(366891.62363194284, curvilinearGrid->GetNode(2, 5).y, tolerance);
    ASSERT_NEAR(366916.02815869858, curvilinearGrid->GetNode(2, 6).y, tolerance);
    ASSERT_NEAR(366948.34436165588, curvilinearGrid->GetNode(2, 7).y, tolerance);
    ASSERT_NEAR(366996.75152488949, curvilinearGrid->GetNode(2, 8).y, tolerance);
}

TEST(CurvilinearGridSmoothing, ComputedDirectionalSmooth_OnMDrirection_ShouldSmoothGrid)
{
    // Set-up
    auto curvilinearGrid = MakeSmallCurvilinearGrid();

    CurvilinearGridSmoothing curvilinearGridSmoothing(*curvilinearGrid, 10);
    curvilinearGridSmoothing.SetLine({80143, 367041}, {80333, 366553});
    curvilinearGridSmoothing.SetBlock({80199, 366749}, {80480, 366869});

    // Execute
    curvilinearGridSmoothing.ComputeDirectional();

    // Assert
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(80053.996925399639, curvilinearGrid->GetNode(1, 0).x, tolerance);
    ASSERT_NEAR(80144.732940478571, curvilinearGrid->GetNode(1, 1).x, tolerance);
    ASSERT_NEAR(80222.990051062123, curvilinearGrid->GetNode(1, 2).x, tolerance);
    ASSERT_NEAR(80314.427829118606, curvilinearGrid->GetNode(1, 3).x, tolerance);
    ASSERT_NEAR(80414.057391982613, curvilinearGrid->GetNode(1, 4).x, tolerance);
    ASSERT_NEAR(80505.096476712482, curvilinearGrid->GetNode(1, 5).x, tolerance);
    ASSERT_NEAR(80595.183339827883, curvilinearGrid->GetNode(1, 6).x, tolerance);
    ASSERT_NEAR(80684.333994102650, curvilinearGrid->GetNode(1, 7).x, tolerance);
    ASSERT_NEAR(80772.567299473958, curvilinearGrid->GetNode(1, 8).x, tolerance);

    ASSERT_NEAR(366827.17869351729, curvilinearGrid->GetNode(1, 0).y, tolerance);
    ASSERT_NEAR(366872.58359778317, curvilinearGrid->GetNode(1, 1).y, tolerance);
    ASSERT_NEAR(366936.38429752010, curvilinearGrid->GetNode(1, 2).y, tolerance);
    ASSERT_NEAR(366981.06765786058, curvilinearGrid->GetNode(1, 3).y, tolerance);
    ASSERT_NEAR(367015.14849423966, curvilinearGrid->GetNode(1, 4).y, tolerance);
    ASSERT_NEAR(367056.48898898275, curvilinearGrid->GetNode(1, 5).y, tolerance);
    ASSERT_NEAR(367099.12347147451, curvilinearGrid->GetNode(1, 6).y, tolerance);
    ASSERT_NEAR(367143.03018172452, curvilinearGrid->GetNode(1, 7).y, tolerance);
    ASSERT_NEAR(367188.18349069095, curvilinearGrid->GetNode(1, 8).y, tolerance);
}

TEST(CurvilinearGridSmoothing, ComputedDirectionalSmooth_OnNDrirection_ShouldSmoothGrid)
{
    // Set-up
    auto curvilinearGrid = MakeSmallCurvilinearGrid();

    CurvilinearGridSmoothing curvilinearGridSmoothing(*curvilinearGrid, 10);
    curvilinearGridSmoothing.SetLine({80199, 366749}, {80480, 366869});
    curvilinearGridSmoothing.SetBlock({80143, 367041}, {80333, 366553});

    // Execute
    curvilinearGridSmoothing.ComputeDirectional();

    // Assert
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(80053.996925399639, curvilinearGrid->GetNode(1, 0).x, tolerance);
    ASSERT_NEAR(80144.692538122108, curvilinearGrid->GetNode(1, 1).x, tolerance);
    ASSERT_NEAR(80234.737077531070, curvilinearGrid->GetNode(1, 2).x, tolerance);
    ASSERT_NEAR(80324.189331719666, curvilinearGrid->GetNode(1, 3).x, tolerance);
    ASSERT_NEAR(80413.125110631052, curvilinearGrid->GetNode(1, 4).x, tolerance);
    ASSERT_NEAR(80505.096476712482, curvilinearGrid->GetNode(1, 5).x, tolerance);
    ASSERT_NEAR(80595.183339827883, curvilinearGrid->GetNode(1, 6).x, tolerance);
    ASSERT_NEAR(80684.333994102650, curvilinearGrid->GetNode(1, 7).x, tolerance);
    ASSERT_NEAR(80772.567299473958, curvilinearGrid->GetNode(1, 8).x, tolerance);

    ASSERT_NEAR(366827.17869351729, curvilinearGrid->GetNode(1, 0).y, tolerance);
    ASSERT_NEAR(366872.56390471710, curvilinearGrid->GetNode(1, 1).y, tolerance);
    ASSERT_NEAR(366919.09182483586, curvilinearGrid->GetNode(1, 2).y, tolerance);
    ASSERT_NEAR(366966.51416593988, curvilinearGrid->GetNode(1, 3).y, tolerance);
    ASSERT_NEAR(367014.64392142039, curvilinearGrid->GetNode(1, 4).y, tolerance);
    ASSERT_NEAR(367056.48898898275, curvilinearGrid->GetNode(1, 5).y, tolerance);
    ASSERT_NEAR(367099.12347147451, curvilinearGrid->GetNode(1, 6).y, tolerance);
    ASSERT_NEAR(367143.03018172452, curvilinearGrid->GetNode(1, 7).y, tolerance);
    ASSERT_NEAR(367188.18349069095, curvilinearGrid->GetNode(1, 8).y, tolerance);
}
