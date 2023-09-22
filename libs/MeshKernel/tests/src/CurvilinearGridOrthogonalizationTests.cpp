#include <gtest/gtest.h>

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridOrthogonalization.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Parameters.hpp>
#include <MeshKernel/Utilities/LinearAlgebra.hpp>
#include <TestUtils/MakeCurvilinearGrids.hpp>

using namespace meshkernel;

TEST(CurvilinearGridOrthogonalization, Compute_OnStronglyNonOrthogonalCurvilinearGrid_ShouldOrthogonalizeGrid)
{
    // Set-up
    lin_alg::Matrix<Point> grid(4, 4);
    grid << Point{0, 0}, Point{0, 10}, Point{0, 20}, Point{0, 30},
        Point{10, 0}, Point{10, 10}, Point{10, 20}, Point{10, 30},
        Point{20, 0}, Point{20, 10}, Point{20, 20}, Point{20, 30},
        Point{30, 0}, Point{30, 10}, Point{30, 20}, Point{30, 30};

    const auto curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(grid, meshkernel::Projection::Cartesian);

    // Move a node, to make the grid strongly non orthogonal
    curvilinearGrid->MoveNode({10.0, 20.0}, {18.0, 12.0});

    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 1;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    meshkernel::CurvilinearGridOrthogonalization curvilinearGridOrthogonalization(curvilinearGrid, orthogonalizationParameters);
    curvilinearGridOrthogonalization.SetBlock({0, 0}, {30, 30});

    // Execute
    auto const orthogonalizedCurvilinearGrid = curvilinearGridOrthogonalization.Compute();

    // Assert the moved nodes has moved towards its original location, making the grid more orthogonal
    constexpr double tolerance = 1e-6;
    ASSERT_NEAR(11.841396536135521, orthogonalizedCurvilinearGrid.m_gridNodes(1, 2).x, tolerance);
    ASSERT_NEAR(18.158586078094562, orthogonalizedCurvilinearGrid.m_gridNodes(1, 2).y, tolerance);
}

TEST(CurvilinearGridOrthogonalization, Compute_OnOrthogonalCurvilinearGrid_ShouldNotModifyGrid)
{
    // Set-up
    lin_alg::Matrix<Point> grid(4, 4);
    grid << Point{0, 0}, Point{0, 10}, Point{0, 20}, Point{0, 30},
        Point{10, 0}, Point{10, 10}, Point{10, 20}, Point{10, 30},
        Point{20, 0}, Point{20, 10}, Point{20, 20}, Point{20, 30},
        Point{30, 0}, Point{30, 10}, Point{30, 20}, Point{30, 30};

    const auto curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(grid, meshkernel::Projection::Cartesian);

    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 1;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    meshkernel::CurvilinearGridOrthogonalization curvilinearGridOrthogonalization(curvilinearGrid, orthogonalizationParameters);
    curvilinearGridOrthogonalization.SetBlock({0, 0}, {30, 30});

    // Execute
    auto const orthogonalizedCurvilinearGrid = curvilinearGridOrthogonalization.Compute();

    // Assert nodes are on the same location because the grid is already orthogonal
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(0.0, orthogonalizedCurvilinearGrid.m_gridNodes(0, 0).x, tolerance);
    ASSERT_NEAR(0.0, orthogonalizedCurvilinearGrid.m_gridNodes(0, 1).x, tolerance);
    ASSERT_NEAR(0.0, orthogonalizedCurvilinearGrid.m_gridNodes(0, 2).x, tolerance);
    ASSERT_NEAR(0.0, orthogonalizedCurvilinearGrid.m_gridNodes(0, 3).x, tolerance);

    ASSERT_NEAR(10.0, orthogonalizedCurvilinearGrid.m_gridNodes(1, 0).x, tolerance);
    ASSERT_NEAR(10.0, orthogonalizedCurvilinearGrid.m_gridNodes(1, 1).x, tolerance);
    ASSERT_NEAR(10.0, orthogonalizedCurvilinearGrid.m_gridNodes(1, 2).x, tolerance);
    ASSERT_NEAR(10.0, orthogonalizedCurvilinearGrid.m_gridNodes(1, 3).x, tolerance);

    ASSERT_NEAR(20.0, orthogonalizedCurvilinearGrid.m_gridNodes(2, 0).x, tolerance);
    ASSERT_NEAR(20.0, orthogonalizedCurvilinearGrid.m_gridNodes(2, 1).x, tolerance);
    ASSERT_NEAR(20.0, orthogonalizedCurvilinearGrid.m_gridNodes(2, 2).x, tolerance);
    ASSERT_NEAR(20.0, orthogonalizedCurvilinearGrid.m_gridNodes(2, 3).x, tolerance);

    ASSERT_NEAR(30.0, orthogonalizedCurvilinearGrid.m_gridNodes(3, 0).x, tolerance);
    ASSERT_NEAR(30.0, orthogonalizedCurvilinearGrid.m_gridNodes(3, 1).x, tolerance);
    ASSERT_NEAR(30.0, orthogonalizedCurvilinearGrid.m_gridNodes(3, 2).x, tolerance);
    ASSERT_NEAR(30.0, orthogonalizedCurvilinearGrid.m_gridNodes(3, 3).x, tolerance);

    ASSERT_NEAR(0.0, orthogonalizedCurvilinearGrid.m_gridNodes(0, 0).y, tolerance);
    ASSERT_NEAR(10.0, orthogonalizedCurvilinearGrid.m_gridNodes(0, 1).y, tolerance);
    ASSERT_NEAR(20.0, orthogonalizedCurvilinearGrid.m_gridNodes(0, 2).y, tolerance);
    ASSERT_NEAR(30.0, orthogonalizedCurvilinearGrid.m_gridNodes(0, 3).y, tolerance);

    ASSERT_NEAR(0.0, orthogonalizedCurvilinearGrid.m_gridNodes(1, 0).y, tolerance);
    ASSERT_NEAR(10.0, orthogonalizedCurvilinearGrid.m_gridNodes(1, 1).y, tolerance);
    ASSERT_NEAR(20.0, orthogonalizedCurvilinearGrid.m_gridNodes(1, 2).y, tolerance);
    ASSERT_NEAR(30.0, orthogonalizedCurvilinearGrid.m_gridNodes(1, 3).y, tolerance);

    ASSERT_NEAR(0.0, orthogonalizedCurvilinearGrid.m_gridNodes(2, 0).y, tolerance);
    ASSERT_NEAR(10.0, orthogonalizedCurvilinearGrid.m_gridNodes(2, 1).y, tolerance);
    ASSERT_NEAR(20.0, orthogonalizedCurvilinearGrid.m_gridNodes(2, 2).y, tolerance);
    ASSERT_NEAR(30.0, orthogonalizedCurvilinearGrid.m_gridNodes(2, 3).y, tolerance);

    ASSERT_NEAR(0.0, orthogonalizedCurvilinearGrid.m_gridNodes(3, 0).y, tolerance);
    ASSERT_NEAR(10.0, orthogonalizedCurvilinearGrid.m_gridNodes(3, 1).y, tolerance);
    ASSERT_NEAR(20.0, orthogonalizedCurvilinearGrid.m_gridNodes(3, 2).y, tolerance);
    ASSERT_NEAR(30.0, orthogonalizedCurvilinearGrid.m_gridNodes(3, 3).y, tolerance);
}

TEST(CurvilinearGridOrthogonalization, Compute_OnONonOrthogonalCurvilinearGrid_ShouldOrthogonalizeGrid)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();

    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    meshkernel::CurvilinearGridOrthogonalization curvilinearGridOrthogonalization(curvilinearGrid, orthogonalizationParameters);
    curvilinearGridOrthogonalization.SetBlock({80154, 366530}, {80610, 367407});
    // Execute
    auto const orthogonalizedCurvilinearGrid = curvilinearGridOrthogonalization.Compute();

    // Assert
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(79983.796374595549, orthogonalizedCurvilinearGrid.m_gridNodes(0, 0).x, tolerance);
    ASSERT_NEAR(80067.930920933941, orthogonalizedCurvilinearGrid.m_gridNodes(0, 1).x, tolerance);
    ASSERT_NEAR(80150.164245250286, orthogonalizedCurvilinearGrid.m_gridNodes(0, 2).x, tolerance);
    ASSERT_NEAR(80231.651392830841, orthogonalizedCurvilinearGrid.m_gridNodes(0, 3).x, tolerance);
    ASSERT_NEAR(80312.801762204603, orthogonalizedCurvilinearGrid.m_gridNodes(0, 4).x, tolerance);
    ASSERT_NEAR(80390.860703898958, orthogonalizedCurvilinearGrid.m_gridNodes(0, 5).x, tolerance);
    ASSERT_NEAR(80466.956631008565, orthogonalizedCurvilinearGrid.m_gridNodes(0, 6).x, tolerance);
    ASSERT_NEAR(80540.510527743099, orthogonalizedCurvilinearGrid.m_gridNodes(0, 7).x, tolerance);
    ASSERT_NEAR(80611.260660513333, orthogonalizedCurvilinearGrid.m_gridNodes(0, 8).x, tolerance);

    ASSERT_NEAR(80056.139417226193, orthogonalizedCurvilinearGrid.m_gridNodes(1, 0).x, tolerance);
    ASSERT_NEAR(80143.392341684361, orthogonalizedCurvilinearGrid.m_gridNodes(1, 1).x, tolerance);
    ASSERT_NEAR(80232.701973556381, orthogonalizedCurvilinearGrid.m_gridNodes(1, 2).x, tolerance);
    ASSERT_NEAR(80323.273756955823, orthogonalizedCurvilinearGrid.m_gridNodes(1, 3).x, tolerance);
    ASSERT_NEAR(80415.126141636414, orthogonalizedCurvilinearGrid.m_gridNodes(1, 4).x, tolerance);
    ASSERT_NEAR(80506.879117047443, orthogonalizedCurvilinearGrid.m_gridNodes(1, 5).x, tolerance);
    ASSERT_NEAR(80595.173874885382, orthogonalizedCurvilinearGrid.m_gridNodes(1, 6).x, tolerance);
    ASSERT_NEAR(80680.093103208215, orthogonalizedCurvilinearGrid.m_gridNodes(1, 7).x, tolerance);
    ASSERT_NEAR(80761.372120477128, orthogonalizedCurvilinearGrid.m_gridNodes(1, 8).x, tolerance);

    ASSERT_NEAR(366936.89538054139, orthogonalizedCurvilinearGrid.m_gridNodes(0, 0).y, tolerance);
    ASSERT_NEAR(366993.16520344437, orthogonalizedCurvilinearGrid.m_gridNodes(0, 1).y, tolerance);
    ASSERT_NEAR(367050.33461904334, orthogonalizedCurvilinearGrid.m_gridNodes(0, 2).y, tolerance);
    ASSERT_NEAR(367108.77222890750, orthogonalizedCurvilinearGrid.m_gridNodes(0, 3).y, tolerance);
    ASSERT_NEAR(367168.86295006215, orthogonalizedCurvilinearGrid.m_gridNodes(0, 4).y, tolerance);
    ASSERT_NEAR(367231.44573248079, orthogonalizedCurvilinearGrid.m_gridNodes(0, 5).y, tolerance);
    ASSERT_NEAR(367292.13654755038, orthogonalizedCurvilinearGrid.m_gridNodes(0, 6).y, tolerance);
    ASSERT_NEAR(367350.89140520856, orthogonalizedCurvilinearGrid.m_gridNodes(0, 7).y, tolerance);
    ASSERT_NEAR(367407.38541954994, orthogonalizedCurvilinearGrid.m_gridNodes(0, 8).y, tolerance);

    ASSERT_NEAR(366823.20595668437, orthogonalizedCurvilinearGrid.m_gridNodes(1, 0).y, tolerance);
    ASSERT_NEAR(366869.41001402331, orthogonalizedCurvilinearGrid.m_gridNodes(1, 1).y, tolerance);
    ASSERT_NEAR(366915.26601702173, orthogonalizedCurvilinearGrid.m_gridNodes(1, 2).y, tolerance);
    ASSERT_NEAR(366961.36595671845, orthogonalizedCurvilinearGrid.m_gridNodes(1, 3).y, tolerance);
    ASSERT_NEAR(367008.20690690872, orthogonalizedCurvilinearGrid.m_gridNodes(1, 4).y, tolerance);
    ASSERT_NEAR(367056.34885035310, orthogonalizedCurvilinearGrid.m_gridNodes(1, 5).y, tolerance);
    ASSERT_NEAR(367105.00611276925, orthogonalizedCurvilinearGrid.m_gridNodes(1, 6).y, tolerance);
    ASSERT_NEAR(367154.47604494903, orthogonalizedCurvilinearGrid.m_gridNodes(1, 7).y, tolerance);
    ASSERT_NEAR(367204.88319783867, orthogonalizedCurvilinearGrid.m_gridNodes(1, 8).y, tolerance);
}

TEST(CurvilinearGridOrthogonalization, Compute_OnONonOrthogonalCurvilinearGridWithMissingElements_ShouldOrthogonalizeGrid)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGridWithMissingFaces();

    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    meshkernel::CurvilinearGridOrthogonalization curvilinearGridOrthogonalization(curvilinearGrid, orthogonalizationParameters);
    curvilinearGridOrthogonalization.SetBlock({80154, 366530}, {80610, 367407});

    // Execute
    auto const orthogonalizedCurvilinearGrid = curvilinearGridOrthogonalization.Compute();

    // Assert
    constexpr double tolerance = 1e-6;

    EXPECT_NEAR(79983.796374595549, orthogonalizedCurvilinearGrid.m_gridNodes(0, 0).x, tolerance);
    EXPECT_NEAR(80069.224277354806, orthogonalizedCurvilinearGrid.m_gridNodes(0, 1).x, tolerance);
    EXPECT_NEAR(80152.817263131525, orthogonalizedCurvilinearGrid.m_gridNodes(0, 2).x, tolerance);
    EXPECT_NEAR(80235.609244143387, orthogonalizedCurvilinearGrid.m_gridNodes(0, 3).x, tolerance);
    EXPECT_NEAR(80317.466702245743, orthogonalizedCurvilinearGrid.m_gridNodes(0, 4).x, tolerance);
    EXPECT_NEAR(80395.149908970227, orthogonalizedCurvilinearGrid.m_gridNodes(0, 5).x, tolerance);
    EXPECT_NEAR(80470.021886679329, orthogonalizedCurvilinearGrid.m_gridNodes(0, 6).x, tolerance);
    EXPECT_NEAR(80542.049673235932, orthogonalizedCurvilinearGrid.m_gridNodes(0, 7).x, tolerance);
    EXPECT_NEAR(80611.260660513333, orthogonalizedCurvilinearGrid.m_gridNodes(0, 8).x, tolerance);

    EXPECT_NEAR(80055.199856352847, orthogonalizedCurvilinearGrid.m_gridNodes(1, 0).x, tolerance);
    EXPECT_NEAR(80143.721071059583, orthogonalizedCurvilinearGrid.m_gridNodes(1, 1).x, tolerance);
    EXPECT_NEAR(80234.545379411153, orthogonalizedCurvilinearGrid.m_gridNodes(1, 2).x, tolerance);
    EXPECT_NEAR(80326.293574499403, orthogonalizedCurvilinearGrid.m_gridNodes(1, 3).x, tolerance);
    EXPECT_NEAR(80418.809428865279, orthogonalizedCurvilinearGrid.m_gridNodes(1, 4).x, tolerance);
    EXPECT_NEAR(80510.299874655640, orthogonalizedCurvilinearGrid.m_gridNodes(1, 5).x, tolerance);
    EXPECT_NEAR(80597.036324757821, orthogonalizedCurvilinearGrid.m_gridNodes(1, 6).x, tolerance);
    EXPECT_NEAR(80680.228401363493, orthogonalizedCurvilinearGrid.m_gridNodes(1, 7).x, tolerance);
    EXPECT_NEAR(80759.822363775238, orthogonalizedCurvilinearGrid.m_gridNodes(1, 8).x, tolerance);

    EXPECT_NEAR(80104.514488846587, orthogonalizedCurvilinearGrid.m_gridNodes(2, 0).x, tolerance);
    EXPECT_NEAR(80197.383420613070, orthogonalizedCurvilinearGrid.m_gridNodes(2, 1).x, tolerance);
    EXPECT_NEAR(80292.449108019704, orthogonalizedCurvilinearGrid.m_gridNodes(2, 2).x, tolerance);
    EXPECT_NEAR(80387.175475807715, orthogonalizedCurvilinearGrid.m_gridNodes(2, 3).x, tolerance);
    EXPECT_NEAR(80480.608251576487, orthogonalizedCurvilinearGrid.m_gridNodes(2, 4).x, tolerance);
    EXPECT_NEAR(80583.120607369667, orthogonalizedCurvilinearGrid.m_gridNodes(2, 5).x, tolerance);
    EXPECT_NEAR(80682.512780448465, orthogonalizedCurvilinearGrid.m_gridNodes(2, 6).x, tolerance);
    EXPECT_NEAR(80780.894620879248, orthogonalizedCurvilinearGrid.m_gridNodes(2, 7).x, tolerance);
    EXPECT_NEAR(80877.580909293247, orthogonalizedCurvilinearGrid.m_gridNodes(2, 8).x, tolerance);

    EXPECT_NEAR(366936.89538054139, orthogonalizedCurvilinearGrid.m_gridNodes(0, 0).y, tolerance);
    EXPECT_NEAR(366994.04511051433, orthogonalizedCurvilinearGrid.m_gridNodes(0, 1).y, tolerance);
    EXPECT_NEAR(367052.21760805714, orthogonalizedCurvilinearGrid.m_gridNodes(0, 2).y, tolerance);
    EXPECT_NEAR(367111.62432093697, orthogonalizedCurvilinearGrid.m_gridNodes(0, 3).y, tolerance);
    EXPECT_NEAR(367172.50484630122, orthogonalizedCurvilinearGrid.m_gridNodes(0, 4).y, tolerance);
    EXPECT_NEAR(367234.88964184484, orthogonalizedCurvilinearGrid.m_gridNodes(0, 5).y, tolerance);
    EXPECT_NEAR(367294.58040378935, orthogonalizedCurvilinearGrid.m_gridNodes(0, 6).y, tolerance);
    EXPECT_NEAR(367352.12121038162, orthogonalizedCurvilinearGrid.m_gridNodes(0, 7).y, tolerance);
    EXPECT_NEAR(367407.38541954994, orthogonalizedCurvilinearGrid.m_gridNodes(0, 8).y, tolerance);

    EXPECT_NEAR(366824.95581170503, orthogonalizedCurvilinearGrid.m_gridNodes(1, 0).y, tolerance);
    EXPECT_NEAR(366871.97091034410, orthogonalizedCurvilinearGrid.m_gridNodes(1, 1).y, tolerance);
    EXPECT_NEAR(366918.90326042997, orthogonalizedCurvilinearGrid.m_gridNodes(1, 2).y, tolerance);
    EXPECT_NEAR(366966.74273245712, orthogonalizedCurvilinearGrid.m_gridNodes(1, 3).y, tolerance);
    EXPECT_NEAR(367014.95754932362, orthogonalizedCurvilinearGrid.m_gridNodes(1, 4).y, tolerance);
    EXPECT_NEAR(367062.00790467981, orthogonalizedCurvilinearGrid.m_gridNodes(1, 5).y, tolerance);
    EXPECT_NEAR(367109.50790550862, orthogonalizedCurvilinearGrid.m_gridNodes(1, 6).y, tolerance);
    EXPECT_NEAR(367157.80620957806, orthogonalizedCurvilinearGrid.m_gridNodes(1, 7).y, tolerance);
    EXPECT_NEAR(367207.16259613103, orthogonalizedCurvilinearGrid.m_gridNodes(1, 8).y, tolerance);

    EXPECT_NEAR(366718.00768776325, orthogonalizedCurvilinearGrid.m_gridNodes(2, 0).y, tolerance);
    EXPECT_NEAR(366755.34794646013, orthogonalizedCurvilinearGrid.m_gridNodes(2, 1).y, tolerance);
    EXPECT_NEAR(366792.50812354451, orthogonalizedCurvilinearGrid.m_gridNodes(2, 2).y, tolerance);
    EXPECT_NEAR(366832.52615748829, orthogonalizedCurvilinearGrid.m_gridNodes(2, 3).y, tolerance);
    EXPECT_NEAR(366870.46128228144, orthogonalizedCurvilinearGrid.m_gridNodes(2, 4).y, tolerance);
    EXPECT_NEAR(366891.62363194284, orthogonalizedCurvilinearGrid.m_gridNodes(2, 5).y, tolerance);
    EXPECT_NEAR(366923.14004067366, orthogonalizedCurvilinearGrid.m_gridNodes(2, 6).y, tolerance);
    EXPECT_NEAR(366957.61329611664, orthogonalizedCurvilinearGrid.m_gridNodes(2, 7).y, tolerance);
    EXPECT_NEAR(366996.07892524434, orthogonalizedCurvilinearGrid.m_gridNodes(2, 8).y, tolerance);
}

TEST(CurvilinearGridOrthogonalization, SetFrozenLine_OnONonOrthogonalGrid_WithCrossingFrozenLines_ShouldThrowAnStdException)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();

    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    meshkernel::CurvilinearGridOrthogonalization curvilinearGridOrthogonalization(curvilinearGrid, orthogonalizationParameters);
    curvilinearGridOrthogonalization.SetBlock({80154, 366530}, {80610, 367407});
    curvilinearGridOrthogonalization.SetLine({80144, 367046}, {80329, 366550});

    // Execute and assert
    ASSERT_THROW(curvilinearGridOrthogonalization.SetLine({80052, 366824}, {80774, 367186}), std::exception);
}

TEST(CurvilinearGridOrthogonalization, Compute_OnONonOrthogonalCurvilinearGridWithFrozenLines_ShouldOrthogonalizeGrid)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();

    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    meshkernel::CurvilinearGridOrthogonalization curvilinearGridOrthogonalization(curvilinearGrid, orthogonalizationParameters);
    curvilinearGridOrthogonalization.SetBlock({80154, 366530}, {80610, 367407});
    curvilinearGridOrthogonalization.SetLine({80144, 367046}, {80329, 366550});

    // Execute
    auto const orthogonalizedCurvilinearGrid = curvilinearGridOrthogonalization.Compute();

    // Assert
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(79983.796374595549, orthogonalizedCurvilinearGrid.m_gridNodes(0, 0).x, tolerance);
    ASSERT_NEAR(80069.479272425073, orthogonalizedCurvilinearGrid.m_gridNodes(0, 1).x, tolerance);
    ASSERT_NEAR(80153.235772058310, orthogonalizedCurvilinearGrid.m_gridNodes(0, 2).x, tolerance);
    ASSERT_NEAR(80234.211288098682, orthogonalizedCurvilinearGrid.m_gridNodes(0, 3).x, tolerance);
    ASSERT_NEAR(80314.661153602894, orthogonalizedCurvilinearGrid.m_gridNodes(0, 4).x, tolerance);
    ASSERT_NEAR(80392.150195938390, orthogonalizedCurvilinearGrid.m_gridNodes(0, 5).x, tolerance);
    ASSERT_NEAR(80467.771942028790, orthogonalizedCurvilinearGrid.m_gridNodes(0, 6).x, tolerance);
    ASSERT_NEAR(80540.898431866168, orthogonalizedCurvilinearGrid.m_gridNodes(0, 7).x, tolerance);
    ASSERT_NEAR(80611.260660513333, orthogonalizedCurvilinearGrid.m_gridNodes(0, 8).x, tolerance);

    ASSERT_NEAR(80055.196755132944, orthogonalizedCurvilinearGrid.m_gridNodes(1, 0).x, tolerance);
    ASSERT_NEAR(80143.960692327732, orthogonalizedCurvilinearGrid.m_gridNodes(1, 1).x, tolerance);
    ASSERT_NEAR(80234.924195000407, orthogonalizedCurvilinearGrid.m_gridNodes(1, 2).x, tolerance);
    ASSERT_NEAR(80324.849820521486, orthogonalizedCurvilinearGrid.m_gridNodes(1, 3).x, tolerance);
    ASSERT_NEAR(80416.240608373060, orthogonalizedCurvilinearGrid.m_gridNodes(1, 4).x, tolerance);
    ASSERT_NEAR(80507.596543850144, orthogonalizedCurvilinearGrid.m_gridNodes(1, 5).x, tolerance);
    ASSERT_NEAR(80595.526594976516, orthogonalizedCurvilinearGrid.m_gridNodes(1, 6).x, tolerance);
    ASSERT_NEAR(80680.087189144266, orthogonalizedCurvilinearGrid.m_gridNodes(1, 7).x, tolerance);
    ASSERT_NEAR(80760.998582117099, orthogonalizedCurvilinearGrid.m_gridNodes(1, 8).x, tolerance);

    ASSERT_NEAR(366936.89538054139, orthogonalizedCurvilinearGrid.m_gridNodes(0, 0).y, tolerance);
    ASSERT_NEAR(366994.21866791911, orthogonalizedCurvilinearGrid.m_gridNodes(0, 1).y, tolerance);
    ASSERT_NEAR(367052.51483841456, orthogonalizedCurvilinearGrid.m_gridNodes(0, 2).y, tolerance);
    ASSERT_NEAR(367110.61675055756, orthogonalizedCurvilinearGrid.m_gridNodes(0, 3).y, tolerance);
    ASSERT_NEAR(367170.31164987158, orthogonalizedCurvilinearGrid.m_gridNodes(0, 4).y, tolerance);
    ASSERT_NEAR(367232.48168405943, orthogonalizedCurvilinearGrid.m_gridNodes(0, 5).y, tolerance);
    ASSERT_NEAR(367292.78650072071, orthogonalizedCurvilinearGrid.m_gridNodes(0, 6).y, tolerance);
    ASSERT_NEAR(367351.20135266299, orthogonalizedCurvilinearGrid.m_gridNodes(0, 7).y, tolerance);
    ASSERT_NEAR(367407.38541954994, orthogonalizedCurvilinearGrid.m_gridNodes(0, 8).y, tolerance);

    ASSERT_NEAR(366824.96156769694, orthogonalizedCurvilinearGrid.m_gridNodes(1, 0).y, tolerance);
    ASSERT_NEAR(366872.12740536616, orthogonalizedCurvilinearGrid.m_gridNodes(1, 1).y, tolerance);
    ASSERT_NEAR(366919.18816119258, orthogonalizedCurvilinearGrid.m_gridNodes(1, 2).y, tolerance);
    ASSERT_NEAR(366964.42115637776, orthogonalizedCurvilinearGrid.m_gridNodes(1, 3).y, tolerance);
    ASSERT_NEAR(367010.41102564143, orthogonalizedCurvilinearGrid.m_gridNodes(1, 4).y, tolerance);
    ASSERT_NEAR(367057.91817534604, orthogonalizedCurvilinearGrid.m_gridNodes(1, 5).y, tolerance);
    ASSERT_NEAR(367106.12547266396, orthogonalizedCurvilinearGrid.m_gridNodes(1, 6).y, tolerance);
    ASSERT_NEAR(367155.26906854485, orthogonalizedCurvilinearGrid.m_gridNodes(1, 7).y, tolerance);
    ASSERT_NEAR(367205.43327878905, orthogonalizedCurvilinearGrid.m_gridNodes(1, 8).y, tolerance);
}
