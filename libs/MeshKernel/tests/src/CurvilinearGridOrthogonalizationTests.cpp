//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2025.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation version 3.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// contact: delft3d.support@deltares.nl
// Stichting Deltares
// P.O. Box 177
// 2600 MH Delft, The Netherlands
//
// All indications and logos of, and references to, "Delft3D" and "Deltares"
// are registered trademarks of Stichting Deltares, and remain the property of
// Stichting Deltares. All rights reserved.
//
//------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridOrthogonalization.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Parameters.hpp>
#include <MeshKernel/Utilities/LinearAlgebra.hpp>
#include <MeshKernel/Utilities/Utilities.hpp>
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

    CurvilinearGrid curvilinearGrid(grid, meshkernel::Projection::cartesian);

    // Move a node, to make the grid strongly non orthogonal
    [[maybe_unused]] auto dummyUnusedAction = curvilinearGrid.MoveNode(meshkernel::Point(10.0, 20.0), meshkernel::Point(18.0, 12.0));

    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 1;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    meshkernel::CurvilinearGridOrthogonalization curvilinearGridOrthogonalization(curvilinearGrid, orthogonalizationParameters);
    curvilinearGridOrthogonalization.SetBlock({0, 0}, {30, 30});

    // Execute
    [[maybe_unused]] auto dummyUndoAction = curvilinearGridOrthogonalization.Compute();

    // Assert the moved nodes has moved towards its original location, making the grid more orthogonal
    constexpr double tolerance = 1e-6;
    ASSERT_NEAR(11.841396536135521, curvilinearGrid.GetNode(1, 2).x, tolerance);
    ASSERT_NEAR(18.158586078094562, curvilinearGrid.GetNode(1, 2).y, tolerance);
}

TEST(CurvilinearGridOrthogonalization, Compute_OnOrthogonalCurvilinearGrid_ShouldNotModifyGrid)
{
    // Set-up
    lin_alg::Matrix<Point> grid(5, 5);
    grid << Point{0, 0}, Point{0, 10}, Point{0, 20}, Point{0, 30}, Point{0, 40},
        Point{10, 0}, Point{10, 10}, Point{10, 20}, Point{10, 30}, Point{10, 40},
        Point{20, 0}, Point{20, 10}, Point{20, 20}, Point{20, 30}, Point{20, 40},
        Point{30, 0}, Point{30, 10}, Point{30, 20}, Point{30, 30}, Point{30, 40},
        Point{40, 0}, Point{40, 10}, Point{40, 20}, Point{40, 30}, Point{40, 40};

    meshkernel::CurvilinearGrid curvilinearGrid(grid, meshkernel::Projection::cartesian);

    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 1;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    meshkernel::CurvilinearGridOrthogonalization curvilinearGridOrthogonalization(curvilinearGrid, orthogonalizationParameters);
    curvilinearGridOrthogonalization.SetBlock({0, 0}, {30, 30});
    curvilinearGridOrthogonalization.SetLine({20.0, 0.0}, {20.0, 30.0});

    // Execute
    [[maybe_unused]] auto dummyUndoAction = curvilinearGridOrthogonalization.Compute();

    // Assert nodes are on the same location because the grid is already orthogonal
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(0.0, curvilinearGrid.GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(0.0, curvilinearGrid.GetNode(0, 1).x, tolerance);
    ASSERT_NEAR(0.0, curvilinearGrid.GetNode(0, 2).x, tolerance);
    ASSERT_NEAR(0.0, curvilinearGrid.GetNode(0, 3).x, tolerance);
    ASSERT_NEAR(0.0, curvilinearGrid.GetNode(0, 4).x, tolerance);

    ASSERT_NEAR(10.0, curvilinearGrid.GetNode(1, 0).x, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid.GetNode(1, 1).x, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid.GetNode(1, 2).x, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid.GetNode(1, 3).x, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid.GetNode(1, 4).x, tolerance);

    ASSERT_NEAR(20.0, curvilinearGrid.GetNode(2, 0).x, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid.GetNode(2, 1).x, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid.GetNode(2, 2).x, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid.GetNode(2, 3).x, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid.GetNode(2, 4).x, tolerance);

    ASSERT_NEAR(30.0, curvilinearGrid.GetNode(3, 0).x, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid.GetNode(3, 1).x, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid.GetNode(3, 2).x, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid.GetNode(3, 3).x, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid.GetNode(3, 4).x, tolerance);

    ASSERT_NEAR(0.0, curvilinearGrid.GetNode(0, 0).y, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid.GetNode(0, 1).y, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid.GetNode(0, 2).y, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid.GetNode(0, 3).y, tolerance);
    ASSERT_NEAR(40.0, curvilinearGrid.GetNode(0, 4).y, tolerance);

    ASSERT_NEAR(0.0, curvilinearGrid.GetNode(1, 0).y, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid.GetNode(1, 1).y, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid.GetNode(1, 2).y, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid.GetNode(1, 3).y, tolerance);
    ASSERT_NEAR(40.0, curvilinearGrid.GetNode(1, 4).y, tolerance);

    ASSERT_NEAR(0.0, curvilinearGrid.GetNode(2, 0).y, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid.GetNode(2, 1).y, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid.GetNode(2, 2).y, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid.GetNode(2, 3).y, tolerance);
    ASSERT_NEAR(40.0, curvilinearGrid.GetNode(2, 4).y, tolerance);

    ASSERT_NEAR(0.0, curvilinearGrid.GetNode(3, 0).y, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid.GetNode(3, 1).y, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid.GetNode(3, 2).y, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid.GetNode(3, 3).y, tolerance);
    ASSERT_NEAR(40.0, curvilinearGrid.GetNode(3, 4).y, tolerance);
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
    meshkernel::CurvilinearGridOrthogonalization curvilinearGridOrthogonalization(*curvilinearGrid, orthogonalizationParameters);
    curvilinearGridOrthogonalization.SetBlock({80154, 366530}, {80610, 367407});
    // Execute
    [[maybe_unused]] auto dummyUndoAction = curvilinearGridOrthogonalization.Compute();

    // Assert
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(79983.796374595549, curvilinearGrid->GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(80067.930920933941, curvilinearGrid->GetNode(0, 1).x, tolerance);
    ASSERT_NEAR(80150.164245250286, curvilinearGrid->GetNode(0, 2).x, tolerance);
    ASSERT_NEAR(80231.651392830841, curvilinearGrid->GetNode(0, 3).x, tolerance);
    ASSERT_NEAR(80312.801762204603, curvilinearGrid->GetNode(0, 4).x, tolerance);
    ASSERT_NEAR(80390.860703898958, curvilinearGrid->GetNode(0, 5).x, tolerance);
    ASSERT_NEAR(80466.956631008565, curvilinearGrid->GetNode(0, 6).x, tolerance);
    ASSERT_NEAR(80540.510527743099, curvilinearGrid->GetNode(0, 7).x, tolerance);
    ASSERT_NEAR(80611.260660513333, curvilinearGrid->GetNode(0, 8).x, tolerance);

    ASSERT_NEAR(80056.139417226193, curvilinearGrid->GetNode(1, 0).x, tolerance);
    ASSERT_NEAR(80143.392341684361, curvilinearGrid->GetNode(1, 1).x, tolerance);
    ASSERT_NEAR(80232.701973556381, curvilinearGrid->GetNode(1, 2).x, tolerance);
    ASSERT_NEAR(80323.273756955823, curvilinearGrid->GetNode(1, 3).x, tolerance);
    ASSERT_NEAR(80415.126141636414, curvilinearGrid->GetNode(1, 4).x, tolerance);
    ASSERT_NEAR(80506.879117047443, curvilinearGrid->GetNode(1, 5).x, tolerance);
    ASSERT_NEAR(80595.173874885382, curvilinearGrid->GetNode(1, 6).x, tolerance);
    ASSERT_NEAR(80680.093103208215, curvilinearGrid->GetNode(1, 7).x, tolerance);
    ASSERT_NEAR(80761.372120477128, curvilinearGrid->GetNode(1, 8).x, tolerance);

    ASSERT_NEAR(366936.89538054139, curvilinearGrid->GetNode(0, 0).y, tolerance);
    ASSERT_NEAR(366993.16520344437, curvilinearGrid->GetNode(0, 1).y, tolerance);
    ASSERT_NEAR(367050.33461904334, curvilinearGrid->GetNode(0, 2).y, tolerance);
    ASSERT_NEAR(367108.77222890750, curvilinearGrid->GetNode(0, 3).y, tolerance);
    ASSERT_NEAR(367168.86295006215, curvilinearGrid->GetNode(0, 4).y, tolerance);
    ASSERT_NEAR(367231.44573248079, curvilinearGrid->GetNode(0, 5).y, tolerance);
    ASSERT_NEAR(367292.13654755038, curvilinearGrid->GetNode(0, 6).y, tolerance);
    ASSERT_NEAR(367350.89140520856, curvilinearGrid->GetNode(0, 7).y, tolerance);
    ASSERT_NEAR(367407.38541954994, curvilinearGrid->GetNode(0, 8).y, tolerance);

    ASSERT_NEAR(366823.20595668437, curvilinearGrid->GetNode(1, 0).y, tolerance);
    ASSERT_NEAR(366869.41001402331, curvilinearGrid->GetNode(1, 1).y, tolerance);
    ASSERT_NEAR(366915.26601702173, curvilinearGrid->GetNode(1, 2).y, tolerance);
    ASSERT_NEAR(366961.36595671845, curvilinearGrid->GetNode(1, 3).y, tolerance);
    ASSERT_NEAR(367008.20690690872, curvilinearGrid->GetNode(1, 4).y, tolerance);
    ASSERT_NEAR(367056.34885035310, curvilinearGrid->GetNode(1, 5).y, tolerance);
    ASSERT_NEAR(367105.00611276925, curvilinearGrid->GetNode(1, 6).y, tolerance);
    ASSERT_NEAR(367154.47604494903, curvilinearGrid->GetNode(1, 7).y, tolerance);
    ASSERT_NEAR(367204.88319783867, curvilinearGrid->GetNode(1, 8).y, tolerance);
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
    meshkernel::CurvilinearGridOrthogonalization curvilinearGridOrthogonalization(*curvilinearGrid, orthogonalizationParameters);
    curvilinearGridOrthogonalization.SetBlock({80154, 366530}, {80610, 367407});

    // Execute
    [[maybe_unused]] auto dummyUndoAction = curvilinearGridOrthogonalization.Compute();

    // Assert
    constexpr double tolerance = 1e-6;
    EXPECT_NEAR(79983.796374595549, curvilinearGrid->GetNode(0, 0).x, tolerance);
    EXPECT_NEAR(80069.224277354806, curvilinearGrid->GetNode(0, 1).x, tolerance);
    EXPECT_NEAR(80152.817263131525, curvilinearGrid->GetNode(0, 2).x, tolerance);
    EXPECT_NEAR(80235.609244143387, curvilinearGrid->GetNode(0, 3).x, tolerance);
    EXPECT_NEAR(80317.466702245743, curvilinearGrid->GetNode(0, 4).x, tolerance);
    EXPECT_NEAR(80395.149908970227, curvilinearGrid->GetNode(0, 5).x, tolerance);
    EXPECT_NEAR(80470.021886679329, curvilinearGrid->GetNode(0, 6).x, tolerance);
    EXPECT_NEAR(80542.049673235932, curvilinearGrid->GetNode(0, 7).x, tolerance);
    EXPECT_NEAR(80611.260660513333, curvilinearGrid->GetNode(0, 8).x, tolerance);

    EXPECT_NEAR(80055.199856352847, curvilinearGrid->GetNode(1, 0).x, tolerance);
    EXPECT_NEAR(80143.721071059583, curvilinearGrid->GetNode(1, 1).x, tolerance);
    EXPECT_NEAR(80234.545379411153, curvilinearGrid->GetNode(1, 2).x, tolerance);
    EXPECT_NEAR(80326.293574499403, curvilinearGrid->GetNode(1, 3).x, tolerance);
    EXPECT_NEAR(80418.809428865279, curvilinearGrid->GetNode(1, 4).x, tolerance);
    EXPECT_NEAR(80510.299874655640, curvilinearGrid->GetNode(1, 5).x, tolerance);
    EXPECT_NEAR(80597.036324757821, curvilinearGrid->GetNode(1, 6).x, tolerance);
    EXPECT_NEAR(80680.228401363493, curvilinearGrid->GetNode(1, 7).x, tolerance);
    EXPECT_NEAR(80759.822363775238, curvilinearGrid->GetNode(1, 8).x, tolerance);

    EXPECT_NEAR(80104.514488846587, curvilinearGrid->GetNode(2, 0).x, tolerance);
    EXPECT_NEAR(80197.383420613070, curvilinearGrid->GetNode(2, 1).x, tolerance);
    EXPECT_NEAR(80292.449108019704, curvilinearGrid->GetNode(2, 2).x, tolerance);
    EXPECT_NEAR(80387.175475807715, curvilinearGrid->GetNode(2, 3).x, tolerance);
    EXPECT_NEAR(80480.608251576487, curvilinearGrid->GetNode(2, 4).x, tolerance);
    EXPECT_NEAR(80583.120607369667, curvilinearGrid->GetNode(2, 5).x, tolerance);
    EXPECT_NEAR(80682.512780448465, curvilinearGrid->GetNode(2, 6).x, tolerance);
    EXPECT_NEAR(80780.894620879248, curvilinearGrid->GetNode(2, 7).x, tolerance);
    EXPECT_NEAR(80877.580909293247, curvilinearGrid->GetNode(2, 8).x, tolerance);

    EXPECT_NEAR(366936.89538054139, curvilinearGrid->GetNode(0, 0).y, tolerance);
    EXPECT_NEAR(366994.04511051433, curvilinearGrid->GetNode(0, 1).y, tolerance);
    EXPECT_NEAR(367052.21760805714, curvilinearGrid->GetNode(0, 2).y, tolerance);
    EXPECT_NEAR(367111.62432093697, curvilinearGrid->GetNode(0, 3).y, tolerance);
    EXPECT_NEAR(367172.50484630122, curvilinearGrid->GetNode(0, 4).y, tolerance);
    EXPECT_NEAR(367234.88964184484, curvilinearGrid->GetNode(0, 5).y, tolerance);
    EXPECT_NEAR(367294.58040378935, curvilinearGrid->GetNode(0, 6).y, tolerance);
    EXPECT_NEAR(367352.12121038162, curvilinearGrid->GetNode(0, 7).y, tolerance);
    EXPECT_NEAR(367407.38541954994, curvilinearGrid->GetNode(0, 8).y, tolerance);

    EXPECT_NEAR(366824.95581170503, curvilinearGrid->GetNode(1, 0).y, tolerance);
    EXPECT_NEAR(366871.97091034410, curvilinearGrid->GetNode(1, 1).y, tolerance);
    EXPECT_NEAR(366918.90326042997, curvilinearGrid->GetNode(1, 2).y, tolerance);
    EXPECT_NEAR(366966.74273245712, curvilinearGrid->GetNode(1, 3).y, tolerance);
    EXPECT_NEAR(367014.95754932362, curvilinearGrid->GetNode(1, 4).y, tolerance);
    EXPECT_NEAR(367062.00790467981, curvilinearGrid->GetNode(1, 5).y, tolerance);
    EXPECT_NEAR(367109.50790550862, curvilinearGrid->GetNode(1, 6).y, tolerance);
    EXPECT_NEAR(367157.80620957806, curvilinearGrid->GetNode(1, 7).y, tolerance);
    EXPECT_NEAR(367207.16259613103, curvilinearGrid->GetNode(1, 8).y, tolerance);

    EXPECT_NEAR(366718.00768776325, curvilinearGrid->GetNode(2, 0).y, tolerance);
    EXPECT_NEAR(366755.34794646013, curvilinearGrid->GetNode(2, 1).y, tolerance);
    EXPECT_NEAR(366792.50812354451, curvilinearGrid->GetNode(2, 2).y, tolerance);
    EXPECT_NEAR(366832.52615748829, curvilinearGrid->GetNode(2, 3).y, tolerance);
    EXPECT_NEAR(366870.46128228144, curvilinearGrid->GetNode(2, 4).y, tolerance);
    EXPECT_NEAR(366891.62363194284, curvilinearGrid->GetNode(2, 5).y, tolerance);
    EXPECT_NEAR(366923.14004067366, curvilinearGrid->GetNode(2, 6).y, tolerance);
    EXPECT_NEAR(366957.61329611664, curvilinearGrid->GetNode(2, 7).y, tolerance);
    EXPECT_NEAR(366996.07892524434, curvilinearGrid->GetNode(2, 8).y, tolerance);
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
    meshkernel::CurvilinearGridOrthogonalization curvilinearGridOrthogonalization(*curvilinearGrid, orthogonalizationParameters);
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
    meshkernel::CurvilinearGridOrthogonalization curvilinearGridOrthogonalization(*curvilinearGrid, orthogonalizationParameters);
    curvilinearGridOrthogonalization.SetBlock({80154, 366530}, {80610, 367407});
    curvilinearGridOrthogonalization.SetLine({80144, 367046}, {80329, 366550});

    // Execute
    [[maybe_unused]] auto dummyUndoAction = curvilinearGridOrthogonalization.Compute();

    // Assert
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(79983.796374595549, curvilinearGrid->GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(80069.479272425073, curvilinearGrid->GetNode(0, 1).x, tolerance);
    ASSERT_NEAR(80153.235772058310, curvilinearGrid->GetNode(0, 2).x, tolerance);
    ASSERT_NEAR(80234.211288098682, curvilinearGrid->GetNode(0, 3).x, tolerance);
    ASSERT_NEAR(80314.661153602894, curvilinearGrid->GetNode(0, 4).x, tolerance);
    ASSERT_NEAR(80392.150195938390, curvilinearGrid->GetNode(0, 5).x, tolerance);
    ASSERT_NEAR(80467.771942028790, curvilinearGrid->GetNode(0, 6).x, tolerance);
    ASSERT_NEAR(80540.898431866168, curvilinearGrid->GetNode(0, 7).x, tolerance);
    ASSERT_NEAR(80611.260660513333, curvilinearGrid->GetNode(0, 8).x, tolerance);

    ASSERT_NEAR(80055.196755132944, curvilinearGrid->GetNode(1, 0).x, tolerance);
    ASSERT_NEAR(80143.960692327732, curvilinearGrid->GetNode(1, 1).x, tolerance);
    ASSERT_NEAR(80234.924195000407, curvilinearGrid->GetNode(1, 2).x, tolerance);
    ASSERT_NEAR(80324.849820521486, curvilinearGrid->GetNode(1, 3).x, tolerance);
    ASSERT_NEAR(80416.240608373060, curvilinearGrid->GetNode(1, 4).x, tolerance);
    ASSERT_NEAR(80507.596543850144, curvilinearGrid->GetNode(1, 5).x, tolerance);
    ASSERT_NEAR(80595.526594976516, curvilinearGrid->GetNode(1, 6).x, tolerance);
    ASSERT_NEAR(80680.087189144266, curvilinearGrid->GetNode(1, 7).x, tolerance);
    ASSERT_NEAR(80760.998582117099, curvilinearGrid->GetNode(1, 8).x, tolerance);

    ASSERT_NEAR(366936.89538054139, curvilinearGrid->GetNode(0, 0).y, tolerance);
    ASSERT_NEAR(366994.21866791911, curvilinearGrid->GetNode(0, 1).y, tolerance);
    ASSERT_NEAR(367052.51483841456, curvilinearGrid->GetNode(0, 2).y, tolerance);
    ASSERT_NEAR(367110.61675055756, curvilinearGrid->GetNode(0, 3).y, tolerance);
    ASSERT_NEAR(367170.31164987158, curvilinearGrid->GetNode(0, 4).y, tolerance);
    ASSERT_NEAR(367232.48168405943, curvilinearGrid->GetNode(0, 5).y, tolerance);
    ASSERT_NEAR(367292.78650072071, curvilinearGrid->GetNode(0, 6).y, tolerance);
    ASSERT_NEAR(367351.20135266299, curvilinearGrid->GetNode(0, 7).y, tolerance);
    ASSERT_NEAR(367407.38541954994, curvilinearGrid->GetNode(0, 8).y, tolerance);

    ASSERT_NEAR(366824.96156769694, curvilinearGrid->GetNode(1, 0).y, tolerance);
    ASSERT_NEAR(366872.12740536616, curvilinearGrid->GetNode(1, 1).y, tolerance);
    ASSERT_NEAR(366919.18816119258, curvilinearGrid->GetNode(1, 2).y, tolerance);
    ASSERT_NEAR(366964.42115637776, curvilinearGrid->GetNode(1, 3).y, tolerance);
    ASSERT_NEAR(367010.41102564143, curvilinearGrid->GetNode(1, 4).y, tolerance);
    ASSERT_NEAR(367057.91817534604, curvilinearGrid->GetNode(1, 5).y, tolerance);
    ASSERT_NEAR(367106.12547266396, curvilinearGrid->GetNode(1, 6).y, tolerance);
    ASSERT_NEAR(367155.26906854485, curvilinearGrid->GetNode(1, 7).y, tolerance);
    ASSERT_NEAR(367205.43327878905, curvilinearGrid->GetNode(1, 8).y, tolerance);
}

TEST(CurvilinearGridOrthogonalization, SetFrozenLine_ShouldFreezeLines)
{
    double deltaX = 10.0;
    double deltaY = 10.0;

    size_t sizeX = 15;
    size_t sizeY = 15;

    // Set-up
    const auto curvilinearGrid = MakeCurvilinearGridRand(0.0, 0.0, deltaX, deltaY, sizeX, sizeY, 0.7, true);

    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 5;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;

    meshkernel::CurvilinearGridOrthogonalization orthogonalisation(*curvilinearGrid, orthogonalizationParameters);

    Point blockLL = curvilinearGrid->GetNode(2, 2);
    Point blockUR = curvilinearGrid->GetNode(12, 12);

    const UInt line1IndexY = 10;
    const UInt line1StartIndex = 0;
    const UInt line1EndIndex = 14;

    Point line1Start = curvilinearGrid->GetNode(line1StartIndex, line1IndexY);
    Point line1End = curvilinearGrid->GetNode(line1EndIndex, line1IndexY);

    const UInt line2IndexY = 10;
    const UInt line2StartIndex = 4;
    const UInt line2EndIndex = 7;

    Point line2Start = curvilinearGrid->GetNode(line2StartIndex, line2IndexY);
    Point line2End = curvilinearGrid->GetNode(line2EndIndex, line2IndexY);

    std::vector<Point> originalLine1Points(line1EndIndex - line1StartIndex + 1);
    std::vector<Point> originalLine2Points(line2EndIndex - line2StartIndex + 1);

    // Collect line points before orthogonalising
    for (UInt i = line1StartIndex; i < line1EndIndex + 1; ++i)
    {
        Point p = curvilinearGrid->GetNode(i, line1IndexY);
        originalLine1Points[i - line1StartIndex] = p;
    }

    for (UInt i = line2StartIndex; i < line2EndIndex + 1; ++i)
    {
        Point p = curvilinearGrid->GetNode(i, line2IndexY);
        originalLine2Points[i - line2StartIndex] = p;
    }

    orthogonalisation.SetLine(line1Start, line1End);
    orthogonalisation.SetLine(line2Start, line2End);

    orthogonalisation.SetBlock(blockLL, blockUR);
    [[maybe_unused]] auto undo = orthogonalisation.Compute();

    const double tolerance = 1.0e-10;

    for (UInt i = line2StartIndex; i < line2EndIndex + 1; ++i)
    {
        EXPECT_NEAR(originalLine2Points[i - line2StartIndex].x, curvilinearGrid->GetNode(i, line2IndexY).x, tolerance);
        EXPECT_NEAR(originalLine2Points[i - line2StartIndex].y, curvilinearGrid->GetNode(i, line2IndexY).y, tolerance);
    }

    for (UInt i = line1StartIndex; i < line1EndIndex + 1; ++i)
    {
        EXPECT_NEAR(originalLine1Points[i].x, curvilinearGrid->GetNode(i, line1IndexY).x, tolerance);
        EXPECT_NEAR(originalLine1Points[i].y, curvilinearGrid->GetNode(i, line1IndexY).y, tolerance);
    }
}

TEST(CurvilinearGridOrthogonalization, Compute_CurvilinearGrid_ShouldOrthogonaliseTopAndRight)
{
    const std::vector<double> random{-0.368462, -0.0413499, -0.281041, 0.178865, 0.434693, 0.0194164,
                                     -0.465428, 0.0297002, -0.492302, -0.433158, 0.186773, 0.430436,
                                     0.0269288, 0.153919, 0.201191, 0.262198, -0.452535, -0.171766,
                                     0.25641, -0.134661, 0.48255, 0.253356, -0.427314, 0.384707,
                                     -0.0635886, -0.0222682, -0.225093, -0.333493, 0.397656, -0.439436,
                                     0.00452289, -0.180967, -0.00602331, -0.409267, -0.426251, -0.115858,
                                     0.413817, -0.0355542, -0.449916, 0.270205, -0.374635, 0.188455, 0.129543};

    double deltaX = 10.0;
    double deltaY = 10.0;

    size_t sizeX = 15;
    size_t sizeY = 15;

    // Set-up
    const auto curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, deltaX, deltaY, sizeX, sizeY);

    // displace grid
    size_t randomCounter = 0;

    auto randomCount = [random, &randomCounter]() mutable
    {
        if (randomCounter == random.size() - 1)
        {
            randomCounter = 0;
        }
        else
        {
            ++randomCounter;
        }

        return randomCounter;
    };

    for (UInt i = 0; i < curvilinearGrid->NumN(); ++i)
    {
        for (UInt j = 0; j < curvilinearGrid->NumM(); ++j)
        {
            double xDisplacement = random[randomCount()] * deltaX;
            double yDisplacement = random[randomCount()] * deltaY;

            curvilinearGrid->GetNode(i, j) += meshkernel::Vector(xDisplacement, yDisplacement);
        }
    }

    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 5;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;

    meshkernel::CurvilinearGridOrthogonalization orthogonalisation(*curvilinearGrid, orthogonalizationParameters);

    UInt bottomLeftIndex = 2;
    UInt topRightIndex = 12;

    Point blockLL = curvilinearGrid->GetNode(bottomLeftIndex, bottomLeftIndex);
    Point blockUR = curvilinearGrid->GetNode(topRightIndex, topRightIndex);

    // Values are not obtained analytically.
    std::vector<double> expectedRightLineX{119.777318, 120.269288, 120.476098807693,
                                           120.167692063759, 120.888155557653, 121.097530069614,
                                           120.704047794346, 119.647257737105, 120.042828014369,
                                           121.756316095996, 122.040604579151, 123.995906412066,
                                           123.031463355581, 123.97656, 122.62198};

    std::vector<double> expectedRightLineY{-2.25093, 11.53919, 17.4214056392148,
                                           29.6912639017348, 39.1864645374443, 45.9879848812012,
                                           62.5038337578463, 73.1581215544347, 81.686469185601,
                                           88.2137897170903, 97.9273961886399, 111.063380173574,
                                           119.272887605696, 125.60564, 135.47465};

    std::vector<double> expectedTopLineX{-1.71766, 8.65339, 16.3087849390702,
                                         31.6757555235772, 39.3153282473397, 50.0873598308451,
                                         56.1169500200071, 72.4679819996561, 83.6122397705629,
                                         86.6085676466389, 98.0347406659406, 107.331890054778,
                                         123.031463355581, 126.31538, 137.18959};

    std::vector<double> expectedTopLineY{122.5641, 124.8255, 122.64834412699,
                                         124.925222917871, 123.169763072201, 121.954057676841,
                                         121.115177221152, 120.909806224416, 120.695014213092,
                                         120.61516654566, 120.182359945175, 119.580290123743,
                                         119.272887605696, 119.586501, 121.78865};

    orthogonalisation.SetBlock(blockLL, blockUR);
    [[maybe_unused]] auto undo = orthogonalisation.Compute();

    const double tolerance = 1.0e-8;

    std::cout.precision(15);

    for (UInt i = 0; i < sizeY; ++i)
    {
        Point p = curvilinearGrid->GetNode(i, topRightIndex);
        EXPECT_NEAR(expectedRightLineX[i], p.x, tolerance);
        EXPECT_NEAR(expectedRightLineY[i], p.y, tolerance);
    }

    for (UInt i = 0; i < sizeX; ++i)
    {
        Point p = curvilinearGrid->GetNode(topRightIndex, i);
        EXPECT_NEAR(expectedTopLineX[i], p.x, tolerance);
        EXPECT_NEAR(expectedTopLineY[i], p.y, tolerance);
    }
}
