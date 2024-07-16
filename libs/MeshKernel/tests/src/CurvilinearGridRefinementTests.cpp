#include <gtest/gtest.h>

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridFullRefinement.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridRefinement.hpp>
#include <MeshKernel/Entities.hpp>

using namespace meshkernel;

TEST(CurvilinearGridRefinement, Compute_OnCurvilinearGrid_ShouldRefine)
{
    // Set-up
    lin_alg::Matrix<Point> grid(4, 4);

    grid << Point{0, 0}, Point{0, 10}, Point{0, 20}, Point{0, 30},
        Point{10, 0}, Point{10, 10}, Point{10, 20}, Point{10, 30},
        Point{20, 0}, Point{20, 10}, Point{20, 20}, Point{20, 30},
        Point{30, 0}, Point{30, 10}, Point{30, 20}, Point{30, 30};

    CurvilinearGrid curvilinearGrid(grid, Projection::cartesian);
    CurvilinearGridRefinement curvilinearGridRefinement(curvilinearGrid, 10);
    curvilinearGridRefinement.SetBlock({10, 20}, {20, 20});

    // Execute
    [[maybe_unused]] auto dummyUndoAction = curvilinearGridRefinement.Compute();

    // Assert
    ASSERT_EQ(4, curvilinearGrid.NumM());
    ASSERT_EQ(13, curvilinearGrid.NumN());

    constexpr double tolerance = 1e-12;
    ASSERT_NEAR(0.0, curvilinearGrid.GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid.GetNode(1, 0).x, tolerance);

    ASSERT_NEAR(11.0, curvilinearGrid.GetNode(2, 0).x, tolerance);
    ASSERT_NEAR(12.0, curvilinearGrid.GetNode(3, 0).x, tolerance);
    ASSERT_NEAR(13.0, curvilinearGrid.GetNode(4, 0).x, tolerance);
    ASSERT_NEAR(14.0, curvilinearGrid.GetNode(5, 0).x, tolerance);
    ASSERT_NEAR(15.0, curvilinearGrid.GetNode(6, 0).x, tolerance);
    ASSERT_NEAR(16.0, curvilinearGrid.GetNode(7, 0).x, tolerance);
    ASSERT_NEAR(17.0, curvilinearGrid.GetNode(8, 0).x, tolerance);
    ASSERT_NEAR(18.0, curvilinearGrid.GetNode(9, 0).x, tolerance);
    ASSERT_NEAR(19.0, curvilinearGrid.GetNode(10, 0).x, tolerance);

    ASSERT_NEAR(20.0, curvilinearGrid.GetNode(11, 0).x, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid.GetNode(12, 0).x, tolerance);

    ASSERT_NEAR(0.0, curvilinearGrid.GetNode(0, 0).y, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid.GetNode(0, 1).y, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid.GetNode(0, 2).y, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid.GetNode(0, 3).y, tolerance);
}

TEST(CurvilinearGridRefinement, Compute_OnCurvilinearGridWithMissingFaces_ShouldRefine)
{
    // Set-up
    lin_alg::Matrix<Point> grid(6, 4);

    grid << Point{0, 0}, Point{0, 10}, Point{0, 20}, Point{0, 30},
        Point{10, 0}, Point{10, 10}, Point{10, 20}, Point{10, 30},
        Point{}, Point{}, Point{20, 20}, Point{20, 30},
        Point{}, Point{}, Point{30, 20}, Point{30, 30},
        Point{40, 0}, Point{40, 10}, Point{40, 20}, Point{40, 30},
        Point{50, 0}, Point{50, 10}, Point{50, 20}, Point{50, 30};

    CurvilinearGrid curvilinearGrid(grid, Projection::cartesian);
    CurvilinearGridRefinement curvilinearGridRefinement(curvilinearGrid, 10);
    curvilinearGridRefinement.SetBlock({10, 20}, {20, 20});

    // Execute
    [[maybe_unused]] auto dummyUndoAction = curvilinearGridRefinement.Compute();

    // Assert
    ASSERT_EQ(4, curvilinearGrid.NumM());
    ASSERT_EQ(15, curvilinearGrid.NumN());

    constexpr double tolerance = 1e-12;

    // vertical gridline 0
    ASSERT_NEAR(0.0, curvilinearGrid.GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(0.0, curvilinearGrid.GetNode(0, 1).x, tolerance);
    ASSERT_NEAR(0.0, curvilinearGrid.GetNode(0, 2).x, tolerance);
    ASSERT_NEAR(0.0, curvilinearGrid.GetNode(0, 3).x, tolerance);

    ASSERT_NEAR(0.0, curvilinearGrid.GetNode(0, 0).y, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid.GetNode(0, 1).y, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid.GetNode(0, 2).y, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid.GetNode(0, 3).y, tolerance);

    // vertical gridline 2
    ASSERT_NEAR(constants::missing::doubleValue, curvilinearGrid.GetNode(2, 0).x, tolerance);
    ASSERT_NEAR(constants::missing::doubleValue, curvilinearGrid.GetNode(2, 1).x, tolerance);
    ASSERT_NEAR(11.0, curvilinearGrid.GetNode(2, 2).x, tolerance);
    ASSERT_NEAR(11.0, curvilinearGrid.GetNode(2, 3).x, tolerance);

    ASSERT_NEAR(constants::missing::doubleValue, curvilinearGrid.GetNode(3, 0).y, tolerance);
    ASSERT_NEAR(constants::missing::doubleValue, curvilinearGrid.GetNode(3, 1).y, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid.GetNode(3, 2).y, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid.GetNode(3, 3).y, tolerance);

    // vertical gridline 10
    ASSERT_NEAR(constants::missing::doubleValue, curvilinearGrid.GetNode(10, 0).x, tolerance);
    ASSERT_NEAR(constants::missing::doubleValue, curvilinearGrid.GetNode(10, 1).x, tolerance);
    ASSERT_NEAR(19.0, curvilinearGrid.GetNode(10, 2).x, tolerance);
    ASSERT_NEAR(19.0, curvilinearGrid.GetNode(10, 3).x, tolerance);

    ASSERT_NEAR(constants::missing::doubleValue, curvilinearGrid.GetNode(10, 0).y, tolerance);
    ASSERT_NEAR(constants::missing::doubleValue, curvilinearGrid.GetNode(10, 1).y, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid.GetNode(10, 2).y, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid.GetNode(10, 3).y, tolerance);

    // vertical gridline 11
    ASSERT_NEAR(constants::missing::doubleValue, curvilinearGrid.GetNode(11, 0).x, tolerance);
    ASSERT_NEAR(constants::missing::doubleValue, curvilinearGrid.GetNode(11, 1).x, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid.GetNode(11, 2).x, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid.GetNode(11, 3).x, tolerance);

    ASSERT_NEAR(constants::missing::doubleValue, curvilinearGrid.GetNode(11, 0).y, tolerance);
    ASSERT_NEAR(constants::missing::doubleValue, curvilinearGrid.GetNode(11, 1).y, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid.GetNode(11, 2).y, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid.GetNode(11, 3).y, tolerance);
}

TEST(CurvilinearGridRefinement, FullRefinementNonUniform)
{
    // Set-up
    lin_alg::Matrix<Point> gridPoints(4, 5);

    gridPoints << Point{0.0, 100.0}, Point{12.0, 100.0}, Point{25.0, 100.0}, Point{50.0, 100.0}, Point{80.0, 100.0},
        Point{0.0, 105.0}, Point{12.0, 105.0}, Point{25.0, 105.0}, Point{50.0, 105.0}, Point{80.0, 105.0},
        Point{0.0, 118.0}, Point{12.0, 118.0}, Point{25.0, 118.0}, Point{50.0, 118.0}, Point{80.0, 118.0},
        Point{0.0, 135.0}, Point{12.0, 135.0}, Point{25.0, 135.0}, Point{50.0, 135.0}, Point{80.0, 135.0};

    std::vector<double> expectedX{0.0, 3.1171875, 6.1875, 9.1640625, 12.0, 14.7109375, 17.5625, 20.8828125, 25.0,
                                  30.1484375, 36.1875, 42.8828125, 50.0, 57.3359375, 64.8125, 72.3828125, 80.0};

    std::vector<double> expectedY{100.0, 101.8, 105.0, 110.6, 118.0, 126.3, 135.0};

    CurvilinearGrid curvilinearGrid(gridPoints, Projection::cartesian);
    CurvilinearGridFullRefinement curvilinearGridRefinement;

    // Refine mesh four times in m-direction and twice in n-direction
    auto undoRefinement = curvilinearGridRefinement.Compute(curvilinearGrid, 4, 2);

    ASSERT_EQ(static_cast<meshkernel::UInt>(expectedX.size()), curvilinearGrid.NumM());
    ASSERT_EQ(static_cast<meshkernel::UInt>(expectedY.size()), curvilinearGrid.NumN());

    constexpr double tolerance = 1.0e-10;

    for (UInt i = 0; i < curvilinearGrid.NumN(); ++i)
    {
        for (UInt j = 0; j < curvilinearGrid.NumM(); ++j)
        {
            EXPECT_NEAR(curvilinearGrid.GetNode(i, j).x, expectedX[j], tolerance);
            EXPECT_NEAR(curvilinearGrid.GetNode(i, j).y, expectedY[i], tolerance);
        }
    }

    //--------------------------------
    // Now test undo
    undoRefinement->Restore();

    ASSERT_EQ(static_cast<meshkernel::UInt>(gridPoints.cols()), curvilinearGrid.NumM());
    ASSERT_EQ(static_cast<meshkernel::UInt>(gridPoints.rows()), curvilinearGrid.NumN());

    for (UInt i = 0; i < curvilinearGrid.NumN(); ++i)
    {
        for (UInt j = 0; j < curvilinearGrid.NumM(); ++j)
        {
            EXPECT_NEAR(curvilinearGrid.GetNode(i, j).x, gridPoints(i, j).x, tolerance);
            EXPECT_NEAR(curvilinearGrid.GetNode(i, j).y, gridPoints(i, j).y, tolerance);
        }
    }
}

TEST(CurvilinearGridRefinement, FullRefinementWithFactor1)
{
    // No refinement should occur in this test

    lin_alg::Matrix<Point> grid(8, 10);

    for (Eigen::Index i = 0; i < grid.rows(); ++i)
    {
        for (Eigen::Index j = 0; j < grid.cols(); ++j)
        {
            grid(i, j) = Point(static_cast<double>(j) * 10.0, static_cast<double>(i) * 10.0);
        }
    }

    CurvilinearGrid curvilinearGrid(grid, Projection::cartesian);
    CurvilinearGridFullRefinement curvilinearGridRefinement;

    const lin_alg::Matrix<Point> originalPoints = curvilinearGrid.GetNodes();

    auto undoRefinement = curvilinearGridRefinement.Compute(curvilinearGrid, 1, 1);

    // No refinement should take place, and so no undo necessary
    EXPECT_EQ(undoRefinement, nullptr);

    ASSERT_EQ(static_cast<meshkernel::UInt>(originalPoints.cols()), curvilinearGrid.NumM());
    ASSERT_EQ(static_cast<meshkernel::UInt>(originalPoints.rows()), curvilinearGrid.NumN());

    for (UInt i = 0; i < curvilinearGrid.NumN(); ++i)
    {
        for (UInt j = 0; j < curvilinearGrid.NumM(); ++j)
        {
            EXPECT_EQ(curvilinearGrid.GetNode(i, j).x, originalPoints(i, j).x);
            EXPECT_EQ(curvilinearGrid.GetNode(i, j).y, originalPoints(i, j).y);
        }
    }
}

TEST(CurvilinearGridRefinement, FullRefinementWithMissingRegion)
{
    constexpr double deltaX = 10.0;
    constexpr double deltaY = 12.0;

    lin_alg::Matrix<Point> gridNodes(6, 10);

    // Defined regions of invalid nodes
    auto region1 = [](const UInt i, const UInt j)
    { return (i == 0 || i == 1) && (3 <= j && j <= 6); };
    auto region2 = [](const UInt i, const UInt j)
    { return (i == 1 || i == 2) && (4 <= j && j <= 5); };

    for (int i = 0; i < gridNodes.rows(); ++i)
    {
        for (int j = 0; j < gridNodes.cols(); ++j)
        {
            if (region1(i, j) || region2(i, j))
            {
                gridNodes(i, j) = Point(constants::missing::doubleValue, constants::missing::doubleValue);
            }
            else
            {
                gridNodes(i, j) = Point(static_cast<double>(j) * deltaX, static_cast<double>(i) * deltaY);
            }
        }
    }

    CurvilinearGrid curvilinearGrid(gridNodes, Projection::cartesian);
    CurvilinearGridFullRefinement curvilinearGridRefinement;

    const UInt mRefinement = 2;
    const UInt nRefinement = 3;

    const UInt expectedSizeM = mRefinement * (static_cast<meshkernel::UInt>(gridNodes.cols()) - 1) + 1;
    const UInt expectedSizeN = nRefinement * (static_cast<meshkernel::UInt>(gridNodes.rows()) - 1) + 1;

    // Refine mesh
    [[maybe_unused]] auto undo = curvilinearGridRefinement.Compute(curvilinearGrid, mRefinement, nRefinement);

    ASSERT_EQ(expectedSizeM, curvilinearGrid.NumM());
    ASSERT_EQ(expectedSizeN, curvilinearGrid.NumN());

    constexpr double tolerance = 1.0e-10;

    for (UInt i = 0; i < curvilinearGrid.NumN(); ++i)
    {
        const double expectedY = static_cast<double>(i) * deltaY / static_cast<double>(nRefinement);

        for (UInt j = 0; j < curvilinearGrid.NumM(); ++j)
        {
            const double expectedX = static_cast<double>(j) * deltaX / static_cast<double>(mRefinement);

            if (curvilinearGrid.GetNode(i, j).IsValid())
            {
                EXPECT_NEAR(curvilinearGrid.GetNode(i, j).x, expectedX, tolerance);
                EXPECT_NEAR(curvilinearGrid.GetNode(i, j).y, expectedY, tolerance);
            }
        }
    }
}

TEST(CurvilinearGridRefinement, IncorrectFullRefinementParameters)
{
    lin_alg::Matrix<Point> grid(4, 5);

    grid << Point{0, 0}, Point{10, 0}, Point{20, 0}, Point{30, 0}, Point{40, 0},
        Point{0, 10}, Point{10, 10}, Point{20, 10}, Point{30, 10}, Point{40, 10},
        Point{0, 20}, Point{10, 20}, Point{20, 20}, Point{30, 20}, Point{40, 20},
        Point{0, 30}, Point{10, 30}, Point{20, 30}, Point{30, 30}, Point{40, 30};

    CurvilinearGrid curvilinearGrid(grid, Projection::cartesian);
    CurvilinearGridFullRefinement curvilinearGridRefinement;

    // Should raise expcetion due to incorrect m- and/or n-refinement factors
    EXPECT_THROW([[maybe_unused]] auto undo = curvilinearGridRefinement.Compute(curvilinearGrid, 0, 2), meshkernel::ConstraintError);
    EXPECT_THROW([[maybe_unused]] auto undo = curvilinearGridRefinement.Compute(curvilinearGrid, 3, 0), meshkernel::ConstraintError);
    EXPECT_THROW([[maybe_unused]] auto undo = curvilinearGridRefinement.Compute(curvilinearGrid, 0, 0), meshkernel::ConstraintError);
    EXPECT_THROW([[maybe_unused]] auto undo = curvilinearGridRefinement.Compute(curvilinearGrid, constants::missing::uintValue, 0), meshkernel::ConstraintError);
    EXPECT_THROW([[maybe_unused]] auto undo = curvilinearGridRefinement.Compute(curvilinearGrid, 0, constants::missing::uintValue), meshkernel::ConstraintError);
    EXPECT_THROW([[maybe_unused]] auto undo = curvilinearGridRefinement.Compute(curvilinearGrid, constants::missing::uintValue, constants::missing::uintValue), meshkernel::ConstraintError);
}
