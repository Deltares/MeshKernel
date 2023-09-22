#include <gtest/gtest.h>

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridDeRefinement.hpp>
#include <MeshKernel/Entities.hpp>

using namespace meshkernel;

TEST(CurvilinearGridDeRefinement, Compute_OnCurvilinearGrid_ShouldDeRefineVerticalGridLines)
{
    // Set-up
    lin_alg::Matrix<Point> grid(5, 4);
    grid << Point{0, 0}, Point{0, 10}, Point{0, 20}, Point{0, 30},
        Point{10, 0}, Point{10, 10}, Point{10, 20}, Point{10, 30},
        Point{15, 0}, Point{15, 10}, Point{15, 20}, Point{15, 30}, // this vertical grid line will be removed
        Point{20, 0}, Point{20, 10}, Point{20, 20}, Point{20, 30},
        Point{30, 0}, Point{30, 10}, Point{30, 20}, Point{30, 30};

    const auto curvilinearGrid = std::make_shared<CurvilinearGrid>(grid, Projection::Type::Cartesian);
    CurvilinearGridDeRefinement curvilinearGridDeRefinement(curvilinearGrid);
    curvilinearGridDeRefinement.SetBlock({10, 20}, {20, 20});

    // Execute
    const auto derefinedGrid = curvilinearGridDeRefinement.Compute();

    // Assert (the vertical line at x=15 is removed)
    ASSERT_EQ(4, derefinedGrid.m_numM);
    ASSERT_EQ(4, derefinedGrid.m_numN);
}

TEST(CurvilinearGridDeRefinement, Compute_OnCurvilinearGridWithMissingFaces_ShouldDeRefineVerticalGridLines)
{
    // Set-up
    lin_alg::Matrix<Point> grid(9, 4);
    grid << Point{0, 0}, Point{0, 10}, Point{0, 20}, Point{0, 30},
        Point{10, 0}, Point{10, 10}, Point{10, 20}, Point{10, 30},
        Point{}, Point{}, Point{11, 20}, Point{11, 30},
        Point{}, Point{}, Point{12, 20}, Point{12, 30},
        Point{}, Point{}, Point{13, 20}, Point{13, 30},
        Point{}, Point{}, Point{20, 20}, Point{20, 30},
        Point{}, Point{}, Point{30, 20}, Point{30, 30},
        Point{40, 0}, Point{40, 10}, Point{40, 20}, Point{40, 30},
        Point{50, 0}, Point{50, 10}, Point{50, 20}, Point{50, 30};

    const auto curvilinearGrid = std::make_shared<CurvilinearGrid>(grid, Projection::Type::Cartesian);
    CurvilinearGridDeRefinement curvilinearGridDeRefinement(curvilinearGrid);
    curvilinearGridDeRefinement.SetBlock({10, 20}, {20, 20});

    // Execute

    const auto derefinedGrid = curvilinearGridDeRefinement.Compute();

    // Assert
    ASSERT_EQ(6, derefinedGrid.m_numM);
    ASSERT_EQ(4, derefinedGrid.m_numN);
}

TEST(CurvilinearGridDeRefinement, Compute_OnCurvilinearGrid_ShouldDeRefineHorizontalGridLines)
{
    // Set-up
    lin_alg::Matrix<Point> grid(4, 5);
    grid << Point{0, 0}, Point{0, 10}, Point{0, 11}, Point{0, 20}, Point{0, 30},
        Point{10, 0}, Point{10, 10}, Point{10, 11}, Point{10, 20}, Point{10, 30},
        Point{20, 0}, Point{20, 10}, Point{20, 11}, Point{20, 20}, Point{20, 30},
        Point{30, 0}, Point{30, 10}, Point{30, 11}, Point{30, 20}, Point{30, 30};

    auto curvilinearGrid = std::make_shared<CurvilinearGrid>(grid, Projection::Type::Cartesian);
    CurvilinearGridDeRefinement curvilinearGridDeRefinement(curvilinearGrid);
    curvilinearGridDeRefinement.SetBlock({10, 10}, {10, 20});

    // Execute
    const auto derefinedGrid = curvilinearGridDeRefinement.Compute();

    // Assert (the vertical line at x=15 is removed)
    ASSERT_EQ(4, derefinedGrid.m_numM);
    ASSERT_EQ(4, derefinedGrid.m_numN);
}