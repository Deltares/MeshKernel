#include <gtest/gtest.h>

#include <MeshKernel/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGridDeRefinement.hpp>
#include <MeshKernel/Entities.hpp>

TEST(CurvilinearGridDeRefinement, Compute_OnCurvilinearGrid_ShouldDeRefineVerticalGridLines)
{
    // Set-up
    std::vector<std::vector<meshkernel::Point>> grid{
        {{0, 0}, {0, 10}, {0, 20}, {0, 30}},
        {{10, 0}, {10, 10}, {10, 20}, {10, 30}},
        {{15, 0}, {15, 10}, {15, 20}, {15, 30}}, // this vertical grid line will be removed
        {{20, 0}, {20, 10}, {20, 20}, {20, 30}},
        {{30, 0}, {30, 10}, {30, 20}, {30, 30}}};

    const auto curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(std::move(grid), meshkernel::Projection::cartesian);
    meshkernel::CurvilinearGridDeRefinement curvilinearGridDeRefinement(curvilinearGrid, {10, 20}, {20, 20});

    // Execute
    const auto derefinedGrid = curvilinearGridDeRefinement.Compute();

    // Assert (the vertical line at x=15 is removed)
    ASSERT_EQ(4, derefinedGrid.m_numM);
    ASSERT_EQ(4, derefinedGrid.m_numN);
}

TEST(CurvilinearGridDeRefinement, Compute_OnCurvilinearGridWithMissingFaces_ShouldDeRefineVerticalGridLines)
{
    // Set-up
    std::vector<std::vector<meshkernel::Point>> grid{
        {{0, 0}, {0, 10}, {0, 20}, {0, 30}},
        {{10, 0}, {10, 10}, {10, 20}, {10, 30}},
        {{}, {}, {11, 20}, {11, 30}},
        {{}, {}, {12, 20}, {12, 30}},
        {{}, {}, {13, 20}, {13, 30}},
        {{}, {}, {20, 20}, {20, 30}},
        {{}, {}, {30, 20}, {30, 30}},
        {{40, 0}, {40, 10}, {40, 20}, {40, 30}},
        {{50, 0}, {50, 10}, {50, 20}, {50, 30}}};

    const auto curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(std::move(grid), meshkernel::Projection::cartesian);
    meshkernel::CurvilinearGridDeRefinement curvilinearGridDeRefinement(curvilinearGrid, {10, 20}, {20, 20});

    // Execute
    const auto derefinedGrid = curvilinearGridDeRefinement.Compute();

    // Assert
    ASSERT_EQ(6, derefinedGrid.m_numM);
    ASSERT_EQ(4, derefinedGrid.m_numN);
}

TEST(CurvilinearGridDeRefinement, Compute_OnCurvilinearGrid_ShouldDeRefineHorizontalGridLines)
{
    // Set-up
    std::vector<std::vector<meshkernel::Point>> grid{
        {{0, 0}, {0, 10}, {0, 11}, {0, 20}, {0, 30}},
        {{10, 0}, {10, 10}, {10, 11}, {10, 20}, {10, 30}},
        {{20, 0}, {20, 10}, {20, 11}, {20, 20}, {20, 30}},
        {{30, 0}, {30, 10}, {30, 11}, {30, 20}, {30, 30}}};

    const auto curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(std::move(grid), meshkernel::Projection::cartesian);
    meshkernel::CurvilinearGridDeRefinement curvilinearGridDeRefinement(curvilinearGrid, {10, 10}, {10, 20});

    // Execute
    const auto derefinedGrid = curvilinearGridDeRefinement.Compute();

    // Assert (the vertical line at x=15 is removed)
    ASSERT_EQ(4, derefinedGrid.m_numM);
    ASSERT_EQ(4, derefinedGrid.m_numN);
}