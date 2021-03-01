#include <gtest/gtest.h>

#include <MeshKernel/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGridRefinement.hpp>
#include <MeshKernel/Entities.hpp>

TEST(CurvilinearGridRefinement, Compute_OnCurvilinearGrid_ShouldRefine)
{
    // Set-up
    std::vector<std::vector<meshkernel::Point>> grid{
        {{0, 0}, {0, 10}, {0, 20}, {0, 30}},
        {{10, 0}, {10, 10}, {10, 20}, {10, 30}},
        {{20, 0}, {20, 10}, {20, 20}, {20, 30}},
        {{30, 0}, {30, 10}, {30, 20}, {30, 30}}};

    const auto curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(grid, meshkernel::Projection::cartesian);
    meshkernel::CurvilinearGridRefinement curvilinearGridRefinement(curvilinearGrid, {10, 20}, {20, 20}, 10, 10);

    // Execute
    curvilinearGridRefinement.Compute();

    // Assert
    ASSERT_EQ(13, curvilinearGrid->m_numM);
    ASSERT_EQ(4, curvilinearGrid->m_numN);

    constexpr double tolerance = 1e-12;
    ASSERT_NEAR(0.0, curvilinearGrid->m_nodes[0][0].x, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid->m_nodes[1][0].x, tolerance);

    ASSERT_NEAR(11.0, curvilinearGrid->m_nodes[2][0].x, tolerance);
    ASSERT_NEAR(12.0, curvilinearGrid->m_nodes[3][0].x, tolerance);
    ASSERT_NEAR(13.0, curvilinearGrid->m_nodes[4][0].x, tolerance);
    ASSERT_NEAR(14.0, curvilinearGrid->m_nodes[5][0].x, tolerance);
    ASSERT_NEAR(15.0, curvilinearGrid->m_nodes[6][0].x, tolerance);
    ASSERT_NEAR(16.0, curvilinearGrid->m_nodes[7][0].x, tolerance);
    ASSERT_NEAR(17.0, curvilinearGrid->m_nodes[8][0].x, tolerance);
    ASSERT_NEAR(18.0, curvilinearGrid->m_nodes[9][0].x, tolerance);
    ASSERT_NEAR(19.0, curvilinearGrid->m_nodes[10][0].x, tolerance);

    ASSERT_NEAR(20.0, curvilinearGrid->m_nodes[11][0].x, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid->m_nodes[12][0].x, tolerance);

    ASSERT_NEAR(0.0, curvilinearGrid->m_nodes[0][0].y, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid->m_nodes[0][1].y, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid->m_nodes[0][2].y, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid->m_nodes[0][3].y, tolerance);
}

TEST(CurvilinearGridRefinement, Compute_OnCurvilinearGridWithMissingFaces_ShouldRefine)
{
    // Set-up
    std::vector<std::vector<meshkernel::Point>> grid{
        {{0, 0}, {0, 10}, {0, 20}, {0, 30}},
        {{10, 0}, {10, 10}, {10, 20}, {10, 30}},
        {{}, {}, {20, 20}, {20, 30}},
        {{}, {}, {30, 20}, {30, 30}},
        {{40, 0}, {40, 10}, {40, 20}, {40, 30}},
        {{50, 0}, {50, 10}, {50, 20}, {50, 30}}};

    const auto curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(grid, meshkernel::Projection::cartesian);
    meshkernel::CurvilinearGridRefinement curvilinearGridRefinement(curvilinearGrid, {10, 20}, {20, 20}, 10, 10);

    // Execute
    curvilinearGridRefinement.Compute();

    // Assert
    ASSERT_EQ(15, curvilinearGrid->m_numM);
    ASSERT_EQ(4, curvilinearGrid->m_numN);

    constexpr double tolerance = 1e-12;

    // vertical gridline 0
    ASSERT_NEAR(0.0, curvilinearGrid->m_nodes[0][0].x, tolerance);
    ASSERT_NEAR(0.0, curvilinearGrid->m_nodes[0][1].x, tolerance);
    ASSERT_NEAR(0.0, curvilinearGrid->m_nodes[0][2].x, tolerance);
    ASSERT_NEAR(0.0, curvilinearGrid->m_nodes[0][3].x, tolerance);

    ASSERT_NEAR(0.0, curvilinearGrid->m_nodes[0][0].y, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid->m_nodes[0][1].y, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid->m_nodes[0][2].y, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid->m_nodes[0][3].y, tolerance);

    // vertical gridline 2
    ASSERT_NEAR(meshkernel::doubleMissingValue, curvilinearGrid->m_nodes[2][0].x, tolerance);
    ASSERT_NEAR(meshkernel::doubleMissingValue, curvilinearGrid->m_nodes[2][1].x, tolerance);
    ASSERT_NEAR(11.0, curvilinearGrid->m_nodes[2][2].x, tolerance);
    ASSERT_NEAR(11.0, curvilinearGrid->m_nodes[2][3].x, tolerance);

    ASSERT_NEAR(meshkernel::doubleMissingValue, curvilinearGrid->m_nodes[3][0].y, tolerance);
    ASSERT_NEAR(meshkernel::doubleMissingValue, curvilinearGrid->m_nodes[3][1].y, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid->m_nodes[3][2].y, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid->m_nodes[3][3].y, tolerance);

    //vertical gridline 10
    ASSERT_NEAR(meshkernel::doubleMissingValue, curvilinearGrid->m_nodes[10][0].x, tolerance);
    ASSERT_NEAR(meshkernel::doubleMissingValue, curvilinearGrid->m_nodes[10][1].x, tolerance);
    ASSERT_NEAR(19.0, curvilinearGrid->m_nodes[10][2].x, tolerance);
    ASSERT_NEAR(19.0, curvilinearGrid->m_nodes[10][3].x, tolerance);

    ASSERT_NEAR(meshkernel::doubleMissingValue, curvilinearGrid->m_nodes[10][0].y, tolerance);
    ASSERT_NEAR(meshkernel::doubleMissingValue, curvilinearGrid->m_nodes[10][1].y, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid->m_nodes[10][2].y, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid->m_nodes[10][3].y, tolerance);

    //vertical gridline 11
    ASSERT_NEAR(meshkernel::doubleMissingValue, curvilinearGrid->m_nodes[11][0].x, tolerance);
    ASSERT_NEAR(meshkernel::doubleMissingValue, curvilinearGrid->m_nodes[11][1].x, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid->m_nodes[11][2].x, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid->m_nodes[11][3].x, tolerance);

    ASSERT_NEAR(meshkernel::doubleMissingValue, curvilinearGrid->m_nodes[11][0].y, tolerance);
    ASSERT_NEAR(meshkernel::doubleMissingValue, curvilinearGrid->m_nodes[11][1].y, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid->m_nodes[11][2].y, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid->m_nodes[11][3].y, tolerance);
}
