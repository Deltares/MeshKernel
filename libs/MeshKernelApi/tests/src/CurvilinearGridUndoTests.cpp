#include <algorithm>
#include <gtest/gtest.h>
#include <random>

#include "MeshKernel/Parameters.hpp"
#include "MeshKernelApi/Mesh2D.hpp"
#include "MeshKernelApi/MeshKernel.hpp"

#include "TestUtils/MakeMeshes.hpp"

TEST(CurvilinearGridUndoTests, DeleteNode)
{
    // Prepare
    int meshKernelId;
    auto errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::MakeGridParameters makeGridParameters;

    makeGridParameters.num_columns = 9;
    makeGridParameters.num_rows = 9;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = 1.0;
    makeGridParameters.block_size_y = 1.0;

    // Execute
    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> node_x(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> node_y(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = node_x.data();
    curvilinearGrid.node_y = node_y.data();
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Make copy of node values.
    std::vector<double> originalNodeX(node_x);
    std::vector<double> originalNodeY(node_y);

    const int deletedNodeIndex = 55;
    constexpr double NullValue = -999.0;
    constexpr double tolerance = 1.0e-10;

    errorCode = meshkernelapi::mkernel_curvilinear_delete_node(meshKernelId, 5.0, 5.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    for (int i = 0; i < curvilinearGrid.num_m * curvilinearGrid.num_n; ++i)
    {
        if (i == deletedNodeIndex)
        {
            EXPECT_NEAR(curvilinearGrid.node_x[i], NullValue, tolerance);
            EXPECT_NEAR(curvilinearGrid.node_y[i], NullValue, tolerance);
        }
        else
        {
            EXPECT_NEAR(curvilinearGrid.node_x[i], originalNodeX[i], tolerance);
            EXPECT_NEAR(curvilinearGrid.node_y[i], originalNodeY[i], tolerance);
        }
    }

    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    bool didUndoOfDeleteNode = false;
    int undoId = meshkernel::constants::missing::intValue;

    errorCode = meshkernelapi::mkernel_undo_state(didUndoOfDeleteNode, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndoOfDeleteNode);
    ASSERT_EQ(meshKernelId, undoId);

    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    for (int i = 0; i < curvilinearGrid.num_m * curvilinearGrid.num_n; ++i)
    {
        EXPECT_NEAR(curvilinearGrid.node_x[i], originalNodeX[i], tolerance);
        EXPECT_NEAR(curvilinearGrid.node_y[i], originalNodeY[i], tolerance);
    }
}

TEST(CurvilinearGridUndoTests, MoveNode)
{
    // Prepare
    int meshKernelId;
    auto errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::MakeGridParameters makeGridParameters;

    makeGridParameters.num_columns = 9;
    makeGridParameters.num_rows = 9;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = 1.0;
    makeGridParameters.block_size_y = 1.0;

    // Execute
    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> node_x(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> node_y(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = node_x.data();
    curvilinearGrid.node_y = node_y.data();
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Make copy of node values.
    std::vector<double> originalNodeX(node_x);
    std::vector<double> originalNodeY(node_y);

    const int movedNodeIndex = 55;
    constexpr double tolerance = 1.0e-10;

    errorCode = meshkernelapi::mkernel_curvilinear_move_node(meshKernelId, 5.0, 5.0, 5.5, 5.25);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    for (int i = 0; i < curvilinearGrid.num_m * curvilinearGrid.num_n; ++i)
    {
        if (i == movedNodeIndex)
        {
            EXPECT_NEAR(curvilinearGrid.node_x[i], 5.5, tolerance);
            EXPECT_NEAR(curvilinearGrid.node_y[i], 5.25, tolerance);
        }
        else
        {
            EXPECT_NEAR(curvilinearGrid.node_x[i], originalNodeX[i], tolerance);
            EXPECT_NEAR(curvilinearGrid.node_y[i], originalNodeY[i], tolerance);
        }
    }

    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    bool didUndoOfDeleteNode = false;
    int undoId = meshkernel::constants::missing::intValue;

    errorCode = meshkernelapi::mkernel_undo_state(didUndoOfDeleteNode, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndoOfDeleteNode);
    ASSERT_EQ(meshKernelId, undoId);

    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    for (int i = 0; i < curvilinearGrid.num_m * curvilinearGrid.num_n; ++i)
    {
        EXPECT_NEAR(curvilinearGrid.node_x[i], originalNodeX[i], tolerance);
        EXPECT_NEAR(curvilinearGrid.node_y[i], originalNodeY[i], tolerance);
    }
}

TEST(CurvilinearGridUndoTests, Smoothing)
{
    int meshKernelId;
    auto errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::MakeGridParameters makeGridParameters;

    constexpr double delta = 1.0;

    makeGridParameters.num_columns = 30;
    makeGridParameters.num_rows = 30;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = delta;
    makeGridParameters.block_size_y = delta;

    // Create a uniform random distribution in (0.01 .. 0.25) * delta
    // lower bound non-zero to ensure simple check works
    std::uniform_real_distribution<double> distribution(0.01 * delta, 0.25 * delta);
    std::default_random_engine engine;

    // Generate curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> node_x(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> node_y(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = node_x.data();
    curvilinearGrid.node_y = node_y.data();
    // Get the nodal values
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    int count = 0;

    // Add some perturbation to the mesh before smoothing.
    for (int i = 0; i < curvilinearGrid.num_n; ++i)
    {
        for (int j = 0; j < curvilinearGrid.num_m; ++j)
        {
            curvilinearGrid.node_x[count] += distribution(engine);
            curvilinearGrid.node_y[count] += distribution(engine);
            ++count;
        }
    }

    // Make copy of node values.
    std::vector<double> originalNodeX(node_x);
    std::vector<double> originalNodeY(node_y);

    // Copy the node values back
    errorCode = meshkernelapi::mkernel_curvilinear_set(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    const int lowerBoundIndex = 10;
    const int upperBoundIndex = 20;

    constexpr double lowerBoundValue = static_cast<double>(lowerBoundIndex);
    constexpr double upperBoundValue = static_cast<double>(upperBoundIndex);

    // apply smoothing
    errorCode = meshkernelapi::mkernel_curvilinear_smoothing(meshKernelId, 1 /* smoothingIterations */,
                                                             lowerBoundValue, lowerBoundValue, upperBoundValue, upperBoundValue);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Get the current state of the curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Expect the nodes of the smoothed region to not be the same as the original nodes
    for (int i = lowerBoundIndex; i < upperBoundIndex; ++i)
    {
        for (int j = lowerBoundIndex; j < upperBoundIndex; ++j)
        {
            int pos = curvilinearGrid.num_n * j + i;
            EXPECT_NE(curvilinearGrid.node_x[pos], originalNodeX[pos]) << "position: " << i << "  " << j;
            EXPECT_NE(curvilinearGrid.node_y[pos], originalNodeY[pos]) << "position: " << i << "  " << j;
        }
    }

    bool didUndoOfSmoothing = false;
    int undoId = meshkernel::constants::missing::intValue;

    // Undo the smoothing
    errorCode = meshkernelapi::mkernel_undo_state(didUndoOfSmoothing, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndoOfSmoothing);
    ASSERT_EQ(meshKernelId, undoId);

    // Get the current state of the curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    constexpr double tolerance = 1.0e-10;

    // After undo of smoothing, compare all the nodes to check that they are the same as the original nodes.
    for (int i = 0; i < curvilinearGrid.num_m * curvilinearGrid.num_n; ++i)
    {
        EXPECT_NEAR(curvilinearGrid.node_x[i], originalNodeX[i], tolerance);
        EXPECT_NEAR(curvilinearGrid.node_y[i], originalNodeY[i], tolerance);
    }
}

#if 0
TEST(CurvilinearGridUndoTests, SmoothingDirectional)
{

    // Prepare
    int meshKernelId;
    auto errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::MakeGridParameters makeGridParameters;

    constexpr double delta = 1.0;

    makeGridParameters.num_columns = 30;
    makeGridParameters.num_rows = 30;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = delta;
    makeGridParameters.block_size_y = delta;

    // Create a uniform distribution in 0 .. 1.
    std::uniform_real_distribution<double> distribution(0.0, 0.25 * delta);
    std::default_random_engine engine;

    // Execute
    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> node_x(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> node_y(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = node_x.data();
    curvilinearGrid.node_y = node_y.data();
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    int count = 0;

    // Add some perturbation to the mesh before smoothing.
    for (int i = 0; i < curvilinearGrid.num_n; ++i)
    {
        for (int j = 0; j < curvilinearGrid.num_m; ++j)
        {
            curvilinearGrid.node_x[count] += distribution(engine);
            curvilinearGrid.node_y[count] += distribution(engine);
            ++count;
        }
    }

    // Make copy of node values.
    std::vector<double> originalNodeX(node_x);
    std::vector<double> originalNodeY(node_y);

    errorCode = meshkernelapi::mkernel_curvilinear_set(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_smoothing_directional(meshKernelId, 1, 10.0, 10.0, 30.0, 30.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    bool didUndoOfSmoothing = false;
    int undoId = meshkernel::constants::missing::intValue;

    errorCode = meshkernelapi::mkernel_undo_state(didUndoOfSmoothing, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndoOfSmoothing);
    ASSERT_EQ(meshKernelId, undoId);

    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    constexpr double tolerance = 1.0e-10;

    // After undo of smoothing, compare all the nodes to check that they are the same as the original nodes.
    for (int i = 0; i < curvilinearGrid.num_m * curvilinearGrid.num_n; ++i)
    {
        EXPECT_NEAR(curvilinearGrid.node_x[i], originalNodeX[i], tolerance);
        EXPECT_NEAR(curvilinearGrid.node_y[i], originalNodeY[i], tolerance);
    }
}
#endif

TEST(CurvilinearGridUndoTests, DeleteInterior)
{
    int meshKernelId;
    auto errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::MakeGridParameters makeGridParameters;

    constexpr double delta = 1.0;

    makeGridParameters.num_columns = 30;
    makeGridParameters.num_rows = 30;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = delta;
    makeGridParameters.block_size_y = delta;

    // Create a uniform random distribution in (0.01 .. 0.25) * delta
    // lower bound non-zero to ensure simple check works
    std::uniform_real_distribution<double> distribution(0.01 * delta, 0.25 * delta);
    std::default_random_engine engine;

    // Generate curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> node_x(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> node_y(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = node_x.data();
    curvilinearGrid.node_y = node_y.data();
    // Get the nodal values
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Make copy of node values.
    std::vector<double> originalNodeX(node_x);
    std::vector<double> originalNodeY(node_y);

    // Copy the node values back
    errorCode = meshkernelapi::mkernel_curvilinear_set(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    const int lowerBoundIndex = 10;
    const int upperBoundIndex = 20;

    constexpr double lowerBoundValue = static_cast<double>(lowerBoundIndex);
    constexpr double upperBoundValue = static_cast<double>(upperBoundIndex);

    // delete interior block
    errorCode = meshkernelapi::mkernel_curvilinear_delete_interior(meshKernelId,
                                                                   lowerBoundValue, lowerBoundValue,
                                                                   upperBoundValue, upperBoundValue);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Get the current state of the curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    constexpr double NullValue = -999.0;

    // Expect the nodes of the smoothed region to not be the same as the original nodes
    for (int i = lowerBoundIndex + 1; i < upperBoundIndex; ++i)
    {
        for (int j = lowerBoundIndex + 1; j < upperBoundIndex; ++j)
        {
            int pos = curvilinearGrid.num_n * j + i;
            EXPECT_EQ(curvilinearGrid.node_x[pos], NullValue) << "position: " << i << "  " << j << "  " << pos;
            EXPECT_EQ(curvilinearGrid.node_y[pos], NullValue) << "position: " << i << "  " << j << "  " << pos;
        }
    }

    bool didUndoOfDeleteInterior = false;
    int undoId = meshkernel::constants::missing::intValue;

    // Undo the delete interior block
    errorCode = meshkernelapi::mkernel_undo_state(didUndoOfDeleteInterior, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndoOfDeleteInterior);
    ASSERT_EQ(meshKernelId, undoId);

    // Get the current state of the curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // After undo of smoothing, compare all the nodes to check that they are the same as the original nodes.
    for (int i = 0; i < curvilinearGrid.num_m * curvilinearGrid.num_n; ++i)
    {
        EXPECT_NE(curvilinearGrid.node_x[i], NullValue);
        EXPECT_NE(curvilinearGrid.node_y[i], NullValue);
    }
}

TEST(CurvilinearGridUndoTests, DeleteExterior)
{
    int meshKernelId;
    auto errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::MakeGridParameters makeGridParameters;

    constexpr double delta = 1.0;

    makeGridParameters.num_columns = 30;
    makeGridParameters.num_rows = 30;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = delta;
    makeGridParameters.block_size_y = delta;

    // Generate curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> node_x(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> node_y(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = node_x.data();
    curvilinearGrid.node_y = node_y.data();
    // Get the nodal values
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Make copy of node values.
    std::vector<double> originalNodeX(node_x);
    std::vector<double> originalNodeY(node_y);

    // Copy the node values back
    errorCode = meshkernelapi::mkernel_curvilinear_set(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    const int lowerBoundIndex = 10;
    const int upperBoundIndex = 20;

    constexpr double lowerBoundValue = static_cast<double>(lowerBoundIndex);
    constexpr double upperBoundValue = static_cast<double>(upperBoundIndex);

    // delete exterior block
    errorCode = meshkernelapi::mkernel_curvilinear_delete_exterior(meshKernelId,
                                                                   lowerBoundValue, lowerBoundValue,
                                                                   upperBoundValue, upperBoundValue);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Get the current state of the curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    constexpr double NullValue = -999.0;

    // Check the block has been invalidated
    auto checkInvalidBlock = [=](int lowerRow, int upperRow, int lowerCol, int upperCol)
    {
        for (int i = lowerRow; i < upperRow; ++i)
        {
            for (int j = lowerCol; j < upperCol; ++j)
            {
                int pos = curvilinearGrid.num_n * j + i;
                EXPECT_EQ(curvilinearGrid.node_x[pos], NullValue) << "position: " << i << "  " << j << "  " << pos;
                EXPECT_EQ(curvilinearGrid.node_y[pos], NullValue) << "position: " << i << "  " << j << "  " << pos;
            }
        }
    };

    // Check the exterior block has been invalidated
    checkInvalidBlock(0, lowerBoundIndex, 0, makeGridParameters.num_columns);
    checkInvalidBlock(upperBoundIndex + 1, makeGridParameters.num_rows, 0, makeGridParameters.num_columns);
    checkInvalidBlock(lowerBoundIndex, upperBoundIndex, 0, lowerBoundIndex);
    checkInvalidBlock(lowerBoundIndex, upperBoundIndex, upperBoundIndex + 1, makeGridParameters.num_columns);

    bool didUndoOfDeleteExterior = false;
    int undoId = meshkernel::constants::missing::intValue;

    // Undo the delete exterior
    errorCode = meshkernelapi::mkernel_undo_state(didUndoOfDeleteExterior, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndoOfDeleteExterior);
    ASSERT_EQ(meshKernelId, undoId);

    // Get the current state of the curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // After undo of smoothing, compare all the nodes to check that they are the same as the original nodes.
    for (int i = 0; i < curvilinearGrid.num_m * curvilinearGrid.num_n; ++i)
    {
        EXPECT_NE(curvilinearGrid.node_x[i], NullValue);
        EXPECT_NE(curvilinearGrid.node_y[i], NullValue);
    }
}

TEST(CurvilinearGridUndoTests, Refine)
{
    int meshKernelId;
    auto errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::MakeGridParameters makeGridParameters;

    constexpr double delta = 1.0;

    makeGridParameters.num_columns = 30;
    makeGridParameters.num_rows = 30;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = delta;
    makeGridParameters.block_size_y = delta;

    // Create a uniform random distribution in (0.01 .. 0.25) * delta
    // lower bound non-zero to ensure simple check works
    std::uniform_real_distribution<double> distribution(0.01 * delta, 0.25 * delta);
    std::default_random_engine engine;

    // Generate curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> node_x(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> node_y(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = node_x.data();
    curvilinearGrid.node_y = node_y.data();
    // Get the nodal values
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Make copy of node values.
    std::vector<double> originalNodeX(node_x);
    std::vector<double> originalNodeY(node_y);

    errorCode = meshkernelapi::mkernel_curvilinear_refine(meshKernelId,
                                                          10.0, 20.0,
                                                          20.0, 20.0,
                                                          2 /* refinement */);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid refinedGrid{};

    // Get the current state of the curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, refinedGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Same number of rows
    EXPECT_EQ(refinedGrid.num_n, curvilinearGrid.num_n);
    // There should be more columns in the mesh
    EXPECT_GT(refinedGrid.num_m, curvilinearGrid.num_m);

    bool didUndoOfRefinement = false;
    int undoId = meshkernel::constants::missing::intValue;

    // Undo the refinement
    errorCode = meshkernelapi::mkernel_undo_state(didUndoOfRefinement, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndoOfRefinement);
    ASSERT_EQ(meshKernelId, undoId);

    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, refinedGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    EXPECT_EQ(refinedGrid.num_n, curvilinearGrid.num_n);
    EXPECT_EQ(refinedGrid.num_m, curvilinearGrid.num_m);

    // Get the current state of the curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    constexpr double tolerance = 1.0e-10;

    // After undo of smoothing, compare all the nodes to check that they are the same as the original nodes.
    for (int i = 0; i < curvilinearGrid.num_m * curvilinearGrid.num_n; ++i)
    {
        EXPECT_NEAR(curvilinearGrid.node_x[i], originalNodeX[i], tolerance);
        EXPECT_NEAR(curvilinearGrid.node_y[i], originalNodeY[i], tolerance);
    }
}

TEST(CurvilinearGridUndoTests, Derefine)
{
    int meshKernelId;
    auto errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::MakeGridParameters makeGridParameters;

    constexpr double delta = 1.0;

    makeGridParameters.num_columns = 30;
    makeGridParameters.num_rows = 30;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = delta;
    makeGridParameters.block_size_y = delta;

    // Create a uniform random distribution in (0.01 .. 0.25) * delta
    // lower bound non-zero to ensure simple check works
    std::uniform_real_distribution<double> distribution(0.01 * delta, 0.25 * delta);
    std::default_random_engine engine;

    // Generate curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> node_x(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> node_y(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = node_x.data();
    curvilinearGrid.node_y = node_y.data();
    // Get the nodal values
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Make copy of node values.
    std::vector<double> originalNodeX(node_x);
    std::vector<double> originalNodeY(node_y);

    errorCode = meshkernelapi::mkernel_curvilinear_derefine(meshKernelId,
                                                            10.0, 20.0,
                                                            20.0, 20.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid refinedGrid{};

    // Get the current state of the curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, refinedGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Same number of rows
    EXPECT_EQ(refinedGrid.num_n, curvilinearGrid.num_n);
    // There should be fewer columns in the mesh
    EXPECT_LT(refinedGrid.num_m, curvilinearGrid.num_m);

    bool didUndoOfRefinement = false;
    int undoId = meshkernel::constants::missing::intValue;

    // Undo the refinement
    errorCode = meshkernelapi::mkernel_undo_state(didUndoOfRefinement, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndoOfRefinement);
    ASSERT_EQ(meshKernelId, undoId);

    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, refinedGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    EXPECT_EQ(refinedGrid.num_n, curvilinearGrid.num_n);
    EXPECT_EQ(refinedGrid.num_m, curvilinearGrid.num_m);

    // Get the current state of the curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    constexpr double tolerance = 1.0e-10;

    // After undo of smoothing, compare all the nodes to check that they are the same as the original nodes.
    for (int i = 0; i < curvilinearGrid.num_m * curvilinearGrid.num_n; ++i)
    {
        EXPECT_NEAR(curvilinearGrid.node_x[i], originalNodeX[i], tolerance);
        EXPECT_NEAR(curvilinearGrid.node_y[i], originalNodeY[i], tolerance);
    }
}

TEST(CurvilinearGridUndoTests, InsertFace)
{
    int meshKernelId;
    auto errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::MakeGridParameters makeGridParameters;

    constexpr double delta = 1.0;

    makeGridParameters.num_columns = 30;
    makeGridParameters.num_rows = 30;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = delta;
    makeGridParameters.block_size_y = delta;

    // Create a uniform random distribution in (0.01 .. 0.25) * delta
    // lower bound non-zero to ensure simple check works
    std::uniform_real_distribution<double> distribution(0.01 * delta, 0.25 * delta);
    std::default_random_engine engine;

    // Generate curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> node_x(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> node_y(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = node_x.data();
    curvilinearGrid.node_y = node_y.data();
    // Get the nodal values
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Make copy of node values.
    std::vector<double> originalNodeX(node_x);
    std::vector<double> originalNodeY(node_y);

    errorCode = meshkernelapi::mkernel_curvilinear_insert_face(meshKernelId, 0.5 * delta, 0.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid refinedGrid{};

    // Get the current state of the curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, refinedGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // 1 extra row for the inserted face
    EXPECT_EQ(refinedGrid.num_n, curvilinearGrid.num_n + 1);
    // Same number of columns
    EXPECT_EQ(refinedGrid.num_m, curvilinearGrid.num_m);

    bool didUndoOfRefinement = false;
    int undoId = meshkernel::constants::missing::intValue;

    // Undo the refinement
    errorCode = meshkernelapi::mkernel_undo_state(didUndoOfRefinement, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndoOfRefinement);
    ASSERT_EQ(meshKernelId, undoId);

    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, refinedGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    EXPECT_EQ(refinedGrid.num_n, curvilinearGrid.num_n);
    EXPECT_EQ(refinedGrid.num_m, curvilinearGrid.num_m);

    // Get the current state of the curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    constexpr double tolerance = 1.0e-10;

    // After undo of smoothing, compare all the nodes to check that they are the same as the original nodes.
    for (int i = 0; i < curvilinearGrid.num_m * curvilinearGrid.num_n; ++i)
    {
        EXPECT_NEAR(curvilinearGrid.node_x[i], originalNodeX[i], tolerance);
        EXPECT_NEAR(curvilinearGrid.node_y[i], originalNodeY[i], tolerance);
    }
}

TEST(CurvilinearGridUndoTests, LineShift)
{
    // Prepare
    int meshKernelId;
    auto errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::MakeGridParameters makeGridParameters;

    makeGridParameters.num_columns = 9;
    makeGridParameters.num_rows = 9;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = 1.0;
    makeGridParameters.block_size_y = 1.0;

    // Execute
    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> node_x(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> node_y(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = node_x.data();
    curvilinearGrid.node_y = node_y.data();
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Make copy of node values.
    std::vector<double> originalNodeX(node_x);
    std::vector<double> originalNodeY(node_y);

    constexpr double tolerance = 1.0e-10;

    // Start compound action
    errorCode = meshkernelapi::mkernel_curvilinear_initialize_line_shift(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_set_block_line_shift(meshKernelId, 2.0, 3.0, 10.0, 10.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_set_line_line_shift(meshKernelId, 5.0, 3.0, 5.0, 15.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Creates an undo action
    errorCode = meshkernelapi::mkernel_curvilinear_move_node_line_shift(meshKernelId, 5.0, 5.0, 6.5, 6.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Creates an undo action
    errorCode = meshkernelapi::mkernel_curvilinear_line_shift(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_finalize_line_shift(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    // End compound action

    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    const int movedNodeIndex = 55;

    // Is it necessary to check the rest of the updated grid, should be covered in other api tests.
    EXPECT_NEAR(curvilinearGrid.node_x[movedNodeIndex], 6.5, tolerance);
    EXPECT_NEAR(curvilinearGrid.node_y[movedNodeIndex], 6.0, tolerance);

    //--------------------------------
    // Need to undo two times for each of the undo actions created

    bool didUndoOfLineShift = false;
    bool didUndoOfMoveNode = false;
    int undoId = meshkernel::constants::missing::intValue;

    errorCode = meshkernelapi::mkernel_undo_state(didUndoOfLineShift, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndoOfLineShift);
    ASSERT_EQ(meshKernelId, undoId);

    errorCode = meshkernelapi::mkernel_undo_state(didUndoOfMoveNode, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndoOfMoveNode);
    ASSERT_EQ(meshKernelId, undoId);

    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    for (int i = 0; i < curvilinearGrid.num_m * curvilinearGrid.num_n; ++i)
    {
        EXPECT_NEAR(curvilinearGrid.node_x[i], originalNodeX[i], tolerance);
        EXPECT_NEAR(curvilinearGrid.node_y[i], originalNodeY[i], tolerance);
    }
}

TEST(CurvilinearGridUndoTests, LineMirror)
{
    // Prepare
    int meshKernelId;
    auto errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::MakeGridParameters makeGridParameters;

    constexpr double delta = 1.0;

    makeGridParameters.num_columns = 30;
    makeGridParameters.num_rows = 30;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = delta;
    makeGridParameters.block_size_y = delta;

    // Execute
    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> node_x(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> node_y(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = node_x.data();
    curvilinearGrid.node_y = node_y.data();
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Make copy of node values.
    std::vector<double> originalNodeX(node_x);
    std::vector<double> originalNodeY(node_y);

    constexpr double tolerance = 1.0e-10;

    // Add an extra column on the left of the domain
    errorCode = meshkernelapi::mkernel_curvilinear_line_mirror(meshKernelId, delta, 0.0, 0.0, 0.0, 30.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid extendedGrid{};

    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, extendedGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Same number of rows
    EXPECT_EQ(extendedGrid.num_n, curvilinearGrid.num_n);
    // One extra row
    EXPECT_EQ(extendedGrid.num_m, curvilinearGrid.num_m + 1);

    bool didUndoOfInsertLine = false;
    int undoId = meshkernel::constants::missing::intValue;

    errorCode = meshkernelapi::mkernel_undo_state(didUndoOfInsertLine, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndoOfInsertLine);
    ASSERT_EQ(meshKernelId, undoId);

    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    for (int i = 0; i < curvilinearGrid.num_m * curvilinearGrid.num_n; ++i)
    {
        EXPECT_NEAR(curvilinearGrid.node_x[i], originalNodeX[i], tolerance);
        EXPECT_NEAR(curvilinearGrid.node_y[i], originalNodeY[i], tolerance);
    }
}

TEST(CurvilinearGridUndoTests, LineAttractionRepulsion)
{
    // Prepare
    int meshKernelId;
    auto errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::MakeGridParameters makeGridParameters;

    constexpr double delta = 10.0;

    makeGridParameters.num_columns = 5;
    makeGridParameters.num_rows = 5;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = delta;
    makeGridParameters.block_size_y = delta;

    // Execute
    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> node_x(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> node_y(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = node_x.data();
    curvilinearGrid.node_y = node_y.data();
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Make copy of node values.
    std::vector<double> originalNodeX(node_x);
    std::vector<double> originalNodeY(node_y);

    constexpr double tolerance = 1.0e-10;

    // Add an extra column on the left of the domain
    errorCode = meshkernelapi::mkernel_curvilinear_line_attraction_repulsion(meshKernelId,
                                                                             0.5,
                                                                             30.0, 0.0, 30.0, 50.0,
                                                                             10.0, 20.0, 50.0, 20.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Check some nodes
    ASSERT_NEAR(17.5, curvilinearGrid.node_x[2], tolerance);
    ASSERT_NEAR(0.0, curvilinearGrid.node_y[2], tolerance);
    ASSERT_NEAR(42.5, curvilinearGrid.node_x[4], tolerance);
    ASSERT_NEAR(0.0, curvilinearGrid.node_y[4], tolerance);

    bool didUndoOfInsertLine = false;
    int undoId = meshkernel::constants::missing::intValue;

    errorCode = meshkernelapi::mkernel_undo_state(didUndoOfInsertLine, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndoOfInsertLine);
    ASSERT_EQ(meshKernelId, undoId);

    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    for (int i = 0; i < curvilinearGrid.num_m * curvilinearGrid.num_n; ++i)
    {
        EXPECT_NEAR(curvilinearGrid.node_x[i], originalNodeX[i], tolerance);
        EXPECT_NEAR(curvilinearGrid.node_y[i], originalNodeY[i], tolerance);
    }
}

TEST(CurvilinearGridUndoTests, OrthogonaliseEntireGrid)
{
    int meshKernelId;
    auto errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::MakeGridParameters makeGridParameters;

    constexpr double delta = 1.0;

    makeGridParameters.num_columns = 30;
    makeGridParameters.num_rows = 30;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = delta;
    makeGridParameters.block_size_y = delta;

    // Create a uniform random distribution in (0.01 .. 0.25) * delta
    // lower bound non-zero to ensure simple check works
    std::uniform_real_distribution<double> distribution(0.01 * delta, 0.25 * delta);
    std::default_random_engine engine;

    // Generate curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> node_x(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> node_y(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = node_x.data();
    curvilinearGrid.node_y = node_y.data();
    // Get the nodal values
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Add some perturbation to the mesh before smoothing.
    for (int i = 0; i < curvilinearGrid.num_n * curvilinearGrid.num_m; ++i)
    {
        curvilinearGrid.node_x[i] += distribution(engine);
        curvilinearGrid.node_y[i] += distribution(engine);
    }

    // Make copy of node values.
    std::vector<double> originalNodeX(node_x);
    std::vector<double> originalNodeY(node_y);

    // Copy the node values back
    errorCode = meshkernelapi::mkernel_curvilinear_set(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::OrthogonalizationParameters orthogonalizationParameters{};
    orthogonalizationParameters.outer_iterations = 1;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;

    //--------------------------------

    // set up orthogonalisation
    errorCode = meshkernelapi::mkernel_curvilinear_initialize_orthogonalize(meshKernelId, orthogonalizationParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_set_block_orthogonalize(meshKernelId, 0.0, 0.0, 30.0, 30.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // apply orthogonalisation
    errorCode = meshkernelapi::mkernel_curvilinear_orthogonalize(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // finalise orthogonalisation
    errorCode = meshkernelapi::mkernel_curvilinear_finalize_orthogonalize(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    //--------------------------------

    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Get the current state of the curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    const std::vector<int> cornerNodes{0, 30, 930, 960};

    // Expect the nodes of the orthogonalised mesh to be different from the original nodes
    for (int i = 0; i < curvilinearGrid.num_m * curvilinearGrid.num_n; ++i)
    {
        // Corner nodes are not moved
        if (std::ranges::find(cornerNodes, i) == cornerNodes.end())
        {
            EXPECT_NE(curvilinearGrid.node_x[i], originalNodeX[i]) << "position: " << i;
            EXPECT_NE(curvilinearGrid.node_y[i], originalNodeY[i]) << "position: " << i;
        }
    }

    bool didUndoOfSmoothing = false;
    int undoId = meshkernel::constants::missing::intValue;

    // Undo the orthogonalisation
    errorCode = meshkernelapi::mkernel_undo_state(didUndoOfSmoothing, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndoOfSmoothing);
    ASSERT_EQ(meshKernelId, undoId);

    // Get the current state of the curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // After undo of smoothing, compare all the nodes to check that they are the same as the original nodes.
    for (int i = 0; i < curvilinearGrid.num_m * curvilinearGrid.num_n; ++i)
    {
        EXPECT_EQ(curvilinearGrid.node_x[i], originalNodeX[i]);
        EXPECT_EQ(curvilinearGrid.node_y[i], originalNodeY[i]);
    }
}

TEST(CurvilinearGridUndoTests, RefineAndOrthogonalise)
{
    int meshKernelId;
    auto errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::MakeGridParameters makeGridParameters;

    constexpr double delta = 1.0;

    makeGridParameters.num_columns = 30;
    makeGridParameters.num_rows = 30;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = delta;
    makeGridParameters.block_size_y = delta;

    // Create a uniform random distribution in (0.01 .. 0.25) * delta
    // lower bound non-zero to ensure simple check works
    std::uniform_real_distribution<double> distribution(0.01 * delta, 0.25 * delta);
    std::default_random_engine engine;

    // Generate curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> node_x(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> node_y(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = node_x.data();
    curvilinearGrid.node_y = node_y.data();
    // Get the nodal values
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Add some perturbation to the mesh before smoothing.
    for (int i = 0; i < curvilinearGrid.num_n * curvilinearGrid.num_m; ++i)
    {
        curvilinearGrid.node_x[i] += distribution(engine);
        curvilinearGrid.node_y[i] += distribution(engine);
    }

    // Make copy of node values.
    std::vector<double> originalNodeX(node_x);
    std::vector<double> originalNodeY(node_y);

    // Copy the node values back
    errorCode = meshkernelapi::mkernel_curvilinear_set(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // End of curvilinear grid initialisation
    //--------------------------------

    // Refine grid
    errorCode = meshkernelapi::mkernel_curvilinear_refine(meshKernelId,
                                                          10.0, 20.0,
                                                          20.0, 20.0,
                                                          2 /* refinement */);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    //--------------------------------
    // Perform orthogonalisation on perturbed and refined grid

    meshkernel::OrthogonalizationParameters orthogonalizationParameters{};
    orthogonalizationParameters.outer_iterations = 1;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;

    // set up orthogonalisation
    errorCode = meshkernelapi::mkernel_curvilinear_initialize_orthogonalize(meshKernelId, orthogonalizationParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_set_block_orthogonalize(meshKernelId, 0.0, 0.0, 30.0, 30.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // apply orthogonalisation
    errorCode = meshkernelapi::mkernel_curvilinear_orthogonalize(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // finalise orthogonalisation
    errorCode = meshkernelapi::mkernel_curvilinear_finalize_orthogonalize(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    //--------------------------------
    // Undo orthogonalisation and refinement

    bool didUndo = false;
    int undoId = meshkernel::constants::missing::intValue;

    // Undo the orthogonalisation
    errorCode = meshkernelapi::mkernel_undo_state(didUndo, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndo);
    ASSERT_EQ(meshKernelId, undoId);

    didUndo = false;

    // Undo the orthogonalisation
    errorCode = meshkernelapi::mkernel_undo_state(didUndo, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndo);
    ASSERT_EQ(meshKernelId, undoId);

    //--------------------------------
    // The mesh should now be in the original (perturbed) state

    // Get the current state of the curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Get the current state of the curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // After undo of smoothing, compare all the nodes to check that they are the same as the original nodes.
    for (int i = 0; i < curvilinearGrid.num_m * curvilinearGrid.num_n; ++i)
    {
        EXPECT_EQ(curvilinearGrid.node_x[i], originalNodeX[i]) << "position: " << i;
        EXPECT_EQ(curvilinearGrid.node_y[i], originalNodeY[i]) << "position: " << i;
    }
}

TEST(CurvilinearGridUndoTests, RefineUndoThenOrthogonalise)
{
    int meshKernelId;
    auto errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::MakeGridParameters makeGridParameters;

    constexpr double delta = 1.0;

    makeGridParameters.num_columns = 30;
    makeGridParameters.num_rows = 30;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = delta;
    makeGridParameters.block_size_y = delta;

    // Create a uniform random distribution in (0.01 .. 0.25) * delta
    // lower bound non-zero to ensure simple check works
    std::uniform_real_distribution<double> distribution(0.01 * delta, 0.25 * delta);
    std::default_random_engine engine;

    // Generate curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> node_x(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> node_y(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = node_x.data();
    curvilinearGrid.node_y = node_y.data();
    // Get the nodal values
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Add some perturbation to the mesh before smoothing.
    for (int i = 0; i < curvilinearGrid.num_n * curvilinearGrid.num_m; ++i)
    {
        curvilinearGrid.node_x[i] += distribution(engine);
        curvilinearGrid.node_y[i] += distribution(engine);
    }

    // Make copy of node values.
    std::vector<double> originalNodeX(node_x);
    std::vector<double> originalNodeY(node_y);

    // Copy the node values back
    errorCode = meshkernelapi::mkernel_curvilinear_set(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // End of curvilinear grid initialisation
    //--------------------------------

    // Refine grid
    errorCode = meshkernelapi::mkernel_curvilinear_refine(meshKernelId,
                                                          10.0, 20.0,
                                                          20.0, 20.0,
                                                          2 /* refinement */);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    //--------------------------------
    // Undo the refinement
    bool didUndo = false;
    int undoId = meshkernel::constants::missing::intValue;

    errorCode = meshkernelapi::mkernel_undo_state(didUndo, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndo);
    ASSERT_EQ(meshKernelId, undoId);

    //--------------------------------
    // Perform orthogonalisation on perturbed grid

    meshkernel::OrthogonalizationParameters orthogonalizationParameters{};
    orthogonalizationParameters.outer_iterations = 1;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;

    // set up orthogonalisation
    errorCode = meshkernelapi::mkernel_curvilinear_initialize_orthogonalize(meshKernelId, orthogonalizationParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_set_block_orthogonalize(meshKernelId, 0.0, 0.0, 30.0, 30.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // apply orthogonalisation
    errorCode = meshkernelapi::mkernel_curvilinear_orthogonalize(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // finalise orthogonalisation
    errorCode = meshkernelapi::mkernel_curvilinear_finalize_orthogonalize(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    //--------------------------------
    // Undo orthognalisation

    didUndo = false;

    errorCode = meshkernelapi::mkernel_undo_state(didUndo, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndo);
    ASSERT_EQ(meshKernelId, undoId);

    //--------------------------------
    // The mesh should now be in the original (perturbed) state

    // Get the current state of the curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Get the current state of the curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // After undo of smoothing, compare all the nodes to check that they are the same as the original nodes.
    for (int i = 0; i < curvilinearGrid.num_m * curvilinearGrid.num_n; ++i)
    {
        EXPECT_EQ(curvilinearGrid.node_x[i], originalNodeX[i]) << "position: " << i;
        EXPECT_EQ(curvilinearGrid.node_y[i], originalNodeY[i]) << "position: " << i;
    }
}

TEST(CurvilinearGridUndoTests, InsertFaceUndoThenMirrorLine)
{
    int meshKernelId;
    auto errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::MakeGridParameters makeGridParameters;

    constexpr double delta = 1.0;

    makeGridParameters.num_columns = 30;
    makeGridParameters.num_rows = 30;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = delta;
    makeGridParameters.block_size_y = delta;

    // Create a uniform random distribution in (0.01 .. 0.25) * delta
    // lower bound non-zero to ensure simple check works
    std::uniform_real_distribution<double> distribution(0.01 * delta, 0.25 * delta);
    std::default_random_engine engine;

    // Generate curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> node_x(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> node_y(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = node_x.data();
    curvilinearGrid.node_y = node_y.data();
    // Get the nodal values
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Make copy of node values.
    std::vector<double> originalNodeX(node_x);
    std::vector<double> originalNodeY(node_y);

    // End of curvilinear grid initialisation
    //--------------------------------

    // Insert face
    errorCode = meshkernelapi::mkernel_curvilinear_insert_face(meshKernelId, 0.5 * delta, 0.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid updatedGrid{};

    // Get the current state of the curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, updatedGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // 1 extra row for the inserted face
    EXPECT_EQ(updatedGrid.num_n, curvilinearGrid.num_n + 1);
    // Same number of columns
    EXPECT_EQ(updatedGrid.num_m, curvilinearGrid.num_m);

    // Undo insert face
    bool didUndo = false;
    int undoId = meshkernel::constants::missing::intValue;

    // Undo the refinement
    errorCode = meshkernelapi::mkernel_undo_state(didUndo, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndo);
    ASSERT_EQ(meshKernelId, undoId);

    //-------------------------------

    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, updatedGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Add an extra column on the left of the domain
    errorCode = meshkernelapi::mkernel_curvilinear_line_mirror(meshKernelId, delta, 0.0, 0.0, 0.0, 30.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, updatedGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, updatedGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Same number of rows as original mesh
    EXPECT_EQ(updatedGrid.num_n, curvilinearGrid.num_n);
    // One extra row than the original mesh
    EXPECT_EQ(updatedGrid.num_m, curvilinearGrid.num_m + 1);

    // Undo the mirror line
    errorCode = meshkernelapi::mkernel_undo_state(didUndo, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndo);
    ASSERT_EQ(meshKernelId, undoId);

    //--------------------------------
    // Grid should be in original state

    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, updatedGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    EXPECT_EQ(updatedGrid.num_n, curvilinearGrid.num_n);
    EXPECT_EQ(updatedGrid.num_m, curvilinearGrid.num_m);

    // Get the current state of the curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    constexpr double tolerance = 1.0e-10;

    // After undo of smoothing, compare all the nodes to check that they are the same as the original nodes.
    for (int i = 0; i < curvilinearGrid.num_m * curvilinearGrid.num_n; ++i)
    {
        EXPECT_NEAR(curvilinearGrid.node_x[i], originalNodeX[i], tolerance);
        EXPECT_NEAR(curvilinearGrid.node_y[i], originalNodeY[i], tolerance);
    }
}

TEST(CurvilinearGridUndoTests, MultiStepUndoTest)
{
    // A multi-step undo test
    // 1. Refine grid in x-direciton
    // 2. undo refinement
    // 3. delete interior block
    // 4. undo deletion
    // 5. Refine grid in y-direction
    // 6. undo refinement

    int meshKernelId;
    auto errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::MakeGridParameters makeGridParameters;

    constexpr double delta = 1.0;

    makeGridParameters.num_columns = 30;
    makeGridParameters.num_rows = 30;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = delta;
    makeGridParameters.block_size_y = delta;

    // Create a uniform random distribution in (0.01 .. 0.25) * delta
    // lower bound non-zero to ensure simple check works
    std::uniform_real_distribution<double> distribution(0.01 * delta, 0.25 * delta);
    std::default_random_engine engine;

    // Generate curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> node_x(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> node_y(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = node_x.data();
    curvilinearGrid.node_y = node_y.data();
    // Get the nodal values
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Make copy of node values.
    std::vector<double> originalNodeX(node_x);
    std::vector<double> originalNodeY(node_y);

    meshkernelapi::CurvilinearGrid updatedGrid{};

    // End of curvilinear grid initialisation

    //--------------------------------
    // Refine grid
    errorCode = meshkernelapi::mkernel_curvilinear_refine(meshKernelId,
                                                          10.0, 20.0,
                                                          20.0, 20.0,
                                                          2 /* refinement */);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, updatedGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    EXPECT_EQ(updatedGrid.num_n, curvilinearGrid.num_n);
    EXPECT_EQ(updatedGrid.num_m, curvilinearGrid.num_m + 10);

    // Undo the refinement
    bool didUndo = false;
    int undoId = meshkernel::constants::missing::intValue;

    errorCode = meshkernelapi::mkernel_undo_state(didUndo, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndo);
    ASSERT_EQ(meshKernelId, undoId);

    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, updatedGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    EXPECT_EQ(updatedGrid.num_n, curvilinearGrid.num_n);
    EXPECT_EQ(updatedGrid.num_m, curvilinearGrid.num_m);

    //--------------------------------
    // Delete interior block

    const int lowerBoundIndex = 10;
    const int upperBoundIndex = 20;

    constexpr double lowerBoundValue = static_cast<double>(lowerBoundIndex);
    constexpr double upperBoundValue = static_cast<double>(upperBoundIndex);

    // delete interior block
    errorCode = meshkernelapi::mkernel_curvilinear_delete_interior(meshKernelId,
                                                                   lowerBoundValue, lowerBoundValue,
                                                                   upperBoundValue, upperBoundValue);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    didUndo = false;
    errorCode = meshkernelapi::mkernel_undo_state(didUndo, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndo);
    ASSERT_EQ(meshKernelId, undoId);

    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, updatedGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    EXPECT_EQ(updatedGrid.num_n, curvilinearGrid.num_n);
    EXPECT_EQ(updatedGrid.num_m, curvilinearGrid.num_m);

    //-------------------------------

    // Refine grid
    errorCode = meshkernelapi::mkernel_curvilinear_refine(meshKernelId,
                                                          0.0, 10.0,
                                                          0.0, 20.0,
                                                          2 /* refinement */);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, updatedGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    EXPECT_EQ(updatedGrid.num_n, curvilinearGrid.num_n + 10);
    EXPECT_EQ(updatedGrid.num_m, curvilinearGrid.num_m);

    // Undo the refinement
    didUndo = false;

    errorCode = meshkernelapi::mkernel_undo_state(didUndo, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndo);
    ASSERT_EQ(meshKernelId, undoId);

    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, updatedGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    EXPECT_EQ(updatedGrid.num_n, curvilinearGrid.num_n);
    EXPECT_EQ(updatedGrid.num_m, curvilinearGrid.num_m);

    //--------------------------------
    // Grid should be in original state

    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, updatedGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    EXPECT_EQ(updatedGrid.num_n, curvilinearGrid.num_n);
    EXPECT_EQ(updatedGrid.num_m, curvilinearGrid.num_m);

    // Get the current state of the curvilinear grid
    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    constexpr double tolerance = 1.0e-10;

    // After undo of smoothing, compare all the nodes to check that they are the same as the original nodes.
    for (int i = 0; i < curvilinearGrid.num_m * curvilinearGrid.num_n; ++i)
    {
        EXPECT_NEAR(curvilinearGrid.node_x[i], originalNodeX[i], tolerance);
        EXPECT_NEAR(curvilinearGrid.node_y[i], originalNodeY[i], tolerance);
    }
}
