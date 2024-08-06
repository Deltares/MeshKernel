#include <algorithm>
#include <gtest/gtest.h>
#include <random>

#include "MeshKernel/Parameters.hpp"
#include "MeshKernelApi/Mesh2D.hpp"
#include "MeshKernelApi/MeshKernel.hpp"

#include "TestUtils/MakeMeshes.hpp"

// namespace aliases
namespace mk = meshkernel;
namespace mkapi = meshkernelapi;

int GenerateCurvilinearMesh(const int meshKernelId, const int nodes, const double delta, const double origin)
{
    meshkernel::MakeGridParameters makeGridParameters;

    // num_columns and num_rows indicate number of elements in each direction, so value = nodes - 1
    makeGridParameters.num_columns = nodes - 1;
    makeGridParameters.num_rows = nodes - 1;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = origin;
    makeGridParameters.origin_y = origin;
    makeGridParameters.block_size_x = delta;
    makeGridParameters.block_size_y = delta;

    // Generate curvilinear grid.
    int errorCode = mkapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    return errorCode;
}

bool CompareCurvilinearGrids(const int meshKernelId1, const int meshKernelId2)
{
    mkapi::CurvilinearGrid curvilinearGrid1{};
    int errorCode = mkapi::mkernel_curvilinear_get_dimensions(meshKernelId1, curvilinearGrid1);

    if (mk::ExitCode::Success != errorCode)
    {
        return false;
    }

    std::vector<double> grid1NodesX(curvilinearGrid1.num_m * curvilinearGrid1.num_n);
    std::vector<double> grid1NodesY(curvilinearGrid1.num_m * curvilinearGrid1.num_n);

    curvilinearGrid1.node_x = grid1NodesX.data();
    curvilinearGrid1.node_y = grid1NodesY.data();

    errorCode = mkapi::mkernel_curvilinear_get_data(meshKernelId1, curvilinearGrid1);

    if (mk::ExitCode::Success != errorCode)
    {
        return false;
    }

    mkapi::CurvilinearGrid curvilinearGrid2{};

    errorCode = mkapi::mkernel_curvilinear_get_dimensions(meshKernelId2, curvilinearGrid2);

    if (mk::ExitCode::Success != errorCode)
    {
        return false;
    }

    std::vector<double> grid2NodesX(curvilinearGrid2.num_m * curvilinearGrid2.num_n);
    std::vector<double> grid2NodesY(curvilinearGrid2.num_m * curvilinearGrid2.num_n);

    curvilinearGrid2.node_x = grid2NodesX.data();
    curvilinearGrid2.node_y = grid2NodesY.data();

    if (curvilinearGrid1.num_n != curvilinearGrid2.num_n || curvilinearGrid1.num_m != curvilinearGrid2.num_m)
    {
        return false;
    }

    errorCode = mkapi::mkernel_curvilinear_get_data(meshKernelId2, curvilinearGrid2);

    size_t count = 0;

    for (int i = 0; i < curvilinearGrid1.num_n; ++i)
    {
        for (int j = 0; j < curvilinearGrid1.num_n; ++j)
        {
            if (curvilinearGrid1.node_x[count] != curvilinearGrid2.node_x[count] || curvilinearGrid1.node_y[count] != curvilinearGrid2.node_y[count])
            {
                return false;
            }

            ++count;
        }
    }

    return true;
}

bool CompareUnstructuredGrids(const int meshKernelId1, const int meshKernelId2)
{
    mkapi::Mesh2D grid1{};
    int errorCode = mkapi::mkernel_mesh2d_get_dimensions(meshKernelId1, grid1);

    if (mk::ExitCode::Success != errorCode)
    {
        return false;
    }

    std::vector<double> grid1NodesX(grid1.num_nodes);
    std::vector<double> grid1NodesY(grid1.num_nodes);
    std::vector<int> grid1Edges(2 * grid1.num_edges);

    grid1.node_x = grid1NodesX.data();
    grid1.node_y = grid1NodesY.data();
    grid1.edge_nodes = grid1Edges.data();

    errorCode = mkapi::mkernel_mesh2d_get_node_edge_data(meshKernelId1, grid1);

    if (mk::ExitCode::Success != errorCode)
    {
        return false;
    }

    mkapi::Mesh2D grid2{};

    errorCode = mkapi::mkernel_mesh2d_get_dimensions(meshKernelId2, grid2);

    if (mk::ExitCode::Success != errorCode)
    {
        return false;
    }

    std::vector<double> grid2NodesX(grid2.num_nodes);
    std::vector<double> grid2NodesY(grid2.num_nodes);
    std::vector<int> grid2Edges(2 * grid2.num_edges);

    grid2.node_x = grid2NodesX.data();
    grid2.node_y = grid2NodesY.data();
    grid2.edge_nodes = grid2Edges.data();

    if (grid1.num_nodes != grid2.num_nodes)
    {
        return false;
    }

    errorCode = mkapi::mkernel_mesh2d_get_node_edge_data(meshKernelId2, grid2);

    size_t count = 0;

    for (int i = 0; i < grid1.num_nodes; ++i)
    {
        if (grid1.node_x[count] != grid2.node_x[count] || grid1.node_y[count] != grid2.node_y[count])
        {
            return false;
        }

        ++count;
    }

    count = 0;

    for (int i = 0; i < 2 * grid1.num_nodes; ++i)
    {
        if (grid1.edge_nodes[count] != grid2.edge_nodes[count])
        {
            return false;
        }

        ++count;
    }

    return true;
}

TEST(UndoTests, BasicAllocationDeallocationTest)
{

    // Two mesh kernel ids will be created
    // checks made on their validity at different stages of the test
    // one of the id's will be deallocated
    // further checks on validity

    int meshKernelId1 = mkapi::mkernel_get_null_identifier();
    int meshKernelId2 = mkapi::mkernel_get_null_identifier();
    int errorCode;
    // Initialised with the opposite of the expected value
    bool isValid = true;

    // Clear the meshkernel state and undo stack before starting the test.
    errorCode = mkapi::mkernel_clear_state();
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    // meshKernelId1 should be not valid
    errorCode = mkapi::mkernel_is_valid_state(meshKernelId1, isValid);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_FALSE(isValid);

    errorCode = mkapi::mkernel_allocate_state(0, meshKernelId1);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    // Initialised with the opposite of the expected value
    // meshKernelId1 should now be valid
    isValid = false;
    errorCode = mkapi::mkernel_is_valid_state(meshKernelId1, isValid);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_TRUE(isValid);

    // Initialised with the opposite of the expected value
    // meshKernelId2 should still not be valid
    isValid = true;
    errorCode = mkapi::mkernel_is_valid_state(meshKernelId2, isValid);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_FALSE(isValid);

    errorCode = mkapi::mkernel_allocate_state(0, meshKernelId2);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    // deallocate meshKernelId1
    errorCode = mkapi::mkernel_deallocate_state(meshKernelId1);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    // Initialised with the opposite of the expected value
    // meshKernelId1 should not be valid
    isValid = true;
    errorCode = mkapi::mkernel_is_valid_state(meshKernelId1, isValid);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_FALSE(isValid);

    // Initialised with the opposite of the expected value
    // meshKernelId2 should still be valid
    isValid = false;
    errorCode = mkapi::mkernel_is_valid_state(meshKernelId2, isValid);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_TRUE(isValid);

    bool didUndo = false;
    int undoId = mkapi::mkernel_get_null_identifier();

    // Undo the deallocation
    errorCode = mkapi::mkernel_undo_state(didUndo, undoId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndo);
    EXPECT_EQ(undoId, meshKernelId1);

    // Initialised with the opposite of the expected value
    isValid = false;
    errorCode = mkapi::mkernel_is_valid_state(meshKernelId1, isValid);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_TRUE(isValid);

    // meshKernelId2 should still be valid
    isValid = false;
    errorCode = mkapi::mkernel_is_valid_state(meshKernelId2, isValid);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_TRUE(isValid);

    // Check undo, there should be nothing to undo
    errorCode = mkapi::mkernel_undo_state(didUndo, undoId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_FALSE(didUndo);
    EXPECT_EQ(undoId, mkapi::mkernel_get_null_identifier());

    // Finalise the test, remove all mesh kernel state objects and undo actions
    errorCode = mkapi::mkernel_expunge_state(meshKernelId2);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    errorCode = mkapi::mkernel_expunge_state(meshKernelId1);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
}

TEST(UndoTests, CurvilinearGridManipulationTests)
{
    const int clg1Size = 4;
    const int clg2Size = 8;

    int meshKernelId1 = mkapi::mkernel_get_null_identifier();
    int meshKernelId2 = mkapi::mkernel_get_null_identifier();
    int errorCode;
    // Initialised with the opposite of the expected value
    bool isValid = true;

    // Clear the meshkernel state and undo stack before starting the test.
    errorCode = mkapi::mkernel_clear_state();
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    // meshKernelId1 should be not valid
    errorCode = mkapi::mkernel_is_valid_state(meshKernelId1, isValid);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_FALSE(isValid);

    errorCode = mkapi::mkernel_allocate_state(0, meshKernelId1);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    errorCode = mkapi::mkernel_allocate_state(0, meshKernelId2);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    // Generate the curvilinear grids.
    errorCode = GenerateCurvilinearMesh(meshKernelId2, clg2Size, 1.0, 0.0);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    errorCode = GenerateCurvilinearMesh(meshKernelId1, clg1Size, 1.0, 3.0);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    // Check the sizes of the curvilinear grids.
    mkapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkapi::mkernel_curvilinear_get_dimensions(meshKernelId2, curvilinearGrid);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(curvilinearGrid.num_n, clg2Size);
    EXPECT_EQ(curvilinearGrid.num_m, clg2Size);

    errorCode = mkapi::mkernel_curvilinear_get_dimensions(meshKernelId1, curvilinearGrid);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(curvilinearGrid.num_n, clg1Size);
    EXPECT_EQ(curvilinearGrid.num_m, clg1Size);

    bool didUndo = false;
    int undoId = mkapi::mkernel_get_null_identifier();

    // Undo the curvilinear grid generation for the last grid generated, meshKernelId1
    errorCode = mkapi::mkernel_undo_state(didUndo, undoId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndo);
    // meshKernelId1 was the last mesh generated.
    EXPECT_EQ(undoId, meshKernelId1);

    // Getting the size of this grid should be an error and the dimensions the null value.
    errorCode = mkapi::mkernel_curvilinear_get_dimensions(meshKernelId1, curvilinearGrid);
    ASSERT_EQ(mk::ExitCode::MeshKernelErrorCode, errorCode);
    EXPECT_EQ(curvilinearGrid.num_n, mkapi::mkernel_get_null_identifier());
    EXPECT_EQ(curvilinearGrid.num_m, mkapi::mkernel_get_null_identifier());

    // Undo the undo of the clg generation
    // The clg for meshKernelId1 should be restored
    errorCode = mkapi::mkernel_redo_state(didUndo, undoId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(undoId, meshKernelId1);

    // With the correct size
    errorCode = mkapi::mkernel_curvilinear_get_dimensions(meshKernelId1, curvilinearGrid);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(curvilinearGrid.num_n, clg1Size);
    EXPECT_EQ(curvilinearGrid.num_m, clg1Size);

    // Deallocate clg 2
    errorCode = mkapi::mkernel_deallocate_state(meshKernelId2);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    errorCode = mkapi::mkernel_curvilinear_get_dimensions(meshKernelId2, curvilinearGrid);
    ASSERT_EQ(mk::ExitCode::MeshKernelErrorCode, errorCode);
    EXPECT_EQ(curvilinearGrid.num_n, mkapi::mkernel_get_null_identifier());
    EXPECT_EQ(curvilinearGrid.num_m, mkapi::mkernel_get_null_identifier());

    // Undo the deallocation or meshKernelId2
    errorCode = mkapi::mkernel_undo_state(didUndo, undoId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndo);
    EXPECT_EQ(undoId, meshKernelId2);

    // The grid should now have the correct size
    errorCode = mkapi::mkernel_curvilinear_get_dimensions(meshKernelId2, curvilinearGrid);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(curvilinearGrid.num_n, clg2Size);
    EXPECT_EQ(curvilinearGrid.num_m, clg2Size);

    // Convert the curvilinear grid id2 to an unstructured grid.
    errorCode = mkapi::mkernel_curvilinear_convert_to_mesh2d(meshKernelId2);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    mkapi::Mesh2D mesh2d{};
    errorCode = mkapi::mkernel_mesh2d_get_dimensions(meshKernelId2, mesh2d);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(mesh2d.num_nodes, clg2Size * clg2Size);
    EXPECT_EQ(mesh2d.num_edges, 2 * (clg2Size - 1) * clg2Size);

    // Now that the grid has been converted to an unstructured grid the dimensions of the curvilienar grid are not
    errorCode = mkapi::mkernel_curvilinear_get_dimensions(meshKernelId2, curvilinearGrid);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(curvilinearGrid.num_n, mkapi::mkernel_get_null_identifier());
    EXPECT_EQ(curvilinearGrid.num_m, mkapi::mkernel_get_null_identifier());

    // Undo the conversion of meshKernelId2 from curvilinear to unstructured.
    errorCode = mkapi::mkernel_undo_state(didUndo, undoId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndo);
    EXPECT_EQ(undoId, meshKernelId2);

    // The unstructured grid is no longer available
    errorCode = mkapi::mkernel_mesh2d_get_dimensions(meshKernelId2, mesh2d);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(mesh2d.num_nodes, mkapi::mkernel_get_null_identifier());
    EXPECT_EQ(mesh2d.num_edges, mkapi::mkernel_get_null_identifier());

    // Finalise the test, remove all mesh kernel state objects and undo actions
    errorCode = mkapi::mkernel_expunge_state(meshKernelId2);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    errorCode = mkapi::mkernel_expunge_state(meshKernelId1);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    errorCode = mkapi::mkernel_is_valid_state(meshKernelId1, isValid);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_FALSE(isValid);

    errorCode = mkapi::mkernel_is_valid_state(meshKernelId2, isValid);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_FALSE(isValid);
}

TEST(UndoTests, UnstructuredGrid)
{
    // setup
    //   1. create two curvilinear grids
    //   2. convert CLGs to unstructured grids.
    // test
    //   3. overwrite second grid using set function
    //   4. undo overwriting, check sizes of new mesh
    //   5. update second grid using add function
    //   6. undo grid update, check sizes of new mesh
    //   7. undo again, this time undo the conversion to unstructured.

    const int clg1Size = 8;
    const int clg2Size = 4;
    const double delta = 1.0;

    int meshKernelId1 = mkapi::mkernel_get_null_identifier();
    int meshKernelId2 = mkapi::mkernel_get_null_identifier();

    // Will be identical to meshKernelId1, to enable comparison after manipulation
    int meshKernelId3 = mkapi::mkernel_get_null_identifier();
    // Will be identical to meshKernelId2, to enable comparison after manipulation
    int meshKernelId4 = mkapi::mkernel_get_null_identifier();

    int errorCode;

    // Clear the meshkernel state and undo stack before starting the test.
    errorCode = mkapi::mkernel_clear_state();
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    errorCode = mkapi::mkernel_allocate_state(0, meshKernelId1);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    errorCode = mkapi::mkernel_allocate_state(0, meshKernelId2);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    errorCode = mkapi::mkernel_allocate_state(0, meshKernelId3);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    errorCode = mkapi::mkernel_allocate_state(0, meshKernelId4);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    errorCode = GenerateCurvilinearMesh(meshKernelId3, clg1Size, 0.5 * delta, static_cast<double>(clg2Size - 1) * delta);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    errorCode = GenerateCurvilinearMesh(meshKernelId1, clg1Size, 0.5 * delta, static_cast<double>(clg2Size - 1) * delta);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    // Convert the curvilinear grid id1 to an unstructured grid.
    errorCode = mkapi::mkernel_curvilinear_convert_to_mesh2d(meshKernelId1);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    // Generate the curvilinear grids, later to be converted to unstructured grids
    errorCode = GenerateCurvilinearMesh(meshKernelId4, clg2Size, delta, 0.0);
    errorCode = GenerateCurvilinearMesh(meshKernelId2, clg2Size, delta, 0.0);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    // Convert the curvilinear grid id2 to an unstructured grid.
    errorCode = mkapi::mkernel_curvilinear_convert_to_mesh2d(meshKernelId2);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    //--------------------------------
    // Execute the test.

    std::vector<double> xCoords(static_cast<size_t>(clg1Size * clg1Size));
    std::vector<double> yCoords(static_cast<size_t>(clg1Size * clg1Size));
    std::vector<int> edges(static_cast<size_t>(4 * clg1Size * (clg1Size - 1)));

    mkapi::Mesh2D mesh2d{};

    // Check dimensions before changing mesh.
    errorCode = mkapi::mkernel_mesh2d_get_dimensions(meshKernelId2, mesh2d);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(mesh2d.num_nodes, clg2Size * clg2Size);
    EXPECT_EQ(mesh2d.num_edges, 2 * (clg2Size - 1) * clg2Size);

    errorCode = mkapi::mkernel_mesh2d_get_dimensions(meshKernelId1, mesh2d);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(mesh2d.num_nodes, clg1Size * clg1Size);
    EXPECT_EQ(mesh2d.num_edges, 2 * (clg1Size - 1) * clg1Size);
    mesh2d.node_x = xCoords.data();
    mesh2d.node_y = yCoords.data();
    mesh2d.edge_nodes = edges.data();

    errorCode = mkapi::mkernel_mesh2d_get_node_edge_data(meshKernelId1, mesh2d);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    // Overwrite meshKernelId2 with the contents of mesh2d.
    errorCode = mkapi::mkernel_mesh2d_set(meshKernelId2, mesh2d);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    // Check dimensions after changing mesh.
    errorCode = mkapi::mkernel_mesh2d_get_dimensions(meshKernelId2, mesh2d);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(mesh2d.num_nodes, clg1Size * clg1Size);
    EXPECT_EQ(mesh2d.num_edges, 2 * (clg1Size - 1) * clg1Size);

    bool didUndo = false;
    int undoId = mkapi::mkernel_get_null_identifier();
    int committedCount = 0;
    int restoredCount = 0;

    errorCode = mkapi::mkernel_undo_state_count_for_id(meshKernelId2, committedCount, restoredCount);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    // 1. Creation of CLG, 2. conversion to unstructured, and 3. mesh2d_set
    EXPECT_EQ(committedCount, 3);
    EXPECT_EQ(restoredCount, 0);

    // Undo the mesh2d_set
    errorCode = mkapi::mkernel_undo_state(didUndo, undoId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndo);
    // meshKernelId2 was the last id set.
    EXPECT_EQ(undoId, meshKernelId2);

    // Check dimensions after undo mesh, should be same as before the last mesh2d_set.
    errorCode = mkapi::mkernel_mesh2d_get_dimensions(meshKernelId2, mesh2d);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(mesh2d.num_nodes, clg2Size * clg2Size);
    EXPECT_EQ(mesh2d.num_edges, 2 * (clg2Size - 1) * clg2Size);

    errorCode = mkapi::mkernel_mesh2d_get_dimensions(meshKernelId1, mesh2d);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    errorCode = mkapi::mkernel_mesh2d_get_node_edge_data(meshKernelId1, mesh2d);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    // Set the contents of mesh2d meshKernelId2.
    errorCode = mkapi::mkernel_mesh2d_set(meshKernelId2, mesh2d);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    // Compare grid1 with grid2, should be the same at this point.
    EXPECT_TRUE(CompareUnstructuredGrids(meshKernelId1, meshKernelId2));

    // Undo the mesh2d_set
    errorCode = mkapi::mkernel_undo_state(didUndo, undoId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndo);
    // meshKernelId2 was the last id set.
    EXPECT_EQ(undoId, meshKernelId2);

    // At this point meshKernelId1 should be restored to its original value

    // Add the contents of mesh2d (which are retrieved from meshKernelId1) meshKernelId2.
    errorCode = mkapi::mkernel_mesh2d_add(meshKernelId2, mesh2d);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    // Check dimensions after changing mesh.
    errorCode = mkapi::mkernel_mesh2d_get_dimensions(meshKernelId2, mesh2d);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(mesh2d.num_nodes, clg1Size * clg1Size + clg2Size * clg2Size);
    EXPECT_EQ(mesh2d.num_edges, 2 * (clg1Size - 1) * clg1Size + 2 * (clg2Size - 1) * clg2Size);

    // Undo the mesh2d_add
    errorCode = mkapi::mkernel_undo_state(didUndo, undoId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndo);
    // meshKernelId2 was the last id set.
    EXPECT_EQ(undoId, meshKernelId2);

    // Check dimensions after changing mesh.
    errorCode = mkapi::mkernel_mesh2d_get_dimensions(meshKernelId2, mesh2d);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(mesh2d.num_nodes, clg2Size * clg2Size);
    EXPECT_EQ(mesh2d.num_edges, 2 * (clg2Size - 1) * clg2Size);

    // Undo the conversion to unstructured
    undoId = mkapi::mkernel_get_null_identifier();
    errorCode = mkapi::mkernel_undo_state(didUndo, undoId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndo);
    // meshKernelId2 was the last id set.
    EXPECT_EQ(undoId, meshKernelId2);

    mkapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkapi::mkernel_curvilinear_get_dimensions(meshKernelId2, curvilinearGrid);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(curvilinearGrid.num_n, clg2Size);
    EXPECT_EQ(curvilinearGrid.num_m, clg2Size);

    // Compare the curvilinear grids
    EXPECT_TRUE(CompareCurvilinearGrids(meshKernelId2, meshKernelId4));

    errorCode = mkapi::mkernel_mesh2d_get_dimensions(meshKernelId2, mesh2d);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(mesh2d.num_nodes, mkapi::mkernel_get_null_identifier());
    EXPECT_EQ(mesh2d.num_edges, mkapi::mkernel_get_null_identifier());

    //--------------------------------
    // Finalise the test, remove all mesh kernel state objects and undo actions
    errorCode = mkapi::mkernel_expunge_state(meshKernelId2);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    // The mesh kernel state object and all undo information should have be removed.
    errorCode = mkapi::mkernel_undo_state_count_for_id(meshKernelId2, committedCount, restoredCount);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(committedCount, 0);
    EXPECT_EQ(restoredCount, 0);

    errorCode = mkapi::mkernel_expunge_state(meshKernelId1);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    // The mesh kernel state object and all undo information should have be removed.
    errorCode = mkapi::mkernel_undo_state_count_for_id(meshKernelId1, committedCount, restoredCount);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(committedCount, 0);
    EXPECT_EQ(restoredCount, 0);
}

// TODO
// add tests for
// 1. check that the SetMesh2dApiNodeEdgeData is correct after undoing a mesh change.
//    basically want to know if the correct nodes and edges (valid and invalid) are retrieved.
//
// 2. more combinations of mk api calls with undo.
