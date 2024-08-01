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

    // Execute
    int errorCode = mkapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    return errorCode;
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

    // Undo the curvilinear grid generation for the last grid generated.
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
    // The clg should be restored
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

    // Undo the deallocation
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

    // Undo the conversion of meshKernelId2 from curvilienar to unstructured.
    errorCode = mkapi::mkernel_undo_state(didUndo, undoId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndo);
    EXPECT_EQ(undoId, meshKernelId2);

    // The unstructured grid is no longer available
    errorCode = mkapi::mkernel_mesh2d_get_dimensions(meshKernelId2, mesh2d);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(mesh2d.num_nodes, mkapi::mkernel_get_null_identifier());
    EXPECT_EQ(mesh2d.num_edges, mkapi::mkernel_get_null_identifier());

    // Expunge the grid and its undo data from the api.
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
    //   4. undo grid update, check sizes of new mesh
    //   5. undo again, this time undo the conversion to unstructured.

    const int clg1Size = 8;
    const int clg2Size = 4;
    const double delta = 1.0;

    int meshKernelId1 = mkapi::mkernel_get_null_identifier();
    int meshKernelId2 = mkapi::mkernel_get_null_identifier();
    int errorCode;

    errorCode = mkapi::mkernel_allocate_state(0, meshKernelId1);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    errorCode = mkapi::mkernel_allocate_state(0, meshKernelId2);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    errorCode = GenerateCurvilinearMesh(meshKernelId1, clg1Size, 0.5 * delta, static_cast<double>(clg2Size - 1) * delta);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    // Convert the curvilinear grid id1 to an unstructured grid.
    errorCode = mkapi::mkernel_curvilinear_convert_to_mesh2d(meshKernelId1);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    // Generate the curvilinear grids, later to be converted to unstructured grids
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

    // Add the contents of mesh2d meshKernelId2.
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

    errorCode = mkapi::mkernel_mesh2d_get_dimensions(meshKernelId2, mesh2d);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(mesh2d.num_nodes, mkapi::mkernel_get_null_identifier());
    EXPECT_EQ(mesh2d.num_edges, mkapi::mkernel_get_null_identifier());

    //--------------------------------
    // Finalise the test

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
