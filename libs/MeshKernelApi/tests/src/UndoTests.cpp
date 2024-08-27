#include <algorithm>
#include <gtest/gtest.h>
#include <random>

#include "MeshKernel/Operations.hpp"
#include "MeshKernel/Parameters.hpp"

#include "MeshKernelApi/BoundingBox.hpp"
#include "MeshKernelApi/Mesh2D.hpp"
#include "MeshKernelApi/MeshKernel.hpp"

#include "TestUtils/MakeMeshes.hpp"

// namespace aliases
namespace mk = meshkernel;
namespace mkapi = meshkernelapi;

// Helper function to generate a curvilinear grid for a meshKernelId
int GenerateCurvilinearMesh(const int meshKernelId,
                            const int nodesX, const int nodesY,
                            const double deltaX, const double deltaY,
                            const double originX, const double originY)
{
    meshkernel::MakeGridParameters makeGridParameters;

    // num_columns and num_rows indicate number of elements in each direction, so value = nodes - 1
    makeGridParameters.num_columns = nodesX - 1;
    makeGridParameters.num_rows = nodesY - 1;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = originX;
    makeGridParameters.origin_y = originY;
    makeGridParameters.block_size_x = deltaX;
    makeGridParameters.block_size_y = deltaY;

    // Generate curvilinear grid.
    int errorCode = mkapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    return errorCode;
}

// Helper function to generate a curvilinear grid for a meshKernelId
int GenerateCurvilinearMesh(const int meshKernelId, const int nodes, const double delta, const double origin)
{
    return GenerateCurvilinearMesh(meshKernelId, nodes, nodes, delta, delta, origin, origin);
}

// Helper function to compare two curvilinear grids
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

    if (mk::ExitCode::Success != errorCode)
    {
        return false;
    }

    size_t count = 0;

    for (int i = 0; i < curvilinearGrid1.num_n; ++i)
    {
        for (int j = 0; j < curvilinearGrid1.num_m; ++j)
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

// Helper function to compare two unstructured grids
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

    if (mk::ExitCode::Success != errorCode)
    {
        return false;
    }

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

// Helper function to compare the undo state for a meshKernelId
bool CheckUndoStateCount(const int meshKernelId, const int expectedCommitted, const int expectedRestored)
{
    int committedCount;
    int restoredCount;
    int errorCode = mkapi::mkernel_undo_state_count_for_id(meshKernelId, committedCount, restoredCount);

    return errorCode == mk::ExitCode::Success && committedCount == expectedCommitted && restoredCount == expectedRestored;
}

TEST(UndoTests, BasicAllocationDeallocationTest)
{

    // Two mesh kernel ids will be created
    // checks made on their validity at different stages of the test
    // one of the id's will be deallocated
    // further checks on validity

    int meshKernelId1 = meshkernel::constants::missing::intValue;
    int meshKernelId2 = meshkernel::constants::missing::intValue;
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
    // meshKernelId2 should still be valid
    isValid = false;
    errorCode = mkapi::mkernel_is_valid_state(meshKernelId2, isValid);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_TRUE(isValid);

    bool didUndo = false;
    int undoId = meshkernel::constants::missing::intValue;

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
    EXPECT_EQ(undoId, meshkernel::constants::missing::intValue);

    // Finalise the test, remove all mesh kernel state objects and undo actions
    errorCode = mkapi::mkernel_expunge_state(meshKernelId2);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    errorCode = mkapi::mkernel_expunge_state(meshKernelId1);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
}

TEST(UndoTests, CurvilinearGridManipulationTests)
{
    // set up
    // 1. generate curvilienar grids
    // 2. check sizes of generated CLG's
    // 3. undo generation of mkId1, and check size of CLG's (should be null)
    // 4. deallocate the state for mkId2.

    const int clg1Size = 4;
    const int clg2Size = 8;

    int meshKernelId1 = meshkernel::constants::missing::intValue;
    int meshKernelId2 = meshkernel::constants::missing::intValue;
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
    mkapi::CurvilinearGrid curvilinearGrid0{};
    errorCode = mkapi::mkernel_curvilinear_get_dimensions(meshKernelId2, curvilinearGrid0);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(curvilinearGrid0.num_n, clg2Size);
    EXPECT_EQ(curvilinearGrid0.num_m, clg2Size);

    mkapi::CurvilinearGrid curvilinearGrid1{};
    errorCode = mkapi::mkernel_curvilinear_get_dimensions(meshKernelId1, curvilinearGrid1);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(curvilinearGrid1.num_n, clg1Size);
    EXPECT_EQ(curvilinearGrid1.num_m, clg1Size);

    bool didUndo = false;
    int undoId = meshkernel::constants::missing::intValue;

    // Undo the curvilinear grid generation for the last grid generated, meshKernelId1
    errorCode = mkapi::mkernel_undo_state(didUndo, undoId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndo);
    // meshKernelId1 was the last mesh generated.
    EXPECT_EQ(undoId, meshKernelId1);

    // Getting the size of this grid should be an error and the dimensions the null value.
    mkapi::CurvilinearGrid curvilinearGrid2{};
    errorCode = mkapi::mkernel_curvilinear_get_dimensions(meshKernelId1, curvilinearGrid2);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(curvilinearGrid2.num_n, 0);
    EXPECT_EQ(curvilinearGrid2.num_m, 0);

    // Undo the undo of the clg generation
    // The clg for meshKernelId1 should be restored
    errorCode = mkapi::mkernel_redo_state(didUndo, undoId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(undoId, meshKernelId1);

    // With the correct size
    mkapi::CurvilinearGrid curvilinearGrid3{};
    errorCode = mkapi::mkernel_curvilinear_get_dimensions(meshKernelId1, curvilinearGrid3);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(curvilinearGrid3.num_n, clg1Size);
    EXPECT_EQ(curvilinearGrid3.num_m, clg1Size);

    // Deallocate clg 2
    errorCode = mkapi::mkernel_deallocate_state(meshKernelId2);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    mkapi::CurvilinearGrid curvilinearGrid4{};
    errorCode = mkapi::mkernel_curvilinear_get_dimensions(meshKernelId2, curvilinearGrid4);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(curvilinearGrid4.num_n, 0);
    EXPECT_EQ(curvilinearGrid4.num_m, 0);

    // Check size of meshKernelId1 is not touched
    mkapi::CurvilinearGrid curvilinear5{};
    errorCode = mkapi::mkernel_curvilinear_get_dimensions(meshKernelId1, curvilinear5);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(curvilinear5.num_n, clg1Size);
    EXPECT_EQ(curvilinear5.num_m, clg1Size);

    // Undo the deallocation or meshKernelId2
    errorCode = mkapi::mkernel_undo_state(didUndo, undoId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndo);
    EXPECT_EQ(undoId, meshKernelId2);

    // The grid should now have the correct size
    mkapi::CurvilinearGrid curvilinearGrid6{};
    errorCode = mkapi::mkernel_curvilinear_get_dimensions(meshKernelId2, curvilinearGrid6);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(curvilinearGrid6.num_n, clg2Size);
    EXPECT_EQ(curvilinearGrid6.num_m, clg2Size);

    // Convert the curvilinear grid id2 to an unstructured grid.
    errorCode = mkapi::mkernel_curvilinear_convert_to_mesh2d(meshKernelId2);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    // Check size of meshKernelId1 is not touched
    mkapi::CurvilinearGrid curvilinearGrid7{};
    errorCode = mkapi::mkernel_curvilinear_get_dimensions(meshKernelId1, curvilinearGrid7);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(curvilinearGrid7.num_n, clg1Size);
    EXPECT_EQ(curvilinearGrid7.num_m, clg1Size);

    mkapi::Mesh2D mesh2d0{};
    errorCode = mkapi::mkernel_mesh2d_get_dimensions(meshKernelId2, mesh2d0);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(mesh2d0.num_nodes, clg2Size * clg2Size);
    EXPECT_EQ(mesh2d0.num_edges, 2 * (clg2Size - 1) * clg2Size);

    // Now that the grid has been converted to an unstructured grid the dimensions of the curvilienar grid are not
    mkapi::CurvilinearGrid curvilinearGrid8{};
    errorCode = mkapi::mkernel_curvilinear_get_dimensions(meshKernelId2, curvilinearGrid8);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(curvilinearGrid8.num_n, 0);
    EXPECT_EQ(curvilinearGrid8.num_m, 0);

    // Undo the conversion of meshKernelId2 from curvilinear to unstructured.
    errorCode = mkapi::mkernel_undo_state(didUndo, undoId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndo);
    EXPECT_EQ(undoId, meshKernelId2);

    // The unstructured grid is no longer available
    mkapi::Mesh2D mesh2d1{};
    errorCode = mkapi::mkernel_mesh2d_get_dimensions(meshKernelId2, mesh2d1);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(mesh2d1.num_nodes, 0);
    EXPECT_EQ(mesh2d1.num_edges, 0);

    // Check size of meshKernelId1 is not touched
    mkapi::CurvilinearGrid curvilinearGrid9{};
    errorCode = mkapi::mkernel_curvilinear_get_dimensions(meshKernelId1, curvilinearGrid9);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(curvilinearGrid9.num_n, clg1Size);
    EXPECT_EQ(curvilinearGrid9.num_m, clg1Size);

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

    int meshKernelId1 = meshkernel::constants::missing::intValue;
    int meshKernelId2 = meshkernel::constants::missing::intValue;

    // Will be identical to meshKernelId1, to enable comparison after manipulation
    int meshKernelId3 = meshkernel::constants::missing::intValue;
    // Will be identical to meshKernelId2, to enable comparison after manipulation
    int meshKernelId4 = meshkernel::constants::missing::intValue;

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
    errorCode = GenerateCurvilinearMesh(meshKernelId2, clg2Size, delta, 0.0);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    // Generate an identical copy of meshKernelId2 to enable easy comparing later
    errorCode = GenerateCurvilinearMesh(meshKernelId4, clg2Size, delta, 0.0);
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
    int undoId = meshkernel::constants::missing::intValue;

    // 1. Creation of CLG, 2. conversion to unstructured, and 3. mesh2d_set
    EXPECT_TRUE(CheckUndoStateCount(meshKernelId2, 3, 0));

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
    undoId = meshkernel::constants::missing::intValue;
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
    EXPECT_EQ(mesh2d.num_nodes, 0);
    EXPECT_EQ(mesh2d.num_edges, 0);

    //--------------------------------
    // Finalise the test, remove all mesh kernel state objects and undo actions
    errorCode = mkapi::mkernel_expunge_state(meshKernelId2);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    // The mesh kernel state object and all undo information should have be removed.
    EXPECT_TRUE(CheckUndoStateCount(meshKernelId2, 0, 0));

    errorCode = mkapi::mkernel_expunge_state(meshKernelId1);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    // The mesh kernel state object and all undo information should have be removed.
    EXPECT_TRUE(CheckUndoStateCount(meshKernelId1, 0, 0));
}

TEST(UndoTests, UnstructuredGridConnection)
{
    // setup
    //   1. create curvilinear grids
    //   2. convert CLG's to unstructured grids
    // test
    //   3. insert node into grid with mkId1 (1 undo action), connect surrounding corner nodes to this new node (4 undo actions)
    //   4. undo the connecting of the nodes and the node insertion (this is to create holes in the node and edge lists)
    //   5. get the mesh data for mkId1
    //   6. connect this mesh data to mkId2, and check result
    //   7. undo connection, check the nodal values
    //   8. undo conversion, mkId2 should now be as it were originally
    //   9. redo conversion and connection, check the nodal values again
    // finalise
    //  10. remove all data from the api.

    const int clg1SizeX = 4;
    const int clg1SizeY = 4;
    const int clg2SizeX = 3;
    const int clg2SizeY = 9;
    const double clg1DeltaX = 4.0;
    const double clg1DeltaY = 4.0;
    const double clg2DeltaX = 1.0;
    const double clg2DeltaY = 1.0;

    const double fraction = 0.03125;

    const double clg1OriginX = 0.0;
    const double clg1OriginY = 0.0;
    const double clg2OriginX = static_cast<double>(clg1SizeX - 1) * clg1DeltaX + fraction * clg2DeltaX;
    const double clg2OriginY = 0.0;

    int meshKernelId1 = meshkernel::constants::missing::intValue;
    int meshKernelId2 = meshkernel::constants::missing::intValue;
    // Will contain the same grid as meshKernelId2 so that it can be compared later
    int meshKernelId3 = meshkernel::constants::missing::intValue;
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

    // Generate the curvilinear grids, later to be converted to unstructured grids
    errorCode = GenerateCurvilinearMesh(meshKernelId1, clg1SizeX, clg1SizeY, clg1DeltaX, clg1DeltaY, clg1OriginX, clg1OriginY);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    errorCode = GenerateCurvilinearMesh(meshKernelId2, clg2SizeX, clg2SizeY, clg2DeltaX, clg2DeltaY, clg2OriginX, clg2OriginY);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    errorCode = GenerateCurvilinearMesh(meshKernelId3, clg2SizeX, clg2SizeY, clg2DeltaX, clg2DeltaY, clg2OriginX, clg2OriginY);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    // Convert the curvilinear grid id1 to an unstructured grid.
    errorCode = mkapi::mkernel_curvilinear_convert_to_mesh2d(meshKernelId1);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    // Convert the curvilinear grid id2 to an unstructured grid.
    errorCode = mkapi::mkernel_curvilinear_convert_to_mesh2d(meshKernelId2);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    //--------------------------------
    // Insert node in middle of an element and connect the surrounding 4 nodes to it.
    // The undo these actions.
    // This is to create some unused data in the node and edge lists.

    // insert node to mesh1
    int newNodeId;
    errorCode = mkapi::mkernel_mesh2d_insert_node(meshKernelId1, 0.5, 0.5, newNodeId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    int node1 = meshkernel::constants::missing::intValue;
    int node2 = meshkernel::constants::missing::intValue;
    int node3 = meshkernel::constants::missing::intValue;
    int node4 = meshkernel::constants::missing::intValue;

    mkapi::BoundingBox boundingBox{-0.1, -0.1, 1.1, 1.1};

    errorCode = mkapi::mkernel_mesh2d_get_location_index(meshKernelId1, 0.0, 0.0,
                                                         1 /*Location::Node*/,
                                                         boundingBox,
                                                         node1);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    errorCode = mkapi::mkernel_mesh2d_get_location_index(meshKernelId1, 1.0, 0.0,
                                                         1 /*Location::Node*/,
                                                         boundingBox,
                                                         node2);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    errorCode = mkapi::mkernel_mesh2d_get_location_index(meshKernelId1, 1.0, 1.0,
                                                         1 /*Location::Node*/,
                                                         boundingBox,
                                                         node3);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    errorCode = mkapi::mkernel_mesh2d_get_location_index(meshKernelId1, 0.0, 1.0,
                                                         1 /*Location::Node*/,
                                                         boundingBox,
                                                         node4);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    int edge1 = meshkernel::constants::missing::intValue;
    int edge2 = meshkernel::constants::missing::intValue;
    int edge3 = meshkernel::constants::missing::intValue;
    int edge4 = meshkernel::constants::missing::intValue;

    errorCode = mkapi::mkernel_mesh2d_insert_edge(meshKernelId1, node1, newNodeId, edge1);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    errorCode = mkapi::mkernel_mesh2d_insert_edge(meshKernelId1, node2, newNodeId, edge2);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    errorCode = mkapi::mkernel_mesh2d_insert_edge(meshKernelId1, node3, newNodeId, edge3);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    errorCode = mkapi::mkernel_mesh2d_insert_edge(meshKernelId1, node4, newNodeId, edge4);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    // 7 undo actions
    // 1 creation, 1 conversion, 1 insert node and 4 insert edges
    EXPECT_TRUE(CheckUndoStateCount(meshKernelId1, 7, 0));

    // now, 5 undos to remove. The adding of 4 edges and 1 inserted node.
    for (int i = 1; i <= 5; ++i)
    {
        bool didUndo = false;
        int undoId = meshkernel::constants::missing::intValue;

        // Undo the deallocation
        errorCode = mkapi::mkernel_undo_state(didUndo, undoId);
        ASSERT_EQ(mk::ExitCode::Success, errorCode);
        EXPECT_TRUE(didUndo);
        EXPECT_EQ(undoId, meshKernelId1);
    }

    EXPECT_TRUE(CheckUndoStateCount(meshKernelId1, 2, 5));
    EXPECT_TRUE(CheckUndoStateCount(meshKernelId2, 2, 0));

    //--------------------------------
    // Now the nodes and edges of meshKernelId1, should contain holes from the undoing above.
    // Connect the meshes, overwriting meshKernelId2 and check the connected mesh is correct.

    mkapi::Mesh2D mesh2d{};
    errorCode = mkapi::mkernel_mesh2d_get_dimensions(meshKernelId1, mesh2d);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    std::vector<double> xNodes(mesh2d.num_nodes);
    std::vector<double> yNodes(mesh2d.num_nodes);
    std::vector<int> edges(2 * mesh2d.num_edges);

    mesh2d.node_x = xNodes.data();
    mesh2d.node_y = yNodes.data();
    mesh2d.edge_nodes = edges.data();

    errorCode = mkapi::mkernel_mesh2d_get_node_edge_data(meshKernelId1, mesh2d);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    // Connect the meshes.
    errorCode = mkapi::mkernel_mesh2d_connect_meshes(meshKernelId2, mesh2d, 0.1);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    errorCode = mkapi::mkernel_mesh2d_get_dimensions(meshKernelId2, mesh2d);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    std::vector<double> xNodesConnected(mesh2d.num_nodes);
    std::vector<double> yNodesConnected(mesh2d.num_nodes);
    std::vector<int> edgesConnected(2 * mesh2d.num_edges);

    mesh2d.node_x = xNodesConnected.data();
    mesh2d.node_y = yNodesConnected.data();
    mesh2d.edge_nodes = edgesConnected.data();

    errorCode = mkapi::mkernel_mesh2d_get_node_edge_data(meshKernelId2, mesh2d);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    const double nullValue = mk::constants::missing::doubleValue;

    // The null values are from the three common nodes of the two meshes being connected.
    std::vector<double> expectedNodesConnectedX{12.0, 13.0, 14.0, 12.0, 13.0, 14.0, 12.0, 13.0, 14.0, 12.0, 13.0, 14.0, 12.0, 13.0, 14.0,
                                                12.0, 13.0, 14.0, 12.0, 13.0, 14.0, 12.0, 13.0, 14.0, 12.0, 13.0, 14.0, 0.0, 4.0, 8.0,
                                                nullValue, 0.0, 4.0, 8.0, nullValue, 0.0, 4.0, 8.0, nullValue, 0.0, 4.0, 8.0, 12.0, 8.0, 8.0};

    // Shift the nodes of subdomain 1 slightly to match the generated mesh
    for (size_t i = 0; i < 27; ++i)
    {
        expectedNodesConnectedX[i] += fraction * clg2DeltaX;
    }

    std::vector<double> expectedNodesConnectedY{0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 4.0, 4.0, 4.0,
                                                5.0, 5.0, 5.0, 6.0, 6.0, 6.0, 7.0, 7.0, 7.0, 8.0, 8.0, 8.0, 0.0, 0.0, 0.0,
                                                nullValue, 4.0, 4.0, 4.0, nullValue, 8.0, 8.0, 8.0, nullValue, 12.0, 12.0, 12.0, 12.0, 2.0, 6.0};

    ASSERT_EQ(mesh2d.num_nodes, 45);
    ASSERT_EQ(mesh2d.num_edges, 84);

    for (size_t i = 0; i < static_cast<size_t>(mesh2d.num_nodes); ++i)
    {
        EXPECT_EQ(expectedNodesConnectedX[i], xNodesConnected[i]);
        EXPECT_EQ(expectedNodesConnectedY[i], yNodesConnected[i]);
    }

    //--------------------------------
    // Undo the connection of the meshes

    bool didUndo = false;
    int undoId = meshkernel::constants::missing::intValue;

    EXPECT_TRUE(CheckUndoStateCount(meshKernelId1, 2, 5));
    EXPECT_TRUE(CheckUndoStateCount(meshKernelId2, 3, 0));

    errorCode = mkapi::mkernel_undo_state(didUndo, undoId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndo);
    EXPECT_EQ(undoId, meshKernelId2);

    std::vector<double> xNodesDisconnected(mesh2d.num_nodes);
    std::vector<double> yNodesDisconnected(mesh2d.num_nodes);
    std::vector<int> edgesDisconnected(2 * mesh2d.num_edges);

    errorCode = mkapi::mkernel_mesh2d_get_dimensions(meshKernelId2, mesh2d);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    mesh2d.node_x = xNodesDisconnected.data();
    mesh2d.node_y = yNodesDisconnected.data();
    mesh2d.edge_nodes = edgesDisconnected.data();

    errorCode = mkapi::mkernel_mesh2d_get_node_edge_data(meshKernelId2, mesh2d);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    ASSERT_EQ(mesh2d.num_nodes, clg2SizeX * clg2SizeY /* 27 */);
    ASSERT_EQ(mesh2d.num_edges, (clg2SizeX - 1) * clg2SizeY + clg2SizeX * (clg2SizeY - 1) /* 42 */);

    std::vector<double> expectedNodesDisconnectedX{12.0, 13.0, 14.0, 12.0, 13.0, 14.0, 12.0, 13.0, 14.0,
                                                   12.0, 13.0, 14.0, 12.0, 13.0, 14.0, 12.0, 13.0, 14.0,
                                                   12.0, 13.0, 14.0, 12.0, 13.0, 14.0, 12.0, 13.0, 14.0};

    std::vector<double> expectedNodesDisconnectedY{0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0,
                                                   3.0, 3.0, 3.0, 4.0, 4.0, 4.0, 5.0, 5.0, 5.0,
                                                   6.0, 6.0, 6.0, 7.0, 7.0, 7.0, 8.0, 8.0, 8.0};

    // Shift the nodes of subdomain 1 slightly to match the generated mesh
    std::ranges::for_each(expectedNodesDisconnectedX, [fraction, clg2DeltaX](double& value)
                          { value += fraction * clg2DeltaX; });

    for (size_t i = 0; i < static_cast<size_t>(mesh2d.num_nodes); ++i)
    {
        EXPECT_EQ(expectedNodesDisconnectedX[i], xNodesDisconnected[i]);
        EXPECT_EQ(expectedNodesDisconnectedY[i], yNodesDisconnected[i]);
    }

    didUndo = false;
    undoId = meshkernel::constants::missing::intValue;

    // Undo mesh conversion of curvilinear to mesh2d
    errorCode = mkapi::mkernel_undo_state(didUndo, undoId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndo);
    EXPECT_EQ(undoId, meshKernelId2);

    EXPECT_TRUE(CheckUndoStateCount(meshKernelId1, 2, 5));
    EXPECT_TRUE(CheckUndoStateCount(meshKernelId2, 1, 2));

    // Check that the curvilienar grid match the original grid.
    EXPECT_TRUE(CompareCurvilinearGrids(meshKernelId2, meshKernelId3));

    didUndo = false;
    undoId = meshkernel::constants::missing::intValue;

    // Redo mesh conversion of curvilinear to mesh2d
    errorCode = mkapi::mkernel_redo_state(didUndo, undoId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndo);
    EXPECT_EQ(undoId, meshKernelId2);

    // Redo connecting meshes
    errorCode = mkapi::mkernel_redo_state(didUndo, undoId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndo);
    EXPECT_EQ(undoId, meshKernelId2);

    EXPECT_TRUE(CheckUndoStateCount(meshKernelId1, 2, 5));
    EXPECT_TRUE(CheckUndoStateCount(meshKernelId2, 3, 0));

    // Check the dimensions are correct
    errorCode = mkapi::mkernel_mesh2d_get_dimensions(meshKernelId2, mesh2d);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    ASSERT_EQ(mesh2d.num_nodes, static_cast<int>(xNodesConnected.size()));
    ASSERT_EQ(2 * mesh2d.num_edges, static_cast<int>(edgesConnected.size()));

    mesh2d.node_x = xNodesConnected.data();
    mesh2d.node_y = yNodesConnected.data();
    mesh2d.edge_nodes = edgesConnected.data();

    // Now check the data is correct
    errorCode = mkapi::mkernel_mesh2d_get_node_edge_data(meshKernelId2, mesh2d);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    // The nodes should be same as above directly after connecting the meshes.
    for (size_t i = 0; i < static_cast<size_t>(mesh2d.num_nodes); ++i)
    {
        EXPECT_EQ(expectedNodesConnectedX[i], xNodesConnected[i]);
        EXPECT_EQ(expectedNodesConnectedY[i], yNodesConnected[i]);
    }

    //--------------------------------
    // Finalise the test, remove all mesh kernel state objects and undo actions
    errorCode = mkapi::mkernel_expunge_state(meshKernelId2);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    // The mesh kernel state object and all undo information should have be removed.
    EXPECT_TRUE(CheckUndoStateCount(meshKernelId2, 0, 0));

    errorCode = mkapi::mkernel_expunge_state(meshKernelId1);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    // The mesh kernel state object and all undo information should have be removed.
    EXPECT_TRUE(CheckUndoStateCount(meshKernelId1, 0, 0));
}

TEST(UndoTests, SetUndoStackSize)
{

    // Two mesh kernel ids will be created
    // checks made on their validity at different stages of the test
    // one of the id's will be deallocated
    // further checks on validity

    int meshKernelId1 = meshkernel::constants::missing::intValue;
    int meshKernelId2 = meshkernel::constants::missing::intValue;
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

    errorCode = mkapi::mkernel_set_undo_size(0);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

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
    // meshKernelId2 should still be valid
    isValid = false;
    errorCode = mkapi::mkernel_is_valid_state(meshKernelId2, isValid);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_TRUE(isValid);

    int commitedSize = -1;
    int restoredSize = -1;

    errorCode = mkapi::mkernel_undo_state_count(commitedSize, restoredSize);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_EQ(commitedSize, 0);
    EXPECT_EQ(restoredSize, 0);

    bool didUndo = false;
    int undoId = meshkernel::constants::missing::intValue;

    // Undo the deallocation
    errorCode = mkapi::mkernel_undo_state(didUndo, undoId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_FALSE(didUndo);
    EXPECT_EQ(undoId, meshkernel::constants::missing::intValue);
}
