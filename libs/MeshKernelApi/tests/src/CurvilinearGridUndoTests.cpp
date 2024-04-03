#include <gtest/gtest.h>

#include "MeshKernel/Parameters.hpp"
#include "MeshKernelApi/Mesh2D.hpp"
#include "MeshKernelApi/MeshKernel.hpp"

#include "CartesianApiTestFixture.hpp"
#include "TestUtils/MakeCurvilinearGrids.hpp"
#include "TestUtils/MakeMeshes.hpp"

TEST (CurvilinearGridUndoTests, DeleteNode)
{

    // Prepare
    int meshKernelId;
    auto errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);


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


    errorCode = meshkernelapi::mkernel_curvilinear_delete_node (meshKernelId, 5.0, 5.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    for (int i = 0; i < curvilinearGrid.num_m * curvilinearGrid.num_n; ++i)
    {
        if (i == deletedNodeIndex)
        {
            EXPECT_NEAR (curvilinearGrid.node_x [i], NullValue, tolerance);
            EXPECT_NEAR (curvilinearGrid.node_y [i], NullValue, tolerance);
        }
        else
        {
            EXPECT_NEAR (curvilinearGrid.node_x [i], originalNodeX [i], tolerance);
            EXPECT_NEAR (curvilinearGrid.node_y [i], originalNodeY [i], tolerance);
        }
    }


    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

}
