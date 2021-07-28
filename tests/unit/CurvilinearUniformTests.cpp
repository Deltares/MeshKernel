#include <gtest/gtest.h>

#include <MeshKernel/CurvilinearGrid/CurvilinearGridCreateUniform.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernelApi/MakeMeshParameters.hpp>
#include <TestUtils/MakeCurvilinearGrids.hpp>

TEST(CurvilinearGrid, CurvilinearGridCreateUniform_WithPolygon_ShouldComputeCurvilinearGrid)
{
    // Setup
    std::vector<meshkernel::Point> polygonNodes{{0.5, 2.5},
                                                {2.5, 0.5},
                                                {5.5, 3.0},
                                                {3.5, 5.0},
                                                {0.5, 2.5}};

    const auto polygons = std::make_shared<meshkernel::Polygons>(polygonNodes, meshkernel::Projection::cartesian);

    meshkernelapi::MakeMeshParameters makeMeshParameters;
    makeMeshParameters.grid_type = 0;
    makeMeshParameters.angle = 0.0;
    makeMeshParameters.origin_x = 0.0;
    makeMeshParameters.origin_y = 0.0;
    makeMeshParameters.num_columns = 3;
    makeMeshParameters.num_rows = 3;
    makeMeshParameters.block_size_x = 1.0;
    makeMeshParameters.block_size_y = 1.0;

    // Execution
    meshkernel::CurvilinearGridCreateUniform const curvilinearGridCreateUniform(makeMeshParameters, meshkernel::Projection::cartesian);
    const auto curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(curvilinearGridCreateUniform.Compute(polygons, 0));

    // Assert, also invalid nodes and edges are included in the curvilinear grid
    auto const numValidNodes = CurvilinearGridCountValidNodes(curvilinearGrid);
    ASSERT_EQ(10, numValidNodes);
}

TEST(CurvilinearGrid, MakeCurvilinearInPolygonSpherical)
{
    // Setup
    std::vector<meshkernel::Point> polygonNodes{{302.002502, 472.130371},
                                                {144.501526, 253.128174},
                                                {368.752930, 112.876755},
                                                {707.755005, 358.879242},
                                                {301.252502, 471.380371},
                                                {302.002502, 472.130371}};

    const auto polygons = std::make_shared<meshkernel::Polygons>(polygonNodes, meshkernel::Projection::spherical);

    meshkernelapi::MakeMeshParameters makeMeshParameters;
    makeMeshParameters.grid_type = 0;
    makeMeshParameters.angle = 0.0;
    makeMeshParameters.origin_x = 0.0;
    makeMeshParameters.origin_y = 0.0;
    makeMeshParameters.num_columns = 3;
    makeMeshParameters.num_rows = 3;
    makeMeshParameters.block_size_x = 5000000.0; //resolution in meters (when using spherical coordinates distances are usually much larger)
    makeMeshParameters.block_size_y = 5000000.0;

    // Execution: function not producing grid points (points gets transformed in meters, therfore everything is outside)
    meshkernel::CurvilinearGridCreateUniform const curvilinearGridCreateUniform(makeMeshParameters, meshkernel::Projection::spherical);
    const auto curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(curvilinearGridCreateUniform.Compute(polygons, 0));

    // Assert
    auto const numValidNodes = CurvilinearGridCountValidNodes(curvilinearGrid);
    ASSERT_EQ(0, numValidNodes);
}

TEST(CurvilinearGrid, MakeCurvilinearInEmptyPolygonSpherical)
{
    //1 Setup

    meshkernelapi::MakeMeshParameters makeMeshParameters;
    makeMeshParameters.grid_type = 0;
    makeMeshParameters.angle = 0.0;
    makeMeshParameters.origin_x = 0.0;
    makeMeshParameters.origin_y = 0.0;
    makeMeshParameters.num_columns = 3;
    makeMeshParameters.num_rows = 3;
    makeMeshParameters.block_size_x = 5000000.0; //resolution in meters (when using spherical coordinates distances are usually much larger)
    makeMeshParameters.block_size_y = 5000000.0;

    // 2 Execution
    meshkernel::CurvilinearGridCreateUniform const curvilinearGridCreateUniform(makeMeshParameters, meshkernel::Projection::spherical);
    const auto [nodes, edges, gridIndices] = curvilinearGridCreateUniform.Compute().ConvertCurvilinearToNodesAndEdges();

    meshkernel::Mesh2D mesh(edges, nodes, meshkernel::Projection::spherical);

    // 3 Assert
    ASSERT_EQ(24, mesh.GetNumEdges());
    ASSERT_EQ(16, mesh.GetNumNodes());

    // x coordinate
    ASSERT_EQ(0.0, mesh.m_nodes[0].x);
    ASSERT_EQ(makeMeshParameters.block_size_x, mesh.m_nodes[1].x);
    ASSERT_EQ(makeMeshParameters.block_size_x * 2, mesh.m_nodes[2].x);
    ASSERT_EQ(makeMeshParameters.block_size_x * 3, mesh.m_nodes[3].x);

    ASSERT_EQ(0.0, mesh.m_nodes[4].x);
    ASSERT_EQ(makeMeshParameters.block_size_x, mesh.m_nodes[5].x);
    ASSERT_EQ(makeMeshParameters.block_size_x * 2, mesh.m_nodes[6].x);
    ASSERT_EQ(makeMeshParameters.block_size_x * 3, mesh.m_nodes[7].x);

    ASSERT_EQ(0.0, mesh.m_nodes[8].x);
    ASSERT_EQ(makeMeshParameters.block_size_x, mesh.m_nodes[9].x);
    ASSERT_EQ(makeMeshParameters.block_size_x * 2, mesh.m_nodes[10].x);
    ASSERT_EQ(makeMeshParameters.block_size_x * 3, mesh.m_nodes[11].x);

    // y coordinate
    ASSERT_EQ(0.0, mesh.m_nodes[0].y);
    ASSERT_EQ(0.0, mesh.m_nodes[1].y);
    ASSERT_EQ(0.0, mesh.m_nodes[2].y);
    ASSERT_EQ(0.0, mesh.m_nodes[3].y);

    ASSERT_EQ(makeMeshParameters.block_size_x, mesh.m_nodes[4].y);
    ASSERT_EQ(makeMeshParameters.block_size_x, mesh.m_nodes[5].y);
    ASSERT_EQ(makeMeshParameters.block_size_x, mesh.m_nodes[6].y);
    ASSERT_EQ(makeMeshParameters.block_size_x, mesh.m_nodes[7].y);

    ASSERT_EQ(3830222.2156113400, mesh.m_nodes[8].y);
    ASSERT_EQ(3830222.2156113400, mesh.m_nodes[9].y);
    ASSERT_EQ(3830222.2156113400, mesh.m_nodes[10].y);
    ASSERT_EQ(3830222.2156113400, mesh.m_nodes[11].y);
}

TEST(CurvilinearGrid, InsertFace_OnBottomLeft_ShouldInsertFace)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();

    // Execution
    curvilinearGrid->InsertFace({80009.0, 366937.0});

    // Assert the new coordinates
    constexpr double tolerance = 1e-12;
    ASSERT_NEAR(79913.595823791460, curvilinearGrid->m_gridNodes[0][0].x, tolerance);
    ASSERT_NEAR(79985.899758176762, curvilinearGrid->m_gridNodes[0][1].x, tolerance);
    ASSERT_NEAR(-999.00000000000000, curvilinearGrid->m_gridNodes[0][2].x, tolerance);

    ASSERT_NEAR(367046.61206756550, curvilinearGrid->m_gridNodes[0][0].y, tolerance);
    ASSERT_NEAR(367110.19327267300, curvilinearGrid->m_gridNodes[0][1].y, tolerance);
    ASSERT_NEAR(-999.0000000000000, curvilinearGrid->m_gridNodes[0][2].y, tolerance);
}

TEST(CurvilinearGrid, InsertFace_OnBottomRight_ShouldInsertFace)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();

    // Execution
    curvilinearGrid->InsertFace({80166.0, 366544.0});

    // Assert the new coordinates
    constexpr double tolerance = 1e-12;
    ASSERT_NEAR(80176.237835892549, curvilinearGrid->m_gridNodes[5][0].x, tolerance);
    ASSERT_NEAR(80259.193944680519, curvilinearGrid->m_gridNodes[5][1].x, tolerance);
    ASSERT_NEAR(-999.0000000000000, curvilinearGrid->m_gridNodes[5][2].x, tolerance);

    ASSERT_NEAR(366433.98982542212, curvilinearGrid->m_gridNodes[5][0].y, tolerance);
    ASSERT_NEAR(366433.75796857959, curvilinearGrid->m_gridNodes[5][1].y, tolerance);
    ASSERT_NEAR(-999.0000000000000, curvilinearGrid->m_gridNodes[5][2].y, tolerance);
}

TEST(CurvilinearGrid, InsertFace_OnTopLeft_ShouldInsertFace)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();

    // Execution
    curvilinearGrid->InsertFace({80612.0, 367407.0});

    // Assert the new coordinates
    constexpr double tolerance = 1e-12;
    ASSERT_NEAR(-999.0000000000000, curvilinearGrid->m_gridNodes[0][6].x, tolerance);
    ASSERT_NEAR(80385.178875081983, curvilinearGrid->m_gridNodes[0][7].x, tolerance);
    ASSERT_NEAR(80449.954021552709, curvilinearGrid->m_gridNodes[0][8].x, tolerance);

    ASSERT_NEAR(-999.0000000000000, curvilinearGrid->m_gridNodes[0][6].y, tolerance);
    ASSERT_NEAR(367549.55640742677, curvilinearGrid->m_gridNodes[0][7].y, tolerance);
    ASSERT_NEAR(367626.58734840894, curvilinearGrid->m_gridNodes[0][8].y, tolerance);
}

TEST(CurvilinearGrid, InsertFace_OnTopRight_ShouldInsertFace)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();

    // Execution
    curvilinearGrid->InsertFace({80870.0, 366541.0});

    // Assert the new coordinates
    constexpr double tolerance = 1e-12;
    ASSERT_NEAR(-999.0000000000000, curvilinearGrid->m_gridNodes[5][6].x, tolerance);
    ASSERT_NEAR(80846.719054016474, curvilinearGrid->m_gridNodes[5][7].x, tolerance);
    ASSERT_NEAR(80959.766315635425, curvilinearGrid->m_gridNodes[5][8].x, tolerance);

    ASSERT_NEAR(-999.0000000000000, curvilinearGrid->m_gridNodes[5][6].y, tolerance);
    ASSERT_NEAR(366346.42989900074, curvilinearGrid->m_gridNodes[5][7].y, tolerance);
    ASSERT_NEAR(366327.01674911042, curvilinearGrid->m_gridNodes[5][8].y, tolerance);
}

TEST(CurvilinearGrid, InsertFace_OnGridWithHoles_ShouldInsertFace)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGridWithMissingFaces();

    // Execution
    curvilinearGrid->InsertFace({80398.0, 366854.0});

    // Assert the new coordinates
    constexpr double tolerance = 1e-12;
    ASSERT_NEAR(80133.930743945661, curvilinearGrid->m_gridNodes[3][0].x, tolerance);
    ASSERT_NEAR(80226.579454943669, curvilinearGrid->m_gridNodes[3][1].x, tolerance);
    ASSERT_NEAR(80320.396380977647, curvilinearGrid->m_gridNodes[3][2].x, tolerance);
    ASSERT_NEAR(80447.108241179725, curvilinearGrid->m_gridNodes[3][3].x, tolerance);
    ASSERT_NEAR(80545.663179203286, curvilinearGrid->m_gridNodes[3][4].x, tolerance);
    ASSERT_NEAR(80623.151266424975, curvilinearGrid->m_gridNodes[3][5].x, tolerance);
    ASSERT_NEAR(80735.800935924854, curvilinearGrid->m_gridNodes[3][6].x, tolerance);
    ASSERT_NEAR(80848.959456291268, curvilinearGrid->m_gridNodes[3][7].x, tolerance);
    ASSERT_NEAR(80962.132385367571, curvilinearGrid->m_gridNodes[3][8].x, tolerance);

    ASSERT_NEAR(366629.99913221149, curvilinearGrid->m_gridNodes[3][0].y, tolerance);
    ASSERT_NEAR(366656.00122676318, curvilinearGrid->m_gridNodes[3][1].y, tolerance);
    ASSERT_NEAR(366679.12590409513, curvilinearGrid->m_gridNodes[3][2].y, tolerance);
    ASSERT_NEAR(366697.14301043766, curvilinearGrid->m_gridNodes[3][3].y, tolerance);
    ASSERT_NEAR(366725.32926280121, curvilinearGrid->m_gridNodes[3][4].y, tolerance);
    ASSERT_NEAR(366718.43748113938, curvilinearGrid->m_gridNodes[3][5].y, tolerance);
    ASSERT_NEAR(366717.43216295895, curvilinearGrid->m_gridNodes[3][6].y, tolerance);
    ASSERT_NEAR(366716.74631036510, curvilinearGrid->m_gridNodes[3][7].y, tolerance);
    ASSERT_NEAR(366718.10475241451, curvilinearGrid->m_gridNodes[3][8].y, tolerance);
}

TEST(CurvilinearGrid, DeleteNode_OnUniformGrid_ShouldDeleteNode)
{
    // Prepare
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();

    // Execute
    curvilinearGrid->DeleteNode({80398.0, 366854.0});

    // The number of nodes was 45 now is 44
    auto const numValidNodes = CurvilinearGridCountValidNodes(curvilinearGrid);
    ASSERT_EQ(numValidNodes, 44);
}
