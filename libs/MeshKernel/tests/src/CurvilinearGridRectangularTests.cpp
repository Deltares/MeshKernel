#include <gtest/gtest.h>

#include <MeshKernel/CurvilinearGrid/CurvilinearGridRectangular.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Polygons.hpp>
#include <TestUtils/MakeCurvilinearGrids.hpp>

using namespace meshkernel;

TEST(CurvilinearGridUniform, CurvilinearGridRectangular_WithPolygon_ShouldComputeCurvilinearGrid)
{
    // Setup
    std::vector<Point> polygonNodes{{0.5, 2.5},
                                    {2.5, 0.5},
                                    {5.5, 3.0},
                                    {3.5, 5.0},
                                    {0.5, 2.5}};

    const auto polygons = std::make_shared<Polygons>(polygonNodes, Projection::cartesian);

    const double angle = 0.0;
    const double blockSizeX = 1.0;
    const double blockSizeY = 1.0;

    // Execution
    CurvilinearGridRectangular const curvilinearGridRectangular(Projection::cartesian);
    const auto curvilinearGrid = curvilinearGridRectangular.Compute(angle,
                                                              blockSizeX,
                                                              blockSizeY,
                                                              polygons,
                                                              0);

    // Assert, also invalid nodes and edges are included in the curvilinear grid
    auto const numValidNodes = CurvilinearGridCountValidNodes(curvilinearGrid);
    ASSERT_EQ(9, numValidNodes);
}

TEST(CurvilinearGridUniform, MakeCurvilinearInPolygonSpherical)
{
    // Setup
    std::vector<Point> polygonNodes{{302.002502, 472.130371},
                                    {144.501526, 253.128174},
                                    {368.752930, 112.876755},
                                    {707.755005, 358.879242},
                                    {301.252502, 471.380371},
                                    {302.002502, 472.130371}};

    const auto polygons = std::make_shared<Polygons>(polygonNodes, Projection::spherical);

    const double angle = 0.0;
    const double blockSizeX = 500.0; // resolution in degree
    const double blockSizeY = 500.0;

    // Execution: function not producing grid points because too large block size
    CurvilinearGridRectangular curvilinearGridCreateRectangular(Projection::spherical);
    const auto grid = curvilinearGridCreateRectangular.Compute(angle,
                                                               blockSizeX,
                                                               blockSizeY,
                                                               polygons,
                                                               0);

    // Assert
    auto const numValidNodes = CurvilinearGridCountValidNodes(grid);
    ASSERT_EQ(0, numValidNodes);
}

TEST(CurvilinearGridUniform, MakeCurvilinearInEmptyPolygonSpherical)
{
    // 1 Setup
    const double angle = 0.0;
    const double originX = 0.0;
    const double originY = 89.0;
    const int numColumns = 3;
    const int numRows = 10;
    const double blockSizeX = 0.1; // resolution in meters (when using spherical coordinates distances are usually much larger)
    const double blockSizeY = 0.1;

    // 2 Execution
    CurvilinearGridRectangular const curvilinearGridCreateRectangular(Projection::spherical);
    const auto [nodes, edges, gridIndices] = curvilinearGridCreateRectangular.Compute(numColumns,
                                                                                      numRows,
                                                                                      originX,
                                                                                      originY,
                                                                                      angle,
                                                                                      blockSizeX,
                                                                                      blockSizeY)
                                                 .ConvertCurvilinearToNodesAndEdges();

    Mesh2D mesh(edges, nodes, Projection::spherical);

    // 3 Assert
    ASSERT_EQ(10, mesh.GetNumEdges());
    ASSERT_EQ(8, mesh.GetNumNodes());

    // x coordinate
    ASSERT_EQ(0.0, mesh.m_nodes[0].x);
    ASSERT_EQ(blockSizeX, mesh.m_nodes[1].x);
    ASSERT_EQ(blockSizeX * 2, mesh.m_nodes[2].x);
    ASSERT_EQ(blockSizeX * 3, mesh.m_nodes[3].x);

    ASSERT_EQ(0.0, mesh.m_nodes[4].x);
    ASSERT_EQ(blockSizeX, mesh.m_nodes[5].x);
    ASSERT_EQ(blockSizeX * 2, mesh.m_nodes[6].x);
    ASSERT_EQ(blockSizeX * 3, mesh.m_nodes[7].x);

    // y coordinate
    ASSERT_EQ(89.0, mesh.m_nodes[0].y);
    ASSERT_EQ(89.0, mesh.m_nodes[1].y);
    ASSERT_EQ(89.0, mesh.m_nodes[2].y);
    ASSERT_EQ(89.0, mesh.m_nodes[3].y);

    ASSERT_EQ(90.0, mesh.m_nodes[4].y);
    ASSERT_EQ(90.0, mesh.m_nodes[5].y);
    ASSERT_EQ(90.0, mesh.m_nodes[6].y);
    ASSERT_EQ(90.0, mesh.m_nodes[7].y);
}

TEST(CurvilinearGridUniformCurvilinearGridUniform, InsertFace_OnBottomLeft_ShouldInsertFace)
{
    // Set-up
    auto curvilinearGrid = MakeSmallCurvilinearGrid();

    // Execution
    curvilinearGrid.InsertFace({80009.0, 366937.0});

    // Assert the new coordinates
    constexpr double tolerance = 1e-12;
    ASSERT_NEAR(79913.595823791460, curvilinearGrid.m_gridNodes(0, 0).x, tolerance);
    ASSERT_NEAR(79985.899758176762, curvilinearGrid.m_gridNodes(0, 1).x, tolerance);
    ASSERT_NEAR(-999.00000000000000, curvilinearGrid.m_gridNodes(0, 2).x, tolerance);

    ASSERT_NEAR(367046.61206756550, curvilinearGrid.m_gridNodes(0, 0).y, tolerance);
    ASSERT_NEAR(367110.19327267300, curvilinearGrid.m_gridNodes(0, 1).y, tolerance);
    ASSERT_NEAR(-999.0000000000000, curvilinearGrid.m_gridNodes(0, 2).y, tolerance);
}

TEST(CurvilinearGridUniform, InsertFace_OnBottomRight_ShouldInsertFace)
{
    // Set-up
    auto curvilinearGrid = MakeSmallCurvilinearGrid();

    // Execution
    curvilinearGrid.InsertFace({80166.0, 366544.0});

    // Assert the new coordinates
    constexpr double tolerance = 1e-12;
    ASSERT_NEAR(80176.237835892549, curvilinearGrid.m_gridNodes(5, 0).x, tolerance);
    ASSERT_NEAR(80259.193944680519, curvilinearGrid.m_gridNodes(5, 1).x, tolerance);
    ASSERT_NEAR(-999.0000000000000, curvilinearGrid.m_gridNodes(5, 2).x, tolerance);

    ASSERT_NEAR(366433.98982542212, curvilinearGrid.m_gridNodes(5, 0).y, tolerance);
    ASSERT_NEAR(366433.75796857959, curvilinearGrid.m_gridNodes(5, 1).y, tolerance);
    ASSERT_NEAR(-999.0000000000000, curvilinearGrid.m_gridNodes(5, 2).y, tolerance);
}

TEST(CurvilinearGridUniform, InsertFace_OnTopLeft_ShouldInsertFace)
{
    // Set-up
    auto curvilinearGrid = MakeSmallCurvilinearGrid();

    // Execution
    curvilinearGrid.InsertFace({80612.0, 367407.0});

    // Assert the new coordinates
    constexpr double tolerance = 1e-12;
    ASSERT_NEAR(-999.0000000000000, curvilinearGrid.m_gridNodes(0, 6).x, tolerance);
    ASSERT_NEAR(80385.178875081983, curvilinearGrid.m_gridNodes(0, 7).x, tolerance);
    ASSERT_NEAR(80449.954021552709, curvilinearGrid.m_gridNodes(0, 8).x, tolerance);

    ASSERT_NEAR(-999.0000000000000, curvilinearGrid.m_gridNodes(0, 6).y, tolerance);
    ASSERT_NEAR(367549.55640742677, curvilinearGrid.m_gridNodes(0, 7).y, tolerance);
    ASSERT_NEAR(367626.58734840894, curvilinearGrid.m_gridNodes(0, 8).y, tolerance);
}

TEST(CurvilinearGridUniform, InsertFace_OnTopRight_ShouldInsertFace)
{
    // Set-up
    auto curvilinearGrid = MakeSmallCurvilinearGrid();

    // Execution
    curvilinearGrid.InsertFace({80870.0, 366541.0});

    // Assert the new coordinates
    constexpr double tolerance = 1e-12;
    ASSERT_NEAR(-999.0000000000000, curvilinearGrid.m_gridNodes(5, 6).x, tolerance);
    ASSERT_NEAR(80846.719054016474, curvilinearGrid.m_gridNodes(5, 7).x, tolerance);
    ASSERT_NEAR(80959.766315635425, curvilinearGrid.m_gridNodes(5, 8).x, tolerance);

    ASSERT_NEAR(-999.0000000000000, curvilinearGrid.m_gridNodes(5, 6).y, tolerance);
    ASSERT_NEAR(366346.42989900074, curvilinearGrid.m_gridNodes(5, 7).y, tolerance);
    ASSERT_NEAR(366327.01674911042, curvilinearGrid.m_gridNodes(5, 8).y, tolerance);
}

TEST(CurvilinearGridUniform, InsertFace_OnGridWithHoles_ShouldInsertFace)
{
    // Set-up
    auto curvilinearGrid = MakeSmallCurvilinearGridWithMissingFaces();

    // Execution
    curvilinearGrid.InsertFace({80398.0, 366854.0});

    // Assert the new coordinates
    constexpr double tolerance = 1e-12;
    ASSERT_NEAR(80133.930743945661, curvilinearGrid.m_gridNodes(3, 0).x, tolerance);
    ASSERT_NEAR(80226.579454943669, curvilinearGrid.m_gridNodes(3, 1).x, tolerance);
    ASSERT_NEAR(80320.396380977647, curvilinearGrid.m_gridNodes(3, 2).x, tolerance);
    ASSERT_NEAR(80447.108241179725, curvilinearGrid.m_gridNodes(3, 3).x, tolerance);
    ASSERT_NEAR(80545.663179203286, curvilinearGrid.m_gridNodes(3, 4).x, tolerance);
    ASSERT_NEAR(80623.151266424975, curvilinearGrid.m_gridNodes(3, 5).x, tolerance);
    ASSERT_NEAR(80735.800935924854, curvilinearGrid.m_gridNodes(3, 6).x, tolerance);
    ASSERT_NEAR(80848.959456291268, curvilinearGrid.m_gridNodes(3, 7).x, tolerance);
    ASSERT_NEAR(80962.132385367571, curvilinearGrid.m_gridNodes(3, 8).x, tolerance);

    ASSERT_NEAR(366629.99913221149, curvilinearGrid.m_gridNodes(3, 0).y, tolerance);
    ASSERT_NEAR(366656.00122676318, curvilinearGrid.m_gridNodes(3, 1).y, tolerance);
    ASSERT_NEAR(366679.12590409513, curvilinearGrid.m_gridNodes(3, 2).y, tolerance);
    ASSERT_NEAR(366697.14301043766, curvilinearGrid.m_gridNodes(3, 3).y, tolerance);
    ASSERT_NEAR(366725.32926280121, curvilinearGrid.m_gridNodes(3, 4).y, tolerance);
    ASSERT_NEAR(366718.43748113938, curvilinearGrid.m_gridNodes(3, 5).y, tolerance);
    ASSERT_NEAR(366717.43216295895, curvilinearGrid.m_gridNodes(3, 6).y, tolerance);
    ASSERT_NEAR(366716.74631036510, curvilinearGrid.m_gridNodes(3, 7).y, tolerance);
    ASSERT_NEAR(366718.10475241451, curvilinearGrid.m_gridNodes(3, 8).y, tolerance);
}

TEST(CurvilinearGridUniform, DeleteNode_OnUniformGrid_ShouldDeleteNode)
{
    // Prepare
     auto curvilinearGrid = MakeSmallCurvilinearGrid();

    // Execute
    curvilinearGrid.DeleteNode({80398.0, 366854.0});

    // The number of nodes was 45 now is 44
    auto const numValidNodes = CurvilinearGridCountValidNodes(curvilinearGrid);
    ASSERT_EQ(numValidNodes, 44);
}

void TestDeleteInteriorNodes(meshkernel::CurvilinearGrid curvilinearGrid,
                             const meshkernel::CurvilinearGridNodeIndices first,
                             const meshkernel::CurvilinearGridNodeIndices second)
{

    // Check first, that all nodes are valid
    for (meshkernel::UInt i = 0; i < curvilinearGrid.m_numN; ++i)
    {
        for (meshkernel::UInt j = 0; j < curvilinearGrid.m_numM; ++j)
        {
            EXPECT_TRUE(curvilinearGrid.GetNode(i, j).IsValid());
        }
    }

    meshkernel::UInt lowerLimitI = std::min(first.m_n, second.m_n) + 1;
    meshkernel::UInt upperLimitI = std::max(first.m_n, second.m_n) - 1;

    meshkernel::UInt lowerLimitJ = std::min(first.m_m, second.m_m) + 1;
    meshkernel::UInt upperLimitJ = std::max(first.m_m, second.m_m) - 1;

    meshkernel::UInt expectedInvalidated = (upperLimitI - lowerLimitI + 1) * (upperLimitJ - lowerLimitJ + 1);
    const auto initialSize = static_cast<meshkernel::UInt>(CurvilinearGridCountValidNodes(curvilinearGrid));

    // Delete the nodes interior to a block
    curvilinearGrid.DeleteInterior(first, second);

    auto inRange = [](const meshkernel::UInt v, const meshkernel::UInt l, const meshkernel::UInt u)
    { return l <= v && v <= u; };

    EXPECT_EQ(initialSize - expectedInvalidated, CurvilinearGridCountValidNodes(curvilinearGrid));

    // Check that these nodes have been set to invalid.
    for (meshkernel::UInt i = 0; i < curvilinearGrid.m_numN; ++i)
    {
        for (meshkernel::UInt j = 0; j < curvilinearGrid.m_numM; ++j)
        {
            if (inRange(i, lowerLimitI, upperLimitI) && inRange(j, lowerLimitJ, upperLimitJ))
            {
                EXPECT_FALSE(curvilinearGrid.GetNode(i, j).IsValid()) << "node should be false: " << i << "  " << j;
            }
            else
            {
                EXPECT_TRUE(curvilinearGrid.GetNode(i, j).IsValid()) << "node should be true: " << i << "  " << j;
            }
        }
    }
}

TEST(CurvilinearGridUniform, DeleteInteriorNodesTest)
{
    // Basic, testing of setting nodes inside a box to invalid
    meshkernel::UInt nx = 10;
    meshkernel::UInt ny = 10;
    auto curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);
    TestDeleteInteriorNodes(curvilinearGrid, {1, 1}, {4, 4});

    // Reset the mesh
    curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);
    TestDeleteInteriorNodes(curvilinearGrid, {2, 1}, {5, 4});

    // Reset the mesh
    curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);
    TestDeleteInteriorNodes(curvilinearGrid, {4, 3}, {7, 8});
}

TEST(CurvilinearGridUniform, DeleteInteriorNodesReverseTest)
{
    // testing of setting nodes inside a box to invalid, with lower and upper reversed

    meshkernel::UInt nx = 10;
    meshkernel::UInt ny = 10;

    // Prepare
    auto curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);
    TestDeleteInteriorNodes(curvilinearGrid, {5, 6}, {1, 2});

    // Reset the mesh
    curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);
    TestDeleteInteriorNodes(curvilinearGrid, {5, 6}, {0, 4});
}

TEST(CurvilinearGridUniform, DeleteInteriorNodesMixedTest)
{
    // testing of setting nodes inside a box to invalid, with lower and upper reversed for any of i and j index

    meshkernel::UInt nx = 100;
    meshkernel::UInt ny = 100;

    auto curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);
    TestDeleteInteriorNodes(curvilinearGrid, {5, 1}, {1, 6});

    // Reset grid
    curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);
    TestDeleteInteriorNodes(curvilinearGrid, {1, 6}, {5, 2});
}

TEST(CurvilinearGridUniform, DeleteInteriorNodesFailureTest)
{
    // testing of setting nodes inside a box to invalid with invalid or out of range indices.

    meshkernel::UInt nx = 10;
    meshkernel::UInt ny = 10;

    // Prepare
    auto curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);

    EXPECT_THROW(curvilinearGrid.DeleteInterior({1, meshkernel::constants::missing::uintValue}, {nx, ny}), meshkernel::ConstraintError);
    EXPECT_THROW(curvilinearGrid.DeleteInterior({1, 1}, {meshkernel::constants::missing::uintValue, ny}), meshkernel::ConstraintError);
    EXPECT_THROW(curvilinearGrid.DeleteInterior({1, 1}, {nx, ny}), meshkernel::ConstraintError);
    EXPECT_THROW(curvilinearGrid.DeleteInterior({nx, 1}, {4, 4}), meshkernel::ConstraintError);
}

void TestDeleteExteriorNodes(meshkernel::CurvilinearGrid& curvilinearGrid,
                             const meshkernel::CurvilinearGridNodeIndices first,
                             const meshkernel::CurvilinearGridNodeIndices second)
{
    // Check first, that all nodes are valid
    for (meshkernel::UInt i = 0; i < curvilinearGrid.m_numN; ++i)
    {
        for (meshkernel::UInt j = 0; j < curvilinearGrid.m_numM; ++j)
        {
            EXPECT_TRUE(curvilinearGrid.GetNode(i, j).IsValid());
        }
    }

    meshkernel::UInt lowerLimitI = std::min(first.m_n, second.m_n);
    meshkernel::UInt upperLimitI = std::max(first.m_n, second.m_n);

    meshkernel::UInt lowerLimitJ = std::min(first.m_m, second.m_m);
    meshkernel::UInt upperLimitJ = std::max(first.m_m, second.m_m);

    meshkernel::UInt expectedValid = (upperLimitI - lowerLimitI + 1) * (upperLimitJ - lowerLimitJ + 1);

    // Delete the nodes outside of a block
    curvilinearGrid.DeleteExterior(first, second);

    auto inRange = [](const meshkernel::UInt v, const meshkernel::UInt l, const meshkernel::UInt u)
    { return l <= v && v <= u; };

    EXPECT_EQ(expectedValid, CurvilinearGridCountValidNodes(curvilinearGrid));

    // Check that these exterior nodes have been set to invalid.
    for (meshkernel::UInt i = 0; i < curvilinearGrid.m_numN; ++i)
    {
        for (meshkernel::UInt j = 0; j < curvilinearGrid.m_numM; ++j)
        {
            if (inRange(i, lowerLimitI, upperLimitI) && inRange(j, lowerLimitJ, upperLimitJ))
            {
                EXPECT_TRUE(curvilinearGrid.GetNode(i, j).IsValid()) << "node should be true: " << i << "  " << j;
            }
            else
            {
                EXPECT_FALSE(curvilinearGrid.GetNode(i, j).IsValid()) << "node should be false: " << i << "  " << j;
            }
        }
    }
}

TEST(CurvilinearGridUniform, DeleteExteriorNodesTest)
{
    // Basic, testing of setting nodes inside a box to invalid
    meshkernel::UInt nx = 10;
    meshkernel::UInt ny = 10;
    auto curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);
    TestDeleteExteriorNodes(curvilinearGrid, {1, 1}, {4, 4});

    // Reset the mesh
    curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);
    TestDeleteExteriorNodes(curvilinearGrid, {2, 1}, {5, 4});

    // Reset the mesh
    curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);
    TestDeleteExteriorNodes(curvilinearGrid, {4, 3}, {7, 8});
}

TEST(CurvilinearGridUniform, DeleteExteriorNodesReverseTest)
{
    // testing of setting nodes inside a box to invalid, with lower and upper reversed

    meshkernel::UInt nx = 10;
    meshkernel::UInt ny = 10;

    // Prepare
    auto curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);
    TestDeleteExteriorNodes(curvilinearGrid, {5, 6}, {1, 2});

    // Reset the mesh
    curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);
    TestDeleteExteriorNodes(curvilinearGrid, {5, 6}, {0, 4});
}

TEST(CurvilinearGridUniform, DeleteExteriorNodesMixedTest)
{
    // testing of setting nodes inside a box to invalid, with lower and upper reversed for any of i and j index

    meshkernel::UInt nx = 100;
    meshkernel::UInt ny = 100;

    auto curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);
    TestDeleteExteriorNodes(curvilinearGrid, {5, 1}, {1, 6});

    // Reset grid
    curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);
    TestDeleteExteriorNodes(curvilinearGrid, {1, 6}, {5, 2});
}

TEST(CurvilinearGridUniform, DeleteExteriorNodesFailureTest)
{
    // testing of setting nodes inside a box to invalid with invalid or out of range indices.

    meshkernel::UInt nx = 10;
    meshkernel::UInt ny = 10;

    // Prepare
    auto curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);

    EXPECT_THROW(curvilinearGrid.DeleteExterior({1, meshkernel::constants::missing::uintValue}, {nx, ny}), meshkernel::ConstraintError);
    EXPECT_THROW(curvilinearGrid.DeleteExterior({1, 1}, {meshkernel::constants::missing::uintValue, ny}), meshkernel::ConstraintError);
    EXPECT_THROW(curvilinearGrid.DeleteExterior({1, 1}, {nx, ny}), meshkernel::ConstraintError);
    EXPECT_THROW(curvilinearGrid.DeleteExterior({nx, 1}, {4, 4}), meshkernel::ConstraintError);
}
