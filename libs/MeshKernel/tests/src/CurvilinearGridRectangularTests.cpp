#include <gtest/gtest.h>

#include "MeshKernel/CurvilinearGrid/CurvilinearGridDeleteExterior.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridDeleteInterior.hpp"
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
    auto const numValidNodes = CurvilinearGridCountValidNodes(*curvilinearGrid);
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
    auto const numValidNodes = CurvilinearGridCountValidNodes(*grid);
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
                                                 ->ConvertCurvilinearToNodesAndEdges();

    Mesh2D mesh(edges, nodes, Projection::spherical);

    // 3 Assert
    ASSERT_EQ(10, mesh.GetNumEdges());
    ASSERT_EQ(8, mesh.GetNumNodes());

    // x coordinate
    ASSERT_EQ(0.0, mesh.Node(0).x);
    ASSERT_EQ(blockSizeX, mesh.Node(1).x);
    ASSERT_EQ(blockSizeX * 2, mesh.Node(2).x);
    ASSERT_EQ(blockSizeX * 3, mesh.Node(3).x);

    ASSERT_EQ(0.0, mesh.Node(4).x);
    ASSERT_EQ(blockSizeX, mesh.Node(5).x);
    ASSERT_EQ(blockSizeX * 2, mesh.Node(6).x);
    ASSERT_EQ(blockSizeX * 3, mesh.Node(7).x);

    // y coordinate
    ASSERT_EQ(89.0, mesh.Node(0).y);
    ASSERT_EQ(89.0, mesh.Node(1).y);
    ASSERT_EQ(89.0, mesh.Node(2).y);
    ASSERT_EQ(89.0, mesh.Node(3).y);

    ASSERT_EQ(90.0, mesh.Node(4).y);
    ASSERT_EQ(90.0, mesh.Node(5).y);
    ASSERT_EQ(90.0, mesh.Node(6).y);
    ASSERT_EQ(90.0, mesh.Node(7).y);
}

TEST(CurvilinearGridUniformCurvilinearGridUniform, InsertFace_OnBottomLeft_ShouldInsertFace)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();

    // Execution
    curvilinearGrid->InsertFace({80009.0, 366937.0});

    // Assert the new coordinates
    constexpr double tolerance = 1e-12;
    ASSERT_NEAR(79913.595823791460, curvilinearGrid->GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(79985.899758176762, curvilinearGrid->GetNode(0, 1).x, tolerance);
    ASSERT_NEAR(-999.00000000000000, curvilinearGrid->GetNode(0, 2).x, tolerance);

    ASSERT_NEAR(367046.61206756550, curvilinearGrid->GetNode(0, 0).y, tolerance);
    ASSERT_NEAR(367110.19327267300, curvilinearGrid->GetNode(0, 1).y, tolerance);
    ASSERT_NEAR(-999.0000000000000, curvilinearGrid->GetNode(0, 2).y, tolerance);
}

TEST(CurvilinearGridUniform, InsertFace_OnBottomRight_ShouldInsertFace)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();

    // Execution
    curvilinearGrid->InsertFace({80166.0, 366544.0});

    // Assert the new coordinates
    constexpr double tolerance = 1e-12;
    ASSERT_NEAR(80176.237835892549, curvilinearGrid->GetNode(5, 0).x, tolerance);
    ASSERT_NEAR(80259.193944680519, curvilinearGrid->GetNode(5, 1).x, tolerance);
    ASSERT_NEAR(-999.0000000000000, curvilinearGrid->GetNode(5, 2).x, tolerance);

    ASSERT_NEAR(366433.98982542212, curvilinearGrid->GetNode(5, 0).y, tolerance);
    ASSERT_NEAR(366433.75796857959, curvilinearGrid->GetNode(5, 1).y, tolerance);
    ASSERT_NEAR(-999.0000000000000, curvilinearGrid->GetNode(5, 2).y, tolerance);
}

TEST(CurvilinearGridUniform, InsertFace_OnTopLeft_ShouldInsertFace)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();

    // Execution
    curvilinearGrid->InsertFace({80612.0, 367407.0});

    // Assert the new coordinates
    constexpr double tolerance = 1e-12;
    ASSERT_NEAR(-999.0000000000000, curvilinearGrid->GetNode(0, 6).x, tolerance);
    ASSERT_NEAR(80385.178875081983, curvilinearGrid->GetNode(0, 7).x, tolerance);
    ASSERT_NEAR(80449.954021552709, curvilinearGrid->GetNode(0, 8).x, tolerance);

    ASSERT_NEAR(-999.0000000000000, curvilinearGrid->GetNode(0, 6).y, tolerance);
    ASSERT_NEAR(367549.55640742677, curvilinearGrid->GetNode(0, 7).y, tolerance);
    ASSERT_NEAR(367626.58734840894, curvilinearGrid->GetNode(0, 8).y, tolerance);
}

TEST(CurvilinearGridUniform, InsertFace_OnTopRight_ShouldInsertFace)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();

    // Execution
    curvilinearGrid->InsertFace({80870.0, 366541.0});

    // Assert the new coordinates
    constexpr double tolerance = 1e-12;
    ASSERT_NEAR(-999.0000000000000, curvilinearGrid->GetNode(5, 6).x, tolerance);
    ASSERT_NEAR(80846.719054016474, curvilinearGrid->GetNode(5, 7).x, tolerance);
    ASSERT_NEAR(80959.766315635425, curvilinearGrid->GetNode(5, 8).x, tolerance);

    ASSERT_NEAR(-999.0000000000000, curvilinearGrid->GetNode(5, 6).y, tolerance);
    ASSERT_NEAR(366346.42989900074, curvilinearGrid->GetNode(5, 7).y, tolerance);
    ASSERT_NEAR(366327.01674911042, curvilinearGrid->GetNode(5, 8).y, tolerance);
}

TEST(CurvilinearGridUniform, InsertFace_OnGridWithHoles_ShouldInsertFace)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGridWithMissingFaces();

    // Execution
    curvilinearGrid->InsertFace({80398.0, 366854.0});

    // Assert the new coordinates
    constexpr double tolerance = 1e-12;
    ASSERT_NEAR(80133.930743945661, curvilinearGrid->GetNode(3, 0).x, tolerance);
    ASSERT_NEAR(80226.579454943669, curvilinearGrid->GetNode(3, 1).x, tolerance);
    ASSERT_NEAR(80320.396380977647, curvilinearGrid->GetNode(3, 2).x, tolerance);
    ASSERT_NEAR(80447.108241179725, curvilinearGrid->GetNode(3, 3).x, tolerance);
    ASSERT_NEAR(80545.663179203286, curvilinearGrid->GetNode(3, 4).x, tolerance);
    ASSERT_NEAR(80623.151266424975, curvilinearGrid->GetNode(3, 5).x, tolerance);
    ASSERT_NEAR(80735.800935924854, curvilinearGrid->GetNode(3, 6).x, tolerance);
    ASSERT_NEAR(80848.959456291268, curvilinearGrid->GetNode(3, 7).x, tolerance);
    ASSERT_NEAR(80962.132385367571, curvilinearGrid->GetNode(3, 8).x, tolerance);

    ASSERT_NEAR(366629.99913221149, curvilinearGrid->GetNode(3, 0).y, tolerance);
    ASSERT_NEAR(366656.00122676318, curvilinearGrid->GetNode(3, 1).y, tolerance);
    ASSERT_NEAR(366679.12590409513, curvilinearGrid->GetNode(3, 2).y, tolerance);
    ASSERT_NEAR(366697.14301043766, curvilinearGrid->GetNode(3, 3).y, tolerance);
    ASSERT_NEAR(366725.32926280121, curvilinearGrid->GetNode(3, 4).y, tolerance);
    ASSERT_NEAR(366718.43748113938, curvilinearGrid->GetNode(3, 5).y, tolerance);
    ASSERT_NEAR(366717.43216295895, curvilinearGrid->GetNode(3, 6).y, tolerance);
    ASSERT_NEAR(366716.74631036510, curvilinearGrid->GetNode(3, 7).y, tolerance);
    ASSERT_NEAR(366718.10475241451, curvilinearGrid->GetNode(3, 8).y, tolerance);
}

TEST(CurvilinearGridUniform, DeleteNode_OnUniformGrid_ShouldDeleteNode)
{
    // Prepare
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();

    // Execute
    curvilinearGrid->DeleteNode({80398.0, 366854.0});

    // The number of nodes was 45 now is 44
    auto const numValidNodes = CurvilinearGridCountValidNodes(*curvilinearGrid);
    ASSERT_EQ(numValidNodes, 44);
}

void TestDeleteInteriorNodes(meshkernel::CurvilinearGrid& curvilinearGrid,
                             const meshkernel::CurvilinearGridNodeIndices first,
                             const meshkernel::CurvilinearGridNodeIndices second)
{

    // Check first, that all nodes are valid
    for (meshkernel::UInt i = 0; i < curvilinearGrid.NumN(); ++i)
    {
        for (meshkernel::UInt j = 0; j < curvilinearGrid.NumM(); ++j)
        {
            EXPECT_TRUE(curvilinearGrid.GetNode(i, j).IsValid());
        }
    }

    meshkernel::UInt lowerLimitI = std::min(first.m_n, second.m_n) + 1;
    meshkernel::UInt upperLimitI = std::max(first.m_n, second.m_n) - 1;

    meshkernel::UInt lowerLimitJ = std::min(first.m_m, second.m_m) + 1;
    meshkernel::UInt upperLimitJ = std::max(first.m_m, second.m_m) - 1;

    meshkernel::UInt expectedInvalidated = (upperLimitI - lowerLimitI + 1) * (upperLimitJ - lowerLimitJ + 1);
    const auto initialSize = static_cast<UInt>(CurvilinearGridCountValidNodes(curvilinearGrid));
    CurvilinearGridDeleteInterior curvilinearGridDeleteInterior(curvilinearGrid);
    curvilinearGridDeleteInterior.m_lowerLeft = {lowerLimitI - 1, lowerLimitJ - 1};
    curvilinearGridDeleteInterior.m_upperRight = {upperLimitI + 1, upperLimitJ + 1};

    // Delete the nodes interior to a block
    curvilinearGridDeleteInterior.Compute();

    auto inRange = [](const meshkernel::UInt v, const meshkernel::UInt l, const meshkernel::UInt u)
    { return l <= v && v <= u; };

    EXPECT_EQ(initialSize - expectedInvalidated, CurvilinearGridCountValidNodes(curvilinearGrid));

    // Check that these nodes have been set to invalid.
    for (meshkernel::UInt i = 0; i < curvilinearGrid.NumN(); ++i)
    {
        for (meshkernel::UInt j = 0; j < curvilinearGrid.NumM(); ++j)
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
    TestDeleteInteriorNodes(*curvilinearGrid, {1, 1}, {4, 4});

    // Reset the mesh
    curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);
    TestDeleteInteriorNodes(*curvilinearGrid, {2, 1}, {5, 4});

    // Reset the mesh
    curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);
    TestDeleteInteriorNodes(*curvilinearGrid, {4, 3}, {7, 8});
}

TEST(CurvilinearGridUniform, DeleteInteriorNodesReverseTest)
{
    // testing of setting nodes inside a box to invalid, with lower and upper reversed

    meshkernel::UInt nx = 10;
    meshkernel::UInt ny = 10;

    // Prepare
    auto curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);
    TestDeleteInteriorNodes(*curvilinearGrid, {5, 6}, {1, 2});

    // Reset the mesh
    curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);
    TestDeleteInteriorNodes(*curvilinearGrid, {5, 6}, {0, 4});
}

TEST(CurvilinearGridUniform, DeleteInteriorNodesMixedTest)
{
    // testing of setting nodes inside a box to invalid, with lower and upper reversed for any of i and j index

    meshkernel::UInt nx = 100;
    meshkernel::UInt ny = 100;

    auto curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);
    TestDeleteInteriorNodes(*curvilinearGrid, {1, 5}, {6, 1});

    // Reset grid
    curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);
    TestDeleteInteriorNodes(*curvilinearGrid, {6, 1}, {2, 5});
}

TEST(CurvilinearGridUniform, DeleteInteriorNodesFailureTest)
{
    // testing of setting nodes inside a box to invalid with invalid or out of range indices.

    meshkernel::UInt nx = 10;
    meshkernel::UInt ny = 10;

    // Prepare
    auto curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);
    CurvilinearGridDeleteInterior curvilinearGridDeleteInterior(*curvilinearGrid);

    EXPECT_THROW(curvilinearGrid->ComputeBlockFromCornerPoints(CurvilinearGridNodeIndices{1, meshkernel::constants::missing::uintValue}, CurvilinearGridNodeIndices{nx, ny}), meshkernel::ConstraintError);

    EXPECT_THROW(curvilinearGrid->ComputeBlockFromCornerPoints(CurvilinearGridNodeIndices{1, 1}, CurvilinearGridNodeIndices{meshkernel::constants::missing::uintValue, ny}), meshkernel::ConstraintError);

    EXPECT_THROW(curvilinearGrid->ComputeBlockFromCornerPoints(CurvilinearGridNodeIndices{1, 1}, CurvilinearGridNodeIndices{nx, ny}), meshkernel::ConstraintError);

    EXPECT_THROW(curvilinearGrid->ComputeBlockFromCornerPoints(CurvilinearGridNodeIndices{10, 1}, CurvilinearGridNodeIndices{4, 4}), meshkernel::ConstraintError);
}

void TestDeleteExteriorNodes(meshkernel::CurvilinearGrid& curvilinearGrid,
                             const meshkernel::CurvilinearGridNodeIndices first,
                             const meshkernel::CurvilinearGridNodeIndices second)
{
    // Check first, that all nodes are valid
    for (meshkernel::UInt i = 0; i < curvilinearGrid.NumN(); ++i)
    {
        for (meshkernel::UInt j = 0; j < curvilinearGrid.NumM(); ++j)
        {
            EXPECT_TRUE(curvilinearGrid.GetNode(i, j).IsValid());
        }
    }

    meshkernel::UInt lowerLimitI = std::min(first.m_n, second.m_n);
    meshkernel::UInt upperLimitI = std::max(first.m_n, second.m_n);

    meshkernel::UInt lowerLimitJ = std::min(first.m_m, second.m_m);
    meshkernel::UInt upperLimitJ = std::max(first.m_m, second.m_m);

    meshkernel::UInt expectedValid = (upperLimitI - lowerLimitI + 1) * (upperLimitJ - lowerLimitJ + 1);

    CurvilinearGridDeleteExterior curvilinearGridDeleteExterior(curvilinearGrid);
    curvilinearGridDeleteExterior.m_lowerLeft = {lowerLimitI, lowerLimitJ};
    curvilinearGridDeleteExterior.m_upperRight = {upperLimitI, upperLimitJ};

    // Delete the nodes outside of a block
    curvilinearGridDeleteExterior.Compute();

    auto inRange = [](const meshkernel::UInt v, const meshkernel::UInt l, const meshkernel::UInt u)
    { return l <= v && v <= u; };

    EXPECT_EQ(expectedValid, CurvilinearGridCountValidNodes(curvilinearGrid));

    // Check that these exterior nodes have been set to invalid.
    for (meshkernel::UInt i = 0; i < curvilinearGrid.NumN(); ++i)
    {
        for (meshkernel::UInt j = 0; j < curvilinearGrid.NumM(); ++j)
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
    TestDeleteExteriorNodes(*curvilinearGrid, {1, 1}, {4, 4});

    // Reset the mesh
    curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);
    TestDeleteExteriorNodes(*curvilinearGrid, {1, 2}, {4, 5});

    // Reset the mesh
    curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);
    TestDeleteExteriorNodes(*curvilinearGrid, {3, 4}, {8, 7});
}

TEST(CurvilinearGridUniform, DeleteExteriorNodesReverseTest)
{
    // testing of setting nodes inside a box to invalid, with lower and upper reversed

    meshkernel::UInt nx = 10;
    meshkernel::UInt ny = 10;

    // Prepare
    auto curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);
    TestDeleteExteriorNodes(*curvilinearGrid, {6, 5}, {2, 1});

    // Reset the mesh
    curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);
    TestDeleteExteriorNodes(*curvilinearGrid, {6, 5}, {4, 0});
}

TEST(CurvilinearGridUniform, DeleteExteriorNodesMixedTest)
{
    // testing of setting nodes inside a box to invalid, with lower and upper reversed for any of i and j index

    meshkernel::UInt nx = 100;
    meshkernel::UInt ny = 100;

    auto curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);
    TestDeleteExteriorNodes(*curvilinearGrid, {1, 5}, {6, 1});

    // Reset grid
    curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);
    TestDeleteExteriorNodes(*curvilinearGrid, {6, 1}, {2, 5});
}

TEST(CurvilinearGridUniform, DeleteExteriorNodesFailureTest)
{
    // testing of setting nodes inside a box to invalid with invalid or out of range indices.

    UInt nx = 10;
    UInt ny = 10;

    // Prepare
    const auto curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, nx, ny);
    CurvilinearGridDeleteExterior curvilinearGridDeleteExterior(*curvilinearGrid);

    EXPECT_THROW(curvilinearGrid->ComputeBlockFromCornerPoints(CurvilinearGridNodeIndices{1, meshkernel::constants::missing::uintValue}, CurvilinearGridNodeIndices{nx, ny}), meshkernel::ConstraintError);

    EXPECT_THROW(curvilinearGrid->ComputeBlockFromCornerPoints(CurvilinearGridNodeIndices{1, 1}, CurvilinearGridNodeIndices{meshkernel::constants::missing::uintValue, ny}), meshkernel::ConstraintError);

    EXPECT_THROW(curvilinearGrid->ComputeBlockFromCornerPoints(CurvilinearGridNodeIndices{1, 1}, CurvilinearGridNodeIndices{nx, ny}), meshkernel::ConstraintError);

    EXPECT_THROW(curvilinearGrid->ComputeBlockFromCornerPoints(CurvilinearGridNodeIndices{nx, 1}, CurvilinearGridNodeIndices{4, 4}), meshkernel::ConstraintError);
}

TEST(CurvilinearGridUniform, NumM_ReturnsNumberOfNodeColumns)
{
    const int numM = 7;
    const int numN = 5;
    const auto subject = MakeCurvilinearGrid(0., 0., 1., 1., numM, numN);
    EXPECT_EQ(7, subject->NumM());
}

TEST(CurvilinearGridUniform, NumM_ReturnsNumberOfRowColumns)
{
    const int numM = 7;
    const int numN = 5;
    const auto subject = MakeCurvilinearGrid(0., 0., 1., 1., numM, numN);
    EXPECT_EQ(5, subject->NumN());
}

TEST(CurvilinearGridUniform, GetNodes_ReturnsCopyOfNodeMatrix)
{
    // NOTE: the storage is transposed, i.e. nodes.rows() corresponds with the number of node columns as indicated by numM and v.v.
    const int numM = 7;
    const int numN = 5;
    const auto subject = MakeCurvilinearGrid(0., 0., 1., 1., numM, numN);
    const auto nodes = subject->GetNodes();
    EXPECT_EQ(numN, nodes.rows());
    EXPECT_EQ(numM, nodes.cols());
}

TEST(CurvilinearGridUniform, GetRowVector_ReturnsVectorOfLengthNumM)
{
    const auto subject = MakeCurvilinearGrid(0., 0., 1., 1., 2, 3);
    EXPECT_EQ(subject->NumN(), subject->GetNodeVectorAtN(1).size());
}

TEST(CurvilinearGridUniform, GetColumnVector_ReturnsVectorOfLengthNumN)
{
    const auto subject = MakeCurvilinearGrid(0., 0., 1., 1., 2, 3);
    EXPECT_EQ(subject->NumN(), subject->GetNodeVectorAtN(1).size());
}

TEST(CurvilinearGridUniform, ConvertCurvilinearToNodesAndEdges_ReturnsSerializedNodes)
{
    // note: this is a characterization test to document the current behavior of MakeCurvilinearGrid,
    // the constructor of CurvilinearGrid taking a matrix of points, and
    // CurvilinearGrid::ConvertCurvilinearToNodesAndEdges
    // A grid of nx=3 and ny=2 returns node coordinates with 3 x-coordinates and 2 y-coordinates!!!
    const auto grid = MakeCurvilinearGrid(2.0, 1.0, 2.0, 1.0, 3, 2);

    EXPECT_EQ(3, grid->NumM());
    EXPECT_EQ(2, grid->NumN());

    const auto [nodes, edges, gridIndices] = grid->ConvertCurvilinearToNodesAndEdges();
    const std::vector<Point> expected_nodes = {{2., 1.}, {4., 1.}, {6., 1.}, {2., 2.}, {4., 2.}, {6., 2.}};

    EXPECT_EQ(expected_nodes.size(), nodes.size());
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        EXPECT_TRUE(meshkernel::IsEqual(expected_nodes[i], nodes[i], 1e-4))
            << "#" << i << ": "
            << "(" << expected_nodes[i].x << "," << expected_nodes[i].y << ") != "
            << "(" << nodes[i].x << "," << nodes[i].y << ")";
    }
}

TEST(CurvilinearGridUniform, ConvertCurvilinearToNodesAndEdges_ReturnsSerializedNodesWithInvalidNode)
{
    // note: this is a characterization test to document the current behavior of MakeCurvilinearGrid,
    // the constructor of CurvilinearGrid taking a matrix of points, and
    // CurvilinearGrid::ConvertCurvilinearToNodesAndEdges
    // A grid of nx=3 and ny=2 returns node coordinates with 3 y-coordinates and 2 x-coordinates!!!
    const auto grid = MakeCurvilinearGrid(2.0, 1.0, 2.0, 1.0, 3, 2);

    EXPECT_EQ(3, grid->NumM());
    EXPECT_EQ(2, grid->NumN());

    grid->GetNode(1, 2).SetInvalid();

    const auto [nodes, edges, gridIndices] = grid->ConvertCurvilinearToNodesAndEdges();
    const std::vector<Point> expected_nodes = {{2., 1.}, {4., 1.}, {6., 1.}, {2., 2.}, {4., 2.}, {-999., -999.}};

    EXPECT_EQ(expected_nodes.size(), nodes.size());
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        EXPECT_TRUE(meshkernel::IsEqual(expected_nodes[i], nodes[i], 1e-4))
            << "#" << i << ": "
            << "(" << expected_nodes[i].x << "," << expected_nodes[i].y << ") != "
            << "(" << nodes[i].x << "," << nodes[i].y << ")";
    }
}

TEST(CurvilinearGridUniform, ConvertCurvilinearToNodesAndEdges_ReturnsSerializedEdges)
{
    // note: this is a characterization test to document the current behavior of MakeCurvilinearGrid,
    // the constructor of CurvilinearGrid taking a matrix of points, and
    // CurvilinearGrid::ConvertCurvilinearToNodesAndEdges
    // A grid of nx=3 and ny=2 returns node coordinates with 3 x-coordinates and 2 y-coordinates
    const auto grid = MakeCurvilinearGrid(2.0, 1.0, 2.0, 1.0, 3, 2);

    EXPECT_EQ(3, grid->NumM());
    EXPECT_EQ(2, grid->NumN());

    const auto [nodes, edges, gridIndices] = grid->ConvertCurvilinearToNodesAndEdges();
    const std::vector<Edge> expected_edges = {{{0u, 3u}, {1u, 4u}, {2u, 5u}, {0u, 1u}, {1u, 2u}, {3u, 4u}, {4u, 5u}}};

    EXPECT_EQ(expected_edges.size(), edges.size());
    for (size_t i = 0; i < edges.size(); ++i)
    {
        EXPECT_EQ(expected_edges[i], edges[i])
            << "#" << i << ": "
            << "[" << expected_edges[i].first << "," << expected_edges[i].second << "] != "
            << "[" << edges[i].first << "," << edges[i].second << "]";
    }
}

TEST(CurvilinearGridUniform, ConvertCurvilinearToNodesAndEdges_ReturnsSerializedEdgesAlsoForEdgesWithInvalidNode)
{
    // note: this is a characterization test to document the current behavior of MakeCurvilinearGrid,
    // the constructor of CurvilinearGrid taking a matrix of points, and
    // CurvilinearGrid::ConvertCurvilinearToNodesAndEdges
    // A grid of nx=3 and ny=2 returns node coordinates with 3 x-coordinates and 2 y-coordinates
    const auto grid = MakeCurvilinearGrid(2.0, 1.0, 2.0, 1.0, 3, 2);

    EXPECT_EQ(3, grid->NumM());
    EXPECT_EQ(2, grid->NumN());

    grid->GetNode(1, 0).SetInvalid();

    const auto [nodes, edges, gridIndices] = grid->ConvertCurvilinearToNodesAndEdges();
    const std::vector<Edge> expected_edges = {{{0u, 3u}, {1u, 4u}, {2u, 5u}, {0u, 1u}, {1u, 2u}, {3u, 4u}, {4u, 5u}}};

    EXPECT_EQ(expected_edges.size(), edges.size());
    for (size_t i = 0; i < edges.size(); ++i)
    {
        EXPECT_EQ(expected_edges[i], edges[i])
            << "#" << i << ": "
            << "[" << expected_edges[i].first << "," << expected_edges[i].second << "] != "
            << "[" << edges[i].first << "," << edges[i].second << "]";
    }
}
