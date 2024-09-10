//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2024.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation version 3.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// contact: delft3d.support@deltares.nl
// Stichting Deltares
// P.O. Box 177
// 2600 MH Delft, The Netherlands
//
// All indications and logos of, and references to, "Delft3D" and "Deltares"
// are registered trademarks of Stichting Deltares, and remain the property of
// Stichting Deltares. All rights reserved.
//
//------------------------------------------------------------------------------

#include <numbers>
#include <random>
#include <utility>

#include <gtest/gtest.h>

#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridCurvature.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridDeRefinement.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridDeleteExterior.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridDeleteInterior.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridLineAttractionRepulsion.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridLineMirror.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridLineShift.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridOrthogonalization.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridRefinement.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridSmoothing.hpp"
#include "MeshKernel/Operations.hpp"
#include "MeshKernel/UndoActions/UndoActionStack.hpp"
#include "MeshKernel/Utilities/LinearAlgebra.hpp"
#include "TestUtils/MakeCurvilinearGrids.hpp"

namespace mk = meshkernel;

TEST(CurvilinearBasicTests, SimpleAddGridLineAtBoundary)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 2.0;

    constexpr size_t nx = 27;
    constexpr size_t ny = 13;

    std::unique_ptr<mk::CurvilinearGrid> grid = MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);
    grid->ComputeGridNodeTypes();

    // Nodes in the original mesh
    const std::vector<mk::Point> originalPoints = grid->ComputeNodes();
    bool addedLine = false;
    std::unique_ptr<mk::UndoAction> undoAction;

    EXPECT_EQ(grid->NumN(), ny);
    EXPECT_EQ(grid->NumM(), nx);

    // Add grid line to bottom of domain
    std::tie(addedLine, undoAction) = grid->AddGridLineAtBoundary({0, 0}, {0, 1});
    EXPECT_EQ(grid->NumN(), ny + 1);
    EXPECT_EQ(grid->NumM(), nx);

    // There should be no start/end offset changes
    EXPECT_EQ(grid->StartOffset().m_n, 0);
    EXPECT_EQ(grid->StartOffset().m_m, 0);
    EXPECT_EQ(grid->EndOffset().m_n, 0);
    EXPECT_EQ(grid->EndOffset().m_m, 0);

    undoAction->Restore();

    EXPECT_EQ(grid->NumN(), ny);
    EXPECT_EQ(grid->NumM(), nx);

    // There should be no start/end offset changes
    EXPECT_EQ(grid->StartOffset().m_n, 1);
    EXPECT_EQ(grid->StartOffset().m_m, 0);
    EXPECT_EQ(grid->EndOffset().m_n, 0);
    EXPECT_EQ(grid->EndOffset().m_m, 0);
}

TEST(CurvilinearBasicTests, SimpleDeleteInterior)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 1.0;

    constexpr size_t nx = 10;
    constexpr size_t ny = 10;

    std::unique_ptr<mk::CurvilinearGrid> grid = MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);
    grid->ComputeGridNodeTypes();

    mk::CurvilinearGridDeleteInterior deleteInterior(*grid);
    deleteInterior.m_lowerLeft = {3, 2};
    deleteInterior.m_upperRight = {6, 7};
    auto deleteAction = deleteInterior.Compute();

    auto validityCheck = [](const mk::CurvilinearGrid& mesh, const mk::CurvilinearGridDeleteInterior& interior, const bool checkDeletedRegion)
    {
        for (mk::UInt r = 0; r < ny; ++r)
        {
            bool invalidRow = checkDeletedRegion && interior.m_lowerLeft.m_n < r && r < interior.m_upperRight.m_n;

            for (mk::UInt c = 0; c < nx; ++c)
            {
                bool invalidCol = checkDeletedRegion && interior.m_lowerLeft.m_m < c && c < interior.m_upperRight.m_m;

                if (invalidRow && invalidCol)
                {
                    EXPECT_FALSE(mesh.GetNode(r, c).IsValid()) << "pos " << r << "  " << c;
                }
                else
                {
                    EXPECT_TRUE(mesh.GetNode(r, c).IsValid()) << "pos " << r << "  " << c;
                }
            }
        }
    };

    validityCheck(*grid, deleteInterior, true /* checkDeletedRegion */);

    //--------------------------------
    // Now undo the deletion of the block

    deleteAction->Restore();
    validityCheck(*grid, deleteInterior, false /* checkDeletedRegion */);

    //--------------------------------
    // Redo the interior bock delete

    deleteAction->Commit();
    validityCheck(*grid, deleteInterior, true /* checkDeletedRegion */);
}

TEST(CurvilinearBasicTests, SimpleDeleteExterior)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 1.0;

    constexpr size_t nx = 50;
    constexpr size_t ny = 50;

    mk::UndoActionStack undoActions;

    std::unique_ptr<mk::CurvilinearGrid> grid = MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);

    auto validityCheck = [](const mk::CurvilinearGrid& mesh, const mk::CurvilinearGridDeleteExterior& exterior, const bool checkDeletedRegion)
    {
        for (mk::UInt r = 0; r < ny; ++r)
        {
            bool invalidRow = checkDeletedRegion && (r < exterior.m_lowerLeft.m_n || exterior.m_upperRight.m_n < r);

            for (mk::UInt c = 0; c < nx; ++c)
            {
                bool invalidCol = checkDeletedRegion && (c < exterior.m_lowerLeft.m_m || exterior.m_upperRight.m_m < c);

                if (invalidRow || invalidCol)
                {
                    EXPECT_FALSE(mesh.GetNode(r, c).IsValid()) << "pos " << r << "  " << c;
                }
                else
                {
                    EXPECT_TRUE(mesh.GetNode(r, c).IsValid()) << "pos " << r << "  " << c;
                }
            }
        }
    };

    //--------------------------------

    mk::CurvilinearGridDeleteExterior deleteExterior(*grid);

    deleteExterior.m_lowerLeft = {20, 15};
    deleteExterior.m_upperRight = {30, 30};

    //--------------------------------

    auto deleteAction = deleteExterior.Compute();

    validityCheck(*grid, deleteExterior, true /* checkDeletedRegion */);

    //--------------------------------
    // undo the deletion of the exterior

    deleteAction->Restore();
    validityCheck(*grid, deleteExterior, false /* checkDeletedRegion */);

    //--------------------------------
    // Redo the exterior delete

    deleteAction->Commit();
    validityCheck(*grid, deleteExterior, true /* checkDeletedRegion */);
}

TEST(CurvilinearBasicTests, AddGridLinesAllRound)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 1.0;

    constexpr size_t nx = 10;
    constexpr size_t ny = 10;

    mk::UndoActionStack undoActions;

    std::unique_ptr<mk::CurvilinearGrid> grid = MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);
    // Nodes in the original mesh
    const std::vector<mk::Point> originalPoints = grid->ComputeNodes();
    grid->ComputeGridNodeTypes();

    //--------------------------------

    mk::CurvilinearGridLineMirror lineMirror(*grid, deltaX);

    EXPECT_EQ(grid->NumN(), ny);
    EXPECT_EQ(grid->NumM(), nx);

    // Bottom
    lineMirror.m_lines.push_back(mk::CurvilinearGridLine({0, 0}, {0, nx - 1}));
    undoActions.Add(lineMirror.Compute());

    EXPECT_EQ(grid->NumN(), ny + 1);
    EXPECT_EQ(grid->NumM(), nx);

    // There should be no start/end offset changes
    EXPECT_EQ(grid->StartOffset().m_n, 0);
    EXPECT_EQ(grid->StartOffset().m_m, 0);
    EXPECT_EQ(grid->EndOffset().m_n, 0);
    EXPECT_EQ(grid->EndOffset().m_m, 0);

    //--------------------------------

    // Left
    lineMirror.m_lines[0] = mk::CurvilinearGridLine({0, 0}, {ny, 0});
    undoActions.Add(lineMirror.Compute());

    EXPECT_EQ(grid->NumN(), ny + 1);
    EXPECT_EQ(grid->NumM(), nx + 1);

    // There should be no start/end offset changes
    EXPECT_EQ(grid->StartOffset().m_n, 0);
    EXPECT_EQ(grid->StartOffset().m_m, 0);
    EXPECT_EQ(grid->EndOffset().m_n, 0);
    EXPECT_EQ(grid->EndOffset().m_m, 0);

    //--------------------------------

    // Top
    lineMirror.m_lines[0] = mk::CurvilinearGridLine({ny, 0}, {ny, nx});
    undoActions.Add(lineMirror.Compute());

    EXPECT_EQ(grid->NumN(), ny + 2);
    EXPECT_EQ(grid->NumM(), nx + 1);

    // There should be no start/end offset changes
    EXPECT_EQ(grid->StartOffset().m_n, 0);
    EXPECT_EQ(grid->StartOffset().m_m, 0);
    EXPECT_EQ(grid->EndOffset().m_n, 0);
    EXPECT_EQ(grid->EndOffset().m_m, 0);

    //--------------------------------

    // Right
    lineMirror.m_lines[0] = mk::CurvilinearGridLine({0, nx}, {ny + 1, nx});
    undoActions.Add(lineMirror.Compute());

    EXPECT_EQ(grid->NumN(), ny + 2);
    EXPECT_EQ(grid->NumM(), nx + 2);

    // There should be no start/end offset changes
    EXPECT_EQ(grid->StartOffset().m_n, 0);
    EXPECT_EQ(grid->StartOffset().m_m, 0);
    EXPECT_EQ(grid->EndOffset().m_n, 0);
    EXPECT_EQ(grid->EndOffset().m_m, 0);

    //--------------------------------

    constexpr double tolerance = 1.0e-14;

    for (mk::UInt i = 0; i < grid->NumN(); ++i)
    {
        EXPECT_NEAR(grid->GetNode(i, 0).x, -1.0, tolerance);
        EXPECT_NEAR(grid->GetNode(i, 0).y, originY - deltaY + static_cast<double>(i) * deltaY, tolerance);
    }

    for (mk::UInt i = 0; i < grid->NumN(); ++i)
    {
        EXPECT_NEAR(grid->GetNode(i, nx + 1).x, -deltaX + static_cast<double>(nx + 1) * deltaX, tolerance);
        EXPECT_NEAR(grid->GetNode(i, nx + 1).y, originY - deltaY + static_cast<double>(i) * deltaY, tolerance);
    }

    for (mk::UInt i = 0; i < grid->NumM(); ++i)
    {
        EXPECT_NEAR(grid->GetNode(0, i).x, originX - deltaX + static_cast<double>(i) * deltaX, tolerance);
        EXPECT_NEAR(grid->GetNode(0, i).y, -deltaY, tolerance);
    }

    for (mk::UInt i = 0; i < grid->NumM(); ++i)
    {
        EXPECT_NEAR(grid->GetNode(ny + 1, i).x, originX - deltaX + static_cast<double>(i) * deltaX, tolerance);
        EXPECT_NEAR(grid->GetNode(ny + 1, i).y, originY + static_cast<double>(ny) * deltaY, tolerance);
    }

    // Now undo the adding of the bondary grid lines checking that the
    // size of the mesh and the start/end offsets are correct.
    //
    undoActions.Undo();

    EXPECT_EQ(grid->NumN(), ny + 2);
    EXPECT_EQ(grid->NumM(), nx + 1);

    // There should be no start/end offset changes
    EXPECT_EQ(grid->StartOffset().m_n, 0);
    EXPECT_EQ(grid->StartOffset().m_m, 0);
    EXPECT_EQ(grid->EndOffset().m_n, 0);
    EXPECT_EQ(grid->EndOffset().m_m, 1);

    //--------------------------------
    undoActions.Undo();

    EXPECT_EQ(grid->NumN(), ny + 1);
    EXPECT_EQ(grid->NumM(), nx + 1);

    // There should be no start/end offset changes
    EXPECT_EQ(grid->StartOffset().m_n, 0);
    EXPECT_EQ(grid->StartOffset().m_m, 0);
    EXPECT_EQ(grid->EndOffset().m_n, 1);
    EXPECT_EQ(grid->EndOffset().m_m, 1);

    //--------------------------------
    undoActions.Undo();

    EXPECT_EQ(grid->NumN(), ny + 1);
    EXPECT_EQ(grid->NumM(), nx);

    // There should be no start/end offset changes
    EXPECT_EQ(grid->StartOffset().m_n, 0);
    EXPECT_EQ(grid->StartOffset().m_m, 1);
    EXPECT_EQ(grid->EndOffset().m_n, 1);
    EXPECT_EQ(grid->EndOffset().m_m, 1);

    //--------------------------------
    undoActions.Undo();

    EXPECT_EQ(grid->NumN(), ny);
    EXPECT_EQ(grid->NumM(), nx);

    // There should be no start/end offset changes
    EXPECT_EQ(grid->StartOffset().m_n, 1);
    EXPECT_EQ(grid->StartOffset().m_m, 1);
    EXPECT_EQ(grid->EndOffset().m_n, 1);
    EXPECT_EQ(grid->EndOffset().m_m, 1);

    //--------------------------------
    // Check that the nodes can be accessed with the correct index.
    // and not the indices without the offset.

    mk::UInt count = 0;

    ASSERT_EQ(originalPoints.size(), grid->NumN() * grid->NumM());

    for (mk::UInt r = 0; r < grid->NumN(); ++r)
    {
        for (mk::UInt c = 0; c < grid->NumM(); ++c)
        {
            EXPECT_EQ(originalPoints[count].x, grid->GetNode(r, c).x);
            EXPECT_EQ(originalPoints[count].y, grid->GetNode(r, c).y);
            ++count;
        }
    }
}

TEST(CurvilinearBasicTests, GridSmoothing)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 1.0;

    constexpr size_t nx = 30;
    constexpr size_t ny = 30;

    std::unique_ptr<mk::CurvilinearGrid> grid = MakeCurvilinearGridRand(originX, originY, deltaX, deltaY, nx, ny, 0.4);

    const lin_alg::Matrix<mk::Point> originalPoints = grid->GetNodes();

    mk::CurvilinearGridSmoothing curvilinearGridSmoothing(*grid, 10);

    curvilinearGridSmoothing.SetBlock({ny * deltaX / 4.0, ny * deltaX / 4.0}, {3.0 * ny * deltaX / 4.0, 3.0 * ny * deltaX / 4.0});

    std::unique_ptr<mk::UndoAction> undoAction = curvilinearGridSmoothing.Compute();

    undoAction->Restore();

    constexpr double tolerance = 1.0e-12;

    ASSERT_EQ(originalPoints.rows(), grid->NumN());
    ASSERT_EQ(originalPoints.cols(), grid->NumM());

    for (mk::UInt r = 0; r < grid->NumN(); ++r)
    {
        for (mk::UInt c = 0; c < grid->NumM(); ++c)
        {
            EXPECT_NEAR(originalPoints(r, c).x, grid->GetNode(r, c).x, tolerance);
            EXPECT_NEAR(originalPoints(r, c).y, grid->GetNode(r, c).y, tolerance);
        }
    }
}

TEST(CurvilinearBasicTests, GridOrthogonalisation)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 1.0;

    constexpr size_t nx = 30;
    constexpr size_t ny = 30;

    std::unique_ptr<mk::CurvilinearGrid> grid = MakeCurvilinearGridRand(originX, originY, deltaX, deltaY, nx, ny, 0.4, false);
    grid->ComputeGridNodeTypes();

    const lin_alg::Matrix<mk::Point> originalPoints = grid->GetNodes();

    mk::OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 1;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;

    meshkernel::CurvilinearGridOrthogonalization curvilinearGridOrthogonalization(*grid, orthogonalizationParameters);

    curvilinearGridOrthogonalization.SetBlock({-10.0, -10.0}, {2.0 * nx * deltaX, 2.0 * nx * deltaX});

    std::unique_ptr<mk::UndoAction> undoAction = curvilinearGridOrthogonalization.Compute();

    undoAction->Restore();

    constexpr double tolerance = 1.0e-12;

    ASSERT_EQ(originalPoints.rows(), grid->NumN());
    ASSERT_EQ(originalPoints.cols(), grid->NumM());

    for (mk::UInt r = 0; r < grid->NumN(); ++r)
    {
        for (mk::UInt c = 0; c < grid->NumM(); ++c)
        {
            EXPECT_NEAR(originalPoints(r, c).x, grid->GetNode(r, c).x, tolerance) << "pos " << r << "  " << c;
            EXPECT_NEAR(originalPoints(r, c).y, grid->GetNode(r, c).y, tolerance) << "pos " << r << "  " << c;
        }
    }
}

TEST(CurvilinearBasicTests, Derefinement)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 1.0;

    constexpr size_t nx = 30;
    constexpr size_t ny = 30;

    std::unique_ptr<mk::CurvilinearGrid> grid = MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);
    grid->ComputeGridNodeTypes();

    const lin_alg::Matrix<mk::Point> originalPoints = grid->GetNodes();

    mk::CurvilinearGridDeRefinement deRefinement(*grid);

    deRefinement.SetBlock({10.0, 20.0}, {20.0, 20.0});

    std::unique_ptr<mk::UndoAction> undoAction = deRefinement.Compute();

    // After de-refinement the size of the grid in the y-direction should be the same
    EXPECT_TRUE(originalPoints.rows() == grid->NumN());
    // After de-refinement the size of the grid in the x-direction should be less
    EXPECT_TRUE(originalPoints.cols() > grid->NumM());

    undoAction->Restore();

    // After undoing the de-refinement the size of the grid in the y-direction should be the same
    ASSERT_TRUE(originalPoints.rows() == grid->NumN());
    // After undoing the de-refinement the size of the grid in the x-direction should be the same
    ASSERT_TRUE(originalPoints.cols() == grid->NumM());

    constexpr double tolerance = 1.0e-10;

    for (mk::UInt r = 0; r < grid->NumN(); ++r)
    {
        for (mk::UInt c = 0; c < grid->NumM(); ++c)
        {
            EXPECT_NEAR(originalPoints(r, c).x, grid->GetNode(r, c).x, tolerance);
            EXPECT_NEAR(originalPoints(r, c).y, grid->GetNode(r, c).y, tolerance);
        }
    }
}

TEST(CurvilinearBasicTests, UndoLineAttractor)
{
    auto grid = MakeSmallCurvilinearGrid();

    grid->ComputeGridNodeTypes();
    const std::vector<mk::Point> originalPoints = grid->ComputeNodes();

    mk::CurvilinearGridLineAttractionRepulsion lineAttractor(*grid, 0.5);

    lineAttractor.SetLine({80266.8, 367104.0}, {80419.3, 366566.2});
    lineAttractor.SetBlock({80198.2, 366750.6}, {80583.1, 366889.8});

    std::unique_ptr<mk::UndoAction> undoAction = lineAttractor.Compute();
    grid->ComputeGridNodeTypes();

    // Asserts
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(80178.014482303217, grid->GetNode(0, 2).x, tolerance);
    ASSERT_NEAR(80266.910680413363, grid->GetNode(1, 2).x, tolerance);
    ASSERT_NEAR(80322.584162464715, grid->GetNode(2, 2).x, tolerance);
    ASSERT_NEAR(80350.500795549306, grid->GetNode(3, 2).x, tolerance);
    ASSERT_NEAR(80362.879671417410, grid->GetNode(4, 2).x, tolerance);

    ASSERT_NEAR(367069.60110549850, grid->GetNode(0, 2).y, tolerance);
    ASSERT_NEAR(366937.57246542675, grid->GetNode(1, 2).y, tolerance);
    ASSERT_NEAR(366803.23746104678, grid->GetNode(2, 2).y, tolerance);
    ASSERT_NEAR(366683.98469820933, grid->GetNode(3, 2).y, tolerance);
    ASSERT_NEAR(366555.11052078847, grid->GetNode(4, 2).y, tolerance);

    undoAction->Restore();

    ASSERT_EQ(originalPoints.size(), grid->GetNumNodes());

    for (mk::UInt i = 0; i < originalPoints.size(); ++i)
    {
        EXPECT_NEAR(originalPoints[i].x, grid->Node(i).x, tolerance);
        EXPECT_NEAR(originalPoints[i].y, grid->Node(i).y, tolerance);
    }
}

TEST(CurvilinearBasicTests, UndoLineShift)
{
    // Set-up
    const auto grid = MakeSmallCurvilinearGrid();
    const std::vector<mk::Point> originalPoints = grid->ComputeNodes();

    meshkernel::CurvilinearGridLineShift curvilinearLineShift(*grid);

    curvilinearLineShift.SetLine({79982.0, 366934.0}, {80155.0, 366530.0});
    curvilinearLineShift.SetBlock({80108.0, 366707.0}, {80291.0, 366792.0});
    // The line shift is made by combining two actions
    // The first is to move the node.
    auto undoMoveNode = curvilinearLineShift.MoveNode({79982.0, 366934.0}, {79872.0, 366876.0});

    // The second actio is to shift the line
    auto undoLineShift = curvilinearLineShift.Compute();

    // Asserts
    const double tolerance = 1e-6;

    ASSERT_NEAR(79872.000000000000, grid->GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(80010.039799507853, grid->GetNode(0, 1).x, tolerance);
    ASSERT_NEAR(80145.970831448722, grid->GetNode(0, 2).x, tolerance);
    ASSERT_NEAR(80225.900042018140, grid->GetNode(0, 3).x, tolerance);
    ASSERT_NEAR(80305.243756829266, grid->GetNode(0, 4).x, tolerance);
    ASSERT_NEAR(80381.747982750283, grid->GetNode(0, 5).x, tolerance);
    ASSERT_NEAR(80458.252208671300, grid->GetNode(0, 6).x, tolerance);
    ASSERT_NEAR(80534.756434592317, grid->GetNode(0, 7).x, tolerance);
    ASSERT_NEAR(80611.260660513333, grid->GetNode(0, 8).x, tolerance);

    ASSERT_NEAR(79970.149644452977, grid->GetNode(1, 0).x, tolerance);
    ASSERT_NEAR(80103.062377666603, grid->GetNode(1, 1).x, tolerance);
    ASSERT_NEAR(80234.924195000407, grid->GetNode(1, 2).x, tolerance);
    ASSERT_NEAR(80324.671765221428, grid->GetNode(1, 3).x, tolerance);
    ASSERT_NEAR(80414.057391982613, grid->GetNode(1, 4).x, tolerance);
    ASSERT_NEAR(80505.096476712482, grid->GetNode(1, 5).x, tolerance);
    ASSERT_NEAR(80595.183339827883, grid->GetNode(1, 6).x, tolerance);
    ASSERT_NEAR(80684.333994102650, grid->GetNode(1, 7).x, tolerance);
    ASSERT_NEAR(80772.567299473958, grid->GetNode(1, 8).x, tolerance);

    ASSERT_NEAR(366876.00000000000, grid->GetNode(0, 0).y, tolerance);
    ASSERT_NEAR(366959.82623907487, grid->GetNode(0, 1).y, tolerance);
    ASSERT_NEAR(367047.36276461056, grid->GetNode(0, 2).y, tolerance);
    ASSERT_NEAR(367104.62934968271, grid->GetNode(0, 3).y, tolerance);
    ASSERT_NEAR(367163.01691965276, grid->GetNode(0, 4).y, tolerance);
    ASSERT_NEAR(367224.10904462705, grid->GetNode(0, 5).y, tolerance);
    ASSERT_NEAR(367285.20116960135, grid->GetNode(0, 6).y, tolerance);
    ASSERT_NEAR(367346.29329457565, grid->GetNode(0, 7).y, tolerance);
    ASSERT_NEAR(367407.38541954994, grid->GetNode(0, 8).y, tolerance);

    ASSERT_NEAR(366781.50715811126, grid->GetNode(1, 0).y, tolerance);
    ASSERT_NEAR(366849.28921837400, grid->GetNode(1, 1).y, tolerance);
    ASSERT_NEAR(366919.18816119258, grid->GetNode(1, 2).y, tolerance);
    ASSERT_NEAR(366966.76979594346, grid->GetNode(1, 3).y, tolerance);
    ASSERT_NEAR(367015.14849423966, grid->GetNode(1, 4).y, tolerance);
    ASSERT_NEAR(367056.48898898275, grid->GetNode(1, 5).y, tolerance);
    ASSERT_NEAR(367099.12347147451, grid->GetNode(1, 6).y, tolerance);
    ASSERT_NEAR(367143.03018172452, grid->GetNode(1, 7).y, tolerance);
    ASSERT_NEAR(367188.18349069095, grid->GetNode(1, 8).y, tolerance);

    // Need to undo both actions
    undoLineShift->Restore();
    undoMoveNode->Restore();

    ASSERT_EQ(originalPoints.size(), grid->GetNumNodes());

    for (mk::UInt i = 0; i < originalPoints.size(); ++i)
    {
        EXPECT_NEAR(originalPoints[i].x, grid->Node(i).x, tolerance);
        EXPECT_NEAR(originalPoints[i].y, grid->Node(i).y, tolerance);
    }
}

TEST(CurvilinearBasicTests, AnotherTest11)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 1.0;

    constexpr size_t nx = 30;
    constexpr size_t ny = 30;

    std::unique_ptr<mk::CurvilinearGrid> grid = MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);
    const std::vector<mk::Point> originalPoints = grid->ComputeNodes();

    meshkernel::CurvilinearGridLineShift lineShift(*grid);

    lineShift.SetLine({5.0, 3.0}, {5.0, 15.0});
    lineShift.SetBlock({2.0, 3.0}, {10.0, 10.0});
    std::unique_ptr<mk::UndoAction> undoMoveNode = lineShift.MoveNode({5.0, 5.0}, {6.5, 6.0});

    // Execute
    std::unique_ptr<mk::UndoAction> undoLineShift = lineShift.Compute();

    //--------------------------------
    // Now undo the actions

    undoLineShift->Restore();
    undoMoveNode->Restore();

    ASSERT_EQ(originalPoints.size(), grid->GetNumNodes());

    const double tolerance = 1e-12;

    for (mk::UInt i = 0; i < originalPoints.size(); ++i)
    {
        EXPECT_NEAR(originalPoints[i].x, grid->Node(i).x, tolerance);
        EXPECT_NEAR(originalPoints[i].y, grid->Node(i).y, tolerance);
    }
}

TEST(CurvilinearBasicTests, CompoundTest)
{
    // Tests combining a number of actions:
    //
    // line shift
    // mesh refinement (x-direction)
    // add boundary grid lines
    // delete interior block
    // more add boundary grid lines
    // mesh refinement (y-direction)
    // more add boundary grid lines
    // delete another interior block
    // undo all of the above, the grid nodes should be the same as the original
    // then redo all actions, the grid nodes should be the same as the nodes were before the undo, but after all the actions were computed.

    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 1.0;

    constexpr size_t nx = 30;
    constexpr size_t ny = 30;

    mk::UndoActionStack undoActions;
    std::unique_ptr<mk::CurvilinearGrid> grid = MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);
    grid->ComputeGridNodeTypes();

    // Nodes in the original mesh
    const std::vector<mk::Point> originalPoints = grid->ComputeNodes();

    meshkernel::CurvilinearGridLineShift lineShift(*grid);

    lineShift.SetLine({5.0, 3.0}, {5.0, 15.0});
    lineShift.SetBlock({2.0, 3.0}, {10.0, 10.0});

    undoActions.Add(lineShift.Compute({5.0, 5.0}, {6.5, 6.0}));

    mk::CurvilinearGridRefinement refinement(*grid, 2);
    refinement.SetBlock({10.0, 20.0}, {20.0, 20.0});

    undoActions.Add(refinement.Compute());

    //--------------------------------

    mk::CurvilinearGridLineMirror lineMirror(*grid, deltaX);

    // Bottom
    lineMirror.m_lines.push_back(mk::CurvilinearGridLine({0, 0}, {0, nx - 1 + 10}));
    undoActions.Add(lineMirror.Compute());
    grid->ComputeGridNodeTypes();

    //--------------------------------

    // Bottom
    lineMirror.m_lines[0] = mk::CurvilinearGridLine({0, 0}, {0, nx - 1 + 10});
    undoActions.Add(lineMirror.Compute());
    grid->ComputeGridNodeTypes();

    //--------------------------------

    mk::CurvilinearGridDeleteInterior deleteInterior(*grid);
    deleteInterior.m_lowerLeft = {7, 4};
    deleteInterior.m_upperRight = {16, 20};
    auto deleteAction2 = deleteInterior.Compute();
    grid->ComputeGridNodeTypes();
    undoActions.Add(std::move(deleteAction2));

    //--------------------------------

    // Left
    lineMirror.m_lines[0] = mk::CurvilinearGridLine({0, 0}, {ny + 1, 0});
    undoActions.Add(lineMirror.Compute());
    grid->ComputeGridNodeTypes();

    //--------------------------------

    // Top
    lineMirror.m_lines[0] = mk::CurvilinearGridLine({ny + 1, 5}, {ny + 1, nx + 7});
    undoActions.Add(lineMirror.Compute());
    grid->ComputeGridNodeTypes();

    //--------------------------------

    // Right
    lineMirror.m_lines[0] = mk::CurvilinearGridLine({3, nx + 10}, {ny - 4, nx + 10});
    undoActions.Add(lineMirror.Compute());
    grid->ComputeGridNodeTypes();

    //--------------------------------

    // Right
    lineMirror.m_lines[0] = mk::CurvilinearGridLine({3, nx + 11}, {ny - 4, nx + 11});
    undoActions.Add(lineMirror.Compute());
    grid->ComputeGridNodeTypes();

    //--------------------------------

    // Right
    lineMirror.m_lines[0] = mk::CurvilinearGridLine({3, nx + 12}, {ny - 4, nx + 12});
    undoActions.Add(lineMirror.Compute());
    grid->ComputeGridNodeTypes();

    //--------------------------------

    mk::CurvilinearGridRefinement refinement2(*grid, 2);
    refinement2.SetBlock({0.0, 10.0}, {0.0, 20.0});

    undoActions.Add(refinement2.Compute());

    //--------------------------------

    // Right
    lineMirror.m_lines[0] = mk::CurvilinearGridLine({3, nx + 13}, {ny - 4, nx + 13});
    undoActions.Add(lineMirror.Compute());
    grid->ComputeGridNodeTypes();

    //--------------------------------

    // Right
    lineMirror.m_lines[0] = mk::CurvilinearGridLine({3, nx + 14}, {ny - 4, nx + 14});
    undoActions.Add(lineMirror.Compute());
    grid->ComputeGridNodeTypes();

    //--------------------------------

    // Right
    lineMirror.m_lines[0] = mk::CurvilinearGridLine({3, nx + 15}, {ny - 4, nx + 15});
    undoActions.Add(lineMirror.Compute());
    grid->ComputeGridNodeTypes();

    //--------------------------------

    // Right
    lineMirror.m_lines[0] = mk::CurvilinearGridLine({3, nx + 16}, {ny - 4, nx + 16});
    undoActions.Add(lineMirror.Compute());
    grid->ComputeGridNodeTypes();

    //--------------------------------

    // Right
    lineMirror.m_lines[0] = mk::CurvilinearGridLine({3, nx + 17}, {ny - 4, nx + 17});
    undoActions.Add(lineMirror.Compute());
    grid->ComputeGridNodeTypes();

    //--------------------------------

    mk::CurvilinearGridDeleteInterior deleteInterior2(*grid);
    deleteInterior2.m_lowerLeft = {9, 25};
    deleteInterior2.m_upperRight = {23, 42};
    auto deleteAction3 = deleteInterior2.Compute();
    grid->ComputeGridNodeTypes();
    undoActions.Add(std::move(deleteAction3));

    //--------------------------------

    // Nodes in the grid after all actions
    const std::vector<mk::Point> refinedPoints = grid->ComputeNodes();

    bool didUndo;

    do
    {
        auto undoOption = undoActions.Undo();
        didUndo = static_cast<bool>(undoOption);
    } while (didUndo);

    constexpr double tolerance = 1.0e-12;

    // Points should be same as inthe original mesh after all actions have bene undone
    mk::UInt index = 0;

    for (mk::UInt i = 0; i < grid->GetNumNodes(); ++i)
    {
        if (grid->Node(i).IsValid())
        {
            EXPECT_NEAR(originalPoints[index].x, grid->Node(i).x, tolerance);
            EXPECT_NEAR(originalPoints[index].y, grid->Node(i).y, tolerance);
            ++index;
        }
    }

    bool didRedo;

    do
    {
        auto redoOption = undoActions.Commit();
        didRedo = static_cast<bool>(redoOption);
    } while (didRedo);

    // Points should be same as in the refined mesh after all actions have bene redone
    for (mk::UInt i = 0; i < grid->GetNumNodes(); ++i)
    {
        EXPECT_NEAR(refinedPoints[i].x, grid->Node(i).x, tolerance);
        EXPECT_NEAR(refinedPoints[i].y, grid->Node(i).y, tolerance);
    }
}

TEST(CurvilinearBasicTests, InsertMultipleFacesAlongAllBoundaryEdges)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 1.0;

    constexpr size_t nx = 11;
    constexpr size_t ny = 11;

    std::unique_ptr<mk::CurvilinearGrid> grid = MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);
    grid->ComputeGridNodeTypes();

    EXPECT_EQ(grid->NumN(), ny);
    EXPECT_EQ(grid->NumM(), nx);

    // Add element to bottom (south) of the grid
    // The number of rows in the array of points should be increased by 1
    [[maybe_unused]] auto undoAction1b = grid->InsertFace({0.5, 0.0});
    EXPECT_EQ(grid->NumN(), ny + 1);
    EXPECT_EQ(grid->NumM(), nx);

    // Add another element to bottom (south) of the grid
    // The number of rows in the array of points should not be increased
    [[maybe_unused]] auto undoAction2b = grid->InsertFace({3.5, 0.0});
    EXPECT_EQ(grid->NumN(), ny + 1);
    EXPECT_EQ(grid->NumM(), nx);

    // Add another element to bottom (south) of the grid
    // The number of rows in the array of points should not be increased
    [[maybe_unused]] auto undoAction3b = grid->InsertFace({2.5, 0.0});
    EXPECT_EQ(grid->NumN(), ny + 1);
    EXPECT_EQ(grid->NumM(), nx);

    //--------------------------------

    // Add element to left (west) of the grid
    // The number of columns in the array of points should be increased by 1
    [[maybe_unused]] auto undoAction1l = grid->InsertFace({0.0, 0.5});
    EXPECT_EQ(grid->NumN(), ny + 1);
    EXPECT_EQ(grid->NumM(), nx + 1);

    // Add another element to left (west) of the grid
    // The number of columns in the array of points should not be increased
    [[maybe_unused]] auto undoAction2l = grid->InsertFace({0.0, 3.5});
    EXPECT_EQ(grid->NumN(), ny + 1);
    EXPECT_EQ(grid->NumM(), nx + 1);

    // Add another element to left (west) of the grid
    // The number of columns in the array of points should not be increased
    [[maybe_unused]] auto undoAction3l = grid->InsertFace({0.0, 2.5});
    EXPECT_EQ(grid->NumN(), ny + 1);
    EXPECT_EQ(grid->NumM(), nx + 1);

    //--------------------------------

    // Add element to top (north) of the grid
    // The number of rows in the array of points should be increased by another 1
    [[maybe_unused]] auto undoAction1t = grid->InsertFace({0.5, 10.0});
    EXPECT_EQ(grid->NumN(), ny + 2);
    EXPECT_EQ(grid->NumM(), nx + 1);

    // Add another element to top (north) of the grid
    // The number of rows in the array of points should not be increased
    [[maybe_unused]] auto undoAction2t = grid->InsertFace({3.5, 10.0});
    EXPECT_EQ(grid->NumN(), ny + 2);
    EXPECT_EQ(grid->NumM(), nx + 1);

    // Add another element to top (north) of the grid
    // The number of rows in the array of points should not be increased
    [[maybe_unused]] auto undoAction3t = grid->InsertFace({2.5, 10.0});
    EXPECT_EQ(grid->NumN(), ny + 2);
    EXPECT_EQ(grid->NumM(), nx + 1);

    //--------------------------------

    // Add element to right (east) of the grid
    // The number of columns in the array of points should be increased by another 1
    [[maybe_unused]] auto undoAction1r = grid->InsertFace({10.0, 0.5});
    EXPECT_EQ(grid->NumN(), ny + 2);
    EXPECT_EQ(grid->NumM(), nx + 2);

    // Add another element to right (east) of the grid
    // The number of columns in the array of points should not be increased
    [[maybe_unused]] auto undoAction2r = grid->InsertFace({10.0, 3.5});
    EXPECT_EQ(grid->NumN(), ny + 2);
    EXPECT_EQ(grid->NumM(), nx + 2);

    // Add another element to right (east) of the grid
    // The number of columns in the array of points should not be increased
    [[maybe_unused]] auto undoAction3r = grid->InsertFace({10.0, 2.5});
    EXPECT_EQ(grid->NumN(), ny + 2);
    EXPECT_EQ(grid->NumM(), nx + 2);

    //--------------------------------
}
