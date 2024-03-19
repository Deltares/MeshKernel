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

#include <gtest/gtest.h>

#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridCurvature.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridDeleteInterior.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Operations.hpp"
#include "MeshKernel/UndoActions/UndoActionStack.hpp"
#include "MeshKernel/Utilities/LinearAlgebra.hpp"
#include "TestUtils/MakeCurvilinearGrids.hpp"

namespace mk = meshkernel;

namespace basic
{

    std::unique_ptr<meshkernel::CurvilinearGrid> MakeCurvilinearGrid(double originX, double originY, double deltaX, double deltaY, size_t nx, size_t ny)
    {
        double y = originY;

        lin_alg::Matrix<meshkernel::Point> points(ny, nx);

        for (size_t n = 0; n < ny; ++n)
        {
            double x = originX;
            double deltaXValue = deltaX;

            for (size_t m = 0; m < nx; ++m)
            {
                points(n, m) = meshkernel::Point(x, y);
                x += deltaXValue;
            }

            y += deltaY;
        }

        auto grid = std::make_unique<meshkernel::CurvilinearGrid>(points, meshkernel::Projection::cartesian);
        grid->SetFlatCopies();
        grid->ComputeGridNodeTypes();
        return grid;
    }

} // namespace basic

TEST(CurvilinearBasicTests, AddGridLineAtBoundary)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 2.0;

    constexpr size_t nx = 27;
    constexpr size_t ny = 13;

    std::unique_ptr<mk::CurvilinearGrid> grid = basic::MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);
    bool addedLine = false;

    std::cout << "pnts: " << grid->GetNode(0, 0).x << ", " << grid->GetNode(0, 0).y
              << "  -- " << grid->GetNode(1, 0).x << ", " << grid->GetNode(1, 0).y
              << std::endl;

    std::cout << " size of grid: " << grid->NumM() << "  " << grid->NumN() << "  " << std::boolalpha << addedLine
              << "  --- " << grid->FullNumM() << "  " << grid->FullNumN() << "  -- "
              << grid->m_startOffset.m_m << ", " << grid->m_startOffset.m_n << "  -- "
              << std::endl;

    grid->AddGridLineAtBoundary({0, 0}, {0, 1});
    grid->SetFlatCopies();

    for (mk::UInt i = 0; i < grid->NumN(); ++i)
    {
        grid->GetNode(0, i).x = grid->GetNode(1, i).x - deltaX;
        grid->GetNode(0, i).y = grid->GetNode(1, i).y;
        std::cout << " pnt: " << grid->GetNode(0, i).x << "  " << grid->GetNode(0, i).y << " -- " << grid->GetNode(1, i).x << "  " << grid->GetNode(1, i).y << std::endl;
    }

    grid->ComputeGridNodeTypes();

    std::cout << "pnts: " << grid->GetNode(0, 0).x << ", " << grid->GetNode(0, 0).y
              << "  -- " << grid->GetNode(1, 0).x << ", " << grid->GetNode(1, 0).y
              << std::endl;

    std::cout << " size of grid: " << grid->NumM() << "  " << grid->NumN() << "  " << std::boolalpha << addedLine
              << "  --- " << grid->FullNumM() << "  " << grid->FullNumN() << "  -- "
              << grid->m_startOffset.m_m << ", " << grid->m_startOffset.m_n << "  -- "
              << std::endl;

    // grid->AddGridLineAtBoundary({0, 26}, {1, 26});
    // grid->ComputeGridNodeTypes();
    // std::cout << " size of grid: " << grid->NumM() << "  " << grid->NumN() << "  " << std::boolalpha << addedLine << "  --- " << grid->FullNumM() << "  " << grid->FullNumN() << std::endl;
    // addedLine = grid->AddGridLineAtBoundary({0, 0}, {0, 1});
    // grid->ComputeGridNodeTypes();
    // std::cout << " size of grid: " << grid->NumM() << "  " << grid->NumN() << "  " << std::boolalpha << addedLine << "  --- " << grid->FullNumM() << "  " << grid->FullNumN() << std::endl;

    grid->DeleteGridLineAtBoundary({0, 0}, {1, 0});
    grid->SetFlatCopies();
    std::cout << "pnts: " << grid->GetNode(0, 0).x << ", " << grid->GetNode(0, 0).y << "  -- "
              << grid->GetNode(1, 0).x << ", " << grid->GetNode(1, 0).y << "  -- "
              << std::endl;
    std::cout << " size of grid: " << grid->NumM() << "  " << grid->NumN() << "  " << std::boolalpha << addedLine
              << "  --- " << grid->FullNumM() << "  " << grid->FullNumN() << "  -- "
              << grid->m_startOffset.m_m << ", " << grid->m_startOffset.m_n << "  -- "
              << std::endl;

    grid->DeleteGridLineAtBoundary({0, 0}, {0, 1});
    grid->SetFlatCopies();
    grid->ComputeGridNodeTypes();
    std::cout << "pnts: " << grid->GetNode(0, 0).x << ", " << grid->GetNode(0, 0).y << "  -- "
              << grid->GetNode(1, 0).x << ", " << grid->GetNode(1, 0).y << "  -- "
              << std::endl;
    std::cout << " size of grid: " << grid->NumM() << "  " << grid->NumN() << "  " << std::boolalpha << addedLine
              << "  --- " << grid->FullNumM() << "  " << grid->FullNumN() << "  -- "
              << grid->m_startOffset.m_m << ", " << grid->m_startOffset.m_n << "  -- "
              << std::endl;

    grid->DeleteGridLineAtBoundary({0, 0}, {0, 1});
    grid->SetFlatCopies();
    grid->ComputeGridNodeTypes();
    std::cout << "pnts: " << grid->GetNode(0, 0).x << ", " << grid->GetNode(0, 0).y << "  -- "
              << grid->GetNode(1, 0).x << ", " << grid->GetNode(1, 0).y << "  -- "
              << std::endl;
    std::cout << " size of grid: " << grid->NumM() << "  " << grid->NumN() << "  " << std::boolalpha << addedLine
              << "  --- " << grid->FullNumM() << "  " << grid->FullNumN() << "  -- "
              << grid->m_startOffset.m_m << ", " << grid->m_startOffset.m_n << "  -- "
              << std::endl;

    grid->AddGridLineAtBoundary({0, 0}, {1, 0});
    grid->SetFlatCopies();
    // grid->AddGridLineAtBoundary({0, 0}, {0, 1});
    grid->ComputeGridNodeTypes();
    std::cout << " size of grid: " << grid->NumM() << "  " << grid->NumN() << "  " << std::boolalpha << addedLine
              << "  --- " << grid->FullNumM() << "  " << grid->FullNumN() << "  -- "
              << grid->m_startOffset.m_m << ", " << grid->m_startOffset.m_n << "  -- "
              << std::endl;

    grid->AddGridLineAtBoundary({0, 0}, {0, 1});
    grid->SetFlatCopies();
    grid->ComputeGridNodeTypes();
    std::cout << " size of grid: " << grid->NumM() << "  " << grid->NumN() << "  " << std::boolalpha << addedLine
              << "  --- " << grid->FullNumM() << "  " << grid->FullNumN() << "  -- "
              << grid->m_startOffset.m_m << ", " << grid->m_startOffset.m_n << "  -- "
              << std::endl;

    grid->AddGridLineAtBoundary({0, 0}, {0, 1});
    grid->SetFlatCopies();
    grid->ComputeGridNodeTypes();
    std::cout << " size of grid: " << grid->NumM() << "  " << grid->NumN() << "  " << std::boolalpha << addedLine
              << "  --- " << grid->FullNumM() << "  " << grid->FullNumN() << "  -- "
              << grid->m_startOffset.m_m << ", " << grid->m_startOffset.m_n << "  -- "
              << std::endl;

    grid->AddGridLineAtBoundary({0, 0}, {0, 1});

    for (mk::UInt i = 0; i < grid->NumN(); ++i)
    {
        grid->GetNode(0, i).x = grid->GetNode(1, i).x - deltaX;
        grid->GetNode(0, i).y = grid->GetNode(1, i).y;
    }

    grid->SetFlatCopies();
    grid->ComputeGridNodeTypes();
    std::cout << " size of grid: " << grid->NumM() << "  " << grid->NumN() << "  " << std::boolalpha << addedLine
              << "  --- " << grid->FullNumM() << "  " << grid->FullNumN() << "  -- "
              << grid->m_startOffset.m_m << ", " << grid->m_startOffset.m_n << "  -- "
              << std::endl;

    grid->DeleteGridLineAtBoundary({0, 0}, {0, 1});
    grid->SetFlatCopies();
    grid->ComputeGridNodeTypes();
    std::cout << "pnts: " << grid->GetNode(0, 0).x << ", " << grid->GetNode(0, 0).y << "  -- "
              << grid->GetNode(1, 0).x << ", " << grid->GetNode(1, 0).y << "  -- "
              << std::endl;
    std::cout << " size of grid: " << grid->NumM() << "  " << grid->NumN() << "  " << std::boolalpha << addedLine
              << "  --- " << grid->FullNumM() << "  " << grid->FullNumN() << "  -- "
              << grid->m_startOffset.m_m << ", " << grid->m_startOffset.m_n << "  -- "
              << std::endl;

    grid->InsertFace({-1.1, 1.0});

    for (mk::UInt i = 0; i < grid->NumN(); ++i)
    {
        grid->GetNode(0, i).x = grid->GetNode(1, i).x - deltaX;
        grid->GetNode(0, i).y = grid->GetNode(1, i).y;
        std::cout << " pnt: " << grid->GetNode(0, i).x << "  " << grid->GetNode(0, i).y << " -- " << grid->GetNode(1, i).x << "  " << grid->GetNode(1, i).y << std::endl;
    }

    grid->SetFlatCopies();
    grid->ComputeGridNodeTypes();
    std::cout << "pnts: " << grid->GetNode(0, 0).x << ", " << grid->GetNode(0, 0).y << "  -- "
              << grid->GetNode(1, 0).x << ", " << grid->GetNode(1, 0).y << "  -- "
              << std::endl;
    std::cout << " size of grid: " << grid->NumM() << "  " << grid->NumN() << "  " << std::boolalpha << addedLine
              << "  --- " << grid->FullNumM() << "  " << grid->FullNumN() << "  -- "
              << grid->m_startOffset.m_m << ", " << grid->m_startOffset.m_n << "  -- "
              << std::endl;

    grid->DeleteGridLineAtBoundary({0, 0}, {0, 1});

    grid->SetFlatCopies();
    grid->ComputeGridNodeTypes();
    std::cout << "pnts: " << grid->GetNode(0, 0).x << ", " << grid->GetNode(0, 0).y << "  -- "
              << grid->GetNode(1, 0).x << ", " << grid->GetNode(1, 0).y << "  -- "
              << std::endl;
    std::cout << " size of grid: " << grid->NumM() << "  " << grid->NumN() << "  " << std::boolalpha << addedLine
              << "  --- " << grid->FullNumM() << "  " << grid->FullNumN() << "  -- "
              << grid->m_startOffset.m_m << ", " << grid->m_startOffset.m_n << "  -- "
              << std::endl;
}

TEST(CurvilinearBasicTests, AddGridLineAtBoundaryUndo)
{
    constexpr double originX = 1.0;
    constexpr double originY = 1.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 1.0;

    constexpr size_t nx = 8;
    constexpr size_t ny = 8;

    // constexpr size_t nx = 27;
    // constexpr size_t ny = 13;

    std::unique_ptr<mk::CurvilinearGrid> grid = basic::MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);
    bool addedLine = false;
    mk::UndoActionStack undoActions;

    grid->print();

    std::cout << " size of grid: " << grid->NumM() << "  " << grid->NumN() << "  " << std::boolalpha << addedLine
              << "  --- " << grid->FullNumM() << "  " << grid->FullNumN() << "  -- "
              << grid->m_startOffset.m_m << ", " << grid->m_startOffset.m_n << "  -- "
              << std::endl;

    // auto [addedAnotherLine, undoAction] = grid->AddGridLineAtBoundary({0, 7}, {1, 7}); ///< UP but on right
    auto [addedAnotherLine, undoAction] = grid->AddGridLineAtBoundary({7, 0}, {7, 1}); ///< Right but on up
    // auto [addedAnotherLine, undoAction] = grid->AddGridLineAtBoundary({0, 0}, {0, 1}); ///< Left but on bottom
    // auto [addedAnotherLine, undoAction] = grid->AddGridLineAtBoundary({0, 0}, {1, 0}); < Bottom but on left

    grid->SetFlatCopies();
    undoAction->Print();
    undoActions.Add(std::move(undoAction));

    // // up but right
    // for (mk::UInt i = 0; i < grid->NumN(); ++i)
    // {
    //     grid->GetNode(i, 8).x = grid->GetNode(i, 7).x + deltaX;
    //     grid->GetNode(i, 8).y = grid->GetNode(i, 7).y;
    // }

    // right but up
    for (mk::UInt i = 0; i < grid->NumM(); ++i)
    {
        grid->GetNode(ny, i).x = grid->GetNode(ny - 1, i).x;
        grid->GetNode(ny, i).y = grid->GetNode(ny - 1, i).y + deltaY;
    }

    // // For left but bottom
    // for (mk::UInt i = 0; i < grid->NumM(); ++i)
    // {
    //     grid->GetNode(0, i).x = grid->GetNode(1, i).x;
    //     grid->GetNode(0, i).y = grid->GetNode(1, i).y - deltaY;
    // }

    // // For bottom but left
    // for (mk::UInt i = 0; i < grid->NumN(); ++i)
    // {
    //     grid->GetNode(i, 0).x = grid->GetNode(i, 1).x - deltaX;
    //     grid->GetNode(i, 0).y = grid->GetNode(i, 1).y;
    // }

    grid->print();

    undoActions.Undo();

    grid->print();
    std::cout << "undoAction->Restore size of grid: " << grid->NumM() << "  " << grid->NumN() << "  " << std::boolalpha << addedLine
              << "  --- " << grid->FullNumM() << "  " << grid->FullNumN() << "  -- "
              << grid->m_startOffset.m_m << ", " << grid->m_startOffset.m_n << "  -- "
              << std::endl;

    undoActions.Commit();
    std::cout << "undoAction->Commit   size of grid: " << grid->NumM() << "  " << grid->NumN() << "  " << std::boolalpha << addedLine
              << "  --- " << grid->FullNumM() << "  " << grid->FullNumN() << "  -- "
              << grid->m_startOffset.m_m << ", " << grid->m_startOffset.m_n << "  -- "
              << std::endl;

    auto [addedAnotherLine2, undoAction2] = grid->AddGridLineAtBoundary({0, 0}, {0, 1});

    // For left but bottom
    for (mk::UInt i = 0; i < grid->NumM(); ++i)
    {
        grid->GetNode(0, i).x = grid->GetNode(1, i).x;
        grid->GetNode(0, i).y = grid->GetNode(1, i).y - deltaY;
    }

    grid->SetFlatCopies();

    grid->print();

    undoAction2->Print();
    undoActions.Add(std::move(undoAction2));

    undoActions.Undo();
    std::cout << "undoAction2->Restore size of grid: " << grid->NumM() << "  " << grid->NumN() << "  " << std::boolalpha << addedLine
              << "  --- " << grid->FullNumM() << "  " << grid->FullNumN() << "  -- "
              << grid->m_startOffset.m_m << ", " << grid->m_startOffset.m_n << "  -- "
              << std::endl;

    undoActions.Commit();
    std::cout << "undoAction->Restore  size of grid: " << grid->NumM() << "  " << grid->NumN() << "  " << std::boolalpha << addedLine
              << "  --- " << grid->FullNumM() << "  " << grid->FullNumN() << "  -- "
              << grid->m_startOffset.m_m << ", " << grid->m_startOffset.m_n << "  -- "
              << std::endl;
}


TEST(CurvilinearBasicTests, AnotherTest)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 1.0;

    constexpr size_t nx = 10;
    constexpr size_t ny = 10;

    std::unique_ptr<mk::CurvilinearGrid> grid = basic::MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);

    mk::CurvilinearGridDeleteInterior deleteInterior(*grid);
    deleteInterior.m_lowerLeft = {3,2};
    deleteInterior.m_upperRight = {6,7};
    auto deleteAction = deleteInterior.Compute2 ();
    // undoAction->Restore ();

    auto addGridLineAction = grid->AddGridLineAtBoundary({0, 0}, {0, 1});

    for (mk::UInt i = 0; i < grid->NumM(); ++i)
    {
        grid->GetNode(0, i).x = grid->GetNode(1, i).x;
        grid->GetNode(0, i).y = grid->GetNode(1, i).y - deltaY;
        std::cout << " add point: "<< grid->GetNode(0, i).x << ", " << grid->GetNode(0, i).y << std::endl;
    }

    // undoAction->Commit ();

    grid->SetFlatCopies ();
    grid->DeleteInvalidNodesAndEdges ();

    grid->printGraph ();



}

TEST(CurvilinearBasicTests, AnotherTest2)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 1.0;

    constexpr size_t nx = 10;
    constexpr size_t ny = 10;

    std::unique_ptr<mk::CurvilinearGrid> grid = basic::MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);

    auto addGridLineAction = grid->AddGridLineAtBoundary({0, 0}, {0, 1});

    for (mk::UInt i = 0; i < grid->NumM(); ++i)
    {
        grid->GetNode(0, i).x = grid->GetNode(1, i).x;
        grid->GetNode(0, i).y = grid->GetNode(1, i).y - deltaY;
        std::cout << " add point: "<< grid->GetNode(0, i).x << ", " << grid->GetNode(0, i).y << std::endl;
    }

    mk::CurvilinearGridDeleteInterior deleteInterior(*grid);
    deleteInterior.m_lowerLeft = {3,2};
    deleteInterior.m_upperRight = {6,7};
    auto deleteAction = deleteInterior.Compute2 ();
    // undoAction->Restore ();


    // undoAction->Commit ();

    grid->SetFlatCopies ();
    grid->DeleteInvalidNodesAndEdges ();

    grid->printGraph ();



}

TEST(CurvilinearBasicTests, AnotherTest3)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 1.0;

    constexpr size_t nx = 10;
    constexpr size_t ny = 10;

    mk::UndoActionStack undoActions;

    std::unique_ptr<mk::CurvilinearGrid> grid = basic::MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);

    //--------------------------------

    auto [lineAdded1, addGridLineAction1] = grid->AddGridLineAtBoundary({0, 0}, {0, 1});
    grid->ComputeGridNodeTypes ();
    undoActions.Add (std::move (addGridLineAction1));

    for (mk::UInt i = 0; i < grid->NumM(); ++i)
    {
        grid->GetNode(0, i).x = grid->GetNode(1, i).x;
        grid->GetNode(0, i).y = grid->GetNode(1, i).y - deltaY;
        std::cout << " add point: "<< grid->GetNode(0, i).x << ", " << grid->GetNode(0, i).y << std::endl;
    }

    grid->ComputeGridNodeTypes ();

    //--------------------------------

    mk::CurvilinearGridDeleteInterior deleteInterior(*grid);
    deleteInterior.m_lowerLeft = {3,2};
    deleteInterior.m_upperRight = {6,7};
    auto deleteAction1 = deleteInterior.Compute2 ();
    undoActions.Add (std::move (deleteAction1));

    //--------------------------------

    auto [lineAdded2, addGridLineAction2] = grid->AddGridLineAtBoundary({0, 0}, {0, 1});
    undoActions.Add (std::move (addGridLineAction2));

    for (mk::UInt i = 0; i < grid->NumM(); ++i)
    {
        grid->GetNode(0, i).x = grid->GetNode(1, i).x;
        grid->GetNode(0, i).y = grid->GetNode(1, i).y - deltaY;
        std::cout << " add point: "<< grid->GetNode(0, i).x << ", " << grid->GetNode(0, i).y << std::endl;
    }

    grid->ComputeGridNodeTypes ();

    //--------------------------------

    auto [lineAdded3, addGridLineAction3] = grid->AddGridLineAtBoundary({0, 0}, {1, 0});
    undoActions.Add (std::move (addGridLineAction3));

    for (mk::UInt i = 0; i < grid->NumN(); ++i)
    {
        grid->GetNode(i, 0).x = grid->GetNode(i, 1).x - deltaX;
        grid->GetNode(i, 0).y = grid->GetNode(i, 1).y;
        std::cout << " add point: "<< grid->GetNode(i, 0).x << ", " << grid->GetNode(i, 0).y << std::endl;
    }

    grid->ComputeGridNodeTypes ();

    //--------------------------------

    // auto [lineAdded4, addGridLineAction4] = grid->AddGridLineAtBoundary({0, 9}, {1, 9});
    auto [lineAdded4, addGridLineAction4] = grid->AddGridLineAtBoundary({0, 10}, {1, 10});
    undoActions.Add (std::move (addGridLineAction4));

    for (mk::UInt i = 0; i < grid->NumN(); ++i)
    {
        grid->GetNode(i, 11).x = grid->GetNode(i, 10).x - deltaX;
        grid->GetNode(i, 11).y = grid->GetNode(i, 10).y;
        std::cout << " add point: "<< grid->GetNode(i, 11).x << ", " << grid->GetNode(i, 11).y << std::endl;
    }

    grid->ComputeGridNodeTypes ();

    //--------------------------------

    deleteInterior.m_lowerLeft = {2,3};
    deleteInterior.m_upperRight = {9,6};
    auto deleteAction2 = deleteInterior.Compute2 ();
    grid->ComputeGridNodeTypes ();
    undoActions.Add (std::move (deleteAction2));

    //--------------------------------

    undoActions.Undo (); // Undo last delete interior action
    undoActions.Undo (); // Undo second add grid line
    undoActions.Undo (); // Undo last add grid line
    undoActions.Undo (); // Undo first delete interior
    undoActions.Undo (); // Undo first add grid line

    undoActions.Commit (); // Redo first add grid line
    undoActions.Commit (); // Redo first delete interior
    undoActions.Commit (); //
    undoActions.Commit (); //


    grid->SetFlatCopies ();
    grid->DeleteInvalidNodesAndEdges ();

    grid->printGraph ();



}
