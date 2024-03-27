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

    std::unique_ptr<meshkernel::CurvilinearGrid> MakeCurvilinearGridRand(double originX, double originY, double deltaX, double deltaY, size_t nx, size_t ny)
    {
        double y = originY;

        std::uniform_real_distribution<double> dist(0.0, 1.0);
        std::default_random_engine re;

        lin_alg::Matrix<meshkernel::Point> points(ny, nx);

        for (size_t n = 0; n < ny; ++n)
        {
            double x = originX;
            double deltaXValue = deltaX;

            for (size_t m = 0; m < nx; ++m)
            {
                points(n, m) = meshkernel::Point(x, y) + mk::Vector(dist(re) * 0.4 * deltaX, dist(re) * 0.4 * deltaY);
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

TEST(CurvilinearBasicTests, DISABLED_AddGridLineAtBoundary)
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

    auto undoAction = grid->InsertFace({-1.1, 1.0});

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

TEST(CurvilinearBasicTests, DISABLED_AddGridLineAtBoundaryUndo)
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

TEST(CurvilinearBasicTests, DISABLED_AnotherTest)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 1.0;

    constexpr size_t nx = 10;
    constexpr size_t ny = 10;

    std::unique_ptr<mk::CurvilinearGrid> grid = basic::MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);

    mk::CurvilinearGridDeleteInterior deleteInterior(*grid);
    deleteInterior.m_lowerLeft = {3, 2};
    deleteInterior.m_upperRight = {6, 7};
    auto deleteAction = deleteInterior.Compute();
    // undoAction->Restore ();

    auto addGridLineAction = grid->AddGridLineAtBoundary({0, 0}, {0, 1});

    for (mk::UInt i = 0; i < grid->NumM(); ++i)
    {
        grid->GetNode(0, i).x = grid->GetNode(1, i).x;
        grid->GetNode(0, i).y = grid->GetNode(1, i).y - deltaY;
        std::cout << " add point: " << grid->GetNode(0, i).x << ", " << grid->GetNode(0, i).y << std::endl;
    }

    // undoAction->Commit ();

    grid->SetFlatCopies();
    grid->DeleteInvalidNodesAndEdges();

    grid->printGraph();
}

TEST(CurvilinearBasicTests, DISABLED_AnotherTest2)
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
        std::cout << " add point: " << grid->GetNode(0, i).x << ", " << grid->GetNode(0, i).y << std::endl;
    }

    mk::CurvilinearGridDeleteInterior deleteInterior(*grid);
    deleteInterior.m_lowerLeft = {3, 2};
    deleteInterior.m_upperRight = {6, 7};
    auto deleteAction = deleteInterior.Compute();
    // undoAction->Restore ();

    // undoAction->Commit ();

    grid->SetFlatCopies();
    grid->DeleteInvalidNodesAndEdges();

    grid->printGraph();
}

TEST(CurvilinearBasicTests, DISABLED_AnotherTest3)
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
    grid->ComputeGridNodeTypes();
    undoActions.Add(std::move(addGridLineAction1));

    for (mk::UInt i = 0; i < grid->NumM(); ++i)
    {
        grid->GetNode(0, i).x = grid->GetNode(1, i).x;
        grid->GetNode(0, i).y = grid->GetNode(1, i).y - deltaY;
        std::cout << " add point: " << grid->GetNode(0, i).x << ", " << grid->GetNode(0, i).y << std::endl;
    }

    grid->ComputeGridNodeTypes();

    //--------------------------------

    mk::CurvilinearGridDeleteInterior deleteInterior(*grid);
    deleteInterior.m_lowerLeft = {3, 2};
    deleteInterior.m_upperRight = {6, 7};
    auto deleteAction1 = deleteInterior.Compute();
    undoActions.Add(std::move(deleteAction1));

    //--------------------------------

    auto [lineAdded2, addGridLineAction2] = grid->AddGridLineAtBoundary({0, 0}, {0, 1});
    undoActions.Add(std::move(addGridLineAction2));

    for (mk::UInt i = 0; i < grid->NumM(); ++i)
    {
        grid->GetNode(0, i).x = grid->GetNode(1, i).x;
        grid->GetNode(0, i).y = grid->GetNode(1, i).y - deltaY;
        std::cout << " add point: " << grid->GetNode(0, i).x << ", " << grid->GetNode(0, i).y << std::endl;
    }

    grid->ComputeGridNodeTypes();

    //--------------------------------

    auto [lineAdded3, addGridLineAction3] = grid->AddGridLineAtBoundary({0, 0}, {1, 0});
    undoActions.Add(std::move(addGridLineAction3));

    for (mk::UInt i = 0; i < grid->NumN(); ++i)
    {
        grid->GetNode(i, 0).x = grid->GetNode(i, 1).x - deltaX;
        grid->GetNode(i, 0).y = grid->GetNode(i, 1).y;
        std::cout << " add point: " << grid->GetNode(i, 0).x << ", " << grid->GetNode(i, 0).y << std::endl;
    }

    grid->ComputeGridNodeTypes();

    //--------------------------------

    // auto [lineAdded4, addGridLineAction4] = grid->AddGridLineAtBoundary({0, 9}, {1, 9});
    auto [lineAdded4, addGridLineAction4] = grid->AddGridLineAtBoundary({0, 10}, {1, 10});
    undoActions.Add(std::move(addGridLineAction4));

    for (mk::UInt i = 0; i < grid->NumN(); ++i)
    {
        grid->GetNode(i, 11).x = grid->GetNode(i, 10).x + deltaX;
        grid->GetNode(i, 11).y = grid->GetNode(i, 10).y;
        std::cout << " add point: " << grid->GetNode(i, 11).x << ", " << grid->GetNode(i, 11).y << std::endl;
    }

    grid->ComputeGridNodeTypes();

    //--------------------------------

    deleteInterior.m_lowerLeft = {2, 4};
    deleteInterior.m_upperRight = {9, 6};
    auto deleteAction2 = deleteInterior.Compute();
    grid->ComputeGridNodeTypes();
    undoActions.Add(std::move(deleteAction2));

    //--------------------------------

    undoActions.Undo(); // Undo last delete interior action
    undoActions.Undo(); // Undo add grid line (0,10)(1,10)
    undoActions.Undo(); // Undo add grid line (0,0)(1,0)

    std::cout << "undo action stack: " << undoActions.Size() << "  " << undoActions.CommittedSize() << "  " << undoActions.RestoredSize() << std::endl;

    grid->ComputeGridNodeTypes();

    //--------------------------------

    // auto [lineAdded4, addGridLineAction4] = grid->AddGridLineAtBoundary({0, 9}, {1, 9});
    auto [lineAdded5, addGridLineAction5] = grid->AddGridLineAtBoundary({0, 9}, {1, 9});
    undoActions.Add(std::move(addGridLineAction5));
    std::cout << "undo action stack: " << undoActions.Size() << "  " << undoActions.CommittedSize() << "  " << undoActions.RestoredSize() << std::endl;

    for (mk::UInt i = 0; i < grid->NumN(); ++i)
    {
        grid->GetNode(i, 10).x = grid->GetNode(i, 9).x + deltaX;
        grid->GetNode(i, 10).y = grid->GetNode(i, 9).y;
        std::cout << " add point: " << grid->GetNode(i, 10).x << ", " << grid->GetNode(i, 10).y << std::endl;
    }

    grid->ComputeGridNodeTypes();

    //--------------------------------

    auto [lineAdded6, addGridLineAction6] = grid->AddGridLineAtBoundary({0, 0}, {1, 0});
    undoActions.Add(std::move(addGridLineAction6));

    for (mk::UInt i = 0; i < grid->NumN(); ++i)
    {
        grid->GetNode(i, 0).x = grid->GetNode(i, 1).x - deltaX;
        grid->GetNode(i, 0).y = grid->GetNode(i, 1).y;
        std::cout << " add point: " << grid->GetNode(i, 0).x << ", " << grid->GetNode(i, 0).y << std::endl;
    }

    grid->ComputeGridNodeTypes();

    // undoActions.Undo(); //
    // undoActions.Undo(); //
    // undoActions.Undo(); // Undo first add grid line

    // undoActions.Commit(); // Redo first add grid line
    // undoActions.Commit(); // Redo first delete interior
    // undoActions.Commit(); //
    // undoActions.Commit(); //
    // undoActions.Commit(); //
    // undoActions.Commit(); //

    grid->SetFlatCopies();
    grid->DeleteInvalidNodesAndEdges();

    grid->printGraph();
}

TEST(CurvilinearBasicTests, DISABLED_AnotherTest4)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 1.0;

    constexpr size_t nx = 50;
    constexpr size_t ny = 50;

    mk::UndoActionStack undoActions;

    std::unique_ptr<mk::CurvilinearGrid> grid = basic::MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);

    //--------------------------------

    auto [lineAdded1, addGridLineAction1] = grid->AddGridLineAtBoundary({0, 0}, {0, 1});
    grid->ComputeGridNodeTypes();
    undoActions.Add(std::move(addGridLineAction1));

    for (mk::UInt i = 0; i < grid->NumM(); ++i)
    {
        grid->GetNode(0, i).x = grid->GetNode(1, i).x;
        grid->GetNode(0, i).y = grid->GetNode(1, i).y - deltaY;
        std::cout << " add point: " << grid->GetNode(0, i).x << ", " << grid->GetNode(0, i).y << std::endl;
    }

    grid->ComputeGridNodeTypes();

    //--------------------------------

    mk::CurvilinearGridDeleteInterior deleteInterior(*grid);
    mk::CurvilinearGridDeleteExterior deleteExterior(*grid);

    deleteInterior.m_lowerLeft = {20, 15};
    deleteInterior.m_upperRight = {30, 30};
    auto deleteAction1 = deleteInterior.Compute();
    undoActions.Add(std::move(deleteAction1));

    //--------------------------------

    auto [lineAdded2, addGridLineAction2] = grid->AddGridLineAtBoundary({0, 0}, {0, 1});
    undoActions.Add(std::move(addGridLineAction2));

    for (mk::UInt i = 0; i < grid->NumM(); ++i)
    {
        grid->GetNode(0, i).x = grid->GetNode(1, i).x;
        grid->GetNode(0, i).y = grid->GetNode(1, i).y - deltaY;
        std::cout << " add point: " << grid->GetNode(0, i).x << ", " << grid->GetNode(0, i).y << std::endl;
    }

    grid->ComputeGridNodeTypes();

    //--------------------------------

    auto [lineAdded3, addGridLineAction3] = grid->AddGridLineAtBoundary({0, 0}, {1, 0});
    undoActions.Add(std::move(addGridLineAction3));

    for (mk::UInt i = 0; i < grid->NumN(); ++i)
    {
        grid->GetNode(i, 0).x = grid->GetNode(i, 1).x - deltaX;
        grid->GetNode(i, 0).y = grid->GetNode(i, 1).y;
        std::cout << " add point: " << grid->GetNode(i, 0).x << ", " << grid->GetNode(i, 0).y << std::endl;
    }

    grid->ComputeGridNodeTypes();

    //--------------------------------

    // auto [lineAdded4, addGridLineAction4] = grid->AddGridLineAtBoundary({0, 9}, {1, 9});
    auto [lineAdded4, addGridLineAction4] = grid->AddGridLineAtBoundary({0, nx}, {1, nx});
    undoActions.Add(std::move(addGridLineAction4));

    for (mk::UInt i = 0; i < grid->NumN(); ++i)
    {
        grid->GetNode(i, nx + 1).x = grid->GetNode(i, nx).x + deltaX;
        grid->GetNode(i, nx + 1).y = grid->GetNode(i, nx).y;
        std::cout << " add point: " << grid->GetNode(i, nx + 1).x << ", " << grid->GetNode(i, nx + 1).y << std::endl;
    }

    grid->ComputeGridNodeTypes();

    // //--------------------------------

    deleteExterior.m_lowerLeft = {5, 5};
    deleteExterior.m_upperRight = {43, 45};
    auto deleteAction3 = deleteExterior.Compute();
    undoActions.Add(std::move(deleteAction3));

    undoActions.Undo();

    // deleteInterior.m_lowerLeft = {2, 4};
    // deleteInterior.m_upperRight = {9, 6};
    // auto deleteAction2 = deleteInterior.Compute();
    // grid->ComputeGridNodeTypes();
    // undoActions.Add(std::move(deleteAction2));

    // //--------------------------------

    // undoActions.Undo(); // Undo last delete interior action
    // undoActions.Undo(); // Undo add grid line (0,10)(1,10)
    // undoActions.Undo(); // Undo add grid line (0,0)(1,0)

    // std::cout << "undo action stack: " << undoActions.Size() << "  " << undoActions.CommittedSize() << "  " << undoActions.RestoredSize() << std::endl;

    // grid->ComputeGridNodeTypes();

    // //--------------------------------

    // // auto [lineAdded4, addGridLineAction4] = grid->AddGridLineAtBoundary({0, 9}, {1, 9});
    // auto [lineAdded5, addGridLineAction5] = grid->AddGridLineAtBoundary({0, 9}, {1, 9});
    // undoActions.Add(std::move(addGridLineAction5));
    // std::cout << "undo action stack: " << undoActions.Size() << "  " << undoActions.CommittedSize() << "  " << undoActions.RestoredSize() << std::endl;

    // for (mk::UInt i = 0; i < grid->NumN(); ++i)
    // {
    //     grid->GetNode(i, 10).x = grid->GetNode(i, 9).x + deltaX;
    //     grid->GetNode(i, 10).y = grid->GetNode(i, 9).y;
    //     std::cout << " add point: " << grid->GetNode(i, 10).x << ", " << grid->GetNode(i, 10).y << std::endl;
    // }

    // grid->ComputeGridNodeTypes();

    // //--------------------------------

    // auto [lineAdded6, addGridLineAction6] = grid->AddGridLineAtBoundary({0, 0}, {1, 0});
    // undoActions.Add(std::move(addGridLineAction6));

    // for (mk::UInt i = 0; i < grid->NumN(); ++i)
    // {
    //     grid->GetNode(i, 0).x = grid->GetNode(i, 1).x - deltaX;
    //     grid->GetNode(i, 0).y = grid->GetNode(i, 1).y;
    //     std::cout << " add point: " << grid->GetNode(i, 0).x << ", " << grid->GetNode(i, 0).y << std::endl;
    // }

    // grid->ComputeGridNodeTypes();

    // // undoActions.Undo(); //
    // // undoActions.Undo(); //
    // // undoActions.Undo(); // Undo first add grid line

    // // undoActions.Commit(); // Redo first add grid line
    // // undoActions.Commit(); // Redo first delete interior
    // // undoActions.Commit(); //
    // // undoActions.Commit(); //
    // // undoActions.Commit(); //
    // // undoActions.Commit(); //

    grid->SetFlatCopies();
    grid->DeleteInvalidNodesAndEdges();

    grid->printGraph();
}

TEST(CurvilinearBasicTests, DISABLED_AnotherTest5)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 1.0;

    constexpr size_t nx = 10;
    constexpr size_t ny = 10;

    mk::UndoActionStack undoActions;

    std::unique_ptr<mk::CurvilinearGrid> grid = basic::MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);

    mk::CurvilinearGridDeleteInterior deleteInterior(*grid);
    // mk::CurvilinearGridDeleteExterior deleteExterior(*grid);

    //--------------------------------

    mk::CurvilinearGridLineMirror lineMirror(*grid, deltaX);

    // Bottom
    lineMirror.m_lines.push_back(mk::CurvilinearGridLine({0, 0}, {0, nx - 1}));
    undoActions.Add(lineMirror.Compute());
    grid->ComputeGridNodeTypes();

    //--------------------------------

    // Bottom
    lineMirror.m_lines[0] = mk::CurvilinearGridLine({0, 0}, {0, nx - 1});
    undoActions.Add(lineMirror.Compute());
    grid->ComputeGridNodeTypes();

    //--------------------------------

    deleteInterior.m_lowerLeft = {2, 4};
    deleteInterior.m_upperRight = {9, 6};
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
    lineMirror.m_lines[0] = mk::CurvilinearGridLine({ny + 1, 0}, {ny + 1, nx});
    undoActions.Add(lineMirror.Compute());
    grid->ComputeGridNodeTypes();

    //--------------------------------

    // Right
    lineMirror.m_lines[0] = mk::CurvilinearGridLine({1, nx}, {ny + 1, nx});
    undoActions.Add(lineMirror.Compute());
    grid->ComputeGridNodeTypes();

    //--------------------------------

    undoActions.Undo();
    undoActions.Undo();
    undoActions.Undo();
    undoActions.Undo();

    undoActions.Commit();
    undoActions.Commit();
    undoActions.Commit();
    undoActions.Commit();

    // undoActions.Undo();

    // undoActions.Undo(); // Undo add grid line (0,10)(1,10)
    // undoActions.Undo(); // Undo add grid line (0,0)(1,0)

    // std::cout << "undo action stack: " << undoActions.Size() << "  " << undoActions.CommittedSize() << "  " << undoActions.RestoredSize() << std::endl;

    // grid->ComputeGridNodeTypes();

    // //--------------------------------

    // // auto [lineAdded4, addGridLineAction4] = grid->AddGridLineAtBoundary({0, 9}, {1, 9});
    // auto [lineAdded5, addGridLineAction5] = grid->AddGridLineAtBoundary({0, 9}, {1, 9});
    // undoActions.Add(std::move(addGridLineAction5));
    // std::cout << "undo action stack: " << undoActions.Size() << "  " << undoActions.CommittedSize() << "  " << undoActions.RestoredSize() << std::endl;

    // for (mk::UInt i = 0; i < grid->NumN(); ++i)
    // {
    //     grid->GetNode(i, 10).x = grid->GetNode(i, 9).x + deltaX;
    //     grid->GetNode(i, 10).y = grid->GetNode(i, 9).y;
    //     std::cout << " add point: " << grid->GetNode(i, 10).x << ", " << grid->GetNode(i, 10).y << std::endl;
    // }

    // grid->ComputeGridNodeTypes();

    // //--------------------------------

    // auto [lineAdded6, addGridLineAction6] = grid->AddGridLineAtBoundary({0, 0}, {1, 0});
    // undoActions.Add(std::move(addGridLineAction6));

    // for (mk::UInt i = 0; i < grid->NumN(); ++i)
    // {
    //     grid->GetNode(i, 0).x = grid->GetNode(i, 1).x - deltaX;
    //     grid->GetNode(i, 0).y = grid->GetNode(i, 1).y;
    //     std::cout << " add point: " << grid->GetNode(i, 0).x << ", " << grid->GetNode(i, 0).y << std::endl;
    // }

    // grid->ComputeGridNodeTypes();

    // // undoActions.Undo(); //
    // // undoActions.Undo(); //
    // // undoActions.Undo(); // Undo first add grid line

    // // undoActions.Commit(); // Redo first add grid line
    // // undoActions.Commit(); // Redo first delete interior
    // // undoActions.Commit(); //
    // // undoActions.Commit(); //
    // // undoActions.Commit(); //
    // // undoActions.Commit(); //

    grid->SetFlatCopies();
    grid->DeleteInvalidNodesAndEdges();

    grid->printGraph();
}

TEST(CurvilinearBasicTests, DISABLED_AnotherTest6)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 1.0;

    constexpr size_t nx = 30;
    constexpr size_t ny = 30;

    std::unique_ptr<mk::CurvilinearGrid> grid = basic::MakeCurvilinearGridRand(originX, originY, deltaX, deltaY, nx, ny);

    grid->ComputeGridNodeTypes();
    grid->SetFlatCopies();
    std::vector<mk::Point> points = grid->Nodes();

    mk::CurvilinearGridSmoothing curvilinearGridSmoothing(*grid, 10);

    curvilinearGridSmoothing.SetBlock({ny * deltaX / 4.0, ny * deltaX / 4.0}, {3.0 * ny * deltaX / 4.0, 3.0 * ny * deltaX / 4.0});

    std::unique_ptr<mk::UndoAction> undoAction = curvilinearGridSmoothing.Compute();

    // undoAction->Restore();
    grid->SetFlatCopies();

    grid->printGraph();

    // for (mk::UInt i = 0; i < points.size(); ++i)
    // {
    //     std::cout << " nodes "
    //               << points[i].x << " ++ " << grid->Node(i).x << " == " << points[i].x - grid->Node(i).x << " ----- "
    //               << points[i].y << " ++ " << grid->Node(i).y << " == " << points[i].y - grid->Node(i).y
    //               << std::endl;
    // }
}

TEST(CurvilinearBasicTests, DISABLED_AnotherTest7)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 1.0;

    constexpr size_t nx = 30;
    constexpr size_t ny = 30;

    std::unique_ptr<mk::CurvilinearGrid> grid = basic::MakeCurvilinearGridRand(originX, originY, deltaX, deltaY, nx, ny);

    grid->ComputeGridNodeTypes();
    grid->SetFlatCopies();
    std::vector<mk::Point> points = grid->Nodes();

    mk::OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 1;
    orthogonalizationParameters.boundary_iterations = 1;
    orthogonalizationParameters.inner_iterations = 1;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;

    meshkernel::CurvilinearGridOrthogonalization curvilinearGridOrthogonalization(*grid, orthogonalizationParameters);

    curvilinearGridOrthogonalization.SetBlock({-10.0, -10.0}, {2.0 * nx * deltaX, 2.0 * nx * deltaX});
    // curvilinearGridOrthogonalization.SetBlock({ny * deltaX / 4.0, ny * deltaX / 4.0}, {3.0 * ny * deltaX / 4.0, 3.0 * ny * deltaX / 4.0});

    std::unique_ptr<mk::UndoAction> undoAction = curvilinearGridOrthogonalization.Compute();

    // undoAction->Restore();
    grid->SetFlatCopies();

    grid->printGraph();

    // for (mk::UInt i = 0; i < points.size(); ++i)
    // {
    //     std::cout << " nodes "
    //               << points[i].x << " ++ " << grid->Node(i).x << " == " << points[i].x - grid->Node(i).x << " ----- "
    //               << points[i].y << " ++ " << grid->Node(i).y << " == " << points[i].y - grid->Node(i).y
    //               << std::endl;
    // }
}

TEST(CurvilinearBasicTests, DISABLED_AnotherTest8)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 1.0;

    constexpr size_t nx = 30;
    constexpr size_t ny = 30;

    std::unique_ptr<mk::CurvilinearGrid> grid = basic::MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);

    grid->ComputeGridNodeTypes();
    grid->SetFlatCopies();
    std::vector<mk::Point> points = grid->Nodes();

    mk::CurvilinearGridDeRefinement deRefinement(*grid);

    // refinement.SetBlock({-10.0, -10.0}, {2.0 * nx * deltaX, 2.0 * nx * deltaX});
    // refinement.SetBlock({ny * deltaX / 4.0, 3.0 * ny * deltaX / 4.0}, {3.0 * ny * deltaX / 4.0, 3.0 * ny * deltaX / 4.0});
    deRefinement.SetBlock({10.0, 20.0}, {20.0, 20.0});

    std::unique_ptr<mk::UndoAction> undoActionDeRef = deRefinement.Compute();
    grid->ComputeGridNodeTypes();
    grid->SetFlatCopies();

    mk::CurvilinearGridRefinement refinement(*grid, 3);
    refinement.SetBlock({10.0, 20.0}, {20.0, 20.0});

    std::unique_ptr<mk::UndoAction> undoActionRef = refinement.Compute();
    grid->ComputeGridNodeTypes();
    grid->SetFlatCopies();

    undoActionRef->Restore();
    undoActionDeRef->Restore();

    undoActionDeRef->Commit();
    undoActionRef->Commit();

    grid->SetFlatCopies();
    grid->printGraph();

    // for (mk::UInt i = 0; i < points.size(); ++i)
    // {
    //     std::cout << " nodes "
    //               << points[i].x << " ++ " << grid->Node(i).x << " == " << points[i].x - grid->Node(i).x << " ----- "
    //               << points[i].y << " ++ " << grid->Node(i).y << " == " << points[i].y - grid->Node(i).y
    //               << std::endl;
    // }
}

TEST(CurvilinearBasicTests, DISABLED_AnotherTest9)
{
    [[maybe_unused]] constexpr double originX = 0.0;
    [[maybe_unused]] constexpr double originY = 0.0;

    [[maybe_unused]] constexpr double deltaX = 1.0;
    [[maybe_unused]] constexpr double deltaY = 1.0;

    [[maybe_unused]] constexpr size_t nx = 30;
    [[maybe_unused]] constexpr size_t ny = 30;

    auto grid = MakeSmallCurvilinearGrid();
    // std::unique_ptr<mk::CurvilinearGrid> grid = basic::MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);

    grid->ComputeGridNodeTypes();
    grid->SetFlatCopies();
    std::vector<mk::Point> points = grid->Nodes();

    mk::CurvilinearGridLineAttractionRepulsion lineAttractor(*grid, 0.5);

    // mk::CurvilinearGridRefinement refinement(*grid, 3);
    // refinement.SetBlock({10.0, 20.0}, {20.0, 20.0});

    // std::unique_ptr<mk::UndoAction> undoActionRef = refinement.Compute();

    // 8.022590004201814008E+04, 3.671046293496827129E+05

    lineAttractor.SetLine({80266.8, 367104.0}, {80419.3, 366566.2});
    lineAttractor.SetBlock({80198.2, 366750.6}, {80583.1, 366889.8});

    // // lineAttractor.SetBlock({-10.0, -10.0}, {2.0 * nx * deltaX, 2.0 * nx * deltaX});
    // // lineAttractor.SetBlock({ny * deltaX / 4.0, 3.0 * ny * deltaX / 4.0}, {3.0 * ny * deltaX / 4.0, 3.0 * ny * deltaX / 4.0});
    // lineAttractor.SetLine({20.0, 10.0}, {10.0, 10.0});
    // lineAttractor.SetBlock({10.0, 20.0}, {20.0, 20.0});

    std::unique_ptr<mk::UndoAction> undoAction = lineAttractor.Compute();
    undoAction->Restore();
    grid->ComputeGridNodeTypes();
    grid->SetFlatCopies();

    grid->printGraph();

    // for (mk::UInt i = 0; i < points.size(); ++i)
    // {
    //     std::cout << " nodes "
    //               << points[i].x << " ++ " << grid->Node(i).x << " == " << points[i].x - grid->Node(i).x << " ----- "
    //               << points[i].y << " ++ " << grid->Node(i).y << " == " << points[i].y - grid->Node(i).y
    //               << std::endl;
    // }
}

TEST(CurvilinearBasicTests, DISABLED_AnotherTest10)
{
    std::cout.precision(18);
    // Set-up
    const auto grid = MakeSmallCurvilinearGrid();
    meshkernel::CurvilinearGridLineShift curvilinearLineShift(*grid);

    // grid->printGraph();

    curvilinearLineShift.SetLine({79982.0, 366934.0}, {80155.0, 366530.0});
    curvilinearLineShift.SetBlock({80108.0, 366707.0}, {80291.0, 366792.0});
    auto undoMoveNode = curvilinearLineShift.MoveNode({79982.0, 366934.0}, {79872.0, 366876.0});

    // curvilinearLineShift.SetLine({79982.0, 366934.0}, {80155.0, 366530.0});
    // curvilinearLineShift.SetBlock({80108.0, 366707.0}, {80291.0, 366792.0});

    // // Move two nodes
    // curvilinearLineShift.MoveNode({79982.0, 366934.0}, {79872.0, 366876.0});
    // curvilinearLineShift.MoveNode({80053.0, 366823.0}, {79932.0, 366773.0});

    // Execute
    auto undoLineShift = curvilinearLineShift.Compute();
    grid->ComputeGridNodeTypes();
    grid->SetFlatCopies();

    // Asserts
    // const double tolerance = 1e-6;

    // ASSERT_NEAR(79872.000000000000, grid->GetNode(0, 0).x, tolerance);
    // ASSERT_NEAR(80010.039799507853, grid->GetNode(0, 1).x, tolerance);
    // ASSERT_NEAR(80145.970831448722, grid->GetNode(0, 2).x, tolerance);
    // ASSERT_NEAR(80225.900042018140, grid->GetNode(0, 3).x, tolerance);
    // ASSERT_NEAR(80305.243756829266, grid->GetNode(0, 4).x, tolerance);
    // ASSERT_NEAR(80381.747982750283, grid->GetNode(0, 5).x, tolerance);
    // ASSERT_NEAR(80458.252208671300, grid->GetNode(0, 6).x, tolerance);
    // ASSERT_NEAR(80534.756434592317, grid->GetNode(0, 7).x, tolerance);
    // ASSERT_NEAR(80611.260660513333, grid->GetNode(0, 8).x, tolerance);

    // ASSERT_NEAR(79970.149644452977, grid->GetNode(1, 0).x, tolerance);
    // ASSERT_NEAR(80103.062377666603, grid->GetNode(1, 1).x, tolerance);
    // ASSERT_NEAR(80234.924195000407, grid->GetNode(1, 2).x, tolerance);
    // ASSERT_NEAR(80324.671765221428, grid->GetNode(1, 3).x, tolerance);
    // ASSERT_NEAR(80414.057391982613, grid->GetNode(1, 4).x, tolerance);
    // ASSERT_NEAR(80505.096476712482, grid->GetNode(1, 5).x, tolerance);
    // ASSERT_NEAR(80595.183339827883, grid->GetNode(1, 6).x, tolerance);
    // ASSERT_NEAR(80684.333994102650, grid->GetNode(1, 7).x, tolerance);
    // ASSERT_NEAR(80772.567299473958, grid->GetNode(1, 8).x, tolerance);

    // ASSERT_NEAR(366876.00000000000, grid->GetNode(0, 0).y, tolerance);
    // ASSERT_NEAR(366959.82623907487, grid->GetNode(0, 1).y, tolerance);
    // ASSERT_NEAR(367047.36276461056, grid->GetNode(0, 2).y, tolerance);
    // ASSERT_NEAR(367104.62934968271, grid->GetNode(0, 3).y, tolerance);
    // ASSERT_NEAR(367163.01691965276, grid->GetNode(0, 4).y, tolerance);
    // ASSERT_NEAR(367224.10904462705, grid->GetNode(0, 5).y, tolerance);
    // ASSERT_NEAR(367285.20116960135, grid->GetNode(0, 6).y, tolerance);
    // ASSERT_NEAR(367346.29329457565, grid->GetNode(0, 7).y, tolerance);
    // ASSERT_NEAR(367407.38541954994, grid->GetNode(0, 8).y, tolerance);

    // ASSERT_NEAR(366781.50715811126, grid->GetNode(1, 0).y, tolerance);
    // ASSERT_NEAR(366849.28921837400, grid->GetNode(1, 1).y, tolerance);
    // ASSERT_NEAR(366919.18816119258, grid->GetNode(1, 2).y, tolerance);
    // ASSERT_NEAR(366966.76979594346, grid->GetNode(1, 3).y, tolerance);
    // ASSERT_NEAR(367015.14849423966, grid->GetNode(1, 4).y, tolerance);
    // ASSERT_NEAR(367056.48898898275, grid->GetNode(1, 5).y, tolerance);
    // ASSERT_NEAR(367099.12347147451, grid->GetNode(1, 6).y, tolerance);
    // ASSERT_NEAR(367143.03018172452, grid->GetNode(1, 7).y, tolerance);
    // ASSERT_NEAR(367188.18349069095, grid->GetNode(1, 8).y, tolerance);

    grid->printGraph();
}

TEST(CurvilinearBasicTests, DISABLED_AnotherTest11)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 1.0;

    constexpr size_t nx = 30;
    constexpr size_t ny = 30;

    std::unique_ptr<mk::CurvilinearGrid> grid = basic::MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);

    meshkernel::CurvilinearGridLineShift lineShift(*grid);

    // grid->printGraph();

    lineShift.SetLine({5.0, 3.0}, {5.0, 15.0});
    lineShift.SetBlock({2.0, 3.0}, {10.0, 10.0});
    std::unique_ptr<mk::UndoAction> undoMoveNode = lineShift.MoveNode({5.0, 5.0}, {6.5, 6.0});

    // lineShift.SetLine({79982.0, 366934.0}, {80155.0, 366530.0});
    // lineShift.SetBlock({80108.0, 366707.0}, {80291.0, 366792.0});

    // // Move two nodes
    // lineShift.MoveNode({79982.0, 366934.0}, {79872.0, 366876.0});
    // lineShift.MoveNode({80053.0, 366823.0}, {79932.0, 366773.0});

    // Execute
    std::unique_ptr<mk::UndoAction> undoAction = lineShift.Compute();

    undoAction->Restore();
    undoMoveNode->Restore();

    grid->ComputeGridNodeTypes();
    grid->SetFlatCopies();
    grid->printGraph();
}

TEST(CurvilinearBasicTests, DISABLED_AnotherTest12)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 1.0;

    constexpr size_t nx = 30;
    constexpr size_t ny = 30;

    std::unique_ptr<mk::CurvilinearGrid> grid = basic::MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);

    auto undoAction1 = grid->InsertFace ({0.5, 30.0});
    undoAction1->Restore ();

    auto undoAction2 = grid->InsertFace ({1.5, 30.0});
    // undoAction2->Restore ();

    grid->ComputeGridNodeTypes();
    grid->SetFlatCopies();
    grid->printGraph();
}

TEST(CurvilinearBasicTests, DISABLED_CompoundTest1)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 1.0;

    constexpr size_t nx = 30;
    constexpr size_t ny = 30;

    mk::UndoActionStack undoActions;
    std::unique_ptr<mk::CurvilinearGrid> grid = basic::MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);

    meshkernel::CurvilinearGridLineShift lineShift(*grid);

    // grid->printGraph();

    lineShift.SetLine({5.0, 3.0}, {5.0, 15.0});
    lineShift.SetBlock({2.0, 3.0}, {10.0, 10.0});

    undoActions.Add(lineShift.MoveNode({5.0, 5.0}, {6.5, 6.0}));
    undoActions.Add(lineShift.Compute());

    mk::CurvilinearGridRefinement refinement(*grid, 2);
    refinement.SetBlock({10.0, 20.0}, {20.0, 20.0});

    undoActions.Add(refinement.Compute());

    //--------------------------------

    mk::CurvilinearGridLineMirror lineMirror(*grid, deltaX);

    // Bottom
    lineMirror.m_lines.push_back(mk::CurvilinearGridLine({0, 0}, {0, nx - 1 + 10}));
    undoActions.Add(lineMirror.Compute());
    grid->ComputeGridNodeTypes();

#if 1

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

#endif

    while (undoActions.Undo())
    {
        // Nothing else to do
    }

    while (undoActions.Commit())
    {
        // Nothing else to do
    }

    // undoActions.Undo();
    // undoActions.Undo();
    // undoActions.Undo();
    // undoActions.Undo();
    // undoActions.Undo();
    // undoActions.Undo();
    // undoActions.Undo();
    // undoActions.Undo();
    // undoActions.Undo();

    grid->ComputeGridNodeTypes();
    grid->SetFlatCopies();
    grid->printGraph();
}
