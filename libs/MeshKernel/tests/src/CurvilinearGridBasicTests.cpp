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
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Operations.hpp"
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
