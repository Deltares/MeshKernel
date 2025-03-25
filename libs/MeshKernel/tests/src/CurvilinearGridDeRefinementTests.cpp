//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2025.
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

#include "MeshKernel/CurvilinearGrid/CurvilinearGridRefinement.hpp"

#include <gtest/gtest.h>

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridDeRefinement.hpp>
#include <MeshKernel/Entities.hpp>

using namespace meshkernel;

TEST(CurvilinearGridDeRefinement, Compute_OnCurvilinearGrid_ShouldDeRefineVerticalGridLines)
{
    // Set-up
    lin_alg::Matrix<Point> grid(4, 5);
    grid << Point{0, 0}, Point{10, 0}, Point{15, 0}, Point{20, 0}, Point{30, 0},
        Point{0, 10}, Point{10, 10}, Point{15, 10}, Point{20, 10}, Point{30, 10},
        Point{0, 20}, Point{10, 20}, Point{15, 20}, Point{20, 20}, Point{30, 20},
        Point{0, 30}, Point{10, 30}, Point{15, 30}, Point{20, 30}, Point{30, 30};

    CurvilinearGrid curvilinearGrid(grid, Projection::cartesian);
    CurvilinearGridDeRefinement curvilinearGridDeRefinement(curvilinearGrid, 2);
    curvilinearGridDeRefinement.SetBlock(Point{10, 20}, Point{20, 20});

    // Execute
    [[maybe_unused]] auto dummyUndoAction = curvilinearGridDeRefinement.Compute();

    // Assert (the vertical line at x=15 is removed)
    ASSERT_EQ(4, curvilinearGrid.NumM());
    ASSERT_EQ(4, curvilinearGrid.NumN());
}

TEST(CurvilinearGridDeRefinement, Compute_OnCurvilinearGridWithMissingFaces_ShouldDeRefineVerticalGridLines)
{
    // Set-up
    lin_alg::Matrix<Point> grid(4, 9);

    grid << Point{0, 0}, Point{10, 0}, Point{}, Point{}, Point{}, Point{}, Point{}, Point{40, 0}, Point{50, 0},
        Point{0, 10}, Point{10, 10}, Point{}, Point{}, Point{}, Point{}, Point{}, Point{40, 10}, Point{50, 10},
        Point{0, 20}, Point{10, 20}, Point{11, 20}, Point{12, 20}, Point{13, 20}, Point{20, 20}, Point{30, 20}, Point{40, 20}, Point{50, 20},
        Point{0, 30}, Point{10, 30}, Point{11, 30}, Point{12, 30}, Point{13, 30}, Point{20, 30}, Point{30, 30}, Point{40, 30}, Point{50, 30};

    CurvilinearGrid curvilinearGrid(grid, Projection::cartesian);
    CurvilinearGridDeRefinement curvilinearGridDeRefinement(curvilinearGrid, 4);
    curvilinearGridDeRefinement.SetBlock(Point{10, 20}, Point{20, 20});

    // Execute
    [[maybe_unused]] auto dummyUndoAction = curvilinearGridDeRefinement.Compute();

    // Assert
    ASSERT_EQ(6, curvilinearGrid.NumM());
    ASSERT_EQ(4, curvilinearGrid.NumN());
}

TEST(CurvilinearGridDeRefinement, Compute_OnCurvilinearGrid_ShouldDeRefineHorizontalGridLines)
{
    // Set-up
    lin_alg::Matrix<Point> grid(5, 4);

    grid << Point{0, 0}, Point{10, 0}, Point{20, 0}, Point{30, 0},
        Point{0, 10}, Point{10, 10}, Point{20, 10}, Point{30, 10},
        Point{0, 11}, Point{10, 11}, Point{20, 11}, Point{30, 11},
        Point{0, 20}, Point{10, 20}, Point{20, 20}, Point{30, 20},
        Point{0, 30}, Point{10, 30}, Point{20, 30}, Point{30, 30};

    CurvilinearGrid curvilinearGrid(grid, Projection::cartesian);
    CurvilinearGridDeRefinement curvilinearGridDeRefinement(curvilinearGrid, 2);
    curvilinearGridDeRefinement.SetBlock(Point{10, 10}, Point{10, 20});

    // Execute
    [[maybe_unused]] auto dummyUndoAction = curvilinearGridDeRefinement.Compute();

    // Assert (the vertical line at x=15 is removed)
    ASSERT_EQ(4, curvilinearGrid.NumN());
    ASSERT_EQ(4, curvilinearGrid.NumM());
}

TEST(CurvilinearGridDeRefinement, Compute_OnRefinedCurvilinearGridWithLargeDerefinementFactor_ShouldDeRefineVerticalGridLines)
{
    // Set-up
    lin_alg::Matrix<Point> grid(10, 10);

    // Fill the grid
    for (int y = 0; y < grid.rows(); ++y)
    {
        for (int x = 0; x < grid.cols(); ++x)
        {
            grid(y, x) = Point{x * 10.0, y * 10.0};
        }
    }

    CurvilinearGrid curvilinearGrid(grid, Projection::cartesian);
    CurvilinearGridRefinement curvilinearGridRefinement(curvilinearGrid, 5);
    curvilinearGridRefinement.SetBlock(Point{20, 0}, Point{80, 0});
    [[maybe_unused]] auto dummyUndoAction = curvilinearGridRefinement.Compute();

    CurvilinearGridDeRefinement curvilinearGridDeRefinement(curvilinearGrid, 100);
    curvilinearGridDeRefinement.SetBlock(Point{20, 0}, Point{80, 0});

    // Execute
    dummyUndoAction = curvilinearGridDeRefinement.Compute();

    // Assert, given the large de-refinement factor all lines between are removed
    ASSERT_EQ(10, curvilinearGrid.NumN());
    ASSERT_EQ(5, curvilinearGrid.NumM());
}

TEST(CurvilinearGridDeRefinement, Compute_OnRefinedCurvilinearGridWithSubsequentDerefinements_ShouldDeRefineVerticallGridLinesAndKeepRightBoundary)
{
    // Set-up
    lin_alg::Matrix<Point> grid(11, 11);

    // Fill the grid
    for (int y = 0; y < grid.rows(); ++y)
    {
        for (int x = 0; x < grid.cols(); ++x)
        {
            grid(y, x) = Point{x * 10.0, y * 10.0};
        }
    }

    CurvilinearGrid curvilinearGrid(grid, Projection::cartesian);

    for (int iter = 0; iter < 4; ++iter)
    {
        CurvilinearGridDeRefinement curvilinearGridDeRefinement(curvilinearGrid, 2);
        curvilinearGridDeRefinement.SetBlock(Point{0, 0}, Point{100, 0});
        [[maybe_unused]] auto dummyUndoAction = curvilinearGridDeRefinement.Compute();
    }

    // Assert, given the large de-refinement factor all lines between are removed
    ASSERT_EQ(11, curvilinearGrid.NumN());
    ASSERT_EQ(2, curvilinearGrid.NumM());

    constexpr double tolerance = 1e-6;
    ASSERT_NEAR(0, curvilinearGrid.GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(100.0, curvilinearGrid.GetNode(0, 1).x, tolerance);
}
