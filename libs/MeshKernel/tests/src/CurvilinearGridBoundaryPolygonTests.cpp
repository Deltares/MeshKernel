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

#include "MakeCurvilinearGrids.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"

#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>

using namespace meshkernel;

TEST(CurvilinearGridBoundaryPolygon, ComputeBoundaryToPolygon_OnValidCurvilinearGrid_ShouldComputeBoundaryPolygon)

{
    // Prepare
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 1.0;

    constexpr UInt nx = 3;
    constexpr UInt ny = 3;

    std::unique_ptr<CurvilinearGrid> grid = MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);
    CurvilinearGridNodeIndices lowerLeft(0, 0);
    CurvilinearGridNodeIndices upperRight(ny - 1, nx - 1);

    // Execute

    const auto boundaryPolygon = grid->ComputeBoundaryPolygons(lowerLeft, upperRight);

    // Assert
    std::vector<Point> expectedBoundaryPolygon{{0.0, 0.0},
                                               {1.0, 0.0},
                                               {2.0, 0.0},
                                               {2.0, 1.0},
                                               {2.0, 2.0},
                                               {1.0, 2.0},
                                               {0.0, 2.0},
                                               {0.0, 1.0},
                                               {0.0, 0.0}};
    ASSERT_THAT(boundaryPolygon, ::testing::ContainerEq(expectedBoundaryPolygon));
}

TEST(CurvilinearGridBoundaryPolygon, ComputeBoundaryToPolygon_OnValidCurvilinearGridWithInvalidNodes_ShouldComputeBoundaryPolygon)

{
    // Prepare
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 1.0;

    constexpr UInt nx = 4;
    constexpr UInt ny = 4;

    std::unique_ptr<CurvilinearGrid> grid = MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);
    auto undoAction = grid->DeleteNode({3.0, 0.0});
    undoAction = grid->DeleteNode({3.0, 3.0});
    CurvilinearGridNodeIndices lowerLeft(0, 0);
    CurvilinearGridNodeIndices upperRight(ny - 1, nx - 1);

    // Execute
    const auto boundaryPolygon = grid->ComputeBoundaryPolygons(lowerLeft, upperRight);

    // Assert
    std::vector<Point> expectedBoundaryPolygon{{0.0, 0.0},
                                               {1.0, 0.0},
                                               {2.0, 0.0},
                                               {2.0, 1.0},
                                               {3.0, 1.0},
                                               {3.0, 2.0},
                                               {2.0, 2.0},
                                               {2.0, 3.0},
                                               {1.0, 3.0},
                                               {0.0, 3.0},
                                               {0.0, 2.0},
                                               {0.0, 1.0},
                                               {0.0, 0.0}};

    ASSERT_THAT(boundaryPolygon, ::testing::ContainerEq(expectedBoundaryPolygon));
}

TEST(CurvilinearGridBoundaryPolygon, ComputeBoundaryToPolygon_OnValidCurvilinearGridWithDisconnectedCurvilinearGrid_ShouldComputeBoundaryPolygons)

{
    // Prepare
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 1.0;

    constexpr UInt nx = 9;
    constexpr UInt ny = 3;

    std::unique_ptr<CurvilinearGrid> grid = MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);
    auto undoAction = grid->DeleteNode({3.0, 1.0});
    undoAction = grid->DeleteNode({4.0, 1.0});
    CurvilinearGridNodeIndices lowerLeft(0, 0);
    CurvilinearGridNodeIndices upperRight(ny - 1, nx - 1);

    // Execute
    const auto boundaryPolygon = grid->ComputeBoundaryPolygons(lowerLeft, upperRight);

    // Assert
    std::vector<Point> expectedBoundaryPolygon{{0.0, 0.0},
                                               {1.0, 0.0},
                                               {2.0, 0.0},
                                               {2.0, 1.0},
                                               {2.0, 2.0},
                                               {1.0, 2.0},
                                               {0.0, 2.0},
                                               {0.0, 1.0},
                                               {0.0, 0.0},
                                               {-999.0, -999.0},
                                               {5.0, 0.0},
                                               {6.0, 0.0},
                                               {7.0, 0.0},
                                               {8.0, 0.0},
                                               {8.0, 1.0},
                                               {8.0, 2.0},
                                               {7.0, 2.0},
                                               {6.0, 2.0},
                                               {5.0, 2.0},
                                               {5.0, 1.0},
                                               {5.0, 0.0}};

    ASSERT_THAT(boundaryPolygon, ::testing::ContainerEq(expectedBoundaryPolygon));
}

TEST(CurvilinearGridBoundaryPolygon, ComputeBoundaryToPolygon_OnValidCurvilinearGridWithDisconnectedCurvilinearGridAndSelection_ShouldComputeBoundaryPolygons)

{
    // Prepare
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 1.0;

    constexpr UInt nx = 9;
    constexpr UInt ny = 3;

    std::unique_ptr<CurvilinearGrid> grid = MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);
    auto undoAction = grid->DeleteNode({3.0, 1.0});
    undoAction = grid->DeleteNode({4.0, 1.0});
    CurvilinearGridNodeIndices lowerLeft(0, 0);
    CurvilinearGridNodeIndices upperRight(1, nx - 1); // not the entire clg

    // Execute
    const auto boundaryPolygon = grid->ComputeBoundaryPolygons(lowerLeft, upperRight);

    // Assert
    std::vector<Point> expectedBoundaryPolygon{{0.0, 0.0},
                                               {1.0, 0.0},
                                               {2.0, 0.0},
                                               {2.0, 1.0},
                                               {1.0, 1.0},
                                               {0.0, 1.0},
                                               {0.0, 0.0},
                                               {-999.0, -999.0},
                                               {5.0, 0.0},
                                               {6.0, 0.0},
                                               {7.0, 0.0},
                                               {8.0, 0.0},
                                               {8.0, 1.0},
                                               {7.0, 1.0},
                                               {6.0, 1.0},
                                               {5.0, 1.0},
                                               {5.0, 0.0}};

    ASSERT_THAT(boundaryPolygon, ::testing::ContainerEq(expectedBoundaryPolygon));
}
