//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2021.
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

#pragma once
#include <memory>

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernelApi/CurvilinearGrid.hpp>

/// @brief Counts the number of valid nodes in a meshkernelapi::CurvilinearGrid instance
/// @param curvilinearGrid The input curvilinear grid
/// @return The number of valid nodes
size_t CurvilinearGridCountValidNodes(meshkernelapi::CurvilinearGrid const& curvilinearGrid);

/// @brief Counts the number of valid nodes in a meshkernel::CurvilinearGrid instance
/// @param The input curvilinear grid
/// @return The number of valid nodes
size_t CurvilinearGridCountValidNodes(meshkernel::CurvilinearGrid const& curvilinearGrid);

/// @brief Makes a small, real world curvi grid
/// See tests/CurvilinearGrids/Orthogonalization/small_curvi_grid.png for a plot of the grid
/// @return A pointer to a curvilinear grid
std::unique_ptr<meshkernel::CurvilinearGrid> MakeSmallCurvilinearGrid();

/// @brief Makes a small, real world curvi grid, with missing faces
/// @return A pointer to a curvilinear grid
std::unique_ptr<meshkernel::CurvilinearGrid> MakeSmallCurvilinearGridWithMissingFaces();

/// @brief Makes a curvilinear grid
/// @return A pointer to a curvilinear grid
std::unique_ptr<meshkernel::CurvilinearGrid> MakeCurvilinearGrid(double originX, double originY, double deltaX, double deltaY, size_t nx, size_t ny);
