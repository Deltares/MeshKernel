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

#include <vector>

#include <MeshKernel/Entities.hpp>

namespace meshkernel
{
    class Polygons;

    /// @brief A class representing a curvilinear grid
    class CurvilinearGrid
    {

    public:
        /// @brief Default constructor
        /// @returns
        CurvilinearGrid() = default;

        /// @brief Create a new curvilinear grid
        /// @param[in] m Number of columns (horizontal direction)
        /// @param[in] n Number of rows (vertical direction)
        CurvilinearGrid(size_t m, size_t n)
        {
            m_grid.resize(m + 1, std::vector<Point>(n + 1, {doubleMissingValue, doubleMissingValue}));
        }

        /// @brief Deletes a curvilinear grid inside a polygon
        /// @param[in] The polygons
        /// @param[in] The index of the polygon to use for deletion
        void Delete(Polygons const& polygons, size_t polygonIndex);

        /// @brief Sets the point to the curvilinear grid
        /// @param[in] grid Input grid points
        CurvilinearGrid(const std::vector<std::vector<Point>>& grid)
        {
            CurvilinearGrid(grid.size(), grid[0].size());
            m_grid = grid;
        }

        std::vector<std::vector<Point>> m_grid; ///< Member variable storing the grid
    };
} // namespace meshkernel
