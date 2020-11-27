//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2020.
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
    /// @brief A class representing a curvilinear grid
    class CurvilinearGrid
    {

    public:
        /// @brief Create a new curvilinear grid
        /// @param[in] m Number of lines
        /// @param[in] n Number of points per line
        void Set(int m, int n)
        {

            int mMax = m + 1;
            int nMax = n + 1;

            m_grid.resize(mMax, std::vector<Point>(nMax, {doubleMissingValue, doubleMissingValue}));
        }

        /// @brief Assign point to the curvilinear grid
        /// @param[in] grid Input grid
        void Set(const std::vector<std::vector<Point>>& grid)
        {
            m_grid = grid;
        }

        std::vector<std::vector<Point>> m_grid; ///< Member variable storing the grid
    };
} // namespace meshkernel
