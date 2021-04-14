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

#include <MeshKernel/Entities.hpp>
#include <memory>

namespace meshkernel
{
    class CurvilinearGrid;

    /// @brief A class implementing the curvilinear grid de-refinement algorithm.
    /// A segment is defined by the first and second point.
    /// The grid lines crossed by the segment are eliminated from the curvilinear grid.
    class CurvilinearGridDeRefinement
    {
    public:
        /// @brief Class constructor
        /// @param[in] grid The input curvilinear grid
        /// @param[in] firstPoint The first node of the segment defining the derefinement zone
        /// @param[in] secondPoint The second node of the segment defining the derefinement zone
        CurvilinearGridDeRefinement(std::shared_ptr<CurvilinearGrid> grid,
                                    const Point& firstPoint,
                                    const Point& secondPoint);

        /// @brief Refine the curvilinear grid
        CurvilinearGrid Compute() const;

    private:
        std::shared_ptr<CurvilinearGrid> m_grid; ///< A pointer to the curvilinear grid to modify
        Point m_firstPoint;                      ///< The first node of the segment defining the derefinement zone
        Point m_secondPoint;                     ///< The second node of the segment defining the derefinement zone
    };
} // namespace meshkernel