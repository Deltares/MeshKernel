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

#include <MeshKernel/CurvilinearGridAlgorithm.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Splines.hpp>

namespace meshkernel
{
    class CurvilinearGrid;

    /// @brief A class implementing the curvilinear grid refinement algorithm
    class CurvilinearGridRefinement : public CurvilinearGridAlgorithm
    {
    public:
        /// @brief Class constructor
        ///
        /// \p firstPoint and \p secondPoint must lie on the same gridline
        /// @param[in] grid The input curvilinear grid
        /// @param[in] refinement  The number of refinement lines between the points set by SetBlock()
        CurvilinearGridRefinement(const std::shared_ptr<CurvilinearGrid>& grid, size_t refinement);

        /// @brief Refine the curvilinear grid
        CurvilinearGrid Compute() override;

    private:
        size_t m_refinement; ///< The selected number of refinement lines
        Splines m_splines;   ///< An instance of the spline class storing the individual grid lines as splines
    };
} // namespace meshkernel
