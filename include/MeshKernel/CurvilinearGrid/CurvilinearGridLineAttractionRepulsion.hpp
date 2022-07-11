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
#include <MeshKernel/CurvilinearGrid/CurvilinearGridAlgorithm.hpp>

namespace meshkernel
{

    /// @brief A class implementing the curvilinear grid line shift. Line is shifted on the new locations.
    /// and the displacement gets distributed on the influence zone, set with SetBlock
    ///
    /// This option provides the possibility to fit the curvilinear grid’s edges to a land boundary.
    class CurvilinearGridLineAttractionRepulsion : public CurvilinearGridAlgorithm
    {
    public:
        /// @brief Class constructor
        /// @param[in] grid             The input curvilinear grid
        /// @param[in] attractionFactor The input attraction/repulsion factor (positive attracts lines/negative repulses lines)
        CurvilinearGridLineAttractionRepulsion(std::shared_ptr<CurvilinearGrid> grid, double attractionFactor);

        /// @brief Computes a new curvilinear grid with line attracted towards \ref m_lines[0]
        /// @return The shifted curvilinear grid
        CurvilinearGrid Compute() override;

    private:
        std::shared_ptr<CurvilinearGrid> m_originalGrid; ///< A pointer to the original grid, used for computation of quantities
        double m_attractionFactor;                       ///< The attraction factor. If 0.0 no line will be attracted to m_lines[0]
    };
} // namespace meshkernel
