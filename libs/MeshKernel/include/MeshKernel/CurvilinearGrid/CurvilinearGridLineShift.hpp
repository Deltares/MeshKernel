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
#include <MeshKernel/Entities.hpp>

namespace meshkernel
{

    /// @brief A class implementing the curvilinear grid line shift. Line is shifted on the new locations.
    /// and the displacement gets distributed on the influence zone, set with SetBlock
    ///
    /// This option provides the possibility to fit the curvilinear grid’s edges to a land boundary.
    class CurvilinearGridLineShift : public CurvilinearGridAlgorithm
    {
    public:
        /// @brief Class constructor
        /// @param[in] grid The input curvilinear grid
        CurvilinearGridLineShift(CurvilinearGrid& grid);

        /// @brief Computes a new curvilinear grid with the line shift
        void Compute() override;

        /// @brief Moves a node from one position to another
        /// @param[in] fromPoint The input position, the closest node on the \ref m_grid grid will be used
        /// @param[in] toPoint The coordinates of the new position
        void MoveNode(Point const& fromPoint, Point const& toPoint);

    private:
        /// @brief Distribute the displacement around the node on the influence zone.
        /// @param[in] node The node to account for. The displacement around this not is calculated subtracting \ref m_grid to \ref m_originalGrid
        void TransformGrid(CurvilinearGridNodeIndices const& node);

        /// @brief Transform the displacement around a node to local or global (TOLOCL)
        /// @param[in] displacement The displacement to transform.
        /// @param[in] node The node position
        /// @param[in] toLocal A boolean to indicate whatever to transform the displacement to local grid (True) or to global grid (false)
        /// @return The new displacement
        Point TransformDisplacement(Point const& displacement, CurvilinearGridNodeIndices const& node, bool toLocal) const;

        CurvilinearGrid m_originalGrid; ///< A pointer to the original grid
    };
} // namespace meshkernel
