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

#include <MeshKernel/CurvilinearGrid/CurvilinearGridLine.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridNodeIndices.hpp>
#include <stddef.h>

namespace meshkernel
{

    /// @brief A struct describing a grid line in the curvilinear grid in terms of node indices
    struct CurvilinearGridLine
    {
        /// @brief The type of grid line, if it is in m or n direction
        enum class GridLineDirection
        {
            MDirection,
            NDirection
        };

        /// @brief CurvilinearGridLine constructor
        /// @param[in] startNode The start node of the grid line
        /// @param[in] endNode The end node of the grid line
        CurvilinearGridLine(CurvilinearGridNodeIndices const& startNode, CurvilinearGridNodeIndices const& endNode);

        /// @brief Inquire if a node is on a grid line
        /// @param[in] node The node to inquire
        /// @return True if the node belongs to the grid line, false otherwise
        [[nodiscard]] bool IsNodeOnLine(CurvilinearGridNodeIndices const& node) const;

        /// @brief Gets the indices of a node on the grid line
        /// @param[in] coordinate The one-dimensional coordinate along the grid line
        /// @return The node indices
        [[nodiscard]] CurvilinearGridNodeIndices GetNodeIndexFromCoordinate(size_t const& coordinate) const;

        /// @brief Inquires if the grid line is an M grid line
        /// @return True if it is an M grid line
        bool IsMGridLine() const { return m_gridLineType == GridLineDirection::MDirection; };

        /// @brief Inquires if the grid line is an N grid line
        /// @return True if it is an N grid line
        bool IsNGridLine() const { return m_gridLineType == GridLineDirection::NDirection; };

        CurvilinearGridNodeIndices m_startNode; ///<The start node of the grid line
        CurvilinearGridNodeIndices m_endNode;   ///<The end node of the grid line
        size_t m_startCoordinate;               ///<The start coordinate. If it is an MDirection, the start m otherwise the start n
        size_t m_endCoordinate;                 ///<The end coordinate. If it is an MDirection, the end m otherwise the end n
        size_t m_constantCoordinate;            ///<The constant coordinate. If it is an MDirection, the n coordinate, otherwise the m coordinate
        GridLineDirection m_gridLineType;       ///<The grid line type
    };
} // namespace meshkernel
