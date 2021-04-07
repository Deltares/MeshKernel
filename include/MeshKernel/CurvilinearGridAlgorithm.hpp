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

#include <MeshKernel/CurvilinearGrid.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Splines.hpp>

namespace meshkernel
{
    class CurvilinearGrid;

    /// @brief A class implementing the curvilinear grid orthogonalization algorithm
    class CurvilinearGridAlgorithm
    {
    public:
        /// @brief Class constructor
        /// @param[in] grid                        The input curvilinear grid
        CurvilinearGridAlgorithm(std::shared_ptr<CurvilinearGrid> grid);

        /// @brief Executes the algorithm, pure virtual
        virtual std::shared_ptr<CurvilinearGrid> Compute() = 0;

        /// @brief Sets a block where the algorithm should execute
        /// @param[in] firstCornerPoint            The first point defining the bounding box
        /// @param[in] secondCornerPoint           The second point defining the bounding box
        void SetBlock(Point const& firstCornerPoint, Point const& secondCornerPoint);

        /// @brief Sets a line, applicable to algorithms such as directional smoothing,
        /// freeze a line in curvilinear grid orthogonalization and line shift
        /// @param[in] firstPoint The point containing the first point of the line to freeze
        /// @param[in] secondPoint The point containing the second point of the line to freeze
        void SetLine(Point const& firstPoint, Point const& secondPoint);

        /// @brief Virtual destructor
        virtual ~CurvilinearGridAlgorithm() = default;

        struct GridLine
        {
            enum class GridLineType
            {
                MGridLine,
                NGridLine
            };

            GridLine(const CurvilinearGrid::NodeIndices& m_start_node, const CurvilinearGrid::NodeIndices& m_end_node)
                : m_startNode(m_start_node),
                  m_endNode(m_end_node)
            {
                // The start-end coordinates of the new line, along m_m or m_n direction
                m_gridLineType = m_startNode.m_m == m_endNode.m_m ? GridLineType::MGridLine : GridLineType::NGridLine;
                m_startCoordinate = m_gridLineType == GridLineType::MGridLine ? m_startNode.m_n : m_startNode.m_m;
                m_endCoordinate = m_gridLineType == GridLineType::MGridLine ? m_endNode.m_n : m_endNode.m_m;
                m_constantCoordinate = m_gridLineType == GridLineType::MGridLine ? m_startNode.m_m : m_startNode.m_n;
            }

            CurvilinearGrid::NodeIndices m_startNode;
            CurvilinearGrid::NodeIndices m_endNode;
            size_t m_startCoordinate;
            size_t m_endCoordinate;
            size_t m_constantCoordinate;
            GridLineType m_gridLineType;
        };

        std::shared_ptr<CurvilinearGrid> m_grid;   ///< A pointer to the curvilinear grid to modify
        std::vector<GridLine> m_lines;             ///< Selected grid lines
        CurvilinearGrid::NodeIndices m_lowerLeft;  ///< The lower left corner of the block
        CurvilinearGrid::NodeIndices m_upperRight; ///< The upper right corner of the block
    };
} // namespace meshkernel
