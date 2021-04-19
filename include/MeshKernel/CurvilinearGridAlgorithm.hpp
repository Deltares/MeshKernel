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
#include <MeshKernel/CurvilinearGridLine.hpp>
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

        /// @brief Executes the algorithm
        virtual CurvilinearGrid Compute() = 0;

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

        CurvilinearGrid m_grid;                    ///< A copy of the original grid, it will be modified by the algorithms
        std::vector<CurvilinearGridLine> m_lines;  ///< Selected grid lines
        CurvilinearGrid::NodeIndices m_lowerLeft;  ///< The lower left corner of the block
        CurvilinearGrid::NodeIndices m_upperRight; ///< The upper right corner of the block
    };
} // namespace meshkernel
