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

#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridLine.hpp"
#include "MeshKernel/Splines.hpp"
#include "MeshKernel/UndoActions/UndoAction.hpp"

namespace meshkernel
{
    class CurvilinearGrid;

    /// @brief A class implementing the curvilinear grid orthogonalization algorithm
    class CurvilinearGridAlgorithm
    {
    public:
        /// @brief Class constructor
        /// @param[in] grid                        The input curvilinear grid
        explicit CurvilinearGridAlgorithm(CurvilinearGrid& grid);

        /// @brief Executes the algorithm
        virtual UndoActionPtr Compute() = 0;

        /// @brief Sets a block where the algorithm should execute
        /// @param[in] firstCornerPoint            The first point defining the bounding box
        /// @param[in] secondCornerPoint           The second point defining the bounding box
        void SetBlock(Point const& firstCornerPoint, Point const& secondCornerPoint);

        /// @brief Sets a block using curvilinear node indices
        /// @param[in] firstCornerPoint            The curvilinear grid index defining the lower left corner of the bounding box
        /// @param[in] secondCornerPoint           The curvilinear grid index defining the upper right corner of the bounding box
        void SetBlock(CurvilinearGridNodeIndices const& firstCornerPoint, CurvilinearGridNodeIndices const& secondCornerPoint);

        /// @brief Sets a line, applicable to algorithms such as directional smoothing,
        /// freeze a line in curvilinear grid orthogonalization and line shift
        /// @param[in] firstPoint The point containing the first point of the line to freeze
        /// @param[in] secondPoint The point containing the second point of the line to freeze
        void SetLine(Point const& firstPoint, Point const& secondPoint);

        /// @brief Fills the frozen nodes mask using the frozen lines
        void ComputeFrozenNodes();

        /// @brief Virtual destructor
        virtual ~CurvilinearGridAlgorithm() = default;

        CurvilinearGrid& m_grid;                  ///< A reference of the grid, modified by the algorithms
        std::vector<CurvilinearGridLine> m_lines; ///< Selected frozen grid lines (need to change)
        CurvilinearGridNodeIndices m_lowerLeft;   ///< The lower left corner of the block
        CurvilinearGridNodeIndices m_upperRight;  ///< The upper right corner of the block
        lin_alg::Matrix<bool> m_isGridNodeFrozen; ///< A mask for setting some of the grid nodes frozen
    };
} // namespace meshkernel
