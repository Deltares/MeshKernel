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
#include <MeshKernel/RTree.hpp>

namespace meshkernel
{
    /// @brief A class representing a curvilinear grid
    class CurvilinearGrid
    {
    public:
        /// @brief Default constructor
        /// @returns
        CurvilinearGrid() = default;

        /// @brief Create a new (empty) curvilinear grid
        /// @param[in] m Number of columns (horizontal direction)
        /// @param[in] n Number of rows (vertical direction)
        CurvilinearGrid(size_t m, size_t n);

        /// @brief Creates a new curvilinear grid from a given set of points
        /// @param[in] grid       The input grid points
        /// @param[in] projection The projection to use
        CurvilinearGrid(const std::vector<std::vector<Point>>& grid, Projection projection);

        /// @brief Builds the node three to find nodes on the curvilinear grid
        void BuildTree();

        /// @brief Get the m and n indices of the closest cu
        std::tuple<int, int> GetNodeIndices(Point point);

        Projection m_projection;                 ///< The projection used
        std::vector<std::vector<Point>> m_nodes; ///< Member variable storing the grid
        RTree m_nodesRTree;                      ///< Spatial R-Tree used to inquire nodes
        size_t m_numM = 0;                       ///< The number of m coordinates (vertical lines)
        size_t m_numN = 0;                       ///< The number of n coordinates (horizontal lines)

    private:
        /// @brief A customized type used for searching the closest m and n indices closest to a point
        struct Node
        {
            double x; ///< The point x coordinate
            double y; ///< The point y coordinate
            size_t m; ///< The point m coordinate
            size_t n; ///< The point n coordinate
        };

        std::vector<Node> m_flattenNodes; ///< The flattened nodes
    };
} // namespace meshkernel
