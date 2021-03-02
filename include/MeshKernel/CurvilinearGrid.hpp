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
#include <MeshKernel/Mesh.hpp>
#include <MeshKernel/RTree.hpp>

namespace meshkernel
{
    /// @brief A class representing a curvilinear grid
    class CurvilinearGrid : public Mesh
    {
    public:
        /// @brief Default constructor
        /// @returns
        CurvilinearGrid() = default;

        /// @brief Creates a new curvilinear grid from a given set of points
        /// @param[in] grid       The input grid points
        /// @param[in] projection The projection to use
        CurvilinearGrid(std::vector<std::vector<Point>> grid, Projection projection);

        /// @brief Converting a curvilinear mesh to a set of nodes, edges and returns the original mapping (gridtonet)
        /// @returns The nodes, the edges, and the original mapping (m and n indices for each node)
        std::tuple<std::vector<Point>, std::vector<Edge>, std::vector<std::pair<size_t, size_t>>> ConvertCurvilinearToNodesAndEdges();

        /// @brief Set internal flat copies of nodes and edges, so the pointer to the first entry is communicated with the front-end
        void SetFlatCopies();

        /// @brief Builds the node three to find nodes on the curvilinear grid
        void BuildTree();

        /// @brief Get the m and n indices of the node closest to the point
        /// @param[in] point       The input grid points
        std::tuple<int, int> GetNodeIndices(Point point);

        size_t m_numM = 0;                           ///< The number of m coordinates (vertical lines)
        size_t m_numN = 0;                           ///< The number of n coordinates (horizontal lines)
        std::vector<std::vector<Point>> m_gridNodes; ///< Member variable storing the grid

    private:
        std::vector<std::pair<size_t, size_t>> m_gridIndices; ///< the original mapping of the flatten nodes in the curvilinear grid
    };
} // namespace meshkernel
