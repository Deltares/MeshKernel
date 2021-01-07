//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2020.
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
#include "Entities.hpp"

#include <vector>

/// \namespace meshkernel
/// @brief Contains the logic of the C++ static library
namespace meshkernel
{
    /// @brief A class describing an unstructured 1d mesh
    class Mesh1D
    {
    public:
        /// @brief Default constructor
        Mesh1D() = default;

        /// @brief Construct a mesh1d starting from the edges and nodes
        /// @param[in] edges The input edges
        /// @param[in] nodes The input nodes
        /// @param[in] projection  The projection to use
        Mesh1D(const std::vector<Edge>& edges,
               const std::vector<Point>& nodes,
               Projection projection);

        /// @brief Inquire if a node is on boundary
        /// @param[in] node The node index
        /// @return If the node is on boundary
        [[nodiscard]] bool IsNodeOnBoundary(size_t node) const { return m_nodesNumEdges[node] == 1; }

        /// @brief Get the number of valid edges
        /// @return The number of valid edges
        [[nodiscard]] auto GetNumEdges() const { return m_numEdges; }

        /// @brief constructs the mesh1d mapping (findcells)
        void FindFaces();

        // nodes
        std::vector<Point> m_nodes;
        std::vector<size_t> m_nodesNumEdges;           ///< For each node, the number of connected edges (nmk)
        std::vector<std::vector<size_t>> m_nodesEdges; ///< For each node, the indices of connected edges (nod%lin)

        // edges
        std::vector<Edge> m_edges;

        // faces
        std::vector<std::vector<int>> m_facesNodes;
        std::vector<Point> m_facesMassCenters;

        // the projection
        Projection m_projection;

    private:
        // counters
        size_t m_numFaces = 0; ///< Number of valid faces (nump)
        size_t m_numNodes = 0; ///< Number of valid nodes in m_nodes
        size_t m_numEdges = 0; ///< Number of valid edges in m_edges
    };
} // namespace meshkernel
