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
    /// @brief A class describing an unstructured mesh (can be both 1d or 2d)
    class Mesh
    {
    public:
        /// @brief Default constructor
        Mesh() = default;

        /// @brief Construct a mesh starting from the edges and nodes
        /// @param[in] edges The input edges
        /// @param[in] nodes The input nodes
        /// @param[in] projection  The projection to use
        Mesh(const std::vector<Edge>& edges,
             const std::vector<Point>& nodes,
             Projection projection);

        /// @brief Inquire if a node is on boundary
        /// @param[in] node The node index
        /// @return If the node is on boundary
        [[nodiscard]] bool IsNodeOnBoundary(size_t node) const { return m_nodesNumEdges[node] == 1; }

        /// @brief Get the number of valid nodes
        /// @return The number of valid node
        [[nodiscard]] auto GetNumNodes() const { return m_numNodes; }

        /// @brief Get the number of valid edges
        /// @return The number of valid edges
        [[nodiscard]] auto GetNumEdges() const { return m_numEdges; }

        /// @brief Get the number of valid faces
        /// @return The number of valid faces
        [[nodiscard]] auto GetNumFaces() const { return m_numFaces; }

        /// @brief Get the number of edges for a face
        /// @param[in] faceIndex The face index
        /// @return The number of valid faces
        [[nodiscard]] auto GetNumFaceEdges(size_t faceIndex) const { return m_numFacesNodes[faceIndex]; }

        /// @brief Get the number of faces an edges shares
        /// @param[in] edgeIndex The edge index
        /// @return The number of faces an edges shares
        [[nodiscard]] auto GetNumEdgesFaces(size_t edgeIndex) const { return m_edgesNumFaces[edgeIndex]; }

        /// @brief Node administration (setnodadmin)
        void NodeAdministration();

        /// @brief Removes all invalid nodes and edges
        void DeleteInvalidNodesAndEdges();

        // nodes
        std::vector<Point> m_nodes;                    ///< The mesh nodes (xk, yk)
        std::vector<std::vector<size_t>> m_nodesEdges; ///< For each node, the indices of connected edges (nod%lin)
        std::vector<size_t> m_nodesNumEdges;           ///< For each node, the number of connected edges (nmk)
        std::vector<int> m_nodeMask;                   ///< The node mask (kc)
        std::vector<std::vector<size_t>> m_nodesNodes; ///< For each node, its neighbours
        std::vector<int> m_nodesTypes;                 ///< The node types (nb)

        // edges
        std::vector<Edge> m_edges;                     ///< The edges, defined as first and second node(kn)
        std::vector<std::vector<size_t>> m_edgesFaces; ///< For each edge, the shared face index (lne)
        std::vector<size_t> m_edgesNumFaces;           ///< For each edge, the number of shared faces(lnn)
        std::vector<double> m_edgeLengths;             ///< The edge lengths
        std::vector<int> m_edgeMask;                   ///< The edge mask (lc)
        std::vector<Point> m_edgesCenters;             ///< The edges centers

        // faces
        std::vector<std::vector<size_t>> m_facesNodes; ///< The nodes composing the faces, in ccw order (netcell%Nod)
        std::vector<size_t> m_numFacesNodes;           ///< The number of nodes composing the face (netcell%N)
        std::vector<std::vector<size_t>> m_facesEdges; ///< The edge indices composing the face (netcell%lin)
        std::vector<Point> m_facesCircumcenters;       ///< The face circumcenters the face circumcenter (xz, yz)
        std::vector<Point> m_facesMassCenters;         ///< The faces centers of mass (xzw, yzw)
        std::vector<double> m_faceArea;                ///< The face area

        Projection m_projection; ///< The projection used

        // counters
        size_t m_numFaces = 0; ///< Number of valid faces (nump)
        size_t m_numNodes = 0; ///< Number of valid nodes in m_nodes
        size_t m_numEdges = 0; ///< Number of valid edges in m_edges
    };
} // namespace meshkernel
