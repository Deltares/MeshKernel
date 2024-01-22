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

#include <array>
#include <vector>

#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/Polygons.hpp"

namespace meshkernel
{
    /// @brief Compute the Casulli refinement for a mesh.
    class CasulliRefinement
    {
    public:
        /// @brief Compute the Casulli refinement for the entire mesh.
        static void Compute(Mesh2D& mesh);

        /// @brief Compute the Casulli refinement for the part of the mesh inside the polygon
        static void Compute(Mesh2D& mesh, const Polygons& polygon);

    private:
        /// @brief Initial size of the edge array
        static constexpr UInt InitialEdgeArraySize = 100;

        /// @brief The maximum number of nodes that a newly created element can have.
        static constexpr UInt MaximumNumberOfNodesInNewlyCreatedElements = 4;

        /// @brief Simple bounded array of 4 UInt values
        /// @typedef LinkNodes
        using LinkNodes = std::array<UInt, 4>;

        /// @brief Initialise the node mask array.
        static std::vector<int> InitialiseNodeMask(const Mesh2D& mesh, const Polygons& polygon);

        /// @brief Create any new nodes that need to be generated.
        static void ComputeNewNodes(Mesh2D& mesh, std::vector<LinkNodes>& newNodes, std::vector<int>& nodeMask);

        /// @brief Connect newly generated nodes
        ///
        /// @param [in, out] mesh The mesh being refined
        /// @param [in] newNodes List of new nodes and connectivity
        /// @param [in] numNodes Number of nodes in original mesh, before refinement.
        /// @param [in] numEdges Number of edges in original mesh, before refinement.
        /// @param [in] numFaces Number of faces in original mesh, before refinement.
        /// @param [in, out] nodeMask Node mask information
        static void LinkNewNodes(Mesh2D& mesh, const std::vector<LinkNodes>& newNodes, const UInt numNodes, const UInt numEdges, const UInt numFaces, std::vector<int>& nodeMask);

        /// @brief Add newly generated nodes to the newNodes list.
        static void StoreNewNode(const Mesh2D& mesh, const UInt nodeId, const UInt link1Index, const UInt link2Index, const UInt newNodeId, std::vector<LinkNodes>& newNodes);

        /// @brief Find elements and nodes that form the patch of elements directly connected to the node.
        static void FindPatchIds(const Mesh2D& mesh,
                                 const UInt currentNode,
                                 std::vector<UInt>& sharedFaces,
                                 std::vector<UInt>& connectedNodes,
                                 std::vector<std::vector<UInt>>& faceNodeMapping);

        /// @brief Delete any unused nodes and performan mesh administration.
        static void Administrate(Mesh2D& mesh, const UInt numNodes, const std::vector<int>& nodeMask);
    };

} // namespace meshkernel
