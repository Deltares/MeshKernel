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
        ///
        /// @param [in, out] mesh Mesh to be refined
        static void Compute(Mesh2D& mesh);

        /// @brief Compute the Casulli refinement for the part of the mesh inside the polygon
        ///
        /// @param [in, out] mesh Mesh to be refined
        /// @param [in] polygon Area within which the mesh will be refined
        static void Compute(Mesh2D& mesh, const Polygons& polygon);

    private:
        ///@brief Indicates status of a node.
        ///
        /// \note the order of the enum values is very important to the working of the refinement algorithm
        enum class NodeMask : char
        {
            NewAssignedNode, //< a new node has been added, current node mask value is any value strictly greater than Unassigned.
            NewGeneralNode,  //< a new node has been added, current node mask value is any assigned value
            Unassigned,      //< Uninitialised state for node mask
            RegisteredNode,  //< Node is to be considered as part of the Casulli refinement
            BoundaryNode,    //< Node lies on the boundary
            CornerNode       //< Node lies at corner of element on the boundary.
        };

        /// @brief Initial size of the edge array
        static constexpr UInt InitialEdgeArraySize = 100;

        /// @brief The maximum number of nodes that a newly created element can have.
        static constexpr UInt MaximumNumberOfNodesInNewlyCreatedElements = 4;

        /// @brief Simple bounded array of 4 UInt values
        /// @typedef EdgeNodes
        using EdgeNodes = std::array<UInt, 4>;

        /// @brief Initialise node mask for boundary nodes
        ///
        /// @param [in] mesh The Mesh
        /// @param [in, out] nodeMask Node mask information
        static void InitialiseBoundaryNodes (const Mesh2D& mesh, std::vector<NodeMask>& nodeMask);

        /// @brief Initialise node mask for corner nodes
        ///
        /// @param [in] mesh The Mesh
        /// @param [in, out] nodeMask Node mask information
        static void InitialiseCornerNodes (const Mesh2D& mesh, std::vector<NodeMask>& nodeMask);

        /// @brief Initialise node mask for face nodes
        ///
        /// @param [in] mesh The Mesh
        /// @param [in, out] nodeMask Node mask information
        static void InitialiseFaceNodes (const Mesh2D& mesh, std::vector<NodeMask>& nodeMask);

        /// @brief Initialise the node mask array.
        ///
        /// @param [in] mesh Mesh used to initialise the node mask
        /// @param [in] polygon Only nodes inside the polygon are to be considered
        static std::vector<NodeMask> InitialiseNodeMask(const Mesh2D& mesh, const Polygons& polygon);

        /// @brief Create any new nodes that need to be generated.
        ///
        /// @param [in, out] mesh The Mesh
        /// @param [in, out] newNodes List of new nodes and connectivity
        /// @param [in, out] nodeMask Node mask information
        static void ComputeNewNodes(Mesh2D& mesh, std::vector<EdgeNodes>& newNodes, std::vector<NodeMask>& nodeMask);

        /// @brief Connect newly generated nodes
        ///
        /// @param [in, out] mesh The mesh being refined
        /// @param [in] newNodes List of new nodes and connectivity
        /// @param [in] numEdges Number of edges in original mesh, before refinement.
        static void ConnectNodes(Mesh2D& mesh, const std::vector<EdgeNodes>& newNodes, const UInt numEdges);

        /// @brief Compute new edges required for refinement
        ///
        /// @param [in, out] mesh The mesh being refined
        /// @param [in] numNodes Number of nodes in original mesh, before refinement.
        /// @param [in] newNodes List of new nodes and connectivity
        /// @param [in, out] nodeMask Node mask information
        static void CreateMissingBoundaryEdges(Mesh2D& mesh, const UInt numNodes, const std::vector<EdgeNodes>& newNodes, std::vector<NodeMask>& nodeMask);

        /// @brief Compute new nodes on faces
        ///
        /// @param [in, out] mesh The mesh being refined
        /// @param [in, out] newNodes List of new nodes and connectivity
        /// @param [in, out] nodeMask Node mask information
        static void ComputeNewFaceNodes(Mesh2D& mesh, std::vector<EdgeNodes>& newNodes, std::vector<NodeMask>& nodeMask);

        /// @brief Compute new nodes on edges
        ///
        /// @param [in, out] mesh The mesh being refined
        /// @param [in] numEdges Number of edges in original mesh, before refinement.
        /// @param [in, out] newNodes List of new nodes and connectivity
        /// @param [in, out] nodeMask Node mask information
        static void ComputeNewEdgeNodes(Mesh2D& mesh, const UInt numEdges, std::vector<EdgeNodes>& newNodes, std::vector<NodeMask>& nodeMask);

        /// @brief Connect edges
        ///
        /// @param [in, out] mesh The mesh being refined
        /// @param [in] currentNode The node being connected
        /// @param [in] newNodes List of new nodes and connectivity
        /// @param [out] edgeCount Number of edges created
        /// @param [in, out] newEdges Identifiers of edges created
        static void ConnectEdges(Mesh2D& mesh,
                                 const UInt currentNode,
                                 const std::vector<EdgeNodes>& newNodes,
                                 UInt& edgeCount,
                                 std::vector<UInt>& newEdges);

        /// @brief Connect face node
        ///
        /// @param [in, out] mesh The mesh being refined
        /// @param [in] currentFace The face being connected
        /// @param [in, out] newNodes List of new nodes and connectivity
        /// @param [in, out] nodeMask Node mask information
        static void ConnectFaceNodes(Mesh2D& mesh,
                                     const UInt currentFace,
                                     const std::vector<EdgeNodes>& newNodes,
                                     std::vector<NodeMask>& nodeMask);

        /// @brief Connect newly generated nodes
        ///
        /// @param [in, out] mesh The mesh being refined
        /// @param [in] newNodes List of new nodes and connectivity
        /// @param [in] numNodes Number of nodes in original mesh, before refinement.
        /// @param [in] numEdges Number of edges in original mesh, before refinement.
        /// @param [in] numFaces Number of faces in original mesh, before refinement.
        /// @param [in, out] nodeMask Node mask information
        static void ConnectNewNodes(Mesh2D& mesh, const std::vector<EdgeNodes>& newNodes, const UInt numNodes, const UInt numEdges, const UInt numFaces, std::vector<NodeMask>& nodeMask);

        /// @brief Add newly generated nodes to the newNodes list.
        ///
        /// @param [in] mesh The Mesh
        /// @param [in] nodeId Node shared by edge1Index and edge2Index
        /// @param [in] edge1Index First of two edges used to determine face
        /// @param [in] edge2Index Second of two edges used to determine face
        /// @param [in] newNodeId Identifier of the new node
        /// @param [in, out] newNodes List of new nodes and connectivity
        static void StoreNewNode(const Mesh2D& mesh, const UInt nodeId, const UInt edge1Index, const UInt edge2Index, const UInt newNodeId, std::vector<EdgeNodes>& newNodes);

        /// @brief Find elements and nodes that form the patch of elements directly connected to the node.
        ///
        /// @param [in] mesh The mesh
        /// @param [in] currentNode Node identifier
        /// @param [out] sharedFaces Identifiers of all faces directly connected to the node
        /// @param [out] connectedNodes Identifiers of the nodes of all faces connected to the node
        /// @param [out] faceNodeMapping Mapping from node index to the position in connectedNodes list.
        static void FindPatchIds(const Mesh2D& mesh,
                                 const UInt currentNode,
                                 std::vector<UInt>& sharedFaces,
                                 std::vector<UInt>& connectedNodes,
                                 std::vector<std::vector<UInt>>& faceNodeMapping);

        /// @brief Delete any unused nodes and perform a mesh administration.
        ///
        /// @param [in, out] mesh The mesh
        /// @param [in] numNodes The original number of nodes in the mesh
        /// @param [in] nodeMask Node mask information
        static void Administrate(Mesh2D& mesh, const UInt numNodes, const std::vector<NodeMask>& nodeMask);
    };

} // namespace meshkernel
