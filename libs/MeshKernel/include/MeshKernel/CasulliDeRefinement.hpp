//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2024.
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
#include "MeshKernel/UndoActions/UndoAction.hpp"

namespace meshkernel
{
    /// @brief Compute the Casulli de-refinement for a mesh.
    class CasulliDeRefinement
    {
    public:
        /// @brief Compute the Casulli de-refinement for the entire mesh.
        ///
        /// @param [in, out] mesh Mesh to be de-refined
        /* [[nodiscard]] */ static std::unique_ptr<meshkernel::UndoAction> Compute(Mesh2D& mesh);

        /// @brief Compute the Casulli de-refinement for the part of the mesh inside the polygon
        ///
        /// @param [in, out] mesh Mesh to be de-refined
        /// @param [in] polygon Area within which the mesh will be de-refined
        /* [[nodiscard]] */ static std::unique_ptr<meshkernel::UndoAction> Compute(Mesh2D& mesh, const Polygons& polygon);

    private:
        // WTF
        enum class ElementMask
        {
            A = 1,         //< front, 'A' cell (used to be node, delete it):  1
            B = 2,         //< front, 'B' cell (used to be link, keep it):    2
            C = 3,         //< 'C' cell (used to be cell, keep it):           3
            NotA = -1,     //< not in front, 'A' cell:                       -1
            NotB = -2,     //< not in front, 'B' cell:                       -2
            Unassigned = 0 //<                                                0
        };

        // TODO how much of this is shared with the Casulli refinement?

        ///@brief Indicates status of a node.
        ///
        /// \note the order of the enum values is very important to the working of the refinement algorithm
        enum class NodeMask : char
        {
            NewAssignedNode, //< a new node has been added, current node mask value is any value strictly greater than Unassigned.
            NewGeneralNode,  //< a new node has been added, current node mask value is any assigned value
            Unassigned,      //< Uninitialised state for node mask
            RegisteredNode,  //< Node is to be considered as part of the Casulli refinement/de-refinement
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
        static void InitialiseBoundaryNodes(const Mesh2D& mesh, std::vector<NodeMask>& nodeMask);

        /// @brief Initialise node mask for corner nodes
        ///
        /// @param [in] mesh The Mesh
        /// @param [in, out] nodeMask Node mask information
        static void InitialiseCornerNodes(const Mesh2D& mesh, std::vector<NodeMask>& nodeMask);

        /// @brief Initialise node mask for face nodes
        ///
        /// @param [in] mesh The Mesh
        /// @param [in, out] nodeMask Node mask information
        static void InitialiseFaceNodes(const Mesh2D& mesh, std::vector<NodeMask>& nodeMask);

        /// @brief Initialise the node mask array.
        ///
        /// @param [in] mesh Mesh used to initialise the node mask
        /// @param [in] polygon Only nodes inside the polygon are to be considered
        static std::vector<NodeMask> InitialiseNodeMask(const Mesh2D& mesh, const Polygons& polygon);

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

        static bool ElementIsSeed(const Mesh2D& mesh, const Polygons& polygon, const UInt face);

        static UInt FindElementSeedIndex(const Mesh2D& mesh, const Polygons& polygon);

        static void FindSurroundingCells(const Mesh2D& mesh,
                                         const Polygons& polygon [[maybe_unused]],
                                         const UInt kCell,
                                         const UInt nMax,
                                         UInt& nDirect, UInt& nIndirect,
                                         std::vector<UInt>& kDirect,
                                         std::vector<UInt>& kIndirect,
                                         std::vector<std::array<UInt, 2>>& kne);

        static void DoDeRefinement(const Mesh2D& mesh, const Polygons& polygon);

        static void UpdateFrontList(const Mesh& mesh, const std::vector<UInt>& frontList, std::vector<UInt>& frontListCopy, const UInt kNew);

    };

} // namespace meshkernel
