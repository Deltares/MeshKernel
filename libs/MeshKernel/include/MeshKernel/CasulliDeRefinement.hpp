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
        std::unique_ptr<meshkernel::UndoAction> Compute(Mesh2D& mesh);

        /// @brief Compute the Casulli de-refinement for the part of the mesh inside the polygon
        ///
        /// @param [in, out] mesh Mesh to be de-refined
        /// @param [in] polygon Area within which the mesh will be de-refined
        std::unique_ptr<meshkernel::UndoAction> Compute(Mesh2D& mesh, const Polygons& polygon);

        /// @brief Compute centres of elements to be deleted.
        ///
        /// Requires that the element mass-centres are pre-computed.
        std::vector<Point> ElementsToDelete(const Mesh2D& mesh, const Polygons& polygon);

    private:
        /// @brief Element mask values used in the de-refinement procedure.
        ///
        /// Enumeration comments are from the original Fortran code.
        enum class ElementMask
        {
            A = 1,         //< front, 'A' cell (used to be node, delete it):  1
            B = 2,         //< front, 'B' cell (used to be link, keep it):    2
            C = 3,         //< 'C' cell (used to be cell, keep it):           3
            NotA = -1,     //< not in front, 'A' cell:                       -1
            NotB = -2,     //< not in front, 'B' cell:                       -2
            Unassigned = 0 //< not assigned a value                           0
        };

        /// @brief Indicate if the element can be a seed element or not.
        bool ElementIsSeed(const Mesh2D& mesh,
                           const std::vector<int>& nodeTypes,
                           const Polygons& polygon,
                           const UInt element);

        /// @brief Find the seed element id to start the mesh de-refinement.
        UInt FindElementSeedIndex(const Mesh2D& mesh,
                                  const std::vector<int>& nodeTypes,
                                  const Polygons& polygon);

        /// @brief Find all elements that are connected along edges to elementId.
        void FindDirectlyConnectedCells(const Mesh2D& mesh,
                                        const UInt elementId,
                                        std::vector<UInt>& connected);

        /// @brief Find all elements that are connected by nodes to elementId.
        void FindIndirectlyConnectedCells(const Mesh2D& mesh,
                                          const UInt elementId,
                                          const std::vector<UInt>& directlyConnected,
                                          std::vector<UInt>& indirectlyConnected);

        /// @brief Find element id's
        void FindAdjacentCells(const Mesh2D& mesh,
                               const std::vector<UInt>& directlyConnected,
                               const std::vector<UInt>& indirectlyConnected,
                               std::vector<std::array<int, 2>>& kne);

        /// @brief Find the elements that are connected to the elementId.
        void FindSurroundingCells(const Mesh2D& mesh,
                                  const Polygons& polygon [[maybe_unused]],
                                  const UInt elementId,
                                  std::vector<UInt>& directlyConnected,
                                  std::vector<UInt>& indirectlyConnected,
                                  std::vector<std::array<int, 2>>& kne);

        /// @brief Initialise the element mask.
        std::vector<ElementMask> InitialiseElementMask(const Mesh2D& mesh,
                                                       const std::vector<int>& nodeTypes,
                                                       const Polygons& polygon);

        /// \brief Determine if the element can be deleted from the mesh or not.
        bool ElementCannotBeDeleted(const Mesh2D& mesh,
                                    const std::vector<int>& nodeTypes,
                                    const Polygons& polygon,
                                    const UInt elementId);

        /// @brief Compute coordinates of the new node.
        Point ComputeNewNodeCoordinates(const Mesh2D& mesh,
                                        const std::vector<int>& nodeTypes,
                                        const UInt nodeId);

        /// @brief
        void UpdateDirectlyConnectedElements(Mesh2D& mesh,
                                             const UInt elementId,
                                             const std::vector<UInt>& directlyConnected,
                                             const std::vector<std::array<int, 2>>& kne);

        /// @brief Get the most significant node type for all nodes of the element.
        int GetNodeCode(const Mesh2D& mesh,
                        const std::vector<int>& nodeTypes,
                        const UInt elementId);

        /// @brief Add element id to the list of id's
        ///
        /// only added is it is not already on the list and the element is a quadrilateral
        void AddElementToList(const Mesh& mesh, const std::vector<UInt>& frontList, std::vector<UInt>& frontListCopy, const UInt kNew);

        /// @brief Redirect nodes of connected cells, deactivate polygons of degree smaller than three
        void RedirectNodesOfConnectedElements(Mesh2D& mesh, const UInt elementId, const UInt nodeId, const std::vector<UInt>& indirectlyConnected);

        /// @brief Removes nodes from the boundary that will not be part of the de-refined mesh.
        void RemoveUnwantedBoundaryNodes(Mesh2D& mesh,
                                         const std::vector<int>& nodeTypes,
                                         const Polygons& polygon,
                                         const std::vector<UInt>& indirectlyConnected);

        /// @brief Delete an element
        void DeleteElement(Mesh2D& mesh,
                           std::vector<int>& nodeTypes,
                           const Polygons& polygon,
                           const UInt elementId,
                           const std::vector<UInt>& directlyConnected,
                           const std::vector<UInt>& indirectlyConnected,
                           const std::vector<std::array<int, 2>>& kne);

        /// @brief Clean up the edge
        void CleanUpEdge(Mesh2D& mesh, const UInt edgeId);

        /// @brief Find the id of the shared node for two edges.
        ///
        /// If no node is shared then return null value
        UInt FindCommonNode(const Mesh2D& mesh, const UInt edgeId1, const UInt edgeId2);

        /// @brief Determine if the element is convex
        bool IsElementConvex(const Mesh2D& mesh, const UInt cell);

        /// @brief Do the Casullu de-refinement
        void DoDeRefinement(Mesh2D& mesh, const Polygons& polygon);

        /// @brief Compute the mesh node types.
        ///
        /// Uses the m_nodeTypes has been generated in the mesh.
        std::vector<int> ComputeNodeTypes(const Mesh2D& mesh, const Polygons& polygon);
    };

} // namespace meshkernel
