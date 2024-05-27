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
        [[nodiscard]] std::unique_ptr<UndoAction> Compute(Mesh2D& mesh) const;

        /// @brief Compute the Casulli de-refinement for the part of the mesh inside the polygon
        ///
        /// @param [in, out] mesh Mesh to be de-refined
        /// @param [in] polygon Area within which the mesh will be de-refined
        [[nodiscard]] std::unique_ptr<UndoAction> Compute(Mesh2D& mesh, const Polygons& polygon) const;

        /// @brief Compute centres of elements to be deleted.
        ///
        /// Requires that the element mass-centres are pre-computed.
        std::vector<Point> ElementsToDelete(const Mesh2D& mesh, const Polygons& polygon) const;

    private:
        /// @brief Maximum number of iterations allowed when initialising the element mask
        static const UInt maxIterationCount = 1000;

        /// @brief Maximum reserve size for arrays used in the de-refinement
        static const UInt maximumSize = 1000;

        /// @brief Label for the element in the de-refined grid.
        ///
        /// The prefix, e.g. WasNode, indicates the original mesh entity around/from which the de-refined element was constructed.
        /// First and after suffix, indicates the order in which the elements are to be processed.
        /// Enumeration values and comments are from the original Fortran code.
        enum class ElementType
        {
            WasNodeFirst = 1,  //< front, 'A' cell (used to be node, delete it):  1
            WasEdgeFirst = 2,  //< front, 'B' cell (used to be link, keep it):    2
            WasCell = 3,       //< 'C' cell (used to be cell, keep it):           3
            WasNodeAfter = -1, //< not in front, 'A' cell:                       -1
            WasEdgeAfter = -2, //< not in front, 'B' cell:                       -2
            Unassigned = 0     //< not assigned a value                           0
        };

        /// @brief Indicate if the element can be a seed element or not.
        bool ElementIsSeed(const Mesh2D& mesh,
                           const std::vector<int>& nodeTypes,
                           const UInt element) const;

        /// @brief Find the seed element id to start the mesh de-refinement.
        UInt FindElementSeedIndex(const Mesh2D& mesh,
                                  const std::vector<int>& nodeTypes) const;

        /// @brief Find all elements that are connected along edges to elementId.
        void FindDirectlyConnectedCells(const Mesh2D& mesh,
                                        const UInt elementId,
                                        std::vector<UInt>& connected) const;

        /// @brief Find all elements that are connected by nodes to elementId.
        void FindIndirectlyConnectedCells(const Mesh2D& mesh,
                                          const UInt elementId,
                                          const std::vector<UInt>& directlyConnected,
                                          std::vector<UInt>& indirectlyConnected) const;

        /// @brief Find element id's
        void FindAdjacentCells(const Mesh2D& mesh,
                               const std::vector<UInt>& directlyConnected,
                               const std::vector<UInt>& indirectlyConnected,
                               std::vector<std::array<int, 2>>& kne) const;

        /// @brief Find the elements that are connected to the elementId.
        void FindSurroundingCells(const Mesh2D& mesh,
                                  const UInt elementId,
                                  std::vector<UInt>& directlyConnected,
                                  std::vector<UInt>& indirectlyConnected,
                                  std::vector<std::array<int, 2>>& kne) const;

        /// @brief Initialise the element mask.
        std::vector<ElementType> InitialiseElementType(const Mesh2D& mesh,
                                                       const std::vector<int>& nodeTypes) const;

        /// \brief Determine if the element can be deleted from the mesh or not.
        bool ElementCannotBeDeleted(const Mesh2D& mesh,
                                    const std::vector<int>& nodeTypes,
                                    const Polygons& polygon,
                                    const UInt elementId) const;

        /// @brief Compute coordinates of the new node.
        Point ComputeNewNodeCoordinates(const Mesh2D& mesh,
                                        const std::vector<int>& nodeTypes,
                                        const UInt nodeId) const;

        /// @brief Update the mesh members for the mesh description and connectivity.
        [[nodiscard]] bool UpdateDirectlyConnectedElements(Mesh2D& mesh,
                                                           const UInt elementId,
                                                           const std::vector<UInt>& directlyConnected,
                                                           const std::vector<std::array<int, 2>>& kne) const;

        /// @brief Update the mesh members for the mesh description and connectivity for triangle elements
        bool UpdateDirectlyConnectedTriangleElements(Mesh2D& mesh,
                                                     const UInt index,
                                                     const UInt connectedElementId,
                                                     const std::vector<std::array<int, 2>>& kne) const;

        /// @brief Update the mesh members for the mesh description and connectivity for non-triangle elements
        ///
        /// That is, element with 4 or more edges.
        void UpdateDirectlyConnectedNonTriangleElements(Mesh2D& mesh,
                                                        const UInt index,
                                                        const UInt elementId,
                                                        const UInt connectedElementId) const;

        /// @brief Get the most significant node type for all nodes of the element.
        int GetNodeCode(const Mesh2D& mesh,
                        const std::vector<int>& nodeTypes,
                        const UInt elementId) const;

        /// @brief Add element id to the list of id's
        ///
        /// only added is it is not already on the list and the element is a quadrilateral
        void AddElementToList(const Mesh& mesh, const std::vector<UInt>& frontList, std::vector<UInt>& frontListCopy, const UInt elementId) const;

        /// @brief Redirect nodes of connected cells, deactivate polygons of degree smaller than three
        void RedirectNodesOfConnectedElements(Mesh2D& mesh, const UInt elementId, const UInt nodeId, const std::vector<UInt>& indirectlyConnected) const;

        /// @brief Removes nodes from the boundary that will not be part of the de-refined mesh.
        [[nodiscard]] bool RemoveUnwantedBoundaryNodes(Mesh2D& mesh,
                                                       const std::vector<int>& nodeTypes,
                                                       const Polygons& polygon,
                                                       const std::vector<UInt>& indirectlyConnected) const;

        /// @brief Delete an element
        [[nodiscard]] bool DeleteElement(Mesh2D& mesh,
                                         std::vector<int>& nodeTypes,
                                         const Polygons& polygon,
                                         const UInt elementId,
                                         const std::vector<UInt>& directlyConnected,
                                         const std::vector<UInt>& indirectlyConnected,
                                         const std::vector<std::array<int, 2>>& kne) const;

        /// @brief Clean up the edge
        ///
        /// @returns Indicates if the cleanp-up was successful or not
        [[nodiscard]] bool CleanUpEdge(Mesh2D& mesh, const UInt edgeId) const;

        /// @brief Do the Casullu de-refinement
        [[nodiscard]] bool DoDeRefinement(Mesh2D& mesh, const Polygons& polygon) const;

        /// @brief Compute the mesh node types.
        ///
        /// Uses the m_nodeTypes has been generated in the mesh.
        std::vector<int> ComputeNodeTypes(const Mesh2D& mesh, const Polygons& polygon) const;
    };

} // namespace meshkernel
