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
    /// @brief Compute the Casulli de-refinement for a mesh (derefine_mesh).
    class CasulliDeRefinement
    {
    public:
        /// @brief Compute the Casulli de-refinement for the entire mesh.
        ///
        /// @param [in, out] mesh Mesh to be de-refined
        [[nodiscard]] static std::unique_ptr<UndoAction> Compute(Mesh2D& mesh);

        /// @brief Compute the Casulli de-refinement for the part of the mesh inside the polygon
        ///
        /// @param [in, out] mesh Mesh to be de-refined
        /// @param [in] polygon Area within which the mesh will be de-refined
        [[nodiscard]] static std::unique_ptr<UndoAction> Compute(Mesh2D& mesh, const Polygons& polygon);

        /// @brief Compute centres of elements to be deleted.
        ///
        /// Requires that the element mass-centres are pre-computed.
        static std::vector<Point> ElementsToDelete(const Mesh2D& mesh, const Polygons& polygon);

    private:
        /// @brief Maximum number of iterations allowed when initialising the element mask
        static constexpr UInt maxIterationCount = 1000;

        /// @brief Maximum reserve size for arrays used in the de-refinement
        static constexpr UInt maximumSize = 1000;

        /// @brief Label for the element in the de-refined grid.
        ///
        /// The prefix, e.g. WasNode, indicates the original mesh entity around/from which the de-refined element was constructed.
        /// First and after suffix, indicates the order in which the elements are to be processed.
        /// Enumeration values and comments are from the original Fortran code.
        enum class ElementType
        {
            WasNodeFirst = 1,  //< front, 'A' face (used to be node, delete it):  1
            WasEdgeFirst = 2,  //< front, 'B' face (used to be edge, keep it):    2
            WasFace = 3,       //< 'C' face (used to be face, keep it):           3
            WasNodeAfter = -1, //< not in front, 'A' face:                       -1
            WasEdgeAfter = -2, //< not in front, 'B' face:                       -2
            Unassigned = 0     //< not assigned a value                           0
        };

        /// @brief Indicates what should happen
        enum class ResultIndicator
        {
            ReturnFromFunction, //< Return (false) from the current function
            BreakInnerLoop,     //< Break the inner loop
            ContinueInnerLoop   //< Continue in the inner loop
        };

        /// @brief Indicate if the element can be a seed element or not.
        static bool ElementIsSeed(const Mesh2D& mesh,
                                  const std::vector<int>& nodeTypes,
                                  const UInt element);

        /// @brief Find the seed element id to start the mesh de-refinement.
        static UInt FindElementSeedIndex(const Mesh2D& mesh,
                                         const std::vector<int>& nodeTypes);

        /// @brief Find all elements that are connected along edges to elementId.
        static void FindDirectlyConnectedFaces(const Mesh2D& mesh,
                                               const UInt elementId,
                                               std::vector<UInt>& connected);

        /// @brief Find all elements that are connected by nodes to elementId.
        static void FindIndirectlyConnectedFaces(const Mesh2D& mesh,
                                                 const UInt elementId,
                                                 const std::vector<UInt>& directlyConnected,
                                                 std::vector<UInt>& indirectlyConnected);

        /// @brief Assign directly connected element indices
        static void AssignDirectlyConnected(const std::vector<UInt>& directlyConnected,
                                            std::array<int, 2>& edgeFaces,
                                            UInt& neighbouringElementId);

        /// @brief Assign indirectly connected element indices
        static void AssignIndirectlyConnected(const std::vector<UInt>& indirectlyConnected,
                                              std::array<int, 2>& edgeFaces,
                                              const UInt neighbouringElementId);

        /// @brief Find element id's
        static void FindAdjacentFaces(const Mesh2D& mesh,
                                      const std::vector<UInt>& directlyConnected,
                                      const std::vector<UInt>& indirectlyConnected,
                                      std::vector<std::array<int, 2>>& edgeFaces);

        /// @brief Find the elements that are connected to the elementId.
        static void FindSurroundingFaces(const Mesh2D& mesh,
                                         const UInt elementId,
                                         std::vector<UInt>& directlyConnected,
                                         std::vector<UInt>& indirectlyConnected,
                                         std::vector<std::array<int, 2>>& edgeFaces);

        /// @brief Update face mask for directly connected elements
        static void UpdateFaceMaskDirectlyConnectedNodeFirst(const std::vector<UInt>& directlyConnected,
                                                             const Mesh2D& mesh,
                                                             const std::vector<UInt>& frontIndex,
                                                             std::vector<UInt>& frontIndexCopy,
                                                             std::vector<ElementType>& faceMask);

        /// @brief Update face mask for indirectly connected elements
        static void UpdateFaceMaskIndirectlyConnectedNodeFirst(const std::vector<UInt>& directlyConnected,
                                                               const Mesh2D& mesh,
                                                               std::vector<ElementType>& faceMask);

        /// @brief Update face mask for directly connected elements
        static void UpdateFaceMaskDirectlyConnectedEdgeFirst(const std::vector<UInt>& directlyConnected,
                                                             const Mesh2D& mesh,
                                                             const std::vector<UInt>& frontIndex,
                                                             std::vector<UInt>& frontIndexCopy,
                                                             std::vector<ElementType>& faceMask);

        /// @brief Update face mask for indirectly connected elements
        static void UpdateFaceMaskIndirectlyConnectedEdgeFirst(const std::vector<UInt>& indirectlyConnected,
                                                               const Mesh2D& mesh,
                                                               const std::vector<UInt>& frontIndex,
                                                               std::vector<UInt>& frontIndexCopy,
                                                               std::vector<ElementType>& faceMask);

        /// @brief Initialise the element mask.
        static std::vector<ElementType> InitialiseElementType(const Mesh2D& mesh,
                                                              const std::vector<int>& nodeTypes);

        /// \brief Determine if the element can be deleted from the mesh or not.
        static bool ElementCannotBeDeleted(const Mesh2D& mesh,
                                           const std::vector<int>& nodeTypes,
                                           const Polygons& polygon,
                                           const UInt elementId);

        /// @brief Compute coordinates of the new node.
        static Point ComputeNewNodeCoordinates(const Mesh2D& mesh,
                                               const std::vector<int>& nodeTypes,
                                               const UInt nodeId);

        /// @brief Get the element index
        static UInt GetElementIndex(const std::array<int, 2>& edgeFaces,
                                    const UInt index);

        /// @brief Find a common edge between elements
        static std::tuple<UInt, UInt> FindCommonEdge(Mesh2D& mesh,
                                                     const UInt leftElementId,
                                                     const UInt rightElementId,
                                                     const UInt connectedElementId);

        /// @brief Update the mesh members for the mesh description and connectivity.
        [[nodiscard]] static bool UpdateDirectlyConnectedElements(Mesh2D& mesh,
                                                                  const UInt elementId,
                                                                  const std::vector<UInt>& directlyConnected,
                                                                  const std::vector<std::array<int, 2>>& edgeFaces);

        /// @brief Update the mesh members for the mesh description and connectivity for triangle elements
        static bool UpdateDirectlyConnectedTriangleElements(Mesh2D& mesh,
                                                            const UInt index,
                                                            const UInt connectedElementId,
                                                            const std::vector<std::array<int, 2>>& edgeFaces);

        /// @brief Update the mesh members for the mesh description and connectivity for non-triangle elements
        ///
        /// That is, element with 4 or more edges.
        static void UpdateDirectlyConnectedNonTriangleElements(Mesh2D& mesh,
                                                               const UInt index,
                                                               const UInt elementId,
                                                               const UInt connectedElementId);

        /// @brief Get the most significant node type for all nodes of the element.
        static int GetNodeCode(const Mesh2D& mesh,
                               const std::vector<int>& nodeTypes,
                               const UInt elementId);

        /// @brief Add element id to the list of id's
        ///
        /// only added is it is not already on the list and the element is a quadrilateral
        static void AddElementToList(const Mesh& mesh,
                                     const std::vector<UInt>& frontList,
                                     std::vector<UInt>& frontListCopy,
                                     const UInt elementId);

        /// @brief Redirect nodes of connected faces, deactivate polygons of degree smaller than three
        static void RedirectNodesOfConnectedElements(Mesh2D& mesh,
                                                     const UInt elementId,
                                                     const UInt nodeId,
                                                     const std::vector<UInt>& indirectlyConnected);

        /// @brief Remove a boundary node or edge
        static ResultIndicator RemoveBoundaryNodeAndEdge(Mesh2D& mesh,
                                                         const Polygons& polygon,
                                                         const std::vector<int>& nodeTypes,
                                                         const UInt faceNodeIndex,
                                                         const UInt connectedElementId,
                                                         const UInt edgeId,
                                                         const UInt previousEdgeId,
                                                         const UInt nodeId);

        /// @brief Removes nodes from the boundary that will not be part of the de-refined mesh.
        [[nodiscard]] static bool RemoveUnwantedBoundaryNodes(Mesh2D& mesh,
                                                              const std::vector<int>& nodeTypes,
                                                              const Polygons& polygon,
                                                              const std::vector<UInt>& indirectlyConnected);

        /// @brief Delete an element
        [[nodiscard]] static bool DeleteElement(Mesh2D& mesh,
                                                std::vector<int>& nodeTypes,
                                                const Polygons& polygon,
                                                const UInt elementId,
                                                const std::vector<UInt>& directlyConnected,
                                                const std::vector<UInt>& indirectlyConnected,
                                                const std::vector<std::array<int, 2>>& edgeFaces);

        /// @brief Clean up the edge
        ///
        /// @returns Indicates if the cleanp-up was successful or not
        [[nodiscard]] static bool CleanUpEdge(Mesh2D& mesh, const UInt edgeId);

        /// @brief Do the Casullu de-refinement
        [[nodiscard]] static bool DoDeRefinement(Mesh2D& mesh, const Polygons& polygon);

        /// @brief Compute the mesh node types.
        ///
        /// Uses the m_nodeTypes has been generated in the mesh.
        static std::vector<int> ComputeNodeTypes(const Mesh2D& mesh, const Polygons& polygon);
    };

} // namespace meshkernel
