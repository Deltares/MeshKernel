//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2023.
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
#include <utility>
#include <vector>

#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/Point.hpp"

namespace meshkernel
{

    /// @brief Connects grids across element faces.
    ///
    /// Can be connected with upto 5 elements or 4 hanging nodes along an edge.
    /// @note Currently only works for quadrilateral elements.
    // Functionality re-implemented from connectcurvilinearquadsddtype.f90
    class ConnectMeshes final
    {
    public:
        /// @brief The default value of the fraction of the edge length used to determine if edges are adjacent.
        static constexpr double DefaultSeparationFraction = 0.4;

        /// @brief Connect grids.
        ///
        /// @param [in,out] mesh The mesh
        /// @param [in] separationFraction The fraction of the shortest edge to use when determining neighbour edge closeness
        void Compute(Mesh2D& mesh, const double separationFraction = DefaultSeparationFraction) const;

    private:
        /// @brief The maximum number of hanging nodes along a single element edge
        static constexpr UInt m_maximumNumberOfIrregularElementsAlongEdge = 5;

        /// @brief Indicates if node should be merged or not
        enum class MergeIndicator
        {
            Initial,
            DoMerge,
            DoNotMerge
        };

        /// @brief Bounded array of integer values.
        using BoundedIntegerArray = std::array<UInt, m_maximumNumberOfIrregularElementsAlongEdge>;

        /// @brief A pair of node-id's that can be merged into a single node.
        using NodesToMerge = std::pair<UInt, UInt>;

        /// @brief Contains irregular edge information.
        struct IrregularEdgeInfo
        {
            ///@brief Hanging node indices
            BoundedIntegerArray hangingNodes{constants::missing::uintValue,
                                             constants::missing::uintValue,
                                             constants::missing::uintValue,
                                             constants::missing::uintValue,
                                             constants::missing::uintValue};
            /// @brief Number of hanging nodes along single irregular edge.
            UInt edgeCount = 0;
            /// @brief Start node of the irregular edge
            UInt startNode = constants::missing::uintValue;
            /// @brief End node of the irregular edge
            UInt endNode = constants::missing::uintValue;
        };

        /// @brief Array of IrregularEdgeInfo structs.
        using IrregularEdgeInfoArray = std::vector<IrregularEdgeInfo>;

        /// @brief Determine if the edges are adjacent, if so then return the start and end points (adjacent.f90)
        ///
        /// @param [in] mesh The mesh
        /// @param [in] separationFraction The fraction of the shortest edge to use when determining neighbour edge closeness
        /// @param [in] edge1 One of the edges for adjacency check
        /// @param [in] edge2 Another of the edges for adjacency check
        /// @param [out] areAdjacent Indicates is edge1 and edge2 are adjacent
        /// @param [out] startNode End point nodes, if not nullvalue then node is hanging node
        /// @param [out] endNode End point nodes, if not nullvalue then node is hanging node
        void AreEdgesAdjacent(const Mesh2D& mesh,
                              const double separationFraction,
                              const UInt edge1,
                              const UInt edge2,
                              bool& areAdjacent,
                              UInt& startNode,
                              UInt& endNode) const;

        /// @brief Find all quadrilateral elements that do no have a neighbour across any of edges.
        ///
        /// @param [in] mesh The mesh
        /// @param [in,out] elementsOnDomainBoundary List of elements that do not have neighbour
        /// @param [in,out] edgesOnDomainBoundary List of edges that do have elements on one side only
        void GetQuadrilateralElementsOnDomainBoundary(const Mesh2D& mesh,
                                                      std::vector<UInt>& elementsOnDomainBoundary,
                                                      std::vector<UInt>& edgesOnDomainBoundary) const;

        /// @brief Get list of node id's ordered with distance from given point.
        ///
        /// @param [in] mesh The mesh
        /// @param [in] nodeIndices List of nodes to have distance computed
        /// @param [in] numberOfNodes Number of nodes to have distance computed
        /// @param [in] point Start point from which node distances are to be computed
        /// @param [out] nearestNeighbours List of nearest nodes in order of distance from point
        void GetOrderedDistanceFromPoint(const Mesh2D& mesh,
                                         const std::vector<UInt>& nodeIndices,
                                         const UInt numberOfNodes,
                                         const Point& point,
                                         BoundedIntegerArray& nearestNeighbours) const;

        /// @brief Merge coincident nodes
        ///
        /// @param [in,out] mesh The mesh
        /// @param [in] nodesToMerge List of nodes to be merged
        /// @param [in,out] mergeIndicator Indicates if node needs to be merged.
        void MergeNodes(Mesh2D& mesh, const std::vector<NodesToMerge>& nodesToMerge, std::vector<MergeIndicator>& mergeIndicator) const;

        /// @brief Free one hanging node along an irregular edge.
        ///
        /// @brief [in] mesh The mesh
        /// @brief [in] hangingNodes List of hanging nodes for edge
        /// @brief [in] startNode End point of regular edge, to which the hanging node will be connected
        /// @brief [in] endNode Other end point of regular edge, to which the hanging node will be connected
        void FreeOneHangingNode(Mesh2D& mesh,
                                const BoundedIntegerArray& hangingNodes,
                                const UInt startNode,
                                const UInt endNode) const;

        /// @brief Free two hanging nodes along an irregular edge.
        ///
        /// @brief [in] mesh The mesh
        /// @brief [in] faceId The element with irregular edge, required to get next adjacent element
        /// @brief [in] edgeId Edge along opposite side of irregular edge, required to get next adjacent element
        /// @brief [in] hangingNodes List of hanging nodes for edge
        /// @brief [in] startNode End point of regular edge, to which the hanging nodes will be connected
        /// @brief [in] endNode Other end point of regular edge, to which the hanging nodes will be connected
        void FreeTwoHangingNodes(Mesh2D& mesh,
                                 const UInt faceId,
                                 const UInt edgeId,
                                 const BoundedIntegerArray& hangingNodes,
                                 const UInt startNode,
                                 const UInt endNode) const;

        /// @brief Free three hanging nodes along an irregular edge.
        ///
        /// @brief [in] mesh The mesh
        /// @brief [in] faceId The element with irregular edge, required to get next adjacent element
        /// @brief [in] edgeId Edge along opposite side of irregular edge, required to get next adjacent element
        /// @brief [in] hangingNodes List of hanging nodes for edge
        /// @brief [in] startNode End point of regular edge, to which the hanging nodes will be connected
        /// @brief [in] endNode Other end point of regular edge, to which the hanging nodes will be connected
        void FreeThreeHangingNodes(Mesh2D& mesh,
                                   const UInt faceId,
                                   const UInt edgeId,
                                   const BoundedIntegerArray& hangingNodes,
                                   const UInt startNode,
                                   const UInt endNode) const;

        /// @brief Free four hanging nodes along an irregular edge.
        ///
        /// @brief [in] mesh The mesh
        /// @brief [in] faceId The element with irregular edge, required to get next adjacent element
        /// @brief [in] edgeId Edge along opposite side of irregular edge, required to get next adjacent element
        /// @brief [in] hangingNodes List of hanging nodes for edge
        /// @brief [in] startNode End point of regular edge, to which the hanging nodes will be connected
        /// @brief [in] endNode Other end point of regular edge, to which the hanging nodes will be connected
        void FreeFourHangingNodes(Mesh2D& mesh,
                                  const UInt faceId,
                                  const UInt edgeId,
                                  const BoundedIntegerArray& hangingNodes,
                                  const UInt startNode,
                                  const UInt endNode) const;

        /// @brief Free any hanging nodes along an irregular edge.
        ///
        /// @brief [in,out] mesh The mesh
        /// @brief [in] numberOfHangingNodes Number of hanging nodes along irregular edge
        /// @brief [in] hangingNodesOnEdge Id's of nodes aong irregular edge
        /// @brief [in] faceId The element with irregular edge, required to get next adjacent element
        /// @brief [in] boundaryEdge The irregular edge
        /// @brief [in] boundaryNode End point of erregular edge, required to order the hanging nodes
        /// @brief [in] edgeId Edge along opposite side of irregular edge, required to get next adjacent element
        void FreeHangingNodes(Mesh2D& mesh,
                              const UInt numberOfHangingNodes,
                              const std::vector<UInt>& hangingNodesOnEdge,
                              const UInt faceId,
                              const Edge& boundaryEdge,
                              const Point& boundaryNode,
                              const UInt edgeId) const;

        /// @brief Find and retain any hanging node id's
        ///
        /// @param [in] mesh The mesh
        /// @param [in] separationFraction The fraction of the shortest edge to use when determining neighbour edge closeness
        /// @param [in] edgesOnDomainBoundary List of edges along domain boundary, more specifically edges with only a single element attached
        /// @param [in,out] irregularEdges List of irregular edges with hanging nodes
        void GatherHangingNodeIds(const Mesh2D& mesh,
                                  const double separationFraction,
                                  const std::vector<UInt>& edgesOnDomainBoundary,
                                  IrregularEdgeInfoArray& irregularEdges) const;

        /// @brief Gather all the nodes that need to be merged.
        ///
        /// @param [in] startNode End point of edge
        /// @param [in] endNode Other end point of edge
        /// @param [in] boundaryEdge Edge with possible coinciding nodes
        /// @param [in,out] nodesToMerge List of nodes to be merged
        /// @param [in,out] mergeIndicator List of indicators, indicating if node has been processed and should be merged or not
        ///
        /// Before merging there may be 2 or more nodes that are at the same point.
        void GatherNodesToMerge(const UInt startNode,
                                const UInt endNode,
                                const Edge& boundaryEdge,
                                std::vector<NodesToMerge>& nodesToMerge,
                                std::vector<MergeIndicator>& mergeIndicator) const;

        /// @brief Gather hanging nodes along the irregular edge.
        ///
        /// @param [in] primaryStartNode End point of irregular edge
        /// @param [in] primaryEndNode Other end point of irregular edge
        /// @param [in] irregularEdge Edge with hanging nodes
        /// @param [in,out] hangingNodesOnEdge List of hanging node along edge
        /// @param [in,out] numberOfHangingNodes Number of hanging nodes along edge
        /// @param [in,out] mergeIndicator List of indicators, indicating if node has been processed and should be merged or not
        void GatherHangingNodes(const UInt primaryStartNode,
                                const UInt primaryEndNode,
                                const Edge& irregularEdge,
                                std::vector<UInt>& hangingNodesOnEdge,
                                UInt& numberOfHangingNodes,
                                std::vector<MergeIndicator>& mergeIndicator) const;
    };

} // namespace meshkernel
