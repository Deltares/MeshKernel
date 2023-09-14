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

    /// @brief Connects curvilinear grids across element faces.
    ///
    /// Can be connected with upto 5 elements or 4 hanging nodes along an edge.
    /// @note Currently only works for quadrilateral elements.
    class ConnectCurvilinearGrids final
    {
    public:
        /// @brief Connect curvilinear grids (connectcurvilinearquadsddtype.f90)
        void Compute(Mesh2D& mesh) const;

    private:
        /// @brief The maximum number of hanging nodes along a single element edge
        static constexpr UInt m_maximumNumberOfIrregularElementsAlongEdge = 5;

        /// @brief Bounded array or integer values.
        using BoundedIntegerArray = std::array<UInt, m_maximumNumberOfIrregularElementsAlongEdge>;

        /// @brief Determine if the edges are adjacent, if so then return the start and end points (adjacent.f90)
        void AreEdgesAdjacent(const Mesh2D& mesh, const UInt edge1, const UInt edge2, bool& areAdjacent, UInt& startNode, UInt& endNode) const;

        /// @brief Find all quadrilateral elements that do no have a neighbour across any of edges.
        void GetQuadrilateralElementsOnDomainBoundary(const Mesh2D& mesh, std::vector<UInt>& elementsOnDomainBoundary, std::vector<UInt>& edgesOnDomainBoundary) const;

        /// @brief Get list of node id's ordered with distance from given point.
        void GetOrderedDistanceFromPoint(const Mesh2D& mesh,
                                         const std::vector<UInt>& nodeIndices,
                                         const UInt numberOfNodes,
                                         const Point& point,
                                         BoundedIntegerArray& nearestNeighbours) const;

        /// @brief Merge coincident nodes
        void MergeNodes(Mesh2D& mesh, const std::vector<std::pair<UInt, UInt>>& nodesToMerge, std::vector<int>& mergeIndicator) const;

        /// @brief Free one hanging node along an irregular edge.
        void FreeOneHangingNode(Mesh2D& mesh, const BoundedIntegerArray& hangingNodes, const UInt startNode, const UInt endNode) const;

        /// @brief Free two hanging nodes along an irregular edge.
        void FreeTwoHangingNodes(Mesh2D& mesh, const UInt faceId, const UInt edgeId, const BoundedIntegerArray& hangingNodes, const UInt startNode, const UInt endNode) const;

        /// @brief Free three hanging nodes along an irregular edge.
        void FreeThreeHangingNodes(Mesh2D& mesh, const UInt faceId, const UInt edgeId, const BoundedIntegerArray& hangingNodes, const UInt startNode, const UInt endNode) const;

        /// @brief Free four hanging nodes along an irregular edge.
        void FreeFourHangingNodes(Mesh2D& mesh, const UInt faceId, const UInt edgeId, const BoundedIntegerArray& hangingNodes, const UInt startNode, const UInt endNode) const;

        /// @brief Free any hanging nodes along an irregular edge.
        void FreeHangingNodes(Mesh2D& mesh,
                              const UInt numberOfHangingNodes,
                              const std::vector<UInt>& hangingNodesOnEdge,
                              const UInt faceId,
                              const Edge& boundaryEdge,
                              const Point& boundaryNode,
                              const UInt edgeId) const;

        /// @brief Find and retain any hanging node id's
        void GatherHangingNodeIds(const Mesh2D& mesh,
                                  const std::vector<UInt>& edgesOnDomainBoundary,
                                  std::vector<BoundedIntegerArray>& irregularEdge,
                                  std::vector<UInt>& irregularEdgeCount,
                                  std::vector<UInt>& irregularStartNodes,
                                  std::vector<UInt>& irregularEndNodes) const;

        /// @brief Gather all the nodes that need to be merged.
        ///
        /// Before merging there may be 2 or more nodes that are at the same point.
        void GatherNodesToMerge(const UInt startNode,
                                const UInt endNode,
                                const Edge& boundaryEdge,
                                std::vector<std::pair<UInt, UInt>>& nodesToMerge,
                                std::vector<int>& mergeIndicator) const;

        /// @brief Gather hanging nodes along the irregular edge.
        void GatherHangingNodes(const UInt primaryStartNode,
                                const UInt primaryEndNode,
                                const Edge& irregularEdge,
                                std::vector<UInt>& hangingNodesOnEdge,
                                UInt& numberOfHangingNodes,
                                std::vector<int>& mergeIndicator) const;
    };

} // namespace meshkernel
