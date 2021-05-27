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
#include <memory>

namespace meshkernel
{
    // Forward declarations
    class Mesh2D;
    class LandBoundaries;

    /// @brief A class used to improve mesh connectivity.
    ///
    /// The edges are flipped in order to reduce the number of edges connected to a node
    /// The optimal number of edges to a node is six.
    /// If we additionally set `triangulateFaces`to true, that results in triangles of 60° in each angle and therefore a nearly ideal mesh.
    /// An additional option `projectToLandBoundary` defines whether we want project to the land boundary.
    class FlipEdges
    {
    public:
        /// @brief Constructor
        /// @param[in] mesh                  The input mesh
        /// @param[in] landBoundary          The land boundary
        /// @param[in] triangulateFaces      Whether to triangulate all faces or not
        /// @param[in] projectToLandBoundary Whether to project to land boundaries or not
        /// @returns If the method succeeded
        FlipEdges(std::shared_ptr<Mesh2D> mesh,
                  std::shared_ptr<LandBoundaries> landBoundary,
                  bool triangulateFaces,
                  bool projectToLandBoundary);

        /// @brief Flip the edges
        void Compute() const;

    private:
        /// @brief Computes the change in topology functional and gets the nodes involved (comp_ntopo)
        /// @param[in]  edge      The current edge
        /// @param[out] nodeLeft  The node at the left side of the edge
        /// @param[out] nodeRight The node at the left side of the edge
        /// @return topologyFunctional The computed functional
        int ComputeTopologyFunctional(size_t edge,
                                      size_t& nodeLeft,
                                      size_t& nodeRight) const;

        /// @brief Determine the optimal number of connected nodes for each node (nmk_opt)
        /// @param nodeIndex
        /// @returns Optimal number of connected nodes
        [[nodiscard]] size_t OptimalNumberOfConnectedNodes(size_t nodeIndex) const;

        /// @brief Compute the difference with the optimal number of edges by counting the numbers of edges that
        ///        connect nodes firstNode and secondNode, and are on the land boundary path (comp_nnow)
        /// @returns Difference form optimum
        [[nodiscard]] int DifferenceFromOptimum(size_t nodeIndex, size_t firstNode, size_t secondNode) const;

        /// @brief Remove a connected edge from a node
        /// @param[in] edgeIndex The index of the edge to remove
        /// @param[in] nodeIndex The index of the node to process
        void DeleteEdgeFromNode(size_t edgeIndex, size_t nodeIndex) const;

        std::shared_ptr<Mesh2D> m_mesh;                   ///< A pointer to the 2D mesh
        std::shared_ptr<LandBoundaries> m_landBoundaries; ///< A pointer to the land boundaries

        bool m_triangulateFaces = false;      ///< Whether to triangulate faces
        bool m_projectToLandBoundary = false; ///< Whether to project to land boundary
    };

} // namespace meshkernel
