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

namespace meshkernel
{
    // Forward declarations
    class Mesh;
    class LandBoundaries;

    class FlipEdges
    {
    public:
        /// <summary>
        /// Ctor
        /// </summary>
        /// <param name="mesh">The input mesh</param>
        /// <param name="landBoundary">The land boundary</param>
        /// <param name="triangulateFaces">Option to triangulate all faces or not</param>
        /// <param name="projectToLandBoundary">Option to project to land boundaries or not</param>
        /// <returns>If the method succeeded</returns>
        FlipEdges(std::shared_ptr<Mesh> mesh,
                  std::shared_ptr<LandBoundaries> landBoundary,
                  bool triangulateFaces,
                  bool projectToLandBoundary);

        /// @brief Flip the edges
        void Compute() const;

    private:
        /// @brief Computes the change in topology functional and gets the nodes involved (comp_ntopo)
        /// @param[in] edge The current edge
        /// @param[out] nodeLeft The node at the left side of the edge
        /// @param[out] nodeRight The node at the left side of the edge
        /// @param[out] topologyFunctional The computed functional
        void ComputeTopologyFunctional(int edge,
                                       int& nodeLeft,
                                       int& nodeRight,
                                       int& topologyFunctional) const;

        /// @brief Determine the optimal number of connected nodes for each node (nmk_opt)
        /// @param nodeIndex
        /// @returns Optimal number of connected nodes
        [[nodiscard]] int OptimalNumberOfConnectedNodes(int nodeIndex) const;

        /// @brief Compute the difference with the optimal number of edges by counting the numbers of edges that
        ///        connect nodes firstNode and secondNode, and are on the land boundary path (comp_nnow)
        /// @returns Difference form optimum
        [[nodiscard]] int DifferenceFromOptimum(int nodeIndex, int firstNode, int secondNode) const;

        /// @brief Remove a connected edge from a node
        /// @param[in] edgeIndex The index of the edge to remove
        /// @param[in] nodeIndex The index of the node to process
        void DeleteEdgeFromNode(int edgeIndex, int nodeIndex) const;

        std::shared_ptr<Mesh> m_mesh;                     // A pointer to mesh
        std::shared_ptr<LandBoundaries> m_landBoundaries; // A pointer to land boundaries

        bool m_triangulateFaces = false;
        bool m_projectToLandBoundary = false;
    };

} // namespace meshkernel
