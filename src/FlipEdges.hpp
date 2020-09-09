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

#include <vector>
#include "Entities.hpp"


namespace MeshKernel
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
        /// <returns></returns>
        FlipEdges();

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

        /// <summary>
        /// Flip the edges
        /// </summary>
        /// <returns>If the method succeeded</returns>
        bool Compute();

    private:

        /// <summary>
        /// Transform non-triangular faces in triangular faces
        /// </summary>
        /// <returns>If the method succeeded</returns>
        bool TriangulateFaces();

        /// <summary>
        /// Computes the change in topology functional and gets the nodes involved (comp_ntopo)
        /// </summary>
        /// <param name="edge">The current edge</param>
        /// <param name="nodeLeft">The node at the left side of the edge</param>
        /// <param name="nodeRight">The node at the left side of the edge</param>
        /// <param name="topologyFunctional">The computed functional</param>
        /// <returns>If the method succeeded</returns>
        bool ComputeTopologyFunctional( int edge,
                                 int& nodeLeft,
                                 int& nodeRight,
                                 int& topologyFunctional) const;

        /// <summary>
        /// Determine the optimal number of connected nodes for each node (nmk_opt)
        /// </summary>
        /// <param name="nodeIndex"></param>
        /// <returns>If the method succeeded</returns>
        int OptimalNumberOfConnectedNodes(int nodeIndex) const;

        /// <summary>
        /// Compute the difference with the optimal number of edges by counting the numbers of edges that 
        /// connect nodes firstNode and secondNode, and are on the land boundary path (comp_nnow)
        /// </summary>
        /// <returns>If the method succeeded</returns>
        int DifferenceFromOptimum(int nodeIndex, int firstNode, int secondNode) const;

        /// <summary>
        /// Remove a connected edge from a node
        /// </summary>
        /// <param name="edgeIndex">The index of the edge to remove</param>
        /// <param name="nodeIndex">The index of the node to process</param>
        /// <returns>If the method succeeded</returns>
        bool DeleteEdgeFromNode(int edgeIndex, int nodeIndex) const;

        std::shared_ptr<Mesh> m_mesh;                                      // A pointer to mesh
        std::shared_ptr<LandBoundaries> m_landBoundaries;                  // A pointer to land boundaries

        bool m_triangulateFaces = false;
        bool m_projectToLandBoundary = false;

    };

}


