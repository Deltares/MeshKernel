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

namespace meshkernel
{
    struct TriangulationWrapper
    {
        enum class TriangulationOptions
        {
            TriangulatePoints = 1,                // generate Delaunay triangulation from input nodes
            GeneratePoints = 2,                   // generate internal nodes in polygon that produce a Delaunay triangulation
            TriangulatePointsAndGenerateFaces = 3 // generate Delaunay triangulation from input nodes with m_faceEdges and m_edgeNodes
        };

        /// @brief
        /// @tparam T A type that contains x and y fields
        /// @param inputNodes The input points
        /// @param numPoints The number of input points
        /// @param triangulationOption Triangulation option, see \ref TriangulationOptions
        /// @param averageTriangleArea An estimation of the average area of triangles (required for option 2)
        /// @param estimatedNumberOfTriangles An estimation of the average number of triangles (required for option 2)
        template <typename T>
        void Compute(const std::vector<T>& inputNodes,
                     int numPoints,
                     TriangulationOptions triangulationOption,
                     double averageTriangleArea,
                     int estimatedNumberOfTriangles); //requires IsCoordinate<T>;

        std::vector<Point> m_nodes;
        std::vector<std::vector<int>> m_faceNodes;
        std::vector<std::vector<int>> m_faceEdges;
        std::vector<std::vector<int>> m_edgeNodes;
        std::vector<std::vector<int>> m_edgesFaces;

        int m_numEdges;
        int m_numNodes;
        int m_numFaces;
    };

} // namespace meshkernel