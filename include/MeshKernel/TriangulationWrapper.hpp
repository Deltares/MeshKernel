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

namespace meshkernel
{
    extern "C"
    {
        void Triangulation(int* jatri, double* xs, double* ys, int* ns, int* indx, int* numtri, int* edgeidx, int* numedge, int* triedge, double* xs3, double* ys3, int* ns3, double* trisize);
    }

    struct Point;
    struct Sample;
    struct TriangulationWrapper
    {
        enum class TriangulationOptions
        {
            TriangulatePoints = 1,                // generate Delaunay triangulation from input nodes
            GeneratePoints = 2,                   // generate internal nodes in polygon that produce a Delaunay triangulation
            TriangulatePointsAndGenerateFaces = 3 // generate Delaunay triangulation from input nodes with m_faceEdges and m_edgeNodes
        };

        std::vector<Point> m_nodes;
        std::vector<std::vector<int>> m_faceNodes;
        std::vector<std::vector<int>> m_faceEdges;
        std::vector<std::vector<int>> m_edgeNodes;
        std::vector<std::vector<int>> m_edgesFaces;

        int m_numEdges;
        int m_numNodes;
        int m_numFaces;

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
                     int estimatedNumberOfTriangles)
        {
            // TODO: move implementation to source files as soon as we can use modules (c++20)

            std::vector<double> xLocalPolygon(inputNodes.size());
            std::vector<double> yLocalPolygon(inputNodes.size());
            for (int i = 0; i < inputNodes.size(); ++i)
            {
                xLocalPolygon[i] = inputNodes[i].x;
                yLocalPolygon[i] = inputNodes[i].y;
            }

            m_numFaces = -1;
            m_numEdges = 0;
            m_numNodes = 0;
            auto intTriangulationOption = static_cast<int>(triangulationOption);

            std::vector<int> faceNodesFlat;
            std::vector<int> edgeNodesFlat;
            std::vector<int> faceEdgesFlat;
            std::vector<double> xNodesFlat;
            std::vector<double> yNodesFlat;

            if (estimatedNumberOfTriangles <= 0)
            {
                estimatedNumberOfTriangles = static_cast<int>(inputNodes.size()) * 6 + 10;
            }

            // If the number of estimated triangles is not sufficient, triangulation must be repeated
            while (m_numFaces < 0)
            {
                m_numFaces = estimatedNumberOfTriangles;
                faceNodesFlat.resize(int(estimatedNumberOfTriangles) * 3);
                edgeNodesFlat.resize(int(estimatedNumberOfTriangles) * 2);
                faceEdgesFlat.resize(int(estimatedNumberOfTriangles) * 3);
                xNodesFlat.resize(int(estimatedNumberOfTriangles) * 3, doubleMissingValue);
                yNodesFlat.resize(int(estimatedNumberOfTriangles) * 3, doubleMissingValue);
                Triangulation(&intTriangulationOption,
                              &xLocalPolygon[0],
                              &yLocalPolygon[0],
                              &numPoints,
                              &faceNodesFlat[0], // INDX
                              &m_numFaces,
                              &edgeNodesFlat[0], // EDGEINDX
                              &m_numEdges,
                              &faceEdgesFlat[0], // TRIEDGE
                              &xNodesFlat[0],
                              &yNodesFlat[0],
                              &m_numNodes,
                              &averageTriangleArea);
                if (estimatedNumberOfTriangles)
                {
                    estimatedNumberOfTriangles = -m_numFaces;
                }
            }

            // Create nodes
            m_nodes.resize(m_numNodes);
            for (int i = 0; i < m_numNodes; ++i)
            {
                m_nodes[i] = {xNodesFlat[i], yNodesFlat[i]};
            }

            // Create m_faceNodes
            m_faceNodes.resize(m_numFaces, std::vector<int>(3, intMissingValue));
            m_faceEdges.resize(m_numFaces, std::vector<int>(3, intMissingValue));
            int faceCounter = 0;
            for (int f = 0; f < m_numFaces; ++f)
            {
                m_faceNodes[f][0] = faceNodesFlat[faceCounter] - 1;
                m_faceEdges[f][0] = faceEdgesFlat[faceCounter] - 1;
                faceCounter++;
                m_faceNodes[f][1] = faceNodesFlat[faceCounter] - 1;
                m_faceEdges[f][1] = faceEdgesFlat[faceCounter] - 1;
                faceCounter++;
                m_faceNodes[f][2] = faceNodesFlat[faceCounter] - 1;
                m_faceEdges[f][2] = faceEdgesFlat[faceCounter] - 1;
                faceCounter++;
            }

            // Create edges
            if (m_numEdges <= 0)
            {
                return;
            }

            m_edgeNodes.resize(m_numEdges, std::vector<int>(2, intMissingValue));
            int edgeCounter = 0;
            for (int e = 0; e < m_numEdges; ++e)
            {
                m_edgeNodes[e][0] = edgeNodesFlat[edgeCounter] - 1;
                edgeCounter++;
                m_edgeNodes[e][1] = edgeNodesFlat[edgeCounter] - 1;
                edgeCounter++;
            }

            m_edgesFaces.resize(m_numEdges, std::vector<int>(2, intMissingValue));
            edgeCounter = 0;
            for (int f = 0; f < m_numFaces; ++f)
            {

                for (int n = 0; n < 3; ++n)
                {
                    auto const edge = faceEdgesFlat[edgeCounter] - 1;
                    edgeCounter++;
                    // For each edge, the shared face index
                    if (m_edgesFaces[edge][0] == intMissingValue)
                    {
                        m_edgesFaces[edge][0] = f;
                    }
                    else
                    {
                        m_edgesFaces[edge][1] = f;
                    }
                }
            }
        }
    };

} // namespace meshkernel
