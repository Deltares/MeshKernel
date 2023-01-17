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

namespace meshkernel
{
    extern "C"
    {
        /// @brief Function of the Triangle library
        ///
        /// \see https://www.cs.cmu.edu/~quake/triangle.html
        void Triangulation(int* jatri, double* xs, double* ys, int* ns, int* indx, int* numtri, int* edgeidx, int* numedge, int* triedge, double* xs3, double* ys3, int* ns3, double* trisize);
    }

    struct Point;
    struct Sample;
    /// @brief Wrapper around the Triangle library
    ///
    /// \see https://www.cs.cmu.edu/~quake/triangle.html
    struct TriangulationWrapper
    {
        /// @brief Enumerator describing all triangulation options
        enum class TriangulationOptions
        {
            TriangulatePoints = 1,                ///< generate Delaunay triangulation from input nodes
            GeneratePoints = 2,                   ///< generate internal nodes in polygon that produce a Delaunay triangulation
            TriangulatePointsAndGenerateFaces = 3 ///< generate Delaunay triangulation from input nodes with m_faceEdges and m_edgeNodes
        };

        std::vector<Point> m_nodes;                    ///< Nodes
        std::vector<std::vector<size_t>> m_faceNodes;  ///< Face nodes
        std::vector<std::vector<size_t>> m_faceEdges;  ///< Face edges
        std::vector<std::vector<size_t>> m_edgeNodes;  ///< Edge nodes
        std::vector<std::vector<size_t>> m_edgesFaces; ///< Edge faces

        size_t m_numEdges; ///< Number of edges
        size_t m_numNodes; ///< Number of nodes
        size_t m_numFaces; ///< Number of faces

        /// @brief
        /// @tparam T A type that contains x and y fields
        /// @param inputNodes The number of input points
        /// @param triangulationOption Triangulation option, see \ref TriangulationOptions
        /// @param averageTriangleArea An estimation of the average area of triangles (required for option 2)
        /// @param estimatedNumberOfTriangles An estimation of the average number of triangles (required for option 2)
        template <typename T>
        void Compute(const std::vector<T>& inputNodes,
                     TriangulationOptions triangulationOption,
                     double averageTriangleArea,
                     size_t estimatedNumberOfTriangles)
        {

            std::vector<double> xLocalPolygon(inputNodes.size());
            std::vector<double> yLocalPolygon(inputNodes.size());
            for (size_t i = 0; i < inputNodes.size(); ++i)
            {
                xLocalPolygon[i] = inputNodes[i].x;
                yLocalPolygon[i] = inputNodes[i].y;
            }

            int numFaces = -1;
            int numEdges = 0;
            int numNodes = 0;

            int numInputNodes = triangulationOption == TriangulationOptions::TriangulatePointsAndGenerateFaces ? static_cast<int>(inputNodes.size()) : static_cast<int>(inputNodes.size() - 1);
            auto intTriangulationOption = static_cast<int>(triangulationOption);

            std::vector<int> faceNodesFlat;
            std::vector<int> edgeNodesFlat;
            std::vector<int> faceEdgesFlat;
            std::vector<double> xNodesFlat;
            std::vector<double> yNodesFlat;

            if (estimatedNumberOfTriangles == 0)
            {
                estimatedNumberOfTriangles = inputNodes.size() * 6 + 10;
            }

            // If the number of estimated triangles is not sufficient, triangulation must be repeated
            while (numFaces < 0)
            {
                numFaces = static_cast<int>(estimatedNumberOfTriangles);
                faceNodesFlat.resize(numFaces * 3);
                edgeNodesFlat.resize(numFaces * 2);
                faceEdgesFlat.resize(numFaces * 3);
                xNodesFlat.resize(numFaces * 3, doubleMissingValue);
                yNodesFlat.resize(numFaces * 3, doubleMissingValue);
                Triangulation(&intTriangulationOption,
                              &xLocalPolygon[0],
                              &yLocalPolygon[0],
                              &numInputNodes,
                              &faceNodesFlat[0], // INDX
                              &numFaces,
                              &edgeNodesFlat[0], // EDGEINDX
                              &numEdges,
                              &faceEdgesFlat[0], // TRIEDGE
                              &xNodesFlat[0],
                              &yNodesFlat[0],
                              &numNodes,
                              &averageTriangleArea);
                if (estimatedNumberOfTriangles > 0)
                {
                    estimatedNumberOfTriangles = -numFaces;
                }
            }

            m_numFaces = numFaces <= 0 ? static_cast<size_t>(0) : static_cast<size_t>(numFaces);
            m_numEdges = static_cast<size_t>(numEdges);
            m_numNodes = static_cast<size_t>(numNodes);

            // Create nodes
            m_nodes.resize(m_numNodes);
            for (size_t i = 0; i < m_numNodes; ++i)
            {
                m_nodes[i] = {xNodesFlat[i], yNodesFlat[i]};
            }

            // Create m_faceNodes
            ResizeAndFill2DVector(m_faceNodes, m_numFaces, 3, true, sizetMissingValue);
            ResizeAndFill2DVector(m_faceEdges, m_numFaces, 3, true, sizetMissingValue);
            size_t faceCounter = 0;
            for (size_t f = 0; f < m_numFaces; ++f)
            {
                m_faceNodes[f][0] = static_cast<size_t>(faceNodesFlat[faceCounter] - 1);
                m_faceEdges[f][0] = static_cast<size_t>(faceEdgesFlat[faceCounter] - 1);
                faceCounter++;
                m_faceNodes[f][1] = static_cast<size_t>(faceNodesFlat[faceCounter] - 1);
                m_faceEdges[f][1] = static_cast<size_t>(faceEdgesFlat[faceCounter] - 1);
                faceCounter++;
                m_faceNodes[f][2] = static_cast<size_t>(faceNodesFlat[faceCounter] - 1);
                m_faceEdges[f][2] = static_cast<size_t>(faceEdgesFlat[faceCounter] - 1);
                faceCounter++;
            }

            // Create edges
            if (m_numEdges == 0)
            {
                return;
            }

            ResizeAndFill2DVector(m_edgeNodes, m_numEdges, 2, true, sizetMissingValue);
            size_t edgeCounter = 0;
            for (size_t e = 0; e < m_numEdges; ++e)
            {
                m_edgeNodes[e][0] = static_cast<size_t>(edgeNodesFlat[edgeCounter] - 1);
                edgeCounter++;
                m_edgeNodes[e][1] = static_cast<size_t>(edgeNodesFlat[edgeCounter] - 1);
                edgeCounter++;
            }

            ResizeAndFill2DVector(m_edgesFaces, m_numEdges, 2, true, sizetMissingValue);
            edgeCounter = 0;
            for (size_t f = 0; f < m_numFaces; ++f)
            {

                for (size_t n = 0; n < numNodesInTriangle; ++n)
                {
                    auto const edge = static_cast<size_t>(faceEdgesFlat[edgeCounter] - 1);
                    edgeCounter++;
                    // For each edge, the shared face index
                    if (m_edgesFaces[edge][0] == sizetMissingValue)
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
