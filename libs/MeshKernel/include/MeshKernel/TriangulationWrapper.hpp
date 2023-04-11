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

#include "MeshKernel/Constants.hpp"

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

        /// @brief Compute the triangulation
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

            m_numFaces = -1;
            m_numEdges = 0;
            m_numNodes = 0;

            int numInputNodes = static_cast<int>(inputNodes.size());
            auto intTriangulationOption = static_cast<int>(triangulationOption);

            if (estimatedNumberOfTriangles == 0)
            {
                estimatedNumberOfTriangles = inputNodes.size() * 6 + 10;
            }

            // If the number of estimated triangles is not sufficient, triangulation must be repeated
            while (m_numFaces < 0)
            {
                m_numFaces = static_cast<int>(estimatedNumberOfTriangles);

                m_faceNodesFlat.resize(estimatedNumberOfTriangles * 3);
                std::fill(m_faceNodesFlat.begin(), m_faceNodesFlat.end(), 0);

                m_edgeNodesFlat.resize(estimatedNumberOfTriangles * 2);
                std::fill(m_edgeNodesFlat.begin(), m_edgeNodesFlat.end(), 0);

                m_faceEdgesFlat.resize(estimatedNumberOfTriangles * 3);
                std::fill(m_faceEdgesFlat.begin(), m_faceEdgesFlat.end(), 0);

                m_xCoordFlat.resize(estimatedNumberOfTriangles * 3, constants::missing::doubleValue);
                std::fill(m_xCoordFlat.begin(), m_xCoordFlat.end(), 0.0);

                m_yCoordFlat.resize(estimatedNumberOfTriangles * 3, constants::missing::doubleValue);
                std::fill(m_yCoordFlat.begin(), m_yCoordFlat.end(), 0.0);

                Triangulation(&intTriangulationOption,
                              xLocalPolygon.data(),
                              yLocalPolygon.data(),
                              &numInputNodes,
                              m_faceNodesFlat.data(), // INDX
                              &m_numFaces,
                              m_edgeNodesFlat.data(), // EDGEINDX
                              &m_numEdges,
                              m_faceEdgesFlat.data(), // TRIEDGE
                              m_xCoordFlat.data(),
                              m_yCoordFlat.data(),
                              &m_numNodes,
                              &averageTriangleArea);
                if (estimatedNumberOfTriangles > 0)
                {
                    estimatedNumberOfTriangles = -m_numFaces;
                }
            }
        }

        /// @brief Build the internal triangulation from the flat triangulation
        void BuildTriangulation();

        /// @brief Gets the number of triangulated edges
        /// @return The number of triangulated edges
        [[nodiscard]] const auto& GetNumEdges() const
        {
            return m_numEdges;
        }

        /// @brief Gets the number of triangulated nodes
        /// @return The number of triangulated nodes
        [[nodiscard]] const auto& GetNumNodes() const
        {
            return m_numNodes;
        }

        /// @brief Gets the number of triangulated faces
        /// @return The number of triangulated faces
        [[nodiscard]] const auto& GetNumFaces() const
        {
            return m_numFaces;
        }

        /// @brief Gets the triangulated nodes
        /// @return The triangulated nodes
        [[nodiscard]] const auto& GetNodes() const
        {
            return m_nodes;
        }

        /// @brief Gets the nodes of a triangulated face
        /// @param faceIndex The face index
        /// @return The triangulated nodes
        [[nodiscard]] const auto& GetFaceNodes(const size_t& faceIndex) const
        {
            return m_faceNodes[faceIndex];
        }

        /// @brief Retrieves the face node
        /// @param faceIndex The index of the face to retrieve the node from
        /// @param nodeIndex The index of the node to retrieve
        /// @return const reference to the node with the specified index for the specified face
        [[nodiscard]] const auto& GetFaceNode(const size_t& faceIndex, const size_t& nodeIndex) const
        {
            return m_faceNodes[faceIndex][nodeIndex];
        }

        /// @brief Retrieves the face edge
        /// @param faceIndex The index of the face to retrieve the edge from
        /// @param edgeIndex The index of the edge to retrieve
        /// @return const reference to the edge with the specified index for the specified face
        [[nodiscard]] const auto& GetFaceEdge(const size_t& faceIndex, const size_t& edgeIndex) const
        {
            return m_faceEdges[faceIndex][edgeIndex];
        }

        /// @brief Retrieves the edge node
        /// @param edgeIndex The index of the edge to retrieve the node from
        /// @param nodeIndex The index of the node to retrieve
        /// @return const reference to the node with the specified index for the specified face
        [[nodiscard]] const auto& GetEdgeNode(const size_t& edgeIndex, const size_t& nodeIndex) const
        {
            return m_edgeNodes[edgeIndex][nodeIndex];
        }

        /// @brief Retrieves the edge face
        /// @param edgeIndex The index of the edge to retrieve the node from
        /// @param faceIndex The index of the face to retrieve
        /// @return const reference to the edge with the specified index for the specified face
        [[nodiscard]] const auto& GetEdgeFace(const size_t& edgeIndex, const size_t& faceIndex) const
        {
            return m_edgesFaces[edgeIndex][faceIndex];
        }

        /// @brief Retrieves the x coordinate of a triangulated node
        /// @param nodeIndex The index of the node to retrieve
        /// @return const reference to the x coordinate
        [[nodiscard]] const auto& GetXCoord(const size_t& nodeIndex) const
        {
            return m_xCoordFlat[nodeIndex];
        }

        /// @brief Retrieves the y coordinate of a triangulated node
        /// @param nodeIndex The index of the node to retrieve
        /// @return const reference to the y coordinate
        [[nodiscard]] const auto& GetYCoord(const size_t& nodeIndex) const
        {
            return m_yCoordFlat[nodeIndex];
        }

    private:

        std::vector<int> m_faceNodesFlat; ///< Face nodes flat array passed to the triangulation library
        std::vector<int> m_edgeNodesFlat; ///< Edge nodes flat array passed to the triangulation library
        std::vector<int> m_faceEdgesFlat; ///< Face edges flat array passed to the triangulation library
        std::vector<double> m_xCoordFlat; ///< x coordinates flat array passed to the triangulation library
        std::vector<double> m_yCoordFlat; ///< y coordinates flat array passed to the triangulation library
        int m_numNodes{0};                ///< Initial number of triangulated nodes
        int m_numEdges{0};                ///< Initial number of triangulated edges
        int m_numFaces{0};                ///< Initial number of triangulated faces

        std::vector<Point> m_nodes;                    ///< Reconstructed vector of nodes
        std::vector<std::vector<size_t>> m_faceNodes;  ///< Reconstructed vector of face nodes
        std::vector<std::vector<size_t>> m_faceEdges;  ///< Reconstructed vector of face edges
        std::vector<std::vector<size_t>> m_edgeNodes;  ///< Reconstructed vector of edge nodes
        std::vector<std::vector<size_t>> m_edgesFaces; ///< Reconstructed vector of edge faces
    };

} // namespace meshkernel
