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

            numFaces = -1;
            numEdges = 0;
            numNodes = 0;

            int numInputNodes = static_cast<int>(inputNodes.size());
            auto intTriangulationOption = static_cast<int>(triangulationOption);

            if (estimatedNumberOfTriangles == 0)
            {
                estimatedNumberOfTriangles = inputNodes.size() * 6 + 10;
            }

            // If the number of estimated triangles is not sufficient, triangulation must be repeated
            while (numFaces < 0)
            {
                numFaces = static_cast<int>(estimatedNumberOfTriangles);

                faceNodesFlat.resize(estimatedNumberOfTriangles * 3);
                std::fill(faceNodesFlat.begin(), faceNodesFlat.end(), 0);

                edgeNodesFlat.resize(estimatedNumberOfTriangles * 2);
                std::fill(edgeNodesFlat.begin(), edgeNodesFlat.end(), 0);

                faceEdgesFlat.resize(estimatedNumberOfTriangles * 3);
                std::fill(faceEdgesFlat.begin(), faceEdgesFlat.end(), 0);

                xCoordFlat.resize(estimatedNumberOfTriangles * 3, constants::missing::doubleValue);
                std::fill(xCoordFlat.begin(), xCoordFlat.end(), 0.0);

                yCoordFlat.resize(estimatedNumberOfTriangles * 3, constants::missing::doubleValue);
                std::fill(yCoordFlat.begin(), yCoordFlat.end(), 0.0);

                Triangulation(&intTriangulationOption,
                              &xLocalPolygon[0],
                              &yLocalPolygon[0],
                              &numInputNodes,
                              &faceNodesFlat[0], // INDX
                              &numFaces,
                              &edgeNodesFlat[0], // EDGEINDX
                              &numEdges,
                              &faceEdgesFlat[0], // TRIEDGE
                              &xCoordFlat[0],
                              &yCoordFlat[0],
                              &numNodes,
                              &averageTriangleArea);
                if (estimatedNumberOfTriangles > 0)
                {
                    estimatedNumberOfTriangles = -numFaces;
                }
            }
        }

        /// @brief Build internal mapping after the flat triangulation structures are filled
        void BuildTriangulation()
        {

            numFaces = numFaces <= 0 ? 0 : numFaces;

            // Create nodes
            m_nodes.resize(numNodes);
            for (auto i = 0; i < numNodes; ++i)
            {
                m_nodes[i] = {xCoordFlat[i], yCoordFlat[i]};
            }

            // Create m_faceNodes
            ResizeAndFill2DVector(m_faceNodes, numFaces, 3, true, constants::missing::sizetValue);
            ResizeAndFill2DVector(m_faceEdges, numFaces, 3, true, constants::missing::sizetValue);
            size_t faceCounter = 0;
            for (size_t f = 0; f < numFaces; ++f)
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
            if (numEdges == 0)
            {
                return;
            }

            ResizeAndFill2DVector(m_edgeNodes, numEdges, 2, true, constants::missing::sizetValue);
            size_t edgeCounter = 0;
            for (size_t e = 0; e < numEdges; ++e)
            {
                m_edgeNodes[e][0] = static_cast<size_t>(edgeNodesFlat[edgeCounter] - 1);
                edgeCounter++;
                m_edgeNodes[e][1] = static_cast<size_t>(edgeNodesFlat[edgeCounter] - 1);
                edgeCounter++;
            }

            ResizeAndFill2DVector(m_edgesFaces, numEdges, 2, true, constants::missing::sizetValue);
            edgeCounter = 0;
            for (size_t f = 0; f < numFaces; ++f)
            {

                for (size_t n = 0; n < Mesh::m_numNodesInTriangle; ++n)
                {
                    auto const edge = static_cast<size_t>(faceEdgesFlat[edgeCounter] - 1);
                    edgeCounter++;
                    // For each edge, the shared face index
                    if (m_edgesFaces[edge][0] == constants::missing::sizetValue)
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

        /// @brief Gets the number of triangulated edges
        /// @return The number of triangulated edges
        [[nodiscard]] const auto& GetNumEdges() const
        {
            return numEdges;
        }

        /// @brief Gets the number of triangulated nodes
        /// @return The number of triangulated nodes
        [[nodiscard]] const auto& GetNumNodes() const
        {
            return numNodes;
        }

        /// @brief Gets the number of triangulated faces
        /// @return The number of triangulated faces
        [[nodiscard]] const auto& GetNumFaces() const
        {
            return numFaces;
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
            return xCoordFlat[nodeIndex];
        }

        /// @brief Retrieves the y coordinate of a triangulated node
        /// @param nodeIndex The index of the node to retrieve
        /// @return const reference to the y coordinate
        [[nodiscard]] const auto& GetYCoord(const size_t& nodeIndex) const
        {
            return yCoordFlat[nodeIndex];
        }

    private:
        std::vector<int> faceNodesFlat; ///< Face nodes flat array passed to the triangulation library
        std::vector<int> edgeNodesFlat; ///< Edge nodes flat array passed to the triangulation library
        std::vector<int> faceEdgesFlat; ///< Face edges flat array passed to the triangulation library
        std::vector<double> xCoordFlat; ///< x coordinates flat array passed to the triangulation library
        std::vector<double> yCoordFlat; ///< y coordinates flat array passed to the triangulation library
        int numNodes{0};                ///< Initial number of triangulated nodes
        int numEdges{0};                ///< Initial number of triangulated edges
        int numFaces{0};                ///< Initial number of triangulated faces

        std::vector<Point> m_nodes;                    ///< Reconstructed vector of nodes
        std::vector<std::vector<size_t>> m_faceNodes;  ///< Reconstructed vector of face nodes
        std::vector<std::vector<size_t>> m_faceEdges;  ///< Reconstructed vector of face edges
        std::vector<std::vector<size_t>> m_edgeNodes;  ///< Reconstructed vector of edge nodes
        std::vector<std::vector<size_t>> m_edgesFaces; ///< Reconstructed vector of edge faces
    };

} // namespace meshkernel
