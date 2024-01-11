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
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Point.hpp"
#include "MeshKernel/PolygonalEnclosure.hpp"

#include <concepts>

namespace meshkernel
{
    extern "C"
    {
        /// @brief Function of the Triangle library
        ///
        /// \see https://www.cs.cmu.edu/~quake/triangle.html
        void Triangulation(int jatri,
                           double const* const xs,
                           double const* const ys,
                           int ns,
                           int* const indx,
                           int* const numtri,
                           int* const edgeidx,
                           int* const numedge,
                           int* const triedge,
                           double* const xs3,
                           double* const ys3,
                           int* const ns3,
                           double trisize);
    }

    class Sample;

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
        template <std::derived_from<Point> T>
        void Compute(const std::vector<T>& inputNodes,
                     TriangulationOptions triangulationOption,
                     double averageTriangleArea,
                     UInt estimatedNumberOfTriangles)
        {
            if (inputNodes.empty())
            {
                throw ConstraintError("The sample is empty.");
            }
            std::vector<double> xLocalPolygon(inputNodes.size());
            std::vector<double> yLocalPolygon(inputNodes.size());
            for (UInt i = 0; i < inputNodes.size(); ++i)
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
                estimatedNumberOfTriangles = static_cast<UInt>(inputNodes.size()) * 6 + 10;
            }

            // If the number of estimated triangles is not sufficient, triangulation must be repeated
            while (m_numFaces < 0)
            {
                m_numFaces = static_cast<int>(estimatedNumberOfTriangles);

                m_faceNodesFlat.resize(estimatedNumberOfTriangles * 3);
                std::ranges::fill(m_faceNodesFlat, 0);

                m_edgeNodesFlat.resize(estimatedNumberOfTriangles * 2);
                std::ranges::fill(m_edgeNodesFlat, 0);

                m_faceEdgesFlat.resize(estimatedNumberOfTriangles * 3);
                std::ranges::fill(m_faceEdgesFlat, 0);

                m_xCoordFlat.resize(estimatedNumberOfTriangles * 3, constants::missing::doubleValue);
                std::ranges::fill(m_xCoordFlat, 0.0);

                m_yCoordFlat.resize(estimatedNumberOfTriangles * 3, constants::missing::doubleValue);
                std::ranges::fill(m_yCoordFlat, 0.0);

                Triangulation(intTriangulationOption,
                              xLocalPolygon.data(),
                              yLocalPolygon.data(),
                              numInputNodes,
                              m_faceNodesFlat.data(), // INDX
                              &m_numFaces,
                              m_edgeNodesFlat.data(), // EDGEINDX
                              &m_numEdges,
                              m_faceEdgesFlat.data(), // TRIEDGE
                              m_xCoordFlat.data(),
                              m_yCoordFlat.data(),
                              &m_numNodes,
                              averageTriangleArea);
                if (estimatedNumberOfTriangles > 0)
                {
                    estimatedNumberOfTriangles = -m_numFaces;
                }
            }
        }

        /// @brief From the set of computed points, select those that are contained within the enclosure
        std::vector<Point> SelectNodes(const PolygonalEnclosure& enclosure) const;

        /// @brief Build the internal triangulation from the flat triangulation
        void BuildTriangulation();

        /// @brief Gets the number of triangulated edges
        /// @return The number of triangulated edges
        [[nodiscard]] int GetNumEdges() const
        {
            return m_numEdges;
        }

        /// @brief Gets the number of triangulated nodes
        /// @return The number of triangulated nodes
        [[nodiscard]] int GetNumNodes() const
        {
            return m_numNodes;
        }

        /// @brief Gets the number of triangulated faces
        /// @return The number of triangulated faces
        [[nodiscard]] int GetNumFaces() const
        {
            return m_numFaces;
        }

        /// @brief Gets the triangulated nodes
        /// @return The triangulated nodes
        [[nodiscard]] const std::vector<Point>& GetNodes() const
        {
            return m_nodes;
        }

        /// @brief Gets the nodes of a triangulated face
        /// @param faceIndex The face index
        /// @return The triangulated nodes
        [[nodiscard]] const std::vector<UInt>& GetFaceNodes(const UInt faceIndex) const
        {
            return m_faceNodes[faceIndex];
        }

        /// @brief Retrieves the face node
        /// @param faceIndex The index of the face to retrieve the node from
        /// @param nodeIndex The index of the node to retrieve
        /// @return const reference to the node with the specified index for the specified face
        [[nodiscard]] UInt GetFaceNode(const UInt faceIndex, const UInt nodeIndex) const
        {
            return m_faceNodes[faceIndex][nodeIndex];
        }

        /// @brief Retrieves the face edge
        /// @param faceIndex The index of the face to retrieve the edge from
        /// @param edgeIndex The index of the edge to retrieve
        /// @return const reference to the edge with the specified index for the specified face
        [[nodiscard]] UInt GetFaceEdge(const UInt faceIndex, const UInt edgeIndex) const
        {
            return m_faceEdges[faceIndex][edgeIndex];
        }

        /// @brief Retrieves the edge node
        /// @param edgeIndex The index of the edge to retrieve the node from
        /// @param nodeIndex The index of the node to retrieve
        /// @return const reference to the node with the specified index for the specified face
        [[nodiscard]] UInt GetEdgeNode(const UInt edgeIndex, const UInt nodeIndex) const
        {
            return m_edgeNodes[edgeIndex][nodeIndex];
        }

        /// @brief Retrieves the edge face
        /// @param edgeIndex The index of the edge to retrieve the node from
        /// @param faceIndex The index of the face to retrieve
        /// @return const reference to the edge with the specified index for the specified face
        [[nodiscard]] UInt GetEdgeFace(const UInt edgeIndex, const UInt faceIndex) const
        {
            return m_edgesFaces[edgeIndex][faceIndex];
        }

        /// @brief Retrieves the x coordinate of a triangulated node
        /// @param nodeIndex The index of the node to retrieve
        /// @return const reference to the x coordinate
        [[nodiscard]] double GetXCoord(const UInt nodeIndex) const
        {
            return m_xCoordFlat[nodeIndex];
        }

        /// @brief Retrieves the y coordinate of a triangulated node
        /// @param nodeIndex The index of the node to retrieve
        /// @return const reference to the y coordinate
        [[nodiscard]] double GetYCoord(const UInt nodeIndex) const
        {
            return m_yCoordFlat[nodeIndex];
        }

        /// @brief Retrieves the (x,y) coordinate of a triangulated node
        /// @param nodeIndex The index of the node to retrieve
        /// @return Point
        Point GetCoord(const UInt nodeIndex) const
        {
            return Point(m_xCoordFlat[nodeIndex], m_yCoordFlat[nodeIndex]);
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

        std::vector<Point> m_nodes;                  ///< Reconstructed vector of nodes
        std::vector<std::vector<UInt>> m_faceNodes;  ///< Reconstructed vector of face nodes
        std::vector<std::vector<UInt>> m_faceEdges;  ///< Reconstructed vector of face edges
        std::vector<std::vector<UInt>> m_edgeNodes;  ///< Reconstructed vector of edge nodes
        std::vector<std::vector<UInt>> m_edgesFaces; ///< Reconstructed vector of edge faces
    };

} // namespace meshkernel
