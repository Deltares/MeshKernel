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
#include <queue>

#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/Point.hpp"

using namespace meshkernel::constants;

namespace meshkernel
{
    /// @brief Construct a curvilinear grid from an unstructured mesh
    class Mesh2DToCurvilinear
    {
    public:
        /// @brief Constructor
        explicit Mesh2DToCurvilinear(Mesh2D& mesh) : m_mesh(mesh)
        {
            if (mesh.GetNumNodes() <= 0)
            {
                throw AlgorithmError("Mesh with no nodes");
            }

            mesh.Administrate();

            if (mesh.GetNumFaces() <= 0)
            {
                throw AlgorithmError("Mesh with no faces");
            }
        }

        void Compute(const Point& point)
        {
            m_mesh.SearchNearestLocation(point, Location::Faces);
            if (m_mesh.GetNumLocations(Location::Faces) <= 0)
            {
                return;
            }

            const auto frontFace = m_mesh.GetLocationsIndices(0, Location::Faces);
            if (m_mesh.GetNumFaceEdges(frontFace) != geometric::numNodesInQuadrilateral)
            {
                return;
            }

            // Check the seed point is contained
            std::vector<Point> polygonPoints;
            for (UInt n = 0; n < geometric::numNodesInQuadrilateral; ++n)
            {
                const auto node = m_mesh.m_facesNodes[frontFace][n];
                polygonPoints.emplace_back(m_mesh.Node(node));
            }
            const auto node = m_mesh.m_facesNodes[frontFace][0];
            polygonPoints.emplace_back(m_mesh.Node(node));

            Polygon polygon(polygonPoints, m_mesh.m_projection);
            if (!polygon.Contains(point))
            {
                return;
            }

            // build local coordinate system
            const auto numNodes = m_mesh.GetNumNodes();
            m_i = std::vector(numNodes, missing::intValue);
            m_j = std::vector(numNodes, missing::intValue);

            const auto firstEdge = m_mesh.m_facesEdges[frontFace][0];
            const auto secondEdge = m_mesh.m_facesEdges[frontFace][1];
            const auto thirdEdge = m_mesh.m_facesEdges[frontFace][2];
            const auto fourthEdge = m_mesh.m_facesEdges[frontFace][3];

            const auto firstNodeIndex = m_mesh.FindCommonNode(firstEdge, secondEdge);
            m_i[firstNodeIndex] = 0;
            m_j[firstNodeIndex] = 0;

            const auto secondNodeIndex = m_mesh.FindCommonNode(secondEdge, thirdEdge);
            m_i[secondNodeIndex] = 1;
            m_j[secondNodeIndex] = 0;

            const auto thirdNodeIndex = m_mesh.FindCommonNode(thirdEdge, fourthEdge);
            m_i[thirdNodeIndex] = 1;
            m_j[thirdNodeIndex] = 1;

            const auto fourthNodeIndex = m_mesh.FindCommonNode(fourthEdge, firstEdge);
            m_i[fourthNodeIndex] = 0;
            m_j[fourthNodeIndex] = 1;

            const auto numFaces = m_mesh.GetNumFaces();
            std::vector visitedFace(numFaces, false);

            std::queue<UInt> q;
            q.push(frontFace);

            while (!q.empty())
            {
                const auto face = q.front();
                q.pop();
                if (visitedFace[face])
                {
                    continue;
                }
                visitedFace[face] = true;

                if (m_mesh.GetNumFaceEdges(face) != geometric::numNodesInQuadrilateral)
                {
                    continue;
                }

                const auto& faceIndices = m_mesh.m_facesNodes[face];

                const auto node0 = faceIndices[0];
                const auto node1 = faceIndices[1];
                const auto node2 = faceIndices[2];
                const auto node3 = faceIndices[3];

                const std::vector localI{m_i[node0], m_i[node1], m_i[node2], m_i[node3]};
                const std::vector localJ{m_j[node0], m_j[node1], m_j[node2], m_j[node3]};
                const std::vector localNodes{node0, node1, node2, node3};

                auto minI = *std::ranges::min_element(localI);
                auto minJ = *std::ranges::min_element(localJ);

                // build a local mapping
                Eigen::Matrix<UInt, 2, 2> matrix;
                for (UInt i = 0; i < localI.size(); ++i)
                {
                    if (localI[i] == minI && localJ[i] == minJ)
                    {
                        matrix(0, 0) = localNodes[i];
                    }
                    if (localI[i] == minI + 1 && localJ[i] == minJ)
                    {
                        matrix(1, 0) = localNodes[i];
                    }
                    if (localI[i] == minI && localJ[i] == minJ + 1)
                    {
                        matrix(0, 1) = localNodes[i];
                    }
                    if (localI[i] == minI + 1 && localJ[i] == minJ + 1)
                    {
                        matrix(1, 1) = localNodes[i];
                    }
                }

                for (UInt d = 0; d < numDirections; ++d)
                {
                    const auto newFaceIndex = computeNeighbourFace(face, d, matrix, visitedFace);
                    if (newFaceIndex != missing::uintValue)
                    {
                        q.push(newFaceIndex);
                    }
                }
            }
        }

    private:
        UInt computeNeighbourFace(const UInt face,
                                  const UInt d,
                                  const Eigen::Matrix2i& matrix,
                                  std::vector<bool>& visitedFace)
        {
            const auto firstNode = matrix(m_nodeFrom[d][0], m_nodeFrom[d][1]);
            const auto secondNode = matrix(m_nodeTo[d][0], m_nodeTo[d][1]);

            // find the edge index
            const auto edgeIndex = m_mesh.FindEdge(firstNode, secondNode);
            if (edgeIndex == missing::uintValue)
            {
                return missing::uintValue;
            }

            // this edge belongs only to the current face
            if (m_mesh.m_edgesNumFaces[edgeIndex] < 2)
            {
                return missing::uintValue;
            }

            const auto newFace = face == m_mesh.m_edgesFaces[edgeIndex][0] ? m_mesh.m_edgesFaces[edgeIndex][1] : m_mesh.m_edgesFaces[edgeIndex][0];

            if (visitedFace[newFace])
            {
                return missing::uintValue;
            }

            if (m_mesh.GetNumFaceEdges(newFace) != geometric::numNodesInQuadrilateral)
            {
                return missing::uintValue;
            }

            int edgeIndexInNewFace = 0;
            for (UInt e = 0; e < m_mesh.GetNumFaceEdges(newFace); ++e)
            {
                if (m_mesh.m_facesEdges[newFace][e] == edgeIndex)
                {
                    edgeIndexInNewFace = e;
                    break;
                }
            }
            auto nextEdgeIndexInNewFace = edgeIndexInNewFace + 1;
            nextEdgeIndexInNewFace = nextEdgeIndexInNewFace > 3 ? 0 : nextEdgeIndexInNewFace;
            const auto nextEdgeInNewFace = m_mesh.m_facesEdges[newFace][nextEdgeIndexInNewFace];
            const auto firstCommonNode = m_mesh.FindCommonNode(edgeIndex, nextEdgeInNewFace);
            const auto i_firstCommonNode = m_i[firstCommonNode] + m_directionsDeltas[d][0];
            const auto j_firstCommonNode = m_j[firstCommonNode] + m_directionsDeltas[d][1];

            auto previousEdgeIndexInNewFace = edgeIndexInNewFace - 1;
            previousEdgeIndexInNewFace = previousEdgeIndexInNewFace < 0 ? 3 : previousEdgeIndexInNewFace;
            const auto previousEdgeInNewFace = m_mesh.m_facesEdges[newFace][previousEdgeIndexInNewFace];
            const auto secondCommonNode = m_mesh.FindCommonNode(edgeIndex, previousEdgeInNewFace);
            const auto i_secondCommonNode = m_i[secondCommonNode] + m_directionsDeltas[d][0];
            const auto j_secondCommonNode = m_j[secondCommonNode] + m_directionsDeltas[d][1];

            const auto invalid = m_i[firstCommonNode] != missing::intValue && m_i[firstCommonNode] != i_firstCommonNode ||
                                 m_j[firstCommonNode] != missing::intValue && m_j[firstCommonNode] != j_firstCommonNode ||
                                 m_i[secondCommonNode] != missing::intValue && m_i[secondCommonNode] != i_secondCommonNode ||
                                 m_j[secondCommonNode] != missing::intValue && m_j[secondCommonNode] != j_secondCommonNode;

            if (!invalid)
            {
                m_i[firstCommonNode] = i_firstCommonNode;
                m_j[firstCommonNode] = j_firstCommonNode;
                m_i[secondCommonNode] = i_secondCommonNode;
                m_j[secondCommonNode] = j_secondCommonNode;
                return newFace;
            }
            return missing::uintValue;
        }

        void ComputeMatrix()
        {
        }

        Mesh2D& m_mesh; ///< The mesh where the edges should be found

        std::vector<int> m_i;
        std::vector<int> m_j;
        int const numDirections = 4;
        std::array<int, 4> m_dirI{0, 1, 1, 0};
        std::array<int, 4> m_dirJ{1, 0, 1, 1};

        std::array<std::array<int, 2>, 4> m_nodeFrom = {{{0, 0},
                                                         {0, 0},
                                                         {1, 0},
                                                         {1, 1}}};

        std::array<std::array<int, 2>, 4> m_nodeTo = {{{0, 0},
                                                       {0, 0},
                                                       {1, 0},
                                                       {1, 1}}};

        std::array<std::array<int, 2>, 4> m_directionsDeltas = {{{-1, 0},
                                                                 {0, -1},
                                                                 {1, 0},
                                                                 {0, 1}}};
    };

} // namespace meshkernel
