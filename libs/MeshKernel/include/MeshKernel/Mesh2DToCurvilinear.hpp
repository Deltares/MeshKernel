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

#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"
#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/Point.hpp"
#include "Utilities/LinearAlgebra.hpp"

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

        std::unique_ptr<CurvilinearGrid> Compute(const Point& point)
        {
            m_mesh.SearchNearestLocation(point, Location::Faces);
            if (m_mesh.GetNumLocations(Location::Faces) <= 0)
            {
                throw AlgorithmError("No valid face found close to the initial point.");
            }

            const auto initialFace = m_mesh.GetLocationsIndices(0, Location::Faces);
            if (m_mesh.GetNumFaceEdges(initialFace) != geometric::numNodesInQuadrilateral)
            {
                throw AlgorithmError("The initial face is not a quadrilateral");
            }

            // Check the seed point is contained
            std::vector<Point> polygonPoints;
            for (UInt n = 0; n < geometric::numNodesInQuadrilateral; ++n)
            {
                const auto node = m_mesh.m_facesNodes[initialFace][n];
                polygonPoints.emplace_back(m_mesh.Node(node));
            }
            const auto node = m_mesh.m_facesNodes[initialFace][0];
            polygonPoints.emplace_back(m_mesh.Node(node));

            Polygon polygon(polygonPoints, m_mesh.m_projection);
            if (!polygon.Contains(point))
            {
                throw AlgorithmError("The initial face does not contain the starting point");
            }

            // build local coordinate system
            const auto numNodes = m_mesh.GetNumNodes();
            m_i = std::vector(numNodes, missing::intValue);
            m_j = std::vector(numNodes, missing::intValue);

            const auto firstEdge = m_mesh.m_facesEdges[initialFace][0];
            const auto secondEdge = m_mesh.m_facesEdges[initialFace][1];
            const auto thirdEdge = m_mesh.m_facesEdges[initialFace][2];
            const auto fourthEdge = m_mesh.m_facesEdges[initialFace][3];

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
            q.push(initialFace);

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

                const auto localNodeMapping = ComputeLocalNodeMapping(face);
                for (UInt d = 0; d < m_numDirections; ++d)
                {
                    const auto newFaceIndex = ComputeNeighbouringFaces(face, localNodeMapping, d, visitedFace);
                    if (newFaceIndex != missing::uintValue)
                    {
                        q.push(newFaceIndex);
                    }
                }
            }

            const auto matrix = ComputeMatrix();
            return std::make_unique<CurvilinearGrid>(matrix, m_mesh.m_projection);
        }

    private:
        Eigen::Matrix<UInt, 2, 2> ComputeLocalNodeMapping(UInt face)
        {
            const auto& faceIndices = m_mesh.m_facesNodes[face];

            const auto node0 = faceIndices[0];
            const auto node1 = faceIndices[1];
            const auto node2 = faceIndices[2];
            const auto node3 = faceIndices[3];

            const std::vector localI{m_i[node0], m_i[node1], m_i[node2], m_i[node3]};
            const std::vector localJ{m_j[node0], m_j[node1], m_j[node2], m_j[node3]};
            const std::vector localNodes{node0, node1, node2, node3};

            const auto minI = *std::ranges::min_element(localI);
            const auto minJ = *std::ranges::min_element(localJ);

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
            return matrix;
        }

        UInt ComputeNeighbouringFaces(const UInt face,
                                      const Eigen::Matrix<UInt, 2, 2>& localNodeMapping,
                                      const UInt d,
                                      std::vector<bool>& visitedFace)
        {

            const auto firstNode = localNodeMapping(m_nodeFrom[d][0], m_nodeFrom[d][1]);
            const auto secondNode = localNodeMapping(m_nodeTo[d][0], m_nodeTo[d][1]);

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

        lin_alg::Matrix<Point> ComputeMatrix()
        {

            int minI = n_maxNumRowsColumns;
            int minJ = n_maxNumRowsColumns;
            const auto numNodes = m_mesh.GetNumNodes();
            for (UInt n = 0; n < numNodes; ++n)
            {
                if (m_i[n] != missing::intValue)
                {
                    minI = std::min(minI, m_i[n]);
                }
                if (m_j[n] != missing::intValue)
                {
                    minJ = std::min(minJ, m_j[n]);
                }
            }

            int maxI = -n_maxNumRowsColumns;
            int maxJ = -n_maxNumRowsColumns;
            for (UInt n = 0; n < numNodes; ++n)
            {
                if (m_i[n] != missing::intValue)
                {
                    m_i[n] = m_i[n] + minI;
                    maxI = std::max(maxI, m_i[n]);
                }
                if (m_j[n] != missing::intValue)
                {
                    m_j[n] = m_j[n] + minJ;
                    maxJ = std::max(maxJ, m_j[n]);
                }
            }

            lin_alg::Matrix<Point> result(maxI + 1, maxJ + 1);
            for (UInt n = 0; n < numNodes; ++n)
            {
                const auto i = m_i[n];
                const auto j = m_j[n];
                if (i != missing::intValue && j != missing::intValue)
                {
                    result(i, j) = m_mesh.Node(n);
                }
                else
                {
                    result(i, j) = Point(missing::doubleValue,
                                         missing::doubleValue);
                }
            }
            return result;
        }

        Mesh2D& m_mesh; ///< The mesh where the edges should be found

        std::vector<int> m_i;
        std::vector<int> m_j;
        UInt const m_numDirections = 4;

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

        int n_maxNumRowsColumns = 1000000;
    };

} // namespace meshkernel
