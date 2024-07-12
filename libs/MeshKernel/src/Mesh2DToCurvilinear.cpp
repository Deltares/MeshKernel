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

#include <queue>

#include "MeshKernel/Mesh2DToCurvilinear.hpp"

using namespace meshkernel;
using namespace constants;

Mesh2DToCurvilinear::Mesh2DToCurvilinear(Mesh2D& mesh) : m_mesh(mesh)
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

std::unique_ptr<CurvilinearGrid> Mesh2DToCurvilinear::Compute(const Point& point)
{
    // 1. Find the face index
    m_mesh.BuildTree(Location::Faces);
    const auto initialFaceIndex = m_mesh.FindLocationIndex(point, Location::Faces);
    if (m_mesh.GetNumFaceEdges(initialFaceIndex) != geometric::numNodesInQuadrilateral)
    {
        throw AlgorithmError("The initial face is not a quadrilateral");
    }

    // 2. Check if the point is inside the face
    std::vector<Point> polygonPoints;
    for (UInt n = 0; n < geometric::numNodesInQuadrilateral; ++n)
    {
        const auto node = m_mesh.m_facesNodes[initialFaceIndex][n];
        polygonPoints.emplace_back(m_mesh.Node(node));
    }
    const auto node = m_mesh.m_facesNodes[initialFaceIndex][0];
    polygonPoints.emplace_back(m_mesh.Node(node));

    Polygon polygon(polygonPoints, m_mesh.m_projection);
    if (!polygon.Contains(point))
    {
        throw AlgorithmError("The initial face does not contain the starting point");
    }

    // 3. Build the local coordinate system
    const auto numNodes = m_mesh.GetNumNodes();
    m_i = std::vector(numNodes, missing::intValue);
    m_j = std::vector(numNodes, missing::intValue);

    const auto firstEdge = m_mesh.m_facesEdges[initialFaceIndex][0];
    const auto secondEdge = m_mesh.m_facesEdges[initialFaceIndex][1];
    const auto thirdEdge = m_mesh.m_facesEdges[initialFaceIndex][2];
    const auto fourthEdge = m_mesh.m_facesEdges[initialFaceIndex][3];

    m_mapping.resize(-3, -3, 3, 3);

    const auto firstNodeIndex = m_mesh.FindCommonNode(firstEdge, secondEdge);
    m_j[firstNodeIndex] = 0;
    m_i[firstNodeIndex] = 0;

    m_mapping.setValue(0, 0, firstNodeIndex);

    const auto secondNodeIndex = m_mesh.FindCommonNode(secondEdge, thirdEdge);
    m_j[secondNodeIndex] = 0;
    m_i[secondNodeIndex] = 1;
    m_mapping.setValue(0, 1, secondNodeIndex);

    const auto thirdNodeIndex = m_mesh.FindCommonNode(thirdEdge, fourthEdge);
    m_i[thirdNodeIndex] = 1;
    m_j[thirdNodeIndex] = 1;
    m_mapping.setValue(1, 1, thirdNodeIndex);

    const auto fourthNodeIndex = m_mesh.FindCommonNode(fourthEdge, firstEdge);
    m_j[fourthNodeIndex] = 1;
    m_i[fourthNodeIndex] = 0;

    m_mapping.setValue(1, 0, fourthNodeIndex);

    // 4. Grow the front using the breath first search algorithm
    const auto numFaces = m_mesh.GetNumFaces();
    std::vector visitedFace(numFaces, false);
    std::queue<UInt> q;
    q.push(initialFaceIndex);
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
        for (auto d = 0u; d < geometric::numNodesInQuadrilateral; ++d)
        {
            const auto newFaceIndex = ComputeNeighbouringFaceNodes(face, localNodeMapping, d, visitedFace);
            if (newFaceIndex != missing::uintValue)
            {
                q.push(newFaceIndex);
            }
        }
    }
    return std::make_unique<CurvilinearGrid>(ComputeCurvilinearMatrix(), m_mesh.m_projection);
}

Eigen::Matrix<UInt, 2, 2> Mesh2DToCurvilinear::ComputeLocalNodeMapping(UInt face) const
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

UInt Mesh2DToCurvilinear::ComputeNeighbouringFaceNodes(const UInt face,
                                                       const Eigen::Matrix<UInt, 2, 2>& localNodeMapping,
                                                       const UInt d,
                                                       const std::vector<bool>& visitedFace)
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
    for (UInt e = 0u; e < geometric::numNodesInQuadrilateral; ++e)
    {
        if (m_mesh.m_facesEdges[newFace][e] == edgeIndex)
        {
            edgeIndexInNewFace = static_cast<int>(e);
            break;
        }
    }
    auto nextEdgeIndexInNewFace = edgeIndexInNewFace + 1;
    nextEdgeIndexInNewFace = nextEdgeIndexInNewFace == geometric::numNodesInQuadrilateral ? 0 : nextEdgeIndexInNewFace;
    const auto nextEdgeInNewFace = m_mesh.m_facesEdges[newFace][nextEdgeIndexInNewFace];
    const auto firstCommonNode = m_mesh.FindCommonNode(edgeIndex, nextEdgeInNewFace);
    const auto firstOtherNode = OtherNodeOfEdge(m_mesh.GetEdge(nextEdgeInNewFace), firstCommonNode);
    const auto iFirstOtherNode = m_i[firstCommonNode] + m_directionsDeltas[d][0];
    const auto jFirstOtherNode = m_j[firstCommonNode] + m_directionsDeltas[d][1];

    auto previousEdgeIndexInNewFace = edgeIndexInNewFace - 1;
    previousEdgeIndexInNewFace = previousEdgeIndexInNewFace == -1 ? geometric::numNodesInQuadrilateral - 1 : previousEdgeIndexInNewFace;
    const auto previousEdgeInNewFace = m_mesh.m_facesEdges[newFace][previousEdgeIndexInNewFace];
    const auto secondCommonNode = m_mesh.FindCommonNode(edgeIndex, previousEdgeInNewFace);
    const auto secondOtherNode = OtherNodeOfEdge(m_mesh.GetEdge(previousEdgeInNewFace), secondCommonNode);
    const auto iSecondCommonNode = m_i[secondCommonNode] + m_directionsDeltas[d][0];
    const auto jSecondCommonNode = m_j[secondCommonNode] + m_directionsDeltas[d][1];

    const auto invalid = (m_i[firstOtherNode] != missing::intValue && m_i[firstOtherNode] != iFirstOtherNode) ||
                         (m_j[firstOtherNode] != missing::intValue && m_j[firstOtherNode] != jFirstOtherNode) ||
                         (m_i[secondOtherNode] != missing::intValue && m_i[secondOtherNode] != iSecondCommonNode) ||
                         (m_j[secondOtherNode] != missing::intValue && m_j[secondOtherNode] != jSecondCommonNode);

    if (invalid)
    {
        return missing::uintValue;
    }
    if (!IsConnectionValid(firstOtherNode, iFirstOtherNode, jFirstOtherNode))
    {
        return missing::uintValue;
    }
    if (!IsConnectionValid(secondOtherNode, iSecondCommonNode, jSecondCommonNode))
    {
        return missing::uintValue;
    }

    m_i[firstOtherNode] = iFirstOtherNode;
    m_j[firstOtherNode] = jFirstOtherNode;
    m_i[secondOtherNode] = iSecondCommonNode;
    m_j[secondOtherNode] = jSecondCommonNode;
    m_mapping.setValue(jFirstOtherNode, iFirstOtherNode, firstOtherNode);
    m_mapping.setValue(jSecondCommonNode, iSecondCommonNode, secondOtherNode);
    return newFace;
}

bool Mesh2DToCurvilinear::IsConnectionValid(const UInt candidateNode, const int iCandidate, const int jCandidate)
{
    const int iLeft = iCandidate - 1;
    const int iRight = iCandidate + 1;
    const int jBottom = jCandidate - 1;
    const int jUp = jCandidate + 1;
    m_mapping.resize(jBottom, iLeft, jUp, iRight);

    if (m_mapping.IsValid(jCandidate, iLeft) && !CheckGridLine(m_mapping.getValue(jCandidate, iLeft), candidateNode))
    {
        return false;
    }

    if (m_mapping.IsValid(jCandidate, iRight) && !CheckGridLine(m_mapping.getValue(jCandidate, iRight), candidateNode))
    {
        return false;
    }

    if (m_mapping.IsValid(jBottom, iCandidate) && !CheckGridLine(m_mapping.getValue(jBottom, iCandidate), candidateNode))
    {
        return false;
    }

    if (m_mapping.IsValid(jUp, iCandidate) && !CheckGridLine(m_mapping.getValue(jUp, iCandidate), candidateNode))
    {
        return false;
    }

    return true;
}

bool Mesh2DToCurvilinear::CheckGridLine(const UInt validNode, const UInt candidateNode) const
{
    bool valid = false;
    for (auto e = 0u; e < m_mesh.m_nodesEdges[candidateNode].size(); ++e)
    {
        const auto edgeIndex = m_mesh.m_nodesEdges[candidateNode][e];
        const auto otherNode = OtherNodeOfEdge(m_mesh.GetEdge(edgeIndex), candidateNode);

        bool doCheck = m_mesh.GetNumFaceEdges(m_mesh.m_edgesFaces[edgeIndex][0]) == geometric::numNodesInQuadrilateral;

        if (m_mesh.m_edgesNumFaces[edgeIndex] == 2)
        {
            doCheck = doCheck || m_mesh.GetNumFaceEdges(m_mesh.m_edgesFaces[edgeIndex][1]) == geometric::numNodesInQuadrilateral;
        }

        if (doCheck && otherNode == validNode)
        {
            valid = true;
            break;
        }
    }

    return valid;
}

lin_alg::Matrix<Point> Mesh2DToCurvilinear::ComputeCurvilinearMatrix()
{
    auto validI = m_i | std::views::filter([](const auto& x)
                                           { return x != missing::intValue; });
    const auto [minI, maxI] = std::ranges::minmax_element(validI);

    auto validJ = m_j | std::views::filter([](const auto& x)
                                           { return x != missing::intValue; });
    const auto [minJ, maxJ] = std::ranges::minmax_element(validJ);

    const auto numM = *maxI - *minI + 1;
    const auto numN = *maxJ - *minJ + 1;

    lin_alg::Matrix<Point> result(numN, numM);

    for (auto n = 0u; n < m_mesh.GetNumNodes(); ++n)
    {
        if (m_j[n] != missing::intValue && m_i[n] != missing::intValue)
        {
            const int j = m_j[n] - *minJ;
            const int i = m_i[n] - *minI;
            result(j, i) = m_mesh.Node(n);
        }
    }

    return result;
}
