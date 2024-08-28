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

#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/FlipEdges.hpp>
#include <MeshKernel/LandBoundaries.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/UndoActions/CompoundUndoAction.hpp>

using meshkernel::FlipEdges;
using meshkernel::Mesh2D;

FlipEdges::FlipEdges(Mesh2D& mesh,
                     LandBoundaries& landBoundary,
                     bool triangulateFaces,
                     bool projectToLandBoundary) : m_mesh(mesh),
                                                   m_landBoundaries(landBoundary),
                                                   m_triangulateFaces(triangulateFaces),
                                                   m_projectToLandBoundary(projectToLandBoundary)
{
    if (m_projectToLandBoundary)
    {
        m_landBoundaries.FindNearestMeshBoundary(LandBoundaries::ProjectToLandBoundaryOption::WholeMesh);
    }
}

std::unique_ptr<meshkernel::UndoAction> FlipEdges::Compute() const
{
    std::unique_ptr<CompoundUndoAction> action = CompoundUndoAction::Create();

    m_mesh.Administrate(action.get());

    if (m_triangulateFaces)
    {
        action->Add(m_mesh.TriangulateFaces());
        m_mesh.Administrate(action.get());
    }

    const UInt MaxIter = 10;
    const auto numEdges = m_mesh.GetNumEdges();
    UInt numFlippedEdges = constants::missing::uintValue;

    for (UInt iteration = 0; iteration < MaxIter; ++iteration)
    {
        if (numFlippedEdges == 0)
        {
            break;
        }
        numFlippedEdges = 0;

        for (UInt e = 0; e < static_cast<UInt>(numEdges); e++)
        {

            if (m_mesh.IsEdgeOnBoundary(e))
            {
                continue;
            }

            // Triangles only
            auto const leftFace = m_mesh.m_edgesFaces[e][0];
            auto const rightFace = m_mesh.m_edgesFaces[e][1];

            if (leftFace == constants::missing::uintValue || rightFace == constants::missing::uintValue)
            {
                continue;
            }

            const auto NumEdgesLeftFace = m_mesh.GetNumFaceEdges(leftFace);
            const auto NumEdgesRightFace = m_mesh.GetNumFaceEdges(rightFace);
            if (NumEdgesLeftFace != constants::geometric::numNodesInTriangle ||
                NumEdgesRightFace != constants::geometric::numNodesInTriangle)
            {
                return action;
            }

            UInt nodeLeft = constants::missing::uintValue;
            UInt nodeRight = constants::missing::uintValue;
            const auto topologyFunctional = ComputeTopologyFunctional(e, nodeLeft, nodeRight);

            if (topologyFunctional >= 0)
            {
                continue;
            }

            // Check if nodes have been masked
            auto const firstNode = m_mesh.GetEdge(e).first;
            auto const secondNode = m_mesh.GetEdge(e).second;

            // Check if the quadrilateral composed by the two adjacent triangles is concave,
            // in which case the diagonals crosses
            const auto [areEdgesCrossing,
                        intersection,
                        crossProduct,
                        firstRatio,
                        secondRatio] = AreSegmentsCrossing(m_mesh.Node(firstNode),
                                                           m_mesh.Node(secondNode),
                                                           m_mesh.Node(nodeLeft),
                                                           m_mesh.Node(nodeRight),
                                                           false,
                                                           m_mesh.m_projection);

            if (!areEdgesCrossing)
            {
                continue;
            }

            // Flip the edges
            action->Add(m_mesh.ResetEdge(e, {nodeLeft, nodeRight}));
            numFlippedEdges++;

            // Find the other edges
            UInt firstEdgeLeftFace = constants::missing::uintValue;
            UInt firstEdgeRightFace = constants::missing::uintValue;
            UInt secondEdgeLeftFace = constants::missing::uintValue;
            UInt secondEdgeRightFace = constants::missing::uintValue;
            for (UInt i = 0; i < NumEdgesLeftFace; i++)
            {
                const auto edgeIndex = m_mesh.m_facesEdges[leftFace][i];
                if (edgeIndex == e)
                {
                    continue;
                }
                const auto first = m_mesh.GetEdge(edgeIndex).first;
                const auto second = m_mesh.GetEdge(edgeIndex).second;

                if (first == firstNode || second == firstNode)
                {
                    firstEdgeLeftFace = edgeIndex;
                }
                if (first == secondNode || second == secondNode)
                {
                    secondEdgeLeftFace = edgeIndex;
                }
            }

            for (UInt i = 0; i < NumEdgesRightFace; i++)
            {
                const auto edgeIndex = m_mesh.m_facesEdges[rightFace][i];
                if (edgeIndex == e)
                {
                    continue;
                }
                const auto first = m_mesh.GetEdge(edgeIndex).first;
                const auto second = m_mesh.GetEdge(edgeIndex).second;

                if (first == firstNode || second == firstNode)
                {
                    firstEdgeRightFace = edgeIndex;
                }
                if (first == secondNode || second == secondNode)
                {
                    secondEdgeRightFace = edgeIndex;
                }
            }

            // Change face orientation
            m_mesh.m_facesNodes[leftFace][0] = nodeLeft;
            m_mesh.m_facesNodes[leftFace][1] = nodeRight;
            m_mesh.m_facesNodes[leftFace][2] = firstNode;

            m_mesh.m_facesEdges[leftFace][0] = e;
            m_mesh.m_facesEdges[leftFace][1] = firstEdgeRightFace;
            m_mesh.m_facesEdges[leftFace][2] = firstEdgeLeftFace;

            m_mesh.m_facesNodes[rightFace][0] = nodeLeft;
            m_mesh.m_facesNodes[rightFace][1] = nodeRight;
            m_mesh.m_facesNodes[rightFace][2] = secondNode;

            m_mesh.m_facesEdges[rightFace][0] = e;
            m_mesh.m_facesEdges[rightFace][1] = secondEdgeRightFace;
            m_mesh.m_facesEdges[rightFace][2] = secondEdgeLeftFace;

            if (m_mesh.m_edgesFaces[firstEdgeRightFace][0] == rightFace)
            {
                m_mesh.m_edgesFaces[firstEdgeRightFace][0] = leftFace;
            }
            else
            {
                m_mesh.m_edgesFaces[firstEdgeRightFace][1] = leftFace;
            }

            if (m_mesh.m_edgesFaces[secondEdgeLeftFace][0] == leftFace)
            {
                m_mesh.m_edgesFaces[secondEdgeLeftFace][0] = rightFace;
            }
            else
            {
                m_mesh.m_edgesFaces[secondEdgeLeftFace][1] = rightFace;
            }

            m_mesh.m_nodesNumEdges[firstNode] = m_mesh.m_nodesNumEdges[firstNode] - 1;
            m_mesh.m_nodesNumEdges[secondNode] = m_mesh.m_nodesNumEdges[secondNode] - 1;
            m_mesh.m_nodesNumEdges[nodeLeft] = m_mesh.m_nodesNumEdges[nodeLeft] + 1;
            m_mesh.m_nodesNumEdges[nodeRight] = m_mesh.m_nodesNumEdges[nodeRight] + 1;

            // Delete edge from m_mesh.m_nodesEdges[firstNode]
            DeleteEdgeFromNode(e, firstNode);

            // Delete edge from m_mesh.m_nodesEdges[secondNode]
            DeleteEdgeFromNode(e, secondNode);

            // Add edge to m_mesh.m_nodesEdges[kl]
            m_mesh.m_nodesEdges[nodeLeft].resize(m_mesh.m_nodesNumEdges[nodeLeft]);
            m_mesh.m_nodesEdges[nodeLeft].back() = e;
            m_mesh.SortEdgesInCounterClockWiseOrder(nodeLeft, nodeLeft);

            // Add edge to m_mesh.m_nodesEdges[kr]
            m_mesh.m_nodesEdges[nodeRight].resize(m_mesh.m_nodesNumEdges[nodeRight]);
            m_mesh.m_nodesEdges[nodeRight].back() = e;
            m_mesh.SortEdgesInCounterClockWiseOrder(nodeRight, nodeRight);
        }
    }

    if (numFlippedEdges != 0)
    {
        throw AlgorithmError("FlipEdges::Compute: Could not complete, there are still edges left to be flipped.");
    }

    // Perform mesh administration
    m_mesh.Administrate(action.get());
    return action;
}

void FlipEdges::DeleteEdgeFromNode(UInt edge, UInt firstNode) const
{
    // Update node, delete edge from m_mesh.m_nodesEdges[firstNode]
    UInt kk = 0;
    while (m_mesh.m_nodesEdges[firstNode][kk] != edge && kk < m_mesh.m_nodesNumEdges[firstNode])
    {
        kk = kk + 1;
    }
    if (m_mesh.m_nodesEdges[firstNode][kk] != edge)
    {
        throw std::invalid_argument("FlipEdges::DeleteEdgeFromNode: The edge does not match the given node.");
    }

    UInt count = 0;
    for (UInt i = 0; i < m_mesh.m_nodesNumEdges[firstNode] + 1; i++)
    {
        if (i + 1 <= kk || i > kk)
        {
            m_mesh.m_nodesEdges[firstNode][count] = m_mesh.m_nodesEdges[firstNode][i];
            count++;
        }
    }

    m_mesh.m_nodesEdges[firstNode].resize(m_mesh.m_nodesNumEdges[firstNode]);
}

int FlipEdges::ComputeTopologyFunctional(UInt edge,
                                         UInt& nodeLeft,
                                         UInt& nodeRight) const
{
    const int largeTopologyFunctionalValue = 1000;

    if (m_mesh.IsEdgeOnBoundary(edge))
    {
        return largeTopologyFunctionalValue;
    }

    const auto firstNode = m_mesh.GetEdge(edge).first;
    const auto secondNode = m_mesh.GetEdge(edge).second;
    const auto faceL = m_mesh.m_edgesFaces[edge][0];
    const auto faceR = m_mesh.m_edgesFaces[edge][1];
    const auto NumEdgesLeftFace = m_mesh.GetNumFaceEdges(faceL);
    const auto NumEdgesRightFace = m_mesh.GetNumFaceEdges(faceR);

    if (NumEdgesLeftFace != constants::geometric::numNodesInTriangle ||
        NumEdgesRightFace != constants::geometric::numNodesInTriangle)
    {
        return largeTopologyFunctionalValue;
    }

    // Find the nodes that are connected to both k1 and k
    UInt sumIndicesLeftFace = 0;
    UInt sumIndicesRightFace = 0;
    for (auto i = 0; i < 3; i++)
    {
        sumIndicesLeftFace += m_mesh.m_facesNodes[faceL][i];
        sumIndicesRightFace += m_mesh.m_facesNodes[faceR][i];
    }

    nodeLeft = sumIndicesLeftFace - firstNode - secondNode;
    nodeRight = sumIndicesRightFace - firstNode - secondNode;

    if (nodeLeft == constants::missing::uintValue || nodeRight == constants::missing::uintValue)
    {
        return largeTopologyFunctionalValue;
    }

    // Check that kl is part of faceL
    bool nodeFound = false;
    for (UInt i = 0; i < NumEdgesLeftFace; i++)
    {
        if (m_mesh.m_facesNodes[faceL][i] == nodeLeft)
        {
            nodeFound = true;
            break;
        }
    }

    if (!nodeFound)
    {
        return largeTopologyFunctionalValue;
    }

    // Check that kr is part of faceR
    nodeFound = false;
    for (UInt i = 0; i < NumEdgesRightFace; i++)
    {
        if (m_mesh.m_facesNodes[faceR][i] == nodeRight)
        {
            nodeFound = true;
            break;
        }
    }

    if (!nodeFound)
    {
        return largeTopologyFunctionalValue;
    }

    //  Compute the change in functional
    const auto n1 = static_cast<int>(m_mesh.m_nodesNumEdges[firstNode]) - static_cast<int>(OptimalNumberOfConnectedNodes(firstNode));
    const auto n2 = static_cast<int>(m_mesh.m_nodesNumEdges[secondNode]) - static_cast<int>(OptimalNumberOfConnectedNodes(secondNode));
    auto nL = static_cast<int>(m_mesh.m_nodesNumEdges[nodeLeft]) - static_cast<int>(OptimalNumberOfConnectedNodes(nodeLeft));
    auto nR = static_cast<int>(m_mesh.m_nodesNumEdges[nodeRight]) - static_cast<int>(OptimalNumberOfConnectedNodes(nodeRight));

    if (m_projectToLandBoundary && m_landBoundaries.GetNumNodes() > 0)
    {
        if (m_landBoundaries.m_meshNodesLandBoundarySegments[firstNode] != constants::missing::uintValue && m_landBoundaries.m_meshNodesLandBoundarySegments[secondNode] != constants::missing::uintValue)
        {
            // Edge is associated with a land boundary -> keep the edge
            return largeTopologyFunctionalValue;
        }

        const auto n1L = DifferenceFromOptimum(firstNode, secondNode, nodeLeft);
        const auto n1R = DifferenceFromOptimum(firstNode, secondNode, nodeRight);

        const auto n2R = DifferenceFromOptimum(secondNode, firstNode, nodeLeft);
        const auto n2L = DifferenceFromOptimum(secondNode, firstNode, nodeRight);

        nL = DifferenceFromOptimum(nodeLeft, firstNode, secondNode);
        nR = DifferenceFromOptimum(nodeRight, firstNode, secondNode);

        const auto topologyFunctional = (n1L - 1) * (n1L - 1) +
                                        (n1R - 1) * (n1R - 1) +
                                        (n2L - 1) * (n2L - 1) +
                                        (n2R - 1) * (n2R - 1) +
                                        2 * ((nL + 1) * (nL + 1) + (nR + 1) * (nR + 1)) -
                                        (n1L * n1L + n1R * n1R + n2L * n2L + n2R * n2R + 2 * (nL * nL + nR * nR));
        return topologyFunctional;
    }

    const auto topologyFunctional = (n1 - 1) * (n1 - 1) +
                                    (n2 - 1) * (n2 - 1) +
                                    (nL + 1) * (nL + 1) +
                                    (nR + 1) * (nR + 1) -
                                    (n1 * n1 + n2 * n2 + nL * nL + nR * nR);

    return topologyFunctional;
}

int FlipEdges::DifferenceFromOptimum(UInt nodeIndex, UInt firstNode, UInt secondNode) const
{
    if (m_landBoundaries.m_meshNodesLandBoundarySegments[nodeIndex] == constants::missing::uintValue)
    {
        return static_cast<int>(m_mesh.m_nodesNumEdges[nodeIndex]) - static_cast<int>(OptimalNumberOfConnectedNodes(nodeIndex));
    }

    // Connected edges needs to be counterclockwise
    const auto sign = sgn(crossProduct(m_mesh.Node(nodeIndex), m_mesh.Node(firstNode), m_mesh.Node(firstNode), m_mesh.Node(secondNode), m_mesh.m_projection));
    const auto isClockWise = sign < 0 ? true : false;
    if (isClockWise)
    {
        const auto firstNodeTemp = firstNode;
        firstNode = secondNode;
        secondNode = firstNodeTemp;
    }

    // Find the first edge connecting firstNode
    UInt edgeIndexConnectingFirstNode = constants::missing::uintValue;
    for (UInt i = 0; i < m_mesh.m_nodesNumEdges[nodeIndex]; i++)
    {
        const auto edgeIndex = m_mesh.m_nodesEdges[nodeIndex][i];

        if (m_mesh.GetEdge(edgeIndex).first == firstNode || m_mesh.GetEdge(edgeIndex).second == firstNode)
        {
            edgeIndexConnectingFirstNode = i;
            break;
        }
    }

    if (edgeIndexConnectingFirstNode == constants::missing::uintValue)
    {
        return 0;
    }

    // Find the first edge connecting secondNode
    UInt edgeIndexConnectingSecondNode = constants::missing::uintValue;
    for (UInt i = 0; i < m_mesh.m_nodesNumEdges[nodeIndex]; i++)
    {
        const auto edgeIndex = m_mesh.m_nodesEdges[nodeIndex][i];

        if (m_mesh.GetEdge(edgeIndex).first == secondNode || m_mesh.GetEdge(edgeIndex).second == secondNode)
        {
            edgeIndexConnectingSecondNode = i;
            break;
        }
    }

    if (edgeIndexConnectingSecondNode == constants::missing::uintValue)
    {
        return 0;
    }

    // Count the numbers of edges clockwise from the one connecting indexFirstNode
    // that are not in a land or mesh boundary path
    auto currentEdgeIndexInNodeEdges = edgeIndexConnectingFirstNode;
    auto edgeIndex = m_mesh.m_nodesEdges[nodeIndex][currentEdgeIndexInNodeEdges];
    auto otherNode = OtherNodeOfEdge(m_mesh.GetEdge(edgeIndex), nodeIndex);

    UInt num = 1;
    while (m_landBoundaries.m_meshNodesLandBoundarySegments[otherNode] == constants::missing::uintValue &&
           !m_mesh.IsEdgeOnBoundary(edgeIndex) &&
           currentEdgeIndexInNodeEdges != edgeIndexConnectingSecondNode)
    {
        currentEdgeIndexInNodeEdges = NextCircularBackwardIndex(currentEdgeIndexInNodeEdges, m_mesh.m_nodesNumEdges[nodeIndex]);
        edgeIndex = m_mesh.m_nodesEdges[nodeIndex][currentEdgeIndexInNodeEdges];
        otherNode = OtherNodeOfEdge(m_mesh.GetEdge(edgeIndex), nodeIndex);
        num++;
    }

    UInt firstEdgeInPathIndex = constants::missing::uintValue;
    if (m_landBoundaries.m_meshNodesLandBoundarySegments[otherNode] != constants::missing::uintValue ||
        m_mesh.IsEdgeOnBoundary(edgeIndex))
    {
        firstEdgeInPathIndex = edgeIndex;
    }

    // If not all edges are visited, count counterclockwise from the one connecting indexSecondNode
    UInt secondEdgeInPathIndex = constants::missing::uintValue;
    if (currentEdgeIndexInNodeEdges != edgeIndexConnectingSecondNode)
    {
        currentEdgeIndexInNodeEdges = edgeIndexConnectingSecondNode;
        edgeIndex = m_mesh.m_nodesEdges[nodeIndex][currentEdgeIndexInNodeEdges];
        otherNode = OtherNodeOfEdge(m_mesh.GetEdge(edgeIndex), nodeIndex);
        num = num + 1;
        while (m_landBoundaries.m_meshNodesLandBoundarySegments[otherNode] == constants::missing::uintValue &&
               !m_mesh.IsEdgeOnBoundary(edgeIndex) &&
               currentEdgeIndexInNodeEdges != edgeIndexConnectingFirstNode &&
               edgeIndex != firstEdgeInPathIndex)
        {
            currentEdgeIndexInNodeEdges = NextCircularForwardIndex(currentEdgeIndexInNodeEdges, m_mesh.m_nodesNumEdges[nodeIndex]);
            edgeIndex = m_mesh.m_nodesEdges[nodeIndex][currentEdgeIndexInNodeEdges];
            otherNode = OtherNodeOfEdge(m_mesh.GetEdge(edgeIndex), nodeIndex);

            if (currentEdgeIndexInNodeEdges != edgeIndexConnectingFirstNode && edgeIndex != firstEdgeInPathIndex)
            {
                num++;
            }
        }

        if ((m_landBoundaries.m_meshNodesLandBoundarySegments[otherNode] != constants::missing::uintValue ||
             m_mesh.IsEdgeOnBoundary(edgeIndex)) &&
            edgeIndex != firstEdgeInPathIndex)
        {
            secondEdgeInPathIndex = edgeIndex;
        }
    }

    // The number of nodes is larger than the connected ones, should not happen
    if (num > m_mesh.m_nodesNumEdges[nodeIndex])
    {
        return 0;
    }

    if (firstEdgeInPathIndex != constants::missing::uintValue && secondEdgeInPathIndex != constants::missing::uintValue)
    {
        // Internal boundary
        return 4;
    }

    return 6;
}

meshkernel::UInt FlipEdges::OptimalNumberOfConnectedNodes(UInt index) const
{
    UInt optimalNumber = 6;
    if (m_mesh.m_nodesTypes[index] == 2)
    {
        optimalNumber = 4;
    }
    if (m_mesh.m_nodesTypes[index] == 3)
    {
        optimalNumber = 3;
    }

    return optimalNumber;
}
