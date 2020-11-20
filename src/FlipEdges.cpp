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
#include <MeshKernel/Operations.cpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh.hpp>
#include <MeshKernel/FlipEdges.hpp>
#include <MeshKernel/LandBoundaries.hpp>
#include <MeshKernel/Exceptions.hpp>

meshkernel::FlipEdges::FlipEdges(std::shared_ptr<Mesh> mesh,
                                 std::shared_ptr<LandBoundaries> landBoundary,
                                 bool triangulateFaces,
                                 bool projectToLandBoundary) : m_mesh(mesh),
                                                               m_landBoundaries(landBoundary),
                                                               m_triangulateFaces(triangulateFaces),
                                                               m_projectToLandBoundary(projectToLandBoundary)
{
    if (m_landBoundaries->GetNumNodes() <= 0)
    {
        m_projectToLandBoundary = false;
    }
    if (m_projectToLandBoundary)
    {
        try
        {
            m_landBoundaries->FindNearestMeshBoundary(4);
        }
        catch (const std::exception&)
        {
            // TODO: log exception: need to rethrow the exception
            m_projectToLandBoundary = false;
        }
    }
};

void meshkernel::FlipEdges::Compute() const
{

    m_mesh->Administrate(Mesh::AdministrationOptions::AdministrateMeshEdgesAndFaces);

    if (m_triangulateFaces)
    {
        m_mesh->TriangulateFaces();
        m_mesh->Administrate(Mesh::AdministrationOptions::AdministrateMeshEdgesAndFaces);
    }

    const int MaxIter = 10;
    const int numEdges = m_mesh->GetNumEdges();
    int numFlippedEdges = intMissingValue;

    for (int iter = 0; iter < MaxIter; iter++)
    {
        if (numFlippedEdges == 0)
        {
            break;
        }
        else
        {
            numFlippedEdges = 0;
        }

        for (int e = 0; e < numEdges; e++)
        {
            // check if nodes have been masked
            auto const firstNode = m_mesh->m_edges[e].first;
            auto const secondNode = m_mesh->m_edges[e].second;

            if (m_mesh->IsEdgeOnBoundary(e))
            {
                continue;
            }

            // triangles only
            auto const leftFace = m_mesh->m_edgesFaces[e][0];
            auto const rightFace = m_mesh->m_edgesFaces[e][1];

            const auto NumEdgesLeftFace = m_mesh->GetNumFaceEdges(leftFace);
            const auto NumEdgesRightFace = m_mesh->GetNumFaceEdges(rightFace);
            if (NumEdgesLeftFace != 3 || NumEdgesRightFace != 3)
            {
                return;
            }

            int nodeLeft = -1;
            int nodeRight = -1;
            int topologyFunctional = 1000;
            ComputeTopologyFunctional(e, nodeLeft, nodeRight, topologyFunctional);

            if (topologyFunctional < 0)
            {
                // Check if the quadrilateral composed by the two adjacent triangles is concave,
                // in which case the diagonals crosses
                Point intersection;
                double crossProduct;
                double firstRatio;
                double secondRatio;
                const auto areEdgesCrossing = AreLinesCrossing(m_mesh->m_nodes[firstNode],
                                                               m_mesh->m_nodes[secondNode],
                                                               m_mesh->m_nodes[nodeLeft],
                                                               m_mesh->m_nodes[nodeRight],
                                                               false,
                                                               intersection,
                                                               crossProduct,
                                                               firstRatio,
                                                               secondRatio,
                                                               m_mesh->m_projection);

                if (!areEdgesCrossing)
                {
                    continue;
                }

                // Flip the edges
                m_mesh->m_edges[e].first = nodeLeft;
                m_mesh->m_edges[e].second = nodeRight;
                numFlippedEdges++;

                // Find the other edges
                int firstEdgeLeftFace;
                int firstEdgeRightFace;
                int secondEdgeLeftFace;
                int secondEdgeRightFace;
                for (int i = 0; i < NumEdgesLeftFace; i++)
                {
                    int edgeIndex = m_mesh->m_facesEdges[leftFace][i];
                    if (edgeIndex == e)
                    {
                        continue;
                    }
                    const int first = m_mesh->m_edges[edgeIndex].first;
                    const int second = m_mesh->m_edges[edgeIndex].second;

                    if (first == firstNode || second == firstNode)
                    {
                        firstEdgeLeftFace = edgeIndex;
                    }
                    if (first == secondNode || second == secondNode)
                    {
                        secondEdgeLeftFace = edgeIndex;
                    }
                }

                for (int i = 0; i < NumEdgesRightFace; i++)
                {
                    int edgeIndex = m_mesh->m_facesEdges[rightFace][i];
                    if (edgeIndex == e)
                    {
                        continue;
                    }
                    const int first = m_mesh->m_edges[edgeIndex].first;
                    const int second = m_mesh->m_edges[edgeIndex].second;

                    if (first == firstNode || second == firstNode)
                    {
                        firstEdgeRightFace = edgeIndex;
                    }
                    if (first == secondNode || second == secondNode)
                    {
                        secondEdgeRightFace = edgeIndex;
                    }
                }

                // change face orientation
                m_mesh->m_facesNodes[leftFace][0] = nodeLeft;
                m_mesh->m_facesNodes[leftFace][1] = nodeRight;
                m_mesh->m_facesNodes[leftFace][2] = firstNode;

                m_mesh->m_facesEdges[leftFace][0] = e;
                m_mesh->m_facesEdges[leftFace][1] = firstEdgeRightFace;
                m_mesh->m_facesEdges[leftFace][2] = firstEdgeLeftFace;

                m_mesh->m_facesNodes[rightFace][0] = nodeLeft;
                m_mesh->m_facesNodes[rightFace][1] = nodeRight;
                m_mesh->m_facesNodes[rightFace][2] = secondNode;

                m_mesh->m_facesEdges[rightFace][0] = e;
                m_mesh->m_facesEdges[rightFace][1] = secondEdgeRightFace;
                m_mesh->m_facesEdges[rightFace][2] = secondEdgeLeftFace;

                if (m_mesh->m_edgesFaces[firstEdgeRightFace][0] == rightFace)
                {
                    m_mesh->m_edgesFaces[firstEdgeRightFace][0] = leftFace;
                }
                else
                {
                    m_mesh->m_edgesFaces[firstEdgeRightFace][1] = leftFace;
                }

                if (m_mesh->m_edgesFaces[secondEdgeLeftFace][0] == leftFace)
                {
                    m_mesh->m_edgesFaces[secondEdgeLeftFace][0] = rightFace;
                }
                else
                {
                    m_mesh->m_edgesFaces[secondEdgeLeftFace][1] = rightFace;
                }

                m_mesh->m_nodesNumEdges[firstNode] = m_mesh->m_nodesNumEdges[firstNode] - 1;
                m_mesh->m_nodesNumEdges[secondNode] = m_mesh->m_nodesNumEdges[secondNode] - 1;
                m_mesh->m_nodesNumEdges[nodeLeft] = m_mesh->m_nodesNumEdges[nodeLeft] + 1;
                m_mesh->m_nodesNumEdges[nodeRight] = m_mesh->m_nodesNumEdges[nodeRight] + 1;

                // Delete edge from m_mesh->m_nodesEdges[firstNode]
                DeleteEdgeFromNode(e, firstNode);

                // Delete edge from m_mesh->m_nodesEdges[secondNode]
                DeleteEdgeFromNode(e, secondNode);

                // Add edge to m_mesh->m_nodesEdges[kl]
                ResizeVectorIfNeeded(m_mesh->m_nodesNumEdges[nodeLeft], m_mesh->m_nodesEdges[nodeLeft]);
                for (int i = 0; i < m_mesh->m_nodesNumEdges[nodeLeft] - 1; i++)
                {
                    m_mesh->m_nodesEdges[nodeLeft][i] = m_mesh->m_nodesEdges[nodeLeft][i];
                }
                m_mesh->m_nodesEdges[nodeLeft][m_mesh->m_nodesNumEdges[nodeLeft] - 1] = e;
                m_mesh->SortEdgesInCounterClockWiseOrder(nodeLeft);

                // Add edge to m_mesh->m_nodesEdges[kr]
                ResizeVectorIfNeeded(m_mesh->m_nodesNumEdges[nodeRight], m_mesh->m_nodesEdges[nodeRight]);
                for (int i = 0; i < m_mesh->m_nodesNumEdges[nodeRight] - 1; i++)
                {
                    m_mesh->m_nodesEdges[nodeRight][i] = m_mesh->m_nodesEdges[nodeRight][i];
                }
                m_mesh->m_nodesEdges[nodeRight][m_mesh->m_nodesNumEdges[nodeRight] - 1] = e;
                m_mesh->SortEdgesInCounterClockWiseOrder(nodeRight);
            }
        }
    }

    if (numFlippedEdges != 0)
    {
        throw AlgorithmError("FlipEdges::Compute: Could not complete, there are still edges left to be flipped.");
    }

    // Perform mesh administration
    m_mesh->Administrate(Mesh::AdministrationOptions::AdministrateMeshEdgesAndFaces);
}

void meshkernel::FlipEdges::DeleteEdgeFromNode(int edge, int firstNode) const
{
    // Update nod, delete edge from m_mesh->m_nodesEdges[firstNode]
    int kk = 0;
    while (m_mesh->m_nodesEdges[firstNode][kk] != edge &&
           kk < m_mesh->m_nodesNumEdges[firstNode])
    {
        kk = kk + 1;
    }
    if (m_mesh->m_nodesEdges[firstNode][kk] != edge)
    {
        throw std::invalid_argument("FlipEdges::DeleteEdgeFromNode: The edge does not match the given node.");
    }

    int count = 0;
    for (int i = 0; i < m_mesh->m_nodesNumEdges[firstNode] + 1; i++)
    {
        if (i <= kk - 1 || i > kk)
        {
            m_mesh->m_nodesEdges[firstNode][count] = m_mesh->m_nodesEdges[firstNode][i];
            count++;
        }
    }
    ResizeVectorIfNeeded(m_mesh->m_nodesNumEdges[firstNode], m_mesh->m_nodesEdges[firstNode]);
}

void meshkernel::FlipEdges::ComputeTopologyFunctional(int edge,
                                                      int& nodeLeft,
                                                      int& nodeRight,
                                                      int& topologyFunctional) const
{

    if (m_mesh->IsEdgeOnBoundary(edge))
    {
        return;
    }

    const auto firstNode = m_mesh->m_edges[edge].first;
    const auto secondNode = m_mesh->m_edges[edge].second;
    const auto faceL = m_mesh->m_edgesFaces[edge][0];
    const auto faceR = m_mesh->m_edgesFaces[edge][1];
    const auto NumEdgesLeftFace = m_mesh->GetNumFaceEdges(faceL);
    const auto NumEdgesRightFace = m_mesh->GetNumFaceEdges(faceR);

    if (NumEdgesLeftFace != 3 || NumEdgesRightFace != 3)
    {
        return;
    }

    // find the nodes that are connected to both k1 and k
    int sumIndicesLeftFace = 0;
    int sumIndicesRightFace = 0;
    for (int i = 0; i < 3; i++)
    {
        sumIndicesLeftFace += m_mesh->m_facesNodes[faceL][i];
        sumIndicesRightFace += m_mesh->m_facesNodes[faceR][i];
    }

    nodeLeft = sumIndicesLeftFace - firstNode - secondNode;
    nodeRight = sumIndicesRightFace - firstNode - secondNode;

    if (nodeLeft < 0 || nodeRight < 0)
    {
        return;
    }

    // check that kl is part of faceL
    bool nodeFound = false;
    for (int i = 0; i < NumEdgesLeftFace; i++)
    {
        if (m_mesh->m_facesNodes[faceL][i] == nodeLeft)
        {
            nodeFound = true;
            break;
        }
    }

    if (!nodeFound)
    {
        return;
    }

    // check that kr is part of faceR
    nodeFound = false;
    for (int i = 0; i < NumEdgesRightFace; i++)
    {
        if (m_mesh->m_facesNodes[faceR][i] == nodeRight)
        {
            nodeFound = true;
            break;
        }
    }

    if (!nodeFound)
    {
        return;
    }

    //  compute the change in functional
    int n1 = m_mesh->m_nodesNumEdges[firstNode] - OptimalNumberOfConnectedNodes(firstNode);
    int n2 = m_mesh->m_nodesNumEdges[secondNode] - OptimalNumberOfConnectedNodes(secondNode);
    int nL = m_mesh->m_nodesNumEdges[nodeLeft] - OptimalNumberOfConnectedNodes(nodeLeft);
    int nR = m_mesh->m_nodesNumEdges[nodeRight] - OptimalNumberOfConnectedNodes(nodeRight);

    if (m_projectToLandBoundary)
    {
        if (m_landBoundaries->m_meshNodesLandBoundarySegments[firstNode] >= 0 && m_landBoundaries->m_meshNodesLandBoundarySegments[secondNode] >= 0)
        {
            //edge is associated with a land boundary, keep the edge
            topologyFunctional = 1000;
            return;
        }
        else
        {
            const auto n1L = DifferenceFromOptimum(firstNode, secondNode, nodeLeft);
            const auto n1R = DifferenceFromOptimum(firstNode, secondNode, nodeRight);

            const auto n2R = DifferenceFromOptimum(secondNode, firstNode, nodeLeft);
            const auto n2L = DifferenceFromOptimum(secondNode, firstNode, nodeRight);

            nL = DifferenceFromOptimum(nodeLeft, firstNode, secondNode);
            nR = DifferenceFromOptimum(nodeRight, firstNode, secondNode);

            topologyFunctional = (n1L - 1) * (n1L - 1) +
                                 (n1R - 1) * (n1R - 1) +
                                 (n2L - 1) * (n2L - 1) +
                                 (n2R - 1) * (n2R - 1) +
                                 2 * ((nL + 1) * (nL + 1) + (nR + 1) * (nR + 1)) -
                                 (n1L * n1L + n1R * n1R + n2L * n2L + n2R * n2R + 2 * (nL * nL + nR * nR));
        }
    }
    else
    {
        topologyFunctional = (n1 - 1) * (n1 - 1) +
                             (n2 - 1) * (n2 - 1) +
                             (nL + 1) * (nL + 1) +
                             (nR + 1) * (nR + 1) -
                             (n1 * n1 + n2 * n2 + nL * nL + nR * nR);
    }
}

//comp_nnow
int meshkernel::FlipEdges::DifferenceFromOptimum(int nodeIndex, int firstNode, int secondNode) const
{
    if (m_landBoundaries->m_meshNodesLandBoundarySegments[nodeIndex] < 0)
    {
        return m_mesh->m_nodesNumEdges[nodeIndex] - OptimalNumberOfConnectedNodes(nodeIndex);
    }

    // connected edges needs to be counterclockwise
    int sign = TwoSegmentsSign(m_mesh->m_nodes[nodeIndex], m_mesh->m_nodes[firstNode], m_mesh->m_nodes[firstNode], m_mesh->m_nodes[secondNode], m_mesh->m_projection);
    bool isClockWise = sign < 0 ? true : false;
    if (isClockWise)
    {
        int firstNodeTemp = firstNode;
        firstNode = secondNode;
        secondNode = firstNodeTemp;
    }

    // find the first edge connecting firstNode
    int edgeIndexConnectingFirstNode = -1;
    for (int i = 0; i < m_mesh->m_nodesNumEdges[nodeIndex]; i++)
    {
        const auto edgeIndex = m_mesh->m_nodesEdges[nodeIndex][i];

        if (m_mesh->m_edges[edgeIndex].first == firstNode || m_mesh->m_edges[edgeIndex].second == firstNode)
        {
            edgeIndexConnectingFirstNode = i;
            break;
        }
    }

    if (edgeIndexConnectingFirstNode == -1)
    {
        return 0;
    }

    // find the first edge connecting secondNode
    int edgeIndexConnectingSecondNode = -1;
    for (int i = 0; i < m_mesh->m_nodesNumEdges[nodeIndex]; i++)
    {
        const auto edgeIndex = m_mesh->m_nodesEdges[nodeIndex][i];

        if (m_mesh->m_edges[edgeIndex].first == secondNode || m_mesh->m_edges[edgeIndex].second == secondNode)
        {
            edgeIndexConnectingSecondNode = i;
            break;
        }
    }

    if (edgeIndexConnectingSecondNode == -1)
    {
        return 0;
    }

    // count the numbers of edges clockwise from the one connecting indexFirstNode
    // that are not in a land or mesh boundary path
    int currentEdgeIndexInNodeEdges = edgeIndexConnectingFirstNode;
    int edgeIndex = m_mesh->m_nodesEdges[nodeIndex][currentEdgeIndexInNodeEdges];
    int otherNode = m_mesh->m_edges[edgeIndex].first + m_mesh->m_edges[edgeIndex].second - nodeIndex;
    int num = 1;
    while (m_landBoundaries->m_meshNodesLandBoundarySegments[otherNode] < 0 &&
           !m_mesh->IsEdgeOnBoundary(edgeIndex) &&
           currentEdgeIndexInNodeEdges != edgeIndexConnectingSecondNode)
    {
        currentEdgeIndexInNodeEdges = NextCircularBackwardIndex(currentEdgeIndexInNodeEdges, m_mesh->m_nodesNumEdges[nodeIndex]);
        edgeIndex = m_mesh->m_nodesEdges[nodeIndex][currentEdgeIndexInNodeEdges];
        otherNode = m_mesh->m_edges[edgeIndex].first + m_mesh->m_edges[edgeIndex].second - nodeIndex;
        num++;
    }

    int firstEdgeInPathIndex = -1;
    if (m_landBoundaries->m_meshNodesLandBoundarySegments[otherNode] >= 0 ||
        m_mesh->IsEdgeOnBoundary(edgeIndex))
    {
        firstEdgeInPathIndex = edgeIndex;
    }

    // If not all edges are visited, count counterclockwise from the one connecting indexSecondNode
    int secondEdgeInPathIndex = -1;
    if (currentEdgeIndexInNodeEdges != edgeIndexConnectingSecondNode)
    {
        currentEdgeIndexInNodeEdges = edgeIndexConnectingSecondNode;
        edgeIndex = m_mesh->m_nodesEdges[nodeIndex][currentEdgeIndexInNodeEdges];
        otherNode = m_mesh->m_edges[edgeIndex].first + m_mesh->m_edges[edgeIndex].second - nodeIndex;
        num = num + 1;
        while (m_landBoundaries->m_meshNodesLandBoundarySegments[otherNode] < 0 &&
               !m_mesh->IsEdgeOnBoundary(edgeIndex) &&
               currentEdgeIndexInNodeEdges != edgeIndexConnectingFirstNode &&
               edgeIndex != firstEdgeInPathIndex)
        {
            currentEdgeIndexInNodeEdges = NextCircularForwardIndex(currentEdgeIndexInNodeEdges, m_mesh->m_nodesNumEdges[nodeIndex]);
            edgeIndex = m_mesh->m_nodesEdges[nodeIndex][currentEdgeIndexInNodeEdges];
            otherNode = m_mesh->m_edges[edgeIndex].first + m_mesh->m_edges[edgeIndex].second - nodeIndex;

            if (currentEdgeIndexInNodeEdges != edgeIndexConnectingFirstNode && edgeIndex != firstEdgeInPathIndex)
            {
                num++;
            }
        }

        if ((m_landBoundaries->m_meshNodesLandBoundarySegments[otherNode] >= 0 ||
             m_mesh->IsEdgeOnBoundary(edgeIndex)) &&
            edgeIndex != firstEdgeInPathIndex)
        {
            secondEdgeInPathIndex = edgeIndex;
        }
    }

    // the number of nodes is larger than the connected ones, should not happen
    if (num > m_mesh->m_nodesNumEdges[nodeIndex])
    {
        return 0;
    }

    int numopt = 6;
    if (firstEdgeInPathIndex >= 0 && secondEdgeInPathIndex >= 0)
    {
        // internal boundary
        numopt = 4;
    }

    return numopt;
};

int meshkernel::FlipEdges::OptimalNumberOfConnectedNodes(int index) const
{
    int optimalNumber = 6;
    if (m_mesh->m_nodesTypes[index] == 2)
    {
        optimalNumber = 4;
    }
    if (m_mesh->m_nodesTypes[index] == 3)
    {
        optimalNumber = 3;
    }

    return optimalNumber;
}
