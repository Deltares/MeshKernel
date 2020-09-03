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
#include <iostream>
#include <algorithm>
#include "Operations.cpp"
#include "Entities.hpp"
#include "Mesh.hpp"
#include "FlipEdges.hpp"
#include "LandBoundaries.hpp"

GridGeom::FlipEdges::FlipEdges() :
    m_mesh(nullptr),
    m_landBoundaries(nullptr)
{
}

GridGeom::FlipEdges::FlipEdges(Mesh* mesh, LandBoundaries* landBoundary, bool triangulateFaces, bool projectToLandBoundary) :
    m_mesh(mesh),
    m_landBoundaries(landBoundary),
    m_triangulateFaces(triangulateFaces),
    m_projectToLandBoundary(projectToLandBoundary)
{
    m_triangulateFaces = true;
    m_projectToLandBoundary = true;

    if (m_projectToLandBoundary)
    {
        bool successful = m_landBoundaries->FindNearestMeshBoundary(4);
        if (!successful)
        {
            m_projectToLandBoundary = false;
        }
    }
};

bool GridGeom::FlipEdges::Compute()
{

    bool successful = m_mesh->Administrate(Mesh::AdministrationOptions::AdministrateMeshEdgesAndFaces);

    if (m_triangulateFaces)
    {
        successful = successful && TriangulateFaces();
        successful = successful && m_mesh->Administrate(Mesh::AdministrationOptions::AdministrateMeshEdgesAndFaces);
    }

    if (!successful)
    {
        return false;
    }

    const int MaxIter = 10;
    const int numNodes = m_mesh->GetNumNodes();
    const int numEdges = m_mesh->GetNumEdges();
    std::vector<int> nodeMask(numNodes);
    int numFlippedEdges;
    for (int iter = 0; iter < MaxIter; iter++)
    {
        std::fill(nodeMask.begin(), nodeMask.end(), 0);
        numFlippedEdges = 0;

        for (int e = 0; e < numEdges; e++)
        {
            // check if nodes have been masked
            auto const firstNode = m_mesh->m_edges[e].first;
            auto const secondNode = m_mesh->m_edges[e].second;

            if (nodeMask[firstNode] != 0 || nodeMask[secondNode] != 0)
            {
                continue;
            }

            if (m_mesh->m_edgesNumFaces[e] != 2)
            {
                continue;
            }

            // triangles only
            auto const faceL = m_mesh->m_edgesFaces[e][0];
            auto const faceR = m_mesh->m_edgesFaces[e][1];

            const auto NumEdgesLeftFace = m_mesh->GetNumFaceEdges(faceL);
            const auto NumEdgesRightFace = m_mesh->GetNumFaceEdges(faceR);
            if (NumEdgesLeftFace != 3 || NumEdgesRightFace != 3)
            {
                return true;
            }

            int kl;
            int kr;
            int ntopo;

            if (e == 14)
            {
                std::cout << "Debug";
            }

            bool successful = TopologyFunctional(e, kl, kr, ntopo);

            if (nodeMask[kl] != 0 || nodeMask[kr] != 0)
            {
                continue;
            }


            if (ntopo < 0)
            {
                // Check if the quadrilateral composed by the two adjacent triangles is concave,
                // in which case the diagonals crosses
                Point intersection;
                double crossProduct;
                double firstRatio;
                double secondRatio;
                const auto areEdgesCrossing = AreLinesCrossing(m_mesh->m_nodes[firstNode],
                    m_mesh->m_nodes[secondNode],
                    m_mesh->m_nodes[kl],
                    m_mesh->m_nodes[kr],
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
                m_mesh->m_edges[e].first = kl;
                m_mesh->m_edges[e].second = kr;
                numFlippedEdges++;

                // Find the other edges
                int L1L;
                int L1R;
                int L2L;
                int L2R;
                for (int i = 0; i < NumEdgesLeftFace; i++)
                {
                    int edgeIndex = m_mesh->m_facesEdges[faceL][i];
                    if (edgeIndex == e)
                    {
                        continue;
                    }
                    const int firstNodeOtherEdge = m_mesh->m_edges[edgeIndex].first;
                    const int secondNodeOtherEdge = m_mesh->m_edges[edgeIndex].second;

                    if (firstNodeOtherEdge == firstNode || secondNodeOtherEdge == firstNode)
                    {
                        L1L = edgeIndex;
                    }
                    if (firstNodeOtherEdge == secondNode || secondNodeOtherEdge == secondNode)
                    {
                        L2L = edgeIndex;
                    }
                }

                for (int i = 0; i < NumEdgesRightFace; i++)
                {
                    int edgeIndex = m_mesh->m_facesEdges[faceR][i];
                    if (edgeIndex == e)
                    {
                        continue;
                    }
                    const int firstNodeOtherEdge = m_mesh->m_edges[edgeIndex].first;
                    const int secondNodeOtherEdge = m_mesh->m_edges[edgeIndex].second;

                    if (firstNodeOtherEdge == firstNode || secondNodeOtherEdge == firstNode)
                    {
                        L1R = edgeIndex;
                    }
                    if (firstNodeOtherEdge == secondNode || secondNodeOtherEdge == secondNode)
                    {
                        L2R = edgeIndex;
                    }
                }

                // change face orientation
                m_mesh->m_facesNodes[faceL][0] = kl;
                m_mesh->m_facesNodes[faceL][1] = kr;
                m_mesh->m_facesNodes[faceL][2] = firstNode;

                m_mesh->m_facesEdges[faceL][0] = e;
                m_mesh->m_facesEdges[faceL][1] = L1R;
                m_mesh->m_facesEdges[faceL][2] = L1L;

                m_mesh->m_facesNodes[faceR][0] = kl;
                m_mesh->m_facesNodes[faceR][1] = kr;
                m_mesh->m_facesNodes[faceR][2] = secondNode;

                m_mesh->m_facesEdges[faceR][0] = e;
                m_mesh->m_facesEdges[faceR][1] = L2R;
                m_mesh->m_facesEdges[faceR][2] = L2L;


                if (m_mesh->m_edgesFaces[L1R][0] == faceR)
                {
                    m_mesh->m_edgesFaces[L1R][0] = faceL;
                }
                else
                {
                    m_mesh->m_edgesFaces[L1R][1] = faceL;
                }


                if (m_mesh->m_edgesFaces[L2L][0] == faceL)
                {
                    m_mesh->m_edgesFaces[L2L][0] = faceR;
                }
                else
                {
                    m_mesh->m_edgesFaces[L2L][1] = faceR;
                }

                m_mesh->m_nodesNumEdges[firstNode] = m_mesh->m_nodesNumEdges[firstNode] - 1;
                m_mesh->m_nodesNumEdges[secondNode] = m_mesh->m_nodesNumEdges[secondNode] - 1;
                m_mesh->m_nodesNumEdges[kl] = m_mesh->m_nodesNumEdges[kl] + 1;
                m_mesh->m_nodesNumEdges[kr] = m_mesh->m_nodesNumEdges[kr] + 1;

                // Delete edge from m_mesh->m_nodesEdges[firstNode]
                DeleteEdgeFromNode(e, firstNode, secondNode);

                // Delete edge from m_mesh->m_nodesEdges[secondNode]
                DeleteEdgeFromNode(e, secondNode, firstNode);

                // Add edge to m_mesh->m_nodesEdges[kl]
                ResizeVectorIfNeeded(m_mesh->m_nodesNumEdges[kl], m_mesh->m_nodesEdges[kl]);
                for (int i = 0; i < m_mesh->m_nodesNumEdges[kl] - 1; i++)
                {
                    m_mesh->m_nodesEdges[kl][i] = m_mesh->m_nodesEdges[kl][i];
                }
                m_mesh->m_nodesEdges[kl][m_mesh->m_nodesNumEdges[kl]- 1] = e;
                m_mesh->SortEdgesInCounterClockWiseOrder(kl);

                // Add edge to m_mesh->m_nodesEdges[kr]
                ResizeVectorIfNeeded(m_mesh->m_nodesNumEdges[kr], m_mesh->m_nodesEdges[kr]);
                for (int i = 0; i < m_mesh->m_nodesNumEdges[kr] - 1; i++)
                {
                    m_mesh->m_nodesEdges[kr][i] = m_mesh->m_nodesEdges[kr][i];
                }
                m_mesh->m_nodesEdges[kr][m_mesh->m_nodesNumEdges[kr] - 1] = e;
                m_mesh->SortEdgesInCounterClockWiseOrder(kr);

            }
        }
    }

    if (numFlippedEdges > 0)
    {
        return false;
    }

    // Perform mesh administration
    m_mesh->Administrate(Mesh::AdministrationOptions::AdministrateMeshEdgesAndFaces);
    return true;
}

bool GridGeom::FlipEdges::DeleteEdgeFromNode(int edge, int firstNode, int secondNode) const
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
        return false;
    }

    int count = 0;
    for (int i = 0; i < m_mesh->m_nodesNumEdges[firstNode] + 1; i++)
    {
        if (i <= kk - 1)
        {
            m_mesh->m_nodesEdges[firstNode][count] = m_mesh->m_nodesEdges[firstNode][i];
            count++;
        }
        else if (i > kk)
        {
            m_mesh->m_nodesEdges[firstNode][count] = m_mesh->m_nodesEdges[firstNode][i];
            count++;
        }
    }
    ResizeVectorIfNeeded(m_mesh->m_nodesNumEdges[firstNode], m_mesh->m_nodesEdges[firstNode]);

    return true;
}

bool GridGeom::FlipEdges::TopologyFunctional(int e,
                                             int& kl,
                                             int& kr,
                                             int& ntopo) const
{

    if (m_mesh->m_edgesNumFaces[e] != 2)
    {
        return true;
    }

    const auto firstNode = m_mesh->m_edges[e].first;
    const auto secondNode = m_mesh->m_edges[e].second;
    const auto faceL = m_mesh->m_edgesFaces[e][0];
    const auto faceR = m_mesh->m_edgesFaces[e][1];
    const auto NumEdgesLeftFace = m_mesh->GetNumFaceEdges(faceL);
    const auto NumEdgesRightFace = m_mesh->GetNumFaceEdges(faceR);

    if (NumEdgesLeftFace != 3 || NumEdgesRightFace != 3)
    {
        return true;
    }

    // find the nodes that are connected to both k1 and k
    int sumIndexsesLeftFace = 0;
    int sumIndexsesRightFace = 0;
    for (int i = 0; i < 3; i++)
    {
        sumIndexsesLeftFace += m_mesh->m_facesNodes[faceL][i];
        sumIndexsesRightFace += m_mesh->m_facesNodes[faceR][i];
    }

    kl = sumIndexsesLeftFace - firstNode - secondNode;
    kr = sumIndexsesRightFace - firstNode - secondNode;

    if (kl < 0 || kr < 0)
    {
        return true;
    }

    // check that kl is part of faceL
    bool nodeFound = false;
    for (int i = 0; i < NumEdgesLeftFace; i++)
    {
        if (m_mesh->m_facesNodes[faceL][i] == kl)
        {
            nodeFound = true;
            break;
        }
    }

    if (!nodeFound)
    {
        return true;
    }

    // check that kr is part of faceR
    nodeFound = false;
    for (int i = 0; i < NumEdgesRightFace; i++)
    {
        if (m_mesh->m_facesNodes[faceR][i] == kr)
        {
            nodeFound = true;
            break;
        }
    }

    if (!nodeFound)
    {
        return true;
    }

    //  compute the change in functional
    int n1 = m_mesh->m_nodesNumEdges[firstNode] - OptimalNumberOfConnectedNodes(firstNode);
    int n2 = m_mesh->m_nodesNumEdges[secondNode] - OptimalNumberOfConnectedNodes(secondNode);
    int nL = m_mesh->m_nodesNumEdges[kl] - OptimalNumberOfConnectedNodes(kl);
    int nR = m_mesh->m_nodesNumEdges[kr] - OptimalNumberOfConnectedNodes(kr);

    if (m_projectToLandBoundary)
    {
        if (m_landBoundaries->m_meshNodesLandBoundarySegments[firstNode] >= 0 && m_landBoundaries->m_meshNodesLandBoundarySegments[secondNode] >= 0)
        {
            //edge is associated with a land boundary, keep the edge
            ntopo = 1000;
            return true;
        }
        else
        {
            auto n1L = DifferenceFromOptimum(firstNode, secondNode, kl);
            auto n1R = DifferenceFromOptimum(firstNode, secondNode, kr);

            auto n2R = DifferenceFromOptimum(secondNode, firstNode, kl);
            auto n2L = DifferenceFromOptimum(secondNode, firstNode, kr);

            auto nL = DifferenceFromOptimum(kl, firstNode, secondNode);
            auto nR = DifferenceFromOptimum(kr, firstNode, secondNode);

            ntopo = (n1L - 1) * (n1L - 1) +
                (n1R - 1) * (n1R - 1) +
                (n2L - 1) * (n2L - 1) +
                (n2R - 1) * (n2R - 1) +
                2.0 * ((nL + 1) * (nL + 1) + (nR + 1) * (nR + 1)) -
                (n1L * n1L + n1R * n1R + n2L * n2L + n2R * n2R + 2.0 * (nL * nL + nR * nR));

        }
    }
    else
    {
        ntopo = (n1 - 1) * (n1 - 1) +
                (n2 - 1) * (n2 - 1) +
                (nL + 1) * (nL + 1) +
                (nR + 1) * (nR + 1) -
                (n1 * n1 + n2 * n2 + nL * nL + nR * nR);
    }



    return true;
}

//comp_nnow
int GridGeom::FlipEdges::DifferenceFromOptimum(int nodeIndex, int firstNode, int secondNode) const
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
           m_mesh->m_edgesNumFaces[edgeIndex] > 1 &&
           currentEdgeIndexInNodeEdges != edgeIndexConnectingSecondNode)
    {
        currentEdgeIndexInNodeEdges = NextCircularBackwardIndex(currentEdgeIndexInNodeEdges, m_mesh->m_nodesNumEdges[nodeIndex]);
        edgeIndex = m_mesh->m_nodesEdges[nodeIndex][currentEdgeIndexInNodeEdges];
        otherNode = m_mesh->m_edges[edgeIndex].first + m_mesh->m_edges[edgeIndex].second - nodeIndex;
        num++;
    }

    int firstEdgeInPathIndex = -1;
    if (m_landBoundaries->m_meshNodesLandBoundarySegments[otherNode] >= 0 ||
        m_mesh->m_edgesNumFaces[edgeIndex] < 2)
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
        while ( m_landBoundaries->m_meshNodesLandBoundarySegments[otherNode] < 0 &&
                m_mesh->m_edgesNumFaces[edgeIndex] > 1 &&
                currentEdgeIndexInNodeEdges != edgeIndexConnectingFirstNode &&
                edgeIndex != firstEdgeInPathIndex )
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
            m_mesh->m_edgesNumFaces[edgeIndex] < 2) && edgeIndex!= firstEdgeInPathIndex)
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

int GridGeom::FlipEdges::OptimalNumberOfConnectedNodes(int index) const
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

bool GridGeom::FlipEdges::TriangulateFaces()
{
    for (int i = 0; i < m_mesh->GetNumFaces(); i++)
    {
        const auto NumEdges = m_mesh->GetNumFaceEdges(i);

        if (NumEdges < 4)
        {
            continue;
        }

        int indexFirstNode = m_mesh->m_facesNodes[i][0];
        for (int j = 2; j < NumEdges - 1; j++)
        {
            const auto nodeIndex = m_mesh->m_facesNodes[i][j];
            int newEdgeIndex;
            m_mesh->ConnectNodes(indexFirstNode, nodeIndex, newEdgeIndex);
        }
    }
    return true;
}