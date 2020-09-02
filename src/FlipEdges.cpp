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

GridGeom::FlipEdges::FlipEdges(Mesh* mesh, LandBoundaries* landBoundary) :
    m_mesh(mesh),
    m_landBoundaries(landBoundary)
{
    m_triangulateFaces = true;
    m_projectToLandBoundary = true;

};

bool GridGeom::FlipEdges::Compute()
{

    m_mesh->FindFaces();
    bool successful = m_mesh->ClassifyNodes();

    if (m_triangulateFaces)
    {
        successful = TriangulateFaces();
        m_mesh->FindFaces();
        successful = m_mesh->ClassifyNodes();
    }

    if (m_projectToLandBoundary)
    {
        m_landBoundaries->FindNearestMeshBoundary(4);
    }

    const int MaxIter = 10;
    const int numNodes = m_mesh->GetNumNodes();
    const int numEdges = m_mesh->GetNumEdges();
    int k1;
    int k2;
    int kl;
    int kr;
    int faceL;
    int faceR;
    int ntopo;

    std::vector<int> nodeMask(numNodes);

    for (int iter = 0; iter < MaxIter; iter++)
    {
        std::fill(nodeMask.begin(), nodeMask.end(), 0);

        for (int e = 0; e < numEdges; e++)
        {
            bool successful = TopologyFunctional(e,
                k1,
                k2,
                kl,
                kr,
                faceL,
                faceR,
                ntopo);




        }


    }




    return true;
}


bool GridGeom::FlipEdges::TopologyFunctional(int e,
    int& k1,
    int& k2,
    int& kl,
    int& kr,
    int& faceL,
    int& faceR,
    int& ntopo) const
{

    if (m_mesh->m_edgesNumFaces[e] != 2)
    {
        return true;
    }

    k1 = m_mesh->m_edges[e].first;
    k2 = m_mesh->m_edges[e].second;
    faceL = m_mesh->m_edgesFaces[e][0];
    faceR = m_mesh->m_edgesFaces[e][1];

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

    kl = sumIndexsesLeftFace - k1 - k2;
    kr = sumIndexsesRightFace - k1 - k2;

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
    int n1 = m_mesh->m_nodesNumEdges[k1] - OptimalNumberOfConnectedNodes(k1);
    int n2 = m_mesh->m_nodesNumEdges[k2] - OptimalNumberOfConnectedNodes(k2);
    int nL = m_mesh->m_nodesNumEdges[kl] - OptimalNumberOfConnectedNodes(kl);
    int nR = m_mesh->m_nodesNumEdges[kr] - OptimalNumberOfConnectedNodes(kr);

    ntopo = (n1 - 1) * (n1 - 1) +
        (n2 - 1) * (n2 - 1) +
        (nL + 1) * (nL + 1) +
        (nR + 1) * (nR + 1) -
        n1 * n1 + n2 * n2 + nL * nL + nR * nR;

    if (m_landBoundaries)
    {
        if (m_landBoundaries->m_meshNodesLandBoundarySegments[k1] >= 0 && m_landBoundaries->m_meshNodesLandBoundarySegments[k2] >= 0)
        {
            //edge is associated with a land boundary, keep the edge
            ntopo = 1000;
        }
        else
        {


        }

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
    bool isccw = sign < 0 ? true : false;

    if (isccw)
    {
        int firstNodeTemp = firstNode;
        firstNode = secondNode;
        secondNode = firstNodeTemp;
    }

    int optimalNumber = 0;
    int num = 0;
    int n = 0;

    // find the first edge connecting firstNode
    int indexFirstNode = -1;
    int edgeIndex = -1;
    for (int i = 0; i < m_mesh->m_nodesNumEdges[nodeIndex]; i++)
    {
        edgeIndex = m_mesh->m_nodesEdges[nodeIndex][i];

        if (m_mesh->m_edges[edgeIndex].first == firstNode || m_mesh->m_edges[edgeIndex].second == firstNode)
        {
            indexFirstNode = i;
            break;
        }
    }

    if (indexFirstNode == -1) 
    {
        return true;
    }

    if (m_mesh->m_edges[edgeIndex].first != nodeIndex && m_mesh->m_edges[edgeIndex].second != nodeIndex)
    {
        return true;
    }

    // find the first edge connecting secondNode
    int indexSecondNode = -1;
    edgeIndex = -1;
    for (int i = 0; i < m_mesh->m_nodesNumEdges[nodeIndex]; i++)
    {
        edgeIndex = m_mesh->m_nodesEdges[nodeIndex][i];

        if (m_mesh->m_edges[edgeIndex].first == secondNode || m_mesh->m_edges[edgeIndex].second == secondNode)
        {
            indexSecondNode = i;
            break;
        }
    }

    if (indexSecondNode == -1)
    {
        return true;
    }

    if (m_mesh->m_edges[edgeIndex].first != nodeIndex && m_mesh->m_edges[edgeIndex].second != nodeIndex)
    {
        return true;
    }

    // Start
    // count the numbers of edges clockwise from the one connecting k1 
    // that are not in a land/net boundary path
    int kk = indexFirstNode;
    edgeIndex = m_mesh->m_nodesEdges[nodeIndex][indexFirstNode];
    int otherNode = m_mesh->m_edges[edgeIndex].first + m_mesh->m_edges[edgeIndex].second - nodeIndex;
    num = 1;
    while ( m_landBoundaries->m_meshNodesLandBoundarySegments[otherNode]<0 && 
            kk != indexSecondNode &&
            m_mesh->m_edgesNumFaces[edgeIndex] > 1 )
    {
        kk = NextCircularBackwardIndex(kk, m_mesh->m_nodesNumEdges[nodeIndex]);
        edgeIndex = m_mesh->m_nodesEdges[nodeIndex][kk];
        int otherNode = m_mesh->m_edges[edgeIndex].first + m_mesh->m_edges[edgeIndex].second - nodeIndex;
        num++;
    }

    



};

int GridGeom::FlipEdges::OptimalNumberOfConnectedNodes(int index) const
{
    int optimalNumber = 6;
    if (m_mesh->m_nodeMask[index] == 1)
    {
        optimalNumber = 4;
    }
    if (m_mesh->m_nodeMask[index] == 3)
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