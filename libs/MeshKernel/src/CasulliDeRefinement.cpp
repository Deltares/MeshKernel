//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2024.
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

#include <algorithm>
#include <iostream> // REMOVE
using namespace std;

#include "MeshKernel/CasulliDeRefinement.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/UndoActions/AddNodeAction.hpp"
#include "MeshKernel/UndoActions/CompoundUndoAction.hpp"

std::unique_ptr<meshkernel::UndoAction> meshkernel::CasulliDeRefinement::Compute(Mesh2D& mesh)
{
    Polygons emptyPolygon;
    return Compute(mesh, emptyPolygon);
}

std::unique_ptr<meshkernel::UndoAction> meshkernel::CasulliDeRefinement::Compute(Mesh2D& mesh, const Polygons& polygon)
{
    std::vector<EdgeNodes> newNodes(mesh.GetNumEdges(), {constants::missing::uintValue, constants::missing::uintValue, constants::missing::uintValue, constants::missing::uintValue});
    std::vector<NodeMask> nodeMask(InitialiseNodeMask(mesh, polygon));
    std::unique_ptr<CompoundUndoAction> refinementAction = CompoundUndoAction::Create();

    [[maybe_unused]] UInt elementSeedIndex = FindElementSeedIndex(mesh, polygon);

    DoDeRefinement (mesh, polygon);

    // [[maybe_unsued]] const UInt numNodes = mesh.GetNumNodes();
    // [[maybe_unsued]] const UInt numEdges = mesh.GetNumEdges();
    // [[maybe_unsued]] const UInt numFaces = mesh.GetNumFaces();

    // refinementAction->Add(ComputeNewNodes(mesh, newNodes, nodeMask));
    // refinementAction->Add(ConnectNewNodes(mesh, newNodes, numNodes, numEdges, numFaces, nodeMask));
    // refinementAction->Add(Administrate(mesh, numNodes, nodeMask));
    return refinementAction;
}

void meshkernel::CasulliDeRefinement::FindSurroundingCells(const Mesh2D& mesh,
                                                           const Polygons& polygon [[maybe_unused]],
                                                           const UInt kCell,
                                                           const UInt nMax [[maybe_unused]],
                                                           UInt& nDirect,
                                                           UInt& nIndirect,
                                                           std::vector<UInt>& kDirect,
                                                           std::vector<UInt>& kIndirect,
                                                           std::vector<std::array<UInt, 2>>& kne)
{


    // Find directly connected cells

    kDirect.clear ();
    kIndirect.clear ();

    nDirect = 0;
    nIndirect = 0;

    for (UInt kk = 0; kk < mesh.m_numFacesNodes[kCell]; ++kk)
    {
        UInt L = mesh.m_facesEdges[kCell][kk];

        if (mesh.m_edgesNumFaces[L] < 2)
        {
            continue;
        }

        UInt kCell2 = mesh.m_edgesFaces[L][0] + mesh.m_edgesFaces[L][1] - kCell; // The other cell?
        bool alreadyVisited = false;

        for (UInt kkk = 0; kkk < kDirect.size (); ++kkk)
        {
            if (kDirect[kkk] == kCell2)
            {
                alreadyVisited = true;
                break;
            }

        }

        if (alreadyVisited)
        {
            continue;
        }

        kDirect.push_back(kCell2);
    }

    // find the cells indirectly connected cells

    for (UInt kk = 0; kk < mesh.m_numFacesNodes [kCell]; ++kk)
    {
        UInt k1 = mesh.m_facesNodes[kCell][kk];

        for (UInt kkk = 0; kkk < mesh.m_nodesNumEdges [k1]; ++kkk)
        {
            UInt L = mesh.m_nodesEdges [k1][kkk];
            bool isFound = false;

            for (UInt i = 0; i < mesh.m_edgesNumFaces [L]; ++i)
            {
                UInt kCell2 = mesh.m_edgesFaces [L][i];

                if (kCell == kCell2)
                {
                    continue;
                }

                isFound = false;

                for (UInt kkkk = 0; kkkk < kDirect.size (); ++kkkk)
                {
                    if (kCell2 == kDirect [kkkk])
                    {
                        isFound = true;
                        break;
                    }
                }

                if (isFound)
                {
                    continue;
                }

                isFound = false;

                for (UInt kkkk = 0; kkkk < kIndirect.size (); ++kkkk)
                {
                    if (kCell2 == kIndirect [kkkk])
                    {
                        isFound = true;
                        break;
                    }
                }

                if (isFound)
                {
                    continue;
                }

                // Add new cell
                kIndirect.push_back(kCell2);
            }
        }
    }

    // Find the adjacent cells

    for (UInt i = 0; i < kDirect.size (); ++i)
    {
        UInt kcell1 = kDirect [i];

        for (UInt j = 0; j < mesh.m_numFacesNodes [kcell1]; ++j)
        {
            UInt L = mesh.m_facesEdges [kcell1][j];

            if (mesh.m_edgesNumFaces [L] < 2)
            {
                continue;
            }

            UInt kCell2 = mesh.m_edgesFaces [L][0] + mesh.m_edgesFaces[L][1] - kcell1;

            for (UInt kk = 0; kk < kDirect.size (); ++kk)
            {
                if (kDirect[kk] == kCell2)
                {
                    if (kne[i][0] == 0)
                    {
                        kne[i][0] = -kCell2;
                    }
                    else
                    {
                        kne[2][1] = -kCell2;
                    }

                    kCell2 = constants::missing::uintValue;
                }
            }

            if (kCell2 == constants::missing::uintValue)
            {
                continue;
            }

            for (UInt kk = 0; kk < kIndirect.size (); ++kk)
            {
                if (kIndirect[kk] == kCell2)
                {
                    if (kne[i][0] == 0)
                    {
                        kne[i][0] = kCell2;
                    }
                    else
                    {
                        kne[2][1] = kCell2;
                    }

                }
           }
        }

    }
}

bool meshkernel::CasulliDeRefinement::ElementIsSeed(const Mesh2D& mesh, const Polygons& polygon [[maybe_unused]], const UInt face)
{
    bool isFace = true;

    for (UInt i = 0; i < mesh.m_numFacesNodes[face];++i)
    {

        if (mesh.m_nodesTypes[mesh.m_facesNodes [face][i]] == 0)
        {
            isFace = false;
            break;
        }

    }

    // for (UInt i = 0; i < mesh.m_facesNodes[face].size(); ++i)
    // {
    //     if (mesh.m_nodesTypes[mesh.m_facesNodes[face][i]] == 0)
    //     {
    //         isFace = false;
    //         break;
    //     }
    // }

    return isFace;
}

meshkernel::UInt meshkernel::CasulliDeRefinement::FindElementSeedIndex(const Mesh2D& mesh, const Polygons& polygon)
{
    UInt seedIndex = constants::missing::uintValue;

    for (UInt e = 0; e < mesh.Edges().size(); ++e)
    {
        if (mesh.m_edgesNumFaces[e] != 1)
        {
            const Edge& edge = mesh.GetEdge(e);

            if (mesh.m_nodesTypes[edge.first] != 2 || mesh.m_nodesTypes[edge.second] != 2)
            {
                continue;
            }

            UInt k1 = mesh.m_edgesFaces[e][0];

            if (mesh.m_facesNodes[k1].size() != constants::geometric::numNodesInQuadrilateral)
            {
                continue;
            }

            bool isFace = ElementIsSeed(mesh, polygon, k1);

            if (!isFace)
            {
                continue;
            }

            seedIndex = k1;
            break;
        }
    }

    // No seed index found, select the first quadrilateral inside the selecting polygon
    if (seedIndex == constants::missing::uintValue)
    {
        for (UInt face = 0; face < mesh.GetNumFaces(); ++face)
        {

            if (mesh.m_facesNodes[face].size() != constants::geometric::numNodesInQuadrilateral)
            {
                continue;
            }

            bool isFace = ElementIsSeed(mesh, polygon, face);

            if (!isFace)
            {
                continue;
            }

            seedIndex = face;
            break;
        }
    }

    // still node element found, so take the first
    if (seedIndex == constants::missing::uintValue)
    {
        seedIndex = 0;
    }

    return seedIndex;
}

void meshkernel::CasulliDeRefinement::UpdateFrontList(const Mesh& mesh, const std::vector<UInt>& frontList, std::vector<UInt>& frontListCopy, const UInt kNew)
{

    if (kNew != constants::missing::uintValue)
    {

        if (mesh.m_numFacesNodes [kNew] != 4)
        {
            return;
        }

        for (UInt i = 0; i < frontList.size (); ++i)
        {
            if (frontList [i] == kNew)
            {
                return;
            }
        }

        frontListCopy.push_back(kNew);
    }

}


void meshkernel::CasulliDeRefinement::DoDeRefinement(const Mesh2D& mesh, const Polygons& polygon)
{
    [[maybe_unused]] UInt seedElement = FindElementSeedIndex(mesh, polygon);
    [[maybe_unused]] UInt iterationCount = 0;
    [[maybe_unused]] UInt nMax = 10; // fix
    [[maybe_unused]] UInt numFront = 1; // fix
    [[maybe_unused]] UInt maxIterationCount = 10; // fix

    [[maybe_unused]] UInt maxNumFront = 10; // fix

    [[maybe_unused]] UInt nDirect = 0; // fix
    [[maybe_unused]] UInt nIndirect = 0; // fix
    [[maybe_unused]] std::vector<UInt> kDirect; // fix
    [[maybe_unused]] std::vector<UInt> kIndirect; // fix
    [[maybe_unused]] std::vector<std::array<UInt, 2>> kne(nMax); // fix
    [[maybe_unused]] std::vector<UInt> frontIndex; // fix
    [[maybe_unused]] std::vector<UInt> frontIndexCopy; // fix

    cout << "CasulliDeRefinement::DoDeRefinement "<< seedElement << endl;

    std::vector<ElementMask> cellMask(mesh.GetNumFaces(), ElementMask::Unassigned);

    cellMask[seedElement] = ElementMask::A;
    frontIndex.push_back(seedElement);

    while (frontIndex.size () > 0 && iterationCount < maxIterationCount)
    {
        ++iterationCount;
        frontIndexCopy.clear ();

        for (UInt i = 0; i < frontIndex.size (); ++i)
        {
            UInt k = frontIndex[i];
            UInt kOther = constants::missing::uintValue;

            FindSurroundingCells(mesh, polygon, k, nMax, nDirect, nIndirect, kDirect, kIndirect, kne);

            std::cout << " for element: " << k << " --- ";

            for (UInt jj = 0; jj < kDirect.size (); ++jj)
            {
                std::cout << kDirect [jj] << "  ";
            }

            std::cout << " --- ";

            for (UInt jj = 0; jj < kIndirect.size (); ++jj)
            {
                std::cout << kIndirect [jj] << "  ";
            }

            std::cout << std::endl;

            if (cellMask[k] == ElementMask::A)
            {
                std::cout << " cellMask[" << k << "] == ElementMask::A "<< std::endl;

                for (UInt j = 0; j < kDirect.size (); ++j)
                {
                    kOther = kDirect[j];

                    if (mesh.m_numFacesNodes[kOther] != 4)
                    {
                        continue;
                    }

                    if ((cellMask[kOther] != ElementMask::A && cellMask[kOther] != ElementMask::NotA) && (cellMask[kOther] != ElementMask::B && cellMask[kOther] != ElementMask::NotB))
                    {
                        cellMask[kOther] = ElementMask::B;
                        UpdateFrontList(mesh, frontIndex, frontIndexCopy, kOther);
                    }
                }

                for (UInt j = 0; j < kIndirect.size (); ++j)
                {
                    kOther = kIndirect[j];

                    if (mesh.m_numFacesNodes[kOther] != 4)
                    {
                        continue;
                    }

                    if (cellMask[kOther] != ElementMask::C)
                    {
                        cellMask[kOther] = ElementMask::C;
                    }
                }

                cellMask[kOther] = ElementMask::NotA;
            }
            else if (cellMask[k] == ElementMask::B)
            {
                std::cout << " cellMask[" << k << "] == ElementMask::B "<< std::endl;

                for (UInt j = 0; j < kDirect.size (); ++j)
                {
                    kOther = kDirect[j];

                    if (mesh.m_numFacesNodes[kOther] != 4)
                    {
                        continue;
                    }

                    if ((cellMask[kOther] != ElementMask::C) && (cellMask[kOther] != ElementMask::A && cellMask[kOther] != ElementMask::NotA) && (cellMask[kOther] != ElementMask::B && cellMask[kOther] != ElementMask::NotB))
                    {
                        cellMask[kOther] = ElementMask::A;
                        UpdateFrontList(mesh, frontIndex, frontIndexCopy, kOther);
                    }
                }

                for (UInt j = 0; j < kIndirect.size (); ++j)
                {
                    kOther = kIndirect[j];

                    if (mesh.m_numFacesNodes[kOther] != 4)
                    {
                        continue;
                    }

                    if ((cellMask[kOther] != ElementMask::B && cellMask[kOther] != ElementMask::NotB) && (cellMask[kOther] != ElementMask::A && cellMask[kOther] != ElementMask::NotA) && cellMask[kOther] != ElementMask::C)
                    {
                        cellMask[kOther] = ElementMask::B;
                        UpdateFrontList(mesh, frontIndex, frontIndexCopy, kOther);
                    }
                }

                cellMask[k] = ElementMask::NotB;
            }
        }

        frontIndex = frontIndexCopy;
        std::cout << "frontIndex.size " << frontIndex.size () << std::endl;
    }


    for (UInt i = 0; i < cellMask.size (); ++i)
    {
        std::cout << "cellMask "<< i << " = " << int(cellMask[i]) << std::endl;
    }

}

void meshkernel::CasulliDeRefinement::InitialiseBoundaryNodes(const Mesh2D& mesh, std::vector<NodeMask>& nodeMask)
{
    // Find nodes that lie on the boundary of the domain.
    for (UInt i = 0; i < mesh.GetNumEdges(); ++i)
    {
        UInt node1 = mesh.GetEdge(i).first;
        UInt node2 = mesh.GetEdge(i).second;

        if (mesh.m_edgesNumFaces[i] == 1)
        {
            if (nodeMask[node1] != NodeMask::Unassigned)
            {
                nodeMask[node1] = NodeMask::BoundaryNode;
            }

            if (nodeMask[node2] != NodeMask::Unassigned)
            {
                nodeMask[node2] = NodeMask::BoundaryNode;
            }
        }
    }
}

void meshkernel::CasulliDeRefinement::InitialiseCornerNodes(const Mesh2D& mesh, std::vector<NodeMask>& nodeMask)
{
    for (UInt i = 0; i < mesh.GetNumNodes(); ++i)
    {

        for (UInt j = 0; j < mesh.m_nodesNumEdges[i]; ++j)
        {
            UInt edge1 = mesh.m_nodesEdges[i][j];

            if (mesh.m_edgesNumFaces[edge1] != 1)
            {
                continue;
            }

            UInt elementId = mesh.m_edgesFaces[edge1][0];
            UInt nodeCount = mesh.m_numFacesNodes[elementId];

            UInt faceEdgeIndex = 0;
            UInt edge2 = mesh.m_facesEdges[elementId][faceEdgeIndex];

            // Check the loop termination, especially the faceEdgeIndex < nodeCount - 1
            // Perhaps change to for loop checking the condition then break.
            while (((mesh.GetEdge(edge2).first != i && mesh.GetEdge(edge2).second != i) || edge2 == edge1) && faceEdgeIndex < nodeCount - 1)
            {
                ++faceEdgeIndex;
                edge2 = mesh.m_facesEdges[elementId][faceEdgeIndex];
            }

            if (mesh.m_edgesNumFaces[edge2] == 1)
            {
                if (nodeMask[i] > NodeMask::Unassigned)
                {
                    nodeMask[i] = NodeMask::CornerNode;
                    break;
                }
            }
        }
    }

    // Find included corner nodes
    for (UInt i = 0; i < mesh.GetNumNodes(); ++i)
    {
        if (nodeMask[i] != NodeMask::Unassigned && mesh.m_nodesTypes[i] == 3)
        {
            nodeMask[i] = NodeMask::CornerNode;
        }
    }
}

void meshkernel::CasulliDeRefinement::InitialiseFaceNodes(const Mesh2D& mesh, std::vector<NodeMask>& nodeMask)
{

    std::vector<UInt> sharedFaces;
    std::vector<UInt> connectedNodes;
    std::vector<std::vector<UInt>> faceNodeMapping;

    for (UInt i = 0; i < mesh.GetNumNodes(); ++i)
    {
        if (nodeMask[i] == NodeMask::Unassigned)
        {
            continue;
        }

        if (mesh.m_nodesNumEdges[i] > 1)
        {
            FindPatchIds(mesh, i, sharedFaces, connectedNodes, faceNodeMapping);
        }
        else
        {
            nodeMask[i] = NodeMask::Unassigned;
            continue;
        }

        UInt elementCount = 0;

        for (UInt j = 0; j < sharedFaces.size(); ++j)
        {
            if (sharedFaces[j] != constants::missing::uintValue)
            {
                ++elementCount;
            }
        }

        if (elementCount == 0)
        {
            nodeMask[i] = NodeMask::Unassigned;
        }

        if (elementCount < mesh.m_nodesNumEdges[i] - 1 && nodeMask[i] == NodeMask::BoundaryNode)
        {
            nodeMask[i] = NodeMask::CornerNode;
        }

        if (elementCount > Mesh::m_maximumNumberOfEdgesPerNode && nodeMask[i] > NodeMask::Unassigned && nodeMask[i] < NodeMask::BoundaryNode)
        {
            nodeMask[i] = NodeMask::CornerNode;
        }
    }
}

std::vector<meshkernel::CasulliDeRefinement::NodeMask> meshkernel::CasulliDeRefinement::InitialiseNodeMask(const Mesh2D& mesh, const Polygons& polygon)
{
    std::vector<NodeMask> nodeMask(10 * mesh.GetNumNodes(), NodeMask::Unassigned);

    // Find nodes that are inside the polygon.
    // If the polygon is empty then all nodes will be taken into account.

    for (UInt i = 0; i < mesh.GetNumNodes(); ++i)
    {
        auto [containsPoint, pointIndex] = polygon.IsPointInPolygons(mesh.Node(i));

        if (containsPoint)
        {
            nodeMask[i] = NodeMask::RegisteredNode;
        }
    }

    InitialiseBoundaryNodes(mesh, nodeMask);
    InitialiseCornerNodes(mesh, nodeMask);
    InitialiseFaceNodes(mesh, nodeMask);

    return nodeMask;
}

void meshkernel::CasulliDeRefinement::FindPatchIds(const Mesh2D& mesh,
                                                   const UInt currentNode,
                                                   std::vector<UInt>& sharedFaces,
                                                   std::vector<UInt>& connectedNodes,
                                                   std::vector<std::vector<UInt>>& faceNodeMapping)
{
    sharedFaces.clear();
    connectedNodes.clear();
    faceNodeMapping.clear();

    if (currentNode >= mesh.GetNumNodes())
    {
        throw AlgorithmError("Node index out of range: {} >= {}", currentNode, mesh.GetNumNodes());
    }

    if (mesh.m_nodesNumEdges[currentNode] < 2)
    {
        return;
    }

    mesh.FindFacesConnectedToNode(currentNode, sharedFaces);

    // no shared face found
    if (sharedFaces.empty())
    {
        return;
    }

    mesh.GetConnectingNodes(currentNode, connectedNodes);
    mesh.FindNodesSharedByFaces(currentNode, sharedFaces, connectedNodes, faceNodeMapping);
}
