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

    [[maybe_unused]] UInt edgeSeedIndex = FineEdgeSeedIndex(mesh, polygon);

    // [[maybe_unsued]] const UInt numNodes = mesh.GetNumNodes();
    // [[maybe_unsued]] const UInt numEdges = mesh.GetNumEdges();
    // [[maybe_unsued]] const UInt numFaces = mesh.GetNumFaces();

    // refinementAction->Add(ComputeNewNodes(mesh, newNodes, nodeMask));
    // refinementAction->Add(ConnectNewNodes(mesh, newNodes, numNodes, numEdges, numFaces, nodeMask));
    // refinementAction->Add(Administrate(mesh, numNodes, nodeMask));
    return refinementAction;
}

void meshkernel::CasulliDeRefinement::FindSurroundingCells(const Mesh2D& mesh, const Polygons& polygon [[maybe_unused]],
                                                           const UInt kCell, const UInt nMax,
                                                           UInt& nDirect, UInt& nIndirect,
                                                           std::vector<UInt>& kDirect,
                                                           std::vector<UInt>& kIndirect,
                                                           std::vector<Edge>& kne)
{

    // Find directly connected cells

    nDirect = 0;
    nIndirect = 0;
    std::ranges::fill(kDirect, 0);
    std::ranges::fill(kIndirect, 0);

    for (UInt kk = 0; kk < mesh.m_numFacesNodes[kCell]; ++kk)
    {
        UInt L = mesh.m_facesEdges[kk];

        if (mesh.m_edgesNumFaces[L] < 2)
        {
            continue;
        }

        UInt kcell2 = mesh.m_edgesFaces[L].first + mesh.m_edgesFaces[L].second - kCell; // The other cell?
    }
}

bool meshkernel::CasulliDeRefinement::ElementIsSeed(const Mesh2D& mesh, const Polygons& polygon [[maybe_unused]], const UInt face)
{
    bool isFace = true;

    for (UInt i = 0; i < mesh.m_facesNodes[face].size(); ++i)
    {
        if (mesh.m_nodesTypes[mesh.m_facesNodes[face][i]] == 0)
        {
            isFace = false;
            break;
        }
    }

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

            if (mesh.m_nodesTypes[edge.first] != 2 or mesh.m_nodesTypes[edge.second] != 2)
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

void meshkernel::CasulliDeRefinement::DoDeRefinement(const Mesh2D& mesh, const Polygons& polygon)
{
    UInt seedElement = FindElementSeedIndex(mesh, polygon);
    UInt iterationCount = 0;

    std::vector<ElementMask> elementMask(mesh.GetNumFaces(), ElementMask::Unassigned);

    elementMask[seedElement] = ElementMask::A;

    while (numFront > 0 && iterationCount < maxIterationCount)
    {
        ++iterationCount;
        numFrontNew = 0;

        for (UInt i = 1; i <= numFront; ++i)
        {
            UInt k = frontIndex[i];

            FindSurroundElements(k, NMAX, ndirect, ndirect, kDirect, kDirect, kne);

            if (cellMask[k] == ElementMask::A)
            {
                for (UINt j = 1; j <= nDirect; ++j)
                {
                    kOther = kDirect[j];

                    if (mesh.m_numFacesNodes[kOther] != 4)
                    {
                        continue;
                    }

                    if ((cellMask[kOther] != ElementMask::A or cellMask[kOther] != ElementMask::NotA) && (cellMask[kOther] != ElementMask::B or cellMask[kOther] != ElementMask::NotB))
                    {
                        cellMask[kOther] = ElementMask::C;
                        UpdateFrontList(kOther);
                    }
                }

                for (UINt j = 1; j <= nIndirect; ++j)
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
                for (UINt j = 1; j <= nDirect; ++j)
                {
                    kOther = kDirect[j];

                    if (mesh.m_numFacesNodes[kOther] != 4)
                    {
                        continue;
                    }

                    if ((cellMask[kOther] != ElementMask::C) && (cellMask[kOther] != ElementMask::B or cellMask[kOther] != ElementMask::NotB))
                    {
                        cellMask[kOther] = ElementMask::C;
                        UpdateFrontList(kOther);
                    }
                }

                for (UINt j = 1; j <= nIndirect; ++j)
                {
                    kOther = kIndirect[j];

                    if (mesh.m_numFacesNodes[kOther] != 4)
                    {
                        continue;
                    }

                    if (cellMask[kOther] != ElementMask::B && (cellMask[kOther] != ElementMask::A || cellMask[kOther] != ElementMask::NotA) && cellMask[kOther] != ElementMask::C)
                    {
                        cellMask[kOther] = ElementMask::B;
                        UpdateFrontList(kOther);
                    }
                }

                cellMask[kOther] = ElementMask::NotB;
            }
        }
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
