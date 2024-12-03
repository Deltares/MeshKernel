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
#include <tuple>

#include "MeshKernel/CasulliRefinement.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/UndoActions/FullUnstructuredGridUndo.hpp"

std::unique_ptr<meshkernel::UndoAction> meshkernel::CasulliRefinement::Compute(Mesh2D& mesh)
{
    Polygons emptyPolygon;
    return Compute(mesh, emptyPolygon);
}

std::unique_ptr<meshkernel::UndoAction> meshkernel::CasulliRefinement::Compute(Mesh2D& mesh, const Polygons& polygon)
{
    std::vector<EdgeNodes> newNodes(mesh.GetNumEdges(), {constants::missing::uintValue, constants::missing::uintValue, constants::missing::uintValue, constants::missing::uintValue});
    std::vector<NodeMask> nodeMask(InitialiseNodeMask(mesh, polygon));

    std::unique_ptr<FullUnstructuredGridUndo> refinementAction = FullUnstructuredGridUndo::Create(mesh);

    const UInt numNodes = mesh.GetNumNodes();
    const UInt numEdges = mesh.GetNumEdges();
    const UInt numFaces = mesh.GetNumFaces();

    ComputeNewNodes(mesh, newNodes, nodeMask);
    ConnectNewNodes(mesh, newNodes, numNodes, numEdges, numFaces, nodeMask);
    Administrate(mesh, numNodes, nodeMask);
    return refinementAction;
}

std::unique_ptr<meshkernel::UndoAction> meshkernel::CasulliRefinement::Compute(Mesh2D& mesh,
                                                                               const Polygons& polygon,
                                                                               const SampleInterpolator& interpolator,
                                                                               const int propertyId,
                                                                               const MeshRefinementParameters& refinementParameters)
{
    // refinementParameters.max_courant_time
    // refinementParameters.min_edge_size

    std::unique_ptr<FullUnstructuredGridUndo> refinementAction = FullUnstructuredGridUndo::Create(mesh);

    bool refinementRequested = true;
    int iterationCount = 0;

    while (refinementRequested && iterationCount < refinementParameters.max_num_refinement_iterations)
    {
        std::vector<EdgeNodes> newNodes(mesh.GetNumEdges(), {constants::missing::uintValue, constants::missing::uintValue, constants::missing::uintValue, constants::missing::uintValue});
        std::vector<NodeMask> nodeMask(InitialiseDepthBasedNodeMask(mesh, polygon, interpolator, propertyId, refinementParameters, refinementRequested));

        if (!refinementRequested)
        {
            break;
        }

        const UInt numNodes = mesh.GetNumNodes();
        const UInt numEdges = mesh.GetNumEdges();
        const UInt numFaces = mesh.GetNumFaces();

        ComputeNewNodes(mesh, newNodes, nodeMask);
        ConnectNewNodes(mesh, newNodes, numNodes, numEdges, numFaces, nodeMask);
        Administrate(mesh, numNodes, nodeMask);
        mesh.ComputeEdgesCenters();
        mesh.ComputeEdgesLengths();
        std::cout << "--------------------------------" << std::endl;
        ++iterationCount;
    }

    return refinementAction;
}

void meshkernel::CasulliRefinement::InitialiseBoundaryNodes(const Mesh2D& mesh, std::vector<NodeMask>& nodeMask)
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

void meshkernel::CasulliRefinement::InitialiseCornerNodes(const Mesh2D& mesh, std::vector<NodeMask>& nodeMask)
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

void meshkernel::CasulliRefinement::InitialiseFaceNodes(const Mesh2D& mesh, std::vector<NodeMask>& nodeMask)
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

void meshkernel::CasulliRefinement::RefineNodeMaskBasedOnDepths(const Mesh2D& mesh,
                                                                const SampleInterpolator& interpolator,
                                                                const int propertyId,
                                                                const MeshRefinementParameters& refinementParameters,
                                                                std::vector<NodeMask>& nodeMask [[maybe_unused]],
                                                                bool& refinementRequested)
{
    std::vector<double> interpolatedDepth(mesh.m_edgesCenters.size());
    interpolator.Interpolate(propertyId, mesh.m_edgesCenters, interpolatedDepth);
    const double maxDtCourant = refinementParameters.max_courant_time;

    refinementRequested = false;

    for (size_t i = 0; i < mesh.GetNumNodes(); ++i)
    {
        bool refineNode = false;

        std::cout << " waveCourant ";

        for (size_t j = 0; j < mesh.m_nodesEdges[i].size(); ++j)
        {
            UInt edgeId = mesh.m_nodesEdges[i][j];

            // If the current edge is already shorter than the minimum edge size then do not refine further
            if (mesh.m_edgeLengths[edgeId] <= refinementParameters.min_edge_size)
            {
                continue;
            }

            const double celerity = constants::physical::sqrt_gravity * std::sqrt(std::abs(interpolatedDepth[edgeId]));
            const double waveCourant = celerity * maxDtCourant / mesh.m_edgeLengths[edgeId];

            std::cout << waveCourant << "  ";

            if (waveCourant < 1.0)
            {
                refineNode = true;
            }
        }

        std::cout << "   " << std::boolalpha << refineNode << std::endl;

        if (nodeMask[i] == NodeMask::RegisteredNode && !refineNode)
        {
            nodeMask[i] = NodeMask::Unassigned;
        }
        else if (nodeMask[i] == NodeMask::RegisteredNode)
        {
            refinementRequested = true;
        }
    }
}

void meshkernel::CasulliRefinement::RegisterNodesInsidePolygon(const Mesh2D& mesh,
                                                               const Polygons& polygon,
                                                               std::vector<NodeMask>& nodeMask)
{
    for (UInt i = 0; i < mesh.GetNumNodes(); ++i)
    {
        auto [containsPoint, pointIndex] = polygon.IsPointInPolygons(mesh.Node(i));

        if (containsPoint)
        {
            nodeMask[i] = NodeMask::RegisteredNode;
        }
    }
}

std::vector<meshkernel::CasulliRefinement::NodeMask> meshkernel::CasulliRefinement::InitialiseDepthBasedNodeMask(const Mesh2D& mesh,
                                                                                                                 const Polygons& polygon,
                                                                                                                 const SampleInterpolator& interpolator,
                                                                                                                 const int propertyId,
                                                                                                                 const MeshRefinementParameters& refinementParameters,
                                                                                                                 bool& refinementRequested)
{
    std::vector<NodeMask> nodeMask(10 * mesh.GetNumNodes(), NodeMask::Unassigned);

    // Find nodes that are inside the polygon.
    // If the polygon is empty then all nodes will be taken into account.

    RegisterNodesInsidePolygon(mesh, polygon, nodeMask);
    RefineNodeMaskBasedOnDepths(mesh, interpolator, propertyId, refinementParameters, nodeMask, refinementRequested);
    InitialiseBoundaryNodes(mesh, nodeMask);
    InitialiseCornerNodes(mesh, nodeMask);
    InitialiseFaceNodes(mesh, nodeMask);

    return nodeMask;
}

std::vector<meshkernel::CasulliRefinement::NodeMask> meshkernel::CasulliRefinement::InitialiseNodeMask(const Mesh2D& mesh, const Polygons& polygon)
{
    std::vector<NodeMask> nodeMask(10 * mesh.GetNumNodes(), NodeMask::Unassigned);

    // Find nodes that are inside the polygon.
    // If the polygon is empty then all nodes will be taken into account.

    RegisterNodesInsidePolygon(mesh, polygon, nodeMask);
    InitialiseBoundaryNodes(mesh, nodeMask);
    InitialiseCornerNodes(mesh, nodeMask);
    InitialiseFaceNodes(mesh, nodeMask);

    return nodeMask;
}

void meshkernel::CasulliRefinement::Administrate(Mesh2D& mesh, const UInt numNodes, const std::vector<NodeMask>& nodeMask)
{
    // Need check only the original nodes in the mesh, hence use of numNodes.
    for (UInt i = 0; i < numNodes; ++i)
    {
        if (nodeMask[i] > NodeMask::Unassigned && nodeMask[i] < NodeMask::CornerNode)
        {
            // With passing false for the collectUndo, DeleteNode should return a null pointer
            [[maybe_unused]] auto undoAction = mesh.DeleteNode(i, false /* collectUndo info */);
        }
    }

    mesh.DeleteInvalidNodesAndEdges ();
    mesh.Administrate();
}

void meshkernel::CasulliRefinement::FindPatchIds(const Mesh2D& mesh,
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

std::vector<meshkernel::UInt> meshkernel::CasulliRefinement::GetNodesToConnect(const Mesh2D& mesh,
                                                                               const std::vector<NodeMask>& nodeMask,
                                                                               const std::vector<UInt>& newEdges,
                                                                               const std::vector<EdgeNodes>& newNodes,
                                                                               const UInt edgeCount,
                                                                               const UInt nodeIndex)
{
    std::vector<UInt> nodesToConnect(edgeCount, constants::missing::uintValue);

    for (UInt j = 0; j < edgeCount; ++j)
    {
        if (mesh.GetEdge(newEdges[j]).first == nodeIndex && nodeMask[newNodes[newEdges[j]][0]] == NodeMask::NewGeneralNode)
        {
            nodesToConnect[j] = newNodes[newEdges[j]][0];
        }

        if (mesh.GetEdge(newEdges[j]).first == nodeIndex && nodeMask[newNodes[newEdges[j]][2]] == NodeMask::NewGeneralNode)
        {
            nodesToConnect[j] = newNodes[newEdges[j]][2];
        }

        if (mesh.GetEdge(newEdges[j]).second == nodeIndex && nodeMask[newNodes[newEdges[j]][1]] == NodeMask::NewGeneralNode)
        {
            nodesToConnect[j] = newNodes[newEdges[j]][1];
        }

        if (mesh.GetEdge(newEdges[j]).second == nodeIndex && nodeMask[newNodes[newEdges[j]][3]] == NodeMask::NewGeneralNode)
        {
            nodesToConnect[j] = newNodes[newEdges[j]][3];
        }
    }

    return nodesToConnect;
}

void meshkernel::CasulliRefinement::ConnectNodes(Mesh2D& mesh, const std::vector<EdgeNodes>& newNodes, const UInt numEdges)
{
    // make the original-edge based new edges
    for (UInt i = 0; i < numEdges; ++i)
    {
        const UInt node1 = newNodes[i][0];
        const UInt node2 = newNodes[i][1];
        const UInt node3 = newNodes[i][2];
        const UInt node4 = newNodes[i][3];

        // Parallel edges, these are the start-end connections
        if (node1 != constants::missing::uintValue && node2 != constants::missing::uintValue && node1 != node2)
        {
            // With passing false for the collectUndo, ConnectNodes should return a null pointer
            [[maybe_unused]] auto ignore = mesh.ConnectNodes(node1, node2, false /* collectUndo info */);
        }

        if (node3 != constants::missing::uintValue && node4 != constants::missing::uintValue && node3 != node4)
        {
            [[maybe_unused]] auto ignore = mesh.ConnectNodes(node3, node4, false /* collectUndo info */);
        }

        // normal edges, these are the left-right connections
        if (node1 != constants::missing::uintValue && node3 != constants::missing::uintValue && node1 != node3)
        {
            [[maybe_unused]] auto ignore = mesh.ConnectNodes(node1, node3, false /* collectUndo info */);
        }

        if (node2 != constants::missing::uintValue && node4 != constants::missing::uintValue && node2 != node4)
        {
            [[maybe_unused]] auto ignore = mesh.ConnectNodes(node2, node4, false /* collectUndo info */);
        }
    }
}

void meshkernel::CasulliRefinement::ConnectFaceNodes(Mesh2D& mesh, const UInt currentFace, const std::vector<EdgeNodes>& newNodes,
                                                     std::vector<NodeMask>& nodeMask)
{

    // Perhaps change quads to maximum number of edges for any shape
    std::array<UInt, constants::geometric::numNodesInQuadrilateral> oldIndex{constants::missing::uintValue, constants::missing::uintValue,
                                                                             constants::missing::uintValue, constants::missing::uintValue};

    std::array<UInt, constants::geometric::numNodesInQuadrilateral> newIndex{constants::missing::uintValue, constants::missing::uintValue,
                                                                             constants::missing::uintValue, constants::missing::uintValue};
    // find the old and new nodes
    for (UInt j = 0; j < mesh.m_numFacesNodes[currentFace]; ++j)
    {
        const UInt previousIndex = (j == 0 ? (mesh.m_numFacesNodes[currentFace] - 1) : (j - 1));

        const UInt edgeId = mesh.m_facesEdges[currentFace][j];
        const UInt previousEdgeId = mesh.m_facesEdges[currentFace][previousIndex];

        oldIndex[j] = mesh.GetEdge(edgeId).first;
        newIndex[j] = newNodes[edgeId][2];

        if (oldIndex[j] != mesh.GetEdge(previousEdgeId).first && oldIndex[j] != mesh.GetEdge(previousEdgeId).second)
        {
            oldIndex[j] = mesh.GetEdge(edgeId).second;
            newIndex[j] = newNodes[edgeId][1];
        }
    }

    for (UInt j = 0; j < mesh.m_numFacesNodes[currentFace]; ++j)
    {
        const UInt previousIndex = (j == 0 ? (mesh.m_numFacesNodes[currentFace] - 1) : (j - 1));
        const UInt nextIndex = (j == mesh.m_numFacesNodes[currentFace] - 1 ? 0 : (j + 1));
        UInt nextNextIndex;

        if (j == mesh.m_numFacesNodes[currentFace] - 2)
        {
            nextNextIndex = 0;
        }
        else if (j == mesh.m_numFacesNodes[currentFace] - 1)
        {
            nextNextIndex = 1;
        }
        else
        {
            nextNextIndex = j + 2;
        }

        const UInt node1 = newIndex[j];
        const UInt node2 = oldIndex[previousIndex];
        const UInt node3 = oldIndex[nextIndex];
        const UInt node4 = oldIndex[nextNextIndex];

        // only one new node: new diagonal edge connects new node with one old node
        if (nodeMask[node1] < NodeMask::Unassigned && nodeMask[node2] == NodeMask::Unassigned && nodeMask[node3] == NodeMask::Unassigned && nodeMask[node4] == NodeMask::Unassigned)
        {
            // With passing false for the collectUndo, ConnectNodes should return a null pointer
            [[maybe_unused]] auto ignore = mesh.ConnectNodes(node1, node4, false /* collectUndo info */);
            break;
        }

        // only one old node: new diagonal edge connects new nodes only (i.e. perpendicular to previous one)
        if (nodeMask[node1] < NodeMask::Unassigned && nodeMask[node2] > NodeMask::Unassigned && nodeMask[node3] > NodeMask::Unassigned && nodeMask[node4] == NodeMask::Unassigned)
        {
            [[maybe_unused]] auto ignore = mesh.ConnectNodes(newIndex[previousIndex], newIndex[nextIndex], false /* collectUndo info */);
            break;
        }

        // two new and opposing nodes: new diagonal edge connects the new nodes
        if (nodeMask[node1] < NodeMask::Unassigned && nodeMask[node2] == NodeMask::Unassigned && nodeMask[node3] == NodeMask::Unassigned && nodeMask[node4] == NodeMask::RegisteredNode)
        {
            [[maybe_unused]] auto ignore = mesh.ConnectNodes(node1, newIndex[nextNextIndex], false /* collectUndo info */);
            break;
        }
    }
}

void meshkernel::CasulliRefinement::ConnectEdges(Mesh2D& mesh, const UInt currentNode, const std::vector<EdgeNodes>& newNodes, UInt& edgeCount, std::vector<UInt>& newEdges)
{
    std::ranges::fill(newEdges, constants::missing::uintValue);
    edgeCount = 0;

    for (UInt j = 0; j < mesh.m_nodesNumEdges[currentNode]; ++j)
    {
        UInt edgeId = mesh.m_nodesEdges[currentNode][j];

        if (mesh.m_edgesNumFaces[edgeId] == 1)
        {
            if (edgeCount >= newEdges.size())
            {
                newEdges.resize(2 * edgeCount + 1);
            }

            newEdges[edgeCount] = edgeId;
            ++edgeCount;
        }
        else
        {
            if (mesh.GetEdge(edgeId).first == currentNode)
            {
                // With passing false for the collectUndo, ConnectNodes should return a null pointer
                [[maybe_unused]] auto ignore1 = mesh.ConnectNodes(currentNode, newNodes[edgeId][0], false /* collectUndo info */);
                [[maybe_unused]] auto ignore2 = mesh.ConnectNodes(currentNode, newNodes[edgeId][2], false /* collectUndo info */);
            }
            else
            {
                [[maybe_unused]] auto ignore1 = mesh.ConnectNodes(currentNode, newNodes[edgeId][1], false /* collectUndo info */);
                [[maybe_unused]] auto ignore2 = mesh.ConnectNodes(currentNode, newNodes[edgeId][3], false /* collectUndo info */);
            }
        }
    }
}

void meshkernel::CasulliRefinement::ConnectNodes(Mesh2D& mesh,
                                                 const NodeMask nodeMask,
                                                 const std::vector<UInt>& nodesToConnect,
                                                 const UInt edgeCount,
                                                 const UInt nodeIndex)
{

    if (nodeMask != NodeMask::CornerNode)
    {
        if (edgeCount != 2)
        {
            throw AlgorithmError("Incorrect number of edges found: {}", edgeCount);
        }
        else
        {
            if (nodesToConnect[0] != constants::missing::uintValue && nodesToConnect[1] != constants::missing::uintValue && nodesToConnect[0] != nodesToConnect[1])
            {
                [[maybe_unused]] auto [edgeId, connectionAction] = mesh.ConnectNodes(nodesToConnect[0], nodesToConnect[1]);
            }
        }
    }
    else
    {
        for (UInt j = 0; j < edgeCount; ++j)
        {
            if (nodesToConnect[j] != constants::missing::uintValue && nodesToConnect[j] != nodeIndex)
            {
                [[maybe_unused]] auto [edgeId, connectionAction] = mesh.ConnectNodes(nodeIndex, nodesToConnect[j]);
            }
        }
    }
}

void meshkernel::CasulliRefinement::CreateMissingBoundaryEdges(Mesh2D& mesh,
                                                               const UInt numNodes,
                                                               const std::vector<EdgeNodes>& newNodes,
                                                               std::vector<NodeMask>& nodeMask)
{
    std::vector<UInt> newEdges(InitialEdgeArraySize);

    // make the missing boundary edges
    for (UInt i = 0; i < numNodes; ++i)
    {
        if (nodeMask[i] < NodeMask::BoundaryNode)
        {
            // boundary and kept nodes only
            continue;
        }

        UInt edgeCount = 0;
        ConnectEdges(mesh, i, newNodes, edgeCount, newEdges);

        if (edgeCount == 0)
        {
            continue;
        }

        std::vector<UInt> nodesToConnect(GetNodesToConnect(mesh, nodeMask, newEdges, newNodes, edgeCount, i));
        ConnectNodes(mesh, nodeMask[i], nodesToConnect, edgeCount, i);
    }
}

void meshkernel::CasulliRefinement::ConnectNewNodes(Mesh2D& mesh, const std::vector<EdgeNodes>& newNodes, const UInt numNodes, const UInt numEdges, const UInt numFaces, std::vector<NodeMask>& nodeMask)
{
    ConnectNodes(mesh, newNodes, numEdges);

    // create the diagonal edges in quads that connect the new mesh with the old mesh
    for (UInt i = 0; i < numFaces; ++i)
    {

        if (mesh.m_numFacesNodes[i] != constants::geometric::numNodesInQuadrilateral)
        {
            continue;
        }

        bool faceIsActive = true;

        for (UInt j = 0; j < mesh.m_facesNodes[i].size(); ++j)
        {
            if (nodeMask[mesh.m_facesNodes[i][j]] == NodeMask::Unassigned)
            {
                faceIsActive = false;
                break;
            }
        }

        // Check for active nodes
        if (faceIsActive)
        {
            ConnectFaceNodes(mesh, i, newNodes, nodeMask);
        }
    }

    CreateMissingBoundaryEdges(mesh, numNodes, newNodes, nodeMask);

    for (UInt i = 0; i < numNodes; ++i)
    {
        if (nodeMask[i] != NodeMask::CornerNode || mesh.m_nodesNumEdges[i] <= MaximumNumberOfNodesInNewlyCreatedElements)
        {
            continue;
        }

        for (UInt j = 0; j < mesh.m_nodesNumEdges[i]; ++j)
        {
            const UInt edgeId = mesh.m_nodesEdges[i][j];

            if (mesh.GetEdge(edgeId).first == i)
            {
                // With passing false for the collectUndo, ConnectNodes should return a null pointer
                [[maybe_unused]] auto ignore1 = mesh.ConnectNodes(i, newNodes[edgeId][0], false /* collectUndo info */);
                [[maybe_unused]] auto ignore2 = mesh.ConnectNodes(i, newNodes[edgeId][2], false /* collectUndo info */);
            }
            else
            {
                [[maybe_unused]] auto ignore1 = mesh.ConnectNodes(i, newNodes[edgeId][1], false /* collectUndo info */);
                [[maybe_unused]] auto ignore2 = mesh.ConnectNodes(i, newNodes[edgeId][3], false /* collectUndo info */);
            }
        }
    }
}

void meshkernel::CasulliRefinement::ComputeNewFaceNodes(Mesh2D& mesh, std::vector<EdgeNodes>& newNodes, std::vector<NodeMask>& nodeMask)
{
    for (UInt i = 0; i < mesh.GetNumFaces(); ++i)
    {
        const Point elementCentre = mesh.m_facesCircumcenters[i];

        for (UInt j = 0; j < mesh.m_numFacesNodes[i]; ++j)
        {
            const UInt elementNode = mesh.m_facesNodes[i][j];

            UInt firstEdgeId = constants::missing::uintValue;
            UInt secondEdgeId = constants::missing::uintValue;
            UInt newNodeId = constants::missing::uintValue;
            std::unique_ptr<AddNodeAction> nodeInsertionAction;

            for (UInt k = 0; k < mesh.m_facesEdges[i].size(); ++k)
            {
                UInt edgeId = mesh.m_facesEdges[i][k];

                if (mesh.GetEdge(edgeId).first == elementNode || mesh.GetEdge(edgeId).second == elementNode)
                {
                    if (firstEdgeId == constants::missing::uintValue)
                    {
                        firstEdgeId = edgeId;
                    }
                    else
                    {
                        secondEdgeId = edgeId;
                        break;
                    }
                }
            }

            if (firstEdgeId == constants::missing::uintValue || secondEdgeId == constants::missing::uintValue)
            {
                // No edges found
                continue;
            }

            if (nodeMask[elementNode] > NodeMask::Unassigned)
            {
                Point newNode = 0.5 * (elementCentre + mesh.Node(elementNode));

                std::tie(newNodeId, nodeInsertionAction) = mesh.InsertNode(newNode);

                nodeMask[newNodeId] = NodeMask::NewAssignedNode;
            }
            else
            {
                newNodeId = elementNode;
            }

            StoreNewNode(mesh, elementNode, firstEdgeId, secondEdgeId, newNodeId, newNodes);
        }
    }
}

void meshkernel::CasulliRefinement::ComputeNewEdgeNodes(Mesh2D& mesh, const UInt numEdges, std::vector<EdgeNodes>& newNodes, std::vector<NodeMask>& nodeMask)
{

    for (UInt i = 0; i < numEdges; ++i)
    {
        UInt newNodeId = constants::missing::uintValue;
        std::unique_ptr<AddNodeAction> nodeInsertionAction;

        if (mesh.m_edgesNumFaces[i] != 1)
        {
            continue;
        }

        const UInt node1 = mesh.GetEdge(i).first;
        const UInt node2 = mesh.GetEdge(i).second;

        if (node1 == constants::missing::uintValue && node2 == constants::missing::uintValue)
        {
            continue;
        }

        const Point edgeCentre = 0.5 * (mesh.Node(node1) + mesh.Node(node2));

        if (nodeMask[node1] != NodeMask::Unassigned)
        {
            const Point newNode = 0.5 * (edgeCentre + mesh.Node(node1));

            std::tie(newNodeId, nodeInsertionAction) = mesh.InsertNode(newNode);
            nodeMask[newNodeId] = NodeMask::NewGeneralNode;
        }
        else
        {
            newNodeId = node1;
        }

        StoreNewNode(mesh, node1, i, i, newNodeId, newNodes);

        if (nodeMask[node2] != NodeMask::Unassigned)
        {
            const Point newNode = 0.5 * (edgeCentre + mesh.Node(node2));

            std::tie(newNodeId, nodeInsertionAction) = mesh.InsertNode(newNode);
            nodeMask[newNodeId] = NodeMask::NewGeneralNode;
        }
        else
        {
            newNodeId = node2;
        }

        StoreNewNode(mesh, node2, i, i, newNodeId, newNodes);
    }
}

void meshkernel::CasulliRefinement::ComputeNewNodes(Mesh2D& mesh, std::vector<EdgeNodes>& newNodes, std::vector<NodeMask>& nodeMask)
{
    // Keep copy of number of edges in mesh before any nodes are added.
    const UInt numEdges = mesh.GetNumEdges();

    ComputeNewFaceNodes(mesh, newNodes, nodeMask);
    ComputeNewEdgeNodes(mesh, numEdges, newNodes, nodeMask);
}

void meshkernel::CasulliRefinement::StoreNewNode(const Mesh2D& mesh, const UInt nodeId, const UInt edge1Index, const UInt edge2Index, const UInt newNodeId, std::vector<EdgeNodes>& newNodes)
{
    UInt edgeId1 = edge1Index;
    UInt edgeId2 = edge2Index;

    if (edgeId1 != constants::missing::uintValue)
    {
        if (edgeId2 == constants::missing::uintValue)
        {
            edgeId2 = edgeId1;
        }
    }
    else
    {
        if (edgeId2 != constants::missing::uintValue)
        {
            edgeId1 = edgeId2;
        }
        else
        {
            throw AlgorithmError("Node edges specified: {}, {}", edgeId1, edgeId2);
        }
    }

    UInt elementId = mesh.FindCommonFace(edgeId1, edgeId2);

    if (elementId == constants::missing::uintValue)
    {
        throw AlgorithmError("No element found that shares edge: {} and {}", edgeId1, edgeId2);
    }

    UInt lr1 = mesh.IsLeftOrRight(elementId, edgeId1);
    UInt lr2 = mesh.IsLeftOrRight(elementId, edgeId2);

    const UInt se1 = mesh.IsStartOrEnd(edgeId1, nodeId);
    const UInt se2 = mesh.IsStartOrEnd(edgeId2, nodeId);

    if (edgeId1 == edgeId2)
    {
        lr1 = 1 - lr1;
        lr2 = 1 - lr2;
    }

    const UInt iPoint1 = se1 + 2 * (1 - lr1);
    const UInt iPoint2 = se2 + 2 * (1 - lr2);

    if (newNodes[edgeId1][iPoint1] == constants::missing::uintValue)
    {
        newNodes[edgeId1][iPoint1] = newNodeId;
    }

    if (newNodes[edgeId2][iPoint2] == constants::missing::uintValue)
    {
        newNodes[edgeId2][iPoint2] = newNodeId;
    }
}
