#include <algorithm>

#include "MeshKernel/CasulliRefinement.hpp"
#include "MeshKernel/Exceptions.hpp"

void meshkernel::CasulliRefinement::Compute(Mesh2D& mesh)
{
    Polygons emptyPolygon;
    Compute(mesh, emptyPolygon);
}

void meshkernel::CasulliRefinement::Compute(Mesh2D& mesh, const Polygons& polygon)
{
    std::vector<EdgeNodes> newNodes(mesh.GetNumEdges(), {constants::missing::uintValue, constants::missing::uintValue, constants::missing::uintValue, constants::missing::uintValue});
    std::vector<NodeMask> nodeMask(InitialiseNodeMask(mesh, polygon));

    const UInt numNodes = mesh.GetNumNodes();
    const UInt numEdges = mesh.GetNumEdges();
    const UInt numFaces = mesh.GetNumFaces();

    ComputeNewNodes(mesh, newNodes, nodeMask);
    ConnectNewNodes(mesh, newNodes, numNodes, numEdges, numFaces, nodeMask);
    Administrate(mesh, numNodes, nodeMask);
}

std::vector<meshkernel::CasulliRefinement::NodeMask> meshkernel::CasulliRefinement::InitialiseNodeMask(const Mesh2D& mesh, const Polygons& polygon)
{
    std::vector<NodeMask> nodeMask(10 * mesh.GetNumNodes(), NodeMask::Unassigned);

    // Find nodes that are inside the polygon.
    // If the polygon is empty then all nodes will be taken into account.
    for (UInt i = 0; i < mesh.GetNumNodes(); ++i)
    {
        auto [containsPoint, pointIndex] = polygon.IsPointInPolygons(mesh.m_nodes[i]);

        if (containsPoint)
        {
            nodeMask[i] = NodeMask::RegisteredNode;
        }
    }

    // Find nodes that lie on the boundary of the domain.
    for (UInt i = 0; i < mesh.GetNumEdges(); ++i)
    {
        UInt node1 = mesh.m_edges[i].first;
        UInt node2 = mesh.m_edges[i].second;

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
            while (((mesh.m_edges[edge2].first != i && mesh.m_edges[edge2].second != i) || edge2 == edge1) && faceEdgeIndex < nodeCount - 1)
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

    return nodeMask;
}

void meshkernel::CasulliRefinement::Administrate(Mesh2D& mesh, const UInt numNodes, const std::vector<NodeMask>& nodeMask)
{
    // Need check only the original nodes in the mesh, hence use of numNodes.
    for (UInt i = 0; i < numNodes; ++i)
    {
        if (nodeMask[i] > NodeMask::Unassigned && nodeMask[i] < NodeMask::CornerNode)
        {
            mesh.DeleteNode(i);
        }
    }

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

void meshkernel::CasulliRefinement::ConnectNewNodes(Mesh2D& mesh, const std::vector<EdgeNodes>& newNodes, const UInt numNodes, const UInt numEdges, const UInt numFaces, std::vector<NodeMask>& nodeMask)
{
    //  make the original-link based new links
    for (UInt i = 0; i < numEdges; ++i)
    {
        const UInt node1 = newNodes[i][0];
        const UInt node2 = newNodes[i][1];
        const UInt node3 = newNodes[i][2];
        const UInt node4 = newNodes[i][3];

        // Parallel edges, these are the start-end connections
        if (node1 != constants::missing::uintValue && node2 != constants::missing::uintValue && node1 != node2)
        {
            mesh.ConnectNodes(node1, node2);
        }

        if (node3 != constants::missing::uintValue && node4 != constants::missing::uintValue && node3 != node4)
        {
            mesh.ConnectNodes(node3, node4);
        }

        // normal edges, these are the left-right connections
        if (node1 != constants::missing::uintValue && node3 != constants::missing::uintValue && node1 != node3)
        {
            mesh.ConnectNodes(node1, node3);
        }

        if (node2 != constants::missing::uintValue && node4 != constants::missing::uintValue && node2 != node4)
        {
            mesh.ConnectNodes(node2, node4);
        }
    }

    // create the diagonal links in quads that connect the new mesh with the old mesh
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
        if (!faceIsActive)
        {
            continue;
        }

        // Perhaps change quads to maximum number of edges for any shape
        std::array<UInt, constants::geometric::numNodesInQuadrilateral> oldIndex{constants::missing::uintValue, constants::missing::uintValue,
                                                                                 constants::missing::uintValue, constants::missing::uintValue};

        std::array<UInt, constants::geometric::numNodesInQuadrilateral> newIndex{constants::missing::uintValue, constants::missing::uintValue,
                                                                                 constants::missing::uintValue, constants::missing::uintValue};
        // find the old and new nodes
        for (UInt j = 0; j < mesh.m_numFacesNodes[i]; ++j)
        {
            const UInt previousIndex = (j == 0 ? (mesh.m_numFacesNodes[i] - 1) : (j - 1));

            const UInt edgeId = mesh.m_facesEdges[i][j];
            const UInt previousEdgeId = mesh.m_facesEdges[i][previousIndex];

            oldIndex[j] = mesh.m_edges[edgeId].first;
            newIndex[j] = newNodes[edgeId][2];

            if (oldIndex[j] != mesh.m_edges[previousEdgeId].first && oldIndex[j] != mesh.m_edges[previousEdgeId].second)
            {
                oldIndex[j] = mesh.m_edges[edgeId].second;
                newIndex[j] = newNodes[edgeId][1];
            }
        }

        for (UInt j = 0; j < mesh.m_numFacesNodes[i]; ++j)
        {
            const UInt previousIndex = (j == 0 ? (mesh.m_numFacesNodes[i] - 1) : (j - 1));
            const UInt nextIndex = (j == mesh.m_numFacesNodes[i] - 1 ? 0 : (j + 1));
            UInt nextNextIndex;

            if (j == mesh.m_numFacesNodes[i] - 2)
            {
                nextNextIndex = 0;
            }
            else if (j == mesh.m_numFacesNodes[i] - 1)
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

            // only one new node: new diagonal link connects new node with one old node
            if (nodeMask[node1] < NodeMask::Unassigned && nodeMask[node2] == NodeMask::Unassigned && nodeMask[node3] == NodeMask::Unassigned && nodeMask[node4] == NodeMask::Unassigned)
            {
                mesh.ConnectNodes(node1, node4);
                break;
            }

            // only one old node: new diagonal link connects new nodes only (i.e. perpendicular to previous one)
            if (nodeMask[node1] < NodeMask::Unassigned && nodeMask[node2] > NodeMask::Unassigned && nodeMask[node3] > NodeMask::Unassigned && nodeMask[node4] == NodeMask::Unassigned)
            {
                mesh.ConnectNodes(newIndex[previousIndex], newIndex[nextIndex]);
                break;
            }

            // two new and opposing nodes: new diagonal link connects the new nodes
            if (nodeMask[node1] < NodeMask::Unassigned && nodeMask[node2] == NodeMask::Unassigned && nodeMask[node3] == NodeMask::Unassigned && nodeMask[node4] == NodeMask::RegisteredNode)
            {
                mesh.ConnectNodes(node1, newIndex[nextNextIndex]);
                break;
            }
        }
    }

    std::vector<UInt> newEdges(InitialEdgeArraySize);

    // make the missing boundary links
    for (UInt i = 0; i < numNodes; ++i)
    {
        if (nodeMask[i] < NodeMask::BoundaryNode)
        {
            // boundary and kept nodes only
            continue;
        }

        UInt edgeCount = 0;
        std::fill(newEdges.begin(), newEdges.end(), constants::missing::uintValue);

        for (UInt j = 0; j < mesh.m_nodesNumEdges[i]; ++j)
        {
            UInt edgeId = mesh.m_nodesEdges[i][j];

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

                if (mesh.m_edges[edgeId].first == i)
                {
                    mesh.ConnectNodes(i, newNodes[edgeId][0]);
                    mesh.ConnectNodes(i, newNodes[edgeId][2]);
                }
                else
                {
                    mesh.ConnectNodes(i, newNodes[edgeId][1]);
                    mesh.ConnectNodes(i, newNodes[edgeId][3]);
                }
            }
        }

        if (edgeCount == 0)
        {
            continue;
        }

        std::vector<UInt> node(edgeCount, constants::missing::uintValue);

        for (UInt j = 0; j < edgeCount; ++j)
        {
            if (mesh.m_edges[newEdges[j]].first == i && nodeMask[newNodes[newEdges[j]][0]] == NodeMask::NewGeneralNode)
            {
                node[j] = newNodes[newEdges[j]][0];
            }

            if (mesh.m_edges[newEdges[j]].first == i && nodeMask[newNodes[newEdges[j]][2]] == NodeMask::NewGeneralNode)
            {
                node[j] = newNodes[newEdges[j]][2];
            }

            if (mesh.m_edges[newEdges[j]].second == i && nodeMask[newNodes[newEdges[j]][1]] == NodeMask::NewGeneralNode)
            {
                node[j] = newNodes[newEdges[j]][1];
            }

            if (mesh.m_edges[newEdges[j]].second == i && nodeMask[newNodes[newEdges[j]][3]] == NodeMask::NewGeneralNode)
            {
                node[j] = newNodes[newEdges[j]][3];
            }
        }

        if (nodeMask[i] != NodeMask::CornerNode)
        {
            if (edgeCount != 2)
            {
                throw AlgorithmError("Incorrect number of edges found: {}", edgeCount);
            }
            else
            {
                if (node[0] != constants::missing::uintValue && node[1] != constants::missing::uintValue && node[0] != node[1])
                {
                    mesh.ConnectNodes(node[0], node[1]);
                }
            }
        }
        else
        {
            for (UInt j = 0; j < edgeCount; ++j)
            {
                if (node[j] != constants::missing::uintValue && node[j] != i)
                {
                    mesh.ConnectNodes(i, node[j]);
                }
            }
        }
    }

    for (UInt i = 0; i < numNodes; ++i)
    {
        if (nodeMask[i] != NodeMask::CornerNode || mesh.m_nodesNumEdges[i] <= MaximumNumberOfNodesInNewlyCreatedElements)
        {
            continue;
        }

        for (UInt j = 0; j < mesh.m_nodesNumEdges[i]; ++j)
        {
            const UInt edgeId = mesh.m_nodesEdges[i][j];

            if (mesh.m_edges[edgeId].first == i)
            {
                mesh.ConnectNodes(i, newNodes[edgeId][0]);
                mesh.ConnectNodes(i, newNodes[edgeId][2]);
            }
            else
            {
                mesh.ConnectNodes(i, newNodes[edgeId][1]);
                mesh.ConnectNodes(i, newNodes[edgeId][3]);
            }
        }
    }
}

void meshkernel::CasulliRefinement::ComputeNewNodes(Mesh2D& mesh, std::vector<EdgeNodes>& newNodes, std::vector<NodeMask>& nodeMask)
{
    // Keep copy of number of edges in mesh before any nodes are added.
    const UInt numEdges = mesh.GetNumEdges();

    for (UInt i = 0; i < mesh.GetNumFaces(); ++i)
    {
        const Point elementCentre = mesh.m_facesCircumcenters[i];

        for (UInt j = 0; j < mesh.m_numFacesNodes[i]; ++j)
        {
            const UInt elementNode = mesh.m_facesNodes[i][j];

            UInt firstEdgeId = constants::missing::uintValue;
            UInt secondEdgeId = constants::missing::uintValue;
            UInt newNodeId = constants::missing::uintValue;

            for (UInt k = 0; k < mesh.m_facesEdges[i].size(); ++k)
            {
                UInt edgeId = mesh.m_facesEdges[i][k];

                if (mesh.m_edges[edgeId].first == elementNode || mesh.m_edges[edgeId].second == elementNode)
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
                Point newNode = 0.5 * (elementCentre + mesh.m_nodes[elementNode]);

                newNodeId = mesh.InsertNode(newNode);
                nodeMask[newNodeId] = NodeMask::NewAssignedNode;
            }
            else
            {
                newNodeId = elementNode;
            }

            StoreNewNode(mesh, elementNode, firstEdgeId, secondEdgeId, newNodeId, newNodes);
        }
    }

    for (UInt i = 0; i < numEdges; ++i)
    {
        UInt newNodeId = constants::missing::uintValue;

        if (mesh.m_edgesNumFaces[i] != 1)
        {
            continue;
        }

        const UInt node1 = mesh.m_edges[i].first;
        const UInt node2 = mesh.m_edges[i].second;

        if (node1 == constants::missing::uintValue && node2 == constants::missing::uintValue)
        {
            continue;
        }

        const Point edgeCentre = 0.5 * (mesh.m_nodes[node1] + mesh.m_nodes[node2]);

        if (nodeMask[node1] != NodeMask::Unassigned)
        {
            const Point newNode = 0.5 * (edgeCentre + mesh.m_nodes[node1]);

            newNodeId = mesh.InsertNode(newNode);
            nodeMask[newNodeId] = NodeMask::NewGeneralNode;
        }
        else
        {
            newNodeId = node1;
        }

        StoreNewNode(mesh, node1, i, i, newNodeId, newNodes);

        if (nodeMask[node2] != NodeMask::Unassigned)
        {
            const Point newNode = 0.5 * (edgeCentre + mesh.m_nodes[node2]);

            newNodeId = mesh.InsertNode(newNode);
            nodeMask[newNodeId] = NodeMask::NewGeneralNode;
        }
        else
        {
            newNodeId = node2;
        }

        StoreNewNode(mesh, node2, i, i, newNodeId, newNodes);
    }
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
