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
    std::vector<LinkNodes> newNodes(mesh.GetNumEdges(), {constants::missing::uintValue, constants::missing::uintValue, constants::missing::uintValue, constants::missing::uintValue});
    std::vector<int> nodeMask(InitialiseNodeMask(mesh, polygon));

    UInt numNodes = mesh.GetNumNodes();
    UInt numEdges = mesh.GetNumEdges();
    UInt numFaces = mesh.GetNumFaces();

    ComputeNewNodes(mesh, newNodes, nodeMask);
    LinkNewNodes(mesh, newNodes, numNodes, numEdges, numFaces, nodeMask);
    Administrate(mesh, numNodes, nodeMask);
}

std::vector<int> meshkernel::CasulliRefinement::InitialiseNodeMask(const Mesh2D& mesh, const Polygons& polygon)
{
    // Need a better initial size
    std::vector<int> nodeMask(10 * mesh.GetNumNodes(), 0);

    // Find nodes that are inside the polygon.
    // If the polygon is empty then all nodes will be taken into account.
    for (UInt i = 0; i < mesh.GetNumNodes(); ++i)
    {
        auto [containsPoint, pointIndex] = polygon.IsPointInPolygons(mesh.m_nodes[i]);

        if (containsPoint)
        {
            nodeMask[i] = 1;
        }
    }

    // Find nodes that lie on the boundary of the domain.
    for (UInt i = 0; i < mesh.GetNumEdges(); ++i)
    {
        UInt node1 = mesh.m_edges[i].first;
        UInt node2 = mesh.m_edges[i].second;

        if (mesh.m_edgesNumFaces[i] == 1)
        {

            if (nodeMask[node1] != 0)
            {
                nodeMask[node1] = 1234;
            }

            if (nodeMask[node2] != 0)
            {
                nodeMask[node2] = 1234;
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
            while (((mesh.m_edges[edge2].first != i && mesh.m_edges[edge2].second) || edge2 == edge1) && faceEdgeIndex < nodeCount - 1)
            {
                ++faceEdgeIndex;
                edge2 = mesh.m_facesEdges[elementId][faceEdgeIndex];
            }

            if (mesh.m_edgesNumFaces[edge2] == 1)
            {
                if (nodeMask[i] > 0)
                {
                    nodeMask[i] = 1235;
                    break;
                }
            }
        }
    }

    // Find included corner nodes
    for (UInt i = 0; i < mesh.GetNumNodes(); ++i)
    {
        if (nodeMask[i] != 0 && mesh.m_nodesTypes[i] == 3)
        {
            nodeMask[i] = 1235;
        }
    }

    std::vector<UInt> sharedFaces;
    std::vector<UInt> connectedNodes;
    std::vector<std::vector<UInt>> faceNodeMapping;

    for (UInt i = 0; i < mesh.GetNumNodes(); ++i)
    {
        if (nodeMask[i] == 0)
        {
            continue;
        }

        if (mesh.m_nodesNumEdges[i] > 1)
        {
            FindPatchIds(mesh, i, sharedFaces, connectedNodes, faceNodeMapping);
        }
        else
        {
            nodeMask[i] = 0;
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
            nodeMask[i] = 0;
        }

        if (elementCount < mesh.m_nodesNumEdges[i] - 1 && nodeMask[i] == 1234)
        {
            nodeMask[i] = 1235;
        }

        if (elementCount > Mesh::m_maximumNumberOfEdgesPerNode && nodeMask[i] > 0 && nodeMask[i] < 1234)
        {
            nodeMask[i] = 1235;
        }
    }

    return nodeMask;
}

void meshkernel::CasulliRefinement::Administrate(Mesh2D& mesh, const UInt numNodes, const std::vector<int>& nodeMask)
{
    // Need check only the original nodes in the mesh, hence use of numNodes.
    for (UInt i = 0; i < numNodes; ++i)
    {
        if (nodeMask[i] > 0 && nodeMask[i] < 1235)
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

void meshkernel::CasulliRefinement::LinkNewNodes(Mesh2D& mesh, const std::vector<LinkNodes>& newNodes, const UInt numNodes, const UInt numEdges, const UInt numFaces, std::vector<int>& nodeMask)
{
    //  make the original-link based new links
    for (UInt i = 0; i < numEdges; ++i)
    {
        UInt node1 = newNodes[i][0];
        UInt node2 = newNodes[i][1];
        UInt node3 = newNodes[i][2];
        UInt node4 = newNodes[i][3];

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
            if (nodeMask[mesh.m_facesNodes[i][j]] == 0)
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
        // Could replace mesh.m_numFacesNodes[i] with 4
        for (UInt j = 0; j < mesh.m_numFacesNodes[i]; ++j)
        {
            UInt previousIndex = (j == 0 ? (mesh.m_numFacesNodes[i] - 1) : (j - 1));

            UInt edgeId = mesh.m_facesEdges[i][j];
            UInt previousEdgeId = mesh.m_facesEdges[i][previousIndex];

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
            UInt previousIndex = (j == 0 ? (mesh.m_numFacesNodes[i] - 1) : (j - 1));
            UInt nextIndex = (j == mesh.m_numFacesNodes[i] - 1 ? 0 : (j + 1));
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

            UInt node1 = newIndex[j];
            UInt node2 = oldIndex[previousIndex];
            UInt node3 = oldIndex[nextIndex];
            UInt node4 = oldIndex[nextNextIndex];

            // only one new node: new diagonal link connects new node with one old node
            if (nodeMask[node1] < 0 && nodeMask[node2] == 0 && nodeMask[node3] == 0 && nodeMask[node4] == 0)
            {
                mesh.ConnectNodes(node1, node4);
                break;
            }

            // only one old node: new diagonal link connects new nodes only (i.e. perpendicular to previous one)
            if (nodeMask[node1] < 0 && nodeMask[node2] > 0 && nodeMask[node3] > 0 && nodeMask[node4] == 0)
            {
                mesh.ConnectNodes(newIndex[previousIndex], newIndex[nextIndex]);
                break;
            }

            // two new and opposing nodes: new diagonal link connects the new nodes
            if (nodeMask[node1] < 0 && nodeMask[node2] == 0 && nodeMask[node3] == 0 && nodeMask[node4] == 1)
            {
                mesh.ConnectNodes(node1, newIndex[nextNextIndex]);
                break;
            }
        }
    }

    std::vector<UInt> link(InitialEdgeArraySize);

    // make the missing boundary links
    for (UInt i = 0; i < numNodes; ++i)
    {
        if (nodeMask[i] < 1234)
        {
            // boundary and kept nodes only
            continue;
        }

        UInt edgeCount = 0;
        std::fill(link.begin(), link.end(), constants::missing::uintValue);

        for (UInt j = 0; j < mesh.m_nodesNumEdges[i]; ++j)
        {
            UInt edgeId = mesh.m_nodesEdges[i][j];

            if (mesh.m_edgesNumFaces[edgeId] == 1)
            {
                if (edgeCount >= link.size())
                {
                    link.resize(2 * edgeCount + 1);
                }

                link[edgeCount] = edgeId;
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
            if (mesh.m_edges[link[j]].first == i && nodeMask[newNodes[link[j]][0]] == -1)
            {
                node[j] = newNodes[link[j]][0];
            }

            if (mesh.m_edges[link[j]].first == i && nodeMask[newNodes[link[j]][2]] == -1)
            {
                node[j] = newNodes[link[j]][2];
            }

            if (mesh.m_edges[link[j]].second == i && nodeMask[newNodes[link[j]][1]] == -1)
            {
                node[j] = newNodes[link[j]][1];
            }

            if (mesh.m_edges[link[j]].second == i && nodeMask[newNodes[link[j]][3]] == -1)
            {
                node[j] = newNodes[link[j]][3];
            }

            //

            if (mesh.m_edges[link[j]].first == i && nodeMask[newNodes[link[j]][0]] == 1236)
            {
                node[j] = newNodes[link[j]][0];
            }

            if (mesh.m_edges[link[j]].first == i && nodeMask[newNodes[link[j]][2]] == 1236)
            {
                node[j] = newNodes[link[j]][2];
            }

            if (mesh.m_edges[link[j]].second == i && nodeMask[newNodes[link[j]][1]] == 1236)
            {
                node[j] = newNodes[link[j]][1];
            }

            if (mesh.m_edges[link[j]].second == i && nodeMask[newNodes[link[j]][3]] == 1236)
            {
                node[j] = newNodes[link[j]][3];
            }
        }

        if (nodeMask[i] != 1235 && nodeMask[i] != 1236)
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
        if (nodeMask[i] != 1235 || mesh.m_nodesNumEdges[i] <= MaximumNumberOfNodesInNewlyCreatedElements)
        {
            continue;
        }

        for (UInt j = 0; j < mesh.m_nodesNumEdges[i]; ++j)
        {
            UInt edgeId = mesh.m_nodesEdges[i][j];

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

void meshkernel::CasulliRefinement::ComputeNewNodes(Mesh2D& mesh, std::vector<LinkNodes>& newNodes, std::vector<int>& nodeMask)
{
    // Keep copy of number of edges in mesh before any nodes are added.
    UInt numEdges = mesh.GetNumEdges();

    for (UInt i = 0; i < mesh.GetNumFaces(); ++i)
    {
        Point elementCentre = mesh.m_facesCircumcenters[i];

        for (UInt j = 0; j < mesh.m_numFacesNodes[i]; ++j)
        {
            UInt elementNode = mesh.m_facesNodes[i][j];

            UInt firstEdgeId = constants::missing::uintValue;
            UInt secondEdgeId = constants::missing::uintValue;
            UInt newPointIndex = constants::missing::uintValue;

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

            if (nodeMask[elementNode] > 0)
            {
                Point newNode = 0.5 * (elementCentre + mesh.m_nodes[elementNode]);

                newPointIndex = mesh.InsertNode(newNode);
                nodeMask[newPointIndex] = -2;
            }
            else
            {
                newPointIndex = elementNode;
            }

            StoreNewNode(mesh, elementNode, firstEdgeId, secondEdgeId, newPointIndex, newNodes);
        }
    }

    for (UInt i = 0; i < numEdges; ++i)
    {
        UInt newPointIndex = constants::missing::uintValue;

        if (mesh.m_edgesNumFaces[i] != 1)
        {
            continue;
        }

        UInt node1 = mesh.m_edges[i].first;
        UInt node2 = mesh.m_edges[i].second;

        if (node1 == constants::missing::uintValue && node2 == constants::missing::uintValue)
        {
            continue;
        }

        Point edgeCentre = 0.5 * (mesh.m_nodes[node1] + mesh.m_nodes[node2]);

        if (nodeMask[node1] != 0)
        {
            Point newNode = 0.5 * (edgeCentre + mesh.m_nodes[node1]);

            newPointIndex = mesh.InsertNode(newNode);
            nodeMask[newPointIndex] = -1;
        }
        else
        {
            newPointIndex = node1;
        }

        StoreNewNode(mesh, node1, i, i, newPointIndex, newNodes);

        if (nodeMask[node2] != 0)
        {
            Point newNode = 0.5 * (edgeCentre + mesh.m_nodes[node2]);

            newPointIndex = mesh.InsertNode(newNode);
            nodeMask[newPointIndex] = -1;
        }
        else
        {
            newPointIndex = node2;
        }

        StoreNewNode(mesh, node2, i, i, newPointIndex, newNodes);
    }
}

void meshkernel::CasulliRefinement::StoreNewNode(const Mesh2D& mesh, const UInt nodeId, const UInt link1Index, const UInt link2Index, const UInt newPointIndex, std::vector<LinkNodes>& newNodes)
{
    UInt edgeId1 = link1Index;
    UInt edgeId2 = link2Index;

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

    UInt elementId = FindCommon(mesh, edgeId1, edgeId2);

    if (elementId == constants::missing::uintValue)
    {
        throw AlgorithmError("No element found that shares edge: {} and {}", edgeId1, edgeId2);
    }

    UInt lr1 = IsLeftRight(mesh, elementId, edgeId1);
    UInt lr2 = IsLeftRight(mesh, elementId, edgeId2);

    UInt se1 = IsStartEnd(mesh, nodeId, edgeId1);
    UInt se2 = IsStartEnd(mesh, nodeId, edgeId2);

    if (edgeId1 == edgeId2)
    {
        lr1 = 1 - lr1;
        lr2 = 1 - lr2;
    }

    UInt iPoint1 = se1 + 2 * (1 - lr1);
    UInt iPoint2 = se2 + 2 * (1 - lr2);

    if (newNodes[edgeId1][iPoint1] == constants::missing::uintValue)
    {
        newNodes[edgeId1][iPoint1] = newPointIndex;
    }

    if (newNodes[edgeId2][iPoint2] == constants::missing::uintValue)
    {
        newNodes[edgeId2][iPoint2] = newPointIndex;
    }
}

meshkernel::UInt meshkernel::CasulliRefinement::IsStartEnd(const Mesh2D& mesh, const UInt nodeId, const UInt edgeId)
{
    UInt isStartEnd = constants::missing::uintValue;

    if (mesh.m_edges[edgeId].first == nodeId)
    {
        isStartEnd = 0;
    }
    else if (mesh.m_edges[edgeId].second == nodeId)
    {
        isStartEnd = 1;
    }

    return isStartEnd;
}

meshkernel::UInt meshkernel::CasulliRefinement::IsLeftRight(const Mesh2D& mesh, const UInt elementId, const UInt edgeId)
{
    UInt isLeftRight = constants::missing::uintValue;
    UInt edgeIndex = constants::missing::uintValue;
    UInt nextEdgeIndex = constants::missing::uintValue;
    UInt endNodeIndex = mesh.m_edges[edgeId].second;

    for (UInt i = 0; i < mesh.m_facesEdges[elementId].size(); ++i)
    {
        UInt faceEdgeId = mesh.m_facesEdges[elementId][i];

        if (faceEdgeId == edgeId)
        {
            edgeIndex = i;
        }
        else if (mesh.m_edges[faceEdgeId].first == endNodeIndex || mesh.m_edges[faceEdgeId].second == endNodeIndex)
        {
            nextEdgeIndex = i;
        }
    }

    if (edgeIndex == constants::missing::uintValue || nextEdgeIndex == constants::missing::uintValue)
    {
        // EdgeId was not found
        return isLeftRight;
    }

    if (nextEdgeIndex == edgeIndex + 1 || nextEdgeIndex + mesh.m_numFacesNodes[elementId] == edgeIndex + 1)
    {
        isLeftRight = 0;
    }
    else if (edgeIndex == nextEdgeIndex + 1 || edgeIndex + mesh.m_numFacesNodes[elementId] == nextEdgeIndex + 1)
    {
        isLeftRight = 1;
    }

    return isLeftRight;
}

meshkernel::UInt meshkernel::CasulliRefinement::FindCommon(const Mesh2D& mesh, const UInt edge1, const UInt edge2)
{
    UInt commonElement = constants::missing::uintValue;

    for (UInt i = 0; i < mesh.m_edgesNumFaces[edge1]; ++i)
    {
        for (UInt j = 0; j < mesh.m_edgesNumFaces[edge2]; ++j)
        {
            if (mesh.m_edgesFaces[edge1][i] == mesh.m_edgesFaces[edge2][j])
            {
                commonElement = mesh.m_edgesFaces[edge1][i];
                break;
            }
        }

        if (commonElement != constants::missing::uintValue)
        {
            break;
        }
    }

    return commonElement;
}
