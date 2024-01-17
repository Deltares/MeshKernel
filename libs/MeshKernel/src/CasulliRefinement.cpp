#include <algorithm>
#include <iostream>
#include <vector>

#include "MeshKernel/CasulliRefinement.hpp"

void meshkernel::CasulliRefinement::Compute(Mesh2D& mesh, const MeshRefinementParameters& meshRefinementParameters [[maybe_unused]])
{
    std::vector<LinkNodes> newNodes(mesh.GetNumEdges(), {constants::missing::uintValue, constants::missing::uintValue, constants::missing::uintValue, constants::missing::uintValue});
    // Need a better initial size
    std::vector<int> kc(10 * mesh.GetNumNodes(), 1);

    UInt numNodes = mesh.GetNumNodes();
    UInt numEdges = mesh.GetNumEdges();
    UInt numFaces = mesh.GetNumFaces();

    Point clickedPoint(0.0, 0.0);

    Initialise(mesh, kc);

    ComputeNewNodes(mesh, newNodes, kc);
    ComputeNewNodesDirectional(mesh, newNodes, kc, clickedPoint);
    // makelinks
    LinkNewNodes(mesh, newNodes, numNodes, numEdges, numFaces, kc);

    // Need check only the original nodes in the mesh, hence use of numNodes.
    for (UInt i = 0; i < numNodes; ++i)
    {
        if (kc[i] > 0 && kc[i] < 1235)
        {
            mesh.DeleteNode(i);
        }
    }

#if 0
    for (UInt L = 0; L < mesh.GetNumEdges(); ++L)
    // for (UInt L = 0; L < numEdges; ++L)
    {
        UInt node1 = mesh.m_edges[L].first;
        UInt node2 = mesh.m_edges[L].second;

        if (node1 == constants::missing::uintValue || node2 == constants::missing::uintValue)
        {
            continue;
        }

        if ((kc[node1] > 0 && kc[node2] > 0) ||
            (kc[node1] == 1235 && kc[node2] == 1235) ||
            (kc[node1] == 0 && kc[node2] == 1235) ||
            (kc[node1] == 1235 && kc[node2] == 0) ||
            (kc[node1] == 0 && kc[node2] == 1236) ||
            (kc[node1] == 1236 && kc[node2] == 0) ||
            (kc[node1] == -1 && kc[node2] == -1) ||
            (kc[node1] == -2 && kc[node2] == -2))
        {

            if (mesh.m_edgesNumFaces[L] == 0) // || mesh.m_edges[L])
            {
                continue;
            }

            mesh.DeleteEdge(L);
        }
    }
#endif

    mesh.Administrate();
}

void meshkernel::CasulliRefinement::OrthonetAdmin(const Mesh2D& mesh, const UInt currentNode, std::vector<UInt>& sharedFaces, std::vector<UInt>& connectedNodes, std::vector<std::vector<UInt>>& faceNodeMapping)
{
    sharedFaces.clear();
    connectedNodes.clear();
    faceNodeMapping.clear();

    if (currentNode >= mesh.GetNumNodes())
    {
        throw MeshKernelError("Node index out of range");
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

void meshkernel::CasulliRefinement::Initialise(const Mesh2D& mesh [[maybe_unused]], std::vector<int>& kc)
{
    std::fill(kc.begin(), kc.end(), 1);

    for (UInt k = 0; k < mesh.GetNumEdges(); ++k)
    {
        UInt node1 = mesh.m_edges[k].first;
        UInt node2 = mesh.m_edges[k].second;

        if (mesh.m_edgesNumFaces[k] == 1)
        {

            if (kc[node1] != 0)
            {
                kc[node1] = 1234;
            }

            if (kc[node2] != 0)
            {
                kc[node2] = 1234;
            }
        }
    }

    for (UInt k = 0; k < mesh.GetNumNodes(); ++k)
    {

        for (UInt kk = 0; k < mesh.m_nodesNumEdges[k]; ++k)
        {
            UInt link1 = mesh.m_nodesEdges[k][kk];

            if (mesh.m_edgesNumFaces[link1] != 1)
            {
                continue;
            }

            UInt icell = mesh.m_edgesFaces[link1][0];
            UInt N = mesh.m_numFacesNodes[icell];

            UInt kkk = 0;
            UInt link2 = mesh.m_facesEdges[icell][kkk];

            // Check the loop termination, especially the kkk < N - 1
            // Perhaps change to for loop checking the condition then break.
            while (((mesh.m_edges[link2].first != k && mesh.m_edges[link2].second) || link2 == link1) && kkk < N - 1)
            {
                ++kkk;
                link2 = mesh.m_facesEdges[icell][kkk];
            }

            if (mesh.m_edgesNumFaces[link2] == 1)
            {
                if (kc[k] > 0)
                {
                    kc[k] = 1235;
                    break;
                }
            }
        }
    }

    for (UInt k = 0; k < mesh.GetNumNodes(); ++k)
    {
        if (mesh.m_nodesTypes[k] == 3)
        {
            kc[k] = 1235;
        }
    }

    UInt nmkx = 0;
    std::vector<UInt> sharedFaces;
    std::vector<UInt> connectedNodes;
    std::vector<std::vector<UInt>> faceNodeMapping;

    for (UInt k0 = 0; k0 < mesh.GetNumNodes(); ++k0)
    {
        if (kc[k0] == 0)
        {
            continue;
        }

        if (mesh.m_nodesNumEdges[k0] > 1)
        {
            OrthonetAdmin(mesh, k0, sharedFaces, connectedNodes, faceNodeMapping);
        }
        else
        {
            kc[k0] = 0;
            continue;
        }

        UInt ncell = 0;

        for (UInt kk = 0; kk < sharedFaces.size(); ++kk)
        {

            if (sharedFaces[kk] != constants::missing::uintValue)
            {
                ++ncell;
            }
        }

        if (ncell == 0)
        {
            kc[k0] = 0;
        }

        if (ncell < mesh.m_nodesNumEdges[k0] - 1 && kc[k0] == 1234)
        {
            kc[k0] = 1235;
        }

        if (ncell > Mesh::m_maximumNumberOfEdgesPerNode && kc[k0] > 0 && kc[k0] < 1234)
        {
            kc[k0] = 1235;
        }

        nmkx = std::max(nmkx, mesh.m_nodesNumEdges[k0]);
    }
}

void meshkernel::CasulliRefinement::LinkNewNodes(Mesh2D& mesh, const std::vector<LinkNodes>& newNodes, const UInt numNodes, const UInt numEdges, const UInt numFaces, std::vector<int>& kc)
{

    //  make the original-link based new links
    for (UInt i = 0; i < numEdges; ++i)
    {
        UInt k1 = newNodes[i][0];
        UInt k2 = newNodes[i][1];
        UInt k3 = newNodes[i][2];
        UInt k4 = newNodes[i][3];

        // Parallel edges, these are the start-end connections

        if (k1 != constants::missing::uintValue && k2 != constants::missing::uintValue && k1 != k2)
        {
            mesh.ConnectNodes(k1, k2);
        }

        if (k3 != constants::missing::uintValue && k4 != constants::missing::uintValue && k3 != k4)
        {
            mesh.ConnectNodes(k3, k4);
        }

        // normal edges, these are the left-right connections

        if (k1 != constants::missing::uintValue && k3 != constants::missing::uintValue && k1 != k3)
        {
            mesh.ConnectNodes(k1, k3);
        }

        if (k2 != constants::missing::uintValue && k4 != constants::missing::uintValue && k2 != k4)
        {
            mesh.ConnectNodes(k2, k4);
        }
    }

    // create the diagonal links in quads that connect the new mesh with the old mesh
    for (UInt i = 0; i < numFaces; ++i)
    {

        if (mesh.m_numFacesNodes[i] != 4)
        {
            continue;
        }

        // Check for active nodes
        // if (!IsActive (mesh, kc, i) {
        //     continue;
        // }

        // Perhaps change quads to maximum number of edges for any shape
        std::array<UInt, constants::geometric::numNodesInQuadrilateral> oldn{constants::missing::uintValue, constants::missing::uintValue, constants::missing::uintValue, constants::missing::uintValue};
        std::array<UInt, constants::geometric::numNodesInQuadrilateral> newn{constants::missing::uintValue, constants::missing::uintValue, constants::missing::uintValue, constants::missing::uintValue};
        // std::array<UInt, constants::geometric::numNodesInQuadrilateral> oldn{0, 0, 0, 0};
        // std::array<UInt, constants::geometric::numNodesInQuadrilateral> newn{0, 0, 0, 0};

        // find the old and new nodes
        // Could replace mesh.m_numFacesNodes[i] with 4
        for (UInt j = 0; j < mesh.m_numFacesNodes[i]; ++j)
        {
            // rename variable kkm1
            UInt kkm1 = (j == 0 ? (mesh.m_numFacesNodes[i] - 1) : (j - 1));

            UInt L = mesh.m_facesEdges[i][j];
            UInt Lm1 = mesh.m_facesEdges[i][kkm1];

            oldn[j] = mesh.m_edges[L].first;
            newn[j] = newNodes[L][2];

            if (oldn[j] != mesh.m_edges[Lm1].first && oldn[j] != mesh.m_edges[Lm1].second)
            {
                oldn[j] = mesh.m_edges[L].second;
                newn[j] = newNodes[L][1];
            }
        }

        for (UInt j = 0; j < mesh.m_numFacesNodes[i]; ++j)
        {
            UInt kkm1 = (j == 0 ? (mesh.m_numFacesNodes[i] - 1) : (j - 1));
            UInt kkp1 = (j == mesh.m_numFacesNodes[i] - 1 ? 0 : (j + 1));
            UInt kkp2;

            if (j == mesh.m_numFacesNodes[i] - 2)
            {
                kkp2 = 0;
            }
            else if (j == mesh.m_numFacesNodes[i] - 1)
            {
                kkp2 = 1;
            }
            else
            {
                kkp2 = j + 2;
            }

            UInt k1 = newn[j];
            UInt k2 = oldn[kkm1];
            UInt k3 = oldn[kkp1];
            UInt k4 = oldn[kkp2];

            // only one new node: new diagonal link connects new node with one old node
            if (kc[k1] < 0 && kc[k2] == 0 && kc[k3] == 0 && kc[k4] == 0)
            {
                mesh.ConnectNodes(k1, k4);
                break;
            }

            // only one old node: new diagonal link connects new nodes only (i.e. perpendicular to previous one)
            if (kc[k1] < 0 && kc[k2] > 0 && kc[k3] > 0 && kc[k4] == 0)
            {
                mesh.ConnectNodes(newn[kkm1], newn[kkp1]);
                break;
            }

            // two new and opposing nodes: new diagonal link connects the new nodes
            if (kc[k1] < 0 && kc[k2] == 0 && kc[k3] == 0 && kc[k4] == 1)
            {
                mesh.ConnectNodes(k1, newn[kkp2]);
                break;
            }
        }
    }

    // Better number than 100
    std::vector<UInt> link(100);

    // make the missing boundary links

    for (UInt k = 0; k < numNodes; ++k)
    {
        if (kc[k] < 1234)
        {
            // boundary and kept nodes only
            continue;
        }

        UInt numLinks = 0;
        std::fill(link.begin(), link.end(), constants::missing::uintValue);

        for (UInt j = 0; j < mesh.m_nodesNumEdges[k]; ++j)
        {
            UInt L = mesh.m_nodesEdges[k][j];

            if (mesh.m_edgesNumFaces[L] == 1)
            {
                link[numLinks] = L;
                ++numLinks;
            }
            else
            {

                if (mesh.m_edges[L].first == k)
                {
                    mesh.ConnectNodes(k, newNodes[L][0]);
                    mesh.ConnectNodes(k, newNodes[L][2]);
                }
                else
                {
                    mesh.ConnectNodes(k, newNodes[L][1]);
                    mesh.ConnectNodes(k, newNodes[L][3]);
                }
            }
        }

        if (numLinks == 0)
        {
            continue;
        }

        std::vector<UInt> node(numLinks, constants::missing::uintValue);

        for (UInt kk = 0; kk < numLinks; ++kk)
        {
            if (mesh.m_edges[link[kk]].first == k && kc[newNodes[link[kk]][0]] == -1)
            {
                node[kk] = newNodes[link[kk]][0];
            }

            if (mesh.m_edges[link[kk]].first == k && kc[newNodes[link[kk]][2]] == -1)
            {
                node[kk] = newNodes[link[kk]][2];
            }

            if (mesh.m_edges[link[kk]].second == k && kc[newNodes[link[kk]][1]] == -1)
            {
                node[kk] = newNodes[link[kk]][1];
            }

            if (mesh.m_edges[link[kk]].second == k && kc[newNodes[link[kk]][3]] == -1)
            {
                node[kk] = newNodes[link[kk]][3];
            }

            //

            if (mesh.m_edges[link[kk]].first == k && kc[newNodes[link[kk]][0]] == 1236)
            {
                node[kk] = newNodes[link[kk]][0];
            }

            if (mesh.m_edges[link[kk]].first == k && kc[newNodes[link[kk]][2]] == 1236)
            {
                node[kk] = newNodes[link[kk]][2];
            }

            if (mesh.m_edges[link[kk]].second == k && kc[newNodes[link[kk]][1]] == 1236)
            {
                node[kk] = newNodes[link[kk]][1];
            }

            if (mesh.m_edges[link[kk]].second == k && kc[newNodes[link[kk]][3]] == 1236)
            {
                node[kk] = newNodes[link[kk]][3];
            }
        }

        if (kc[k] != 1235 && kc[k] != 1236)
        {
            if (numLinks != 2)
            {
                // Error
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
            for (UInt kk = 0; kk < numLinks; ++kk)
            {
                if (node[kk] != constants::missing::uintValue && node[kk] != k)
                {
                    mesh.ConnectNodes(k, node[kk]);
                }
            }
        }
    }

    for (UInt k = 0; k < numNodes; ++k)
    {

        // if (kc[k] != 1235 || mesh.m_nodesNumEdges[k] <= constants::geometric::numNodesInQuadrilateral)
        if (kc[k] != 1235 || mesh.m_nodesNumEdges[k] < 5)
        {
            continue;
        }

        for (UInt kk = 0; kk < mesh.m_nodesNumEdges[k]; ++kk)
        {
            UInt L = mesh.m_nodesEdges[k][kk];

            if (mesh.m_edges[L].first == k)
            {
                mesh.ConnectNodes(k, newNodes[L][0]);
                mesh.ConnectNodes(k, newNodes[L][2]);
            }
            else
            {
                mesh.ConnectNodes(k, newNodes[L][1]);
                mesh.ConnectNodes(k, newNodes[L][3]);
            }
        }
    }
}

void meshkernel::CasulliRefinement::ComputeNewNodes(Mesh2D& mesh, std::vector<LinkNodes>& newNodes, std::vector<int>& kc)
{
    [[maybe_unused]] UInt numEdges = mesh.GetNumEdges();

    for (UInt i = 0; i < mesh.GetNumFaces(); ++i)
    {
        Point elementCentre = mesh.m_facesCircumcenters[i];

        for (UInt j = 0; j < mesh.m_numFacesNodes[i]; ++j)
        {

            UInt knode = mesh.m_facesNodes[i][j];

            UInt link1 = constants::missing::uintValue;
            UInt link2 = constants::missing::uintValue;
            UInt knew = constants::missing::uintValue;

            for (UInt k = 0; k < mesh.m_facesEdges[i].size(); ++k)
            {
                UInt L = mesh.m_facesEdges[i][k];

                if (mesh.m_edges[L].first == knode || mesh.m_edges[L].second == knode)
                {
                    if (link1 == constants::missing::uintValue)
                    {
                        link1 = L;
                    }
                    else
                    {
                        link2 = L;
                        break;
                    }
                }
            }

            if (link1 == constants::missing::uintValue || link2 == constants::missing::uintValue)
            {
                // No links found
                continue;
            }

            if (kc[knode] > 0)
            {
                Point newNode = 0.5 * (elementCentre + mesh.m_nodes[knode]);
                knew = mesh.InsertNode(newNode);
                kc[knew] = -2;
            }
            else
            {
                knew = knode;
            }

            StoreNewNode(mesh, knode, link1, link2, knew, newNodes);
        }
    }

    for (UInt i = 0; i < numEdges; ++i)
    {
        UInt knew = constants::missing::uintValue;

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

        if (kc[node1] != 0)
        {
            Point newNode = 0.5 * (edgeCentre + mesh.m_nodes[node1]);
            knew = mesh.InsertNode(newNode);
            kc[knew] = -1;
        }
        else
        {
            knew = node1;
        }

        StoreNewNode(mesh, node1, i, i, knew, newNodes);

        if (kc[node2] != 0)
        {
            Point newNode = 0.5 * (edgeCentre + mesh.m_nodes[node2]);
            knew = mesh.InsertNode(newNode);
            kc[knew] = -1;
        }
        else
        {
            knew = node2;
        }

        StoreNewNode(mesh, node2, i, i, knew, newNodes);
    }

    dummy = 1;
}

void meshkernel::CasulliRefinement::StoreNewNode(const Mesh2D& mesh, const UInt knode, const UInt link1, const UInt link2, const UInt knew, std::vector<LinkNodes>& newNodes)
{
    UInt link1Copy = link1;
    UInt link2Copy = link2;

    if (link1Copy != constants::missing::uintValue)
    {
        if (link2Copy == constants::missing::uintValue)
        {
            link2Copy = link1Copy;
        }
    }
    else
    {
        if (link2Copy != constants::missing::uintValue)
        {
            link1Copy = link2Copy;
        }
        else
        {
            // Log error/warning message "no links specified"
            return;
        }
    }

    UInt element = FindCommon(mesh, link1Copy, link2Copy);

    if (element == constants::missing::uintValue)
    {
        // error no elements found
    }

    UInt lr1 = IsLeftRight(mesh, element, link1Copy);
    UInt lr2 = IsLeftRight(mesh, element, link2Copy);

    UInt se1 = IsStartEnd(mesh, knode, link1Copy);
    UInt se2 = IsStartEnd(mesh, knode, link2Copy);

    if (link1Copy == link2Copy)
    {
        lr1 = 1 - lr1;
        lr2 = 1 - lr2;
    }

    UInt iPoint1 = 0 + se1 + 2 * (1 - lr1);
    UInt iPoint2 = 0 + se2 + 2 * (1 - lr2);

    if (newNodes[link1Copy][iPoint1] == constants::missing::uintValue)
    {
        newNodes[link1Copy][iPoint1] = knew;
    }

    if (newNodes[link2Copy][iPoint2] == constants::missing::uintValue)
    {
        newNodes[link2Copy][iPoint2] = knew;
    }
}

meshkernel::UInt meshkernel::CasulliRefinement::IsStartEnd(const Mesh2D& mesh, const UInt nodeId, const UInt edgeId)
{
    UInt isStartEnd = constants::missing::uintValue;

    if (mesh.m_edges[edgeId].first == nodeId)
    {
        isStartEnd = 0;
    }

    if (mesh.m_edges[edgeId].second == nodeId)
    {
        isStartEnd = 1;
    }

    return isStartEnd;
}

meshkernel::UInt meshkernel::CasulliRefinement::IsLeftRight(const Mesh2D& mesh, const UInt elementId, const UInt edgeId)
{
    UInt isLeftRight = constants::missing::uintValue;
    UInt kself = constants::missing::uintValue;
    UInt knext = constants::missing::uintValue;
    UInt kend = mesh.m_edges[edgeId].second;

    for (UInt i = 0; i < mesh.m_facesEdges[elementId].size(); ++i)
    {
        UInt L1 = mesh.m_facesEdges[elementId][i];

        if (L1 == edgeId)
        {
            kself = i;
        }
        else if (mesh.m_edges[L1].first == kend || mesh.m_edges[L1].second == kend)
        {
            knext = i;
        }
    }

    if (kself == constants::missing::uintValue || knext == constants::missing::uintValue)
    {
        return isLeftRight;
    }

    if (knext == kself + 1 || knext + mesh.m_numFacesNodes[elementId] == kself + 1)
    {
        isLeftRight = 0;
    }
    else if (kself == knext + 1 || kself + mesh.m_numFacesNodes[elementId] == knext + 1)
    {
        isLeftRight = 1;
    }

    return isLeftRight;
}

meshkernel::UInt meshkernel::CasulliRefinement::FindCommon(const Mesh2D& mesh, const UInt l1, const UInt l2)
{
    UInt commonElement = constants::missing::uintValue;

    for (UInt i = 0; i < mesh.m_edgesNumFaces[l1]; ++i)
    {
        for (UInt j = 0; j < mesh.m_edgesNumFaces[l2]; ++j)
        {
            if (mesh.m_edgesFaces[l1][i] == mesh.m_edgesFaces[l2][j])
            {
                commonElement = mesh.m_edgesFaces[l1][i];
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

void meshkernel::CasulliRefinement::ComputeNewNodesDirectional(Mesh2D& mesh, std::vector<LinkNodes>& newNodes, std::vector<int>& kc, const Point& clickedPoint)
{
}
