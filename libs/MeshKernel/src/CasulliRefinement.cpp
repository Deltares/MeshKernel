#include <vector>
#include <iostream>

#include "MeshKernel/CasulliRefinement.hpp"

void meshkernel::CasulliRefinement::Compute (Mesh2D& mesh, const MeshRefinementParameters& meshRefinementParameters[[maybe_unused]])
{
    std::vector<LinkNodes> newNodes (mesh.GetNumEdges(), {constants::missing::uintValue, constants::missing::uintValue, constants::missing::uintValue, constants::missing::uintValue});

    ComputeNewNodes (mesh, newNodes);
}

void meshkernel::CasulliRefinement::ComputeNewNodes (Mesh2D& mesh, std::vector<LinkNodes>& newNodes)
{
    std::vector<int> kc(4 * mesh.GetNumNodes(), 1);

    for (UInt i = 0; i < mesh.GetNumFaces(); ++i)
    {
        Point elementCentre = mesh.m_facesCircumcenters[i];

        for (UInt j = 0; j < mesh.m_numFacesNodes[i]; ++j)
        {

            UInt knode = mesh.m_facesNodes [i][j];

            UInt link1 = constants::missing::uintValue;
            UInt link2 = constants::missing::uintValue;
            UInt knew =  constants::missing::uintValue;

            for (UInt k = 0; k < mesh.m_facesEdges[i].size (); ++k)
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

            if (kc[knode] != 0)
            {
                Point newNode = 0.5 * (elementCentre + mesh.m_nodes [knode]);
                std::cout << " new point " << newNode.x << ", " << newNode.y << std::endl;
                knew = mesh.InsertNode (newNode);
                kc[knew] = -2;
            }
            else
            {
                knew = knode;
            }

            StoreNewNode (mesh, knode, link1, link2, knew, newNodes);
        }

        std::cout << std::endl;
    }

    [[maybe_unused]]int dummy = 1;

}

void meshkernel::CasulliRefinement::StoreNewNode (const Mesh2D& mesh, const UInt knode, const UInt link1, const UInt link2, const UInt knew, std::vector<LinkNodes>& newNodes)
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

meshkernel::UInt meshkernel::CasulliRefinement::IsStartEnd (const Mesh2D& mesh, const UInt nodeId, const UInt edgeId)
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

meshkernel::UInt meshkernel::CasulliRefinement::IsLeftRight (const Mesh2D& mesh, const UInt elementId, const UInt edgeId)
{
    UInt isLeftRight = constants::missing::uintValue;
    UInt kself = constants::missing::uintValue;
    UInt knext = constants::missing::uintValue;
    UInt kend = mesh.m_edges [edgeId].second;

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

meshkernel::UInt meshkernel::CasulliRefinement::FindCommon (const Mesh2D& mesh, const UInt l1, const UInt l2)
{
    UInt commonElement = constants::missing::uintValue;

    for (UInt i = 0; i < mesh.m_edgesNumFaces [l1]; ++i)
    {
        for (UInt j = 0; j < mesh.m_edgesNumFaces [l2]; ++j)
        {
            if (mesh.m_edgesFaces [l1][i] == mesh.m_edgesFaces [l2][j])
            {
                commonElement = mesh.m_edgesFaces [l1][i];
                break;
            }
        }
    }

    return commonElement;
}
