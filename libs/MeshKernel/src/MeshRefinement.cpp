//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2021.
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

#include <MeshKernel/AveragingInterpolation.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/MeshRefinement.hpp>
#include <MeshKernel/Operations.hpp>

using meshkernel::Mesh2D;
using meshkernel::MeshRefinement;

MeshRefinement::MeshRefinement(std::shared_ptr<Mesh2D> mesh,
                               std::shared_ptr<MeshInterpolation> interpolant,
                               const MeshRefinementParameters& meshRefinementParameters)
    : m_mesh(mesh),
      m_interpolant(interpolant)
{
    CheckMeshRefinementParameters(meshRefinementParameters);
    m_meshRefinementParameters = meshRefinementParameters;
    m_refinementType = static_cast<MeshRefinementType::Type>(m_meshRefinementParameters.refinement_type);
}

MeshRefinement::MeshRefinement(std::shared_ptr<Mesh2D> mesh,
                               std::shared_ptr<MeshInterpolation> interpolant,
                               const MeshRefinementParameters& meshRefinementParameters,
                               bool useNodalRefinement)
    : m_mesh(mesh),
      m_interpolant(interpolant),
      m_useNodalRefinement(useNodalRefinement)
{
    CheckMeshRefinementParameters(meshRefinementParameters);
    m_meshRefinementParameters = meshRefinementParameters;
    m_refinementType = static_cast<MeshRefinementType::Type>(m_meshRefinementParameters.refinement_type);
}

MeshRefinement::MeshRefinement(std::shared_ptr<Mesh2D> mesh,
                               const Polygons& polygon,
                               const MeshRefinementParameters& meshRefinementParameters)
    : m_mesh(mesh),
      m_polygons(polygon)
{
    CheckMeshRefinementParameters(meshRefinementParameters);
    m_meshRefinementParameters = meshRefinementParameters;
}

void MeshRefinement::Compute()
{
    // administrate mesh once more
    m_mesh->Administrate();

    // all faces and edges refined
    m_faceMask.resize(m_mesh->GetNumFaces(), 1);
    m_edgeMask.resize(m_mesh->GetNumEdges(), -1);

    // get bounding box
    Point lowerLeft{constants::missing::doubleValue, constants::missing::doubleValue};
    Point upperRight{constants::missing::doubleValue, constants::missing::doubleValue};
    if (m_mesh->m_projection == Projection::Spherical)
    {
        const auto boundingBox = BoundingBox(m_mesh->m_nodes);
        lowerLeft = boundingBox.lowerLeft();
        upperRight = boundingBox.upperRight();
    }

    // select the nodes to refine
    auto const isRefinementBasedOnSamples = m_interpolant != nullptr;
    if (!isRefinementBasedOnSamples && m_meshRefinementParameters.refine_intersected == 1)
    {
        const auto edgeMask = m_mesh->MaskEdgesOfFacesInPolygon(m_polygons, false, true);
        m_nodeMask = m_mesh->NodeMaskFromEdgeMask(edgeMask);
    }
    else
    {
        m_nodeMask = m_mesh->NodeMaskFromPolygon(m_polygons, true);
    }

    FindBrotherEdges();

    // set_initial_mask
    ComputeNodeMaskAtPolygonPerimeter();

    auto numFacesAfterRefinement = m_mesh->GetNumFaces();
    // reserve some extra capacity for the node mask
    m_nodeMask.reserve(m_nodeMask.size() * 2);
    for (auto level = 0; level < m_meshRefinementParameters.max_num_refinement_iterations; level++)
    {
        // Compute all edge lengths at once
        ComputeEdgeBelowMinSizeAfterRefinement();

        // computes the edge and face refinement mask from samples
        if (isRefinementBasedOnSamples)
        {
            ComputeRefinementMasksFromSamples();
        }
        else
        {
            ComputeRefinementMaskFromEdgeSize();
        }

        if (level == 0)
        {
            // if one face node is in polygon enable face refinement
            for (UInt f = 0; f < m_mesh->GetNumFaces(); ++f)
            {
                bool activeNodeFound = false;
                for (UInt n = 0; n < m_mesh->GetNumFaceEdges(f); ++n)
                {
                    const auto nodeIndex = m_mesh->m_facesNodes[f][n];
                    if (m_nodeMask[nodeIndex] != 0 && m_nodeMask[nodeIndex] != -2)
                    {
                        activeNodeFound = true;
                        break;
                    }
                }
                if (!activeNodeFound)
                {
                    m_faceMask[f] = 0;
                }
            }
        }
        if (level > 0)
        {
            // if one face node is not in polygon disable refinement
            for (UInt f = 0; f < m_mesh->GetNumFaces(); f++)
            {
                for (UInt n = 0; n < m_mesh->GetNumFaceEdges(f); n++)
                {
                    const auto nodeIndex = m_mesh->m_facesNodes[f][n];
                    if (m_nodeMask[nodeIndex] != 1)
                    {
                        m_faceMask[f] = 0;
                        break;
                    }
                }
            }
        }

        ComputeEdgesRefinementMask();

        ComputeIfFaceShouldBeSplit();

        UInt numFacesToRefine = 0;
        for (UInt f = 0; f < m_mesh->GetNumFaces(); f++)
        {
            if (m_faceMask[f] != 0)
            {
                numFacesToRefine++;
            }
        }
        if (numFacesToRefine == 0)
        {
            break;
        }
        numFacesAfterRefinement = numFacesAfterRefinement * 4;

        // spit the edges
        RefineFacesBySplittingEdges();

        m_mesh->OffsetSphericalCoordinates(lowerLeft.x, upperRight.x);

        m_mesh->Administrate();

        m_faceMask.resize(m_mesh->GetNumFaces());
        m_edgeMask.resize(m_mesh->GetNumEdges());

        FindBrotherEdges();
    }

    // remove isolated hanging nodes and connect if needed
    if (m_meshRefinementParameters.connect_hanging_nodes == 1)
    {
        ConnectHangingNodes();
        m_mesh->Administrate();
    }
}

void MeshRefinement::ComputeRefinementMaskFromEdgeSize()
{
    std::ranges::fill(m_edgeMask, 0);
    std::ranges::fill(m_faceMask, 0);
    for (UInt f = 0; f < m_mesh->GetNumFaces(); ++f)
    {
        for (UInt n = 0; n < m_mesh->GetNumFaceEdges(f); ++n)
        {
            const auto e = m_mesh->m_facesEdges[f][n];
            if (!m_isEdgeBelowMinSizeAfterRefinement[e])
            {
                m_edgeMask[e] = -1;
                m_faceMask[f] = 1;
            }
        }
    }
}

meshkernel::UInt MeshRefinement::DeleteIsolatedHangingnodes()
{

    UInt numRemovedIsolatedHangingNodes = 0;
    for (UInt e = 0; e < m_mesh->GetNumEdges(); ++e)
    {
        const auto brotherEdgeIndex = m_brotherEdges[e];
        if (brotherEdgeIndex == constants::missing::uintValue)
        {
            continue;
        }

        const auto commonNode = m_mesh->FindCommonNode(e, brotherEdgeIndex);
        if (commonNode == constants::missing::uintValue)
        {
            continue;
        }

        if (commonNode > 0 && m_mesh->m_nodesNumEdges[commonNode] == 2)
        {
            for (UInt f = 0; f < m_mesh->m_edgesNumFaces[e]; ++f)
            {
                const auto faceIndex = m_mesh->m_edgesFaces[e][f];

                if (faceIndex != m_mesh->m_edgesFaces[brotherEdgeIndex][0] &&
                    faceIndex != m_mesh->m_edgesFaces[brotherEdgeIndex][std::min(m_mesh->m_edgesNumFaces[brotherEdgeIndex], static_cast<UInt>(1))])
                {
                    throw AlgorithmError("Algorithm error.");
                }

                UInt ee = 0;
                UInt nn = 0;
                for (UInt n = 0; n < m_mesh->GetNumFaceEdges(faceIndex); ++n)
                {
                    const auto edgeIndex = m_mesh->m_facesEdges[faceIndex][n];
                    if (edgeIndex != brotherEdgeIndex)
                    {
                        m_mesh->m_facesEdges[faceIndex][ee] = edgeIndex;
                        ee++;
                    }

                    const auto nodeIndex = m_mesh->m_facesEdges[faceIndex][n];
                    if (nodeIndex != commonNode)
                    {
                        m_mesh->m_facesNodes[faceIndex][nn] = nodeIndex;
                        nn++;
                    }
                }

                m_mesh->m_numFacesNodes[faceIndex] -= 1;

                if (m_mesh->m_numFacesNodes[faceIndex] != ee || m_mesh->m_numFacesNodes[faceIndex] != nn)
                {
                    throw AlgorithmError("Algorithm error.");
                }
            }

            const auto otherNodeIndex = OtherNodeOfEdge(m_mesh->m_edges[brotherEdgeIndex], commonNode);

            // update lin admin
            if (m_mesh->m_edges[e].first == commonNode)
            {
                m_mesh->m_edges[e].first = otherNodeIndex;
            }
            else
            {
                m_mesh->m_edges[e].second = otherNodeIndex;
            }

            // change nod adm of other node
            for (UInt ee = 0; ee < m_mesh->m_nodesNumEdges[otherNodeIndex]; ++ee)
            {
                if (m_mesh->m_nodesEdges[otherNodeIndex][ee] == brotherEdgeIndex)
                {
                    m_mesh->m_nodesEdges[otherNodeIndex][ee] = e;
                    break;
                }
            }

            // delete node
            m_mesh->DeleteNode(commonNode);

            m_mesh->DeleteEdge(brotherEdgeIndex);

            m_brotherEdges[brotherEdgeIndex] = constants::missing::uintValue;

            numRemovedIsolatedHangingNodes++;
        }
    }
    return numRemovedIsolatedHangingNodes;
}

void MeshRefinement::ConnectHangingNodes()
{
    std::vector edgeEndNodeCache(Mesh::m_maximumNumberOfNodesPerFace, constants::missing::uintValue);
    std::vector hangingNodeCache(Mesh::m_maximumNumberOfNodesPerFace, constants::missing::uintValue);

    for (UInt f = 0; f < m_mesh->GetNumFaces(); ++f)
    {
        std::ranges::fill(edgeEndNodeCache, constants::missing::uintValue);
        std::ranges::fill(hangingNodeCache, constants::missing::uintValue);
        const auto numEdges = m_mesh->GetNumFaceEdges(f);
        if (numEdges > Mesh::m_maximumNumberOfNodesPerFace)
        {
            continue;
        }

        UInt numNonHangingNodes = 0;
        for (UInt n = 0; n < numEdges; ++n)
        {
            const auto e = NextCircularBackwardIndex(n, numEdges);
            const auto ee = NextCircularForwardIndex(n, numEdges);

            const auto edgeIndex = m_mesh->m_facesEdges[f][n];
            const auto firstEdgeIndex = m_mesh->m_facesEdges[f][e];
            const auto secondEdgeIndex = m_mesh->m_facesEdges[f][ee];
            if (m_brotherEdges[edgeIndex] == secondEdgeIndex)
            {
                continue;
            }

            if (numNonHangingNodes > Mesh::m_maximumNumberOfNodesPerFace - 1)
            {
                return;
            }

            edgeEndNodeCache[numNonHangingNodes] = m_mesh->FindCommonNode(edgeIndex, secondEdgeIndex);
            if (edgeEndNodeCache[numNonHangingNodes] == constants::missing::uintValue)
            {
                throw AlgorithmError("Could not find common node.");
            }

            if (m_brotherEdges[edgeIndex] == firstEdgeIndex)
            {
                hangingNodeCache[numNonHangingNodes] = m_mesh->FindCommonNode(edgeIndex, firstEdgeIndex);
                if (hangingNodeCache[numNonHangingNodes] == constants::missing::uintValue)
                {
                    throw AlgorithmError("Could not find common node.");
                }
            }
            numNonHangingNodes++;
        }

        const auto numHangingNodes = numEdges - numNonHangingNodes;
        if (numHangingNodes == 0)
            continue;

        // Quads
        if (numNonHangingNodes == Mesh::m_numNodesQuads)
        {
            switch (numHangingNodes)
            {
            case 1: // one hanging node
                for (UInt n = 0; n < numNonHangingNodes; ++n)
                {
                    if (hangingNodeCache[n] == constants::missing::uintValue)
                    {
                        continue;
                    }

                    auto ee = NextCircularBackwardIndex(n, numNonHangingNodes);
                    ee = NextCircularBackwardIndex(ee, numNonHangingNodes);
                    const auto eee = NextCircularForwardIndex(n, numNonHangingNodes);
                    m_mesh->ConnectNodes(edgeEndNodeCache[ee], hangingNodeCache[n]);
                    m_mesh->ConnectNodes(edgeEndNodeCache[eee], hangingNodeCache[n]);

                    break;
                }
                break;
            case 2: // two hanging node
                for (UInt n = 0; n < numNonHangingNodes; ++n)
                {
                    if (hangingNodeCache[n] == constants::missing::uintValue)
                    {
                        continue;
                    }

                    const auto e = NextCircularBackwardIndex(n, numNonHangingNodes);
                    const auto ee = NextCircularForwardIndex(n, numNonHangingNodes);
                    const auto eee = NextCircularForwardIndex(n + 1, numNonHangingNodes);
                    if (hangingNodeCache[e] != constants::missing::uintValue) // left neighbor
                    {
                        m_mesh->ConnectNodes(hangingNodeCache[e], hangingNodeCache[n]);
                        m_mesh->ConnectNodes(hangingNodeCache[n], edgeEndNodeCache[ee]);
                        m_mesh->ConnectNodes(edgeEndNodeCache[ee], hangingNodeCache[e]);
                    }
                    else if (hangingNodeCache[ee] != constants::missing::uintValue) // right neighbor
                    {
                        m_mesh->ConnectNodes(hangingNodeCache[n], hangingNodeCache[ee]);
                        m_mesh->ConnectNodes(hangingNodeCache[ee], edgeEndNodeCache[eee]);
                        m_mesh->ConnectNodes(edgeEndNodeCache[eee], hangingNodeCache[n]);
                    }
                    else if (hangingNodeCache[eee] != constants::missing::uintValue) // hanging nodes must be opposing
                    {
                        m_mesh->ConnectNodes(hangingNodeCache[n], hangingNodeCache[eee]);
                    }
                    break;
                }
                break;
            default:
                break;
            }
        }
        else if (numNonHangingNodes == 3)
        {
            switch (numHangingNodes)
            {
            case 1: // one hanging node
                for (UInt n = 0; n < numNonHangingNodes; ++n)
                {
                    if (hangingNodeCache[n] == constants::missing::uintValue)
                    {
                        continue;
                    }
                    const auto e = NextCircularForwardIndex(n, numNonHangingNodes);
                    m_mesh->ConnectNodes(hangingNodeCache[n], edgeEndNodeCache[e]);
                    break;
                }
                break;
            case 2: // two hanging node
                for (UInt n = 0; n < numNonHangingNodes; ++n)
                {
                    if (hangingNodeCache[n] == constants::missing::uintValue)
                    {
                        continue;
                    }
                    const auto e = NextCircularBackwardIndex(n, numNonHangingNodes);
                    const auto ee = NextCircularForwardIndex(n, numNonHangingNodes);
                    if (hangingNodeCache[e] != constants::missing::uintValue) // left neighbor
                    {
                        m_mesh->ConnectNodes(hangingNodeCache[n], hangingNodeCache[e]);
                    }
                    else
                    {
                        m_mesh->ConnectNodes(hangingNodeCache[n], hangingNodeCache[ee]);
                    }
                    break;
                }
                break;
            default:
                break;
            }
        }
        else
        {
            throw std::invalid_argument("MeshRefinement::connect_hanging_nodes: The number of non-hanging nodes is neither 3 nor 4.");
        }
    }
}

void MeshRefinement::RefineFacesBySplittingEdges()
{
    const auto numEdgesBeforeRefinement = m_mesh->GetNumEdges();

    // Add new nodes where required
    std::vector<UInt> notHangingFaceNodes;
    notHangingFaceNodes.reserve(Mesh::m_maximumNumberOfNodesPerFace);
    std::vector<UInt> nonHangingEdges;
    nonHangingEdges.reserve(Mesh::m_maximumNumberOfNodesPerFace);

    std::vector<Point> facePolygonWithoutHangingNodes;
    facePolygonWithoutHangingNodes.reserve(Mesh::m_maximumNumberOfNodesPerFace);
    std::vector<UInt> localEdgesNumFaces;
    localEdgesNumFaces.reserve(Mesh::m_maximumNumberOfEdgesPerFace);

    for (UInt e = 0; e < m_mesh->GetNumEdges(); e++)
    {
        if (m_edgeMask[e] == 0)
        {
            continue;
        }

        // Compute the center of the edge
        const auto firstNodeIndex = m_mesh->m_edges[e].first;
        const auto secondNodeIndex = m_mesh->m_edges[e].second;
        const auto firstNode = m_mesh->m_nodes[firstNodeIndex];
        const auto secondNode = m_mesh->m_nodes[secondNodeIndex];

        Point middle{(firstNode.x + secondNode.x) * 0.5, (firstNode.y + secondNode.y) * 0.5};
        if (m_mesh->m_projection == Projection::Spherical)
        {

            middle.y = (firstNode.y + secondNode.y) / 2.0;
            if (std::abs(firstNode.x - secondNode.x) > 180.0)
            {
                middle.x += 180.0;
            }

            // fix at poles
            const auto firstNodeAtPole = IsPointOnPole(firstNode);
            const auto secondNodeAtPole = IsPointOnPole(secondNode);
            if (firstNodeAtPole && !secondNodeAtPole)
            {
                middle.x = secondNode.x;
            }
            else if (!firstNodeAtPole && secondNodeAtPole)
            {
                middle.x = firstNode.x;
            }
        }

        const auto newNodeIndex = m_mesh->InsertNode(middle);
        m_edgeMask[e] = static_cast<int>(newNodeIndex);

        // set mask on the new node
        m_nodeMask.emplace_back(1);
        if (m_nodeMask[firstNodeIndex] == 0 && m_nodeMask[secondNodeIndex] == 0)
        {
            m_nodeMask[newNodeIndex] = 0;
        }
        else if (m_nodeMask[firstNodeIndex] != 1 || m_nodeMask[secondNodeIndex] != 1)
        {
            m_nodeMask[newNodeIndex] = -1;
        }
        else
        {
            m_nodeMask[newNodeIndex] = 1;
        }
    }

    for (UInt f = 0; f < m_mesh->GetNumFaces(); f++)
    {
        if (m_faceMask[f] == 0)
        {
            continue;
        }

        const auto numEdges = m_mesh->GetNumFaceEdges(f);
        // check if the parent face is crossed by the enclosing polygon
        bool isParentCrossed = false;
        for (UInt e = 0; e < numEdges; ++e)
        {
            const auto n = m_mesh->m_facesNodes[f][e];
            if (m_nodeMask[n] != 1)
            {
                isParentCrossed = true;
                break;
            }
        }

        m_mesh->ComputeFaceClosedPolygonWithLocalMappings(f, m_polygonNodesCache, m_localNodeIndicesCache, m_globalEdgeIndicesCache);

        UInt numBrotherEdges = 0;
        notHangingFaceNodes.clear();
        nonHangingEdges.clear();

        for (UInt e = 0; e < numEdges; e++)
        {
            const auto firstEdge = NextCircularBackwardIndex(e, numEdges);
            const auto secondEdge = NextCircularForwardIndex(e, numEdges);

            auto mappedEdge = m_localNodeIndicesCache[e];
            const auto edgeIndex = m_mesh->m_facesEdges[f][mappedEdge];

            mappedEdge = m_localNodeIndicesCache[firstEdge];
            const auto firstEdgeIndex = m_mesh->m_facesEdges[f][mappedEdge];

            mappedEdge = m_localNodeIndicesCache[secondEdge];
            const auto secondEdgeIndex = m_mesh->m_facesEdges[f][mappedEdge];

            if (edgeIndex == constants::missing::uintValue)
            {
                continue;
            }

            if (m_brotherEdges[edgeIndex] == secondEdgeIndex && secondEdgeIndex != constants::missing::uintValue)
            {
                numBrotherEdges++;
                const auto newNode = m_mesh->FindCommonNode(edgeIndex, m_brotherEdges[edgeIndex]);
                if (newNode == constants::missing::uintValue)
                {
                    throw AlgorithmError("Could not find common node.");
                }

                notHangingFaceNodes.emplace_back(newNode);
            }
            else if ((m_brotherEdges[edgeIndex] != firstEdgeIndex || m_brotherEdges[edgeIndex] == constants::missing::uintValue) && m_edgeMask[edgeIndex] != 0)
            {
                notHangingFaceNodes.emplace_back(m_edgeMask[edgeIndex]);
            }

            if (notHangingFaceNodes.size() >= Mesh::m_maximumNumberOfNodesPerFace)
            {
                return;
            }

            // check if start of this link is hanging
            if (m_brotherEdges[edgeIndex] != firstEdgeIndex || firstEdgeIndex != constants::missing::uintValue)
            {
                nonHangingEdges.emplace_back(e);
            }
        }

        // compute new center node : circumcenter without hanging nodes for quads, c / g otherwise
        facePolygonWithoutHangingNodes.clear();
        localEdgesNumFaces.clear();
        for (const auto& edge : nonHangingEdges)
        {
            facePolygonWithoutHangingNodes.emplace_back(m_polygonNodesCache[edge]);

            const auto mappedEdge = m_localNodeIndicesCache[edge];
            const auto edgeIndex = m_mesh->m_facesEdges[f][mappedEdge];
            if (edgeIndex != constants::missing::uintValue)
            {
                localEdgesNumFaces.emplace_back(m_mesh->m_edgesNumFaces[edgeIndex]);
            }
            else
            {
                localEdgesNumFaces.emplace_back(1);
            }
        }

        // quads
        Point splittingNode(m_mesh->m_facesMassCenters[f]);
        if (localEdgesNumFaces.size() == Mesh::m_numNodesQuads && m_meshRefinementParameters.use_mass_center_when_refining == 0)
        {
            // close the polygon before computing the face circumcenter
            facePolygonWithoutHangingNodes.emplace_back(facePolygonWithoutHangingNodes.front());
            localEdgesNumFaces.emplace_back(localEdgesNumFaces.front());

            splittingNode = m_mesh->ComputeFaceCircumenter(facePolygonWithoutHangingNodes,
                                                           localEdgesNumFaces);

            if (m_mesh->m_projection == Projection::Spherical)
            {
                auto miny = std::numeric_limits<double>::max();
                auto maxy = std::numeric_limits<double>::lowest();
                for (const auto& node : facePolygonWithoutHangingNodes)
                {
                    miny = std::min(node.y, miny);
                    maxy = std::max(node.y, maxy);
                }

                const auto middlelatitude = (miny + maxy) / 2.0;
                const auto ydiff = maxy - miny;
                if (ydiff > 1e-8)
                {
                    splittingNode.y = miny + 2.0 * (middlelatitude - miny) / ydiff * (splittingNode.y - miny);
                }
            }
        }

        if (localEdgesNumFaces.size() >= Mesh::m_numNodesQuads)
        {
            if (notHangingFaceNodes.size() > 2)
            {
                const auto newNodeIndex = m_mesh->InsertNode(splittingNode);

                for (const auto& notHangingNode : notHangingFaceNodes)
                {
                    m_mesh->ConnectNodes(notHangingNode, newNodeIndex);
                }

                m_nodeMask.emplace_back(1);
                if (isParentCrossed)
                {
                    // inactive nodes in cells crossed by polygon
                    m_nodeMask[newNodeIndex] = -1;
                }
            }
            else if (notHangingFaceNodes.size() == 2)
            {
                m_mesh->ConnectNodes(notHangingFaceNodes[0], notHangingFaceNodes[1]);
            }
        }
        else
        {
            for (UInt n = 0; n < notHangingFaceNodes.size(); ++n)
            {
                const auto nn = NextCircularForwardIndex(n, static_cast<UInt>(notHangingFaceNodes.size()));
                m_mesh->ConnectNodes(notHangingFaceNodes[n], notHangingFaceNodes[nn]);
            }
        }
    }

    // Split original edges
    for (UInt e = 0; e < numEdgesBeforeRefinement; ++e)
    {
        if (m_edgeMask[e] > 0)
        {
            const auto newEdgeIndex = m_mesh->ConnectNodes(m_edgeMask[e], m_mesh->m_edges[e].second);
            m_mesh->m_edges[e].second = m_edgeMask[e];
            m_brotherEdges.resize(m_mesh->GetNumEdges());
            m_brotherEdges[newEdgeIndex] = e;
            m_brotherEdges[e] = newEdgeIndex;
        }
    }
}

void MeshRefinement::ComputeNodeMaskAtPolygonPerimeter()
{
    for (UInt f = 0; f < m_mesh->GetNumFaces(); f++)
    {
        bool crossing = false;
        const auto numnodes = m_mesh->GetNumFaceEdges(f);
        for (UInt n = 0; n < numnodes; n++)
        {
            const auto nodeIndex = m_mesh->m_facesNodes[f][n];
            if (m_nodeMask[nodeIndex] == 0)
            {
                crossing = true;
                break;
            }
        }

        if (crossing)
        {
            m_faceMask[f] = 0;
            for (UInt n = 0; n < numnodes; n++)
            {
                const auto nodeIndex = m_mesh->m_facesNodes[f][n];
                if (m_nodeMask[nodeIndex] == 1)
                {
                    m_nodeMask[nodeIndex] = -2;
                }
            }
        }
    }
}

void MeshRefinement::ComputeRefinementMasksFromSamples()
{
    std::ranges::fill(m_edgeMask, 0);
    std::ranges::fill(m_faceMask, 0);

    m_polygonNodesCache.resize(Mesh::m_maximumNumberOfNodesPerFace + 1);
    m_localNodeIndicesCache.resize(Mesh::m_maximumNumberOfNodesPerFace + 1, constants::missing::uintValue);
    m_globalEdgeIndicesCache.resize(Mesh::m_maximumNumberOfEdgesPerFace + 1, constants::missing::uintValue);
    m_refineEdgeCache.resize(Mesh::m_maximumNumberOfEdgesPerFace, 0);

    // Compute all interpolated values
    m_interpolant->Compute();

    for (UInt f = 0; f < m_mesh->GetNumFaces(); f++)
    {
        FindHangingNodes(f);

        ComputeRefinementMasksFromSamples(f);
    }

    for (auto& edge : m_edgeMask)
    {
        edge = -edge;
    }
    SmoothRefinementMasks();
}

void MeshRefinement::FindHangingNodes(UInt face)
{
    const auto numFaceNodes = m_mesh->GetNumFaceEdges(face);

    if (numFaceNodes > Mesh::m_maximumNumberOfEdgesPerNode)
    {
        throw AlgorithmError("The number of face nodes is greater than the maximum number of edges per node.");
    }

    m_isHangingNodeCache.resize(Mesh::m_maximumNumberOfNodesPerFace);
    m_isHangingEdgeCache.resize(Mesh::m_maximumNumberOfEdgesPerFace);
    std::fill(m_isHangingNodeCache.begin(), m_isHangingNodeCache.end(), false);
    std::fill(m_isHangingEdgeCache.begin(), m_isHangingEdgeCache.end(), false);

    auto kknod = numFaceNodes;
    for (UInt n = 0; n < numFaceNodes; n++)
    {
        const auto edgeIndex = m_mesh->m_facesEdges[face][n];

        // check if the parent edge is in the cell
        if (m_brotherEdges[edgeIndex] != constants::missing::uintValue)
        {
            const auto e = NextCircularBackwardIndex(n, numFaceNodes);
            const auto ee = NextCircularForwardIndex(n, numFaceNodes);
            const auto firstEdgeIndex = m_mesh->m_facesEdges[face][e];
            const auto secondEdgeIndex = m_mesh->m_facesEdges[face][ee];

            UInt commonNode = constants::missing::uintValue;
            if (m_brotherEdges[edgeIndex] == firstEdgeIndex)
            {
                commonNode = m_mesh->FindCommonNode(edgeIndex, firstEdgeIndex);
                if (commonNode == constants::missing::uintValue)
                {
                    throw AlgorithmError("Could not find common node.");
                }
            }
            else if (m_brotherEdges[edgeIndex] == secondEdgeIndex)
            {
                commonNode = m_mesh->FindCommonNode(edgeIndex, secondEdgeIndex);
                if (commonNode == constants::missing::uintValue)
                {
                    throw AlgorithmError("Could not find common node.");
                }
            }

            if (commonNode != constants::missing::uintValue)
            {
                m_isHangingEdgeCache[n] = true;
                for (UInt nn = 0; nn < numFaceNodes; nn++)
                {
                    kknod = NextCircularForwardIndex(kknod, numFaceNodes);

                    if (m_mesh->m_facesNodes[face][kknod] == commonNode && !m_isHangingNodeCache[kknod])
                    {
                        m_isHangingNodeCache[kknod] = true;
                        break;
                    }
                }
            }
        }
    }
}

meshkernel::UInt MeshRefinement::CountHangingNodes() const
{
    meshkernel::UInt result = 0;
    for (const auto& v : m_isHangingNodeCache)
    {
        if (v)
        {
            result++;
        }
    }
    return result;
}

meshkernel::UInt MeshRefinement::CountHangingEdges() const
{
    meshkernel::UInt result = 0;
    for (const auto& v : m_isHangingEdgeCache)
    {
        if (v)
        {
            result++;
        }
    }
    return result;
}

meshkernel::UInt MeshRefinement::CountEdgesToRefine(UInt face) const
{
    const auto numFaceNodes = m_mesh->GetNumFaceEdges(face);

    UInt result = 0;

    for (UInt n = 0; n < numFaceNodes; n++)
    {
        const auto edgeIndex = m_mesh->m_facesEdges[face][n];
        if (m_edgeMask[edgeIndex] != 0)
        {
            result += 1;
        }
    }
    return result;
}

void MeshRefinement::ComputeRefinementMasksFromSamples(UInt face)
{

    if (IsEqual(m_interpolant->GetFaceResult(face), constants::missing::doubleValue))
    {
        return;
    }

    size_t numEdgesToBeRefined = 0;
    std::ranges::fill(m_refineEdgeCache, 0);

    switch (m_refinementType)
    {
    case MeshRefinementType::RefinementLevels:
        if (m_interpolant->GetFaceResult(face) <= 0)
        {
            return;
        }
        for (UInt i = 0; i < m_mesh->GetNumFaceEdges(face); i++)
        {
            numEdgesToBeRefined++;
            m_refineEdgeCache[i] = 1;
        }
        break;
    case MeshRefinementType::WaveCourant:
        if (m_useNodalRefinement)
        {
            ComputeFaceLocationTypes();
        }
        for (size_t e = 0; e < m_mesh->GetNumFaceEdges(face); ++e)
        {
            const auto edge = m_mesh->m_facesEdges[face][e];
            if (m_mesh->m_edgeLengths[edge] < m_mergingDistance)
            {
                numEdgesToBeRefined++;
                continue;
            }
            if (m_isEdgeBelowMinSizeAfterRefinement[edge])
            {
                continue;
            }

            bool doRefinement;
            if (m_useNodalRefinement)
            {
                if (m_faceLocationType[face] == FaceLocation::Land)
                {
                    doRefinement = false;
                }
                else if (m_faceLocationType[face] == FaceLocation::LandWater)
                {
                    doRefinement = true;
                }
                else
                {
                    const auto edgeDepth = std::abs(m_interpolant->GetEdgeResult(edge));
                    doRefinement = IsRefineNeededBasedOnCourantCriteria(edge, edgeDepth);
                }
            }
            else
            {
                const double faceDepth = m_interpolant->GetFaceResult(face);
                doRefinement = IsRefineNeededBasedOnCourantCriteria(edge, faceDepth);
            }

            if (doRefinement)
            {
                numEdgesToBeRefined++;
                m_refineEdgeCache[e] = 1;
            }
        }

        if (numEdgesToBeRefined > 0)
        {
            numEdgesToBeRefined = 0;
            for (UInt i = 0; i < m_mesh->GetNumFaceEdges(face); i++)
            {
                if (m_refineEdgeCache[i] == 1 || m_isHangingNodeCache[i])
                {
                    numEdgesToBeRefined++;
                }
            }
        }

        if (m_meshRefinementParameters.directional_refinement == 0)
        {
            if (numEdgesToBeRefined == m_mesh->GetNumFaceEdges(face))
            {
                for (UInt i = 0; i < m_mesh->GetNumFaceEdges(face); i++)
                {
                    if (!m_isHangingNodeCache[i])
                    {
                        m_refineEdgeCache[i] = 1;
                    }
                }
            }
            else
            {
                numEdgesToBeRefined = 0;
            }
        }
        break;

    default:
        throw AlgorithmError("Invalid refinement type");
    }

    // Compute face and edge masks
    if (numEdgesToBeRefined > 1)
    {
        m_faceMask[face] = 1;

        for (size_t n = 0; n < m_mesh->GetNumFaceEdges(face); n++)
        {
            if (m_refineEdgeCache[n] == 1)
            {
                const auto edgeIndex = m_mesh->m_facesEdges[face][n];
                if (edgeIndex != constants::missing::uintValue)
                {
                    m_edgeMask[edgeIndex] = 1;
                }
            }
        }
    }
}

void MeshRefinement::ComputeFaceLocationTypes()
{
    m_faceLocationType.resize(m_mesh->GetNumFaces());
    std::ranges::fill(m_faceLocationType, FaceLocation::Water);
    for (UInt face = 0; face < m_mesh->GetNumFaces(); face++)
    {
        double maxVal = -std::numeric_limits<double>::max();
        double minVal = std::numeric_limits<double>::max();
        for (UInt e = 0; e < m_mesh->GetNumFaceEdges(face); ++e)
        {
            const auto node = m_mesh->m_facesNodes[face][e];
            const auto val = m_interpolant->GetNodeResult(node);
            maxVal = std::max(maxVal, val);
            minVal = std::min(minVal, val);
        }

        if (minVal > 0.0)
        {
            m_faceLocationType[face] = FaceLocation::Land;
        }
        if (maxVal >= 0.0 && (!IsEqual(minVal, constants::missing::doubleValue) && minVal < 0.0))
        {
            m_faceLocationType[face] = FaceLocation::LandWater;
        }
    }
}

void MeshRefinement::ComputeEdgeBelowMinSizeAfterRefinement()
{
    m_mesh->ComputeEdgesLengths();
    m_isEdgeBelowMinSizeAfterRefinement.resize(m_mesh->GetNumEdges());
    for (UInt e = 0; e < m_mesh->GetNumEdges(); e++)
    {
        const double newEdgeLength = 0.5 * m_mesh->m_edgeLengths[e];
        m_isEdgeBelowMinSizeAfterRefinement[e] = newEdgeLength < m_meshRefinementParameters.min_edge_size;
    }
}

bool MeshRefinement::IsRefineNeededBasedOnCourantCriteria(UInt edge, double depthValues) const
{
    const double maxDtCourant = 120.0;
    const double celerity = constants::physical::sqrt_gravity * std::sqrt(std::abs(depthValues));
    const double waveCourant = celerity * maxDtCourant / m_mesh->m_edgeLengths[edge];
    bool doRefinement = waveCourant < 1.0;
    return doRefinement;
}

void MeshRefinement::ComputeEdgesRefinementMask()
{
    bool repeat = true;
    UInt iter = 0;
    const UInt numMaxIterations = 6;
    std::vector<int> isQuadEdge(Mesh::m_numNodesQuads);
    std::vector<UInt> numOfEdges(Mesh::m_maximumNumberOfEdgesPerFace);

    while (repeat && iter < numMaxIterations)
    {
        repeat = false;
        iter++;

        for (UInt f = 0; f < m_mesh->GetNumFaces(); f++)
        {
            if (m_faceMask[f] == 0)
            {
                continue;
            }

            const auto numHangingNodes = CountHangingEdges();

            const auto numFaceNodes = m_mesh->GetNumFaceEdges(f);

            // non-quads
            const auto numNodesEffective = numFaceNodes - numHangingNodes;
            if (numNodesEffective != Mesh::m_numNodesQuads)
            {
                for (UInt n = 0; n < numFaceNodes; n++)
                {
                    const auto e = NextCircularBackwardIndex(n, numFaceNodes);
                    const auto ee = NextCircularForwardIndex(n, numFaceNodes);

                    const auto edgeIndex = m_mesh->m_facesEdges[f][n];

                    const auto firstEdgeIndex = m_mesh->m_facesEdges[f][e];
                    const auto secondEdgeIndex = m_mesh->m_facesEdges[f][ee];

                    // do not refine edges with an hanging node
                    if (m_brotherEdges[edgeIndex] != firstEdgeIndex && m_brotherEdges[edgeIndex] != secondEdgeIndex)
                    {
                        m_edgeMask[edgeIndex] = 1;
                    }
                }
            }
            if (numNodesEffective == Mesh::m_numNodesQuads)
            {
                // number the links in the cell, links that share a hanging node will have the same number
                UInt num = 0;
                for (UInt n = 0; n < numFaceNodes; n++)
                {
                    const auto edgeIndex = m_mesh->m_facesEdges[f][n];
                    numOfEdges[n] = num;

                    if (m_edgeMask[edgeIndex] != 0)
                    {
                        isQuadEdge[num] = m_edgeMask[edgeIndex];
                    }

                    const auto ee = NextCircularForwardIndex(n, numFaceNodes);

                    const auto secondEdgeIndex = m_mesh->m_facesEdges[f][ee];

                    if (n != numFaceNodes - 1 && m_brotherEdges[edgeIndex] != secondEdgeIndex)
                    {
                        num++;
                    }
                    if (m_brotherEdges[edgeIndex] == secondEdgeIndex)
                    {
                        isQuadEdge[num] = 1;
                    }
                }

                if (num + 1 != Mesh::m_numNodesQuads)
                {
                    throw AlgorithmError("The number the links in the cell is not equal to 3.");
                }

                UInt numEdgesToRefine = 0;
                UInt firstEdgeIndex = 0;
                UInt secondEdgeIndex = 0;
                for (UInt i = 0; i < Mesh::m_numNodesQuads; i++)
                {
                    if (isQuadEdge[i] != 0)
                    {
                        numEdgesToRefine++;
                        if (firstEdgeIndex == 0)
                        {
                            firstEdgeIndex = i;
                        }
                        else if (secondEdgeIndex == 0)
                        {
                            secondEdgeIndex = i;
                        }
                    }
                }

                const auto edgeIndexDifference = secondEdgeIndex - firstEdgeIndex;
                bool refineAllEdges = false;
                if (numEdgesToRefine == 2 && (edgeIndexDifference == 1 || edgeIndexDifference == 3))
                {
                    repeat = true;
                    refineAllEdges = true;
                }

                for (UInt n = 0; n < numFaceNodes; n++)
                {
                    const auto edgeIndex = m_mesh->m_facesEdges[f][n];
                    if (m_edgeMask[edgeIndex] > 0)
                    {
                        continue;
                    }

                    if (refineAllEdges != true && m_edgeMask[edgeIndex] != -1)
                    {
                        continue;
                    }

                    const auto e = NextCircularBackwardIndex(n, numFaceNodes);
                    const auto ee = NextCircularForwardIndex(n, numFaceNodes);

                    if (numOfEdges[n] != numOfEdges[e] && numOfEdges[n] != numOfEdges[ee])
                    {
                        m_edgeMask[edgeIndex] = 1;
                    }
                }
            }
        }
    }

    if (repeat)
    {
        throw AlgorithmError("Solution did not converge.");
    }

    // only keep m_edgeMask = 1, set other values to 0
    for (auto& value : m_edgeMask)
    {
        if (value != 1)
        {
            value = 0;
        }
    }
}

void MeshRefinement::ComputeIfFaceShouldBeSplit()
{
    const UInt maxiter = 1000;
    UInt num = 1;
    UInt iter = 0;
    while (num != 0)
    {
        iter++;
        if (iter > maxiter)
        {
            break;
        }

        num = 0;
        for (UInt f = 0; f < m_mesh->GetNumFaces(); f++)
        {
            if (m_faceMask[f] != 0 && m_faceMask[f] != -1)
            {
                continue;
            }

            const auto numEdgesToRefine = CountEdgesToRefine(f);

            FindHangingNodes(f);
            const auto numHangingEdges = CountHangingEdges();
            const auto numHangingNodes = CountHangingNodes();

            bool isSplittingRequired = false;

            // check if the edge has a brother edge and needs to be refined
            const auto numFaceNodes = m_mesh->GetNumFaceEdges(f);

            if (numFaceNodes > Mesh::m_maximumNumberOfEdgesPerFace)
            {
                return;
            }

            for (UInt n = 0; n < numFaceNodes; n++)
            {
                const auto edgeIndex = m_mesh->m_facesEdges[f][n];
                if (m_isHangingEdgeCache[n] && m_edgeMask[edgeIndex] > 0)
                {
                    isSplittingRequired = true;
                }
            }

            // compute the effective face type
            const auto numNodesEffective = numFaceNodes - static_cast<UInt>(static_cast<double>(numHangingEdges) / 2.0);
            if (2 * (numFaceNodes - numNodesEffective) != numHangingEdges)
            {
                // uneven number of brotherlinks
                // TODO: ADD DOT
            }

            if (numFaceNodes + numEdgesToRefine > Mesh::m_maximumNumberOfEdgesPerFace || // would result in unsupported cells after refinement
                numFaceNodes - numHangingNodes - numEdgesToRefine <= 1 ||                // faces with only one unrefined edge
                numNodesEffective == numEdgesToRefine)                                   // refine all edges
            {
                isSplittingRequired = true;
            }

            if (isSplittingRequired)
            {
                if (m_faceMask[f] != -1)
                {
                    m_faceMask[f] = 2;
                }
                else
                {
                    m_faceMask[f] = -2;
                }

                for (UInt n = 0; n < numFaceNodes; n++)
                {
                    const auto edgeIndex = m_mesh->m_facesEdges[f][n];
                    if (!m_isHangingEdgeCache[n] && m_edgeMask[edgeIndex] == 0)
                    {
                        m_edgeMask[edgeIndex] = 1;
                        num++;
                    }
                    if (iter == maxiter)
                    {
                        // TODO: ADD DOT/MESSAGES
                    }
                }
            }
        }
    }
}

void MeshRefinement::FindBrotherEdges()
{
    m_brotherEdges.resize(m_mesh->GetNumEdges());
    std::ranges::fill(m_brotherEdges, constants::missing::uintValue);

    for (UInt n = 0; n < m_mesh->GetNumNodes(); n++)
    {
        const auto numEdgesNodes = m_mesh->m_nodesNumEdges[n];
        for (UInt e = 0; e < numEdgesNodes; e++)
        {

            const auto firstEdgeIndex = m_mesh->m_nodesEdges[n][e];
            if (m_mesh->GetNumEdgesFaces(firstEdgeIndex) < 1)
            {
                continue;
            }

            const auto ee = NextCircularForwardIndex(e, numEdgesNodes);
            const auto secondEdgeIndex = m_mesh->m_nodesEdges[n][ee];
            if (m_mesh->GetNumEdgesFaces(secondEdgeIndex) < 1)
            {
                continue;
            }

            // both edges should share the same face
            const auto firstEdgeLeftFace = m_mesh->m_edgesFaces[firstEdgeIndex][0];
            const auto firstEdgeRighFace = m_mesh->GetNumEdgesFaces(firstEdgeIndex) == 1 ? firstEdgeLeftFace : m_mesh->m_edgesFaces[firstEdgeIndex][1];
            const auto secondEdgeLeftFace = m_mesh->m_edgesFaces[secondEdgeIndex][0];
            const auto secondEdgeRighFace = m_mesh->GetNumEdgesFaces(secondEdgeIndex) == 1 ? secondEdgeLeftFace : m_mesh->m_edgesFaces[secondEdgeIndex][1];

            if (firstEdgeLeftFace != secondEdgeLeftFace &&
                firstEdgeLeftFace != secondEdgeRighFace &&
                firstEdgeRighFace != secondEdgeLeftFace &&
                firstEdgeRighFace != secondEdgeRighFace)
            {
                continue;
            }

            // check if node k is in the middle
            const auto firstEdgeOtherNode = OtherNodeOfEdge(m_mesh->m_edges[firstEdgeIndex], n);
            const auto secondEdgeOtherNode = OtherNodeOfEdge(m_mesh->m_edges[secondEdgeIndex], n);
            const auto center = ComputeMiddlePointAccountingForPoles(m_mesh->m_nodes[firstEdgeOtherNode], m_mesh->m_nodes[secondEdgeOtherNode], m_mesh->m_projection);

            // compute tolerance
            const auto firstEdgeLength = ComputeDistance(m_mesh->m_nodes[firstEdgeOtherNode], m_mesh->m_nodes[n], m_mesh->m_projection);
            const auto secondEdgeLength = ComputeDistance(m_mesh->m_nodes[secondEdgeOtherNode], m_mesh->m_nodes[n], m_mesh->m_projection);
            const auto minConnectionDistance = 1e-4 * std::max(firstEdgeLength, secondEdgeLength);

            // The center of the two edges coincides with the shared node
            const auto distanceFromCentre = ComputeDistance(center, m_mesh->m_nodes[n], m_mesh->m_projection);
            if (distanceFromCentre < minConnectionDistance)
            {
                m_brotherEdges[firstEdgeIndex] = secondEdgeIndex;
                m_brotherEdges[secondEdgeIndex] = firstEdgeIndex;
            }
        }
    }
}

void MeshRefinement::SmoothRefinementMasks()
{
    if (m_meshRefinementParameters.directional_refinement == 1)
    {
        throw AlgorithmError("Directional refinement cannot be used in combination with smoothing. Please set directional refinement to off!");
    }
    if (m_meshRefinementParameters.smoothing_iterations == 0)
    {
        return;
    }

    std::vector splitEdge(m_edgeMask.size(), false);

    for (int iter = 0; iter < m_meshRefinementParameters.smoothing_iterations; ++iter)
    {
        std::fill(splitEdge.begin(), splitEdge.end(), false);

        for (UInt f = 0; f < m_mesh->GetNumFaces(); ++f)
        {
            if (m_faceMask[f] != 1)
            {
                continue;
            }

            const auto numEdges = m_mesh->GetNumFaceEdges(f);

            for (UInt e = 0; e < numEdges; ++e)
            {
                const auto edgeIndex = m_mesh->m_facesEdges[f][e];
                const auto nextEdgeIndex = NextCircularForwardIndex(e, numEdges);
                const auto previousEdgeIndex = NextCircularBackwardIndex(e, numEdges);
                const auto split = m_brotherEdges[edgeIndex] != m_mesh->m_facesEdges[f][nextEdgeIndex] &&
                                   m_brotherEdges[edgeIndex] != m_mesh->m_facesEdges[f][previousEdgeIndex];

                if (split)
                {
                    splitEdge[edgeIndex] = true;
                }
            }
        }

        // update face refinement mask
        for (UInt f = 0; f < m_mesh->GetNumFaces(); ++f)
        {
            const auto numEdges = m_mesh->GetNumFaceEdges(f);

            for (UInt e = 0; e < numEdges; ++e)
            {
                const auto edgeIndex = m_mesh->m_facesEdges[f][e];
                if (splitEdge[edgeIndex])
                {
                    m_faceMask[f] = 1;
                }
            }
        }

        // update edge refinement mask
        for (UInt f = 0; f < m_mesh->GetNumFaces(); ++f)
        {
            if (m_faceMask[f] != 1)
            {
                continue;
            }

            const auto numEdges = m_mesh->GetNumFaceEdges(f);

            for (UInt e = 0; e < numEdges; ++e)
            {
                const auto edgeIndex = m_mesh->m_facesEdges[f][e];
                const auto nextEdgeIndex = NextCircularForwardIndex(e, numEdges);
                const auto previousEdgeIndex = NextCircularBackwardIndex(e, numEdges);
                const auto split = m_brotherEdges[edgeIndex] != m_mesh->m_facesEdges[f][nextEdgeIndex] &&
                                   m_brotherEdges[edgeIndex] != m_mesh->m_facesEdges[f][previousEdgeIndex];

                if (split)
                {
                    m_edgeMask[edgeIndex] = 1;
                }
            }
        }
    }
}
