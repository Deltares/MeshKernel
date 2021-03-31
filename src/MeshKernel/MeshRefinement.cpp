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

#include <algorithm>
#include <cmath>
#include <numeric>
#include <tuple>
#include <vector>

#include <MeshKernel/AveragingInterpolation.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/MeshRefinement.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/RTree.hpp>

using meshkernel::Mesh2D;
using meshkernel::MeshRefinement;

MeshRefinement::MeshRefinement(std::shared_ptr<Mesh2D> mesh,
                               std::shared_ptr<AveragingInterpolation> averaging,
                               const meshkernelapi::SampleRefineParameters& sampleRefineParameters,
                               const meshkernelapi::InterpolationParameters& interpolationParameters) : m_mesh(mesh),
                                                                                                        m_averaging(averaging),
                                                                                                        m_sampleRefineParameters(sampleRefineParameters),
                                                                                                        m_interpolationParameters(interpolationParameters)
{
    m_refinementType = static_cast<RefinementType>(m_sampleRefineParameters.RefinementType);
};

MeshRefinement::MeshRefinement(std::shared_ptr<Mesh2D> mesh,
                               const Polygons& polygon,
                               const meshkernelapi::InterpolationParameters& interpolationParameters) : m_mesh(mesh),
                                                                                                        m_polygons(polygon),
                                                                                                        m_interpolationParameters(interpolationParameters){};

void MeshRefinement::Compute()
{
    // administrate mesh once more
    m_mesh->Administrate(Mesh2D::AdministrationOption::AdministrateMeshEdgesAndFaces);

    // all faces and edges refined
    m_faceMask.resize(m_mesh->GetNumFaces(), 1);
    m_edgeMask.resize(m_mesh->GetNumEdges(), -1);

    // get bounding box
    Point lowerLeft{doubleMissingValue, doubleMissingValue};
    Point upperRight{doubleMissingValue, doubleMissingValue};
    if (m_mesh->m_projection == Projection::spherical)
    {
        const auto boundingBox = GetBoundingBox(m_mesh->m_nodes);
        lowerLeft = std::get<0>(boundingBox);
        upperRight = std::get<1>(boundingBox);
    }

    // select the nodes to refine
    const auto isRefinementBasedOnSamples = m_averaging == nullptr ? false : true;
    if (!isRefinementBasedOnSamples && m_interpolationParameters.RefineIntersected == 1)
    {
        m_mesh->MaskFaceEdgesInPolygon(m_polygons, false, true);
        m_mesh->ComputeNodeMaskFromEdgeMask();
    }
    else
    {
        m_mesh->MaskNodesInPolygons(m_polygons, true);
    }

    FindBrotherEdges();

    //set_initial_mask
    ComputeNodeMaskAtPolygonPerimeter();

    auto numFacesAfterRefinement = m_mesh->GetNumFaces();
    for (auto level = 0; level < m_interpolationParameters.MaxNumberOfRefinementIterations; level++)
    {
        if (level > 0)
        {
            FindBrotherEdges();
        }

        // Compute all edge lengths at once
        m_mesh->ComputeEdgesLengths();

        const auto numEdgesBeforeRefinement = m_mesh->GetNumEdges();

        // computes the edge and face refinement mask from samples
        if (isRefinementBasedOnSamples)
        {
            ComputeRefinementMasksFromSamples();

            for (auto& edge : m_edgeMask)
            {
                edge = -edge;
            }

            //TODO: implement SmoothEdgeRefinementMask SmoothEdgeRefinementMask();
        }
        else
        {
            std::fill(m_faceMask.begin(), m_faceMask.end(), 1);
            std::fill(m_edgeMask.begin(), m_edgeMask.end(), -1);
        }

        if (level == 0)
        {
            //if one face node is in polygon enable face refinement
            for (auto f = 0; f < m_mesh->GetNumFaces(); ++f)
            {
                bool activeNodeFound = false;
                for (auto n = 0; n < m_mesh->GetNumFaceEdges(f); ++n)
                {
                    const auto nodeIndex = m_mesh->m_facesNodes[f][n];
                    if (m_mesh->m_nodeMask[nodeIndex] != 0 && m_mesh->m_nodeMask[nodeIndex] != -2)
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
            //if one face node is not in polygon disable refinement
            for (auto f = 0; f < m_mesh->GetNumFaces(); f++)
            {
                for (auto n = 0; n < m_mesh->GetNumFaceEdges(f); n++)
                {
                    const auto nodeIndex = m_mesh->m_facesNodes[f][n];
                    if (m_mesh->m_nodeMask[nodeIndex] != 1)
                    {
                        m_faceMask[f] = 0;
                        break;
                    }
                }
            }
        }

        ComputeEdgesRefinementMask();

        ComputeIfFaceShouldBeSplit();

        size_t numFacesToRefine = 0;
        for (auto f = 0; f < m_mesh->GetNumFaces(); f++)
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
        RefineFacesBySplittingEdges(numEdgesBeforeRefinement);

        m_mesh->OffsetSphericalCoordinates(lowerLeft.x, upperRight.x);

        m_mesh->Administrate(Mesh2D::AdministrationOption::AdministrateMeshEdgesAndFaces);

        m_faceMask.resize(m_mesh->GetNumFaces());
        m_edgeMask.resize(m_mesh->GetNumEdges());
    }

    //remove isolated hanging nodes and connect if needed
    if (m_sampleRefineParameters.ConnectHangingNodes == 1)
    {
        ConnectHangingNodes();
        m_mesh->Administrate(Mesh2D::AdministrationOption::AdministrateMeshEdgesAndFaces);
    }
}

size_t MeshRefinement::DeleteIsolatedHangingnodes()
{

    size_t numRemovedIsolatedHangingNodes = 0;
    for (auto e = 0; e < m_mesh->GetNumEdges(); ++e)
    {
        const auto brotherEdgeIndex = m_brotherEdges[e];
        if (brotherEdgeIndex == sizetMissingValue)
        {
            continue;
        }

        const auto commonNode = m_mesh->FindCommonNode(e, brotherEdgeIndex);
        if (commonNode == sizetMissingValue)
        {
            continue;
        }

        if (commonNode > 0 && m_mesh->m_nodesNumEdges[commonNode] == 2)
        {
            for (auto f = 0; f < m_mesh->m_edgesNumFaces[e]; ++f)
            {
                const auto faceIndex = m_mesh->m_edgesFaces[e][f];

                if (faceIndex != m_mesh->m_edgesFaces[brotherEdgeIndex][0] &&
                    faceIndex != m_mesh->m_edgesFaces[brotherEdgeIndex][std::min(m_mesh->m_edgesNumFaces[brotherEdgeIndex], static_cast<size_t>(1))])
                {
                    throw AlgorithmError("MeshRefinement::DeleteIsolatedHangingnodes: Algorithm error.");
                }

                size_t ee = 0;
                size_t nn = 0;
                for (auto n = 0; n < m_mesh->GetNumFaceEdges(faceIndex); ++n)
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
                    throw AlgorithmError("MeshRefinement::DeleteIsolatedHangingnodes: Algorithm error.");
                }
            }

            const auto otherNodeIndex = OtherNodeOfEdge(m_mesh->m_edges[brotherEdgeIndex], commonNode);

            //update lin admin
            if (m_mesh->m_edges[e].first == commonNode)
            {
                m_mesh->m_edges[e].first = otherNodeIndex;
            }
            else
            {
                m_mesh->m_edges[e].second = otherNodeIndex;
            }

            //change nod adm of other node
            for (auto ee = 0; ee < m_mesh->m_nodesNumEdges[otherNodeIndex]; ++ee)
            {
                if (m_mesh->m_nodesEdges[otherNodeIndex][ee] == brotherEdgeIndex)
                {
                    m_mesh->m_nodesEdges[otherNodeIndex][ee] = e;
                    break;
                }
            }

            //delete node
            m_mesh->DeleteNode(commonNode);

            m_mesh->DeleteEdge(brotherEdgeIndex);

            m_brotherEdges[brotherEdgeIndex] = sizetMissingValue;

            numRemovedIsolatedHangingNodes++;
        }
    }
    return numRemovedIsolatedHangingNodes;
}

void MeshRefinement::ConnectHangingNodes()
{
    std::vector<size_t> edgeEndNodeCache(maximumNumberOfNodesPerFace, sizetMissingValue);
    std::vector<size_t> hangingNodeCache(maximumNumberOfNodesPerFace, sizetMissingValue);

    for (auto f = 0; f < m_mesh->GetNumFaces(); ++f)
    {
        std::fill(edgeEndNodeCache.begin(), edgeEndNodeCache.end(), sizetMissingValue);
        std::fill(hangingNodeCache.begin(), hangingNodeCache.end(), sizetMissingValue);
        const auto numEdges = m_mesh->GetNumFaceEdges(f);
        if (numEdges > maximumNumberOfNodesPerFace)
        {
            continue;
        }

        size_t numNonHangingNodes = 0;
        for (auto n = 0; n < numEdges; ++n)
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

            if (numNonHangingNodes > maximumNumberOfNodesPerFace - 1)
            {
                return;
            }

            edgeEndNodeCache[numNonHangingNodes] = m_mesh->FindCommonNode(edgeIndex, secondEdgeIndex);
            if (edgeEndNodeCache[numNonHangingNodes] == sizetMissingValue)
            {
                throw AlgorithmError("MeshRefinement::ConnectHangingNodes: Could not find common node.");
            }

            if (m_brotherEdges[edgeIndex] == firstEdgeIndex)
            {
                hangingNodeCache[numNonHangingNodes] = m_mesh->FindCommonNode(edgeIndex, firstEdgeIndex);
                if (hangingNodeCache[numNonHangingNodes] == sizetMissingValue)
                {
                    throw AlgorithmError("MeshRefinement::ConnectHangingNodes: Could not find common node.");
                }
            }
            numNonHangingNodes++;
        }

        const auto numHangingNodes = numEdges - numNonHangingNodes;
        if (numHangingNodes == 0)
            continue;

        // Quads
        if (numNonHangingNodes == numNodesQuads)
        {
            switch (numHangingNodes)
            {
            case 1: // one hanging node
                for (auto n = 0; n < numNonHangingNodes; ++n)
                {
                    if (hangingNodeCache[n] == sizetMissingValue)
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
                for (auto n = 0; n < numNonHangingNodes; ++n)
                {
                    if (hangingNodeCache[n] == sizetMissingValue)
                    {
                        continue;
                    }

                    const auto e = NextCircularBackwardIndex(n, numNonHangingNodes);
                    const auto ee = NextCircularForwardIndex(n, numNonHangingNodes);
                    const auto eee = NextCircularForwardIndex(n + 1, numNonHangingNodes);
                    if (hangingNodeCache[e] != sizetMissingValue) // left neighbor
                    {
                        m_mesh->ConnectNodes(hangingNodeCache[e], hangingNodeCache[n]);
                        m_mesh->ConnectNodes(hangingNodeCache[n], edgeEndNodeCache[ee]);
                        m_mesh->ConnectNodes(edgeEndNodeCache[ee], hangingNodeCache[e]);
                    }
                    else if (hangingNodeCache[ee] != sizetMissingValue) // right neighbor
                    {
                        m_mesh->ConnectNodes(hangingNodeCache[n], hangingNodeCache[ee]);
                        m_mesh->ConnectNodes(hangingNodeCache[ee], edgeEndNodeCache[eee]);
                        m_mesh->ConnectNodes(edgeEndNodeCache[eee], hangingNodeCache[n]);
                    }
                    else if (hangingNodeCache[eee] != sizetMissingValue) // hanging nodes must be opposing
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
                for (auto n = 0; n < numNonHangingNodes; ++n)
                {
                    if (hangingNodeCache[n] == sizetMissingValue)
                    {
                        continue;
                    }
                    const auto e = NextCircularForwardIndex(n, numNonHangingNodes);
                    m_mesh->ConnectNodes(hangingNodeCache[n], edgeEndNodeCache[e]);
                    break;
                }
                break;
            case 2: // two hanging node
                for (auto n = 0; n < numNonHangingNodes; ++n)
                {
                    if (hangingNodeCache[n] == sizetMissingValue)
                    {
                        continue;
                    }
                    const auto e = NextCircularBackwardIndex(n, numNonHangingNodes);
                    const auto ee = NextCircularForwardIndex(n, numNonHangingNodes);
                    if (hangingNodeCache[e] != sizetMissingValue) // left neighbor
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
            throw std::invalid_argument("MeshRefinement::ConnectHangingNodes: The number of non-hanging nodes is neither 3 nor 4.");
        }
    }
}

void MeshRefinement::RefineFacesBySplittingEdges(size_t numEdgesBeforeRefinement)
{
    //Add new nodes where required
    std::vector<size_t> notHangingFaceNodes;
    notHangingFaceNodes.reserve(maximumNumberOfNodesPerFace);
    std::vector<size_t> nonHangingEdges;
    nonHangingEdges.reserve(maximumNumberOfNodesPerFace);

    std::vector<Point> facePolygonWithoutHangingNodes;
    facePolygonWithoutHangingNodes.reserve(maximumNumberOfNodesPerFace);
    std::vector<size_t> localEdgesNumFaces;
    localEdgesNumFaces.reserve(maximumNumberOfEdgesPerFace);

    for (auto e = 0; e < m_mesh->GetNumEdges(); e++)
    {
        if (m_edgeMask[e] == 0)
        {
            continue;
        }

        //Compute the center of the edge
        const auto firstNodeIndex = m_mesh->m_edges[e].first;
        const auto secondNodeIndex = m_mesh->m_edges[e].second;
        const auto firstNode = m_mesh->m_nodes[firstNodeIndex];
        const auto secondNode = m_mesh->m_nodes[secondNodeIndex];

        Point middle{(firstNode.x + secondNode.x) * 0.5, (firstNode.y + secondNode.y) * 0.5};
        if (m_mesh->m_projection == Projection::spherical)
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
        m_mesh->m_nodeMask[newNodeIndex] = 1;
        if (m_mesh->m_nodeMask[firstNodeIndex] == 0 && m_mesh->m_nodeMask[secondNodeIndex] == 0)
        {
            m_mesh->m_nodeMask[newNodeIndex] = 0;
        }
        else if (m_mesh->m_nodeMask[firstNodeIndex] != 1 || m_mesh->m_nodeMask[secondNodeIndex] != 1)
        {
            m_mesh->m_nodeMask[newNodeIndex] = -1;
        }
        else
        {
            m_mesh->m_nodeMask[newNodeIndex] = 1;
        }
    }

    for (auto f = 0; f < m_mesh->GetNumFaces(); f++)
    {
        if (m_faceMask[f] == 0)
        {
            continue;
        }

        const auto numEdges = m_mesh->GetNumFaceEdges(f);
        // check if the parent face is crossed by the enclosing polygon
        bool isParentCrossed = false;
        for (auto e = 0; e < numEdges; ++e)
        {
            const auto n = m_mesh->m_facesNodes[f][e];
            if (m_mesh->m_nodeMask[n] != 1)
            {
                isParentCrossed = true;
                break;
            }
        }

        m_mesh->ComputeFaceClosedPolygonWithLocalMappings(f, m_polygonNodesCache, m_localNodeIndicesCache, m_globalEdgeIndicesCache);

        size_t numBrotherEdges = 0;
        notHangingFaceNodes.clear();
        nonHangingEdges.clear();

        for (auto e = 0; e < numEdges; e++)
        {
            const auto firstEdge = NextCircularBackwardIndex(e, numEdges);
            const auto secondEdge = NextCircularForwardIndex(e, numEdges);

            auto mappedEdge = m_localNodeIndicesCache[e];
            const auto edgeIndex = m_mesh->m_facesEdges[f][mappedEdge];

            mappedEdge = m_localNodeIndicesCache[firstEdge];
            const auto firstEdgeIndex = m_mesh->m_facesEdges[f][mappedEdge];

            mappedEdge = m_localNodeIndicesCache[secondEdge];
            const auto secondEdgeIndex = m_mesh->m_facesEdges[f][mappedEdge];

            if (edgeIndex == sizetMissingValue)
            {
                continue;
            }

            if (m_brotherEdges[edgeIndex] == secondEdgeIndex && secondEdgeIndex != sizetMissingValue)
            {
                numBrotherEdges++;
                const auto newNode = m_mesh->FindCommonNode(edgeIndex, m_brotherEdges[edgeIndex]);
                if (newNode == sizetMissingValue)
                {
                    throw AlgorithmError("MeshRefinement::RefineFacesBySplittingEdges: Could not find common node.");
                }

                notHangingFaceNodes.emplace_back(newNode);
            }
            else if ((m_brotherEdges[edgeIndex] != firstEdgeIndex || m_brotherEdges[edgeIndex] == sizetMissingValue) && m_edgeMask[edgeIndex] != 0)
            {
                notHangingFaceNodes.emplace_back(m_edgeMask[edgeIndex]);
            }

            if (notHangingFaceNodes.size() >= maximumNumberOfNodesPerFace)
            {
                return;
            }

            // check if start of this link is hanging
            if (m_brotherEdges[edgeIndex] != firstEdgeIndex || firstEdgeIndex != sizetMissingValue)
            {
                nonHangingEdges.emplace_back(e);
            }
        }

        //compute new center node : circumcenter without hanging nodes for quads, c / g otherwise
        facePolygonWithoutHangingNodes.clear();
        localEdgesNumFaces.clear();
        for (const auto& edge : nonHangingEdges)
        {
            facePolygonWithoutHangingNodes.emplace_back(m_polygonNodesCache[edge]);

            const auto mappedEdge = m_localNodeIndicesCache[edge];
            const auto edgeIndex = m_mesh->m_facesEdges[f][mappedEdge];
            if (edgeIndex != sizetMissingValue)
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
        if (localEdgesNumFaces.size() == numNodesQuads && m_interpolationParameters.UseMassCenterWhenRefining == 0)
        {
            // close the polygon before computing the face circumcenter
            facePolygonWithoutHangingNodes.emplace_back(facePolygonWithoutHangingNodes.front());
            localEdgesNumFaces.emplace_back(localEdgesNumFaces.front());

            splittingNode = m_mesh->ComputeFaceCircumenter(facePolygonWithoutHangingNodes,
                                                           localEdgesNumFaces);

            if (m_mesh->m_projection == Projection::spherical)
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

        if (localEdgesNumFaces.size() >= numNodesQuads)
        {
            if (notHangingFaceNodes.size() > 2)
            {
                const auto newNodeIndex = m_mesh->InsertNode(splittingNode);

                for (const auto& notHangingNode : notHangingFaceNodes)
                {
                    m_mesh->ConnectNodes(notHangingNode, newNodeIndex);
                }

                m_mesh->m_nodeMask[newNodeIndex] = 1;
                if (isParentCrossed)
                {
                    //inactive nodes in cells crossed by polygon
                    m_mesh->m_nodeMask[newNodeIndex] = -1;
                }
            }
            else if (notHangingFaceNodes.size() == 2)
            {
                m_mesh->ConnectNodes(notHangingFaceNodes[0], notHangingFaceNodes[1]);
            }
        }
        else
        {
            for (auto n = 0; n < notHangingFaceNodes.size(); ++n)
            {
                const auto nn = NextCircularForwardIndex(n, notHangingFaceNodes.size());
                m_mesh->ConnectNodes(notHangingFaceNodes[n], notHangingFaceNodes[nn]);
            }
        }
    }

    //Split original edges
    for (auto e = 0; e < numEdgesBeforeRefinement; ++e)
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
    for (auto f = 0; f < m_mesh->GetNumFaces(); f++)
    {
        bool crossing = false;
        const auto numnodes = m_mesh->GetNumFaceEdges(f);
        for (auto n = 0; n < numnodes; n++)
        {
            const auto nodeIndex = m_mesh->m_facesNodes[f][n];
            if (m_mesh->m_nodeMask[nodeIndex] == 0)
            {
                crossing = true;
                break;
            }
        }

        if (crossing)
        {
            m_faceMask[f] = 0;
            for (auto n = 0; n < numnodes; n++)
            {
                const auto nodeIndex = m_mesh->m_facesNodes[f][n];
                if (m_mesh->m_nodeMask[nodeIndex] == 1)
                {
                    m_mesh->m_nodeMask[nodeIndex] = -2;
                }
            }
        }
    }
}

void MeshRefinement::ComputeRefinementMasksFromSamples()
{
    std::fill(m_edgeMask.begin(), m_edgeMask.end(), 0);
    std::fill(m_faceMask.begin(), m_faceMask.end(), 0);
    m_polygonNodesCache.resize(maximumNumberOfNodesPerFace + 1);
    m_localNodeIndicesCache.resize(maximumNumberOfNodesPerFace + 1, sizetMissingValue);
    m_globalEdgeIndicesCache.resize(maximumNumberOfEdgesPerFace + 1, sizetMissingValue);
    std::vector<size_t> refineEdgeCache(maximumNumberOfEdgesPerFace);

    // Compute all interpolated values
    m_averaging->Compute();

    for (auto f = 0; f < m_mesh->GetNumFaces(); f++)
    {

        size_t numHangingEdges;
        size_t numHangingNodes;
        size_t numEdgesToRefine;
        FindHangingNodes(f, numHangingEdges, numHangingNodes, numEdgesToRefine);

        std::fill(refineEdgeCache.begin(), refineEdgeCache.end(), 0);
        size_t numEdgesToBeRefined = 0;
        ComputeEdgesRefinementMaskFromSamples(f, refineEdgeCache, numEdgesToBeRefined);

        m_faceMask[f] = 0;
        if (numEdgesToBeRefined > 1)
        {
            m_faceMask[f] = 1;

            for (auto n = 0; n < m_mesh->GetNumFaceEdges(f); n++)
            {
                if (refineEdgeCache[n] == 1)
                {
                    const auto edgeIndex = m_mesh->m_facesEdges[f][n];
                    if (edgeIndex != sizetMissingValue)
                    {
                        m_edgeMask[edgeIndex] = 1;
                    }
                }
            }
        }
    }
};

void MeshRefinement::FindHangingNodes(size_t face,
                                      size_t& numHangingEdges,
                                      size_t& numHangingNodes,
                                      size_t& numEdgesToRefine)
{

    numEdgesToRefine = 0;
    numHangingEdges = 0;
    numHangingNodes = 0;
    const auto numFaceNodes = m_mesh->GetNumFaceEdges(face);

    if (numFaceNodes > maximumNumberOfEdgesPerNode)
    {
        throw AlgorithmError("MeshRefinement::FindHangingNodes: The number of face nodes is greater than the maximum number of edges per node.");
    }

    m_isHangingNodeCache.resize(maximumNumberOfNodesPerFace);
    m_isHangingEdgeCache.resize(maximumNumberOfEdgesPerFace);
    std::fill(m_isHangingNodeCache.begin(), m_isHangingNodeCache.end(), false);
    std::fill(m_isHangingEdgeCache.begin(), m_isHangingEdgeCache.end(), false);

    auto kknod = numFaceNodes;
    for (auto n = 0; n < numFaceNodes; n++)
    {
        const auto edgeIndex = m_mesh->m_facesEdges[face][n];
        if (m_edgeMask[edgeIndex] != 0)
        {
            numEdgesToRefine += 1;
        }

        // check if the parent edge is in the cell
        if (m_brotherEdges[edgeIndex] != sizetMissingValue)
        {
            const auto e = NextCircularBackwardIndex(n, numFaceNodes);
            const auto ee = NextCircularForwardIndex(n, numFaceNodes);
            const auto firstEdgeIndex = m_mesh->m_facesEdges[face][e];
            const auto secondEdgeIndex = m_mesh->m_facesEdges[face][ee];

            size_t commonNode = sizetMissingValue;
            if (m_brotherEdges[edgeIndex] == firstEdgeIndex)
            {
                commonNode = m_mesh->FindCommonNode(edgeIndex, firstEdgeIndex);
                if (commonNode == sizetMissingValue)
                {
                    throw AlgorithmError("MeshRefinement::FindHangingNodes: Could not find common node.");
                }
            }
            else if (m_brotherEdges[edgeIndex] == secondEdgeIndex)
            {
                commonNode = m_mesh->FindCommonNode(edgeIndex, secondEdgeIndex);
                if (commonNode == sizetMissingValue)
                {
                    throw AlgorithmError("MeshRefinement::FindHangingNodes: Could not find common node.");
                }
            }

            if (commonNode != sizetMissingValue)
            {
                m_isHangingEdgeCache[n] = true;
                numHangingEdges++;
                for (auto nn = 0; nn < numFaceNodes; nn++)
                {
                    kknod = NextCircularForwardIndex(kknod, numFaceNodes);

                    if (m_mesh->m_facesNodes[face][kknod] == commonNode && !m_isHangingNodeCache[kknod])
                    {
                        numHangingNodes++;
                        m_isHangingNodeCache[kknod] = true;
                        break;
                    }
                }
            }
        }
    }
}

void MeshRefinement::ComputeEdgesRefinementMaskFromSamples(size_t face,
                                                           std::vector<size_t>& refineEdgeCache,
                                                           size_t& numEdgesToBeRefined)
{
    numEdgesToBeRefined = 0;
    if (m_refinementType == RefinementType::RidgeRefinement)
    {
        throw AlgorithmError("MeshRefinement::ComputeEdgesRefinementMaskFromSamples: This functionality is not implemented yet.");
    }

    const auto refinementValue = m_averaging->GetResults()[face];

    if (m_refinementType == RefinementType::RefinementLevels && refinementValue <= 0)
    {
        return;
    }

    if (IsEqual(refinementValue, doubleMissingValue))
    {
        return;
    }

    // compute all lengths
    const auto numEdges = m_mesh->GetNumFaceEdges(face);
    m_mesh->ComputeFaceClosedPolygonWithLocalMappings(face, m_polygonNodesCache, m_localNodeIndicesCache, m_globalEdgeIndicesCache);
    for (auto i = 0; i < numEdges; i++)
    {
        const auto edgeIndex = m_globalEdgeIndicesCache[i];
        if (m_mesh->m_edgeLengths[edgeIndex] < mergingDistance)
        {
            numEdgesToBeRefined++;
            continue;
        }

        bool doRefinement = false;

        // based on wave courant
        if (m_refinementType == RefinementType::WaveCourant)
        {
            const double newEdgeLength = 0.5 * m_mesh->m_edgeLengths[edgeIndex];
            const double c = std::sqrt(gravity * std::abs(refinementValue));
            const double waveCourant = c * (m_sampleRefineParameters.MinimumCellSize / std::sqrt(gravity)) / m_mesh->m_edgeLengths[edgeIndex];
            doRefinement = waveCourant < 1.0 && std::abs(newEdgeLength - m_sampleRefineParameters.MinimumCellSize) < std::abs(m_mesh->m_edgeLengths[edgeIndex] - m_sampleRefineParameters.MinimumCellSize);
        }

        // based on refinement levels
        if (m_refinementType == RefinementType::RefinementLevels && refinementValue > 0.0)
        {
            doRefinement = true;
        }

        if (doRefinement)
        {
            numEdgesToBeRefined++;
            refineEdgeCache[i] = 1;
        }
    }

    if (numEdgesToBeRefined > 0)
    {
        numEdgesToBeRefined = 0;
        for (auto i = 0; i < numEdges; i++)
        {
            if (refineEdgeCache[i] == 1 || m_isHangingNodeCache[i])
            {
                numEdgesToBeRefined++;
            }
        }
    }

    if (!m_directionalRefinement)
    {
        if (numEdgesToBeRefined == numEdges)
        {
            for (auto i = 0; i < numEdges; i++)
            {
                if (!m_isHangingNodeCache[i])
                {
                    refineEdgeCache[i] = 1;
                }
            }
        }
        else
        {
            numEdgesToBeRefined = 0;
        }
    }
}

void MeshRefinement::ComputeEdgesRefinementMask()
{
    bool repeat = true;
    size_t iter = 0;
    const size_t numMaxIterations = 6;
    std::vector<int> isQuadEdge(numNodesQuads);
    std::vector<size_t> numOfEdges(maximumNumberOfEdgesPerFace);

    while (repeat && iter < numMaxIterations)
    {
        repeat = false;
        iter++;

        for (auto f = 0; f < m_mesh->GetNumFaces(); f++)
        {
            if (m_faceMask[f] == 0)
            {
                continue;
            }

            size_t numHangingEdges;
            size_t numHangingNodes;
            size_t numEdgesToRefine;
            FindHangingNodes(f, numHangingEdges, numHangingNodes, numEdgesToRefine);

            const auto numFaceNodes = m_mesh->GetNumFaceEdges(f);

            // non-quads
            const auto numNodesEffective = numFaceNodes - numHangingNodes;
            if (numNodesEffective != numNodesQuads)
            {
                for (auto n = 0; n < numFaceNodes; n++)
                {
                    const auto e = NextCircularBackwardIndex(n, numFaceNodes);
                    const auto ee = NextCircularForwardIndex(n, numFaceNodes);

                    const auto edgeIndex = m_mesh->m_facesEdges[f][n];

                    const auto firstEdgeIndex = m_mesh->m_facesEdges[f][e];
                    const auto secondEdgeIndex = m_mesh->m_facesEdges[f][ee];

                    //do not refine edges with an hanging node
                    if (m_brotherEdges[edgeIndex] != firstEdgeIndex && m_brotherEdges[edgeIndex] != secondEdgeIndex)
                    {
                        m_edgeMask[edgeIndex] = 1;
                    }
                }
            }
            if (numNodesEffective == numNodesQuads)
            {
                // number the links in the cell, links that share a hanging node will have the same number
                size_t num = 0;
                for (auto n = 0; n < numFaceNodes; n++)
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

                if (num + 1 != numNodesQuads)
                {
                    throw AlgorithmError("MeshRefinement::ComputeEdgesRefinementMask: The number the links in the cell is not equals 3.");
                }

                numEdgesToRefine = 0;
                size_t firstEdgeIndex = 0;
                size_t secondEdgeIndex = 0;
                for (auto i = 0; i < numNodesQuads; i++)
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

                for (auto n = 0; n < numFaceNodes; n++)
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
        throw AlgorithmError("MeshRefinement::ComputeEdgesRefinementMask: Solution did not converge.");
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
    const size_t maxiter = 1000;
    size_t num = 1;
    size_t iter = 0;
    while (num != 0)
    {
        iter++;
        if (iter > maxiter)
        {
            break;
        }

        num = 0;
        for (auto f = 0; f < m_mesh->GetNumFaces(); f++)
        {
            if (m_faceMask[f] != 0 && m_faceMask[f] != -1)
            {
                continue;
            }

            size_t numHangingEdges;
            size_t numHangingNodes;
            size_t numEdgesToRefine;
            FindHangingNodes(f, numHangingEdges, numHangingNodes, numEdgesToRefine);

            bool isSplittingRequired = false;

            // check if the edge has a brother edge and needs to be refined
            const auto numFaceNodes = m_mesh->GetNumFaceEdges(f);

            if (numFaceNodes > maximumNumberOfEdgesPerFace)
            {
                return;
            }

            for (auto n = 0; n < numFaceNodes; n++)
            {
                const auto edgeIndex = m_mesh->m_facesEdges[f][n];
                if (m_isHangingEdgeCache[n] && m_edgeMask[edgeIndex] > 0)
                {
                    isSplittingRequired = true;
                }
            }

            //compute the effective face type
            const auto numNodesEffective = numFaceNodes - static_cast<size_t>(static_cast<double>(numHangingEdges) / 2.0);
            if (2 * (numFaceNodes - numNodesEffective) != numHangingEdges)
            {
                //uneven number of brotherlinks
                // TODO: ADD DOT
            }

            if (numFaceNodes + numEdgesToRefine > maximumNumberOfEdgesPerFace || // would result in unsupported cells after refinement
                numFaceNodes - numHangingNodes - numEdgesToRefine <= 1 ||        // faces with only one unrefined edge
                numNodesEffective == numEdgesToRefine)                           // refine all edges
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

                for (auto n = 0; n < numFaceNodes; n++)
                {
                    const auto edgeIndex = m_mesh->m_facesEdges[f][n];
                    if (!m_isHangingEdgeCache[n] && m_edgeMask[edgeIndex] == 0)
                    {
                        m_edgeMask[edgeIndex] = 1;
                        num++;
                    }
                    if (iter == maxiter)
                    {
                        //TODO: ADD DOT/MESSAGES
                    }
                }
            }
        }
    }
}

void MeshRefinement::FindBrotherEdges()
{
    m_brotherEdges.resize(m_mesh->GetNumEdges());
    std::fill(m_brotherEdges.begin(), m_brotherEdges.end(), sizetMissingValue);

    for (auto n = 0; n < m_mesh->GetNumNodes(); n++)
    {
        const auto numEdgesNodes = m_mesh->m_nodesNumEdges[n];
        for (auto e = 0; e < numEdgesNodes; e++)
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

            //check if node k is in the middle
            const auto firstEdgeOtherNode = OtherNodeOfEdge(m_mesh->m_edges[firstEdgeIndex], n);
            const auto secondEdgeOtherNode = OtherNodeOfEdge(m_mesh->m_edges[secondEdgeIndex], n);

            //compute tolerance
            const auto firstEdgeSquaredLength = ComputeSquaredDistance(m_mesh->m_nodes[firstEdgeOtherNode], m_mesh->m_nodes[n], m_mesh->m_projection);
            const auto secondEdgeSquaredLength = ComputeSquaredDistance(m_mesh->m_nodes[secondEdgeOtherNode], m_mesh->m_nodes[n], m_mesh->m_projection);
            const auto squaredTolerance = 0.0000001 * std::max(firstEdgeSquaredLength, secondEdgeSquaredLength);

            //The center of the two edges coincides with the shared node
            const auto center = ComputeMiddlePointAccountingForPoles(m_mesh->m_nodes[firstEdgeOtherNode], m_mesh->m_nodes[secondEdgeOtherNode], m_mesh->m_projection);
            const auto squaredDistanceFromCentre = ComputeSquaredDistance(center, m_mesh->m_nodes[n], m_mesh->m_projection);
            if (squaredDistanceFromCentre < squaredTolerance)
            {
                m_brotherEdges[firstEdgeIndex] = secondEdgeIndex;
                m_brotherEdges[secondEdgeIndex] = firstEdgeIndex;
            }
        }
    }
}

void MeshRefinement::SmoothEdgeRefinementMask() const
{
    throw AlgorithmError("MeshRefinement::SmoothEdgeRefinementMask: Not implemented yet.");
}
