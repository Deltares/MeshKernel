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

#include <stdexcept>
#include <vector>

#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/LandBoundaries.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Polygons.hpp>

namespace meshkernel
{
    LandBoundaries::LandBoundaries(const std::vector<Point>& landBoundary,
                                   std::shared_ptr<Mesh2D> mesh,
                                   std::shared_ptr<Polygons> polygons) : m_mesh(mesh),
                                                                         m_polygons(polygons)
    {
        if (!landBoundary.empty())
        {
            m_nodes.reserve(10000);
            std::copy(landBoundary.begin(), landBoundary.end(), std::back_inserter(m_nodes));
            m_polygonNodesCache.resize(maximumNumberOfNodesPerFace);
        }
    }

    void LandBoundaries::Administrate()
    {
        if (m_nodes.empty())
        {
            return;
        }

        // Computes the land boundary nodes inside a polygon
        m_nodeFaceIndices = m_mesh->PointFaceIndices(m_nodes);

        // do not consider the landboundary nodes outside the polygon
        std::vector<bool> nodeMask(m_nodes.size() - 1, false);
        for (auto n = 0; n < m_nodes.size() - 1; n++)
        {
            if (!m_nodes[n].IsValid() || !m_nodes[n + 1].IsValid())
            {
                continue;
            }

            const bool firstPointInPolygon = m_polygons->IsPointInPolygon(m_nodes[n], 0);
            const bool secondPointInPolygon = m_polygons->IsPointInPolygon(m_nodes[n + 1], 0);

            if (firstPointInPolygon || secondPointInPolygon)
            {
                nodeMask[n] = true;
            }
        }

        // mesh boundary to polygon
        const std::vector<Point> polygonNodes;
        const auto meshBoundaryPolygon = m_mesh->MeshBoundaryToPolygon(polygonNodes);

        // mask all land boundary nodes close to the mesh boundary (distanceFromMeshNode < minDistance)
        for (auto n = 0; n < m_nodes.size() - 1; n++)
        {
            if (!nodeMask[n] || meshBoundaryPolygon.empty())
            {
                continue;
            }

            Point firstPoint = m_nodes[n];
            Point secondPoint = m_nodes[n + 1];
            bool landBoundaryIsClose = false;

            for (auto nn = 0; nn < meshBoundaryPolygon.size() - 2; nn++)
            {
                Point firstMeshBoundaryNode = meshBoundaryPolygon[nn];
                Point secondMeshBoundaryNode = meshBoundaryPolygon[nn + 1];

                if (!firstMeshBoundaryNode.IsValid() || !secondMeshBoundaryNode.IsValid())
                {
                    continue;
                }

                const double edgeLength = ComputeDistance(firstMeshBoundaryNode, secondMeshBoundaryNode, m_mesh->m_projection);
                const double minDistance = m_closeFactor * edgeLength;

                Point normalPoint;
                double ratio;
                const double distanceFromMeshNode = DistanceFromLine(firstMeshBoundaryNode, firstPoint, secondPoint, m_mesh->m_projection, normalPoint, ratio);

                if (distanceFromMeshNode <= minDistance)
                {
                    landBoundaryIsClose = true;
                    break;
                }
            }

            if (landBoundaryIsClose)
            {
                nodeMask[n] = true;
            }
        }

        // find the start/end node of the land boundaries.
        // Emplace back them in m_validLandBoundaries if the land-boundary segment is close to a mesh node
        // TODO:  Why is splitting in two segments required?
        const auto indices = FindIndices(m_nodes, 0, m_nodes.size(), doubleMissingValue);
        m_validLandBoundaries.reserve(indices.size());
        for (const auto& index : indices)
        {
            // if the start node has a valid mask
            if (nodeMask[index[0]])
            {
                m_validLandBoundaries.emplace_back(index);
            }
        }

        // Generate two segments for closed land boundaries
        const auto numSegmentIndicesBeforeSplitting = m_validLandBoundaries.size();
        for (auto i = 0; i < numSegmentIndicesBeforeSplitting; i++)
        {
            const auto startSegmentIndex = m_validLandBoundaries[i][0];
            const auto endSegmentIndex = m_validLandBoundaries[i][1];
            if (endSegmentIndex - startSegmentIndex > 1)
            {
                const auto split = static_cast<size_t>(startSegmentIndex + (endSegmentIndex - startSegmentIndex) / 2);
                m_validLandBoundaries[i][1] = split;
                m_validLandBoundaries.emplace_back(std::initializer_list<size_t>{split, endSegmentIndex});
            }
        }
    };

    void LandBoundaries::FindNearestMeshBoundary(ProjectToLandBoundaryOption projectToLandBoundaryOption)
    {
        if (m_nodes.empty() || projectToLandBoundaryOption == ProjectToLandBoundaryOption::DoNotProjectToLandBoundary)
        {
            return;
        }

        m_closeFactor = m_closeWholeMeshFactor;
        if (projectToLandBoundaryOption == ProjectToLandBoundaryOption::OuterMeshBoundaryToLandBoundary ||
            projectToLandBoundaryOption == ProjectToLandBoundaryOption::InnerAndOuterMeshBoundaryToLandBoundary)
        {
            m_findOnlyOuterMeshBoundary = true;
            m_closeFactor = m_closeToLandBoundaryFactor;
        }

        Administrate();

        m_nodeMask.resize(m_mesh->GetNumNodes(), sizetMissingValue);
        m_faceMask.resize(m_mesh->GetNumFaces(), false);
        m_edgeMask.resize(m_mesh->GetNumEdges(), sizetMissingValue);
        m_meshNodesLandBoundarySegments.resize(m_mesh->GetNumNodes(), sizetMissingValue);
        m_nodesMinDistances.resize(m_mesh->GetNumNodes(), doubleMissingValue);

        // Loop over the segments of the land boundary and assign each node to the land boundary segment index
        for (auto landBoundarySegment = 0; landBoundarySegment < m_validLandBoundaries.size(); landBoundarySegment++)
        {
            size_t numPaths = 0;
            size_t numRejectedPaths = 0;
            MakePath(landBoundarySegment, numPaths, numRejectedPaths);

            if (numRejectedPaths > 0 && projectToLandBoundaryOption == ProjectToLandBoundaryOption::InnerAndOuterMeshBoundaryToLandBoundary)
            {
                m_findOnlyOuterMeshBoundary = false;
                MakePath(landBoundarySegment, numPaths, numRejectedPaths);
                m_findOnlyOuterMeshBoundary = true;
            }
        }

        // connect the m_mesh nodes
        if (m_findOnlyOuterMeshBoundary)
        {
            std::vector<size_t> connectedNodes;
            for (auto e = 0; e < m_mesh->GetNumEdges(); e++)
            {
                if (!m_mesh->IsEdgeOnBoundary(e))
                {
                    continue;
                }

                AssignLandBoundaryPolylineToMeshNodes(e, true, connectedNodes, 0);
            }
        }
    };

    void LandBoundaries::AssignLandBoundaryPolylineToMeshNodes(size_t edgeIndex, bool initialize, std::vector<size_t>& nodes, size_t numNodes)
    {
        if (m_nodes.empty())
        {
            return;
        }

        std::vector<size_t> nodesLoc;
        size_t numNodesLoc;

        if (initialize)
        {
            if (!m_mesh->IsEdgeOnBoundary(edgeIndex) || m_mesh->m_edges[edgeIndex].first == sizetMissingValue || m_mesh->m_edges[edgeIndex].second == sizetMissingValue)
                throw std::invalid_argument("LandBoundaries::AssignLandBoundaryPolylineToMeshNodes: Cannot not assign segment to mesh nodes.");

            const auto firstMeshNode = m_mesh->m_edges[edgeIndex].first;
            const auto secondMeshNode = m_mesh->m_edges[edgeIndex].second;

            if (m_meshNodesLandBoundarySegments[firstMeshNode] != sizetMissingValue &&
                m_meshNodesLandBoundarySegments[secondMeshNode] == sizetMissingValue &&
                m_nodeMask[firstMeshNode] != sizetMissingValue &&
                m_nodeMask[secondMeshNode] != sizetMissingValue)
            {
                nodesLoc.resize(3);
                nodesLoc[0] = firstMeshNode;
                nodesLoc[1] = secondMeshNode;
                numNodesLoc = 2;
            }
            else if (m_meshNodesLandBoundarySegments[firstMeshNode] == sizetMissingValue &&
                     m_meshNodesLandBoundarySegments[secondMeshNode] != sizetMissingValue &&
                     m_nodeMask[firstMeshNode] != sizetMissingValue &&
                     m_nodeMask[secondMeshNode] != sizetMissingValue)
            {
                nodesLoc.resize(3);
                nodesLoc[0] = secondMeshNode;
                nodesLoc[1] = firstMeshNode;
                numNodesLoc = 2;
            }
            else
            {
                //not a valid edge
                return;
            }
        }
        else
        {
            nodesLoc.resize(numNodes + 1);
            numNodesLoc = numNodes;
            std::copy(nodes.begin(), nodes.end(), nodesLoc.begin());
        }

        const auto maxNodes = *std::max_element(nodesLoc.begin(), nodesLoc.end() - 1);
        if (numNodesLoc > maxNodes)
        {
            return;
        }

        const auto lastVisitedNode = nodesLoc[numNodesLoc - 1];

        for (auto e = 0; e < m_mesh->m_nodesNumEdges[lastVisitedNode]; e++)
        {
            const auto edge = m_mesh->m_nodesEdges[lastVisitedNode][e];

            if (!m_mesh->IsEdgeOnBoundary(edge))
                continue;

            const auto otherNode = OtherNodeOfEdge(m_mesh->m_edges[edge], lastVisitedNode);

            // path stopped
            if (m_nodeMask[otherNode] == sizetMissingValue)
                break;

            // TODO: C++ 20 for(auto& i :  views::reverse(vec))
            bool otherNodeAlreadyVisited = false;
            for (auto n = numNodesLoc - 1; n < numNodesLoc; --n)
            {
                if (otherNode == nodesLoc[n])
                {
                    otherNodeAlreadyVisited = true;
                    break;
                }
            }

            if (otherNodeAlreadyVisited)
                continue;

            nodesLoc[numNodesLoc] = otherNode;

            if (m_meshNodesLandBoundarySegments[otherNode] != sizetMissingValue)
            {
                // Now check landboundary for otherNode
                for (auto n = 1; n < numNodesLoc; n++)
                {
                    const auto meshNode = nodesLoc[n];

                    double minimumDistance;
                    Point pointOnLandBoundary;
                    size_t nearestLandBoundaryNodeIndex = sizetMissingValue;
                    double edgeRatio;
                    NearestLandBoundaryNode(m_mesh->m_projection, m_mesh->m_nodes[meshNode], 0, m_nodes.size(), minimumDistance, pointOnLandBoundary, nearestLandBoundaryNodeIndex, edgeRatio);

                    // find the segment index of the found point
                    size_t landboundarySegmentIndex = std::numeric_limits<size_t>::max();
                    for (auto s = 0; s < m_validLandBoundaries.size(); s++)
                    {
                        if (nearestLandBoundaryNodeIndex >= m_validLandBoundaries[s][0] && nearestLandBoundaryNodeIndex < m_validLandBoundaries[s][1])
                        {
                            landboundarySegmentIndex = s;
                            break;
                        }
                    }

                    if (landboundarySegmentIndex == std::numeric_limits<size_t>::max())
                    {
                        throw AlgorithmError("LandBoundaries::AssignLandBoundaryPolylineToMeshNodes: No segment index found: cannot assign segment to mesh nodes.");
                    }

                    if ((nearestLandBoundaryNodeIndex == m_validLandBoundaries[landboundarySegmentIndex][0] && edgeRatio < 0.0) ||
                        (nearestLandBoundaryNodeIndex == m_validLandBoundaries[landboundarySegmentIndex][1] - 1 && edgeRatio > 1.0))
                    {
                        if (m_addLandboundaries)
                        {
                            AddLandBoundary(nodesLoc, numNodesLoc, lastVisitedNode);
                            m_meshNodesLandBoundarySegments[meshNode] = m_validLandBoundaries.size() - 1; //last added ;and boundary
                        }
                    }
                    else
                    {
                        m_meshNodesLandBoundarySegments[meshNode] = landboundarySegmentIndex;
                    }
                }
            }
            else
            {
                AssignLandBoundaryPolylineToMeshNodes(edge, false, nodesLoc, numNodesLoc + 1);
            }
        }
    }

    void LandBoundaries::AddLandBoundary(const std::vector<size_t>& nodesLoc, size_t numNodesLoc, size_t nodeIndex)
    {
        if (m_nodes.empty())
        {
            return;
        }

        const auto startSegmentIndex = m_meshNodesLandBoundarySegments[nodesLoc[0]];
        const auto endSegmentIndex = m_meshNodesLandBoundarySegments[nodesLoc[numNodesLoc]];

        if (startSegmentIndex == sizetMissingValue || startSegmentIndex >= m_validLandBoundaries.size() ||
            endSegmentIndex == sizetMissingValue || endSegmentIndex >= m_validLandBoundaries.size())
        {
            throw std::invalid_argument("LandBoundaries::AddLandBoundary: Invalid segment index.");
        }

        // find start/end
        auto startNode = m_validLandBoundaries[startSegmentIndex][0];
        auto endNode = m_validLandBoundaries[startSegmentIndex][1];

        Point newNodeLeft;
        if (ComputeSquaredDistance(m_mesh->m_nodes[nodeIndex], m_nodes[startNode], m_mesh->m_projection) <= ComputeSquaredDistance(m_mesh->m_nodes[nodeIndex], m_nodes[endNode], m_mesh->m_projection))
        {
            newNodeLeft = m_nodes[startNode];
        }
        else
        {
            newNodeLeft = m_nodes[endNode];
        }

        Point newNodeRight;
        if (endSegmentIndex == startSegmentIndex)
        {
            newNodeRight = m_nodes[startNode] + m_nodes[endNode] - newNodeLeft;
        }
        else
        {
            // find start/end
            startNode = m_validLandBoundaries[endSegmentIndex][0];
            endNode = m_validLandBoundaries[endSegmentIndex][1];
            if (ComputeSquaredDistance(m_mesh->m_nodes[nodeIndex], m_nodes[startNode], m_mesh->m_projection) <= ComputeSquaredDistance(m_mesh->m_nodes[nodeIndex], m_nodes[endNode], m_mesh->m_projection))
            {
                newNodeRight = m_nodes[startNode];
            }
            else
            {
                newNodeRight = m_nodes[endNode];
            }
        }

        // Update  nodes
        m_nodes.emplace_back(Point{doubleMissingValue, doubleMissingValue});
        m_nodes.emplace_back(newNodeLeft);
        m_nodes.emplace_back(newNodeRight);
        m_nodes.emplace_back(Point{doubleMissingValue, doubleMissingValue});

        // Update segment indices
        m_validLandBoundaries.push_back(std::initializer_list<size_t>{m_nodes.size() - 3, m_nodes.size() - 2});
    }

    void LandBoundaries::MakePath(size_t landBoundaryPolyline,
                                  size_t& numNodesInPath,
                                  size_t& numRejectedNodesInPath)
    {
        if (m_nodes.empty())
        {
            return;
        }

        const auto startLandBoundaryIndex = m_validLandBoundaries[landBoundaryPolyline][0];
        const auto endLandBoundaryIndex = m_validLandBoundaries[landBoundaryPolyline][1];

        if (startLandBoundaryIndex >= m_nodes.size() || startLandBoundaryIndex >= endLandBoundaryIndex)
        {
            throw std::invalid_argument("LandBoundaries::MakePath: Invalid boundary index.");
        }

        // fractional location of the projected outer nodes(min and max) on the land boundary segment
        ComputeMeshNodeMask(landBoundaryPolyline);

        size_t startMeshNode = sizetMissingValue;
        size_t endMeshNode = sizetMissingValue;
        FindStartEndMeshNodes(landBoundaryPolyline, startMeshNode, endMeshNode);

        if (startMeshNode == sizetMissingValue || endMeshNode == sizetMissingValue || startMeshNode == endMeshNode)
        {
            throw AlgorithmError("LandBoundaries::MakePath: Cannot not find valid mesh nodes.");
        }

        std::vector<size_t> connectedNodes;
        ShortestPath(landBoundaryPolyline, startLandBoundaryIndex, endLandBoundaryIndex, startMeshNode, connectedNodes);

        auto lastSegment = m_meshNodesLandBoundarySegments[endMeshNode];
        size_t lastNode = sizetMissingValue;
        auto currentNode = endMeshNode;
        double minDinstanceFromLandBoundary;
        double distanceFromLandBoundary;
        Point nodeOnLandBoundary;
        size_t currentNodeLandBoundaryNodeIndex;
        double currentNodeEdgeRatio;
        bool stopPathSearch; //the path has been temporarily stopped(true) or not (false)
        size_t numConnectedNodes = 0;
        numRejectedNodesInPath = 0;
        numNodesInPath = 0;

        while (true)
        {
            stopPathSearch = true;

            if (m_meshNodesLandBoundarySegments[currentNode] != sizetMissingValue)
            {
                // Multiple boundary segments: take the nearest
                const auto previousLandBoundarySegment = m_meshNodesLandBoundarySegments[currentNode];
                const auto previousStartMeshNode = m_validLandBoundaries[previousLandBoundarySegment][0];
                const auto previousEndMeshNode = m_validLandBoundaries[previousLandBoundarySegment][1];
                double previousMinDistance;

                NearestLandBoundaryNode(m_mesh->m_projection, m_mesh->m_nodes[currentNode], previousStartMeshNode, previousEndMeshNode,
                                        previousMinDistance, nodeOnLandBoundary, currentNodeLandBoundaryNodeIndex, currentNodeEdgeRatio);

                NearestLandBoundaryNode(m_mesh->m_projection, m_mesh->m_nodes[currentNode], startLandBoundaryIndex, endLandBoundaryIndex,
                                        distanceFromLandBoundary, nodeOnLandBoundary, currentNodeLandBoundaryNodeIndex, currentNodeEdgeRatio);

                const auto minDinstanceFromLandBoundaryCurrentNode = m_nodesMinDistances[currentNode];

                if (distanceFromLandBoundary <= previousMinDistance &&
                    distanceFromLandBoundary < m_minDistanceFromLandFactor * minDinstanceFromLandBoundaryCurrentNode)
                {
                    stopPathSearch = false;
                }
            }
            else
            {
                if (IsEqual(m_nodesMinDistances[currentNode], doubleMissingValue))
                {
                    NearestLandBoundaryNode(m_mesh->m_projection, m_mesh->m_nodes[currentNode], 0, m_nodes.size(),
                                            minDinstanceFromLandBoundary, nodeOnLandBoundary, currentNodeLandBoundaryNodeIndex, currentNodeEdgeRatio);
                    m_nodesMinDistances[currentNode] = minDinstanceFromLandBoundary;
                }
                else
                {
                    minDinstanceFromLandBoundary = m_nodesMinDistances[currentNode];
                }

                NearestLandBoundaryNode(m_mesh->m_projection, m_mesh->m_nodes[currentNode], startLandBoundaryIndex, endLandBoundaryIndex,
                                        distanceFromLandBoundary, nodeOnLandBoundary, currentNodeLandBoundaryNodeIndex, currentNodeEdgeRatio);

                if (distanceFromLandBoundary < m_minDistanceFromLandFactor * minDinstanceFromLandBoundary &&
                    (!m_findOnlyOuterMeshBoundary || m_mesh->m_nodesTypes[currentNode] == 2 || m_mesh->m_nodesTypes[currentNode] == 3))
                {
                    stopPathSearch = false;
                }
            }

            if (stopPathSearch)
            {
                if (numConnectedNodes == 1 && lastSegment != sizetMissingValue)
                {
                    m_meshNodesLandBoundarySegments[lastNode] = lastSegment;
                }
                numConnectedNodes = 0;
                numRejectedNodesInPath += 1;
            }
            else
            {
                lastSegment = m_meshNodesLandBoundarySegments[currentNode];
                lastNode = currentNode;

                numNodesInPath += 1;
                numConnectedNodes += 1;

                m_meshNodesLandBoundarySegments[lastNode] = landBoundaryPolyline;
            }

            if (currentNode == startMeshNode)
            {
                break;
            }

            const auto nextEdgeIndex = connectedNodes[currentNode];
            if (nextEdgeIndex == sizetMissingValue || nextEdgeIndex >= m_mesh->GetNumEdges())
            {
                break;
            }

            currentNode = OtherNodeOfEdge(m_mesh->m_edges[nextEdgeIndex], currentNode);

            if (currentNode == sizetMissingValue || currentNode >= m_mesh->GetNumNodes())
            {
                break;
            }
        }

        if (numConnectedNodes == 1)
        {
            m_meshNodesLandBoundarySegments[lastNode] = lastSegment;
        }
    }

    void LandBoundaries::ComputeMeshNodeMask(size_t landBoundaryPolyline)
    {
        if (m_nodes.empty())
        {
            return;
        }

        const auto startLandBoundaryIndex = m_validLandBoundaries[landBoundaryPolyline][0];
        const auto endLandBoundaryIndex = m_validLandBoundaries[landBoundaryPolyline][1];

        // Try to find a face crossed by the current land boundary polyline:
        // 1. One of the land boundary nodes is inside a face
        // 2. Or one of the land boundary segments is crossing a boundary mesh edge
        size_t crossedFaceIndex = sizetMissingValue;
        for (auto i = startLandBoundaryIndex; i < endLandBoundaryIndex; i++)
        {
            crossedFaceIndex = m_nodeFaceIndices[i];
            if (crossedFaceIndex != sizetMissingValue)
            {
                break;
            }

            auto [face, edge] = m_mesh->IsSegmentCrossingABoundaryEdge(m_nodes[i], m_nodes[i + 1]);
            crossedFaceIndex = face;
            if (crossedFaceIndex != sizetMissingValue)
            {
                break;
            }
        }

        std::fill(m_nodeMask.begin(), m_nodeMask.end(), sizetMissingValue);
        if (m_landMask)
        {
            std::fill(m_faceMask.begin(), m_faceMask.end(), false);
            std::fill(m_edgeMask.begin(), m_edgeMask.end(), sizetMissingValue);
            //m_faceMask assumes crossedFace has already been done.
            if (crossedFaceIndex != sizetMissingValue && crossedFaceIndex < m_mesh->GetNumFaces())
            {
                m_faceMask[crossedFaceIndex] = true;
            }

            std::vector<size_t> landBoundaryFaces{crossedFaceIndex};
            MaskMeshFaceMask(landBoundaryFaces,
                             landBoundaryPolyline);

            // Mask all nodes of the masked faces
            for (auto f = 0; f < m_mesh->GetNumFaces(); f++)
            {
                if (m_faceMask[f])
                {
                    for (auto n = 0; n < m_mesh->GetNumFaceEdges(f); n++)
                    {
                        m_nodeMask[m_mesh->m_facesNodes[f][n]] = landBoundaryPolyline;
                    }
                }
            }
        }
        else
        {
            for (auto& e : m_nodeMask)
            {
                e = landBoundaryPolyline;
            }
        }

        for (auto n = 0; n < m_mesh->GetNumNodes(); n++)
        {
            if (m_nodeMask[n] != sizetMissingValue)
            {
                const bool inPolygon = m_polygons->IsPointInPolygon(m_mesh->m_nodes[n], 0);
                if (!inPolygon)
                {
                    m_nodeMask[n] = sizetMissingValue;
                }
            }
        }
    }

    void LandBoundaries::MaskMeshFaceMask(std::vector<size_t>& initialFaces,
                                          size_t landBoundaryPolyline)
    {
        if (m_nodes.empty())
        {
            return;
        }

        std::vector<size_t> nextFaces;
        nextFaces.reserve(initialFaces.size());
        for (const auto& face : initialFaces)
        {
            // no face was crossed by the land boundary: mask boundary faces only
            // these are the faces that are close (up to a certain tolerance) by a land boundary
            if (face == sizetMissingValue)
            {
                for (auto e = 0; e < m_mesh->GetNumEdges(); e++)
                {
                    // only boundary edges are considered
                    if (!m_mesh->IsEdgeOnBoundary(e))
                    {
                        continue;
                    }

                    const auto currentFace = m_mesh->m_edgesFaces[e][0];
                    // already masked
                    if (m_faceMask[currentFace])
                    {
                        continue;
                    }

                    for (const auto& edge : m_mesh->m_facesEdges[currentFace])
                    {
                        const auto landBoundaryNode = IsMeshEdgeCloseToLandBoundaries(edge, landBoundaryPolyline);
                        if (landBoundaryNode != sizetMissingValue)
                        {
                            m_faceMask[currentFace] = true;
                            break;
                        }
                    }
                }
            }
            else
            {
                // face is crossed
                if (m_mesh->GetNumFaces() < numNodesInTriangle)
                {
                    continue;
                }

                for (const auto& currentEdge : m_mesh->m_facesEdges[face])
                {
                    // is a boundary edge, continue
                    if (m_mesh->IsEdgeOnBoundary(currentEdge))
                    {
                        continue;
                    }

                    const auto otherFace = face == m_mesh->m_edgesFaces[currentEdge][0] ? m_mesh->m_edgesFaces[currentEdge][1] : m_mesh->m_edgesFaces[currentEdge][0];

                    // already masked
                    if (m_faceMask[otherFace])
                    {
                        continue;
                    }

                    bool isFaceFound = false;
                    for (const auto& edge : m_mesh->m_facesEdges[otherFace])
                    {
                        if (m_edgeMask[edge] == 1)
                        {
                            // previously visited crossed edge
                            isFaceFound = true;
                            continue;
                        }
                        if (m_edgeMask[edge] == 0)
                        {
                            // previously visited uncrossed edge
                            continue;
                        }

                        // visited edge
                        m_edgeMask[edge] = 0;
                        const auto landBoundaryNode = IsMeshEdgeCloseToLandBoundaries(edge, landBoundaryPolyline);

                        if (landBoundaryNode != sizetMissingValue)
                        {
                            m_edgeMask[edge] = 1;
                            isFaceFound = true;
                        }
                    }

                    m_faceMask[otherFace] = isFaceFound;
                    if (m_faceMask[otherFace])
                    {
                        nextFaces.emplace_back(otherFace);
                    }
                }
            }
        }

        if (!nextFaces.empty())
        {
            MaskMeshFaceMask(nextFaces, landBoundaryPolyline);
        }
    }

    size_t LandBoundaries::IsMeshEdgeCloseToLandBoundaries(size_t meshEdgeIndex,
                                                           size_t landBoundaryPolyline)
    {
        if (m_nodes.empty())
        {
            return sizetMissingValue;
        }

        const auto startLandBoundaryIndex = m_validLandBoundaries[landBoundaryPolyline][0];
        const auto endLandBoundaryIndex = m_validLandBoundaries[landBoundaryPolyline][1];

        size_t landBoundaryNode = sizetMissingValue;
        const auto startNode = std::max(std::min(static_cast<size_t>(0), endLandBoundaryIndex - 1), startLandBoundaryIndex);

        if (m_mesh->m_edges[meshEdgeIndex].first == sizetMissingValue || m_mesh->m_edges[meshEdgeIndex].second == sizetMissingValue)
            return false;

        const auto firstMeshNode = m_mesh->m_nodes[m_mesh->m_edges[meshEdgeIndex].first];
        const auto secondMeshNode = m_mesh->m_nodes[m_mesh->m_edges[meshEdgeIndex].second];

        const double meshEdgeLength = ComputeDistance(firstMeshNode, secondMeshNode, m_mesh->m_projection);

        const double distanceFactor = m_findOnlyOuterMeshBoundary ? m_closeToLandBoundaryFactor : m_closeWholeMeshFactor;
        const double closeDistance = meshEdgeLength * distanceFactor;

        auto currentNode = startNode;
        size_t searchIter = 0;
        int stepNode = 0;
        while (searchIter < 3)
        {
            const double landBoundaryLength = ComputeSquaredDistance(m_nodes[currentNode], m_nodes[currentNode + 1], m_mesh->m_projection);

            if (landBoundaryLength > 0)
            {
                double ratioFirstMeshNode = 0.0;
                double ratioSecondMeshNode = 1.0;
                Point normalPoint;
                const double distanceFromLandBoundaryFirstMeshNode = DistanceFromLine(firstMeshNode,
                                                                                      m_nodes[currentNode],
                                                                                      m_nodes[currentNode + 1],
                                                                                      m_mesh->m_projection,
                                                                                      normalPoint,
                                                                                      ratioFirstMeshNode);

                if (distanceFromLandBoundaryFirstMeshNode < closeDistance)
                {
                    landBoundaryNode = currentNode;
                    // the projection of firstMeshNode is within the segment currentNode / currentNode + 1
                    if (ratioFirstMeshNode >= 0.0 && ratioFirstMeshNode <= 1.0)
                    {
                        break;
                    }
                }
                else
                {
                    // check the second point
                    Point normalPoint;
                    const auto distanceFromLandBoundarySecondMeshNode = DistanceFromLine(secondMeshNode,
                                                                                         m_nodes[currentNode],
                                                                                         m_nodes[currentNode + 1],
                                                                                         m_mesh->m_projection,
                                                                                         normalPoint,
                                                                                         ratioSecondMeshNode);

                    if (distanceFromLandBoundarySecondMeshNode < closeDistance)
                    {
                        landBoundaryNode = currentNode;
                        // the projection of secondMeshNode is within the segment currentNode / currentNode + 1
                        if (ratioSecondMeshNode >= 0.0 && ratioSecondMeshNode <= 1.0)
                            break;
                    }
                }
            }

            // search the next land boundary edge if projection is not within is within the segment currentNode / currentNode + 1
            searchIter = 0;
            while ((searchIter == 0 || currentNode < startLandBoundaryIndex || currentNode > endLandBoundaryIndex - 1) && searchIter < 3)
            {
                searchIter += 1;
                if (stepNode < 0)
                {
                    stepNode = -stepNode + 1;
                }
                else
                {
                    stepNode = -stepNode - 1;
                }
                currentNode = currentNode + stepNode;
            }
        }

        return landBoundaryNode;
    }

    void LandBoundaries::FindStartEndMeshNodes(size_t landBoundaryPolyline,
                                               size_t& startMeshNode,
                                               size_t& endMeshNode)
    {
        if (m_nodes.empty())
        {
            return;
        }

        const auto startLandBoundaryIndex = m_validLandBoundaries[landBoundaryPolyline][0];
        const auto endLandBoundaryIndex = m_validLandBoundaries[landBoundaryPolyline][1];
        const auto leftIndex = endLandBoundaryIndex - 1;
        const auto rightIndex = startLandBoundaryIndex;

        // compute the start and end point of the land boundary respectively
        const auto nextLeftIndex = std::min(leftIndex + 1, endLandBoundaryIndex);
        Point startPoint = {m_nodes[nextLeftIndex].x, m_nodes[nextLeftIndex].y};
        Point endPoint = {m_nodes[rightIndex].x, m_nodes[rightIndex].y};

        const bool isStartPointInsideAPolygon = m_polygons->IsPointInPolygon(startPoint, 0);
        const bool isEndPointInsideAPolygon = m_polygons->IsPointInPolygon(endPoint, 0);

        if (!isStartPointInsideAPolygon)
        {
            startPoint.x = m_nodes[leftIndex + 1].x;
            startPoint.y = m_nodes[leftIndex + 1].y;
        }
        if (!isEndPointInsideAPolygon)
        {
            endPoint.x = m_nodes[rightIndex].x;
            endPoint.y = m_nodes[rightIndex].y;
        }

        // Get the edges that are closest the land boundary
        Point normalPoint;
        double ratio;

        auto minDistStart = std::numeric_limits<double>::max();
        auto minDistEnd = std::numeric_limits<double>::max();
        size_t startEdge = sizetMissingValue;
        size_t endEdge = sizetMissingValue;

        for (auto e = 0; e < m_mesh->GetNumEdges(); e++)
        {
            if (m_mesh->m_edges[e].first == sizetMissingValue || m_mesh->m_edges[e].second == sizetMissingValue)
            {
                continue;
            }

            if (m_nodeMask[m_mesh->m_edges[e].first] == sizetMissingValue || m_nodeMask[m_mesh->m_edges[e].second] == sizetMissingValue)
            {
                continue;
            }

            const double distanceFromFirstMeshNode = DistanceFromLine(startPoint, m_mesh->m_nodes[m_mesh->m_edges[e].first], m_mesh->m_nodes[m_mesh->m_edges[e].second], m_mesh->m_projection, normalPoint, ratio);
            const double distanceFromSecondMeshNode = DistanceFromLine(endPoint, m_mesh->m_nodes[m_mesh->m_edges[e].first], m_mesh->m_nodes[m_mesh->m_edges[e].second], m_mesh->m_projection, normalPoint, ratio);

            if (distanceFromFirstMeshNode < minDistStart)
            {
                startEdge = e;
                minDistStart = distanceFromFirstMeshNode;
            }
            if (distanceFromSecondMeshNode < minDistEnd)
            {
                endEdge = e;
                minDistEnd = distanceFromSecondMeshNode;
            }
        }

        if (startEdge == sizetMissingValue || endEdge == sizetMissingValue)
        {
            throw std::invalid_argument("LandBoundaries::FindStartEndMeshNodes: Cannot find startMeshNode or endMeshNode.");
        }
        FindStartEndMeshNodesFromEdges(startEdge, endEdge, startPoint, endPoint, startMeshNode, endMeshNode);
    }

    void LandBoundaries::FindStartEndMeshNodesFromEdges(size_t startEdge,
                                                        size_t endEdge,
                                                        Point startPoint,
                                                        Point endPoint,
                                                        size_t& startMeshNode,
                                                        size_t& endMeshNode) const
    {
        if (m_nodes.empty())
        {
            return;
        }

        auto firstMeshNodeIndex = m_mesh->m_edges[startEdge].first;
        auto secondMeshNodeIndex = m_mesh->m_edges[startEdge].second;
        auto firstDistance = ComputeSquaredDistance(m_mesh->m_nodes[firstMeshNodeIndex], startPoint, m_mesh->m_projection);
        auto secondDistance = ComputeSquaredDistance(m_mesh->m_nodes[secondMeshNodeIndex], startPoint, m_mesh->m_projection);

        if (firstDistance <= secondDistance)
        {
            startMeshNode = firstMeshNodeIndex;
        }
        else
        {
            startMeshNode = secondMeshNodeIndex;
        }

        firstMeshNodeIndex = m_mesh->m_edges[endEdge].first;
        secondMeshNodeIndex = m_mesh->m_edges[endEdge].second;
        firstDistance = ComputeSquaredDistance(m_mesh->m_nodes[firstMeshNodeIndex], endPoint, m_mesh->m_projection);
        secondDistance = ComputeSquaredDistance(m_mesh->m_nodes[secondMeshNodeIndex], endPoint, m_mesh->m_projection);

        if (firstDistance <= secondDistance)
        {
            endMeshNode = firstMeshNodeIndex;
        }
        else
        {
            endMeshNode = secondMeshNodeIndex;
        }
    }

    void LandBoundaries::ShortestPath(size_t landBoundaryPolyline,
                                      size_t startLandBoundaryIndex,
                                      size_t endLandBoundaryIndex,
                                      size_t startMeshNode,
                                      std::vector<size_t>& connectedNodes)
    {
        if (m_nodes.empty())
        {
            return;
        }

        connectedNodes.resize(m_mesh->GetNumNodes(), sizetMissingValue);
        std::fill(connectedNodes.begin(), connectedNodes.end(), sizetMissingValue);

        // infinite distance for all nodes
        std::vector<double> nodeDistances(m_mesh->GetNumNodes(), std::numeric_limits<double>::max());
        std::vector<bool> isVisited(m_mesh->GetNumNodes(), false);

        auto currentNodeIndex = startMeshNode;
        nodeDistances[startMeshNode] = 0.0;
        while (true)
        {
            isVisited[currentNodeIndex] = true;
            Point currentNode = m_mesh->m_nodes[currentNodeIndex];
            Point currentNodeOnLandBoundary;
            size_t currentNodeLandBoundaryNodeIndex;
            double currentNodeEdgeRatio;
            double currentNodeDistance;
            NearestLandBoundaryNode(m_mesh->m_projection, currentNode, startLandBoundaryIndex, endLandBoundaryIndex,
                                    currentNodeDistance, currentNodeOnLandBoundary, currentNodeLandBoundaryNodeIndex, currentNodeEdgeRatio);

            if (currentNodeLandBoundaryNodeIndex == sizetMissingValue)
                throw AlgorithmError("LandBoundaries::ShortestPath: Cannot compute the nearest node on the land boundary.");

            for (const auto& edgeIndex : m_mesh->m_nodesEdges[currentNodeIndex])
            {
                if (m_mesh->m_edges[edgeIndex].first == sizetMissingValue || m_mesh->m_edges[edgeIndex].second == sizetMissingValue)
                {
                    continue;
                }

                const auto neighbouringNodeIndex = OtherNodeOfEdge(m_mesh->m_edges[edgeIndex], currentNodeIndex);

                if (isVisited[neighbouringNodeIndex])
                {
                    continue;
                }

                const auto neighbouringNode = m_mesh->m_nodes[neighbouringNodeIndex];
                Point neighbouringNodeOnLandBoundary;
                size_t neighbouringNodeLandBoundaryNodeIndex;
                double neighbouringNodeEdgeRatio;
                double neighbouringNodeDistance;
                NearestLandBoundaryNode(m_mesh->m_projection, neighbouringNode, startLandBoundaryIndex, endLandBoundaryIndex,
                                        neighbouringNodeDistance, neighbouringNodeOnLandBoundary, neighbouringNodeLandBoundaryNodeIndex, neighbouringNodeEdgeRatio);

                Point middlePoint{
                    (currentNode.x + neighbouringNode.x) * 0.5,
                    (currentNode.y + neighbouringNode.y) * 0.5};

                Point middlePointOnLandBoundary;
                size_t middlePointLandBoundaryNodeIndex;
                double middlePointEdgeRatio;
                double middlePointDistance;
                NearestLandBoundaryNode(m_mesh->m_projection, middlePoint, startLandBoundaryIndex, endLandBoundaryIndex,
                                        middlePointDistance, middlePointOnLandBoundary, middlePointLandBoundaryNodeIndex, middlePointEdgeRatio);

                double maximumDistance = std::max(currentNodeDistance, neighbouringNodeDistance);

                if (currentNodeLandBoundaryNodeIndex < neighbouringNodeLandBoundaryNodeIndex)
                {
                    for (auto n = currentNodeLandBoundaryNodeIndex + 1; n < neighbouringNodeLandBoundaryNodeIndex; ++n)
                    {
                        double ratio;
                        middlePointDistance = DistanceFromLine(m_nodes[n], currentNode, neighbouringNode, m_mesh->m_projection, middlePointOnLandBoundary, ratio);
                        if (middlePointDistance > maximumDistance)
                            maximumDistance = middlePointDistance;
                    }
                }
                else if (currentNodeLandBoundaryNodeIndex > neighbouringNodeLandBoundaryNodeIndex)
                {
                    for (auto n = neighbouringNodeLandBoundaryNodeIndex + 1; n < currentNodeLandBoundaryNodeIndex; ++n)
                    {
                        double ratio;
                        middlePointDistance = DistanceFromLine(m_nodes[n], currentNode, neighbouringNode, m_mesh->m_projection, middlePointOnLandBoundary, ratio);
                        if (middlePointDistance > maximumDistance)
                            maximumDistance = middlePointDistance;
                    }
                }

                // In case of netboundaries only: set penalty when edge is not on the boundary
                if (m_findOnlyOuterMeshBoundary && !m_mesh->IsEdgeOnBoundary(edgeIndex))
                    maximumDistance = 1e6 * maximumDistance;

                const double edgeLength = ComputeDistance(currentNode, neighbouringNode, m_mesh->m_projection);
                const double correctedDistance = nodeDistances[currentNodeIndex] + edgeLength * maximumDistance;

                if (correctedDistance < nodeDistances[neighbouringNodeIndex])
                {
                    nodeDistances[neighbouringNodeIndex] = correctedDistance;
                    connectedNodes[neighbouringNodeIndex] = edgeIndex;
                }
            }

            // linear search with masking
            currentNodeIndex = 0;
            double minValue = std::numeric_limits<double>::max();
            for (auto n = 0; n < m_mesh->GetNumNodes(); n++)
            {
                if (m_nodeMask[n] == landBoundaryPolyline && !isVisited[n] && nodeDistances[n] < minValue)
                {
                    currentNodeIndex = n;
                    minValue = nodeDistances[n];
                }
            }

            if (currentNodeIndex >= m_mesh->GetNumNodes() ||
                nodeDistances[currentNodeIndex] == std::numeric_limits<double>::max() ||
                isVisited[currentNodeIndex])
            {
                break;
            }
        }
    }

    void LandBoundaries::NearestLandBoundaryNode(const Projection& projection,
                                                 const Point& node,
                                                 size_t startLandBoundaryIndex,
                                                 size_t endLandBoundaryIndex,
                                                 double& minimumDistance,
                                                 Point& pointOnLandBoundary,
                                                 size_t& nearestLandBoundaryNodeIndex,
                                                 double& edgeRatio)
    {
        if (m_nodes.empty())
        {
            return;
        }

        minimumDistance = std::numeric_limits<double>::max();
        nearestLandBoundaryNodeIndex = sizetMissingValue;
        edgeRatio = -1.0;
        pointOnLandBoundary = node;
        for (auto n = startLandBoundaryIndex; n < endLandBoundaryIndex; n++)
        {
            if (!m_nodes[n].IsValid() || !m_nodes[n + 1].IsValid())
            {
                continue;
            }

            Point normalPoint{doubleMissingValue, doubleMissingValue};
            double ratio = 0.0;
            const double distanceFromLandBoundary = DistanceFromLine(node, m_nodes[n], m_nodes[n + 1], projection, normalPoint, ratio);

            if (distanceFromLandBoundary > 0.0 && distanceFromLandBoundary < minimumDistance)
            {
                minimumDistance = distanceFromLandBoundary;
                pointOnLandBoundary = normalPoint;
                nearestLandBoundaryNodeIndex = n;
                edgeRatio = ratio;
            }
        }
    }

    void LandBoundaries::SnapMeshToLandBoundaries()
    {
        if (m_nodes.empty() || m_meshNodesLandBoundarySegments.empty())
        {
            return;
        }

        const auto numNodes = m_mesh->GetNumNodes();
        for (auto n = 0; n < numNodes; ++n)
        {
            if (m_mesh->m_nodesTypes[n] == 1 || m_mesh->m_nodesTypes[n] == 2 || m_mesh->m_nodesTypes[n] == 3)
            {
                const auto meshNodeToLandBoundarySegment = m_meshNodesLandBoundarySegments[n];
                if (meshNodeToLandBoundarySegment == sizetMissingValue)
                {
                    continue;
                }

                double minimumDistance;
                Point pointOnLandBoundary;
                size_t nearestLandBoundaryNodeIndex;
                double edgeRatio;

                NearestLandBoundaryNode(m_mesh->m_projection,
                                        m_mesh->m_nodes[n],
                                        m_validLandBoundaries[meshNodeToLandBoundarySegment][0],
                                        m_validLandBoundaries[meshNodeToLandBoundarySegment][1],
                                        minimumDistance,
                                        pointOnLandBoundary,
                                        nearestLandBoundaryNodeIndex,
                                        edgeRatio);

                m_mesh->m_nodes[n] = pointOnLandBoundary;
            }
        }
    }
}; // namespace meshkernel
