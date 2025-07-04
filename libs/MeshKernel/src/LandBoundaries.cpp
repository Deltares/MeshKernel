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

#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/LandBoundaries.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/UndoActions/CompoundUndoAction.hpp>

using meshkernel::LandBoundaries;
using meshkernel::Mesh2D;

LandBoundaries::LandBoundaries(const std::vector<Point>& boundaryNodes,
                               Mesh2D& mesh) : m_mesh(mesh),
                                               m_polygons(std::vector<Point>(), mesh.m_projection),
                                               m_landBoundary(boundaryNodes)
{
    m_polygonNodesCache.clear();
}

LandBoundaries::LandBoundaries(const std::vector<Point>& boundaryNodes,
                               Mesh2D& mesh,
                               const Polygons& polygons) : m_mesh(mesh),
                                                           m_polygons(polygons),
                                                           m_landBoundary(boundaryNodes)
{
    if (!m_landBoundary.IsEmpty())
    {
        m_polygonNodesCache.resize(constants::geometric::maximumNumberOfNodesPerFace);
    }
}

void LandBoundaries::Administrate()
{
    if (m_landBoundary.IsEmpty())
    {
        return;
    }

    // Computes the land boundary nodes inside a polygon
    m_nodeFaceIndices = m_mesh.PointFaceIndices(m_landBoundary.GetNodes());

    // Do not consider the landboundary nodes outside the polygon
    std::vector<bool> nodeMask(m_landBoundary.GetNodeMask(m_polygons));

    // Mesh boundary to polygon
    const std::vector<Point> polygonNodes;
    const auto meshBoundaryPolygon = m_mesh.ComputeBoundaryPolygons(polygonNodes);

    // Find the start/end node of the land boundaries
    const auto landBoundaryIndices = m_landBoundary.FindPolylineIndices();

    // Emplace start/end node of the land boundaries back in m_validLandBoundaries,
    // if the land boundary segment is close to a mesh node
    m_validLandBoundaries.reserve(landBoundaryIndices.size());
    for (const auto& landBoundaryIndex : landBoundaryIndices)
    {
        // If the start node has a valid mask
        const auto& [startIndex, endIndex] = landBoundaryIndex;

        if (nodeMask[startIndex])
        {
            m_validLandBoundaries.emplace_back(landBoundaryIndex);
        }
    }

    // Split every land boundary in order to get rid of closed land boundaries
    // (some pieces of software cannot handle closed polygons, only polylines)
    const auto numSegmentIndicesBeforeSplitting = static_cast<UInt>(m_validLandBoundaries.size());
    for (UInt i = 0; i < numSegmentIndicesBeforeSplitting; ++i)
    {
        const auto [startSegmentIndex, endSegmentIndex] = m_validLandBoundaries[i];
        if (endSegmentIndex - startSegmentIndex > 1)
        {
            const auto split = startSegmentIndex + (endSegmentIndex - startSegmentIndex) / 2;
            m_validLandBoundaries[i].second = split;
            m_validLandBoundaries.emplace_back(split, endSegmentIndex);
        }
    }
}

void LandBoundaries::FindNearestMeshBoundary(ProjectToLandBoundaryOption projectToLandBoundaryOption)
{
    if (m_landBoundary.IsEmpty() ||
        projectToLandBoundaryOption == ProjectToLandBoundaryOption::DoNotProjectToLandBoundary ||
        projectToLandBoundaryOption == ProjectToLandBoundaryOption::ToOriginalNetBoundary)
    {
        return;
    }

    if (projectToLandBoundaryOption == ProjectToLandBoundaryOption::OuterMeshBoundaryToLandBoundary ||
        projectToLandBoundaryOption == ProjectToLandBoundaryOption::InnerAndOuterMeshBoundaryToLandBoundary)
    {
        m_findOnlyOuterMeshBoundary = true;
    }

    Administrate();

    m_nodeMask.resize(m_mesh.GetNumNodes(), constants::missing::uintValue);
    m_faceMask.resize(m_mesh.GetNumFaces(), false);
    m_edgeMask.resize(m_mesh.GetNumEdges(), constants::missing::uintValue);
    m_meshNodesLandBoundarySegments.resize(m_mesh.GetNumNodes(), constants::missing::uintValue);
    m_nodesMinDistances.resize(m_mesh.GetNumNodes(), constants::missing::doubleValue);

    // Loop over the segments of the land boundary and assign each node to the land boundary segment index
    for (UInt landBoundarySegment = 0; landBoundarySegment < m_validLandBoundaries.size(); ++landBoundarySegment)
    {
        UInt numRejectedNodesInPath = 0;
        UInt numNodesInPath = 0;
        const bool pathFound = MakePath(landBoundarySegment, numNodesInPath, numRejectedNodesInPath);

        if (pathFound && projectToLandBoundaryOption == ProjectToLandBoundaryOption::InnerAndOuterMeshBoundaryToLandBoundary)
        {
            m_findOnlyOuterMeshBoundary = false;
            MakePath(landBoundarySegment, numNodesInPath, numRejectedNodesInPath);
            m_findOnlyOuterMeshBoundary = true;
        }
    }

    // Connect the m_mesh nodes
    if (m_findOnlyOuterMeshBoundary)
    {
        std::vector<UInt> connectedNodes;
        for (UInt e = 0; e < m_mesh.GetNumEdges(); ++e)
        {
            if (!m_mesh.IsEdgeOnBoundary(e))
            {
                continue;
            }

            AssignLandBoundaryPolylineToMeshNodes(e, true, connectedNodes, 0);
        }
    }
}

bool LandBoundaries::InitialiseNodeLocations(const bool initialize,
                                             const UInt edgeIndex,
                                             const std::vector<UInt>& nodes,
                                             const UInt numNodes,
                                             std::vector<UInt>& nodesLoc,
                                             UInt& numNodesLoc) const
{
    if (initialize)
    {

        if (!m_mesh.IsEdgeOnBoundary(edgeIndex) || m_mesh.GetEdge(edgeIndex).first == constants::missing::uintValue || m_mesh.GetEdge(edgeIndex).second == constants::missing::uintValue)
        {
            throw AlgorithmError("Cannot not assign segment to mesh nodes.");
        }

        const auto firstMeshNode = m_mesh.GetEdge(edgeIndex).first;
        const auto secondMeshNode = m_mesh.GetEdge(edgeIndex).second;

        if (m_meshNodesLandBoundarySegments[firstMeshNode] != constants::missing::uintValue &&
            m_meshNodesLandBoundarySegments[secondMeshNode] == constants::missing::uintValue &&
            m_nodeMask[firstMeshNode] != constants::missing::uintValue &&
            m_nodeMask[secondMeshNode] != constants::missing::uintValue)
        {
            nodesLoc.resize(3);
            nodesLoc[0] = firstMeshNode;
            nodesLoc[1] = secondMeshNode;
            numNodesLoc = 2;
        }
        else if (m_meshNodesLandBoundarySegments[firstMeshNode] == constants::missing::uintValue &&
                 m_meshNodesLandBoundarySegments[secondMeshNode] != constants::missing::uintValue &&
                 m_nodeMask[firstMeshNode] != constants::missing::uintValue &&
                 m_nodeMask[secondMeshNode] != constants::missing::uintValue)
        {
            nodesLoc.resize(3);
            nodesLoc[0] = secondMeshNode;
            nodesLoc[1] = firstMeshNode;
            numNodesLoc = 2;
        }
        else
        {
            // not a valid edge
            return false;
        }
    }
    else
    {
        nodesLoc.resize(numNodes + 1);
        numNodesLoc = numNodes;
        std::copy(nodes.begin(), nodes.end(), nodesLoc.begin());
    }

    return true;
}

meshkernel::UInt LandBoundaries::GetSegmentIndex(const UInt nearestLandBoundaryNodeIndex) const
{
    UInt landboundarySegmentIndex = std::numeric_limits<UInt>::max();

    for (UInt s = 0; s < m_validLandBoundaries.size(); s++)
    {
        const auto& [startIndex, endIndex] = m_validLandBoundaries[s];

        if (nearestLandBoundaryNodeIndex >= startIndex && nearestLandBoundaryNodeIndex < endIndex)
        {
            landboundarySegmentIndex = s;
            break;
        }
    }

    return landboundarySegmentIndex;
}

void LandBoundaries::AssignLandBoundaryPolylineToMeshNodes(UInt edgeIndex, bool initialize, const std::vector<UInt>& nodes, UInt numNodes)
{
    if (m_landBoundary.IsEmpty())
    {
        return;
    }

    std::vector<UInt> nodesLoc;
    UInt numNodesLoc;

    if (!InitialiseNodeLocations(initialize, edgeIndex, nodes, numNodes, nodesLoc, numNodesLoc))
    {
        // Invalid edge
        return;
    }

    const auto maxNodes = *std::max_element(nodesLoc.begin(), nodesLoc.end() - 1);
    if (numNodesLoc > maxNodes)
    {
        return;
    }

    const auto lastVisitedNode = nodesLoc[numNodesLoc - 1];

    for (UInt e = 0; e < m_mesh.GetNumNodesEdges(lastVisitedNode); e++)
    {
        const auto edge = m_mesh.m_nodesEdges[lastVisitedNode][e];

        if (!m_mesh.IsEdgeOnBoundary(edge))
            continue;

        const auto otherNode = OtherNodeOfEdge(m_mesh.GetEdge(edge), lastVisitedNode);

        // path stopped
        if (m_nodeMask[otherNode] == constants::missing::uintValue)
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

        if (m_meshNodesLandBoundarySegments[otherNode] != constants::missing::uintValue)
        {
            // Now check landboundary for otherNode
            for (UInt n = 1; n < numNodesLoc; n++)
            {
                const auto meshNode = nodesLoc[n];
                const auto [minimumDistance,
                            pointOnLandBoundary,
                            nearestLandBoundaryNodeIndex,
                            edgeRatio] = NearestLandBoundarySegment(constants::missing::uintValue, m_mesh.Node(meshNode));

                // find the segment index of the found point
                UInt landboundarySegmentIndex = GetSegmentIndex(nearestLandBoundaryNodeIndex);

                if (landboundarySegmentIndex == std::numeric_limits<UInt>::max())
                {
                    throw AlgorithmError("No segment index found: cannot assign segment to mesh nodes.");
                }
                const auto& [startIndex, endIndex] = m_validLandBoundaries[landboundarySegmentIndex];

                if ((nearestLandBoundaryNodeIndex == startIndex && edgeRatio < 0.0) ||
                    (nearestLandBoundaryNodeIndex == endIndex - 1 && edgeRatio > 1.0))
                {
                    if (m_addLandboundaries)
                    {
                        AddLandBoundary(nodesLoc, numNodesLoc, lastVisitedNode);
                        m_meshNodesLandBoundarySegments[meshNode] = static_cast<UInt>(m_validLandBoundaries.size()) - 1; // last added ;and boundary
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

void LandBoundaries::AddLandBoundary(const std::vector<UInt>& nodesLoc, UInt numNodesLoc, UInt nodeIndex)
{
    if (m_landBoundary.IsEmpty())
    {
        return;
    }

    const auto startSegmentIndex = m_meshNodesLandBoundarySegments[nodesLoc[0]];
    const auto endSegmentIndex = m_meshNodesLandBoundarySegments[nodesLoc[numNodesLoc]];

    if (startSegmentIndex == constants::missing::uintValue || startSegmentIndex >= m_validLandBoundaries.size() ||
        endSegmentIndex == constants::missing::uintValue || endSegmentIndex >= m_validLandBoundaries.size())
    {
        throw std::invalid_argument("LandBoundaries::AddLandBoundary: Invalid segment index.");
    }

    // find start/end
    const auto& [startNodeLeftBoundary, endNodeLeftBoundary] = m_validLandBoundaries[startSegmentIndex];

    Point newNodeLeft;

    newNodeLeft = m_landBoundary.ClosestPoint(m_mesh.Node(nodeIndex), startNodeLeftBoundary, endNodeLeftBoundary, m_mesh.m_projection);

    Point newNodeRight;
    if (endSegmentIndex == startSegmentIndex)
    {
        newNodeRight = m_landBoundary.Node(startNodeLeftBoundary) + (m_landBoundary.Node(endNodeLeftBoundary) - newNodeLeft);
    }
    else
    {
        // find start/end
        const auto& [startNodeRightBoundary, endNodeRightBoundary] = m_validLandBoundaries[endSegmentIndex];
        newNodeRight = m_landBoundary.ClosestPoint(m_mesh.Node(nodeIndex), startNodeRightBoundary, endNodeRightBoundary, m_mesh.m_projection);
    }

    // Update land boundary nodes
    m_landBoundary.AddSegment(newNodeLeft, newNodeRight);

    // Update segment indices
    m_validLandBoundaries.emplace_back(std::make_pair<UInt, UInt>(static_cast<UInt>(m_landBoundary.GetNumNodes()) - 3, static_cast<UInt>(m_landBoundary.GetNumNodes()) - 2));
}

bool LandBoundaries::StopPathSearch(const UInt landBoundaryIndex, const UInt currentNode)
{
    bool stopPathSearch = true;

    if (m_meshNodesLandBoundarySegments[currentNode] != constants::missing::uintValue)
    {
        // Multiple boundary segments: take the nearest
        const auto previousLandBoundarySegment = m_meshNodesLandBoundarySegments[currentNode];

        const auto [previousMinDistance,
                    nodeOnPreviousLandBoundary,
                    nodeOnPreviousLandBoundaryNodeIndex,
                    previousLandBoundaryEdgeRatio] = NearestLandBoundarySegment(previousLandBoundarySegment, m_mesh.Node(currentNode));

        const auto [distanceFromLandBoundary,
                    nodeOnLandBoundary,
                    currentNodeLandBoundaryNodeIndex,
                    currentNodeEdgeRatio] = NearestLandBoundarySegment(landBoundaryIndex, m_mesh.Node(currentNode));

        const auto minDinstanceFromLandBoundaryCurrentNode = m_nodesMinDistances[currentNode];

        if (distanceFromLandBoundary <= previousMinDistance &&
            distanceFromLandBoundary < m_minDistanceFromLandFactor * minDinstanceFromLandBoundaryCurrentNode)
        {
            stopPathSearch = false;
        }
    }
    else
    {
        if (IsEqual(m_nodesMinDistances[currentNode], constants::missing::doubleValue))
        {
            const auto [minDinstanceFromLandBoundary,
                        nodeOnLandBoundary,
                        currentNodeLandBoundaryNodeIndex,
                        currentNodeEdgeRatio] = NearestLandBoundarySegment(constants::missing::uintValue, m_mesh.Node(currentNode));

            m_nodesMinDistances[currentNode] = minDinstanceFromLandBoundary;
        }

        const auto [distanceFromLandBoundary,
                    nodeOnLandBoundary,
                    currentNodeLandBoundaryNodeIndex,
                    currentNodeEdgeRatio] = NearestLandBoundarySegment(landBoundaryIndex, m_mesh.Node(currentNode));

        if (distanceFromLandBoundary < m_minDistanceFromLandFactor * m_nodesMinDistances[currentNode] &&
            (!m_findOnlyOuterMeshBoundary || m_mesh.GetNodeType(currentNode) == MeshNodeType::Boundary || m_mesh.GetNodeType(currentNode) == MeshNodeType::Corner))
        {
            stopPathSearch = false;
        }
    }

    return stopPathSearch;
}

bool LandBoundaries::MakePath(UInt landBoundaryIndex, UInt& numNodesInPath, UInt& numRejectedNodesInPath)
{
    numNodesInPath = 0;
    numRejectedNodesInPath = 0;
    if (m_landBoundary.IsEmpty())
    {
        return false;
    }

    const auto& [startIndex, endIndex] = m_validLandBoundaries[landBoundaryIndex];
    if (startIndex >= m_landBoundary.GetNumNodes() || startIndex >= endIndex)
    {
        throw std::invalid_argument("LandBoundaries::MakePath: Invalid boundary index.");
    }

    // Fractional location of the projected outer nodes(min and max) on the land boundary segment
    ComputeMeshNodeMask(landBoundaryIndex);
    const auto [startMeshNode, endMeshNode] = FindStartEndMeshNodesDijkstraAlgorithm(landBoundaryIndex);
    if (startMeshNode == constants::missing::uintValue)
    {
        return false;
    }

    const auto connectedNodeEdges = ShortestPath(landBoundaryIndex, startMeshNode);
    auto lastSegment = m_meshNodesLandBoundarySegments[endMeshNode];
    UInt lastNode = constants::missing::uintValue;
    auto currentNode = endMeshNode;
    UInt numConnectedNodes = 0;

    while (true)
    {
        if (StopPathSearch(landBoundaryIndex, currentNode))
        {
            if (numConnectedNodes == 1 && lastSegment != constants::missing::uintValue)
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

            m_meshNodesLandBoundarySegments[lastNode] = landBoundaryIndex;
        }

        if (currentNode == startMeshNode)
        {
            break;
        }

        const auto nextEdgeIndex = connectedNodeEdges[currentNode];
        if (nextEdgeIndex == constants::missing::uintValue || nextEdgeIndex >= m_mesh.GetNumEdges())
        {
            break;
        }

        currentNode = OtherNodeOfEdge(m_mesh.GetEdge(nextEdgeIndex), currentNode);

        if (currentNode == constants::missing::uintValue || currentNode >= m_mesh.GetNumNodes())
        {
            break;
        }
    }

    if (numConnectedNodes == 1)
    {
        m_meshNodesLandBoundarySegments[lastNode] = lastSegment;
    }

    return true;
}

void LandBoundaries::ComputeMeshNodeMask(UInt landBoundaryIndex)
{
    if (m_landBoundary.IsEmpty())
    {
        return;
    }

    const auto& [startLandBoundaryIndex, endLandBoundaryIndex] = m_validLandBoundaries[landBoundaryIndex];

    // Try to find a face crossed by the current land boundary polyline:
    // 1. One of the land boundary nodes is inside a face
    // 2. Or one of the land boundary segments is crossing a boundary mesh edge
    UInt crossedFaceIndex = constants::missing::uintValue;
    for (auto i = startLandBoundaryIndex; i < endLandBoundaryIndex; i++)
    {
        crossedFaceIndex = m_nodeFaceIndices[i];
        if (crossedFaceIndex != constants::missing::uintValue)
        {
            break;
        }

        auto [face, edge] = m_mesh.IsSegmentCrossingABoundaryEdge(m_landBoundary.Node(i), m_landBoundary.Node(i + 1));
        crossedFaceIndex = face;
        if (crossedFaceIndex != constants::missing::uintValue)
        {
            break;
        }
    }

    std::fill(m_nodeMask.begin(), m_nodeMask.end(), constants::missing::uintValue);
    if (m_landMask)
    {
        std::fill(m_faceMask.begin(), m_faceMask.end(), false);
        std::fill(m_edgeMask.begin(), m_edgeMask.end(), constants::missing::uintValue);
        // m_faceMask assumes crossedFace has already been done.
        if (crossedFaceIndex != constants::missing::uintValue)
        {
            m_faceMask[crossedFaceIndex] = true;
        }

        std::vector<UInt> landBoundaryFaces{crossedFaceIndex};
        MaskMeshFaceMask(landBoundaryIndex, landBoundaryFaces);

        // Mask all nodes of the masked faces
        for (UInt f = 0; f < m_mesh.GetNumFaces(); f++)
        {
            if (m_faceMask[f])
            {
                for (UInt n = 0; n < m_mesh.GetNumFaceEdges(f); n++)
                {
                    m_nodeMask[m_mesh.m_facesNodes[f][n]] = landBoundaryIndex;
                }
            }
        }
    }
    else
    {
        for (auto& e : m_nodeMask)
        {
            e = landBoundaryIndex;
        }
    }

    for (UInt n = 0; n < m_mesh.GetNumNodes(); n++)
    {
        if (m_nodeMask[n] != constants::missing::uintValue)
        {
            const bool inPolygon = m_polygons.IsPointInPolygon(m_mesh.Node(n), 0);
            if (!inPolygon)
            {
                m_nodeMask[n] = constants::missing::uintValue;
            }
        }
    }
}

bool LandBoundaries::ContainsCrossedFace(const UInt landBoundaryIndex, const UInt otherFace)
{
    bool isFaceFound = false;
    for (const auto& edge : m_mesh.m_facesEdges[otherFace])
    {
        if (m_edgeMask[edge] == 1)
        {
            // Previously visited crossed edge
            isFaceFound = true;
            continue;
        }
        if (m_edgeMask[edge] == 0)
        {
            // Previously visited uncrossed edge
            continue;
        }

        // Visited edge
        m_edgeMask[edge] = 0;
        const auto landBoundaryNode = IsMeshEdgeCloseToLandBoundaries(landBoundaryIndex, edge);

        if (landBoundaryNode != constants::missing::uintValue)
        {
            m_edgeMask[edge] = 1;
            isFaceFound = true;
        }
    }

    return isFaceFound;
}

void LandBoundaries::MaskFacesCloseToBoundary(const UInt landBoundaryIndex)
{
    for (UInt e = 0; e < m_mesh.GetNumEdges(); e++)
    {
        // only boundary edges are considered
        if (!m_mesh.IsEdgeOnBoundary(e))
        {
            continue;
        }

        const auto face = m_mesh.m_edgesFaces[e][0];
        // already masked
        if (m_faceMask[face])
        {
            continue;
        }

        for (const auto& edge : m_mesh.m_facesEdges[face])
        {
            const auto landBoundaryNode = IsMeshEdgeCloseToLandBoundaries(landBoundaryIndex, edge);
            if (landBoundaryNode != constants::missing::uintValue)
            {
                m_faceMask[face] = true;
                break;
            }
        }
    }
}

void LandBoundaries::MaskMeshFaceMask(UInt landBoundaryIndex, const std::vector<UInt>& initialFaces)
{
    if (m_landBoundary.IsEmpty())
    {
        return;
    }

    std::vector<UInt> nextFaces;
    nextFaces.reserve(initialFaces.size());
    for (const auto& face : initialFaces)
    {
        // No face was crossed by the land boundary: mask boundary faces only
        // These are the faces that are close (up to a certain tolerance) by a land boundary
        if (face == constants::missing::uintValue)
        {
            MaskFacesCloseToBoundary(landBoundaryIndex);
        }
        else
        {
            // Face is crossed
            if (m_mesh.GetNumFaces() < constants::geometric::numNodesInTriangle)
            {
                continue;
            }

            for (const auto& currentEdge : m_mesh.m_facesEdges[face])
            {
                // If it is a boundary edge, continue
                if (m_mesh.IsEdgeOnBoundary(currentEdge))
                {
                    continue;
                }

                const auto otherFace = face == m_mesh.m_edgesFaces[currentEdge][0] ? m_mesh.m_edgesFaces[currentEdge][1] : m_mesh.m_edgesFaces[currentEdge][0];

                // Already masked
                if (m_faceMask[otherFace])
                {
                    continue;
                }

                m_faceMask[otherFace] = ContainsCrossedFace(landBoundaryIndex, otherFace);

                if (m_faceMask[otherFace])
                {
                    nextFaces.emplace_back(otherFace);
                }
            }
        }
    }

    if (!nextFaces.empty())
    {
        MaskMeshFaceMask(landBoundaryIndex, nextFaces);
    }
}

void LandBoundaries::GetLandBoundaryNode(const double closeDistance,
                                         const Point& firstMeshNode,
                                         const Point& secondMeshNode,
                                         const UInt currentNode,
                                         UInt& landBoundaryNode,
                                         bool& isWithinSegment) const
{
    const double landBoundaryLength = ComputeSquaredDistance(m_landBoundary.Node(currentNode),
                                                             m_landBoundary.Node(currentNode + 1),
                                                             m_mesh.m_projection);

    isWithinSegment = false;
    if (landBoundaryLength > 0.0)
    {

        const auto [distanceFromLandBoundaryFirstMeshNode, normalPoint, ratioFirstMeshNode] = DistanceFromLine(firstMeshNode,
                                                                                                               m_landBoundary.Node(currentNode),
                                                                                                               m_landBoundary.Node(currentNode + 1),
                                                                                                               m_mesh.m_projection);

        if (distanceFromLandBoundaryFirstMeshNode < closeDistance)
        {
            landBoundaryNode = currentNode;
            // The projection of firstMeshNode is within the segment currentNode / currentNode + 1
            if (ratioFirstMeshNode >= 0.0 && ratioFirstMeshNode <= 1.0)
            {
                isWithinSegment = true;
            }
        }
        else
        {
            // Check the second point
            const auto [distanceFromLandBoundarySecondMeshNode, normalPoint, ratioSecondMeshNode] = DistanceFromLine(secondMeshNode,
                                                                                                                     m_landBoundary.Node(currentNode),
                                                                                                                     m_landBoundary.Node(currentNode + 1),
                                                                                                                     m_mesh.m_projection);

            if (distanceFromLandBoundarySecondMeshNode < closeDistance)
            {
                landBoundaryNode = currentNode;
                // The projection of secondMeshNode is within the segment currentNode / currentNode + 1
                if (ratioSecondMeshNode >= 0.0 && ratioSecondMeshNode <= 1.0)
                {
                    isWithinSegment = true;
                }
            }
        }
    }
}

meshkernel::UInt LandBoundaries::IsMeshEdgeCloseToLandBoundaries(UInt landBoundaryIndex, UInt edge)
{
    UInt landBoundaryNode = constants::missing::uintValue;
    if (m_landBoundary.IsEmpty())
    {
        return landBoundaryNode;
    }

    const auto& [startLandBoundaryIndex, endLandBoundaryIndex] = m_validLandBoundaries[landBoundaryIndex];

    const auto startNode = std::max(std::min(static_cast<UInt>(0), endLandBoundaryIndex - 1), startLandBoundaryIndex);
    if (m_mesh.GetEdge(edge).first == constants::missing::uintValue || m_mesh.GetEdge(edge).second == constants::missing::uintValue)
    {
        return landBoundaryNode;
    }

    const auto firstMeshNode = m_mesh.Node(m_mesh.GetEdge(edge).first);
    const auto secondMeshNode = m_mesh.Node(m_mesh.GetEdge(edge).second);

    const double meshEdgeLength = ComputeDistance(firstMeshNode, secondMeshNode, m_mesh.m_projection);
    const double distanceFactor = m_findOnlyOuterMeshBoundary ? m_closeToLandBoundaryFactor : m_closeWholeMeshFactor;

    const double closeDistance = meshEdgeLength * distanceFactor;

    // Now search over the land boundaries
    auto currentNode = startNode;
    UInt searchIterations = 0;
    int stepNode = 0;
    const UInt maximumNumberOfIterations = 3;
    while (searchIterations < maximumNumberOfIterations)
    {
        bool isWithinSegment = false;

        GetLandBoundaryNode(closeDistance, firstMeshNode, secondMeshNode, currentNode, landBoundaryNode, isWithinSegment);
        if (isWithinSegment)
        {
            break;
        }

        // Search the next land boundary edge if projection is not within is within the segment currentNode / currentNode + 1
        searchIterations = 0;
        while ((searchIterations == 0 || currentNode < startLandBoundaryIndex || currentNode > endLandBoundaryIndex - 1) && searchIterations < 3)
        {
            searchIterations += 1;
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

std::tuple<meshkernel::UInt, meshkernel::UInt> LandBoundaries::FindStartEndMeshNodesDijkstraAlgorithm(UInt landBoundaryIndex)
{
    if (m_landBoundary.IsEmpty())
    {
        return {constants::missing::uintValue, constants::missing::uintValue};
    }

    const auto& [startLandBoundaryIndex, endLandBoundaryIndex] = m_validLandBoundaries[landBoundaryIndex];
    const auto leftIndex = endLandBoundaryIndex - 1;
    const auto rightIndex = startLandBoundaryIndex;

    // Compute the start and end point of the land boundary respectively
    const auto nextLeftIndex = std::min(leftIndex + 1, endLandBoundaryIndex);
    const Point startPoint = m_landBoundary.Node(nextLeftIndex);
    const Point endPoint = m_landBoundary.Node(rightIndex);

    // Get the edges that are closest the land boundary
    auto minDistStart = std::numeric_limits<double>::max();
    auto minDistEnd = std::numeric_limits<double>::max();
    UInt startEdge = constants::missing::uintValue;
    UInt endEdge = constants::missing::uintValue;

    for (UInt e = 0; e < m_mesh.GetNumEdges(); e++)
    {
        // If the edge has an invalid node, continue
        if (m_mesh.GetEdge(e).first == constants::missing::uintValue || m_mesh.GetEdge(e).second == constants::missing::uintValue)
        {
            continue;
        }

        // Use only edges with both nodes masked
        if (m_nodeMask[m_mesh.GetEdge(e).first] == constants::missing::uintValue || m_nodeMask[m_mesh.GetEdge(e).second] == constants::missing::uintValue)
        {
            continue;
        }

        const auto [distanceFromFirstMeshNode, normalFirstMeshNode, ratioFirstMeshNode] = DistanceFromLine(startPoint,
                                                                                                           m_mesh.Node(m_mesh.GetEdge(e).first),
                                                                                                           m_mesh.Node(m_mesh.GetEdge(e).second),
                                                                                                           m_mesh.m_projection);

        const auto [distanceFromSecondMeshNode, normalSecondMeshNode, ratioSecondMeshNode] = DistanceFromLine(endPoint,
                                                                                                              m_mesh.Node(m_mesh.GetEdge(e).first),
                                                                                                              m_mesh.Node(m_mesh.GetEdge(e).second),
                                                                                                              m_mesh.m_projection);

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

    if (startEdge == constants::missing::uintValue || endEdge == constants::missing::uintValue)
    {
        return {constants::missing::uintValue, constants::missing::uintValue};
    }

    const auto startMeshNode = FindStartEndMeshNodesFromEdges(startEdge, startPoint);
    const auto endMeshNode = FindStartEndMeshNodesFromEdges(endEdge, endPoint);

    return {startMeshNode, endMeshNode};
}

meshkernel::UInt LandBoundaries::FindStartEndMeshNodesFromEdges(UInt edge, Point point) const
{
    if (m_landBoundary.IsEmpty())
    {
        return constants::missing::uintValue;
    }

    const auto firstMeshNodeIndex = m_mesh.GetEdge(edge).first;
    const auto secondMeshNodeIndex = m_mesh.GetEdge(edge).second;
    const auto firstDistance = ComputeSquaredDistance(m_mesh.Node(firstMeshNodeIndex), point, m_mesh.m_projection);
    const auto secondDistance = ComputeSquaredDistance(m_mesh.Node(secondMeshNodeIndex), point, m_mesh.m_projection);

    if (firstDistance <= secondDistance)
    {
        return firstMeshNodeIndex;
    }
    return secondMeshNodeIndex;
}

std::vector<meshkernel::UInt> LandBoundaries::ShortestPath(UInt landBoundaryIndex,
                                                           UInt startMeshNode)
{
    std::vector<UInt> connectedNodeEdges;
    if (m_landBoundary.IsEmpty())
    {
        return connectedNodeEdges;
    }

    connectedNodeEdges.resize(m_mesh.GetNumNodes(), constants::missing::uintValue);
    std::fill(connectedNodeEdges.begin(), connectedNodeEdges.end(), constants::missing::uintValue);

    // Infinite distance for all nodes
    std::vector<double> nodeDistances(m_mesh.GetNumNodes(), std::numeric_limits<double>::max());
    std::vector<bool> isVisited(m_mesh.GetNumNodes(), false);

    auto currentNodeIndex = startMeshNode;
    nodeDistances[startMeshNode] = 0.0;
    while (true)
    {
        isVisited[currentNodeIndex] = true;
        const Point currentNode = m_mesh.Node(currentNodeIndex);

        const auto [currentNodeDistance,
                    currentNodeOnLandBoundary,
                    currentNodeLandBoundaryNodeIndex,
                    currentNodeEdgeRatio] = NearestLandBoundarySegment(landBoundaryIndex, currentNode);

        if (currentNodeLandBoundaryNodeIndex == constants::missing::uintValue)
        {
            throw AlgorithmError("ShortestPath: Cannot compute the nearest node on the land boundary.");
        }

        for (const auto& edgeIndex : m_mesh.m_nodesEdges[currentNodeIndex])
        {
            if (m_mesh.GetEdge(edgeIndex).first == constants::missing::uintValue || m_mesh.GetEdge(edgeIndex).second == constants::missing::uintValue)
            {
                continue;
            }

            const auto neighbouringNodeIndex = OtherNodeOfEdge(m_mesh.GetEdge(edgeIndex), currentNodeIndex);

            if (isVisited[neighbouringNodeIndex])
            {
                continue;
            }

            const auto neighboringNode = m_mesh.Node(neighbouringNodeIndex);

            const auto [neighboringNodeDistance,
                        neighboringNodeOnLandBoundary,
                        neighboringNodeLandBoundaryNodeIndex,
                        neighboringNodeEdgeRatio] = NearestLandBoundarySegment(landBoundaryIndex, neighboringNode);

            double maximumDistance = std::max(currentNodeDistance, neighboringNodeDistance);

            if (currentNodeLandBoundaryNodeIndex < neighboringNodeLandBoundaryNodeIndex)
            {
                for (auto n = currentNodeLandBoundaryNodeIndex + 1; n < neighboringNodeLandBoundaryNodeIndex; ++n)
                {
                    const auto [middlePointDistance, middlePointOnLandBoundary, ratio] = DistanceFromLine(m_landBoundary.Node(n), currentNode, neighboringNode, m_mesh.m_projection);
                    if (middlePointDistance > maximumDistance)
                    {
                        maximumDistance = middlePointDistance;
                    }
                }
            }
            else if (currentNodeLandBoundaryNodeIndex > neighboringNodeLandBoundaryNodeIndex)
            {
                for (auto n = neighboringNodeLandBoundaryNodeIndex + 1; n < currentNodeLandBoundaryNodeIndex; ++n)
                {
                    const auto [middlePointDistance, middlePointOnLandBoundary, ratio] = DistanceFromLine(m_landBoundary.Node(n), currentNode, neighboringNode, m_mesh.m_projection);
                    if (middlePointDistance > maximumDistance)
                    {
                        maximumDistance = middlePointDistance;
                    }
                }
            }

            // In case of netboundaries only: set penalty when edge is not on the boundary
            if (m_findOnlyOuterMeshBoundary && !m_mesh.IsEdgeOnBoundary(edgeIndex))
            {
                maximumDistance = 1e6 * maximumDistance;
            }

            const double edgeLength = ComputeDistance(currentNode, neighboringNode, m_mesh.m_projection);
            const double correctedDistance = nodeDistances[currentNodeIndex] + edgeLength * maximumDistance;

            if (correctedDistance < nodeDistances[neighbouringNodeIndex])
            {
                nodeDistances[neighbouringNodeIndex] = correctedDistance;
                connectedNodeEdges[neighbouringNodeIndex] = edgeIndex;
            }
        }

        // Linear search with masking
        currentNodeIndex = 0;
        double minValue = std::numeric_limits<double>::max();
        for (UInt n = 0; n < m_mesh.GetNumNodes(); n++)
        {
            if (m_nodeMask[n] == landBoundaryIndex && !isVisited[n] && nodeDistances[n] < minValue)
            {
                currentNodeIndex = n;
                minValue = nodeDistances[n];
            }
        }

        if (currentNodeIndex >= m_mesh.GetNumNodes() ||
            IsEqual(nodeDistances[currentNodeIndex], std::numeric_limits<double>::max()) ||
            isVisited[currentNodeIndex])
        {
            break;
        }
    }

    return connectedNodeEdges;
}

std::tuple<double, meshkernel::Point, meshkernel::UInt, double> LandBoundaries::NearestLandBoundarySegment(UInt segmentIndex, const Point& node) const
{
    double minimumDistance = std::numeric_limits<double>::max();
    Point pointOnLandBoundary = node;
    UInt nearestLandBoundaryNodeIndex = constants::missing::uintValue;
    double edgeRatio = -1.0;

    if (m_landBoundary.IsEmpty())
    {
        return {minimumDistance, pointOnLandBoundary, nearestLandBoundaryNodeIndex, edgeRatio};
    }

    const auto startLandBoundaryIndex = segmentIndex == constants::missing::uintValue ? 0 : m_validLandBoundaries[segmentIndex].first;
    const auto endLandBoundaryIndex = segmentIndex == constants::missing::uintValue ? static_cast<UInt>(m_landBoundary.GetNumNodes() - 1) : m_validLandBoundaries[segmentIndex].second;

    for (auto n = startLandBoundaryIndex; n < endLandBoundaryIndex; n++)
    {
        if (!m_landBoundary.Node(n).IsValid() || !m_landBoundary.Node(n + 1).IsValid())
        {
            continue;
        }

        const auto [distanceFromLandBoundary, normalPoint, ratio] = DistanceFromLine(node, m_landBoundary.Node(n), m_landBoundary.Node(n + 1), m_mesh.m_projection);

        if (distanceFromLandBoundary > 0.0 && distanceFromLandBoundary < minimumDistance)
        {
            minimumDistance = distanceFromLandBoundary;
            pointOnLandBoundary = normalPoint;
            nearestLandBoundaryNodeIndex = n;
            edgeRatio = ratio;
        }
    }

    return {minimumDistance, pointOnLandBoundary, nearestLandBoundaryNodeIndex, edgeRatio};
}

std::unique_ptr<meshkernel::UndoAction> LandBoundaries::SnapMeshToLandBoundaries() const
{
    if (m_landBoundary.IsEmpty() || m_meshNodesLandBoundarySegments.empty())
    {
        return nullptr;
    }

    std::unique_ptr<CompoundUndoAction> action = CompoundUndoAction::Create();

    const auto numNodes = m_mesh.GetNumNodes();
    for (UInt n = 0; n < numNodes; ++n)
    {
        if (m_mesh.GetNodeType(n) == MeshNodeType::Internal || m_mesh.GetNodeType(n) == MeshNodeType::Boundary || m_mesh.GetNodeType(n) == MeshNodeType::Corner)
        {
            const auto meshNodeToLandBoundarySegment = m_meshNodesLandBoundarySegments[n];
            if (meshNodeToLandBoundarySegment == constants::missing::uintValue)
            {
                continue;
            }

            const auto [minimumDistance,
                        pointOnLandBoundary,
                        nearestLandBoundaryNodeIndex,
                        edgeRatio] = NearestLandBoundarySegment(meshNodeToLandBoundarySegment, m_mesh.Node(n));

            action->Add(m_mesh.ResetNode(n, pointOnLandBoundary));
        }
    }

    return action;
}
