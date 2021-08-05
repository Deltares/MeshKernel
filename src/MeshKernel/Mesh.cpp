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

#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/Mesh.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Polygons.hpp>

using meshkernel::Mesh;

Mesh::Mesh(const std::vector<Edge>& edges,
           const std::vector<Point>& nodes,
           Projection projection) : m_nodes(nodes), m_edges(edges), m_projection(projection){};

void Mesh::NodeAdministration()
{
    // assume no duplicated links
    for (auto e = 0; e < GetNumEdges(); e++)
    {
        const auto firstNode = m_edges[e].first;
        const auto secondNode = m_edges[e].second;

        if (firstNode == sizetMissingValue || secondNode == sizetMissingValue)
        {
            continue;
        }

        if (m_nodesNumEdges[firstNode] >= maximumNumberOfEdgesPerNode || m_nodesNumEdges[secondNode] >= maximumNumberOfEdgesPerNode)
        {
            continue;
        }

        // Search for previously connected edges
        auto alreadyAddedEdge = false;
        for (auto i = 0; i < m_nodesNumEdges[firstNode]; ++i)
        {
            if (const auto currentEdge = m_edges[m_nodesEdges[firstNode][i]]; currentEdge.first == secondNode || currentEdge.second == secondNode)
            {
                alreadyAddedEdge = true;
                break;
            }
        }
        if (!alreadyAddedEdge)
        {
            m_nodesEdges[firstNode][m_nodesNumEdges[firstNode]] = e;
            m_nodesNumEdges[firstNode]++;
        }

        // Search for previously connected edges
        alreadyAddedEdge = false;
        for (auto i = 0; i < m_nodesNumEdges[secondNode]; ++i)
        {
            if (const auto currentEdge = m_edges[m_nodesEdges[secondNode][i]]; currentEdge.first == firstNode || currentEdge.second == firstNode)
            {
                alreadyAddedEdge = true;
                break;
            }
        }
        if (!alreadyAddedEdge)
        {
            m_nodesEdges[secondNode][m_nodesNumEdges[secondNode]] = e;
            m_nodesNumEdges[secondNode]++;
        }
    }

    // resize
    for (auto n = 0; n < GetNumNodes(); n++)
    {
        m_nodesEdges[n].resize(m_nodesNumEdges[n]);
    }
}

void Mesh::DeleteInvalidNodesAndEdges()
{

    // Mask nodes connected to valid edges
    std::vector<bool> connectedNodes(m_nodes.size(), false);
    size_t numInvalidEdges = 0;

    for (const auto& [firstNode, secondNode] : m_edges)
    {
        if (firstNode == sizetMissingValue || secondNode == sizetMissingValue)
        {
            numInvalidEdges++;
            continue;
        }
        connectedNodes[firstNode] = true;
        connectedNodes[secondNode] = true;
    }

    // Count all invalid nodes (note: there might be nodes that are not connected to an edge)
    size_t numInvalidNodes = 0;
    for (auto n = 0; n < m_nodes.size(); ++n)
    {
        // invalidate nodes that are not connected
        if (!connectedNodes[n])
        {
            m_nodes[n] = {doubleMissingValue, doubleMissingValue};
        }

        if (!m_nodes[n].IsValid())
        {
            numInvalidNodes++;
        }
    }

    // If nothing to invalidate return
    if (numInvalidEdges == 0 && numInvalidNodes == 0)
    {
        return;
    }

    // Flag invalid nodes
    std::vector<size_t> validNodesIndices(m_nodes.size());
    std::fill(validNodesIndices.begin(), validNodesIndices.end(), sizetMissingValue);
    size_t validIndex = 0;
    for (auto n = 0; n < m_nodes.size(); ++n)
    {
        if (m_nodes[n].IsValid())
        {
            validNodesIndices[n] = validIndex;
            validIndex++;
        }
    }

    // Flag invalid edges
    for (auto& [firstNode, secondNode] : m_edges)
    {
        if (firstNode != sizetMissingValue && secondNode != sizetMissingValue && validNodesIndices[firstNode] != sizetMissingValue && validNodesIndices[secondNode] != sizetMissingValue)
        {
            firstNode = validNodesIndices[firstNode];
            secondNode = validNodesIndices[secondNode];
            continue;
        }

        firstNode = sizetMissingValue;
        secondNode = sizetMissingValue;
    }

    // Remove invalid nodes, without reducing capacity
    const auto endNodeVector = std::remove_if(m_nodes.begin(), m_nodes.end(), [](const Point& n) { return !n.IsValid(); });
    m_nodes.erase(endNodeVector, m_nodes.end());

    // Remove invalid edges, without reducing capacity
    const auto endEdgeVector = std::remove_if(m_edges.begin(), m_edges.end(), [](const Edge& e) { return e.first == sizetMissingValue || e.second == sizetMissingValue; });
    m_edges.erase(endEdgeVector, m_edges.end());
}

void Mesh::MergeTwoNodes(size_t firstNodeIndex, size_t secondNodeIndex)
{
    if (firstNodeIndex >= GetNumNodes() || secondNodeIndex >= GetNumNodes())
    {
        throw std::invalid_argument("Mesh::MergeTwoNodes: Either the first or the second node-index is invalid.");
    }

    auto edgeIndex = FindEdge(firstNodeIndex, secondNodeIndex);
    if (edgeIndex != sizetMissingValue)
    {
        m_edges[edgeIndex].first = sizetMissingValue;
        m_edges[edgeIndex].second = sizetMissingValue;
    }

    // check if there is another edge starting at firstEdgeOtherNode and ending at secondNode
    for (auto n = 0; n < m_nodesNumEdges[firstNodeIndex]; n++)
    {
        const auto firstEdgeIndex = m_nodesEdges[firstNodeIndex][n];
        const auto firstEdge = m_edges[firstEdgeIndex];
        const auto firstEdgeOtherNode = OtherNodeOfEdge(firstEdge, firstNodeIndex);
        if (firstEdgeOtherNode != sizetMissingValue && firstEdgeOtherNode != secondNodeIndex)
        {
            for (auto nn = 0; nn < m_nodesNumEdges[firstEdgeOtherNode]; nn++)
            {
                const auto secondEdgeIndex = m_nodesEdges[firstEdgeOtherNode][nn];
                auto secondEdge = m_edges[secondEdgeIndex];
                const auto secondNodeSecondEdge = OtherNodeOfEdge(secondEdge, firstEdgeOtherNode);
                if (secondNodeSecondEdge == secondNodeIndex)
                {
                    m_edges[secondEdgeIndex].first = sizetMissingValue;
                    m_edges[secondEdgeIndex].second = sizetMissingValue;
                }
            }
        }
    }

    // add all valid edges starting at secondNode
    std::vector<size_t> secondNodeEdges(maximumNumberOfEdgesPerNode, sizetMissingValue);
    size_t numSecondNodeEdges = 0;
    for (auto n = 0; n < m_nodesNumEdges[secondNodeIndex]; n++)
    {
        edgeIndex = m_nodesEdges[secondNodeIndex][n];
        if (m_edges[edgeIndex].first != sizetMissingValue)
        {
            secondNodeEdges[numSecondNodeEdges] = edgeIndex;
            numSecondNodeEdges++;
        }
    }

    // add all valid edges starting at firstNode are assigned to the second node
    for (auto n = 0; n < m_nodesNumEdges[firstNodeIndex]; n++)
    {
        edgeIndex = m_nodesEdges[firstNodeIndex][n];
        if (m_edges[edgeIndex].first != sizetMissingValue)
        {
            secondNodeEdges[numSecondNodeEdges] = edgeIndex;
            if (m_edges[edgeIndex].first == firstNodeIndex)
            {
                m_edges[edgeIndex].first = secondNodeIndex;
            }
            if (m_edges[edgeIndex].second == firstNodeIndex)
            {
                m_edges[edgeIndex].second = secondNodeIndex;
            }
            numSecondNodeEdges++;
        }
    }

    // re-assign edges to second node
    m_nodesEdges[secondNodeIndex] = std::vector<size_t>(secondNodeEdges.begin(), secondNodeEdges.begin() + numSecondNodeEdges);
    m_nodesNumEdges[secondNodeIndex] = numSecondNodeEdges;

    // remove edges to first node
    m_nodesEdges[firstNodeIndex] = std::vector<size_t>(0);
    m_nodesNumEdges[firstNodeIndex] = 0;
    m_nodes[firstNodeIndex] = {doubleMissingValue, doubleMissingValue};

    m_nodesRTreeRequiresUpdate = true;
    m_edgesRTreeRequiresUpdate = true;
}

void Mesh::MergeNodesInPolygon(const Polygons& polygon, double mergingDistance)
{
    // first filter the nodes in polygon
    std::vector<Point> filteredNodes(GetNumNodes());
    std::vector<size_t> originalNodeIndices(GetNumNodes(), sizetMissingValue);
    size_t index = 0;
    for (auto i = 0; i < GetNumNodes(); ++i)
    {
        const bool inPolygon = polygon.IsPointInPolygon(m_nodes[i], 0);
        if (inPolygon)
        {
            filteredNodes[index] = m_nodes[i];
            originalNodeIndices[index] = i;
            index++;
        }
    }
    filteredNodes.resize(index);

    // Update the R-Tree of the mesh nodes
    RTree nodesRtree;
    nodesRtree.BuildTree(filteredNodes);

    // merge the closest nodes
    auto const mergingDistanceSquared = mergingDistance * mergingDistance;
    for (auto i = 0; i < filteredNodes.size(); ++i)
    {
        nodesRtree.SearchPoints(filteredNodes[i], mergingDistanceSquared);

        const auto resultSize = nodesRtree.GetQueryResultSize();
        if (resultSize > 1)
        {
            for (auto j = 0; j < nodesRtree.GetQueryResultSize(); j++)
            {
                const auto nodeIndexInFilteredNodes = nodesRtree.GetQueryResult(j);
                if (nodeIndexInFilteredNodes != i)
                {
                    MergeTwoNodes(originalNodeIndices[i], originalNodeIndices[nodeIndexInFilteredNodes]);
                    nodesRtree.DeleteNode(i);
                }
            }
        }
    }

    AdministrateNodesEdges();
}

size_t Mesh::ConnectNodes(size_t startNode, size_t endNode)
{
    const auto edgeIndex = FindEdge(startNode, endNode);

    // The nodes are already connected
    if (edgeIndex != sizetMissingValue)
        return sizetMissingValue;

    // increment the edges container
    const auto newEdgeIndex = GetNumEdges();
    m_edges.resize(newEdgeIndex + 1);
    m_edges[newEdgeIndex].first = startNode;
    m_edges[newEdgeIndex].second = endNode;

    m_edgesRTreeRequiresUpdate = true;

    return newEdgeIndex;
}

size_t Mesh::InsertNode(const Point& newPoint)
{
    const auto newSize = GetNumNodes() + 1;
    const auto newNodeIndex = GetNumNodes();

    m_nodes.resize(newSize);
    m_nodesNumEdges.resize(newSize);
    m_nodesEdges.resize(newSize);

    m_nodes[newNodeIndex] = newPoint;
    m_nodesNumEdges[newNodeIndex] = 0;

    m_nodesRTreeRequiresUpdate = true;

    return newNodeIndex;
}

void Mesh::DeleteNode(size_t node)
{
    if (node >= GetNumNodes())
    {
        throw std::invalid_argument("Mesh::DeleteNode: The index of the node to be deleted does not exist.");
    }

    for (auto e = 0; e < m_nodesNumEdges[node]; e++)
    {
        const auto edgeIndex = m_nodesEdges[node][e];
        DeleteEdge(edgeIndex);
    }
    m_nodes[node] = {doubleMissingValue, doubleMissingValue};

    m_nodesRTreeRequiresUpdate = true;
}

void Mesh::DeleteEdge(size_t edge)
{
    if (edge == sizetMissingValue)
    {
        throw std::invalid_argument("Mesh::DeleteEdge: The index of the edge to be deleted does not exist.");
    }

    m_edges[edge].first = sizetMissingValue;
    m_edges[edge].second = sizetMissingValue;

    m_edgesRTreeRequiresUpdate = true;
}

void Mesh::ComputeEdgesLengths()
{
    auto const numEdges = GetNumEdges();
    m_edgeLengths.resize(numEdges, doubleMissingValue);
    for (auto e = 0; e < numEdges; e++)
    {
        auto const first = m_edges[e].first;
        auto const second = m_edges[e].second;
        m_edgeLengths[e] = ComputeDistance(m_nodes[first], m_nodes[second], m_projection);
    }
}

void Mesh::ComputeEdgesCenters()
{
    m_edgesCenters = ComputeEdgeCenters(m_nodes, m_edges);
}

size_t Mesh::FindCommonNode(size_t firstEdgeIndex, size_t secondEdgeIndex) const
{
    const auto firstEdgeFirstNode = m_edges[firstEdgeIndex].first;
    const auto firstEdgeEdgeSecondNode = m_edges[firstEdgeIndex].second;

    const auto secondEdgeFirstNode = m_edges[secondEdgeIndex].first;
    const auto secondEdgeSecondNode = m_edges[secondEdgeIndex].second;

    if (firstEdgeFirstNode == sizetMissingValue || firstEdgeEdgeSecondNode == sizetMissingValue || secondEdgeFirstNode == sizetMissingValue || secondEdgeSecondNode == sizetMissingValue)
    {
        throw std::invalid_argument("Mesh::FindCommonNode: At least one of the given edges is invalid.");
    }

    if (firstEdgeFirstNode == secondEdgeFirstNode || firstEdgeFirstNode == secondEdgeSecondNode)
    {
        return firstEdgeFirstNode;
    }
    if (firstEdgeEdgeSecondNode == secondEdgeFirstNode || firstEdgeEdgeSecondNode == secondEdgeSecondNode)
    {
        return firstEdgeEdgeSecondNode;
    }
    return sizetMissingValue;
}

size_t Mesh::FindEdge(size_t firstNodeIndex, size_t secondNodeIndex) const
{
    if (firstNodeIndex == sizetMissingValue || secondNodeIndex == sizetMissingValue)
    {
        throw std::invalid_argument("Mesh::FindEdge: Invalid node index.");
    }

    size_t edgeIndex = sizetMissingValue;
    for (auto n = 0; n < m_nodesNumEdges[firstNodeIndex]; n++)
    {
        const auto localEdgeIndex = m_nodesEdges[firstNodeIndex][n];
        const auto firstEdgeOtherNode = OtherNodeOfEdge(m_edges[localEdgeIndex], firstNodeIndex);
        if (firstEdgeOtherNode == secondNodeIndex)
        {
            edgeIndex = localEdgeIndex;
            break;
        }
    }
    return edgeIndex;
}

size_t Mesh::FindNodeCloseToAPoint(Point const& point, double searchRadius)
{
    if (GetNumNodes() <= 0)
    {
        throw std::invalid_argument("Mesh::FindNodeCloseToAPoint: There are no valid nodes.");
    }

    SearchNearestLocation(point, searchRadius * searchRadius, MeshLocations::Nodes);

    if (GetNumLocations(MeshLocations::Nodes) > 0)
    {
        return GetLocationsIndices(0, MeshLocations::Nodes);
    }

    throw AlgorithmError("Mesh::FindNodeCloseToAPoint: Could not find the node index close to a point.");
}

size_t Mesh::FindNodeCloseToAPoint(Point point, const std::vector<bool>& oneDNodeMask)
{
    if (GetNumNodes() <= 0)
    {
        throw std::invalid_argument("Mesh::FindNodeCloseToAPoint: There are no valid nodes.");
    }

    // create rtree a first time
    if (m_nodesRTree.Empty())
    {
        m_nodesRTree.BuildTree(m_nodes);
        m_nodesRTreeRequiresUpdate = false;
    }

    m_nodesRTree.SearchNearestPoint(point);
    const auto resultSize = m_nodesRTree.GetQueryResultSize();

    // no results found
    if (resultSize <= 0)
    {
        throw AlgorithmError("Mesh::FindNodeCloseToAPoint: query result size <= 0.");
    }

    // resultSize > 0, no node mask applied
    if (oneDNodeMask.empty())
    {
        return m_nodesRTree.GetQueryResult(0);
    }

    // resultSize > 0, a mask is applied
    for (auto index = 0; index < resultSize; ++index)
    {
        const auto nodeIndex = m_nodesRTree.GetQueryResult(index);
        if (oneDNodeMask[nodeIndex])
        {
            return nodeIndex;
        }
    }

    throw AlgorithmError("Mesh::FindNodeCloseToAPoint: Could not find the node index close to a point.");
}

size_t Mesh::FindEdgeCloseToAPoint(Point point)
{
    if (GetNumEdges() == 0)
    {
        throw std::invalid_argument("Mesh::GetNodeIndex: There are no valid edges.");
    }

    SearchNearestLocation(point, MeshLocations::Edges);

    if (GetNumLocations(MeshLocations::Edges) >= 1)
    {
        return GetLocationsIndices(0, MeshLocations::Edges);
    }

    throw AlgorithmError("Mesh::FindEdgeCloseToAPoint: Could not find the closest edge to a point.");
}

void Mesh::MoveNode(Point newPoint, size_t nodeindex)
{
    const Point nodeToMove = m_nodes.at(nodeindex);

    const auto dx = GetDx(nodeToMove, newPoint, m_projection);
    const auto dy = GetDy(nodeToMove, newPoint, m_projection);

    const auto distanceNodeToMoveFromNewPoint = std::sqrt(dx * dx + dy * dy);
    for (auto n = 0; n < GetNumNodes(); ++n)
    {
        const auto nodeDx = GetDx(m_nodes[n], nodeToMove, m_projection);
        const auto nodeDy = GetDy(m_nodes[n], nodeToMove, m_projection);
        const double distanceCurrentNodeFromNewPoint = std::sqrt(nodeDx * nodeDx + nodeDy * nodeDy);

        const auto factor = 0.5 * (1.0 + std::cos(std::min(distanceCurrentNodeFromNewPoint / distanceNodeToMoveFromNewPoint, 1.0) * M_PI));

        m_nodes[n].x += dx * factor;
        m_nodes[n].y += dy * factor;
    }

    m_nodesRTreeRequiresUpdate = true;
    m_edgesRTreeRequiresUpdate = true;
}

bool Mesh::IsFaceOnBoundary(size_t face) const
{

    bool isFaceOnBoundary = false;

    for (auto e = 0; e < GetNumFaceEdges(face); ++e)
    {
        const auto edge = m_facesEdges[face][e];
        if (IsEdgeOnBoundary(edge))
        {
            isFaceOnBoundary = true;
            break;
        }
    }
    return isFaceOnBoundary;
}

void Mesh::SortEdgesInCounterClockWiseOrder(size_t startNode, size_t endNode)
{

    std::vector<double> edgeAngles(maximumNumberOfEdgesPerNode);
    std::vector<std::size_t> indices(maximumNumberOfEdgesPerNode);
    std::vector<size_t> edgeNodeCopy(maximumNumberOfEdgesPerNode);
    for (auto n = startNode; n <= endNode; n++)
    {
        if (!m_nodes[n].IsValid())
        {
            continue;
        }

        double phi0 = 0.0;
        double phi;
        std::fill(edgeAngles.begin(), edgeAngles.end(), 0.0);
        for (auto edgeIndex = 0; edgeIndex < m_nodesNumEdges[n]; edgeIndex++)
        {

            auto firstNode = m_edges[m_nodesEdges[n][edgeIndex]].first;
            auto secondNode = m_edges[m_nodesEdges[n][edgeIndex]].second;
            if (firstNode == sizetMissingValue || secondNode == sizetMissingValue)
            {
                continue;
            }

            if (secondNode == n)
            {
                secondNode = firstNode;
                firstNode = n;
            }

            const auto deltaX = GetDx(m_nodes[secondNode], m_nodes[firstNode], m_projection);
            const auto deltaY = GetDy(m_nodes[secondNode], m_nodes[firstNode], m_projection);
            if (abs(deltaX) < minimumDeltaCoordinate && abs(deltaY) < minimumDeltaCoordinate)
            {
                if (deltaY < 0.0)
                {
                    phi = -M_PI / 2.0;
                }
                else
                {
                    phi = M_PI / 2.0;
                }
            }
            else
            {
                phi = atan2(deltaY, deltaX);
            }

            if (edgeIndex == 0)
            {
                phi0 = phi;
            }

            edgeAngles[edgeIndex] = phi - phi0;
            if (edgeAngles[edgeIndex] < 0.0)
            {
                edgeAngles[edgeIndex] = edgeAngles[edgeIndex] + 2.0 * M_PI;
            }
        }

        // Performing sorting
        indices.resize(m_nodesNumEdges[n]);
        edgeNodeCopy.clear();
        std::copy(m_nodesEdges[n].begin(), m_nodesEdges[n].end(), std::back_inserter(edgeNodeCopy));
        iota(indices.begin(), indices.end(), 0);
        sort(indices.begin(), indices.end(), [&](std::size_t const& i1, std::size_t const& i2) { return edgeAngles[i1] < edgeAngles[i2]; });

        for (std::size_t edgeIndex = 0; edgeIndex < m_nodesNumEdges[n]; edgeIndex++)
        {
            m_nodesEdges[n][edgeIndex] = edgeNodeCopy[indices[edgeIndex]];
        }
    }
}

void Mesh::BuildTree(MeshLocations meshLocation)
{
    if (meshLocation == MeshLocations::Nodes && m_nodesRTree.Empty())
    {
        m_nodesRTree.BuildTree(m_nodes);
        m_nodesRTreeRequiresUpdate = false;
    }

    if (meshLocation == MeshLocations::Edges && m_edgesRTree.Empty())
    {
        ComputeEdgesCenters();
        m_edgesRTree.BuildTree(m_edgesCenters);
        m_edgesRTreeRequiresUpdate = false;
    }

    if (meshLocation == MeshLocations::Faces && m_facesRTree.Empty())
    {
        m_facesRTree.BuildTree(m_facesCircumcenters);
    }
}

void Mesh::SearchNearestLocation(Point point, MeshLocations meshLocation)
{
    BuildTree(meshLocation);
    if (meshLocation == MeshLocations::Nodes)
    {
        m_nodesRTree.SearchNearestPoint(point);
    }

    if (meshLocation == MeshLocations::Edges)
    {
        m_edgesRTree.SearchNearestPoint(point);
    }

    if (meshLocation == MeshLocations::Faces)
    {
        m_facesRTree.SearchNearestPoint(point);
    }
}

void Mesh::SearchNearestLocation(Point point, double squaredRadius, MeshLocations meshLocation)
{
    BuildTree(meshLocation);
    if (meshLocation == MeshLocations::Nodes)
    {
        m_nodesRTree.SearchNearestPoint(point, squaredRadius);
    }

    if (meshLocation == MeshLocations::Edges)
    {
        m_edgesRTree.SearchNearestPoint(point, squaredRadius);
    }

    if (meshLocation == MeshLocations::Faces)
    {
        m_facesRTree.SearchNearestPoint(point, squaredRadius);
    }
}

void Mesh::SearchLocations(Point point, double squaredRadius, MeshLocations meshLocation)
{
    BuildTree(meshLocation);
    if (meshLocation == MeshLocations::Nodes)
    {
        m_nodesRTree.SearchPoints(point, squaredRadius);
    }

    if (meshLocation == MeshLocations::Edges)
    {
        m_edgesRTree.SearchPoints(point, squaredRadius);
    }

    if (meshLocation == MeshLocations::Faces)
    {
        m_facesRTree.SearchPoints(point, squaredRadius);
    }
}

size_t Mesh::GetNumLocations(MeshLocations meshLocation) const
{
    if (meshLocation == MeshLocations::Nodes)
    {
        return m_nodesRTree.GetQueryResultSize();
    }

    if (meshLocation == MeshLocations::Edges)
    {
        return m_edgesRTree.GetQueryResultSize();
    }

    if (meshLocation == MeshLocations::Faces)
    {
        return m_facesRTree.GetQueryResultSize();
    }

    return sizetMissingValue;
}

size_t Mesh::GetLocationsIndices(size_t index, MeshLocations meshLocation)
{
    if (meshLocation == MeshLocations::Nodes)
    {
        return m_nodesRTree.GetQueryResult(index);
    }

    if (meshLocation == MeshLocations::Edges)
    {
        return m_edgesRTree.GetQueryResult(index);
    }

    if (meshLocation == MeshLocations::Faces)
    {
        return m_facesRTree.GetQueryResult(index);
    }

    return sizetMissingValue;
}

void Mesh::AdministrateNodesEdges()
{
    DeleteInvalidNodesAndEdges();

    if (m_nodesRTreeRequiresUpdate && !m_nodesRTree.Empty())
    {
        m_nodesRTree.BuildTree(m_nodes);
        m_nodesRTreeRequiresUpdate = false;
    }

    if (m_edgesRTreeRequiresUpdate && !m_edgesRTree.Empty())
    {
        ComputeEdgesCenters();
        m_edgesRTree.BuildTree(m_edgesCenters);
        m_edgesRTreeRequiresUpdate = false;
    }

    // return if there are no nodes or no edges
    if (m_nodes.empty() || m_edges.empty())
    {
        return;
    }

    m_nodesEdges.resize(m_nodes.size());
    std::fill(m_nodesEdges.begin(), m_nodesEdges.end(), std::vector<size_t>(maximumNumberOfEdgesPerNode, sizetMissingValue));

    m_nodesNumEdges.resize(m_nodes.size());
    std::fill(m_nodesNumEdges.begin(), m_nodesNumEdges.end(), 0);

    NodeAdministration();

    SortEdgesInCounterClockWiseOrder(0, GetNumNodes() - 1);
}

double Mesh::ComputeMaxLengthSurroundingEdges(size_t node)
{

    if (m_edgeLengths.empty())
    {
        ComputeEdgesLengths();
    }

    auto maxEdgeLength = std::numeric_limits<double>::lowest();
    for (auto ee = 0; ee < m_nodesNumEdges[node]; ++ee)
    {
        const auto edge = m_nodesEdges[node][ee];
        maxEdgeLength = std::max(maxEdgeLength, m_edgeLengths[edge]);
    }

    return maxEdgeLength;
}

std::vector<meshkernel::Point> Mesh::ComputeLocations(MeshLocations location) const
{
    std::vector<Point> result;
    if (location == MeshLocations::Nodes)
    {
        result.reserve(GetNumNodes());
        for (const auto& n : m_nodes)
        {
            result.emplace_back(n);
        }
    }
    if (location == MeshLocations::Edges)
    {
        result.reserve(GetNumEdges());
        for (const auto& [firstNode, secondNode] : m_edges)
        {

            if (firstNode != sizetMissingValue && secondNode != sizetMissingValue)
            {
                result.emplace_back((m_nodes[firstNode] + m_nodes[secondNode]) * 0.5);
            }
            else
            {
                result.emplace_back(Point{});
            }
        }
    }
    if (location == MeshLocations::Faces)
    {
        result.reserve(GetNumFaces());
        for (const auto& massCentre : m_facesMassCenters)
        {
            result.emplace_back(massCentre);
        }
    }
    return result;
}

Mesh& Mesh::operator+=(Mesh const& rhs)
{
    if (m_projection != rhs.m_projection || rhs.GetNumNodes() == 0 || rhs.GetNumEdges() == 0)
    {
        throw std::invalid_argument("Mesh2D::operator+=: The two meshes cannot be added.");
    }

    const auto rhsNumNodes = rhs.GetNumNodes();
    const auto rhsNumEdges = rhs.GetNumEdges();

    auto const numNodes = GetNumNodes();
    auto const numEdges = GetNumEdges();
    m_edges.resize(GetNumEdges() + rhsNumEdges);
    m_nodes.resize(GetNumNodes() + rhsNumNodes);

    //copy mesh nodes
    for (auto n = numNodes; n < numNodes + rhsNumNodes; ++n)
    {
        const auto index = n - numNodes;
        m_nodes[n] = rhs.m_nodes[index];
    }

    //copy mesh edges
    for (auto e = numEdges; e < numEdges + rhsNumEdges; ++e)
    {
        const auto index = e - numEdges;
        m_edges[e].first = rhs.m_edges[index].first + numNodes;
        m_edges[e].second = rhs.m_edges[index].second + numNodes;
    }

    m_nodesRTreeRequiresUpdate = true;
    m_edgesRTreeRequiresUpdate = true;

    AdministrateNodesEdges();

    return *this;
}
