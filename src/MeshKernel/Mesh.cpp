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

#pragma once
#include "MeshKernel/Mesh.hpp"

#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/Operations.hpp>

#include <MeshKernel/Entities.hpp>
#include <stdexcept>
#include <vector>

meshkernel::Mesh::Mesh(const std::vector<Edge>& edges,
                       const std::vector<Point>& nodes,
                       Projection projection) : m_nodes(nodes), m_edges(edges), m_projection(projection){};

void meshkernel::Mesh::NodeAdministration()
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

void meshkernel::Mesh::DeleteInvalidNodesAndEdges()
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
        m_numNodes = m_nodes.size();
        m_numEdges = m_edges.size();
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
    m_numNodes = m_nodes.size();

    // Remove invalid edges, without reducing capacity
    const auto endEdgeVector = std::remove_if(m_edges.begin(), m_edges.end(), [](const Edge& e) { return e.first == sizetMissingValue || e.second == sizetMissingValue; });
    m_edges.erase(endEdgeVector, m_edges.end());
    m_numEdges = m_edges.size();
}

void meshkernel::Mesh::MergeTwoNodes(size_t firstNodeIndex, size_t secondNodeIndex)
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

size_t meshkernel::Mesh::ConnectNodes(size_t startNode, size_t endNode)
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
    m_numEdges++;

    m_edgesRTreeRequiresUpdate = true;

    return newEdgeIndex;
}

size_t meshkernel::Mesh::InsertNode(const Point& newPoint)
{
    const auto newSize = GetNumNodes() + 1;
    const auto newNodeIndex = GetNumNodes();

    m_nodes.resize(newSize);
    m_nodeMask.resize(newSize);
    m_nodesNumEdges.resize(newSize);
    m_nodesEdges.resize(newSize);

    m_numNodes++;

    m_nodes[newNodeIndex] = newPoint;
    m_nodeMask[newNodeIndex] = static_cast<int>(newNodeIndex);
    m_nodesNumEdges[newNodeIndex] = 0;

    m_nodesRTreeRequiresUpdate = true;

    return newNodeIndex;
}

void meshkernel::Mesh::DeleteNode(size_t nodeIndex)
{
    if (nodeIndex >= GetNumNodes())
    {
        throw std::invalid_argument("Mesh::DeleteNode: The index of the node to be deleted does not exist.");
    }

    for (auto e = 0; e < m_nodesNumEdges[nodeIndex]; e++)
    {
        const auto edgeIndex = m_nodesEdges[nodeIndex][e];
        DeleteEdge(edgeIndex);
    }
    m_nodes[nodeIndex] = {doubleMissingValue, doubleMissingValue};
    m_numNodes--;

    m_nodesRTreeRequiresUpdate = true;
}

void meshkernel::Mesh::DeleteEdge(size_t edgeIndex)
{
    if (edgeIndex == sizetMissingValue)
    {
        throw std::invalid_argument("Mesh::DeleteEdge: The index of the edge to be deleted does not exist.");
    }

    m_edges[edgeIndex].first = sizetMissingValue;
    m_edges[edgeIndex].second = sizetMissingValue;

    m_edgesRTreeRequiresUpdate = true;
}

void meshkernel::Mesh::ComputeEdgesLengths()
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

void meshkernel::Mesh::ComputeEdgesCenters()
{
    m_edgesCenters = ComputeEdgeCenters(m_nodes, m_edges);
}

size_t meshkernel::Mesh::FindCommonNode(size_t firstEdgeIndex, size_t secondEdgeIndex) const
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

size_t meshkernel::Mesh::FindEdge(size_t firstNodeIndex, size_t secondNodeIndex) const
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

size_t meshkernel::Mesh::FindNodeCloseToAPoint(Point point, double searchRadius)
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

    double const searchRadiusSquared = searchRadius * searchRadius;
    m_nodesRTree.NearestNeighborsOnSquaredDistance(point, searchRadiusSquared);
    const auto resultSize = m_nodesRTree.GetQueryResultSize();

    if (resultSize > 0)
    {
        return m_nodesRTree.GetQueryIndex(0);
    }

    throw AlgorithmError("Mesh::FindNodeCloseToAPoint: Could not find the node index close to a point.");
}

size_t meshkernel::Mesh::FindNodeCloseToAPoint(Point point, const std::vector<bool>& oneDNodeMask)
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

    m_nodesRTree.NearestNeighbors(point);
    const auto resultSize = m_nodesRTree.GetQueryResultSize();

    // no results found
    if (resultSize <= 0)
    {
        throw AlgorithmError("Mesh::FindNodeCloseToAPoint: query result size <= 0.");
    }

    // resultSize > 0, no node mask applied
    if (oneDNodeMask.empty())
    {
        return m_nodesRTree.GetQueryIndex(0);
    }

    // resultSize > 0, a mask is applied
    for (auto index = 0; index < resultSize; ++index)
    {
        const auto nodeIndex = m_nodesRTree.GetQueryIndex(index);
        if (oneDNodeMask[nodeIndex])
        {
            return nodeIndex;
        }
    }

    throw AlgorithmError("Mesh::FindNodeCloseToAPoint: Could not find the node index close to a point.");
}

size_t meshkernel::Mesh::FindEdgeCloseToAPoint(Point point)
{
    if (GetNumEdges() == 0)
    {
        throw std::invalid_argument("Mesh::GetNodeIndex: There are no valid edges.");
    }

    if (m_edgesRTree.Empty())
    {
        ComputeEdgesCenters();
        m_edgesRTree.BuildTree(m_edgesCenters);
        m_edgesRTreeRequiresUpdate = false;
    }

    m_edgesRTree.NearestNeighbors(point);
    auto const resultSize = m_edgesRTree.GetQueryResultSize();
    if (resultSize >= 1)
    {
        const auto edgeIndex = m_edgesRTree.GetQueryIndex(0);
        return edgeIndex;
    }

    throw AlgorithmError("Mesh::FindEdgeCloseToAPoint: Could not find the closest edge to a point.");
}

void meshkernel::Mesh::MoveNode(Point newPoint, size_t nodeindex)
{
    const Point nodeToMove = m_nodes[nodeindex];

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

bool meshkernel::Mesh::IsFaceOnBoundary(size_t face) const
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

void meshkernel::Mesh::SortEdgesInCounterClockWiseOrder(size_t node)
{
    if (!m_nodes[node].IsValid())
    {
        throw std::invalid_argument("Mesh::SortEdgesInCounterClockWiseOrder: Invalid nodes.");
    }

    double phi0 = 0.0;
    double phi;
    std::vector<double> edgeAngles(maximumNumberOfEdgesPerNode);
    std::fill(edgeAngles.begin(), edgeAngles.end(), 0.0);
    for (auto edgeIndex = 0; edgeIndex < m_nodesNumEdges[node]; edgeIndex++)
    {

        auto firstNode = m_edges[m_nodesEdges[node][edgeIndex]].first;
        auto secondNode = m_edges[m_nodesEdges[node][edgeIndex]].second;
        if (firstNode == sizetMissingValue || secondNode == sizetMissingValue)
        {
            continue;
        }

        if (secondNode == node)
        {
            secondNode = firstNode;
            firstNode = node;
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
    std::vector<std::size_t> indices(m_nodesNumEdges[node]);
    std::vector<size_t> edgeNodeCopy{m_nodesEdges[node]};
    iota(indices.begin(), indices.end(), 0);
    sort(indices.begin(), indices.end(), [&](std::size_t i1, std::size_t i2) { return edgeAngles[i1] < edgeAngles[i2]; });

    for (std::size_t edgeIndex = 0; edgeIndex < m_nodesNumEdges[node]; edgeIndex++)
    {
        m_nodesEdges[node][edgeIndex] = edgeNodeCopy[indices[edgeIndex]];
    }
}

void meshkernel::Mesh::AdministrateNodesEdges()
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
    if (m_numNodes == 0 || m_numEdges == 0)
    {
        return;
    }

    m_nodesEdges.resize(m_nodes.size());
    std::fill(m_nodesEdges.begin(), m_nodesEdges.end(), std::vector<size_t>(maximumNumberOfEdgesPerNode, sizetMissingValue));

    m_nodesNumEdges.resize(m_nodes.size());
    std::fill(m_nodesNumEdges.begin(), m_nodesNumEdges.end(), 0);

    NodeAdministration();

    for (auto n = 0; n < GetNumNodes(); n++)
    {
        SortEdgesInCounterClockWiseOrder(n);
    }
}
