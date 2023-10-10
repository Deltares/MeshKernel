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

#include <cmath>
#include <numeric>

#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/Mesh.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Polygons.hpp>

using meshkernel::Mesh;

Mesh::Mesh(const std::vector<Edge>& edges,
           const std::vector<Point>& nodes,
           Projection projection) : m_nodes(nodes), m_edges(edges), m_projection(projection) {}

bool Mesh::NodeAdministration()
{
    // assume no duplicated links
    for (UInt e = 0; e < static_cast<UInt>(GetNumEdges()); e++)
    {
        const auto firstNode = m_edges[e].first;
        const auto secondNode = m_edges[e].second;

        if (firstNode == constants::missing::uintValue || secondNode == constants::missing::uintValue)
        {
            continue;
        }

        if (m_nodesNumEdges[firstNode] >= Mesh::m_maximumNumberOfEdgesPerNode || m_nodesNumEdges[secondNode] >= Mesh::m_maximumNumberOfEdgesPerNode)
        {
            continue;
        }

        // Search for previously connected edges
        auto alreadyAddedEdge = false;
        for (UInt i = 0; i < m_nodesNumEdges[firstNode]; ++i)
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
        for (UInt i = 0; i < m_nodesNumEdges[secondNode]; ++i)
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
    for (UInt n = 0; n < GetNumNodes(); n++)
    {
        m_nodesEdges[n].resize(m_nodesNumEdges[n]);
    }

    UInt quadrilateralCount = 0;

    for (UInt n = 0; n < GetNumNodes(); n++)
    {
        if (m_nodesNumEdges[n] == constants::geometric::numNodesInQuadrilateral)
        {
            // It is assumed that a node of a quadrilateral will have four connected edges.
            // most of the time this assumption is true.
            ++quadrilateralCount;
        }
    }

    return quadrilateralCount > GetNumNodes() / 2;
}

void Mesh::DeleteInvalidNodesAndEdges()
{

    // Mask nodes connected to valid edges
    std::vector<bool> connectedNodes(m_nodes.size(), false);
    UInt numInvalidEdges = 0;

    for (const auto& [firstNode, secondNode] : m_edges)
    {
        if (firstNode == constants::missing::uintValue || secondNode == constants::missing::uintValue)
        {
            numInvalidEdges++;
            continue;
        }
        connectedNodes[firstNode] = true;
        connectedNodes[secondNode] = true;
    }

    // Count all invalid nodes (note: there might be nodes that are not connected to an edge)
    UInt numInvalidNodes = 0;
    for (UInt n = 0; n < m_nodes.size(); ++n)
    {
        // invalidate nodes that are not connected
        if (!connectedNodes[n])
        {
            m_nodes[n] = {constants::missing::doubleValue, constants::missing::doubleValue};
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
    std::vector<UInt> validNodesIndices(m_nodes.size());
    std::ranges::fill(validNodesIndices, constants::missing::uintValue);
    UInt validIndex = 0;
    for (UInt n = 0; n < m_nodes.size(); ++n)
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
        if (firstNode != constants::missing::uintValue && secondNode != constants::missing::uintValue && validNodesIndices[firstNode] != constants::missing::uintValue && validNodesIndices[secondNode] != constants::missing::uintValue)
        {
            firstNode = validNodesIndices[firstNode];
            secondNode = validNodesIndices[secondNode];
            continue;
        }

        firstNode = constants::missing::uintValue;
        secondNode = constants::missing::uintValue;
    }

    // Remove invalid nodes, without reducing capacity
    const auto endNodeVector = std::remove_if(m_nodes.begin(), m_nodes.end(), [](const Point& n)
                                              { return !n.IsValid(); });
    m_nodes.erase(endNodeVector, m_nodes.end());

    // Remove invalid edges, without reducing capacity
    const auto endEdgeVector = std::remove_if(m_edges.begin(), m_edges.end(), [](const Edge& e)
                                              { return e.first == constants::missing::uintValue || e.second == constants::missing::uintValue; });
    m_edges.erase(endEdgeVector, m_edges.end());
}

void Mesh::MergeTwoNodes(UInt firstNodeIndex, UInt secondNodeIndex)
{
    if (firstNodeIndex >= GetNumNodes() || secondNodeIndex >= GetNumNodes())
    {
        throw std::invalid_argument("Mesh::MergeTwoNodes: Either the first or the second node-index is invalid.");
    }

    auto edgeIndex = FindEdge(firstNodeIndex, secondNodeIndex);
    if (edgeIndex != constants::missing::uintValue)
    {
        m_edges[edgeIndex].first = constants::missing::uintValue;
        m_edges[edgeIndex].second = constants::missing::uintValue;
    }

    // check if there is another edge starting at firstEdgeOtherNode and ending at secondNode
    for (UInt n = 0; n < m_nodesNumEdges[firstNodeIndex]; n++)
    {
        const auto firstEdgeIndex = m_nodesEdges[firstNodeIndex][n];
        const auto& firstEdge = m_edges[firstEdgeIndex];
        const auto firstEdgeOtherNode = OtherNodeOfEdge(firstEdge, firstNodeIndex);
        if (firstEdgeOtherNode != constants::missing::uintValue && firstEdgeOtherNode != secondNodeIndex)
        {
            for (UInt nn = 0; nn < m_nodesNumEdges[firstEdgeOtherNode]; nn++)
            {
                const auto secondEdgeIndex = m_nodesEdges[firstEdgeOtherNode][nn];
                auto secondEdge = m_edges[secondEdgeIndex];
                const auto secondNodeSecondEdge = OtherNodeOfEdge(secondEdge, firstEdgeOtherNode);
                if (secondNodeSecondEdge == secondNodeIndex)
                {
                    m_edges[secondEdgeIndex].first = constants::missing::uintValue;
                    m_edges[secondEdgeIndex].second = constants::missing::uintValue;
                }
            }
        }
    }

    // add all valid edges starting at secondNode
    std::vector<UInt> secondNodeEdges(Mesh::m_maximumNumberOfEdgesPerNode, constants::missing::uintValue);
    UInt numSecondNodeEdges = 0;
    for (UInt n = 0; n < m_nodesNumEdges[secondNodeIndex]; n++)
    {
        edgeIndex = m_nodesEdges[secondNodeIndex][n];
        if (m_edges[edgeIndex].first != constants::missing::uintValue)
        {
            secondNodeEdges[numSecondNodeEdges] = edgeIndex;
            numSecondNodeEdges++;
        }
    }

    // add all valid edges starting at firstNode are assigned to the second node
    for (UInt n = 0; n < m_nodesNumEdges[firstNodeIndex]; n++)
    {
        edgeIndex = m_nodesEdges[firstNodeIndex][n];
        if (m_edges[edgeIndex].first != constants::missing::uintValue)
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
    m_nodesEdges[secondNodeIndex] = std::vector<UInt>(secondNodeEdges.begin(), secondNodeEdges.begin() + numSecondNodeEdges);
    m_nodesNumEdges[secondNodeIndex] = numSecondNodeEdges;

    // remove edges to first node
    m_nodesEdges[firstNodeIndex] = std::vector<UInt>(0);
    m_nodesNumEdges[firstNodeIndex] = 0;
    m_nodes[firstNodeIndex] = {constants::missing::doubleValue, constants::missing::doubleValue};

    m_nodesRTreeRequiresUpdate = true;
    m_edgesRTreeRequiresUpdate = true;
}

void Mesh::MergeNodesInPolygon(const Polygons& polygon, double mergingDistance)
{
    // first filter the nodes in polygon
    std::vector<Point> filteredNodes(GetNumNodes());
    std::vector<UInt> originalNodeIndices(GetNumNodes(), constants::missing::uintValue);
    UInt index = 0;
    for (UInt i = 0; i < GetNumNodes(); ++i)
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
    for (UInt i = 0; i < filteredNodes.size(); ++i)
    {
        nodesRtree.SearchPoints(filteredNodes[i], mergingDistanceSquared);

        const auto resultSize = nodesRtree.GetQueryResultSize();
        if (resultSize > 1)
        {
            for (UInt j = 0; j < nodesRtree.GetQueryResultSize(); j++)
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

meshkernel::UInt Mesh::ConnectNodes(UInt startNode, UInt endNode)
{
    const auto edgeIndex = FindEdge(startNode, endNode);

    // The nodes are already connected
    if (edgeIndex != constants::missing::uintValue)
        return constants::missing::uintValue;

    // increment the edges container
    const auto newEdgeIndex = GetNumEdges();
    m_edges.resize(newEdgeIndex + 1);
    m_edges[newEdgeIndex].first = startNode;
    m_edges[newEdgeIndex].second = endNode;

    m_edgesRTreeRequiresUpdate = true;

    return newEdgeIndex;
}

meshkernel::UInt Mesh::InsertNode(const Point& newPoint)
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

void Mesh::DeleteNode(UInt node)
{
    if (node >= GetNumNodes())
    {
        throw std::invalid_argument("Mesh::DeleteNode: The index of the node to be deleted does not exist.");
    }

    for (UInt e = 0; e < m_nodesNumEdges[node]; e++)
    {
        const auto edgeIndex = m_nodesEdges[node][e];
        DeleteEdge(edgeIndex);
    }
    m_nodes[node] = {constants::missing::doubleValue, constants::missing::doubleValue};

    m_nodesRTreeRequiresUpdate = true;
}

void Mesh::DeleteEdge(UInt edge)
{
    if (edge == constants::missing::uintValue)
    {
        throw std::invalid_argument("Mesh::DeleteEdge: The index of the edge to be deleted does not exist.");
    }

    m_edges[edge].first = constants::missing::uintValue;
    m_edges[edge].second = constants::missing::uintValue;

    m_edgesRTreeRequiresUpdate = true;
}

void Mesh::ComputeEdgesLengths()
{
    auto const numEdges = GetNumEdges();
    m_edgeLengths.resize(numEdges, constants::missing::doubleValue);
    for (UInt e = 0; e < numEdges; e++)
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

meshkernel::UInt Mesh::FindCommonNode(UInt firstEdgeIndex, UInt secondEdgeIndex) const
{
    const auto firstEdgeFirstNode = m_edges[firstEdgeIndex].first;
    const auto firstEdgeEdgeSecondNode = m_edges[firstEdgeIndex].second;

    const auto secondEdgeFirstNode = m_edges[secondEdgeIndex].first;
    const auto secondEdgeSecondNode = m_edges[secondEdgeIndex].second;

    if (firstEdgeFirstNode == constants::missing::uintValue || firstEdgeEdgeSecondNode == constants::missing::uintValue || secondEdgeFirstNode == constants::missing::uintValue || secondEdgeSecondNode == constants::missing::uintValue)
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
    return constants::missing::uintValue;
}

meshkernel::UInt Mesh::FindEdge(UInt firstNodeIndex, UInt secondNodeIndex) const
{
    if (firstNodeIndex == constants::missing::uintValue || secondNodeIndex == constants::missing::uintValue)
    {
        throw std::invalid_argument("Mesh::FindEdge: Invalid node index.");
    }

    UInt edgeIndex = constants::missing::uintValue;
    for (UInt n = 0; n < m_nodesNumEdges[firstNodeIndex]; n++)
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

meshkernel::UInt Mesh::FindNodeCloseToAPoint(Point const& point, double searchRadius)
{
    if (GetNumNodes() <= 0)
    {
        throw std::invalid_argument("Mesh::FindNodeCloseToAPoint: There are no valid nodes.");
    }

    SearchNearestLocation(point, searchRadius * searchRadius, Location::Nodes);

    if (GetNumLocations(Location::Nodes) > 0)
    {
        return GetLocationsIndices(0, Location::Nodes);
    }

    throw AlgorithmError("Could not find the node index close to a point.");
}

meshkernel::UInt Mesh::FindNodeCloseToAPoint(Point point, const std::vector<bool>& oneDNodeMask)
{
    if (GetNumNodes() <= 0)
    {
        throw std::invalid_argument("Mesh::FindNodeCloseToAPoint: There are no valid nodes.");
    }

    SearchNearestLocation(point, Location::Nodes);

    if (GetNumLocations(Location::Nodes) <= 0)
    {
        throw AlgorithmError("Query result size <= 0.");
    }

    // resultSize > 0, no node mask applied
    if (oneDNodeMask.empty())
    {
        return GetLocationsIndices(0, Location::Nodes);
    }

    // resultSize > 0, a mask is applied
    for (UInt index = 0; index < GetNumLocations(Location::Nodes); ++index)
    {
        const auto nodeIndex = GetLocationsIndices(index, Location::Nodes);
        if (oneDNodeMask[nodeIndex])
        {
            return nodeIndex;
        }
    }

    throw AlgorithmError("Could not find the node index close to a point.");
}

meshkernel::UInt Mesh::FindEdgeCloseToAPoint(Point point)
{
    if (GetNumEdges() == 0)
    {
        throw std::invalid_argument("Mesh::GetNodeIndex: There are no valid edges.");
    }

    SearchNearestLocation(point, Location::Edges);

    if (GetNumLocations(Location::Edges) >= 1)
    {
        return GetLocationsIndices(0, Location::Edges);
    }

    throw AlgorithmError("Could not find the closest edge to a point.");
}

void Mesh::MoveNode(Point newPoint, UInt nodeindex)
{
    if (nodeindex >= m_nodes.size())
    {
        throw ConstraintError("Invalid node index: {}", nodeindex);
    }

    const Point nodeToMove = m_nodes[nodeindex];
    const auto dx = GetDx(nodeToMove, newPoint, m_projection);
    const auto dy = GetDy(nodeToMove, newPoint, m_projection);

    const double distanceNodeToMoveFromNewPointSquared = dx * dx + dy * dy;
    const double distanceNodeToMoveFromNewPointSquaredInv = 1.0 / distanceNodeToMoveFromNewPointSquared;

    for (UInt n = 0; n < GetNumNodes(); ++n)
    {
        const auto nodeDx = GetDx(m_nodes[n], nodeToMove, m_projection);
        const auto nodeDy = GetDy(m_nodes[n], nodeToMove, m_projection);
        const double distanceCurrentNodeFromNewPointSquared = nodeDx * nodeDx + nodeDy * nodeDy;

        if (distanceCurrentNodeFromNewPointSquared <= distanceNodeToMoveFromNewPointSquared)
        {
            const auto factor = 0.5 * (1.0 + std::cos(std::sqrt(distanceCurrentNodeFromNewPointSquared * distanceNodeToMoveFromNewPointSquaredInv) * M_PI));

            m_nodes[n].x += dx * factor;
            m_nodes[n].y += dy * factor;
        }
    }

    m_nodesRTreeRequiresUpdate = true;
    m_edgesRTreeRequiresUpdate = true;
}

bool Mesh::IsFaceOnBoundary(UInt face) const
{

    bool isFaceOnBoundary = false;

    for (UInt e = 0; e < GetNumFaceEdges(face); ++e)
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

void Mesh::SortEdgesInCounterClockWiseOrder(UInt startNode, UInt endNode)
{

    std::vector<double> edgeAngles(m_maximumNumberOfEdgesPerNode);
    std::vector<UInt> indices(m_maximumNumberOfEdgesPerNode);
    std::vector<UInt> edgeNodeCopy(m_maximumNumberOfEdgesPerNode);
    for (UInt n = startNode; n <= endNode; n++)
    {
        if (!m_nodes[n].IsValid())
        {
            continue;
        }

        double phi0 = 0.0;
        double phi;
        std::ranges::fill(edgeAngles, 0.0);
        for (UInt edgeIndex = 0; edgeIndex < m_nodesNumEdges[n]; edgeIndex++)
        {

            auto firstNode = m_edges[m_nodesEdges[n][edgeIndex]].first;
            auto secondNode = m_edges[m_nodesEdges[n][edgeIndex]].second;
            if (firstNode == constants::missing::uintValue || secondNode == constants::missing::uintValue)
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
            if (std::abs(deltaX) < m_minimumDeltaCoordinate && std::abs(deltaY) < m_minimumDeltaCoordinate)
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
                phi = std::atan2(deltaY, deltaX);
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
        sort(indices.begin(), indices.end(), [&](UInt const& i1, UInt const& i2)
             { return edgeAngles[i1] < edgeAngles[i2]; });

        for (UInt edgeIndex = 0; edgeIndex < m_nodesNumEdges[n]; edgeIndex++)
        {
            m_nodesEdges[n][edgeIndex] = edgeNodeCopy[indices[edgeIndex]];
        }
    }
}

void Mesh::SearchNearestLocation(Point point, Location meshLocation)
{
    switch (meshLocation)
    {
    case Location::Nodes:
        m_nodesRTree.SearchNearestPoint(point);
        break;
    case Location::Edges:
        m_edgesRTree.SearchNearestPoint(point);
        break;
    case Location::Faces:
        m_facesRTree.SearchNearestPoint(point);
        break;
    case Location::Unknown:
    default:
        throw std::runtime_error("Mesh2D::SearchNearestLocation: Mesh location has not been set.");
    }
}

void Mesh::SearchNearestLocation(Point point, double squaredRadius, Location meshLocation)
{
    switch (meshLocation)
    {
    case Location::Faces:
        m_facesRTree.SearchNearestPoint(point, squaredRadius);
        break;
    case Location::Nodes:
        m_nodesRTree.SearchNearestPoint(point, squaredRadius);
        break;
    case Location::Edges:
        m_edgesRTree.SearchNearestPoint(point, squaredRadius);
        break;
    case Location::Unknown:
    default:
        throw std::runtime_error("Mesh2D::SearchNearestLocation: Mesh location has not been set.");
    }
}

void Mesh::SearchLocations(Point point, double squaredRadius, Location meshLocation)
{
    switch (meshLocation)
    {
    case Location::Faces:
        m_facesRTree.SearchPoints(point, squaredRadius);
        break;
    case Location::Nodes:
        m_nodesRTree.SearchPoints(point, squaredRadius);
        break;
    case Location::Edges:
        m_edgesRTree.SearchPoints(point, squaredRadius);
        break;
    case Location::Unknown:
    default:
        throw std::runtime_error("Mesh2D::SearchLocations: Mesh location has not been set.");
    }
}

void Mesh::BuildTree(Location meshLocation)
{
    switch (meshLocation)
    {
    case Location::Faces:
        if (m_facesRTreeRequiresUpdate)
        {
            m_facesRTree.BuildTree(m_facesCircumcenters);
            m_facesRTreeRequiresUpdate = false;
        }
        break;
    case Location::Nodes:
        if (m_nodesRTreeRequiresUpdate)
        {

            m_nodesRTree.BuildTree(m_nodes);
            m_nodesRTreeRequiresUpdate = false;
        }
        break;
    case Location::Edges:
        if (m_edgesRTreeRequiresUpdate)
        {
            ComputeEdgesCenters();
            m_edgesRTree.BuildTree(m_edgesCenters);
            m_edgesRTreeRequiresUpdate = false;
        }
        break;
    case Location::Unknown:
    default:
        throw std::runtime_error("Mesh2D::SearchLocations: Mesh location has not been set.");
    }
}

void Mesh::BuildTree(Location meshLocation, const BoundingBox& boundingBox)
{
    switch (meshLocation)
    {
    case Location::Faces:
        if (m_facesRTreeRequiresUpdate || m_boundingBoxCache != boundingBox)
        {
            m_facesRTree.BuildTree(m_facesCircumcenters, boundingBox);
            m_facesRTreeRequiresUpdate = false;
            m_boundingBoxCache = boundingBox;
        }
        break;
    case Location::Nodes:
        if (m_nodesRTreeRequiresUpdate || m_boundingBoxCache != boundingBox)
        {
            m_nodesRTree.BuildTree(m_nodes, boundingBox);
            m_nodesRTreeRequiresUpdate = false;
            m_boundingBoxCache = boundingBox;
        }
        break;
    case Location::Edges:
        if (m_edgesRTreeRequiresUpdate || m_boundingBoxCache != boundingBox)
        {
            ComputeEdgesCenters();
            m_edgesRTree.BuildTree(m_edgesCenters, boundingBox);
            m_edgesRTreeRequiresUpdate = false;
            m_boundingBoxCache = boundingBox;
        }
        break;
    case Location::Unknown:
    default:
        throw std::runtime_error("Invalid location");
    }
}

meshkernel::UInt Mesh::GetNumLocations(Location meshLocation) const
{
    switch (meshLocation)
    {
    case Location::Faces:
        return m_facesRTree.GetQueryResultSize();
    case Location::Nodes:
        return m_nodesRTree.GetQueryResultSize();
    case Location::Edges:
        return m_edgesRTree.GetQueryResultSize();
    case Location::Unknown:
    default:
        return constants::missing::uintValue;
    }
}

meshkernel::UInt Mesh::GetLocationsIndices(UInt index, Location meshLocation)
{
    switch (meshLocation)
    {
    case Location::Faces:
        return m_facesRTree.GetQueryResult(index);
    case Location::Nodes:
        return m_nodesRTree.GetQueryResult(index);
    case Location::Edges:
        return m_edgesRTree.GetQueryResult(index);
    case Location::Unknown:
    default:
        return constants::missing::uintValue;
    }
}

void Mesh::Administrate()
{
    AdministrateNodesEdges();
}

void Mesh::AdministrateNodesEdges()
{
    DeleteInvalidNodesAndEdges();

    // return if there are no nodes or no edges
    if (m_nodes.empty() || m_edges.empty())
    {
        return;
    }

    m_nodesEdges.resize(m_nodes.size());
    std::ranges::fill(m_nodesEdges, std::vector(m_maximumNumberOfEdgesPerNode, constants::missing::uintValue));

    m_nodesNumEdges.resize(m_nodes.size());
    std::ranges::fill(m_nodesNumEdges, 0);

    NodeAdministration();

    SortEdgesInCounterClockWiseOrder(0, GetNumNodes() - 1);
}

double Mesh::ComputeMaxLengthSurroundingEdges(UInt node)
{

    if (m_edgeLengths.empty())
    {
        ComputeEdgesLengths();
    }

    auto maxEdgeLength = std::numeric_limits<double>::lowest();
    for (UInt ee = 0; ee < m_nodesNumEdges[node]; ++ee)
    {
        const auto edge = m_nodesEdges[node][ee];
        maxEdgeLength = std::max(maxEdgeLength, m_edgeLengths[edge]);
    }

    return maxEdgeLength;
}

std::vector<meshkernel::Point> Mesh::ComputeLocations(Location location) const
{
    std::vector<Point> result;
    if (location == Location::Nodes)
    {
        result.reserve(GetNumNodes());
        for (const auto& n : m_nodes)
        {
            result.emplace_back(n);
        }
    }
    if (location == Location::Edges)
    {
        result.reserve(GetNumEdges());
        for (const auto& [firstNode, secondNode] : m_edges)
        {

            if (firstNode != constants::missing::uintValue && secondNode != constants::missing::uintValue)
            {
                result.emplace_back((m_nodes[firstNode] + m_nodes[secondNode]) * 0.5);
            }
            else
            {
                result.emplace_back(Point{});
            }
        }
    }
    if (location == Location::Faces)
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
    if (m_projection != rhs.m_projection)
    {
        throw std::invalid_argument("Mesh2D::operator+=: The two meshes cannot be added because they have different projections");
    }

    if (rhs.GetNumNodes() == 0 || rhs.GetNumEdges() == 0)
    {
        return *this;
    }

    const auto rhsNumNodes = rhs.GetNumNodes();
    const auto rhsNumEdges = rhs.GetNumEdges();

    auto const numNodes = GetNumNodes();
    auto const numEdges = GetNumEdges();
    m_edges.resize(GetNumEdges() + rhsNumEdges);
    m_nodes.resize(GetNumNodes() + rhsNumNodes);

    // copy mesh nodes
    for (auto n = numNodes; n < numNodes + rhsNumNodes; ++n)
    {
        const auto index = n - numNodes;
        m_nodes[n] = rhs.m_nodes[index];
    }

    // copy mesh edges
    for (auto e = numEdges; e < numEdges + rhsNumEdges; ++e)
    {
        const auto index = e - numEdges;
        m_edges[e].first = rhs.m_edges[index].first + numNodes;
        m_edges[e].second = rhs.m_edges[index].second + numNodes;
    }

    m_nodesRTreeRequiresUpdate = true;
    m_edgesRTreeRequiresUpdate = true;

    Administrate();

    return *this;
}
