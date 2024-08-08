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

#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Mesh.hpp"
#include "MeshKernel/Operations.hpp"
#include "MeshKernel/Polygons.hpp"
#include "MeshKernel/RangeCheck.hpp"
#include "MeshKernel/UndoActions/CompoundUndoAction.hpp"
#include "MeshKernel/Utilities/RTreeFactory.hpp"

using meshkernel::Mesh;

Mesh::Mesh() : Mesh(Projection::cartesian)
{
}

Mesh::Mesh(Projection projection) : m_projection(projection)
{

    m_RTrees.emplace(Location::Nodes, RTreeFactory::Create(m_projection));
    m_RTrees.emplace(Location::Edges, RTreeFactory::Create(m_projection));
    m_RTrees.emplace(Location::Faces, RTreeFactory::Create(m_projection));
}

Mesh::Mesh(const std::vector<Edge>& edges,
           const std::vector<Point>& nodes,
           Projection projection) : m_projection(projection),
                                    m_nodes(nodes),
                                    m_edges(edges)

{
    m_RTrees.emplace(Location::Nodes, RTreeFactory::Create(m_projection));
    m_RTrees.emplace(Location::Edges, RTreeFactory::Create(m_projection));
    m_RTrees.emplace(Location::Faces, RTreeFactory::Create(m_projection));
    DeleteInvalidNodesAndEdges();
}

bool Mesh::NodeAdministration()
{
    // assume no duplicated links
    for (UInt e = 0; e < GetNumEdges(); e++)
    {
        const auto firstNode = m_edges[e].first;
        const auto secondNode = m_edges[e].second;

        if (firstNode == constants::missing::uintValue || secondNode == constants::missing::uintValue)
        {
            continue;
        }

        if (m_nodesNumEdges[firstNode] >= m_maximumNumberOfEdgesPerNode || m_nodesNumEdges[secondNode] >= m_maximumNumberOfEdgesPerNode)
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

void Mesh::FindConnectedNodes(std::vector<bool>& connectedNodes, UInt& numInvalidEdges) const
{
    numInvalidEdges = 0;

    connectedNodes.resize(m_nodes.size());
    std::fill(connectedNodes.begin(), connectedNodes.end(), false);

    for (const auto& [firstNode, secondNode] : m_edges)
    {
        if (firstNode == constants::missing::uintValue ||
            secondNode == constants::missing::uintValue ||
            !m_nodes[firstNode].IsValid() ||
            !m_nodes[secondNode].IsValid())
        {
            numInvalidEdges++;
        }
        else
        {
            connectedNodes[firstNode] = true;
            connectedNodes[secondNode] = true;
        }
    }
}

void Mesh::InvalidateUnConnectedNodes(const std::vector<bool>& connectedNodes, UInt& numInvalidNodes, CompoundUndoAction* undoAction)
{
    numInvalidNodes = 0;

    for (UInt n = 0; n < m_nodes.size(); ++n)
    {
        // invalidate nodes that are not connected
        if (!connectedNodes[n])
        {
            if (undoAction == nullptr)
            {
                m_nodes[n].SetInvalid();
            }
            else
            {
                undoAction->Add(ResetNode(n, {constants::missing::doubleValue, constants::missing::doubleValue}));
            }
        }

        if (!m_nodes[n].IsValid())
        {
            numInvalidNodes++;
        }
    }
}

void Mesh::DeleteInvalidNodesAndEdges()
{

    // Mask nodes connected to valid edges
    std::vector<bool> connectedNodes(m_nodes.size(), false);
    UInt numInvalidEdges = 0;
    // Count all invalid nodes (note: there might be nodes that are not connected to an edge)
    UInt numInvalidNodes = 0;

    FindConnectedNodes(connectedNodes, numInvalidEdges);
    InvalidateUnConnectedNodes(connectedNodes, numInvalidNodes);

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
        if (firstNode != constants::missing::uintValue &&
            secondNode != constants::missing::uintValue &&
            validNodesIndices[firstNode] != constants::missing::uintValue &&
            validNodesIndices[secondNode] != constants::missing::uintValue)
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

void Mesh::SetUnConnectedNodesAndEdgesToInvalid(CompoundUndoAction* undoAction)
{
    // Mask nodes connected to valid edges
    std::vector<bool> connectedNodes(m_nodes.size(), false);
    UInt numInvalidEdges = 0;
    // Count all invalid nodes (note: there might be nodes that are not connected to an edge)
    UInt numInvalidNodes = 0;

    FindConnectedNodes(connectedNodes, numInvalidEdges);
    InvalidateUnConnectedNodes(connectedNodes, numInvalidNodes, undoAction);

    // If there is nothing to invalidate then return
    if (numInvalidEdges == 0 && numInvalidNodes == 0)
    {
        return;
    }

    // Flag invalid nodes
    std::vector<bool> nodeIsValid(m_nodes.size());

    for (UInt n = 0; n < m_nodes.size(); ++n)
    {
        nodeIsValid[n] = m_nodes[n].IsValid();
    }

    // Flag invalid edges
    for (UInt e = 0; e < m_edges.size(); ++e)
    {
        Edge& edge = m_edges[e];

        if (edge.first == constants::missing::uintValue ||
            edge.second == constants::missing::uintValue ||
            !nodeIsValid[edge.first] || !nodeIsValid[edge.second])
        {
            if (undoAction == nullptr)
            {
                edge = {constants::missing::uintValue, constants::missing::uintValue};
            }
            else
            {
                undoAction->Add(ResetEdge(e, {constants::missing::uintValue, constants::missing::uintValue}));
            }
        }
    }
}

std::unique_ptr<meshkernel::UndoAction> Mesh::MergeTwoNodes(UInt firstNodeIndex, UInt secondNodeIndex)
{
    if (firstNodeIndex == constants::missing::uintValue)
    {
        throw ConstraintError("The first node-index has the null value");
    }

    if (secondNodeIndex == constants::missing::uintValue)
    {
        throw ConstraintError("The second node-index has the null value");
    }

    if (firstNodeIndex >= GetNumNodes())
    {
        throw RangeError("The first node-index is out of range: value = {}, number of nodes = {}", firstNodeIndex, GetNumNodes());
    }

    if (secondNodeIndex >= GetNumNodes())
    {
        throw RangeError("The second node-index is out of range: value = {}, number of nodes = {}", secondNodeIndex, GetNumNodes());
    }

    if (firstNodeIndex == secondNodeIndex)
    {
        // Nothing to do if the two nodes have the same index.
        return nullptr;
    }

    std::unique_ptr<CompoundUndoAction> undoAction = CompoundUndoAction::Create();

    auto edgeIndex = FindEdge(firstNodeIndex, secondNodeIndex);

    if (edgeIndex != constants::missing::uintValue)
    {
        undoAction->Add(DeleteEdge(edgeIndex));
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
                    undoAction->Add(DeleteEdge(secondEdgeIndex));
                }
            }
        }
    }

    // add all valid edges starting at secondNode
    std::vector<UInt> secondNodeEdges(m_maximumNumberOfEdgesPerNode, constants::missing::uintValue);
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
                undoAction->Add(ResetEdge(edgeIndex, {secondNodeIndex, m_edges[edgeIndex].second}));
            }
            else if (m_edges[edgeIndex].second == firstNodeIndex)
            {
                undoAction->Add(ResetEdge(edgeIndex, {m_edges[edgeIndex].first, secondNodeIndex}));
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

    // Set the node to be invalid
    undoAction->Add(ResetNode(firstNodeIndex, {constants::missing::doubleValue, constants::missing::doubleValue}));

    m_nodesRTreeRequiresUpdate = true;
    m_edgesRTreeRequiresUpdate = true;
    return undoAction;
}

std::unique_ptr<meshkernel::UndoAction> Mesh::MergeNodesInPolygon(const Polygons& polygon, double mergingDistance)
{

    // first filter the nodes in polygon
    const auto numNodes = GetNumNodes();
    std::vector<Point> filteredNodes(numNodes);
    std::vector<UInt> originalNodeIndices(numNodes, constants::missing::uintValue);
    const auto isNodeInPolygon = IsLocationInPolygon(polygon, Location::Nodes);

    UInt filteredNodeCount = 0;
    for (UInt i = 0; i < numNodes; ++i)
    {
        if (isNodeInPolygon[i])
        {
            filteredNodes[filteredNodeCount] = m_nodes[i];
            originalNodeIndices[filteredNodeCount] = i;
            filteredNodeCount++;
        }
    }

    // no node to merge
    if (filteredNodeCount == 0)
    {
        return nullptr;
    }

    std::unique_ptr<CompoundUndoAction> undoAction = CompoundUndoAction::Create();
    filteredNodes.resize(filteredNodeCount);

    AdministrateNodesEdges();

    // Update the R-Tree of the mesh nodes
    const auto nodesRtree = RTreeFactory::Create(m_projection);
    nodesRtree->BuildTree(filteredNodes);

    // merge the closest nodes
    auto const mergingDistanceSquared = mergingDistance * mergingDistance;
    for (UInt i = 0; i < filteredNodes.size(); ++i)
    {
        nodesRtree->SearchPoints(filteredNodes[i], mergingDistanceSquared);

        const auto resultSize = nodesRtree->GetQueryResultSize();
        if (resultSize > 1)
        {
            for (UInt j = 0; j < nodesRtree->GetQueryResultSize(); j++)
            {
                const auto nodeIndexInFilteredNodes = nodesRtree->GetQueryResult(j);
                if (nodeIndexInFilteredNodes != i)
                {
                    undoAction->Add(MergeTwoNodes(originalNodeIndices[i], originalNodeIndices[nodeIndexInFilteredNodes]));
                    nodesRtree->DeleteNode(i);
                }
            }
        }
    }

    AdministrateNodesEdges();
    return undoAction;
}

std::tuple<meshkernel::UInt, std::unique_ptr<meshkernel::AddNodeAction>> Mesh::InsertNode(const Point& newPoint)
{
    const auto newNodeIndex = GetNumNodes();

    m_nodes.resize(newNodeIndex + 1);
    m_nodesNumEdges.resize(newNodeIndex + 1);
    m_nodesEdges.resize(newNodeIndex + 1);

    std::unique_ptr<AddNodeAction> undoAction = AddNodeAction::Create(*this, newNodeIndex, newPoint);
    CommitAction(*undoAction);

    return {newNodeIndex, std::move(undoAction)};
}

std::tuple<meshkernel::UInt, std::unique_ptr<meshkernel::AddEdgeAction>> Mesh::ConnectNodes(UInt startNode, UInt endNode)
{
    if (FindEdge(startNode, endNode) != constants::missing::uintValue)
    {
        return {constants::missing::uintValue, nullptr};
    }

    // increment the edges container
    const auto newEdgeIndex = GetNumEdges();
    m_edges.resize(newEdgeIndex + 1);

    std::unique_ptr<AddEdgeAction> undoAction = AddEdgeAction::Create(*this, newEdgeIndex, startNode, endNode);
    CommitAction(*undoAction);
    return {newEdgeIndex, std::move(undoAction)};
}

std::unique_ptr<meshkernel::ResetNodeAction> Mesh::ResetNode(const UInt nodeId, const Point& newValue)
{
    if (nodeId >= GetNumNodes())
    {
        throw ConstraintError("The node index, {}, is not in range.", nodeId);
    }

    std::unique_ptr<ResetNodeAction> undoAction = ResetNodeAction::Create(*this, nodeId, m_nodes[nodeId], newValue);
    CommitAction(*undoAction);
    return undoAction;
}

void Mesh::SetNode(const UInt nodeId, const Point& newValue)
{
    if (nodeId >= GetNumNodes())
    {
        throw ConstraintError("The node index, {}, is not in range.", nodeId);
    }

    m_nodes[nodeId] = newValue;
}

std::unique_ptr<meshkernel::DeleteEdgeAction> Mesh::DeleteEdge(UInt edge)
{
    if (edge == constants::missing::uintValue) [[unlikely]]
    {
        throw std::invalid_argument("Mesh::DeleteEdge: The index of the edge to be deleted does not exist.");
    }

    std::unique_ptr<meshkernel::DeleteEdgeAction> undoAction = DeleteEdgeAction::Create(*this, edge, m_edges[edge].first, m_edges[edge].second);

    CommitAction(*undoAction);
    return undoAction;
}

std::unique_ptr<meshkernel::DeleteNodeAction> Mesh::DeleteNode(UInt node)
{
    if (node >= GetNumNodes()) [[unlikely]]
    {
        throw std::invalid_argument("Mesh::DeleteNode: The index of the node to be deleted does not exist.");
    }

    std::unique_ptr<DeleteNodeAction> undoAction = DeleteNodeAction::Create(*this, node, m_nodes[node]);

    for (UInt e = 0; e < m_nodesEdges[node].size(); e++)
    {
        const auto edgeIndex = m_nodesEdges[node][e];
        undoAction->Add(DeleteEdge(edgeIndex));
    }

    CommitAction(*undoAction);
    return undoAction;
}

void Mesh::ComputeEdgesLengths()
{
    auto const numEdges = GetNumEdges();
    m_edgeLengths.resize(numEdges, constants::missing::doubleValue);

    // TODO could be openmp loop
    for (UInt e = 0; e < numEdges; e++)
    {
        auto const first = m_edges[e].first;
        auto const second = m_edges[e].second;

        if (first != constants::missing::uintValue && second != constants::missing::uintValue) [[likely]]
        {
            m_edgeLengths[e] = ComputeDistance(m_nodes[first], m_nodes[second], m_projection);
        }
    }
}

double Mesh::ComputeMinEdgeLength(const Polygons& polygon) const
{
    auto const numEdges = GetNumEdges();
    auto result = std::numeric_limits<double>::max();

    const auto isNodeInPolygon = IsLocationInPolygon(polygon, Location::Nodes);
    for (UInt e = 0; e < numEdges; e++)
    {
        const auto& [firstNode, secondNode] = m_edges[e];
        if (isNodeInPolygon[firstNode] || isNodeInPolygon[secondNode])
        {
            result = std::min(result, m_edgeLengths[e]);
        }
    }
    return result;
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
        throw ConstraintError("Mesh::FindEdge: Invalid node index: first {}, second {}", firstNodeIndex, secondNodeIndex);
    }

    for (UInt n = 0; n < m_nodesNumEdges[firstNodeIndex]; n++)
    {
        const auto edgeIndex = m_nodesEdges[firstNodeIndex][n];
        const auto firstEdgeOtherNode = OtherNodeOfEdge(m_edges[edgeIndex], firstNodeIndex);
        const auto edgeFound = firstEdgeOtherNode == secondNodeIndex;
        if (edgeFound)
        {
            return edgeIndex;
        }
    }

    return constants::missing::uintValue;
}

meshkernel::UInt Mesh::FindEdgeWithLinearSearch(UInt firstNodeIndex, UInt secondNodeIndex) const
{
    for (UInt edgeIndex = 0; edgeIndex < GetNumEdges(); edgeIndex++)
    {
        const auto& [firstNode, secondNode] = m_edges[edgeIndex];
        const auto edgeFound = (firstNode == firstNodeIndex && secondNode == secondNodeIndex) ||
                               (secondNode == firstNodeIndex && firstNode == secondNodeIndex);
        if (edgeFound)
        {
            return edgeIndex;
        }
    }

    return constants::missing::uintValue;
}

meshkernel::UInt Mesh::FindNodeCloseToAPoint(Point const& point, double searchRadius)
{
    if (GetNumNodes() <= 0)
    {
        return constants::missing::uintValue;
    }
    BuildTree(Location::Nodes);

    const auto& rtree = m_RTrees.at(Location::Nodes);

    rtree->SearchNearestPoint(point, searchRadius * searchRadius);

    if (rtree->GetQueryResultSize() > 0)
    {
        return rtree->GetQueryResult(0);
    }

    return constants::missing::uintValue;
}

meshkernel::UInt Mesh::FindLocationIndex(Point point,
                                         Location location,
                                         const std::vector<bool>& locationMask,
                                         const BoundingBox& boundingBox)
{
    BuildTree(location, boundingBox);
    const auto& rtree = m_RTrees.at(location);
    if (rtree->Empty())
    {
        throw AlgorithmError("Empty RTree");
    }

    rtree->SearchNearestPoint(point);
    const auto numLocations = rtree->GetQueryResultSize();

    if (numLocations <= 0)
    {
        throw AlgorithmError("Query result size <= 0.");
    }

    // resultSize > 0, no node mask applied
    if (locationMask.empty())
    {
        return rtree->GetQueryResult(0);
    }

    // resultSize > 0, a mask is applied
    for (UInt index = 0; index < numLocations; ++index)
    {
        const auto locationIndex = rtree->GetQueryResult(index);
        if (locationMask[locationIndex])
        {
            return locationIndex;
        }
    }

    throw AlgorithmError("Could not find a valid location close to a point.");
}

std::unique_ptr<meshkernel::UndoAction> Mesh::MoveNode(Point newPoint, UInt nodeindex)
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

    std::vector<UInt> movedNodeIndex;
    std::vector<Vector> nodeDisplacement;
    movedNodeIndex.reserve(GetNumNodes());
    nodeDisplacement.reserve(GetNumNodes());

    for (UInt n = 0; n < GetNumNodes(); ++n)
    {
        const auto nodeDx = GetDx(m_nodes[n], nodeToMove, m_projection);
        const auto nodeDy = GetDy(m_nodes[n], nodeToMove, m_projection);
        const double distanceCurrentNodeFromNewPointSquared = nodeDx * nodeDx + nodeDy * nodeDy;

        if (distanceCurrentNodeFromNewPointSquared <= distanceNodeToMoveFromNewPointSquared)
        {
            const auto factor = 0.5 * (1.0 + std::cos(std::sqrt(distanceCurrentNodeFromNewPointSquared * distanceNodeToMoveFromNewPointSquaredInv) * M_PI));

            nodeDisplacement.emplace_back(dx * factor, dy * factor);
            movedNodeIndex.emplace_back(n);
        }
    }

    std::unique_ptr<NodeTranslationAction> undoAction = NodeTranslationAction::Create(*this, movedNodeIndex);

    for (UInt i = 0; i < movedNodeIndex.size(); ++i)
    {
        m_nodes[movedNodeIndex[i]] += nodeDisplacement[i];
    }

    m_nodesRTreeRequiresUpdate = true;
    m_edgesRTreeRequiresUpdate = true;
    return undoAction;
}

std::unique_ptr<meshkernel::ResetEdgeAction> Mesh::ResetEdge(UInt edgeId, const Edge& edge)
{
    std::unique_ptr<meshkernel::ResetEdgeAction> undoAction = ResetEdgeAction::Create(*this, edgeId, m_edges[edgeId], edge);
    CommitAction(*undoAction);
    return undoAction;
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

void Mesh::Administrate(CompoundUndoAction* undoAction)
{
    AdministrateNodesEdges(undoAction);
}

void Mesh::AdministrateNodesEdges(CompoundUndoAction* undoAction)
{
    SetUnConnectedNodesAndEdgesToInvalid(undoAction);

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

meshkernel::UInt Mesh::GetLocalFaceNodeIndex(const UInt faceIndex, const UInt nodeIndex) const
{
    UInt faceNodeIndex = constants::missing::uintValue;

    const auto numFaceNodes = GetNumFaceEdges(faceIndex);

    for (UInt n = 0; n < numFaceNodes; ++n)
    {
        if (m_facesNodes[faceIndex][n] == nodeIndex)
        {
            faceNodeIndex = n;
            break;
        }
    }

    return faceNodeIndex;
}

std::vector<bool> Mesh::IsLocationInPolygon(const Polygons& polygon, Location location) const
{
    const auto locations = ComputeLocations(location);
    std::vector<bool> result(locations.size(), false);
    for (UInt i = 0; i < result.size(); ++i)
    {
        result[i] = polygon.IsPointInPolygon(locations[i], 0);
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

std::unique_ptr<meshkernel::UndoAction> Mesh::Join(const Mesh& rhs)
{
    std::unique_ptr<FullUnstructuredGridUndo> joinAction = FullUnstructuredGridUndo::Create(*this);
    *this += rhs;
    return joinAction;
}

meshkernel::UInt Mesh::GetNumValidNodes() const
{
    return static_cast<UInt>(std::ranges::count_if(m_nodes, [](const Point& p)
                                                   { return p.IsValid(); }));
}

meshkernel::UInt Mesh::GetNumValidEdges() const
{
    UInt count = 0;

    for (UInt i = 0; i < m_edges.size(); ++i)
    {
        if (IsValidEdge(i))
        {
            ++count;
        }
    }

    return count;
}

meshkernel::UInt Mesh::GetEdgeIndex(const UInt elementId, const UInt edgeId) const
{
    if (elementId == constants::missing::uintValue || edgeId == constants::missing::uintValue)
    {
        return constants::missing::uintValue;
    }

    if (elementId >= GetNumFaces())
    {
        throw ConstraintError("Element id is greater than the number of elements: {} >= {}", elementId, GetNumFaces());
    }

    if (edgeId >= GetNumEdges())
    {
        throw ConstraintError("edge id is greater than the number of edges: {} >= {}", edgeId, GetNumEdges());
    }

    const std::vector<UInt>& edgeIds = m_facesEdges[elementId];

    for (UInt e = 0; e < edgeIds.size(); ++e)
    {
        if (edgeIds[e] == edgeId)
        {
            return e;
        }
    }

    return constants::missing::uintValue;
}

meshkernel::UInt Mesh::GetNodeIndex(const UInt elementId, const UInt nodeId) const
{
    if (elementId == constants::missing::uintValue || nodeId == constants::missing::uintValue)
    {
        return constants::missing::uintValue;
    }

    if (elementId >= GetNumFaces())
    {
        throw ConstraintError("Element id is greater than the number of elements: {} >= {}", elementId, GetNumFaces());
    }

    if (nodeId >= GetNumValidNodes())
    {
        throw ConstraintError("node id is greater than the number of nodes: {} >= {}", nodeId, GetNumValidNodes());
    }

    const std::vector<UInt>& nodeIds = m_facesNodes[elementId];

    // TODO use Operations::FindIndex when curvilinear grid from splines has been added to master
    // return FindIndex (m_facesNodes[elementId], nodeId);

    for (UInt n = 0; n < nodeIds.size(); ++n)
    {
        if (nodeIds[n] == nodeId)
        {
            return n;
        }
    }

    return constants::missing::uintValue;
}

std::vector<meshkernel::UInt> Mesh::GetValidNodeMapping() const
{
    std::vector<meshkernel::UInt> nodeMap(GetNumNodes());
    UInt count = 0;

    for (UInt i = 0; i < m_nodes.size(); ++i)
    {
        if (m_nodes[i].IsValid())
        {
            nodeMap[count] = i;
            ++count;
        }
    }

    nodeMap.resize(count);
    return nodeMap;
}

std::vector<meshkernel::UInt> Mesh::GetValidEdgeMapping() const
{
    std::vector<meshkernel::UInt> edgeMap(GetNumEdges());
    UInt count = 0;

    for (UInt i = 0; i < m_edges.size(); ++i)
    {
        if (IsValidEdge(i))
        {
            edgeMap[count] = i;
            ++count;
        }
    }

    edgeMap.resize(count);
    return edgeMap;
}

bool Mesh::IsValidEdge(const UInt edgeId) const
{
    if (edgeId >= m_edges.size())
    {
        throw ConstraintError("The edge index is out of bounds. {} >= {}.", edgeId, m_edges.size());
    }

    return m_edges[edgeId].first != constants::missing::uintValue && m_edges[edgeId].second != constants::missing::uintValue &&
           m_nodes[m_edges[edgeId].first].IsValid() && m_nodes[m_edges[edgeId].second].IsValid();
}

void Mesh::BuildTree(Location location, const BoundingBox& boundingBox)
{
    switch (location)
    {
    case Location::Faces:
        if (m_facesRTreeRequiresUpdate || m_boundingBoxCache != boundingBox)
        {
            Administrate();
            m_RTrees.at(Location::Faces)->BuildTree(m_facesCircumcenters, boundingBox);
            m_facesRTreeRequiresUpdate = false;
            m_boundingBoxCache = boundingBox;
        }
        break;
    case Location::Nodes:
        if (m_nodesRTreeRequiresUpdate || m_boundingBoxCache != boundingBox)
        {
            m_RTrees.at(Location::Nodes)->BuildTree(m_nodes, boundingBox);
            m_nodesRTreeRequiresUpdate = false;
            m_boundingBoxCache = boundingBox;
        }
        break;
    case Location::Edges:
        if (m_edgesRTreeRequiresUpdate || m_boundingBoxCache != boundingBox)
        {
            ComputeEdgesCenters();
            m_RTrees.at(Location::Edges)->BuildTree(m_edgesCenters, boundingBox);
            m_edgesRTreeRequiresUpdate = false;
            m_boundingBoxCache = boundingBox;
        }
        break;
    case Location::Unknown:
    default:
        throw std::runtime_error("Invalid location");
    }
}

//--------------------------------

void Mesh::CommitAction(const AddNodeAction& undoAction)
{
    m_nodes[undoAction.NodeId()] = undoAction.Node();
    m_nodesNumEdges[undoAction.NodeId()] = 0;
    m_nodesRTreeRequiresUpdate = true;
}

void Mesh::CommitAction(const AddEdgeAction& undoAction)
{
    m_edges[undoAction.EdgeId()] = undoAction.GetEdge();
    m_edgesRTreeRequiresUpdate = true;
}

void Mesh::CommitAction(const ResetNodeAction& undoAction)
{
    m_nodes[undoAction.NodeId()] = undoAction.UpdatedNode();
    m_nodesRTreeRequiresUpdate = true;
    m_edgesRTreeRequiresUpdate = true;
}

void Mesh::CommitAction(const ResetEdgeAction& undoAction)
{
    m_edges[undoAction.EdgeId()] = undoAction.UpdatedEdge();
    m_edgesRTreeRequiresUpdate = true;
}

void Mesh::CommitAction(const DeleteEdgeAction& undoAction)
{
    m_edges[undoAction.EdgeId()] = {constants::missing::uintValue, constants::missing::uintValue};
    m_edgesRTreeRequiresUpdate = true;
}

void Mesh::CommitAction(const DeleteNodeAction& undoAction)
{
    m_nodes[undoAction.NodeId()] = {constants::missing::doubleValue, constants::missing::doubleValue};
    m_nodesRTreeRequiresUpdate = true;
}

void Mesh::CommitAction(NodeTranslationAction& undoAction)
{
    undoAction.Swap(m_nodes);
}

void Mesh::CommitAction(MeshConversionAction& undoAction)
{
    undoAction.Swap(m_nodes, m_projection);
}

void Mesh::CommitAction(FullUnstructuredGridUndo& undoAction)
{
    undoAction.Swap(m_nodes, m_edges);
    Administrate();
}

void Mesh::RestoreAction(const AddNodeAction& undoAction)
{
    m_nodes[undoAction.NodeId()] = Point(constants::missing::doubleValue, constants::missing::doubleValue);
    m_nodesNumEdges[undoAction.NodeId()] = 0;
    m_nodesRTreeRequiresUpdate = true;
}

void Mesh::RestoreAction(const AddEdgeAction& undoAction)
{
    m_edges[undoAction.EdgeId()] = {constants::missing::uintValue, constants::missing::uintValue};
    m_edgesRTreeRequiresUpdate = true;
}

void Mesh::RestoreAction(const ResetNodeAction& undoAction)
{
    m_nodes[undoAction.NodeId()] = undoAction.InitialNode();
    m_nodesRTreeRequiresUpdate = true;
    m_edgesRTreeRequiresUpdate = true;
}

void Mesh::RestoreAction(const ResetEdgeAction& undoAction)
{
    m_edges[undoAction.EdgeId()] = undoAction.InitialEdge();
    m_edgesRTreeRequiresUpdate = true;
}

void Mesh::RestoreAction(const DeleteEdgeAction& undoAction)
{
    m_edges[undoAction.EdgeId()] = undoAction.GetEdge();
    m_edgesRTreeRequiresUpdate = true;
}

void Mesh::RestoreAction(const DeleteNodeAction& undoAction)
{
    m_nodes[undoAction.NodeId()] = undoAction.Node();
    m_nodesRTreeRequiresUpdate = true;
}

void Mesh::RestoreAction(NodeTranslationAction& undoAction)
{
    undoAction.Swap(m_nodes);
}

void Mesh::RestoreAction(MeshConversionAction& undoAction)
{
    undoAction.Swap(m_nodes, m_projection);
}

void Mesh::RestoreAction(FullUnstructuredGridUndo& undoAction)
{
    undoAction.Swap(m_nodes, m_edges);
    Administrate();
}
