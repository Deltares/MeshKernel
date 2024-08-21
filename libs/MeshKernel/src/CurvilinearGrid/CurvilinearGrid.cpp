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

#include <utility>

#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridLine.hpp"
#include "MeshKernel/CurvilinearGrid/UndoActions/ResetCurvilinearNodeAction.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Operations.hpp"
#include "MeshKernel/Polygons.hpp"
#include "MeshKernel/Utilities/NumericFunctions.hpp"
#include "MeshKernel/Utilities/RTreeFactory.hpp"

using meshkernel::CurvilinearGrid;
using meshkernel::CurvilinearGridNodeIndices;

CurvilinearGrid::CurvilinearGrid(const CurvilinearGrid& grid) : m_projection(Projection::cartesian),
                                                                m_gridNodes(grid.m_gridNodes),
                                                                m_gridFacesMask(grid.m_gridFacesMask),
                                                                m_gridNodesTypes(grid.m_gridNodesTypes),
                                                                m_gridIndices(grid.m_gridIndices)
{
    m_RTrees.emplace(Location::Nodes, RTreeFactory::Create(m_projection));
    m_RTrees.emplace(Location::Edges, RTreeFactory::Create(m_projection));
    m_RTrees.emplace(Location::Faces, RTreeFactory::Create(m_projection));
}

CurvilinearGrid::CurvilinearGrid(CurvilinearGrid&& grid) noexcept : m_projection(Projection::cartesian),
                                                                    m_gridNodes(std::move(grid.m_gridNodes)),
                                                                    m_gridFacesMask(std::move(grid.m_gridFacesMask)),
                                                                    m_gridNodesTypes(std::move(grid.m_gridNodesTypes)),
                                                                    m_gridIndices(std::move(grid.m_gridIndices)),
                                                                    m_RTrees(std::move(grid.m_RTrees)) {}

CurvilinearGrid::CurvilinearGrid(Projection projection) : m_projection(projection)
{
    m_RTrees.emplace(Location::Nodes, RTreeFactory::Create(m_projection));
    m_RTrees.emplace(Location::Edges, RTreeFactory::Create(m_projection));
    m_RTrees.emplace(Location::Faces, RTreeFactory::Create(m_projection));
}

CurvilinearGrid::CurvilinearGrid(lin_alg::Matrix<Point> const& grid, Projection projection) : m_projection(projection)
{
    m_RTrees.emplace(Location::Nodes, RTreeFactory::Create(m_projection));
    m_RTrees.emplace(Location::Edges, RTreeFactory::Create(m_projection));
    m_RTrees.emplace(Location::Faces, RTreeFactory::Create(m_projection));
    SetGridNodes(grid);
}

CurvilinearGrid::CurvilinearGrid(lin_alg::Matrix<Point>&& grid, Projection projection) : m_projection(projection)
{
    m_RTrees.emplace(Location::Nodes, RTreeFactory::Create(m_projection));
    m_RTrees.emplace(Location::Edges, RTreeFactory::Create(m_projection));
    m_RTrees.emplace(Location::Faces, RTreeFactory::Create(m_projection));
    SetGridNodes(std::move(grid));
}

CurvilinearGrid& CurvilinearGrid::operator=(CurvilinearGrid&& copy) noexcept
{
    if (this != &copy)
    {
        m_gridNodes = std::move(copy.m_gridNodes);
        m_gridFacesMask = std::move(copy.m_gridFacesMask);
        m_gridNodesTypes = std::move(copy.m_gridNodesTypes);
        m_gridIndices = std::move(copy.m_gridIndices);
        m_RTrees = std::move(copy.m_RTrees);
        m_edges = std::move(copy.m_edges);

        m_projection = std::exchange(copy.m_projection, Projection::cartesian);
        m_nodesRTreeRequiresUpdate = std::exchange(copy.m_nodesRTreeRequiresUpdate, false);
        m_edgesRTreeRequiresUpdate = std::exchange(copy.m_edgesRTreeRequiresUpdate, false);
        m_facesRTreeRequiresUpdate = std::exchange(copy.m_facesRTreeRequiresUpdate, false);
        m_boundingBoxCache = std::exchange(copy.m_boundingBoxCache, BoundingBox());
        m_startOffset = std::exchange(copy.m_startOffset, CurvilinearGridNodeIndices(0, 0));
        m_endOffset = std::exchange(copy.m_endOffset, CurvilinearGridNodeIndices(0, 0));
    }

    return *this;
}

CurvilinearGrid& CurvilinearGrid::operator=(const CurvilinearGrid& copy)
{
    if (this != &copy)
    {
        m_projection = copy.m_projection;
        m_gridNodes = copy.m_gridNodes;
        m_gridFacesMask = copy.m_gridFacesMask;
        m_gridNodesTypes = copy.m_gridNodesTypes;
        m_gridIndices = copy.m_gridIndices;

        m_nodesRTreeRequiresUpdate = true;
        m_edgesRTreeRequiresUpdate = true;
        m_facesRTreeRequiresUpdate = true;

        m_RTrees.emplace(Location::Nodes, RTreeFactory::Create(m_projection));
        m_RTrees.emplace(Location::Edges, RTreeFactory::Create(m_projection));
        m_RTrees.emplace(Location::Faces, RTreeFactory::Create(m_projection));

        m_boundingBoxCache = copy.m_boundingBoxCache;

        m_edges = copy.m_edges;

        m_startOffset = copy.m_startOffset;
        m_endOffset = copy.m_endOffset;

        SetGridNodes(m_gridNodes);
    }

    return *this;
}

void CurvilinearGrid::SetGridNodes(const lin_alg::Matrix<Point>& gridNodes)
{
    if (gridNodes.rows() <= 1 || gridNodes.cols() <= 1)
    {
        throw std::invalid_argument("CurvilinearGrid::CurvilinearGrid: Invalid curvilinear grid nodes");
    }

    m_gridNodes = gridNodes;

    m_nodesRTreeRequiresUpdate = true;
    m_edgesRTreeRequiresUpdate = true;
    m_facesRTreeRequiresUpdate = true;

    m_gridIndices = ComputeNodeIndices();
}

void CurvilinearGrid::SetGridNodes(lin_alg::Matrix<Point>&& gridNodes)
{
    if (gridNodes.rows() <= 1 || gridNodes.cols() <= 1)
    {
        throw std::invalid_argument("CurvilinearGrid::CurvilinearGrid: Invalid curvilinear grid nodes");
    }

    m_gridNodes = std::move(gridNodes);

    m_nodesRTreeRequiresUpdate = true;
    m_edgesRTreeRequiresUpdate = true;
    m_facesRTreeRequiresUpdate = true;

    m_gridIndices = ComputeNodeIndices();
}

void CurvilinearGrid::Delete(std::shared_ptr<Polygons> polygons, UInt polygonIndex)
{
    // no polygons available
    if (polygons->IsEmpty())
    {
        return;
    }

    // no grid, return
    if (IsEmpty())
    {
        return;
    }

    const auto numN = NumN();
    const auto numM = NumM();

    lin_alg::Matrix<bool> nodeBasedMask(numN, numM);
    nodeBasedMask.fill(false);

    lin_alg::Matrix<bool> faceBasedMask(numN - 1, numM - 1);
    faceBasedMask.fill(true);

    // Mark points inside a polygonIndex
    for (UInt n = 0; n < numN; ++n)
    {
        for (UInt m = 0; m < numM; ++m)
        {
            const auto isInPolygon = polygons->IsPointInPolygon(GetNode(n, m), polygonIndex);

            if (isInPolygon)
            {
                nodeBasedMask(n, m) = true;
            }
        }
    }

    // Mark faces when all nodes are inside
    for (UInt n = 0; n < numN - 1; ++n)
    {
        for (UInt m = 0; m < numM - 1; ++m)
        {
            if (!nodeBasedMask(n, m) ||
                !nodeBasedMask(n + 1, m) ||
                !nodeBasedMask(n, m + 1) ||
                !nodeBasedMask(n + 1, m + 1))
            {
                faceBasedMask(n, m) = false;
            }
        }
    }

    // Mark only the nodes of faces completely included in the polygonIndex
    nodeBasedMask.fill(false);
    for (UInt n = 0; n < numN - 1; ++n)
    {
        for (UInt m = 0; m < numM - 1; ++m)
        {
            if (faceBasedMask(n, m))
            {
                nodeBasedMask(n, m) = true;
                nodeBasedMask(n + 1, m) = true;
                nodeBasedMask(n, m + 1) = true;
                nodeBasedMask(n + 1, m + 1) = true;
            }
        }
    }

    // mark points inside a polygonIndex
    for (UInt n = 0; n < numN; ++n)
    {
        for (UInt m = 0; m < numM; ++m)
        {
            if (!nodeBasedMask(n, m))
            {
                GetNode(n, m).x = constants::missing::doubleValue;
                GetNode(n, m).y = constants::missing::doubleValue;
            }
        }
    }
}

void CurvilinearGrid::BuildTree(Location location, const BoundingBox& boundingBox)
{
    switch (location)
    {
    case Location::Faces:
        if (m_facesRTreeRequiresUpdate || m_boundingBoxCache != boundingBox)
        {
            const auto faceCenters = ComputeFaceCenters();
            m_RTrees.at(Location::Faces)->BuildTree(faceCenters, boundingBox);
            m_facesRTreeRequiresUpdate = false;
            m_boundingBoxCache = boundingBox;
        }
        break;
    case Location::Nodes:
        if (m_nodesRTreeRequiresUpdate || m_boundingBoxCache != boundingBox)
        {
            const auto nodes = ComputeNodes();
            m_RTrees.at(Location::Nodes)->BuildTree(nodes, boundingBox);
            m_nodesRTreeRequiresUpdate = false;
            m_boundingBoxCache = boundingBox;
        }
        break;
    case Location::Edges:
        if (m_edgesRTreeRequiresUpdate || m_boundingBoxCache != boundingBox)
        {
            m_edges = ComputeEdges();
            const auto edgeCenters = ComputeEdgesCenters();
            m_RTrees.at(Location::Edges)->BuildTree(edgeCenters, boundingBox);
            m_edgesRTreeRequiresUpdate = false;
            m_boundingBoxCache = boundingBox;
        }
        break;
    case Location::Unknown:
    default:
        throw std::runtime_error("Invalid location");
    }
}

meshkernel::UInt CurvilinearGrid::FindLocationIndex(Point point,
                                                    Location location,
                                                    const BoundingBox& boundingBox)
{
    BuildTree(location, boundingBox);
    const auto& rtree = m_RTrees.at(location);
    if (rtree->Empty())
    {
        throw AlgorithmError("Empty RTree");
    }

    rtree->SearchNearestPoint(point);

    if (rtree->GetQueryResultSize() <= 0)
    {
        throw AlgorithmError("Query result size <= 0.");
    }

    return rtree->GetQueryResult(0);
}

std::tuple<CurvilinearGridNodeIndices, CurvilinearGridNodeIndices> CurvilinearGrid::GetEdgeNodeIndices(Point const& point)
{
    BuildTree(Location::Edges);

    m_RTrees.at(Location::Edges)->SearchNearestPoint(point);

    if (m_RTrees.at(Location::Edges)->GetQueryResultSize() == 0)
    {
        return {{}, {}};
    }

    const auto edgeIndex = m_RTrees.at(Location::Edges)->GetQueryResult(0);
    auto const firstNode = m_edges[edgeIndex].first;
    auto const secondNode = m_edges[edgeIndex].second;

    return {m_gridIndices[firstNode], m_gridIndices[secondNode]};
}

bool CurvilinearGrid::AreFaceNodesValid(UInt n, UInt m) const
{
    return GetNode(n, m).IsValid() &&
           GetNode(n, m + 1).IsValid() &&
           GetNode(n + 1, m).IsValid() &&
           GetNode(n + 1, m + 1).IsValid();
}

std::tuple<CurvilinearGridNodeIndices, CurvilinearGridNodeIndices> CurvilinearGrid::ComputeBlockFromCornerPoints(Point const& firstCornerPoint, Point const& secondCornerPoint)
{
    // Get the m and n indices from the point coordinates
    auto const firstNodeIndex = FindLocationIndex(firstCornerPoint, Location::Nodes);
    auto const secondNodeIndex = FindLocationIndex(secondCornerPoint, Location::Nodes);

    // Compute bounding box as node indices from corner points
    return ComputeBlockFromCornerPoints(m_gridIndices[firstNodeIndex], m_gridIndices[secondNodeIndex]);
}

std::tuple<CurvilinearGridNodeIndices, CurvilinearGridNodeIndices>
CurvilinearGrid::ComputeBlockFromCornerPoints(const CurvilinearGridNodeIndices& firstNode, const CurvilinearGridNodeIndices& secondNode) const
{
    CurvilinearGridNodeIndices lowerLeft(std::min(firstNode.m_n, secondNode.m_n), std::min(firstNode.m_m, secondNode.m_m));
    CurvilinearGridNodeIndices upperRight(std::max(firstNode.m_n, secondNode.m_n), std::max(firstNode.m_m, secondNode.m_m));

    if (!lowerLeft.IsValid() || !upperRight.IsValid())
    {
        throw ConstraintError("Invalid index: first index - {{{}, {}}}, second index - {{{}, {}}}", lowerLeft.m_m, lowerLeft.m_n, upperRight.m_m, upperRight.m_n);
    }

    if (lowerLeft.m_m >= NumM() || lowerLeft.m_n >= NumN())
    {
        throw ConstraintError("Invalid index: first index {{{}, {}}} not in mesh limits {{{}, {}}}", lowerLeft.m_m, lowerLeft.m_n, NumM(), NumN());
    }

    if (upperRight.m_m >= NumM() || upperRight.m_n >= NumN())
    {
        throw ConstraintError("Invalid index: second index {{{}, {}}} not in mesh limits {{{}, {}}}", upperRight.m_m, upperRight.m_n, NumM(), NumN());
    }

    return {lowerLeft, upperRight};
}

void CurvilinearGrid::ComputeGridFacesMask()
{
    // Flag valid faces
    lin_alg::ResizeAndFillMatrix(m_gridFacesMask, FullNumN() - 1, FullNumM() - 1, false, false);
    std::vector<std::vector<bool>> isFaceMaskValidFlatArray(NumN() - 1, std::vector<bool>(NumM() - 1, false));

    for (UInt n = 0; n < NumN() - 1; ++n)
    {
        for (UInt m = 0; m < NumM() - 1; ++m)
        {
            // Only if all grid nodes of the face are valid, the face is valid
            if (!AreFaceNodesValid(n, m))
            {
                continue;
            }
            isFaceMaskValidFlatArray[n][m] = true;
            IsFaceMaskValid(n, m) = true;
        }
    }
}

void CurvilinearGrid::RemoveInvalidNodes(bool invalidNodesToRemove)
{

    if (!invalidNodesToRemove)
    {
        return;
    }

    // Compute the face mask
    ComputeGridFacesMask();

    std::vector validNodeMask(NumN(), std::vector(NumM(), false));
    for (UInt n = 0; n < NumN() - 1; ++n)
    {
        for (UInt m = 0; m < NumM() - 1; ++m)
        {
            if (IsFaceMaskValid(n, m))
            {
                validNodeMask[n][m] = true;
                validNodeMask[n][m + 1] = true;
                validNodeMask[n + 1][m] = true;
                validNodeMask[n + 1][m + 1] = true;
            }
        }
    }

    invalidNodesToRemove = false;
    for (UInt n = 0; n < NumN(); ++n)
    {
        for (UInt m = 0; m < NumM(); ++m)
        {
            if (!validNodeMask[n][m] && GetNode(n, m).IsValid())
            {
                GetNode(n, m) = {constants::missing::doubleValue, constants::missing::doubleValue};
                invalidNodesToRemove = true;
            }
        }
    }

    RemoveInvalidNodes(invalidNodesToRemove);
}

void CurvilinearGrid::ComputeGridNodeTypes()
{
    RemoveInvalidNodes(true);
    lin_alg::ResizeAndFillMatrix(m_gridNodesTypes, FullNumN(), FullNumM(), false, NodeType::Invalid);

    // Flag faces based on boundaries
    for (UInt n = 0; n < NumN(); ++n)
    {
        for (UInt m = 0; m < NumM(); ++m)
        {

            if (!GetNode(n, m).IsValid())
            {
                continue;
            }

            // Bottom side
            if (m == 0 && n == 0)
            {
                GetNodeType(n, m) = NodeType::BottomLeft;
                continue;
            }
            if (m == 0 && n == NumN() - 1)
            {
                GetNodeType(n, m) = NodeType::BottomRight;
                continue;
            }
            if (m == 0 && !GetNode(n - 1, m).IsValid())
            {
                GetNodeType(n, m) = NodeType::BottomLeft;
                continue;
            }
            if (m == 0 && !GetNode(n + 1, m).IsValid())
            {
                GetNodeType(n, m) = NodeType::BottomRight;
                continue;
            }
            if (m == 0)
            {
                GetNodeType(n, m) = NodeType::Bottom;
                continue;
            }
            // Upper side
            if (m == NumM() - 1 && n == 0)
            {
                GetNodeType(n, m) = NodeType::UpperLeft;
                continue;
            }
            if (m == NumM() - 1 && n == NumN() - 1)
            {
                GetNodeType(n, m) = NodeType::UpperRight;
                continue;
            }
            if (m == NumM() - 1 && !GetNode(n - 1, m).IsValid())
            {
                GetNodeType(n, m) = NodeType::UpperLeft;
                continue;
            }
            if (m == NumM() - 1 && !GetNode(n + 1, m).IsValid())
            {
                GetNodeType(n, m) = NodeType::UpperRight;
                continue;
            }
            if (m == NumM() - 1)
            {
                GetNodeType(n, m) = NodeType::Up;
                continue;
            }
            // Bottom side
            if (n == 0 && !GetNode(n, m - 1).IsValid())
            {
                GetNodeType(n, m) = NodeType::BottomLeft;
                continue;
            }
            if (n == 0 && !GetNode(n, m + 1).IsValid())
            {
                GetNodeType(n, m) = NodeType::UpperRight;
                continue;
            }
            if (n == 0)
            {
                GetNodeType(n, m) = NodeType::Left;
                continue;
            }
            // Upper side
            if (n == NumN() - 1 && !GetNode(n, m - 1).IsValid())
            {
                GetNodeType(n, m) = NodeType::BottomRight;
                continue;
            }
            if (n == NumN() - 1 && !GetNode(n, m + 1).IsValid())
            {
                GetNodeType(n, m) = NodeType::UpperRight;
                continue;
            }
            if (n == NumN() - 1)
            {
                GetNodeType(n, m) = NodeType::Right;
                continue;
            }

            auto const isBottomRightFaceValid = IsFaceMaskValid(n, m - 1);
            auto const isBottomLeftFaceValid = IsFaceMaskValid(n - 1, m - 1);
            auto const isTopRightFaceValid = IsFaceMaskValid(n, m);
            auto const isTopLeftFaceValid = IsFaceMaskValid(n - 1, m);

            if (isTopRightFaceValid &&
                isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                GetNodeType(n, m) = NodeType::InternalValid;
                continue;
            }
            if (!isTopRightFaceValid &&
                isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                GetNodeType(n, m) = NodeType::BottomLeft;
                continue;
            }
            if (isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                GetNodeType(n, m) = NodeType::BottomRight;
                continue;
            }
            if (isTopRightFaceValid &&
                isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                GetNodeType(n, m) = NodeType::UpperRight;
                continue;
            }
            if (isTopRightFaceValid &&
                isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                GetNodeType(n, m) = NodeType::UpperLeft;
                continue;
            }

            if (isTopRightFaceValid &&
                isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                GetNodeType(n, m) = NodeType::Bottom;
                continue;
            }
            if (isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                GetNodeType(n, m) = NodeType::Left;
                continue;
            }

            if (!isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                GetNodeType(n, m) = NodeType::Up;
                continue;
            }

            if (!isTopRightFaceValid &&
                isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                GetNodeType(n, m) = NodeType::Right;
                continue;
            }

            if (isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                GetNodeType(n, m) = NodeType::BottomLeft;
                continue;
            }

            if (!isTopRightFaceValid &&
                isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                GetNodeType(n, m) = NodeType::BottomRight;
                continue;
            }

            if (!isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                GetNodeType(n, m) = NodeType::UpperRight;
                continue;
            }

            if (!isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                GetNodeType(n, m) = NodeType::UpperLeft;
            }
        }
    }
}

meshkernel::UndoActionPtr CurvilinearGrid::InsertFace(Point const& point)
{
    if (!point.IsValid())
    {
        throw std::invalid_argument("CurvilinearGrid::InsertFace: invalid point provided");
    }

    // Gets the indices of the closest edge to the specified point (they are neighbors by construction)
    auto const [firstNode, secondNode] = GetEdgeNodeIndices(point);

    if (!firstNode.IsValid() || !secondNode.IsValid())
    {
        throw std::invalid_argument("CurvilinearGrid::InsertFace: no valid nodes found");
    }

    // Compute the grid node types
    ComputeGridNodeTypes();

    // Add a new edge
    UndoActionPtr undoAction = AddEdge(firstNode, secondNode);

    // Re-compute quantities
    ComputeGridNodeTypes();
    m_gridIndices = ComputeNodeIndices();

    return undoAction;
}

std::tuple<bool, meshkernel::UndoActionPtr> CurvilinearGrid::AddGridLineAtBoundary(CurvilinearGridNodeIndices const& firstNode,
                                                                                   CurvilinearGridNodeIndices const& secondNode)
{

    if (!firstNode.IsValid() || !secondNode.IsValid())
    {
        throw ConstraintError("Invalid indices");
    }

    if (!IsInRange(firstNode.m_m, 0u, NumM()) || !IsInRange(firstNode.m_n, 0u, NumN()))
    {
        throw ConstraintError("First index {{{}, {}}} not in mesh limits {{{}, {}}},  {{{}, {}}}",
                              firstNode.m_n, firstNode.m_m, NumN(), NumM(), m_gridNodes.rows(), m_gridNodes.cols());
    }

    if (!IsInRange(secondNode.m_m, 0u, NumM()) || !IsInRange(secondNode.m_n, 0u, NumN()))
    {
        throw ConstraintError("Second index {{{}, {}}} not in mesh limits {{{}, {}}}", secondNode.m_n, secondNode.m_m, NumN(), NumM());
    }

    // If both nodes are invalid, we can substitute the invalid values. New allocation is not needed.
    bool const areNodesValid = GetNode(firstNode.m_n, firstNode.m_m).IsValid() && GetNode(secondNode.m_n, secondNode.m_m).IsValid();

    // Allocation depends on directions
    bool gridSizeChanged = false;
    auto const gridLineType = GetBoundaryGridLineType(firstNode, secondNode);

    UndoActionPtr undoAction;

    if (areNodesValid)
    {

        if (gridLineType == BoundaryGridLineType::Bottom)
        {
            if (firstNode.m_n == 0 || secondNode.m_n == 0)
            {

                // n-direction
                if (m_startOffset.m_n == 0)
                {
                    lin_alg::InsertRow(m_gridNodes, lin_alg::RowVector<Point>(FullNumM()), 0);
                }
                else
                {
                    m_startOffset.m_n -= 1;
                }

                undoAction = AddGridLineUndoAction::Create(*this, {1, 0}, {0, 0});
                gridSizeChanged = true;
            }
        }

        if (gridLineType == BoundaryGridLineType::Top)
        {
            if (firstNode.m_n == NumN() - 1 || secondNode.m_n == NumN() - 1)
            {
                // n-direction
                if (m_endOffset.m_n == 0)
                {
                    lin_alg::InsertRow(m_gridNodes, lin_alg::RowVector<Point>(FullNumM()), FullNumN());
                }
                else
                {
                    m_endOffset.m_n -= 1;
                }

                undoAction = AddGridLineUndoAction::Create(*this, {0, 0}, {1, 0});
                gridSizeChanged = true;
            }
        }

        if (gridLineType == BoundaryGridLineType::Right)
        {
            if (firstNode.m_m == NumM() - 1 || secondNode.m_m == NumM() - 1)
            {
                // m-direction
                if (m_endOffset.m_m == 0)
                {
                    lin_alg::InsertCol(m_gridNodes, lin_alg::ColVector<Point>(FullNumN()), FullNumM());
                }
                else
                {
                    m_endOffset.m_m -= 1;
                }

                undoAction = AddGridLineUndoAction::Create(*this, {0, 0}, {0, 1});
                gridSizeChanged = true;
            }
        }

        if (gridLineType == BoundaryGridLineType::Left)
        {

            if (firstNode.m_m == 0 || secondNode.m_m == 0)
            {

                // m-direction
                if (m_startOffset.m_m == 0)
                {
                    lin_alg::InsertCol(m_gridNodes, lin_alg::ColVector<Point>(FullNumN()), 0);
                }
                else
                {
                    m_startOffset.m_m -= 1;
                }

                undoAction = AddGridLineUndoAction::Create(*this, {0, 1}, {0, 0});
                gridSizeChanged = true;
            }
        }
    }

    return {gridSizeChanged, std::move(undoAction)};
}

CurvilinearGrid::BoundaryGridLineType CurvilinearGrid::GetBoundaryGridLineType(CurvilinearGridNodeIndices const& firstNode, CurvilinearGridNodeIndices const& secondNode) const
{
    auto const firstNodeType = GetNodeType(firstNode.m_n, firstNode.m_m);
    auto const secondNodeType = GetNodeType(secondNode.m_n, secondNode.m_m);

    if (firstNodeType == NodeType::InternalValid || firstNodeType == NodeType::Invalid ||
        secondNodeType == NodeType::InternalValid || secondNodeType == NodeType::Invalid)
    {
        throw std::invalid_argument("CurvilinearGrid::GetBoundaryGridLineType: Not a boundary grid line");
    }

    if (firstNodeType == NodeType::Bottom && secondNodeType == NodeType::Bottom)
    {
        return BoundaryGridLineType::Left;
    }
    if (firstNodeType == NodeType::Up && secondNodeType == NodeType::Up)
    {
        return BoundaryGridLineType::Right;
    }
    if (firstNodeType == NodeType::Left && secondNodeType == NodeType::Left)
    {
        return BoundaryGridLineType::Bottom;
    }
    if (firstNodeType == NodeType::Right && secondNodeType == NodeType::Right)
    {
        return BoundaryGridLineType::Top;
    }

    // Corner nodes
    if (firstNode.m_m == secondNode.m_m) // Left or Right
    {
        if (firstNode.m_m + 1 < m_gridNodes.cols() && GetNode(firstNode.m_n, firstNode.m_m + 1).IsValid() && GetNode(secondNode.m_n, secondNode.m_m + 1).IsValid())
        {
            return BoundaryGridLineType::Left;
        }
        return BoundaryGridLineType::Right;
    }
    if (firstNode.m_n == secondNode.m_n) // Bottom or top
    {
        if (firstNode.m_n + 1 < m_gridNodes.rows() && GetNode(firstNode.m_n + 1, firstNode.m_m).IsValid() && GetNode(secondNode.m_n + 1, secondNode.m_m).IsValid())
        {
            return BoundaryGridLineType::Bottom;
        }
        return BoundaryGridLineType::Top;
    }

    throw std::invalid_argument("CurvilinearGrid::GetBoundaryGridLineType: Invalid node types");
}

meshkernel::UndoActionPtr CurvilinearGrid::AddEdge(CurvilinearGridNodeIndices const& firstNode, CurvilinearGridNodeIndices const& secondNode)
{

    // Allocate new grid line if needed
    auto const gridLineType = GetBoundaryGridLineType(firstNode, secondNode);
    // Allocation depends on directions

    std::unique_ptr<CompoundUndoAction> undoAddEdge = CompoundUndoAction::Create();

    m_edgesRTreeRequiresUpdate = true;

    if (gridLineType == BoundaryGridLineType::Bottom)
    {
        auto const firstNewNodeCoordinates = GetNode(firstNode.m_n, firstNode.m_m) * 2.0 - GetNode(firstNode.m_n + 1, firstNode.m_m);
        auto const secondNewNodeCoordinates = GetNode(secondNode.m_n, secondNode.m_m) * 2.0 - GetNode(secondNode.m_n + 1, secondNode.m_m);
        auto [isGridLineAdded, addGridLineUndo] = AddGridLineAtBoundary(firstNode, secondNode);

        if (isGridLineAdded)
        {
            undoAddEdge->Add(std::move(addGridLineUndo));
            undoAddEdge->Add(MoveNode(CurvilinearGridNodeIndices(0, firstNode.m_m), firstNewNodeCoordinates));
            undoAddEdge->Add(MoveNode(CurvilinearGridNodeIndices(0, secondNode.m_m), secondNewNodeCoordinates));
        }
        else
        {
            undoAddEdge->Add(MoveNode(CurvilinearGridNodeIndices(firstNode.m_n - 1, firstNode.m_m), firstNewNodeCoordinates));
            undoAddEdge->Add(MoveNode(CurvilinearGridNodeIndices(secondNode.m_n - 1, secondNode.m_m), secondNewNodeCoordinates));
        }
    }

    if (gridLineType == BoundaryGridLineType::Top)
    {
        auto const firstNewNodeCoordinates = GetNode(firstNode.m_n, firstNode.m_m) * 2.0 - GetNode(firstNode.m_n - 1, firstNode.m_m);
        auto const secondNewNodeCoordinates = GetNode(secondNode.m_n, secondNode.m_m) * 2.0 - GetNode(secondNode.m_n - 1, secondNode.m_m);
        auto [isGridLineAdded, addGridLineUndo] = AddGridLineAtBoundary(firstNode, secondNode);

        if (isGridLineAdded)
        {
            undoAddEdge->Add(std::move(addGridLineUndo));
        }

        undoAddEdge->Add(MoveNode(CurvilinearGridNodeIndices(firstNode.m_n + 1, firstNode.m_m), firstNewNodeCoordinates));
        undoAddEdge->Add(MoveNode(CurvilinearGridNodeIndices(secondNode.m_n + 1, secondNode.m_m), secondNewNodeCoordinates));
    }

    if (gridLineType == BoundaryGridLineType::Left)
    {
        auto const firstNewNodeCoordinates = GetNode(firstNode.m_n, firstNode.m_m) * 2.0 - GetNode(firstNode.m_n, firstNode.m_m + 1);
        auto const secondNewNodeCoordinates = GetNode(secondNode.m_n, secondNode.m_m) * 2.0 - GetNode(secondNode.m_n, secondNode.m_m + 1);
        auto [isGridLineAdded, addGridLineUndo] = AddGridLineAtBoundary(firstNode, secondNode);

        if (isGridLineAdded)
        {
            undoAddEdge->Add(std::move(addGridLineUndo));
            // Assign the new coordinates
            undoAddEdge->Add(MoveNode(CurvilinearGridNodeIndices(firstNode.m_n, 0), firstNewNodeCoordinates));
            undoAddEdge->Add(MoveNode(CurvilinearGridNodeIndices(secondNode.m_n, 0), secondNewNodeCoordinates));
        }
        else
        {
            undoAddEdge->Add(MoveNode(CurvilinearGridNodeIndices(firstNode.m_n, firstNode.m_m - 1), firstNewNodeCoordinates));
            undoAddEdge->Add(MoveNode(CurvilinearGridNodeIndices(secondNode.m_n, secondNode.m_m - 1), secondNewNodeCoordinates));
        }
    }

    if (gridLineType == BoundaryGridLineType::Right)
    {
        auto const firstNewNodeCoordinates = GetNode(firstNode.m_n, firstNode.m_m) * 2.0 - GetNode(firstNode.m_n, firstNode.m_m - 1);
        auto const secondNewNodeCoordinates = GetNode(secondNode.m_n, secondNode.m_m) * 2.0 - GetNode(secondNode.m_n, secondNode.m_m - 1);
        auto [isGridLineAdded, addGridLineUndo] = AddGridLineAtBoundary(firstNode, secondNode);

        if (isGridLineAdded)
        {
            undoAddEdge->Add(std::move(addGridLineUndo));
        }

        undoAddEdge->Add(MoveNode(CurvilinearGridNodeIndices(firstNode.m_n, firstNode.m_m + 1), firstNewNodeCoordinates));
        undoAddEdge->Add(MoveNode(CurvilinearGridNodeIndices(secondNode.m_n, secondNode.m_m + 1), secondNewNodeCoordinates));
    }

    if (gridLineType == BoundaryGridLineType::Bottom ||
        gridLineType == BoundaryGridLineType::Right ||
        gridLineType == BoundaryGridLineType::Top ||
        gridLineType == BoundaryGridLineType::Left)
    {
        m_nodesRTreeRequiresUpdate = true;
        m_edgesRTreeRequiresUpdate = true;
        m_facesRTreeRequiresUpdate = true;
    }

    return undoAddEdge;
}

std::tuple<double, double, double>
CurvilinearGrid::ComputeDirectionalSmoothingFactors(CurvilinearGridNodeIndices const& gridpoint,
                                                    const CurvilinearGridNodeIndices& pointOnSmoothingLineIndices,
                                                    const CurvilinearGridNodeIndices& lowerLeftIndices,
                                                    const CurvilinearGridNodeIndices& upperRightIndices)
{
    // horizontal smoothing factor
    // integer
    const auto horizontalDelta = gridpoint.m_m > pointOnSmoothingLineIndices.m_m ? gridpoint.m_m - pointOnSmoothingLineIndices.m_m : pointOnSmoothingLineIndices.m_m - gridpoint.m_m;
    // integer
    const auto maxHorizontalDelta = gridpoint.m_m > pointOnSmoothingLineIndices.m_m ? upperRightIndices.m_m - pointOnSmoothingLineIndices.m_m : pointOnSmoothingLineIndices.m_m - lowerLeftIndices.m_m;
    // double
    const auto mSmoothingFactor = maxHorizontalDelta == 0 ? 1.0 : (1.0 + std::cos(M_PI * static_cast<double>(horizontalDelta) / static_cast<double>(maxHorizontalDelta))) * 0.5;

    // vertical smoothing factor
    const auto verticalDelta = gridpoint.m_n > pointOnSmoothingLineIndices.m_n ? gridpoint.m_n - pointOnSmoothingLineIndices.m_n : pointOnSmoothingLineIndices.m_n - gridpoint.m_n;
    const auto maxVerticalDelta = gridpoint.m_n > pointOnSmoothingLineIndices.m_n ? upperRightIndices.m_n - pointOnSmoothingLineIndices.m_n : pointOnSmoothingLineIndices.m_n - lowerLeftIndices.m_n;
    const auto nSmoothingFactor = maxVerticalDelta == 0 ? 1.0 : (1.0 + std::cos(M_PI * static_cast<double>(verticalDelta) / static_cast<double>(maxVerticalDelta))) * 0.5;

    // mixed smoothing factor
    const auto mixedSmoothingFactor = std::sqrt(nSmoothingFactor * mSmoothingFactor);

    return {mSmoothingFactor, nSmoothingFactor, mixedSmoothingFactor};
}

double CurvilinearGrid::ComputeAverageNodalDistance(CurvilinearGridNodeIndices const& index, CurvilinearGridLine::GridLineDirection direction)
{
    if (index.m_m > NumM() || index.m_n > NumN())
    {
        throw std::invalid_argument("CurvilinearGrid::ComputeAverageNodalDistance: invalid index coordinates");
    }
    if (direction != CurvilinearGridLine::GridLineDirection::MDirection && direction != CurvilinearGridLine::GridLineDirection::NDirection)
    {
        throw std::invalid_argument("CurvilinearGrid::ComputeAverageNodalDistance: invalid direction");
    }

    if (direction == CurvilinearGridLine::GridLineDirection::MDirection)
    {
        double numEdges = 0.0;
        double leftDistance = 0.0;
        double rightDistance = 0.0;
        if (index.m_m > 0 && GetNode(index.m_n, index.m_m - 1).IsValid())
        {
            leftDistance = ComputeDistance(GetNode(index.m_n, index.m_m), GetNode(index.m_n, index.m_m - 1), m_projection);
            numEdges += 1;
        }
        if (index.m_m + 1 < NumM() && GetNode(index.m_n, index.m_m + 1).IsValid())
        {
            rightDistance = ComputeDistance(GetNode(index.m_n, index.m_m), GetNode(index.m_n, index.m_m + 1), m_projection);
            numEdges += 1;
        }
        return numEdges == 0 ? 0.0 : (leftDistance + rightDistance) / numEdges;
    }
    if (direction == CurvilinearGridLine::GridLineDirection::NDirection)
    {
        double numEdges = 0.0;
        double bottomDistance = 0.0;
        double upDistance = 0.0;
        if (index.m_n > 0 && GetNode(index.m_n - 1, index.m_m).IsValid())
        {
            bottomDistance = ComputeDistance(GetNode(index.m_n, index.m_m), GetNode(index.m_n - 1, index.m_m), m_projection);
            numEdges += 1;
        }
        if (index.m_n + 1 < NumN() && GetNode(index.m_n + 1, index.m_m).IsValid())
        {
            upDistance = ComputeDistance(GetNode(index.m_n, index.m_m), GetNode(index.m_n + 1, index.m_m), m_projection);
            numEdges += 1;
        }
        return numEdges == 0 ? 0.0 : (bottomDistance + upDistance) / numEdges;
    }

    throw std::invalid_argument("CurvilinearGrid::ComputeAverageNodalDistance: Invalid direction");
}

meshkernel::Point CurvilinearGrid::TransformDisplacement(Point const& displacement, CurvilinearGridNodeIndices const& node, bool isLocal) const
{
    Point left = GetNode(node.m_n, node.m_m);
    Point right = left;
    if (node.m_n < NumN() - 1 && GetNode(node.m_n + 1, node.m_m).IsValid())
    {
        right = GetNode(node.m_n + 1, node.m_m);
    }
    if (node.m_n > 0 && GetNode(node.m_n - 1, node.m_m).IsValid())
    {
        left = GetNode(node.m_n - 1, node.m_m);
    }

    const auto horizontalDistance = ComputeDistance(right, left, m_projection);
    const auto horizontalDelta = right - left;

    if (isLocal && horizontalDistance > 0.0)
    {
        return {(displacement.x * horizontalDelta.x + displacement.y * horizontalDelta.y) / horizontalDistance,
                (displacement.y * horizontalDelta.x - displacement.x * horizontalDelta.y) / horizontalDistance};
    }
    if (!isLocal && horizontalDistance > 0.0)
    {
        return {(displacement.x * horizontalDelta.x - displacement.y * horizontalDelta.y) / horizontalDistance,
                (displacement.x * horizontalDelta.y + displacement.y * horizontalDelta.x) / horizontalDistance};
    }

    return {0.0, 0.0};
}

meshkernel::UndoActionPtr CurvilinearGrid::DeleteNode(Point const& point)
{
    // Get the m and n indices from the point coordinates
    auto const nodeToDeleteIndex = FindLocationIndex(point, Location::Nodes);
    const auto nodeToDelete = m_gridIndices[nodeToDeleteIndex];

    std::unique_ptr<ResetCurvilinearNodeAction> undoAction;

    if (nodeToDelete.IsValid())
    {
        undoAction = ResetCurvilinearNodeAction::Create(*this, nodeToDelete, GetNode(nodeToDelete),
                                                        {constants::missing::doubleValue, constants::missing::doubleValue},
                                                        true /* recomputeNodeTypes */);

        // Invalidate gridnodes
        GetNode(nodeToDelete.m_n, nodeToDelete.m_m) = {constants::missing::doubleValue, constants::missing::doubleValue};
        // Re-compute quantities
        ComputeGridNodeTypes();
        m_gridIndices = ComputeNodeIndices();
    }

    return undoAction;
}

meshkernel::UndoActionPtr CurvilinearGrid::MoveNode(const CurvilinearGridNodeIndices& nodeIndex, Point const& toPoint)
{
    // Check the node indices are valid
    if (!nodeIndex.IsValid())
    {
        throw std::invalid_argument("CurvilinearGrid::MoveNode node indices not found");
    }

    std::unique_ptr<ResetCurvilinearNodeAction> undoAction = ResetCurvilinearNodeAction::Create(*this, nodeIndex, GetNode(nodeIndex), toPoint,
                                                                                                false /* recomputeNodeTypes */);

    // move fromPoint to toPoint
    GetNode(nodeIndex.m_n, nodeIndex.m_m) = toPoint;
    return undoAction;
}

meshkernel::UndoActionPtr CurvilinearGrid::MoveNode(Point const& fromPoint, Point const& toPoint)
{
    // Get the node indices of fromPoint
    auto const nodeToMoveIndex = FindLocationIndex(fromPoint, Location::Nodes);
    const auto nodeToMove = m_gridIndices[nodeToMoveIndex];
    return MoveNode(nodeToMove, toPoint);
}

meshkernel::BoundingBox CurvilinearGrid::GetBoundingBox() const
{

    Point lowerLeft(std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
    Point upperRight(-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());
    UInt last = static_cast<UInt>(NumM() - 1);

    // Only need to loop over boundary nodes

    // First loop over lower boundary (i,0)
    for (UInt i = 0; i < NumM(); ++i)
    {
        lowerLeft.x = std::min(lowerLeft.x, GetNode(0, i).x);
        lowerLeft.y = std::min(lowerLeft.y, GetNode(0, i).y);
        upperRight.x = std::max(upperRight.x, GetNode(0, i).x);
        upperRight.y = std::max(upperRight.y, GetNode(0, i).y);
    }

    // First loop over right boundary (last,i)
    for (UInt i = 0; i < NumN(); ++i)
    {
        lowerLeft.x = std::min(lowerLeft.x, GetNode(i, last).x);
        lowerLeft.y = std::min(lowerLeft.y, GetNode(i, last).y);
        upperRight.x = std::max(upperRight.x, GetNode(i, last).x);
        upperRight.y = std::max(upperRight.y, GetNode(i, last).y);
    }

    // This assumes that each column has the same number of points
    last = NumN() - 1;

    // First loop over upper boundary (i,last)
    for (UInt i = 0; i < NumM(); ++i)
    {
        lowerLeft.x = std::min(lowerLeft.x, GetNode(last, i).x);
        lowerLeft.y = std::min(lowerLeft.y, GetNode(last, i).y);
        upperRight.x = std::max(upperRight.x, GetNode(last, i).x);
        upperRight.y = std::max(upperRight.y, GetNode(last, i).y);
    }

    // First loop over left boundary (0,i)
    for (UInt i = 0; i < NumN(); ++i)
    {
        lowerLeft.x = std::min(lowerLeft.x, GetNode(i, 0).x);
        lowerLeft.y = std::min(lowerLeft.y, GetNode(i, 0).y);
        upperRight.x = std::max(upperRight.x, GetNode(i, 0).x);
        upperRight.y = std::max(upperRight.y, GetNode(i, 0).y);
    }

    return {lowerLeft, upperRight};
}

void CurvilinearGrid::RestoreAction(const AddGridLineUndoAction& undoAction)
{
    m_startOffset += undoAction.StartOffset();
    m_endOffset += undoAction.EndOffset();

    // The node and edge trees need to be rebuilt
    m_nodesRTreeRequiresUpdate = true;
    m_edgesRTreeRequiresUpdate = true;
    m_gridIndices = ComputeNodeIndices();

    // Since the size of the grid has changed the node-vector needs to be reset
    ComputeGridNodeTypes();
}

void CurvilinearGrid::CommitAction(const AddGridLineUndoAction& undoAction)
{
    m_startOffset -= undoAction.StartOffset();
    m_endOffset -= undoAction.EndOffset();
    m_nodesRTreeRequiresUpdate = true;
    m_edgesRTreeRequiresUpdate = true;
    m_gridIndices = ComputeNodeIndices();
    ComputeGridNodeTypes();
}

void CurvilinearGrid::RestoreAction(CurvilinearGridBlockUndoAction& undoAction)
{
    undoAction.Swap(*this);
    m_nodesRTreeRequiresUpdate = true;
    m_edgesRTreeRequiresUpdate = true;
    ComputeGridNodeTypes();
}

void CurvilinearGrid::CommitAction(CurvilinearGridBlockUndoAction& undoAction)
{
    undoAction.Swap(*this);
    m_nodesRTreeRequiresUpdate = true;
    m_edgesRTreeRequiresUpdate = true;
    ComputeGridNodeTypes();
}

void CurvilinearGrid::RestoreAction(CurvilinearGridRefinementUndoAction& undoAction)
{
    undoAction.Swap(m_gridNodes, m_startOffset, m_endOffset);
    m_nodesRTreeRequiresUpdate = true;
    m_edgesRTreeRequiresUpdate = true;
    m_gridIndices = ComputeNodeIndices();
    ComputeGridNodeTypes();
}

void CurvilinearGrid::CommitAction(CurvilinearGridRefinementUndoAction& undoAction)
{
    undoAction.Swap(m_gridNodes, m_startOffset, m_endOffset);
    m_nodesRTreeRequiresUpdate = true;
    m_edgesRTreeRequiresUpdate = true;
    m_gridIndices = ComputeNodeIndices();
    ComputeGridNodeTypes();
}

void CurvilinearGrid::RestoreAction(const ResetCurvilinearNodeAction& undoAction)
{
    GetNode(undoAction.NodeId()) = undoAction.InitialNode();

    if (undoAction.RecalculateNodeTypes())
    {
        ComputeGridNodeTypes();
    }
}

void CurvilinearGrid::CommitAction(const ResetCurvilinearNodeAction& undoAction)
{
    GetNode(undoAction.NodeId()) = undoAction.UpdatedNode();

    if (undoAction.RecalculateNodeTypes())
    {
        ComputeGridNodeTypes();
    }
}

std::vector<meshkernel::Point> CurvilinearGrid::ComputeNodes() const
{
    if (!IsValid())
    {
        throw AlgorithmError("Invalid curvilinear grid ");
    }

    std::vector<Point> result(NumN() * NumM());
    UInt ind = 0;
    for (UInt n = 0; n < NumN(); n++)
    {
        for (UInt m = 0; m < NumM(); m++)
        {

            result[ind] = GetNode(n, m);
            ind++;
        }
    }
    return result;
}

std::vector<meshkernel::Edge> CurvilinearGrid::ComputeEdges() const
{
    if (!IsValid())
    {
        throw AlgorithmError("Invalid curvilinear grid ");
    }

    const auto numM = NumM();
    const auto numN = NumN();

    std::vector<Edge> result(numM * (numN - 1) +
                             (numM - 1) * numN);

    UInt ind = 0;

    for (UInt n = 0; n < numN - 1; n++)
    {
        for (UInt m = 0; m < numM; m++)
        {
            result[ind].first = numM * n + m;
            result[ind].second = numM * (n + 1) + m;
            ind++;
        }
    }
    for (UInt n = 0; n < numN; n++)
    {
        for (UInt m = 0; m < numM - 1; m++)
        {
            result[ind].first = numM * n + m;
            result[ind].second = numM * n + m + 1;
            ind++;
        }
    }

    return result;
}

std::vector<meshkernel::Point> CurvilinearGrid::ComputeEdgesCenters() const
{
    std::vector<Point> result(GetNumEdges());
    UInt index = 0;

    const auto edgeIndices = ComputeEdgeIndices();

    for (const auto& edgeIndex : edgeIndices)
    {
        const auto first = edgeIndex[0];
        const auto second = edgeIndex[1];
        result[index] = (GetNode(first.m_n, first.m_m) + GetNode(second.m_n, second.m_m)) * 0.5;
        index++;
    }

    return result;
}

std::vector<CurvilinearGridNodeIndices> CurvilinearGrid::ComputeNodeIndices() const
{
    std::vector result(NumN() * NumM(),
                       CurvilinearGridNodeIndices{constants::missing::uintValue,
                                                  constants::missing::uintValue});

    UInt ind = 0;
    for (UInt n = 0; n < NumN(); n++)
    {
        for (UInt m = 0; m < NumM(); m++)
        {

            result[ind] = {n, m};
            ind++;
        }
    }
    return result;
}

std::vector<CurvilinearGrid::CurvilinearEdgeNodeIndices> CurvilinearGrid::ComputeEdgeIndices() const
{

    std::vector<CurvilinearEdgeNodeIndices> result(GetNumEdges());
    UInt index = 0;
    for (UInt n = 0; n < NumN() - 1; n++)
    {
        for (UInt m = 0; m < NumM(); m++)
        {
            result[index][0] = {n, m};
            result[index][1] = {n + 1, m};
            index++;
        }
    }
    for (UInt n = 0; n < NumN(); n++)
    {
        for (UInt m = 0; m < NumM() - 1; m++)
        {
            result[index][0] = {n, m};
            result[index][1] = {n, m + 1};
            index++;
        }
    }
    return result;
}

std::set<CurvilinearGrid::CurvilinearEdge> CurvilinearGrid::ComputeBoundaryEdges(const CurvilinearGridNodeIndices& lowerLeft, const CurvilinearGridNodeIndices& upperRight) const
{
    std::vector<Point> result;
    std::set<CurvilinearEdge> boundaryEdges;
    for (UInt n = lowerLeft.m_n; n < upperRight.m_n; n++)
    {
        for (UInt m = lowerLeft.m_m; m < upperRight.m_m; m++)
        {
            CurvilinearFaceNodeIndices faceIndices;

            faceIndices[0] = {n, m};
            faceIndices[1] = {n, m + 1};
            faceIndices[2] = {n + 1, m + 1};
            faceIndices[3] = {n + 1, m};

            if (!std::ranges::all_of(faceIndices, [this](const CurvilinearGridNodeIndices& curviNode)
                                     { return GetNode(curviNode.m_n, curviNode.m_m).IsValid(); }))
            {
                continue;
            }

            for (UInt i = 0u; i < constants::geometric::numNodesInQuadrilateral; ++i)
            {
                const auto firstCurvilinearNodeIndex = faceIndices[i];
                const auto nextIndex = NextCircularForwardIndex(i, constants::geometric::numNodesInQuadrilateral);
                const auto secondCurvilinearNodeIndex = faceIndices[nextIndex];

                const auto edge = std::make_pair(std::min(firstCurvilinearNodeIndex, secondCurvilinearNodeIndex),
                                                 std::max(firstCurvilinearNodeIndex, secondCurvilinearNodeIndex));

                const auto it = boundaryEdges.find(edge);
                if (it != boundaryEdges.end())
                {
                    boundaryEdges.erase(it);
                }
                else
                {
                    boundaryEdges.insert(edge);
                }
            }
        }
    }

    return boundaryEdges;
}

std::vector<meshkernel::Point> CurvilinearGrid::ComputeBoundaryPolygons(const CurvilinearGridNodeIndices& lowerLeft, const CurvilinearGridNodeIndices& upperRight) const
{
    std::vector<Point> result;

    auto boundaryEdges = ComputeBoundaryEdges(lowerLeft, upperRight);
    if (boundaryEdges.empty())
    {
        return result;
    }

    auto currentEdge = boundaryEdges.begin();
    auto startNode = currentEdge->first;
    auto sharedNode = currentEdge->second;

    // iterate over all boundary edges, until all boundary polygons are found
    while (!boundaryEdges.empty())
    {
        result.push_back(GetNode(startNode.m_n, startNode.m_m));
        // iterate over connected edges, until the start node is found and the current boundary polygon is closed
        while (sharedNode != startNode)
        {
            result.push_back(GetNode(sharedNode.m_n, sharedNode.m_m));

            boundaryEdges.erase(currentEdge);
            currentEdge = std::ranges::find_if(boundaryEdges, [&](const CurvilinearEdge& edge)
                                               { return edge.first == sharedNode ||
                                                        edge.second == sharedNode; });

            sharedNode = currentEdge->first == sharedNode ? currentEdge->second : currentEdge->first;
        }
        result.push_back(GetNode(startNode.m_n, startNode.m_m));
        boundaryEdges.erase(currentEdge);
        if (boundaryEdges.empty())
        {
            return result;
        }

        result.emplace_back(constants::missing::doubleValue, constants::missing::doubleValue);

        // startNode node for the next boundary polygon
        currentEdge = boundaryEdges.begin();
        startNode = currentEdge->first;
        sharedNode = currentEdge->second;
    }

    return result;
}

std::vector<meshkernel::Point> CurvilinearGrid::ComputeFaceCenters() const
{
    std::vector<Point> result;
    result.reserve((NumN() - 1) * (NumM() - 1));

    for (UInt n = 0; n < NumN() - 1; n++)
    {
        for (UInt m = 0; m < NumM() - 1; m++)
        {
            Point massCenter{0.0, 0.0};

            massCenter += GetNode(n, m);
            massCenter += GetNode(n, m + 1);
            massCenter += GetNode(n + 1, m + 1);
            massCenter += GetNode(n + 1, m);

            result.push_back(massCenter * 0.25);
        }
    }
    return result;
}
