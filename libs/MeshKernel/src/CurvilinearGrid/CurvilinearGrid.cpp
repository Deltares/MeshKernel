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

#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridLine.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Operations.hpp"
#include "MeshKernel/Polygons.hpp"
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

void CurvilinearGrid::SetGridNodes(const lin_alg::Matrix<Point>& gridNodes)
{
    m_gridNodes = gridNodes;

    if (!IsValid())
    {
        throw std::invalid_argument("CurvilinearGrid::CurvilinearGrid: Invalid curvilinear grid");
    }

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
    if (lin_alg::MatrixIsEmpty(m_gridNodes))
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
            const auto isInPolygon = polygons->IsPointInPolygon(m_gridNodes(n, m), polygonIndex);
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
                m_gridNodes(n, m).x = constants::missing::doubleValue;
                m_gridNodes(n, m).y = constants::missing::doubleValue;
            }
        }
    }
}

bool CurvilinearGrid::IsValid() const
{
    if (m_gridNodes.size() == 0)
    {
        return false;
    }
    if (NumM() < 2)
    {
        return false;
    }
    if (NumN() < 2)
    {
        return false;
    }

    return true;
}

void CurvilinearGrid::BuildTree(Location location)
{
    switch (location)
    {
    case Location::Faces:
        if (m_facesRTreeRequiresUpdate)
        {
            const auto faceCenters = ComputeFaceCenters();
            m_RTrees.at(Location::Faces)->BuildTree(faceCenters);
            m_facesRTreeRequiresUpdate = false;
        }
        break;
    case Location::Nodes:
        if (m_nodesRTreeRequiresUpdate)
        {
            const auto nodes = ComputeNodes();
            m_RTrees.at(Location::Nodes)->BuildTree(nodes);
            m_nodesRTreeRequiresUpdate = false;
        }
        break;
    case Location::Edges:
        if (m_edgesRTreeRequiresUpdate)
        {
            m_edges = ComputeEdges();
            const auto edgeCenters = ComputeEdgesCenters();
            m_RTrees.at(Location::Edges)->BuildTree(edgeCenters);
            m_edgesRTreeRequiresUpdate = false;
        }
        break;
    case Location::Unknown:
    default:
        throw std::runtime_error("Mesh2D::SearchLocations: Mesh location has not been set.");
    }
}

meshkernel::UInt CurvilinearGrid::FindLocationIndex(Point point,
                                                    Location location,
                                                    const std::vector<bool>& locationMask)
{
    BuildTree(location);
    const auto& rtree = m_RTrees.at(location);
    if (rtree->Empty())
    {
        return constants::missing::uintValue;
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

CurvilinearGridNodeIndices CurvilinearGrid::GetNodeIndices(Point point)
{
    BuildTree(Location::Nodes);

    m_RTrees.at(Location::Nodes)->SearchNearestPoint(point);

    if (m_RTrees.at(Location::Nodes)->GetQueryResultSize() == 0)
    {
        return {constants::missing::uintValue, constants::missing::uintValue};
    }

    const auto nodeIndex = m_RTrees.at(Location::Nodes)->GetQueryResult(0);
    return m_gridIndices[nodeIndex];
}

std::tuple<CurvilinearGridNodeIndices, CurvilinearGridNodeIndices> CurvilinearGrid::GetEdgeNodeIndices(Point const& point)
{
    BuildTree(Location::Edges);

    m_RTrees.at(Location::Edges)->SearchNearestPoint(point);

    if (m_RTrees.at(Location::Edges)->GetQueryResultSize() == 0)
    {
        return {{}, {}};
    }

    const auto nodeIndex = m_RTrees.at(Location::Edges)->GetQueryResult(0);
    auto const firstNode = m_edges[nodeIndex].first;
    auto const secondNode = m_edges[nodeIndex].second;

    return {m_gridIndices[firstNode], m_gridIndices[secondNode]};
}

bool CurvilinearGrid::AreFaceNodesValid(UInt n, UInt m) const
{
    return m_gridNodes(n, m).IsValid() &&
           m_gridNodes(n, m + 1).IsValid() &&
           m_gridNodes(n + 1, m).IsValid() &&
           m_gridNodes(n + 1, m + 1).IsValid();
}

std::tuple<CurvilinearGridNodeIndices, CurvilinearGridNodeIndices> CurvilinearGrid::ComputeBlockFromCornerPoints(Point const& firstCornerPoint, Point const& secondCornerPoint)
{
    // Get the m and n indices from the point coordinates
    auto const firstNode = GetNodeIndices(firstCornerPoint);
    auto const secondNode = GetNodeIndices(secondCornerPoint);

    // Compute bounding box as node indices from corner points
    return ComputeBlockFromCornerPoints(firstNode, secondNode);
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
    lin_alg::ResizeAndFillMatrix(m_gridFacesMask, NumN() - 1, NumM() - 1, false, false);
    for (UInt n = 0; n < NumN() - 1; ++n)
    {
        for (UInt m = 0; m < NumM() - 1; ++m)
        {
            // Only if all grid nodes of the face are valid, the face is valid
            if (!AreFaceNodesValid(n, m))
            {
                continue;
            }
            m_gridFacesMask(n, m) = true;
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
            if (m_gridFacesMask(n, m))
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
            if (!validNodeMask[n][m] && m_gridNodes(n, m).IsValid())
            {
                m_gridNodes(n, m) = {constants::missing::doubleValue, constants::missing::doubleValue};
                invalidNodesToRemove = true;
            }
        }
    }

    RemoveInvalidNodes(invalidNodesToRemove);
}

void CurvilinearGrid::ComputeGridNodeTypes()
{
    RemoveInvalidNodes(true);
    lin_alg::ResizeAndFillMatrix(m_gridNodesTypes, NumN(), NumM(), false, NodeType::Invalid);

    // Flag faces based on boundaries
    for (UInt n = 0; n < NumN(); ++n)
    {
        for (UInt m = 0; m < NumM(); ++m)
        {

            if (!m_gridNodes(n, m).IsValid())
            {
                continue;
            }

            // Bottom side
            if (m == 0 && n == 0)
            {
                m_gridNodesTypes(n, m) = NodeType::BottomLeft;
                continue;
            }
            if (m == 0 && n == NumN() - 1)
            {
                m_gridNodesTypes(n, m) = NodeType::BottomRight;
                continue;
            }
            if (m == 0 && !m_gridNodes(n - 1, m).IsValid())
            {
                m_gridNodesTypes(n, m) = NodeType::BottomLeft;
                continue;
            }
            if (m == 0 && !m_gridNodes(n + 1, m).IsValid())
            {
                m_gridNodesTypes(n, m) = NodeType::BottomRight;
                continue;
            }
            if (m == 0)
            {
                m_gridNodesTypes(n, m) = NodeType::Bottom;
                continue;
            }
            // Upper side
            if (m == NumM() - 1 && n == 0)
            {
                m_gridNodesTypes(n, m) = NodeType::UpperLeft;
                continue;
            }
            if (m == NumM() - 1 && n == NumN() - 1)
            {
                m_gridNodesTypes(n, m) = NodeType::UpperRight;
                continue;
            }
            if (m == NumM() - 1 && !m_gridNodes(n - 1, m).IsValid())
            {
                m_gridNodesTypes(n, m) = NodeType::UpperLeft;
                continue;
            }
            if (m == NumM() - 1 && !m_gridNodes(n + 1, m).IsValid())
            {
                m_gridNodesTypes(n, m) = NodeType::UpperRight;
                continue;
            }
            if (m == NumM() - 1)
            {
                m_gridNodesTypes(n, m) = NodeType::Up;
                continue;
            }
            // Bottom side
            if (n == 0 && !m_gridNodes(n, m - 1).IsValid())
            {
                m_gridNodesTypes(n, m) = NodeType::BottomLeft;
                continue;
            }
            if (n == 0 && !m_gridNodes(n, m + 1).IsValid())
            {
                m_gridNodesTypes(n, m) = NodeType::UpperRight;
                continue;
            }
            if (n == 0)
            {
                m_gridNodesTypes(n, m) = NodeType::Left;
                continue;
            }
            // Upper side
            if (n == NumN() - 1 && !m_gridNodes(n, m - 1).IsValid())
            {
                m_gridNodesTypes(n, m) = NodeType::BottomRight;
                continue;
            }
            if (n == NumN() - 1 && !m_gridNodes(n, m + 1).IsValid())
            {
                m_gridNodesTypes(n, m) = NodeType::UpperRight;
                continue;
            }
            if (n == NumN() - 1)
            {
                m_gridNodesTypes(n, m) = NodeType::Right;
                continue;
            }

            auto const isBottomRightFaceValid = m_gridFacesMask(n, m - 1);
            auto const isBottomLeftFaceValid = m_gridFacesMask(n - 1, m - 1);
            auto const isTopRightFaceValid = m_gridFacesMask(n, m);
            auto const isTopLeftFaceValid = m_gridFacesMask(n - 1, m);

            if (isTopRightFaceValid &&
                isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                m_gridNodesTypes(n, m) = NodeType::InternalValid;
                continue;
            }
            if (!isTopRightFaceValid &&
                isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                m_gridNodesTypes(n, m) = NodeType::BottomLeft;
                continue;
            }
            if (isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                m_gridNodesTypes(n, m) = NodeType::BottomRight;
                continue;
            }
            if (isTopRightFaceValid &&
                isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                m_gridNodesTypes(n, m) = NodeType::UpperRight;
                continue;
            }
            if (isTopRightFaceValid &&
                isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                m_gridNodesTypes(n, m) = NodeType::UpperLeft;
                continue;
            }

            if (isTopRightFaceValid &&
                isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                m_gridNodesTypes(n, m) = NodeType::Bottom;
                continue;
            }
            if (isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                m_gridNodesTypes(n, m) = NodeType::Left;
                continue;
            }

            if (!isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                m_gridNodesTypes(n, m) = NodeType::Up;
                continue;
            }

            if (!isTopRightFaceValid &&
                isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                m_gridNodesTypes(n, m) = NodeType::Right;
                continue;
            }

            if (isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                m_gridNodesTypes(n, m) = NodeType::BottomLeft;
                continue;
            }

            if (!isTopRightFaceValid &&
                isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                m_gridNodesTypes(n, m) = NodeType::BottomRight;
                continue;
            }

            if (!isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                m_gridNodesTypes(n, m) = NodeType::UpperRight;
                continue;
            }

            if (!isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                m_gridNodesTypes(n, m) = NodeType::UpperLeft;
            }
        }
    }
}

void CurvilinearGrid::InsertFace(Point const& point)
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
    AddEdge(firstNode, secondNode);

    // Re-compute quantities
    ComputeGridNodeTypes();
    m_gridIndices = ComputeNodeIndices();
}

bool CurvilinearGrid::AddGridLineAtBoundary(CurvilinearGridNodeIndices const& firstNode,
                                            CurvilinearGridNodeIndices const& secondNode)
{
    // If both nodes are invalid, we can substitute the invalid values. New allocation is not needed.
    bool const areNodesValid = m_gridNodes(firstNode.m_n, firstNode.m_m).IsValid() &&
                               m_gridNodes(secondNode.m_n, secondNode.m_m).IsValid();

    // Allocation depends on directions
    bool gridSizeChanged = false;
    auto const gridLineType = GetBoundaryGridLineType(firstNode, secondNode);

    if (areNodesValid)
    {

        if (gridLineType == BoundaryGridLineType::Left)
        {
            lin_alg::InsertRow(m_gridNodes,
                               lin_alg::RowVector<Point>(NumM()),
                               0);
            gridSizeChanged = true;
        }
        if (gridLineType == BoundaryGridLineType::Right)
        {
            lin_alg::InsertRow(m_gridNodes,
                               lin_alg::RowVector<Point>(NumM()),
                               NumN());
            gridSizeChanged = true;
        }
        if (gridLineType == BoundaryGridLineType::Up)
        {
            lin_alg::InsertCol(m_gridNodes,
                               lin_alg::ColVector<Point>(NumN()),
                               NumM());
            gridSizeChanged = true;
        }
        if (gridLineType == BoundaryGridLineType::Bottom)
        {
            lin_alg::InsertCol(m_gridNodes,
                               lin_alg::ColVector<Point>(NumN()),
                               0);
            gridSizeChanged = true;
        }
    }

    return gridSizeChanged;
}

CurvilinearGrid::BoundaryGridLineType CurvilinearGrid::GetBoundaryGridLineType(CurvilinearGridNodeIndices const& firstNode, CurvilinearGridNodeIndices const& secondNode) const
{
    auto const firstNodeType = m_gridNodesTypes(firstNode.m_n, firstNode.m_m);
    auto const secondNodeType = m_gridNodesTypes(secondNode.m_n, secondNode.m_m);

    if (firstNodeType == NodeType::InternalValid || firstNodeType == NodeType::Invalid ||
        secondNodeType == NodeType::InternalValid || secondNodeType == NodeType::Invalid)
    {
        throw std::invalid_argument("CurvilinearGrid::GetBoundaryGridLineType: Not a boundary grid line");
    }

    if (firstNodeType == NodeType::Bottom || secondNodeType == NodeType::Bottom ||
        (firstNodeType == NodeType::BottomLeft && secondNodeType == NodeType::BottomRight) ||
        (firstNodeType == NodeType::BottomRight && secondNodeType == NodeType::BottomLeft))
    {
        return BoundaryGridLineType::Bottom;
    }
    if (firstNodeType == NodeType::Up || secondNodeType == NodeType::Up ||
        (firstNodeType == NodeType::UpperLeft && secondNodeType == NodeType::UpperRight) ||
        (firstNodeType == NodeType::UpperRight && secondNodeType == NodeType::UpperLeft))
    {
        return BoundaryGridLineType::Up;
    }
    if (firstNodeType == NodeType::Left || secondNodeType == NodeType::Left ||
        (firstNodeType == NodeType::BottomLeft && secondNodeType == NodeType::UpperLeft) ||
        (firstNodeType == NodeType::UpperLeft && secondNodeType == NodeType::BottomLeft))
    {
        return BoundaryGridLineType::Left;
    }
    if (firstNodeType == NodeType::Right || secondNodeType == NodeType::Right ||
        (firstNodeType == NodeType::BottomRight && secondNodeType == NodeType::UpperRight) ||
        (firstNodeType == NodeType::UpperRight && secondNodeType == NodeType::BottomRight))
    {
        return BoundaryGridLineType::Right;
    }

    throw std::invalid_argument("CurvilinearGrid::GetBoundaryGridLineType: Invalid node types");
}

void CurvilinearGrid::AddEdge(CurvilinearGridNodeIndices const& firstNode, CurvilinearGridNodeIndices const& secondNode)
{

    // Allocate new grid line if needed
    auto const gridLineType = GetBoundaryGridLineType(firstNode, secondNode);
    // Allocation depends on directions
    if (gridLineType == BoundaryGridLineType::Left)
    {
        auto const firstNewNodeCoordinates = m_gridNodes(firstNode.m_n, firstNode.m_m) * 2.0 - m_gridNodes(firstNode.m_n + 1, firstNode.m_m);
        auto const secondNewNodeCoordinates = m_gridNodes(secondNode.m_n, secondNode.m_m) * 2.0 - m_gridNodes(secondNode.m_n + 1, secondNode.m_m);
        auto const isGridLineAdded = AddGridLineAtBoundary(firstNode, secondNode);
        if (isGridLineAdded)
        {
            m_gridNodes(0, firstNode.m_m) = firstNewNodeCoordinates;
            m_gridNodes(0, secondNode.m_m) = secondNewNodeCoordinates;
            return;
        }

        m_gridNodes(firstNode.m_n - 1, firstNode.m_m) = firstNewNodeCoordinates;
        m_gridNodes(secondNode.m_n - 1, secondNode.m_m) = secondNewNodeCoordinates;
        return;
    }
    if (gridLineType == BoundaryGridLineType::Right)
    {
        auto const firstNewNodeCoordinates = m_gridNodes(firstNode.m_n, firstNode.m_m) * 2.0 - m_gridNodes(firstNode.m_n - 1, firstNode.m_m);
        auto const secondNewNodeCoordinates = m_gridNodes(secondNode.m_n, secondNode.m_m) * 2.0 - m_gridNodes(secondNode.m_n - 1, secondNode.m_m);
        AddGridLineAtBoundary(firstNode, secondNode);
        m_gridNodes(firstNode.m_n + 1, firstNode.m_m) = firstNewNodeCoordinates;
        m_gridNodes(secondNode.m_n + 1, secondNode.m_m) = secondNewNodeCoordinates;

        return;
    }
    if (gridLineType == BoundaryGridLineType::Bottom)
    {
        auto const firstNewNodeCoordinates = m_gridNodes(firstNode.m_n, firstNode.m_m) * 2.0 - m_gridNodes(firstNode.m_n, firstNode.m_m + 1);
        auto const secondNewNodeCoordinates = m_gridNodes(secondNode.m_n, secondNode.m_m) * 2.0 - m_gridNodes(secondNode.m_n, secondNode.m_m + 1);
        auto const isGridLineAdded = AddGridLineAtBoundary(firstNode, secondNode);
        if (isGridLineAdded)
        {
            // Assign the new coordinates
            m_gridNodes(firstNode.m_n, 0) = firstNewNodeCoordinates;
            m_gridNodes(secondNode.m_n, 0) = secondNewNodeCoordinates;
            return;
        }

        m_gridNodes(firstNode.m_n, firstNode.m_m - 1) = firstNewNodeCoordinates;
        m_gridNodes(secondNode.m_n, secondNode.m_m - 1) = secondNewNodeCoordinates;
    }

    if (gridLineType == BoundaryGridLineType::Up)
    {
        auto const firstNewNodeCoordinates = m_gridNodes(firstNode.m_n, firstNode.m_m) * 2.0 - m_gridNodes(firstNode.m_n, firstNode.m_m - 1);
        auto const secondNewNodeCoordinates = m_gridNodes(secondNode.m_n, secondNode.m_m) * 2.0 - m_gridNodes(secondNode.m_n, secondNode.m_m - 1);
        AddGridLineAtBoundary(firstNode, secondNode);
        m_gridNodes(firstNode.m_n, firstNode.m_m + 1) = firstNewNodeCoordinates;
        m_gridNodes(secondNode.m_n, secondNode.m_m + 1) = secondNewNodeCoordinates;
    }
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
        if (index.m_m > 0 && m_gridNodes(index.m_n, index.m_m - 1).IsValid())
        {
            leftDistance = ComputeDistance(m_gridNodes(index.m_n, index.m_m), m_gridNodes(index.m_n, index.m_m - 1), m_projection);
            numEdges += 1;
        }
        if (index.m_m + 1 < NumM() && m_gridNodes(index.m_n, index.m_m + 1).IsValid())
        {
            rightDistance = ComputeDistance(m_gridNodes(index.m_n, index.m_m), m_gridNodes(index.m_n, index.m_m + 1), m_projection);
            numEdges += 1;
        }
        return numEdges == 0 ? 0.0 : (leftDistance + rightDistance) / numEdges;
    }
    if (direction == CurvilinearGridLine::GridLineDirection::NDirection)
    {
        double numEdges = 0.0;
        double bottomDistance = 0.0;
        double upDistance = 0.0;
        if (index.m_n > 0 && m_gridNodes(index.m_n - 1, index.m_m).IsValid())
        {
            bottomDistance = ComputeDistance(m_gridNodes(index.m_n, index.m_m), m_gridNodes(index.m_n - 1, index.m_m), m_projection);
            numEdges += 1;
        }
        if (index.m_n + 1 < NumN() && m_gridNodes(index.m_n + 1, index.m_m).IsValid())
        {
            upDistance = ComputeDistance(m_gridNodes(index.m_n, index.m_m), m_gridNodes(index.m_n + 1, index.m_m), m_projection);
            numEdges += 1;
        }
        return numEdges == 0 ? 0.0 : (bottomDistance + upDistance) / numEdges;
    }

    throw std::invalid_argument("CurvilinearGrid::ComputeAverageNodalDistance: Invalid direction");
}

meshkernel::Point CurvilinearGrid::TransformDisplacement(Point const& displacement, CurvilinearGridNodeIndices const& node, bool isLocal) const
{
    Point left = m_gridNodes(node.m_n, node.m_m);
    Point right = left;
    if (node.m_n < NumN() - 1 && m_gridNodes(node.m_n + 1, node.m_m).IsValid())
    {
        right = m_gridNodes(node.m_n + 1, node.m_m);
    }
    if (node.m_n > 0 && m_gridNodes(node.m_n - 1, node.m_m).IsValid())
    {
        left = m_gridNodes(node.m_n - 1, node.m_m);
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

void CurvilinearGrid::DeleteNode(Point const& point)
{
    // Get the m and n indices from the point coordinates
    auto const nodeToDelete = GetNodeIndices(point);

    if (nodeToDelete.IsValid())
    {
        // Invalidate gridnodes
        m_gridNodes(nodeToDelete.m_n, nodeToDelete.m_m) = {constants::missing::doubleValue, constants::missing::doubleValue};
        // Re-compute quantities
        ComputeGridNodeTypes();
        m_gridIndices = ComputeNodeIndices();
    }
}

void CurvilinearGrid::MoveNode(Point const& fromPoint, Point const& toPoint)
{
    // Get the node indices of fromPoint
    auto const nodeIndex = GetNodeIndices(fromPoint);

    // Check the node indices are valid
    if (!nodeIndex.IsValid())
    {
        throw std::invalid_argument("CurvilinearGrid::MoveNode node indices not found");
    }

    // move fromPoint to toPoint
    m_gridNodes(nodeIndex.m_n, nodeIndex.m_m) = toPoint;
}

meshkernel::BoundingBox CurvilinearGrid::GetBoundingBox() const
{

    Point lowerLeft(std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
    Point upperRight(-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());
    size_t last = NumM() - 1;

    // Only need to loop over boundary nodes

    // First loop over lower boundary (i,0)
    for (Eigen::Index i = 0; i < NumM(); ++i)
    {
        lowerLeft.x = std::min(lowerLeft.x, m_gridNodes(0, i).x);
        lowerLeft.y = std::min(lowerLeft.y, m_gridNodes(0, i).y);
        upperRight.x = std::max(upperRight.x, m_gridNodes(0, i).x);
        upperRight.y = std::max(upperRight.y, m_gridNodes(0, i).y);
    }

    // First loop over right boundary (last,i)
    for (Eigen::Index i = 0; i < NumN(); ++i)
    {
        lowerLeft.x = std::min(lowerLeft.x, m_gridNodes(i, last).x);
        lowerLeft.y = std::min(lowerLeft.y, m_gridNodes(i, last).y);
        upperRight.x = std::max(upperRight.x, m_gridNodes(i, last).x);
        upperRight.y = std::max(upperRight.y, m_gridNodes(i, last).y);
    }

    // This assumes that each column has the same number of points
    last = NumN() - 1;

    // First loop over upper boundary (i,last)
    for (Eigen::Index i = 0; i < NumM(); ++i)
    {
        lowerLeft.x = std::min(lowerLeft.x, m_gridNodes(last, i).x);
        lowerLeft.y = std::min(lowerLeft.y, m_gridNodes(last, i).y);
        upperRight.x = std::max(upperRight.x, m_gridNodes(last, i).x);
        upperRight.y = std::max(upperRight.y, m_gridNodes(last, i).y);
    }

    // First loop over left boundary (0,i)
    for (Eigen::Index i = 0; i < NumN(); ++i)
    {
        lowerLeft.x = std::min(lowerLeft.x, m_gridNodes(i, 0).x);
        lowerLeft.y = std::min(lowerLeft.y, m_gridNodes(i, 0).y);
        upperRight.x = std::max(upperRight.x, m_gridNodes(i, 0).x);
        upperRight.y = std::max(upperRight.y, m_gridNodes(i, 0).y);
    }

    return BoundingBox(lowerLeft, upperRight);
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

            result[ind] = m_gridNodes(n, m);
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
    if (!IsValid())
    {
        throw AlgorithmError("Invalid curvilinear grid ");
    }

    std::vector<Point> edgesCenters(GetNumEdges());
    UInt index = 0;
    for (UInt n = 0; n < NumN() - 1; n++)
    {
        for (UInt m = 0; m < NumM(); m++)
        {
            edgesCenters[index] = (m_gridNodes(n, m) + m_gridNodes(n + 1, m)) * 0.5;
            index++;
        }
    }
    for (UInt n = 0; n < NumN(); n++)
    {
        for (UInt m = 0; m < NumM() - 1; m++)
        {
            edgesCenters[index] = (m_gridNodes(n, m) + m_gridNodes(n, m + 1)) * 0.5;
            index++;
        }
    }

    return edgesCenters;
}

std::vector<CurvilinearGridNodeIndices> CurvilinearGrid::ComputeNodeIndices() const
{
    if (!IsValid())
    {
        throw AlgorithmError("Invalid curvilinear grid ");
    }

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

std::vector<meshkernel::Point> CurvilinearGrid::ComputeFaceCenters() const
{
    if (!IsValid())
    {
        throw AlgorithmError("Invalid curvilinear grid ");
    }

    std::vector<Point> result((NumM() - 1) * (NumN() - 1));

    UInt index = 0;
    for (UInt n = 0; n < NumN() - 1; n++)
    {
        for (UInt m = 0; m < NumM() - 1; m++)
        {
            result[index] = (m_gridNodes(n, m) +
                             m_gridNodes(n, m + 1) +
                             m_gridNodes(n + 1, m + 1) +
                             m_gridNodes(n + 1, m)) *
                            0.25;
            index++;
        }
    }

    return result;
}
