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

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLine.hpp>
#include <MeshKernel/Operations.hpp>

using meshkernel::CurvilinearGrid;
using meshkernel::CurvilinearGridNodeIndices;

CurvilinearGrid::CurvilinearGrid(std::vector<std::vector<Point>> const& grid, Projection projection) : m_gridNodes(grid)
{
    if (!IsValid())
    {
        throw std::invalid_argument("CurvilinearGrid::CurvilinearGrid: Invalid curvilinear grid ");
    }

    m_projection = projection;
    m_numM = m_gridNodes.size();
    m_numN = m_gridNodes[0].size();

    SetFlatCopies();
}

void CurvilinearGrid::SetFlatCopies()
{
    m_numM = m_gridNodes.size();
    m_numN = m_gridNodes[0].size();
    const auto [nodes, edges, gridIndices] = ConvertCurvilinearToNodesAndEdges();
    m_nodes = nodes;
    m_edges = edges;
    m_gridIndices = gridIndices;
}

std::tuple<std::vector<meshkernel::Point>, std::vector<meshkernel::Edge>, std::vector<CurvilinearGridNodeIndices>> CurvilinearGrid::ConvertCurvilinearToNodesAndEdges()
{
    if (!IsValid())
    {
        throw std::invalid_argument("CurvilinearGrid::ConvertCurvilinearToNodesAndEdges: Invalid curvilinear grid ");
    }

    std::vector<Point> nodes(m_gridNodes.size() * m_gridNodes[0].size());
    std::vector<Edge> edges(m_gridNodes.size() * (m_gridNodes[0].size() - 1) + (m_gridNodes.size() - 1) * m_gridNodes[0].size());
    std::vector<std::vector<size_t>> nodeIndices(m_gridNodes.size(), std::vector<size_t>(m_gridNodes[0].size(), sizetMissingValue));
    std::vector<CurvilinearGridNodeIndices> gridIndices(nodes.size(), CurvilinearGridNodeIndices{sizetMissingValue, sizetMissingValue});

    size_t ind = 0;
    for (size_t m = 0; m < m_gridNodes.size(); m++)
    {
        for (size_t n = 0; n < m_gridNodes[0].size(); n++)
        {
            if (m_gridNodes[m][n].IsValid())
            {
                nodes[ind] = m_gridNodes[m][n];
                nodeIndices[m][n] = ind;
                gridIndices[ind] = {m, n};
                ind++;
            }
        }
    }
    nodes.resize(ind);

    ind = 0;
    for (auto m = 0; m < m_gridNodes.size() - 1; m++)
    {
        for (auto n = 0; n < m_gridNodes[0].size(); n++)
        {
            if (nodeIndices[m][n] != sizetMissingValue && nodeIndices[m + 1][n] != sizetMissingValue)
            {
                edges[ind].first = nodeIndices[m][n];
                edges[ind].second = nodeIndices[m + 1][n];
                ind++;
            }
        }
    }

    for (auto m = 0; m < m_gridNodes.size(); m++)
    {
        for (auto n = 0; n < m_gridNodes[0].size() - 1; n++)
        {
            if (nodeIndices[m][n] != sizetMissingValue && nodeIndices[m][n + 1] != sizetMissingValue)
            {
                edges[ind].first = nodeIndices[m][n];
                edges[ind].second = nodeIndices[m][n + 1];
                ind++;
            }
        }
    }
    edges.resize(ind);

    return {nodes, edges, gridIndices};
}

bool CurvilinearGrid::IsValid() const
{
    if (m_gridNodes.empty())
    {
        return false;
    }
    if (m_gridNodes[0].empty())
    {
        return false;
    }
    if (m_gridNodes.size() < 2)
    {
        return false;
    }
    if (m_gridNodes[0].size() < 2)
    {
        return false;
    }

    return true;
}

CurvilinearGridNodeIndices CurvilinearGrid::GetNodeIndices(Point point)
{
    SearchNearestNeighbors(point, MeshLocations::Nodes);
    if (GetNumNearestNeighbors(MeshLocations::Nodes) == 0)
    {
        return {sizetMissingValue, sizetMissingValue};
    }

    const auto nodeIndex = GetNearestNeighborIndex(0, MeshLocations::Nodes);
    return m_gridIndices[nodeIndex];
}

std::tuple<CurvilinearGridNodeIndices, CurvilinearGridNodeIndices> CurvilinearGrid::GetEdgeNodeIndices(Point const& point)
{
    SearchNearestNeighbors(point, MeshLocations::Edges);
    if (GetNumNearestNeighbors(MeshLocations::Edges) == 0)
    {
        return {{}, {}};
    }

    const auto nodeIndex = GetNearestNeighborIndex(0, MeshLocations::Edges);
    auto const firstNode = m_edges[nodeIndex].first;
    auto const secondNode = m_edges[nodeIndex].second;

    return {m_gridIndices[firstNode], m_gridIndices[secondNode]};
}

bool CurvilinearGrid::IsValidFace(size_t m, size_t n) const
{
    return m_gridNodes[m][n].IsValid() &&
           m_gridNodes[m + 1][n].IsValid() &&
           m_gridNodes[m][n + 1].IsValid() &&
           m_gridNodes[m + 1][n + 1].IsValid();
};

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
    return {{std::min(firstNode.m_m, secondNode.m_m), std::min(firstNode.m_n, secondNode.m_n)},
            {std::max(firstNode.m_m, secondNode.m_m), std::max(firstNode.m_n, secondNode.m_n)}};
}

void CurvilinearGrid::ComputeGridFacesMask()
{
    // Flag valid faces
    m_gridFacesMask.resize(m_numM - 1, std::vector<bool>(m_numN - 1, false));
    for (auto m = 0; m < m_numM - 1; ++m)
    {
        for (auto n = 0; n < m_numN - 1; ++n)
        {
            // Only if all grid nodes of the face are valid, the face is valid
            if (!IsValidFace(m, n))
            {
                continue;
            }
            m_gridFacesMask[m][n] = true;
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

    invalidNodesToRemove = false;
    // Flag nodes not connected to valid faces
    for (auto m = 1; m < m_numM - 1; ++m)
    {
        for (auto n = 1; n < m_numN - 1; ++n)
        {
            if (m_gridNodes[m][n].IsValid() &&
                !m_gridFacesMask[m][n] &&
                !m_gridFacesMask[m - 1][n] &&
                !m_gridFacesMask[m - 1][n - 1] &&
                !m_gridFacesMask[m][n - 1])
            {
                m_gridNodes[m][n] = {doubleMissingValue, doubleMissingValue};
                invalidNodesToRemove = true;
            }
        }
    }

    for (auto m = 1; m < m_numM - 1; ++m)
    {
        if (m_gridNodes[m][0].IsValid() &&
            !m_gridFacesMask[m - 1][0] &&
            !m_gridFacesMask[m][0])
        {
            m_gridNodes[m][0] = {doubleMissingValue, doubleMissingValue};
        }
    }

    for (auto n = 1; n < m_numN - 1; ++n)
    {
        if (m_gridNodes[0][n].IsValid() &&
            !m_gridFacesMask[0][n - 1] &&
            !m_gridFacesMask[0][n])
        {
            m_gridNodes[0][n] = {doubleMissingValue, doubleMissingValue};
        }
    }

    if (m_gridNodes[0][0].IsValid() && !m_gridFacesMask[0][0])
    {
        m_gridNodes[0][0] = {doubleMissingValue, doubleMissingValue};
    }

    RemoveInvalidNodes(invalidNodesToRemove);
}

void CurvilinearGrid::ComputeGridNodeTypes()
{
    RemoveInvalidNodes(true);
    m_gridNodesTypes.resize(m_numM, std::vector<NodeType>(m_numN, NodeType::Invalid));

    // Flag faces based on boundaries
    for (size_t m = 0; m < m_numM; ++m)
    {
        for (size_t n = 0; n < m_numN; ++n)
        {

            if (!m_gridNodes[m][n].IsValid())
            {
                continue;
            }

            // Left side
            if (m == 0 && n == 0)
            {
                m_gridNodesTypes[m][n] = NodeType::BottomLeft;
                continue;
            }
            if (m == 0 && n == m_numN - 1)
            {
                m_gridNodesTypes[m][n] = NodeType::UpperLeft;
                continue;
            }
            if (m == 0 && !m_gridNodes[m][n - 1].IsValid())
            {
                m_gridNodesTypes[m][n] = NodeType::BottomLeft;
                continue;
            }
            if (m == 0 && !m_gridNodes[m][n + 1].IsValid())
            {
                m_gridNodesTypes[m][n] = NodeType::UpperLeft;
                continue;
            }
            if (m == 0)
            {
                m_gridNodesTypes[m][n] = NodeType::Left;
                continue;
            }
            // Right side
            if (m == m_numM - 1 && n == 0)
            {
                m_gridNodesTypes[m][n] = NodeType::BottomRight;
                continue;
            }
            if (m == m_numM - 1 && n == m_numN - 1)
            {
                m_gridNodesTypes[m][n] = NodeType::UpperRight;
                continue;
            }
            if (m == m_numM - 1 && !m_gridNodes[m][n - 1].IsValid())
            {
                m_gridNodesTypes[m][n] = NodeType::BottomRight;
                continue;
            }
            if (m == m_numM - 1 && !m_gridNodes[m][n + 1].IsValid())
            {
                m_gridNodesTypes[m][n] = NodeType::UpperRight;
                continue;
            }
            if (m == m_numM - 1)
            {
                m_gridNodesTypes[m][n] = NodeType::Right;
                continue;
            }
            // Bottom side
            if (n == 0 && !m_gridNodes[m - 1][n].IsValid())
            {
                m_gridNodesTypes[m][n] = NodeType::BottomLeft;
                continue;
            }
            if (n == 0 && !m_gridNodes[m + 1][n].IsValid())
            {
                m_gridNodesTypes[m][n] = NodeType::BottomRight;
                continue;
            }
            if (n == 0)
            {
                m_gridNodesTypes[m][n] = NodeType::Bottom;
                continue;
            }
            // Upper side
            if (n == m_numN - 1 && !m_gridNodes[m - 1][n].IsValid())
            {
                m_gridNodesTypes[m][n] = NodeType::UpperLeft;
                continue;
            }
            if (n == m_numN - 1 && !m_gridNodes[m + 1][n].IsValid())
            {
                m_gridNodesTypes[m][n] = NodeType::UpperRight;
                continue;
            }
            if (n == m_numN - 1)
            {
                m_gridNodesTypes[m][n] = NodeType::Up;
                continue;
            }

            auto const isTopLeftFaceValid = m_gridFacesMask[m - 1][n];
            auto const isTopRightFaceValid = m_gridFacesMask[m][n];
            auto const isBottomLeftFaceValid = m_gridFacesMask[m - 1][n - 1];
            auto const isBottomRightFaceValid = m_gridFacesMask[m][n - 1];

            if (isTopRightFaceValid &&
                isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                m_gridNodesTypes[m][n] = NodeType::InternalValid;
                continue;
            }
            if (!isTopRightFaceValid &&
                isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                m_gridNodesTypes[m][n] = NodeType::BottomLeft;
                continue;
            }
            if (isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                m_gridNodesTypes[m][n] = NodeType::BottomRight;
                continue;
            }
            if (isTopRightFaceValid &&
                isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                m_gridNodesTypes[m][n] = NodeType::UpperRight;
                continue;
            }
            if (isTopRightFaceValid &&
                isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                m_gridNodesTypes[m][n] = NodeType::UpperLeft;
                continue;
            }

            if (isTopRightFaceValid &&
                isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                m_gridNodesTypes[m][n] = NodeType::Bottom;
                continue;
            }
            if (isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                m_gridNodesTypes[m][n] = NodeType::Left;
                continue;
            }

            if (!isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                m_gridNodesTypes[m][n] = NodeType::Up;
                continue;
            }

            if (!isTopRightFaceValid &&
                isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                m_gridNodesTypes[m][n] = NodeType::Right;
                continue;
            }

            if (isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                m_gridNodesTypes[m][n] = NodeType::BottomLeft;
                continue;
            }

            if (!isTopRightFaceValid &&
                isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                m_gridNodesTypes[m][n] = NodeType::BottomRight;
                continue;
            }

            if (!isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                m_gridNodesTypes[m][n] = NodeType::UpperRight;
                continue;
            }

            if (!isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                m_gridNodesTypes[m][n] = NodeType::UpperLeft;
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
    SetFlatCopies();
}

bool CurvilinearGrid::AddGridLineAtBoundary(CurvilinearGridNodeIndices const& firstNode, CurvilinearGridNodeIndices const& secondNode)
{
    // If both nodes are invalid, we can substitute the invalid values. New allocation is not needed.
    bool const areNodesValid = m_gridNodes[firstNode.m_m][firstNode.m_n].IsValid() && m_gridNodes[secondNode.m_m][secondNode.m_n].IsValid();

    // Allocation depends on directions
    bool gridSizeChanged = false;
    auto const gridLineType = GetBoundaryGridLineType(firstNode, secondNode);
    if (gridLineType == BoundaryGridLineType::Left && areNodesValid)
    {
        m_gridNodes.emplace(m_gridNodes.begin(), std::vector<Point>(m_gridNodes[0].size()));
        gridSizeChanged = true;
    }
    if (gridLineType == BoundaryGridLineType::Right && areNodesValid)
    {
        m_gridNodes.emplace_back(std::vector<Point>(m_gridNodes[0].size()));
        gridSizeChanged = true;
    }
    if (gridLineType == BoundaryGridLineType::Up && areNodesValid)
    {
        for (auto& gridNodes : m_gridNodes)
        {
            gridNodes.emplace_back();
        }
        gridSizeChanged = true;
    }
    if (gridLineType == BoundaryGridLineType::Bottom && areNodesValid)
    {
        for (auto& gridNodes : m_gridNodes)
        {
            gridNodes.emplace(gridNodes.begin());
        }
        gridSizeChanged = true;
    }

    return gridSizeChanged;
}

CurvilinearGrid::BoundaryGridLineType CurvilinearGrid::GetBoundaryGridLineType(CurvilinearGridNodeIndices const& firstNode, CurvilinearGridNodeIndices const& secondNode) const
{
    auto const firstNodeType = m_gridNodesTypes[firstNode.m_m][firstNode.m_n];
    auto const secondNodeType = m_gridNodesTypes[secondNode.m_m][secondNode.m_n];

    if (firstNodeType == NodeType::InternalValid || firstNodeType == NodeType::Invalid ||
        secondNodeType == NodeType::InternalValid || secondNodeType == NodeType::Invalid)
    {
        throw std::invalid_argument("CurvilinearGrid::GetBoundaryGridLineType: Not a boundary grid line");
    }

    if (firstNodeType == NodeType::Bottom || secondNodeType == NodeType::Bottom ||
        firstNodeType == NodeType::BottomLeft && secondNodeType == NodeType::BottomRight ||
        firstNodeType == NodeType::BottomRight && secondNodeType == NodeType::BottomLeft)
    {
        return BoundaryGridLineType::Bottom;
    }
    if (firstNodeType == NodeType::Up || secondNodeType == NodeType::Up ||
        firstNodeType == NodeType::UpperLeft && secondNodeType == NodeType::UpperRight ||
        firstNodeType == NodeType::UpperRight && secondNodeType == NodeType::UpperLeft)
    {
        return BoundaryGridLineType::Up;
    }
    if (firstNodeType == NodeType::Left || secondNodeType == NodeType::Left ||
        firstNodeType == NodeType::BottomLeft && secondNodeType == NodeType::UpperLeft ||
        firstNodeType == NodeType::UpperLeft && secondNodeType == NodeType::BottomLeft)
    {
        return BoundaryGridLineType::Left;
    }
    if (firstNodeType == NodeType::Right || secondNodeType == NodeType::Right ||
        firstNodeType == NodeType::BottomRight && secondNodeType == NodeType::UpperRight ||
        firstNodeType == NodeType::UpperRight && secondNodeType == NodeType::BottomRight)
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
        auto const firstNewNodeCoordinates = m_gridNodes[firstNode.m_m][firstNode.m_n] * 2.0 - m_gridNodes[firstNode.m_m + 1][firstNode.m_n];
        auto const secondNewNodeCoordinates = m_gridNodes[secondNode.m_m][secondNode.m_n] * 2.0 - m_gridNodes[secondNode.m_m + 1][secondNode.m_n];
        auto const isGridLineAdded = AddGridLineAtBoundary(firstNode, secondNode);
        if (isGridLineAdded)
        {
            m_gridNodes.front()[firstNode.m_n] = firstNewNodeCoordinates;
            m_gridNodes.front()[secondNode.m_n] = secondNewNodeCoordinates;
            return;
        }

        m_gridNodes[firstNode.m_m - 1][firstNode.m_n] = firstNewNodeCoordinates;
        m_gridNodes[secondNode.m_m - 1][secondNode.m_n] = secondNewNodeCoordinates;
        return;
    }
    if (gridLineType == BoundaryGridLineType::Right)
    {
        auto const firstNewNodeCoordinates = m_gridNodes[firstNode.m_m][firstNode.m_n] * 2.0 - m_gridNodes[firstNode.m_m - 1][firstNode.m_n];
        auto const secondNewNodeCoordinates = m_gridNodes[secondNode.m_m][secondNode.m_n] * 2.0 - m_gridNodes[secondNode.m_m - 1][secondNode.m_n];
        AddGridLineAtBoundary(firstNode, secondNode);
        m_gridNodes[firstNode.m_m + 1][firstNode.m_n] = firstNewNodeCoordinates;
        m_gridNodes[secondNode.m_m + 1][secondNode.m_n] = secondNewNodeCoordinates;

        return;
    }
    if (gridLineType == BoundaryGridLineType::Bottom)
    {
        auto const firstNewNodeCoordinates = m_gridNodes[firstNode.m_m][firstNode.m_n] * 2.0 - m_gridNodes[firstNode.m_m][firstNode.m_n + 1];
        auto const secondNewNodeCoordinates = m_gridNodes[secondNode.m_m][secondNode.m_n] * 2.0 - m_gridNodes[secondNode.m_m][secondNode.m_n + 1];
        auto const isGridLineAdded = AddGridLineAtBoundary(firstNode, secondNode);
        if (isGridLineAdded)
        {
            // Assign the new coordinates
            m_gridNodes[firstNode.m_m].front() = firstNewNodeCoordinates;
            m_gridNodes[secondNode.m_m].front() = secondNewNodeCoordinates;
            return;
        }

        m_gridNodes[firstNode.m_m][firstNode.m_n - 1] = firstNewNodeCoordinates;
        m_gridNodes[secondNode.m_m][secondNode.m_n - 1] = secondNewNodeCoordinates;
    }

    if (gridLineType == BoundaryGridLineType::Up)
    {
        auto const firstNewNodeCoordinates = m_gridNodes[firstNode.m_m][firstNode.m_n] * 2.0 - m_gridNodes[firstNode.m_m][firstNode.m_n - 1];
        auto const secondNewNodeCoordinates = m_gridNodes[secondNode.m_m][secondNode.m_n] * 2.0 - m_gridNodes[secondNode.m_m][secondNode.m_n - 1];
        AddGridLineAtBoundary(firstNode, secondNode);
        m_gridNodes[firstNode.m_m][firstNode.m_n + 1] = firstNewNodeCoordinates;
        m_gridNodes[secondNode.m_m][secondNode.m_n + 1] = secondNewNodeCoordinates;
    }
}

std::tuple<double, double, double>
CurvilinearGrid::ComputeDirectionalSmoothingFactors(CurvilinearGridNodeIndices const& gridpoint,
                                                    const CurvilinearGridNodeIndices& pointOnSmoothingLineIndices,
                                                    const CurvilinearGridNodeIndices& lowerLeftIndices,
                                                    const CurvilinearGridNodeIndices& upperRightIndices)
{
    // horizontal smoothing factor
    const auto horizontalDelta = gridpoint.m_m > pointOnSmoothingLineIndices.m_m ? gridpoint.m_m - pointOnSmoothingLineIndices.m_m : pointOnSmoothingLineIndices.m_m - gridpoint.m_m;
    const auto maxHorizontalDelta = gridpoint.m_m > pointOnSmoothingLineIndices.m_m ? upperRightIndices.m_m - pointOnSmoothingLineIndices.m_m : pointOnSmoothingLineIndices.m_m - lowerLeftIndices.m_m;
    const auto mSmoothingFactor = maxHorizontalDelta == 0 ? 1.0 : (1.0 + std::cos(M_PI * static_cast<double>(horizontalDelta) / static_cast<double>(maxHorizontalDelta))) * 0.5;

    // vertical smoothing factor
    const auto verticalDelta = gridpoint.m_n > pointOnSmoothingLineIndices.m_n ? gridpoint.m_n - pointOnSmoothingLineIndices.m_n : pointOnSmoothingLineIndices.m_n - gridpoint.m_n;
    const auto maxVerticalDelta = gridpoint.m_n > pointOnSmoothingLineIndices.m_n ? upperRightIndices.m_n - pointOnSmoothingLineIndices.m_n : pointOnSmoothingLineIndices.m_n - lowerLeftIndices.m_n;
    const auto nSmoothingFactor = maxVerticalDelta == 0 ? 1.0 : (1.0 + std::cos(M_PI * static_cast<double>(verticalDelta) / static_cast<double>(maxVerticalDelta))) * 0.5;

    // mixed smoothing factor
    const auto mixedSmoothingFactor = std::sqrt(nSmoothingFactor * mSmoothingFactor);

    return {mSmoothingFactor, nSmoothingFactor, mixedSmoothingFactor};
}

CurvilinearGrid CurvilinearGrid::CloneCurvilinearGrid() const
{
    return CurvilinearGrid(m_gridNodes, m_projection);
}

double CurvilinearGrid::ComputeAverageNodalDistance(CurvilinearGridNodeIndices const& index, CurvilinearGridLine::GridLineDirection direction)
{
    if (index.m_m > m_gridNodes.size() || index.m_n > m_gridNodes[0].size())
    {
        throw std::invalid_argument("CurvilinearGrid::ComputeAverageNodalDistance: invalid index coordinates");
    }
    if (direction != CurvilinearGridLine::GridLineDirection::MDirection && direction != CurvilinearGridLine::GridLineDirection::NDirection)
    {
        throw std::invalid_argument("CurvilinearGrid::ComputeAverageNodalDistance: invalid direction");
    }

    if (direction == CurvilinearGridLine::GridLineDirection::MDirection)
    {
        int numEdges = 0.0;
        double leftDistance = 0.0;
        double rightDistance = 0.0;
        if (index.m_m > 0 && m_gridNodes[index.m_m - 1][index.m_n].IsValid())
        {
            leftDistance = ComputeDistance(m_gridNodes[index.m_m][index.m_n], m_gridNodes[index.m_m - 1][index.m_n], m_projection);
            numEdges += 1;
        }
        if (index.m_m + 1 < m_gridNodes.size() && m_gridNodes[index.m_m + 1][index.m_n].IsValid())
        {
            rightDistance = ComputeDistance(m_gridNodes[index.m_m][index.m_n], m_gridNodes[index.m_m + 1][index.m_n], m_projection);
            numEdges += 1;
        }
        return numEdges == 0 ? 0.0 : (leftDistance + rightDistance) / numEdges;
    }
    if (direction == CurvilinearGridLine::GridLineDirection::NDirection)
    {
        int numEdges = 0;
        double bottomDistance = 0.0;
        double upDistance = 0.0;
        if (index.m_n > 0 && m_gridNodes[index.m_m][index.m_n - 1].IsValid())
        {
            bottomDistance = ComputeDistance(m_gridNodes[index.m_m][index.m_n], m_gridNodes[index.m_m][index.m_n - 1], m_projection);
            numEdges += 1;
        }
        if (index.m_n + 1 < m_gridNodes[0].size() && m_gridNodes[index.m_m][index.m_n + 1].IsValid())
        {
            upDistance = ComputeDistance(m_gridNodes[index.m_m][index.m_n], m_gridNodes[index.m_m][index.m_n + 1], m_projection);
            numEdges += 1;
        }
        return numEdges == 0 ? 0.0 : (bottomDistance + upDistance) / numEdges;
    }

    throw std::invalid_argument("CurvilinearGrid::ComputeAverageNodalDistance: Invalid direction");
}

meshkernel::Point CurvilinearGrid::TransformDisplacement(Point const& displacement, CurvilinearGridNodeIndices const& node, bool isLocal) const
{
    Point left = m_gridNodes[node.m_m][node.m_n];
    Point right = left;
    if (node.m_m < m_numM - 1 && m_gridNodes[node.m_m + 1][node.m_n].IsValid())
    {
        right = m_gridNodes[node.m_m + 1][node.m_n];
    }
    if (node.m_m > 0 && m_gridNodes[node.m_m - 1][node.m_n].IsValid())
    {
        left = m_gridNodes[node.m_m - 1][node.m_n];
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
