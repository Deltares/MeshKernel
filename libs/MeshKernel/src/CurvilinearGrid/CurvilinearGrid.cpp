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

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLine.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Polygons.hpp>

using meshkernel::CurvilinearGrid;
using meshkernel::CurvilinearGridNodeIndices;

CurvilinearGrid::CurvilinearGrid(const CurvilinearGrid& grid) : Mesh(grid.m_edges, grid.m_nodes, grid.m_projection),
                                                                m_gridNodes(grid.m_gridNodes),
                                                                m_gridFacesMask(grid.m_gridFacesMask),
                                                                m_gridNodesTypes(grid.m_gridNodesTypes),
                                                                m_gridIndices(grid.m_gridIndices),
                                                                m_numM(grid.NumM()),
                                                                m_numN(grid.NumN())
{
}

CurvilinearGrid::CurvilinearGrid(Projection projection) : Mesh(projection) {}

CurvilinearGrid::CurvilinearGrid(lin_alg::Matrix<Point> const& grid, Projection projection) : Mesh(projection)
{
    SetGridNodes(grid);
}

void CurvilinearGrid::SetGridNodes(const lin_alg::Matrix<Point>& gridNodes)
{
    m_gridNodes = gridNodes;

    if (!IsValid())
    {
        throw std::invalid_argument("CurvilinearGrid::CurvilinearGrid: Invalid curvilinear grid");
    }

    m_numM = static_cast<UInt>(m_gridNodes.rows());
    m_numN = static_cast<UInt>(m_gridNodes.cols());

    m_nodesRTreeRequiresUpdate = true;
    m_edgesRTreeRequiresUpdate = true;
    m_facesRTreeRequiresUpdate = true;

    SetFlatCopies();
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

    const auto numN = static_cast<UInt>(m_gridNodes.rows());
    const auto numM = static_cast<UInt>(m_gridNodes.cols());

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

void CurvilinearGrid::SetFlatCopies()
{
    if (lin_alg::MatrixIsEmpty(m_gridNodes))
    {
        return;
    }

    m_numM = static_cast<UInt>(m_gridNodes.rows());
    m_numN = static_cast<UInt>(m_gridNodes.cols());
    const auto [nodes, edges, gridIndices] = ConvertCurvilinearToNodesAndEdges();
    m_nodes = nodes;
    m_edges = edges;
    m_gridIndices = gridIndices;
}

std::tuple<std::vector<meshkernel::Point>,
           std::vector<meshkernel::Edge>,
           std::vector<CurvilinearGridNodeIndices>>
CurvilinearGrid::ConvertCurvilinearToNodesAndEdges() const
{
    if (!IsValid())
    {
        throw std::invalid_argument("CurvilinearGrid::ConvertCurvilinearToNodesAndEdges: Invalid curvilinear grid ");
    }

    std::vector<Point> nodes(m_gridNodes.rows() * m_gridNodes.cols());
    std::vector<Edge> edges(m_gridNodes.rows() * (m_gridNodes.cols() - 1) +
                            (m_gridNodes.rows() - 1) * m_gridNodes.cols());
    lin_alg::Matrix<UInt> nodeIndices(m_gridNodes.rows(), m_gridNodes.cols());
    nodeIndices.setConstant(constants::missing::uintValue);
    std::vector<CurvilinearGridNodeIndices> gridIndices(nodes.size(),
                                                        CurvilinearGridNodeIndices{constants::missing::uintValue,
                                                                                   constants::missing::uintValue});

    UInt ind = 0;
    for (UInt m = 0; m < m_gridNodes.rows(); m++)
    {
        for (UInt n = 0; n < m_gridNodes.cols(); n++)
        {
            nodes[ind] = m_gridNodes(m, n);
            nodeIndices(m, n) = ind;
            gridIndices[ind] = {m, n};
            ind++;
        }
    }

    ind = 0;
    for (UInt m = 0; m < m_gridNodes.rows() - 1; m++)
    {
        for (UInt n = 0; n < m_gridNodes.cols(); n++)
        {
            if (nodeIndices(m, n) != constants::missing::uintValue &&
                nodeIndices(m + 1, n) != constants::missing::uintValue)
            {
                edges[ind].first = nodeIndices(m, n);
                edges[ind].second = nodeIndices(m + 1, n);
                ind++;
            }
        }
    }

    for (UInt m = 0; m < m_gridNodes.rows(); m++)
    {
        for (UInt n = 0; n < m_gridNodes.cols() - 1; n++)
        {
            if (nodeIndices(m, n) != constants::missing::uintValue &&
                nodeIndices(m, n + 1) != constants::missing::uintValue)
            {
                edges[ind].first = nodeIndices(m, n);
                edges[ind].second = nodeIndices(m, n + 1);
                ind++;
            }
        }
    }
    edges.resize(ind);

    return {nodes, edges, gridIndices};
}

bool CurvilinearGrid::IsValid() const
{
    if (m_gridNodes.size() == 0)
    {
        return false;
    }
    if (m_gridNodes.rows() < 2)
    {
        return false;
    }
    if (m_gridNodes.cols() < 2)
    {
        return false;
    }

    return true;
}

CurvilinearGridNodeIndices CurvilinearGrid::GetNodeIndices(Point point)
{
    BuildTree(Location::Nodes);
    SearchNearestLocation(point, Location::Nodes);
    if (GetNumLocations(Location::Nodes) == 0)
    {
        return {constants::missing::uintValue, constants::missing::uintValue};
    }

    const auto nodeIndex = GetLocationsIndices(0, Location::Nodes);
    return m_gridIndices[nodeIndex];
}

std::tuple<CurvilinearGridNodeIndices, CurvilinearGridNodeIndices> CurvilinearGrid::GetEdgeNodeIndices(Point const& point)
{
    BuildTree(Location::Edges);
    SearchNearestLocation(point, Location::Edges);
    if (GetNumLocations(Location::Edges) == 0)
    {
        return {{}, {}};
    }

    const auto nodeIndex = GetLocationsIndices(0, Location::Edges);
    auto const firstNode = m_edges[nodeIndex].first;
    auto const secondNode = m_edges[nodeIndex].second;

    return {m_gridIndices[firstNode], m_gridIndices[secondNode]};
}

bool CurvilinearGrid::AreFaceNodesValid(UInt m, UInt n) const
{
    return m_gridNodes(m, n).IsValid() &&
           m_gridNodes(m + 1, n).IsValid() &&
           m_gridNodes(m, n + 1).IsValid() &&
           m_gridNodes(m + 1, n + 1).IsValid();
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
    CurvilinearGridNodeIndices lowerLeft(std::min(firstNode.m_m, secondNode.m_m), std::min(firstNode.m_n, secondNode.m_n));
    CurvilinearGridNodeIndices upperRight(std::max(firstNode.m_m, secondNode.m_m), std::max(firstNode.m_n, secondNode.m_n));

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
    lin_alg::ResizeAndFillMatrix(m_gridFacesMask, NumM() - 1, NumN() - 1, false, false);
    for (UInt m = 0; m < NumM() - 1; ++m)
    {
        for (UInt n = 0; n < NumN() - 1; ++n)
        {
            // Only if all grid nodes of the face are valid, the face is valid
            if (!AreFaceNodesValid(m, n))
            {
                continue;
            }
            m_gridFacesMask(m, n) = true;
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
    for (UInt m = 1; m < NumM() - 1; ++m)
    {
        for (UInt n = 1; n < NumN() - 1; ++n)
        {
            if (m_gridNodes(m, n).IsValid() &&
                !m_gridFacesMask(m, n) &&
                !m_gridFacesMask(m - 1, n) &&
                !m_gridFacesMask(m - 1, n - 1) &&
                !m_gridFacesMask(m, n - 1))
            {
                m_gridNodes(m, n) = {constants::missing::doubleValue, constants::missing::doubleValue};
                invalidNodesToRemove = true;
            }
        }
    }

    for (UInt m = 1; m < NumM() - 1; ++m)
    {
        if (m_gridNodes(m, 0).IsValid() &&
            !m_gridFacesMask(m - 1, 0) &&
            !m_gridFacesMask(m, 0))
        {
            m_gridNodes(m, 0) = {constants::missing::doubleValue, constants::missing::doubleValue};
        }
    }

    for (UInt n = 1; n < NumN() - 1; ++n)
    {
        if (m_gridNodes(0, n).IsValid() &&
            !m_gridFacesMask(0, n - 1) &&
            !m_gridFacesMask(0, n))
        {
            m_gridNodes(0, n) = {constants::missing::doubleValue, constants::missing::doubleValue};
        }
    }

    if (m_gridNodes(0, 0).IsValid() && !m_gridFacesMask(0, 0))
    {
        m_gridNodes(0, 0) = {constants::missing::doubleValue, constants::missing::doubleValue};
    }

    RemoveInvalidNodes(invalidNodesToRemove);
}

void CurvilinearGrid::ComputeGridNodeTypes()
{
    RemoveInvalidNodes(true);
    lin_alg::ResizeAndFillMatrix(m_gridNodesTypes, NumM(), NumN(), false, NodeType::Invalid);

    // Flag faces based on boundaries
    for (UInt m = 0; m < NumM(); ++m)
    {
        for (UInt n = 0; n < NumN(); ++n)
        {

            if (!m_gridNodes(m, n).IsValid())
            {
                continue;
            }

            // Left side
            if (m == 0 && n == 0)
            {
                m_gridNodesTypes(m, n) = NodeType::BottomLeft;
                continue;
            }
            if (m == 0 && n == NumN() - 1)
            {
                m_gridNodesTypes(m, n) = NodeType::UpperLeft;
                continue;
            }
            if (m == 0 && !m_gridNodes(m, n - 1).IsValid())
            {
                m_gridNodesTypes(m, n) = NodeType::BottomLeft;
                continue;
            }
            if (m == 0 && !m_gridNodes(m, n + 1).IsValid())
            {
                m_gridNodesTypes(m, n) = NodeType::UpperLeft;
                continue;
            }
            if (m == 0)
            {
                m_gridNodesTypes(m, n) = NodeType::Left;
                continue;
            }
            // Right side
            if (m == NumM() - 1 && n == 0)
            {
                m_gridNodesTypes(m, n) = NodeType::BottomRight;
                continue;
            }
            if (m == NumM() - 1 && n == NumN() - 1)
            {
                m_gridNodesTypes(m, n) = NodeType::UpperRight;
                continue;
            }
            if (m == NumM() - 1 && !m_gridNodes(m, n - 1).IsValid())
            {
                m_gridNodesTypes(m, n) = NodeType::BottomRight;
                continue;
            }
            if (m == NumM() - 1 && !m_gridNodes(m, n + 1).IsValid())
            {
                m_gridNodesTypes(m, n) = NodeType::UpperRight;
                continue;
            }
            if (m == NumM() - 1)
            {
                m_gridNodesTypes(m, n) = NodeType::Right;
                continue;
            }
            // Bottom side
            if (n == 0 && !m_gridNodes(m - 1, n).IsValid())
            {
                m_gridNodesTypes(m, n) = NodeType::BottomLeft;
                continue;
            }
            if (n == 0 && !m_gridNodes(m + 1, n).IsValid())
            {
                m_gridNodesTypes(m, n) = NodeType::BottomRight;
                continue;
            }
            if (n == 0)
            {
                m_gridNodesTypes(m, n) = NodeType::Bottom;
                continue;
            }
            // Upper side
            if (n == NumN() - 1 && !m_gridNodes(m - 1, n).IsValid())
            {
                m_gridNodesTypes(m, n) = NodeType::UpperLeft;
                continue;
            }
            if (n == NumN() - 1 && !m_gridNodes(m + 1, n).IsValid())
            {
                m_gridNodesTypes(m, n) = NodeType::UpperRight;
                continue;
            }
            if (n == NumN() - 1)
            {
                m_gridNodesTypes(m, n) = NodeType::Up;
                continue;
            }

            auto const isTopLeftFaceValid = m_gridFacesMask(m - 1, n);
            auto const isTopRightFaceValid = m_gridFacesMask(m, n);
            auto const isBottomLeftFaceValid = m_gridFacesMask(m - 1, n - 1);
            auto const isBottomRightFaceValid = m_gridFacesMask(m, n - 1);

            if (isTopRightFaceValid &&
                isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                m_gridNodesTypes(m, n) = NodeType::InternalValid;
                continue;
            }
            if (!isTopRightFaceValid &&
                isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                m_gridNodesTypes(m, n) = NodeType::BottomLeft;
                continue;
            }
            if (isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                m_gridNodesTypes(m, n) = NodeType::BottomRight;
                continue;
            }
            if (isTopRightFaceValid &&
                isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                m_gridNodesTypes(m, n) = NodeType::UpperRight;
                continue;
            }
            if (isTopRightFaceValid &&
                isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                m_gridNodesTypes(m, n) = NodeType::UpperLeft;
                continue;
            }

            if (isTopRightFaceValid &&
                isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                m_gridNodesTypes(m, n) = NodeType::Bottom;
                continue;
            }
            if (isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                m_gridNodesTypes(m, n) = NodeType::Left;
                continue;
            }

            if (!isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                m_gridNodesTypes(m, n) = NodeType::Up;
                continue;
            }

            if (!isTopRightFaceValid &&
                isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                m_gridNodesTypes(m, n) = NodeType::Right;
                continue;
            }

            if (isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                m_gridNodesTypes(m, n) = NodeType::BottomLeft;
                continue;
            }

            if (!isTopRightFaceValid &&
                isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                m_gridNodesTypes(m, n) = NodeType::BottomRight;
                continue;
            }

            if (!isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                m_gridNodesTypes(m, n) = NodeType::UpperRight;
                continue;
            }

            if (!isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                m_gridNodesTypes(m, n) = NodeType::UpperLeft;
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

bool CurvilinearGrid::AddGridLineAtBoundary(CurvilinearGridNodeIndices const& firstNode,
                                            CurvilinearGridNodeIndices const& secondNode)
{
    // If both nodes are invalid, we can substitute the invalid values. New allocation is not needed.
    bool const areNodesValid = m_gridNodes(firstNode.m_m, firstNode.m_n).IsValid() &&
                               m_gridNodes(secondNode.m_m, secondNode.m_n).IsValid();

    // Allocation depends on directions
    bool gridSizeChanged = false;
    auto const gridLineType = GetBoundaryGridLineType(firstNode, secondNode);

    if (areNodesValid)
    {

        if (gridLineType == BoundaryGridLineType::Left)
        {
            lin_alg::InsertRow(m_gridNodes,
                               lin_alg::RowVector<Point>(m_gridNodes.cols()),
                               0);
            gridSizeChanged = true;
        }
        if (gridLineType == BoundaryGridLineType::Right)
        {
            lin_alg::InsertRow(m_gridNodes,
                               lin_alg::RowVector<Point>(m_gridNodes.cols()),
                               m_gridNodes.rows());
            gridSizeChanged = true;
        }
        if (gridLineType == BoundaryGridLineType::Up)
        {
            lin_alg::InsertCol(m_gridNodes,
                               lin_alg::ColVector<Point>(m_gridNodes.rows()),
                               m_gridNodes.cols());
            gridSizeChanged = true;
        }
        if (gridLineType == BoundaryGridLineType::Bottom)
        {
            lin_alg::InsertCol(m_gridNodes,
                               lin_alg::ColVector<Point>(m_gridNodes.rows()),
                               0);
            gridSizeChanged = true;
        }
    }

    return gridSizeChanged;
}

CurvilinearGrid::BoundaryGridLineType CurvilinearGrid::GetBoundaryGridLineType(CurvilinearGridNodeIndices const& firstNode, CurvilinearGridNodeIndices const& secondNode) const
{
    auto const firstNodeType = m_gridNodesTypes(firstNode.m_m, firstNode.m_n);
    auto const secondNodeType = m_gridNodesTypes(secondNode.m_m, secondNode.m_n);

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
        auto const firstNewNodeCoordinates = m_gridNodes(firstNode.m_m, firstNode.m_n) * 2.0 - m_gridNodes(firstNode.m_m + 1, firstNode.m_n);
        auto const secondNewNodeCoordinates = m_gridNodes(secondNode.m_m, secondNode.m_n) * 2.0 - m_gridNodes(secondNode.m_m + 1, secondNode.m_n);
        auto const isGridLineAdded = AddGridLineAtBoundary(firstNode, secondNode);
        if (isGridLineAdded)
        {
            m_gridNodes(0, firstNode.m_n) = firstNewNodeCoordinates;
            m_gridNodes(0, secondNode.m_n) = secondNewNodeCoordinates;
            return;
        }

        m_gridNodes(firstNode.m_m - 1, firstNode.m_n) = firstNewNodeCoordinates;
        m_gridNodes(secondNode.m_m - 1, secondNode.m_n) = secondNewNodeCoordinates;
        return;
    }
    if (gridLineType == BoundaryGridLineType::Right)
    {
        auto const firstNewNodeCoordinates = m_gridNodes(firstNode.m_m, firstNode.m_n) * 2.0 - m_gridNodes(firstNode.m_m - 1, firstNode.m_n);
        auto const secondNewNodeCoordinates = m_gridNodes(secondNode.m_m, secondNode.m_n) * 2.0 - m_gridNodes(secondNode.m_m - 1, secondNode.m_n);
        AddGridLineAtBoundary(firstNode, secondNode);
        m_gridNodes(firstNode.m_m + 1, firstNode.m_n) = firstNewNodeCoordinates;
        m_gridNodes(secondNode.m_m + 1, secondNode.m_n) = secondNewNodeCoordinates;

        return;
    }
    if (gridLineType == BoundaryGridLineType::Bottom)
    {
        auto const firstNewNodeCoordinates = m_gridNodes(firstNode.m_m, firstNode.m_n) * 2.0 - m_gridNodes(firstNode.m_m, firstNode.m_n + 1);
        auto const secondNewNodeCoordinates = m_gridNodes(secondNode.m_m, secondNode.m_n) * 2.0 - m_gridNodes(secondNode.m_m, secondNode.m_n + 1);
        auto const isGridLineAdded = AddGridLineAtBoundary(firstNode, secondNode);
        if (isGridLineAdded)
        {
            // Assign the new coordinates
            m_gridNodes(firstNode.m_m, 0) = firstNewNodeCoordinates;
            m_gridNodes(secondNode.m_m, 0) = secondNewNodeCoordinates;
            return;
        }

        m_gridNodes(firstNode.m_m, firstNode.m_n - 1) = firstNewNodeCoordinates;
        m_gridNodes(secondNode.m_m, secondNode.m_n - 1) = secondNewNodeCoordinates;
    }

    if (gridLineType == BoundaryGridLineType::Up)
    {
        auto const firstNewNodeCoordinates = m_gridNodes(firstNode.m_m, firstNode.m_n) * 2.0 - m_gridNodes(firstNode.m_m, firstNode.m_n - 1);
        auto const secondNewNodeCoordinates = m_gridNodes(secondNode.m_m, secondNode.m_n) * 2.0 - m_gridNodes(secondNode.m_m, secondNode.m_n - 1);
        AddGridLineAtBoundary(firstNode, secondNode);
        m_gridNodes(firstNode.m_m, firstNode.m_n + 1) = firstNewNodeCoordinates;
        m_gridNodes(secondNode.m_m, secondNode.m_n + 1) = secondNewNodeCoordinates;
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
    if (index.m_m > m_gridNodes.rows() || index.m_n > m_gridNodes.cols())
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
        if (index.m_m > 0 && m_gridNodes(index.m_m - 1, index.m_n).IsValid())
        {
            leftDistance = ComputeDistance(m_gridNodes(index.m_m, index.m_n), m_gridNodes(index.m_m - 1, index.m_n), m_projection);
            numEdges += 1;
        }
        if (index.m_m + 1 < m_gridNodes.rows() && m_gridNodes(index.m_m + 1, index.m_n).IsValid())
        {
            rightDistance = ComputeDistance(m_gridNodes(index.m_m, index.m_n), m_gridNodes(index.m_m + 1, index.m_n), m_projection);
            numEdges += 1;
        }
        return numEdges == 0 ? 0.0 : (leftDistance + rightDistance) / numEdges;
    }
    if (direction == CurvilinearGridLine::GridLineDirection::NDirection)
    {
        double numEdges = 0.0;
        double bottomDistance = 0.0;
        double upDistance = 0.0;
        if (index.m_n > 0 && m_gridNodes(index.m_m, index.m_n - 1).IsValid())
        {
            bottomDistance = ComputeDistance(m_gridNodes(index.m_m, index.m_n), m_gridNodes(index.m_m, index.m_n - 1), m_projection);
            numEdges += 1;
        }
        if (index.m_n + 1 < m_gridNodes.cols() && m_gridNodes(index.m_m, index.m_n + 1).IsValid())
        {
            upDistance = ComputeDistance(m_gridNodes(index.m_m, index.m_n), m_gridNodes(index.m_m, index.m_n + 1), m_projection);
            numEdges += 1;
        }
        return numEdges == 0 ? 0.0 : (bottomDistance + upDistance) / numEdges;
    }

    throw std::invalid_argument("CurvilinearGrid::ComputeAverageNodalDistance: Invalid direction");
}

meshkernel::Point CurvilinearGrid::TransformDisplacement(Point const& displacement, CurvilinearGridNodeIndices const& node, bool isLocal) const
{
    Point left = m_gridNodes(node.m_m, node.m_n);
    Point right = left;
    if (node.m_m < NumM() - 1 && m_gridNodes(node.m_m + 1, node.m_n).IsValid())
    {
        right = m_gridNodes(node.m_m + 1, node.m_n);
    }
    if (node.m_m > 0 && m_gridNodes(node.m_m - 1, node.m_n).IsValid())
    {
        left = m_gridNodes(node.m_m - 1, node.m_n);
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
        m_gridNodes(nodeToDelete.m_m, nodeToDelete.m_n) = {constants::missing::doubleValue, constants::missing::doubleValue};
        // Re-compute quantities
        ComputeGridNodeTypes();
        SetFlatCopies();
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
    m_gridNodes(nodeIndex.m_m, nodeIndex.m_n) = toPoint;
}

meshkernel::BoundingBox CurvilinearGrid::GetBoundingBox() const
{

    Point lowerLeft(std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
    Point upperRight(-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());
    size_t last = m_gridNodes.rows() - 1;

    // Only need to loop over boundary nodes

    // First loop over lower boundary (i,0)
    for (Eigen::Index i = 0; i < m_gridNodes.rows(); ++i)
    {
        lowerLeft.x = std::min(lowerLeft.x, m_gridNodes(i, 0).x);
        lowerLeft.y = std::min(lowerLeft.y, m_gridNodes(i, 0).y);
        upperRight.x = std::max(upperRight.x, m_gridNodes(i, 0).x);
        upperRight.y = std::max(upperRight.y, m_gridNodes(i, 0).y);
    }

    // First loop over right boundary (last,i)
    for (Eigen::Index i = 0; i < m_gridNodes.cols(); ++i)
    {
        lowerLeft.x = std::min(lowerLeft.x, m_gridNodes(last, i).x);
        lowerLeft.y = std::min(lowerLeft.y, m_gridNodes(last, i).y);
        upperRight.x = std::max(upperRight.x, m_gridNodes(last, i).x);
        upperRight.y = std::max(upperRight.y, m_gridNodes(last, i).y);
    }

    // This assumes that each column has the same number of points
    last = m_gridNodes.cols() - 1;

    // First loop over upper boundary (i,last)
    for (Eigen::Index i = 0; i < m_gridNodes.rows(); ++i)
    {
        lowerLeft.x = std::min(lowerLeft.x, m_gridNodes(i, last).x);
        lowerLeft.y = std::min(lowerLeft.y, m_gridNodes(i, last).y);
        upperRight.x = std::max(upperRight.x, m_gridNodes(i, last).x);
        upperRight.y = std::max(upperRight.y, m_gridNodes(i, last).y);
    }

    // First loop over left boundary (0,i)
    for (Eigen::Index i = 0; i < m_gridNodes.cols(); ++i)
    {
        lowerLeft.x = std::min(lowerLeft.x, m_gridNodes(0, i).x);
        lowerLeft.y = std::min(lowerLeft.y, m_gridNodes(0, i).y);
        upperRight.x = std::max(upperRight.x, m_gridNodes(0, i).x);
        upperRight.y = std::max(upperRight.y, m_gridNodes(0, i).y);
    }

    return BoundingBox(lowerLeft, upperRight);
}
