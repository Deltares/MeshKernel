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
#include <MeshKernel/CurvilinearGrid/ResetCurvilinearNodeAction.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/Utilities/NumericFunctions.hpp>

using meshkernel::CurvilinearGrid;
using meshkernel::CurvilinearGridNodeIndices;

CurvilinearGrid::CurvilinearGrid(const CurvilinearGrid& grid) : Mesh(grid.m_edges, grid.Nodes(), grid.m_projection),
                                                                m_gridNodes(grid.m_gridNodes),
                                                                m_gridFacesMask(grid.m_gridFacesMask),
                                                                m_gridNodesTypes(grid.m_gridNodesTypes),
                                                                m_gridIndices(grid.m_gridIndices)
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

void CurvilinearGrid::SetFlatCopies()
{
    if (IsEmpty())
    {
        return;
    }

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

    std::vector<Point> nodes(NumN() * NumM());
    std::vector<Edge> edges(NumM() * (NumN() - 1) +
                            (NumM() - 1) * NumN());
    lin_alg::Matrix<UInt> nodeIndices(NumN(), NumM());
    nodeIndices.setConstant(constants::missing::uintValue);
    std::vector gridIndices(nodes.size(),
                            CurvilinearGridNodeIndices{constants::missing::uintValue,
                                                       constants::missing::uintValue});

    UInt ind = 0;
    for (UInt n = 0; n < NumN(); n++)
    {
        for (UInt m = 0; m < NumM(); m++)
        {

            nodes[ind] = GetNode(n, m);
            nodeIndices(n, m) = ind;
            gridIndices[ind] = {n, m};
            ind++;
        }
    }

    ind = 0;
    for (UInt n = 0; n < NumN() - 1; n++)
    {
        for (UInt m = 0; m < NumM(); m++)
        {

            edges[ind].first = nodeIndices(n, m);
            edges[ind].second = nodeIndices(n + 1, m);
            ind++;
        }
    }

    for (UInt n = 0; n < NumN(); n++)
    {
        for (UInt m = 0; m < NumM() - 1; m++)
        {
            edges[ind].first = nodeIndices(n, m);
            edges[ind].second = nodeIndices(n, m + 1);
            ind++;
        }
    }
    edges.resize(ind);

    return {nodes, edges, gridIndices};
}

bool CurvilinearGrid::IsValid() const
{
    if (IsEmpty())
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
    lin_alg::ResizeAndFillMatrix(m_gridFacesMask, FullNumN() - 1, FullNumM() - 1, false, false);
    for (UInt n = 0; n < NumN() - 1; ++n)
    {
        for (UInt m = 0; m < NumM() - 1; ++m)
        {
            // Only if all grid nodes of the face are valid, the face is valid
            if (!AreFaceNodesValid(n, m))
            {
                continue;
            }
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

    invalidNodesToRemove = false;
    // Flag nodes not connected to valid faces
    for (UInt n = 1; n < NumN() - 1; ++n)
    {
        for (UInt m = 1; m < NumM() - 1; ++m)
        {
            if (GetNode(n, m).IsValid() &&
                !IsFaceMaskValid(n, m) &&
                !IsFaceMaskValid(n, m - 1) &&
                !IsFaceMaskValid(n - 1, m - 1) &&
                !IsFaceMaskValid(n - 1, m))
            {
                GetNode(n, m) = {constants::missing::doubleValue, constants::missing::doubleValue};
                invalidNodesToRemove = true;
            }
        }
    }

    for (UInt m = 1; m < NumM() - 1; ++m)
    {
        if (GetNode(0, m).IsValid() &&
            !IsFaceMaskValid(0, m - 1) &&
            !IsFaceMaskValid(0, m))
        {
            GetNode(0, m) = {constants::missing::doubleValue, constants::missing::doubleValue};
        }
    }

    for (UInt n = 1; n < NumN() - 1; ++n)
    {
        if (GetNode(n, 0).IsValid() &&
            !IsFaceMaskValid(n - 1, 0) &&
            !IsFaceMaskValid(n, 0))
        {
            GetNode(n, 0) = {constants::missing::doubleValue, constants::missing::doubleValue};
        }
    }

    if (GetNode(0, 0).IsValid() && !IsFaceMaskValid(0, 0))
    {
        GetNode(0, 0) = {constants::missing::doubleValue, constants::missing::doubleValue};
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
    SetFlatCopies();

    return undoAction;
}

void CurvilinearGrid::DeleteGridLineAtBoundary(CurvilinearGridNodeIndices const& firstNode, CurvilinearGridNodeIndices const& secondNode)
{
    bool const areNodesValid = GetNode(firstNode.m_n, firstNode.m_m).IsValid() &&
                               GetNode(secondNode.m_n, secondNode.m_m).IsValid();

    // Allocation depends on directions
    auto const gridLineType = GetBoundaryGridLineType(firstNode, secondNode);

    if (areNodesValid)
    {
        if (gridLineType == BoundaryGridLineType::Left)
        {
            // n-direction
            m_startOffset.m_n += 1;
            // need to reset the m_gridNodes, m_gridFacesMask, m_gridNodesTypes
            // for the deleted row
        }

        if (gridLineType == BoundaryGridLineType::Right)
        {
            // n-direction
            m_endOffset.m_n += 1;
            // need to reset the m_gridNodes, m_gridFacesMask, m_gridNodesTypes
            // for the deleted row
        }

        if (gridLineType == BoundaryGridLineType::Up)
        {
            // m-direction
            m_endOffset.m_m += 1;
            // need to reset the m_gridNodes, m_gridFacesMask, m_gridNodesTypes
            // for the deleted column
        }

        if (gridLineType == BoundaryGridLineType::Bottom)
        {
            // m-direction
            m_startOffset.m_m += 1;
            // need to reset the m_gridNodes, m_gridFacesMask, m_gridNodesTypes
            // for the deleted column
        }
    }
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
    bool const areNodesValid = GetNode(firstNode.m_n, firstNode.m_m).IsValid() &&
                               GetNode(secondNode.m_n, secondNode.m_m).IsValid();

    // Allocation depends on directions
    bool gridSizeChanged = false;
    auto const gridLineType = GetBoundaryGridLineType(firstNode, secondNode);

    UndoActionPtr undoAction;

    if (areNodesValid)
    {

        if (gridLineType == BoundaryGridLineType::Left)
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

        if (gridLineType == BoundaryGridLineType::Right)
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

        if (gridLineType == BoundaryGridLineType::Up)
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

        if (gridLineType == BoundaryGridLineType::Bottom)
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

meshkernel::UndoActionPtr CurvilinearGrid::AddEdge(CurvilinearGridNodeIndices const& firstNode, CurvilinearGridNodeIndices const& secondNode)
{

    // Allocate new grid line if needed
    auto const gridLineType = GetBoundaryGridLineType(firstNode, secondNode);
    // Allocation depends on directions

    std::unique_ptr<CompoundUndoAction> undoAddEdge = CompoundUndoAction::Create ();

    if (gridLineType == BoundaryGridLineType::Left)
    {
        auto const firstNewNodeCoordinates = GetNode(firstNode.m_n, firstNode.m_m) * 2.0 - GetNode(firstNode.m_n + 1, firstNode.m_m);
        auto const secondNewNodeCoordinates = GetNode(secondNode.m_n, secondNode.m_m) * 2.0 - GetNode(secondNode.m_n + 1, secondNode.m_m);
        auto [isGridLineAdded, addGridLineUndo] = AddGridLineAtBoundary(firstNode, secondNode);

        if (isGridLineAdded)
        {
            undoAddEdge->Add (std::move(addGridLineUndo));
            undoAddEdge->Add (MoveNode (CurvilinearGridNodeIndices(0, firstNode.m_m), firstNewNodeCoordinates));
            undoAddEdge->Add (MoveNode (CurvilinearGridNodeIndices(0, secondNode.m_m), secondNewNodeCoordinates));
        }
        else
        {
            undoAddEdge->Add (MoveNode (CurvilinearGridNodeIndices(firstNode.m_n - 1, firstNode.m_m), firstNewNodeCoordinates));
            undoAddEdge->Add (MoveNode (CurvilinearGridNodeIndices(secondNode.m_n - 1, secondNode.m_m), secondNewNodeCoordinates));
        }

        return undoAddEdge;
    }

    if (gridLineType == BoundaryGridLineType::Right)
    {
        auto const firstNewNodeCoordinates = GetNode(firstNode.m_n, firstNode.m_m) * 2.0 - GetNode(firstNode.m_n - 1, firstNode.m_m);
        auto const secondNewNodeCoordinates = GetNode(secondNode.m_n, secondNode.m_m) * 2.0 - GetNode(secondNode.m_n - 1, secondNode.m_m);
        auto [isGridLineAdded, addGridLineUndo] = AddGridLineAtBoundary(firstNode, secondNode);
        undoAddEdge->Add (std::move(addGridLineUndo));
        undoAddEdge->Add (MoveNode (CurvilinearGridNodeIndices(firstNode.m_n + 1, firstNode.m_m), firstNewNodeCoordinates));
        undoAddEdge->Add (MoveNode (CurvilinearGridNodeIndices(secondNode.m_n + 1, secondNode.m_m), secondNewNodeCoordinates));

        return undoAddEdge;
    }

    if (gridLineType == BoundaryGridLineType::Bottom)
    {
        auto const firstNewNodeCoordinates = GetNode(firstNode.m_n, firstNode.m_m) * 2.0 - GetNode(firstNode.m_n, firstNode.m_m + 1);
        auto const secondNewNodeCoordinates = GetNode(secondNode.m_n, secondNode.m_m) * 2.0 - GetNode(secondNode.m_n, secondNode.m_m + 1);
        auto [isGridLineAdded, addGridLineUndo] = AddGridLineAtBoundary(firstNode, secondNode);

        if (isGridLineAdded)
        {
            undoAddEdge->Add (std::move(addGridLineUndo));
            // Assign the new coordinates
            undoAddEdge->Add (MoveNode (CurvilinearGridNodeIndices(firstNode.m_n, 0), firstNewNodeCoordinates));
            undoAddEdge->Add (MoveNode (CurvilinearGridNodeIndices(secondNode.m_n, 0), secondNewNodeCoordinates));
        }
        else
        {
            GetNode(firstNode.m_n, firstNode.m_m - 1) = firstNewNodeCoordinates;
            GetNode(secondNode.m_n, secondNode.m_m - 1) = secondNewNodeCoordinates;
        }

        return undoAddEdge;
    }

    if (gridLineType == BoundaryGridLineType::Up)
    {
        auto const firstNewNodeCoordinates = GetNode(firstNode.m_n, firstNode.m_m) * 2.0 - GetNode(firstNode.m_n, firstNode.m_m - 1);
        auto const secondNewNodeCoordinates = GetNode(secondNode.m_n, secondNode.m_m) * 2.0 - GetNode(secondNode.m_n, secondNode.m_m - 1);
        [[maybe_unused]] auto [isGridLineAdded, addGridLineUndo] = AddGridLineAtBoundary(firstNode, secondNode);
        undoAddEdge->Add (std::move(addGridLineUndo));
        undoAddEdge->Add (MoveNode (CurvilinearGridNodeIndices(firstNode.m_n, firstNode.m_m + 1), firstNewNodeCoordinates));
        undoAddEdge->Add (MoveNode (CurvilinearGridNodeIndices(secondNode.m_n, secondNode.m_m + 1), secondNewNodeCoordinates));
    }

    return  undoAddEdge;
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
    auto const nodeToDelete = GetNodeIndices(point);

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
        SetFlatCopies();
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
    auto const nodeIndex = GetNodeIndices(fromPoint);
    return MoveNode (nodeIndex, toPoint);

    // // Check the node indices are valid
    // if (!nodeIndex.IsValid())
    // {
    //     throw std::invalid_argument("CurvilinearGrid::MoveNode node indices not found");
    // }

    // std::unique_ptr<ResetCurvilinearNodeAction> undoAction = ResetCurvilinearNodeAction::Create(*this, nodeIndex, GetNode(nodeIndex), toPoint,
    //                                                                                             false /* recomputeNodeTypes */);

    // // move fromPoint to toPoint
    // GetNode(nodeIndex.m_n, nodeIndex.m_m) = toPoint;
    // return undoAction;
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
        lowerLeft.x = std::min(lowerLeft.x, GetNode(0, i).x);
        lowerLeft.y = std::min(lowerLeft.y, GetNode(0, i).y);
        upperRight.x = std::max(upperRight.x, GetNode(0, i).x);
        upperRight.y = std::max(upperRight.y, GetNode(0, i).y);
    }

    // First loop over right boundary (last,i)
    for (Eigen::Index i = 0; i < NumN(); ++i)
    {
        lowerLeft.x = std::min(lowerLeft.x, GetNode(i, last).x);
        lowerLeft.y = std::min(lowerLeft.y, GetNode(i, last).y);
        upperRight.x = std::max(upperRight.x, GetNode(i, last).x);
        upperRight.y = std::max(upperRight.y, GetNode(i, last).y);
    }

    // This assumes that each column has the same number of points
    last = NumN() - 1;

    // First loop over upper boundary (i,last)
    for (Eigen::Index i = 0; i < NumM(); ++i)
    {
        lowerLeft.x = std::min(lowerLeft.x, GetNode(last, i).x);
        lowerLeft.y = std::min(lowerLeft.y, GetNode(last, i).y);
        upperRight.x = std::max(upperRight.x, GetNode(last, i).x);
        upperRight.y = std::max(upperRight.y, GetNode(last, i).y);
    }

    // First loop over left boundary (0,i)
    for (Eigen::Index i = 0; i < NumN(); ++i)
    {
        lowerLeft.x = std::min(lowerLeft.x, GetNode(i, 0).x);
        lowerLeft.y = std::min(lowerLeft.y, GetNode(i, 0).y);
        upperRight.x = std::max(upperRight.x, GetNode(i, 0).x);
        upperRight.y = std::max(upperRight.y, GetNode(i, 0).y);
    }

    return BoundingBox(lowerLeft, upperRight);
}

void CurvilinearGrid::print(std::ostream& out)
{
    out << std::endl;
    out << "grid dim: " << NumN() << " x " << NumM() << std::endl;

    for (int row = NumN() - 1; row >= 0; --row)
    {
        for (UInt col = 0; col < NumM(); ++col)
        {
            const Point p = GetNode(static_cast<UInt>(row), col);
            out << "{" << p.x << ", " << p.y << "}  ";
        }

        out << std::endl;
    }

    out << std::endl;
}

void CurvilinearGrid::printGraph(std::ostream& out)
{
    Print(m_nodes, m_edges, out);
}

void CurvilinearGrid::RestoreAction(AddGridLineUndoAction& undoAction)
{
    m_startOffset += undoAction.StartOffset();
    m_endOffset += undoAction.EndOffset();
}

void CurvilinearGrid::CommitAction(AddGridLineUndoAction& undoAction)
{
    m_startOffset -= undoAction.StartOffset();
    m_endOffset -= undoAction.EndOffset();
}

void CurvilinearGrid::RestoreAction(CurvilinearGridBlockUndo& undoAction)
{
    undoAction.Swap(*this);
}

void CurvilinearGrid::CommitAction(CurvilinearGridBlockUndo& undoAction)
{
    undoAction.Swap(*this);
}

void CurvilinearGrid::RestoreAction(CurvilinearGridRefinementUndoAction& undoAction)
{
    undoAction.Swap(m_gridNodes, m_startOffset, m_endOffset);
    ComputeGridNodeTypes();
}

void CurvilinearGrid::CommitAction(CurvilinearGridRefinementUndoAction& undoAction)
{
    undoAction.Swap(m_gridNodes, m_startOffset, m_endOffset);
    ComputeGridNodeTypes();
}

void CurvilinearGrid::RestoreAction(const ResetCurvilinearNodeAction& undoAction)
{
    GetNode(undoAction.NodeId()) = undoAction.InitialNode();

    if (undoAction.RecalculateNodeTypes())
    {
        ComputeGridNodeTypes();
        SetFlatCopies();
    }
}

void CurvilinearGrid::CommitAction(const ResetCurvilinearNodeAction& undoAction)
{
    GetNode(undoAction.NodeId()) = undoAction.UpdatedNode();

    if (undoAction.RecalculateNodeTypes())
    {
        ComputeGridNodeTypes();
        SetFlatCopies();
    }
}
