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
#include <MeshKernel/CurvilinearGrid/CurvilinearGridNodeIndices.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridSmoothing.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridUtilities.hpp>
#include <MeshKernel/CurvilinearGrid/UndoActions/CurvilinearGridBlockUndoAction.hpp>
#include <MeshKernel/Operations.hpp>

using meshkernel::CurvilinearGrid;
using meshkernel::CurvilinearGridSmoothing;

CurvilinearGridSmoothing::CurvilinearGridSmoothing(CurvilinearGrid& grid, UInt smoothingIterations) : CurvilinearGridAlgorithm(grid),
                                                                                                      m_smoothingIterations(smoothingIterations)

{
    // Allocate cache for storing grid nodes values
    lin_alg::ResizeAndFillMatrix(m_gridNodesCache, m_grid.NumN(), m_grid.NumM(), true);

    // Compute the grid node types
    m_grid.ComputeGridNodeTypes();
}

meshkernel::UndoActionPtr CurvilinearGridSmoothing::Compute()
{
    CurvilinearGridNodeIndices upperLimit = m_upperRight;
    ++upperLimit.m_n;
    ++upperLimit.m_m;

    std::unique_ptr<CurvilinearGridBlockUndoAction> undoAction = CurvilinearGridBlockUndoAction::Create(m_grid, m_lowerLeft, upperLimit);

    // Compute the frozen nodes
    ComputeFrozenNodes();

    // Perform smoothing iterations
    for (UInt smoothingIterations = 0; smoothingIterations < m_smoothingIterations; ++smoothingIterations)
    {
        Solve();
    }

    return undoAction;
}

std::unique_ptr<CurvilinearGrid> CurvilinearGridSmoothing::ComputeDirectional(const Point& firstPoint, const Point& secondPoint)
{
    const auto gridLine = m_grid.GetGridLine(firstPoint, secondPoint);

    // Points are coinciding, this no smoothing zone
    if ((gridLine.IsMGridLine() && m_lowerLeft.m_n == m_upperRight.m_n) ||
        (gridLine.IsNGridLine() && m_lowerLeft.m_m == m_upperRight.m_m))
    {
        throw std::invalid_argument("CurvilinearGridSmoothing::Compute The points defining the smoothing area have the same direction of the smoothing line.");
    }

    // Re-compute the smoothing block: the block must coincide with the line end and start nodes
    CurvilinearGridNodeIndices lowerLeftBlock;
    CurvilinearGridNodeIndices upperRightBlock;
    if (gridLine.IsMGridLine())
    {
        lowerLeftBlock = {std::min(m_lowerLeft.m_n, m_upperRight.m_n), gridLine.m_startCoordinate};
        upperRightBlock = {std::max(m_lowerLeft.m_n, m_upperRight.m_n), gridLine.m_endCoordinate};
    }
    else
    {
        lowerLeftBlock = {gridLine.m_startCoordinate, std::min(m_lowerLeft.m_m, m_upperRight.m_m)};
        upperRightBlock = {gridLine.m_endCoordinate, std::max(m_lowerLeft.m_m, m_upperRight.m_m)};
    }
    m_lowerLeft = lowerLeftBlock;
    m_upperRight = upperRightBlock;

    // Compute the frozen nodes
    ComputeFrozenNodes();

    // Perform smoothing iterations
    for (UInt smoothingIterations = 0; smoothingIterations < m_smoothingIterations; ++smoothingIterations)
    {
        SolveDirectional(gridLine);
    }

    return std::make_unique<CurvilinearGrid>(m_grid);
}

std::tuple<meshkernel::Point, meshkernel::Point> CurvilinearGridSmoothing::ComputeGridDelta(const UInt n,
                                                                                            const UInt m,
                                                                                            const CurvilinearGridLine& gridLine) const
{
    Point firstDelta;
    Point secondDelta;

    if (gridLine.IsNGridLine())
    {
        firstDelta = m_gridNodesCache(n, m) - m_gridNodesCache(n - 1, m);
        secondDelta = m_gridNodesCache(n, m) - m_gridNodesCache(n + 1, m);
    }
    else
    {
        firstDelta = m_gridNodesCache(n, m) - m_gridNodesCache(n, m - 1);
        secondDelta = m_gridNodesCache(n, m) - m_gridNodesCache(n, m + 1);
    }

    return {firstDelta, secondDelta};
}

meshkernel::Point CurvilinearGridSmoothing::ComputeSmoothedGridNode(const UInt n,
                                                                    const UInt m,
                                                                    const CurvilinearGridLine& gridLine,
                                                                    const double firstLengthSquared,
                                                                    const double secondLengthSquared) const
{
    const double smoothingFactor = 0.5;

    const auto maxlength = std::max(firstLengthSquared, secondLengthSquared);
    const auto characteristicLength = std::abs(secondLengthSquared - firstLengthSquared) * 0.5;
    const auto [mSmoothing, nSmoothing, mixedSmoothing] = CurvilinearGrid::ComputeDirectionalSmoothingFactors({n, m}, gridLine.m_startNode, m_lowerLeft, m_upperRight);

    Point gridNode;

    if (gridLine.IsNGridLine())
    {
        // smooth along vertical
        const auto a = maxlength < 1e-8 ? 0.5 : mSmoothing * smoothingFactor * characteristicLength / maxlength;
        const auto maxDelta = firstLengthSquared > secondLengthSquared ? m_gridNodesCache(n - 1, m) - m_grid.GetNode(n, m) : m_gridNodesCache(n + 1, m) - m_grid.GetNode(n, m);
        gridNode = m_gridNodesCache(n, m) + maxDelta * a;
    }
    else
    {
        // smooth along horizontal
        const auto a = maxlength < 1e-8 ? 0.5 : nSmoothing * smoothingFactor * characteristicLength / maxlength;
        const auto maxDelta = firstLengthSquared > secondLengthSquared ? m_gridNodesCache(n, m - 1) - m_grid.GetNode(n, m) : m_gridNodesCache(n, m + 1) - m_grid.GetNode(n, m);
        gridNode = m_gridNodesCache(n, m) + maxDelta * a;
    }

    return gridNode;
}

void CurvilinearGridSmoothing::SolveDirectional(const CurvilinearGridLine& gridLine)
{

    // assign current nodal values to the m_gridNodesCache
    m_gridNodesCache = m_grid.GetNodes();

    auto isInvalidValidNode = [&](auto const& n, auto const& m)
    {
        if (gridLine.IsNGridLine())
        {
            return m_grid.GetNodeType(n, m) != NodeType::InternalValid &&
                   m_grid.GetNodeType(n, m) != NodeType::Bottom &&
                   m_grid.GetNodeType(n, m) != NodeType::Up;
        }

        return m_grid.GetNodeType(n, m) != NodeType::InternalValid &&
               m_grid.GetNodeType(n, m) != NodeType::Left &&
               m_grid.GetNodeType(n, m) != NodeType::Right;
    };

    // Apply smoothing
    for (auto n = m_lowerLeft.m_n; n <= m_upperRight.m_n; ++n)
    {
        for (auto m = m_lowerLeft.m_m; m <= m_upperRight.m_m; ++m)
        {
            // No update to perform for frozen nodes
            if (m_isGridNodeFrozen(n, m))
            {
                continue;
            }

            // Apply line smoothing only in internal nodes
            if (isInvalidValidNode(n, m))
            {
                continue;
            }

            // Calculate influence radius
            auto [firstDelta, secondDelta] = ComputeGridDelta(n, m, gridLine);

            const auto firstLengthSquared = lengthSquared(firstDelta);
            const auto secondLengthSquared = lengthSquared(secondDelta);

            m_grid.GetNode(n, m) = ComputeSmoothedGridNode(n, m, gridLine, firstLengthSquared, secondLengthSquared);
        }
    }
}

void CurvilinearGridSmoothing::Solve()
{
    double const a = 0.5;
    double const b = 1.0 - a;

    // assign current nodal values to the m_gridNodesCache
    m_gridNodesCache = m_grid.GetNodes();

    // Apply smoothing
    for (auto n = m_lowerLeft.m_n; n <= m_upperRight.m_n; ++n)
    {
        for (auto m = m_lowerLeft.m_m; m <= m_upperRight.m_m; ++m)
        {
            // No update to perform for frozen nodes
            if (m_isGridNodeFrozen(n, m))
            {
                continue;
            }

            // It is invalid or a corner point, skip smoothing
            if (m_grid.GetNodeType(n, m) == NodeType::Invalid ||
                m_grid.GetNodeType(n, m) == NodeType::BottomLeft ||
                m_grid.GetNodeType(n, m) == NodeType::UpperLeft ||
                m_grid.GetNodeType(n, m) == NodeType::BottomRight ||
                m_grid.GetNodeType(n, m) == NodeType::UpperRight)
            {
                continue;
            }

            // Compute new position based on a smoothing operator
            if (m_grid.GetNodeType(n, m) == NodeType::InternalValid)
            {
                const auto val = m_gridNodesCache(n, m) * a + (m_gridNodesCache(n - 1, m) + m_gridNodesCache(n + 1, m)) * 0.25 * b +
                                 (m_gridNodesCache(n, m - 1) + m_gridNodesCache(n, m + 1)) * 0.25 * b;

                m_grid.GetNode(n, m) = val;
                continue;
            }

            // For the point on the boundaries first computed the new position
            Point newNodePosition;
            if (m_grid.GetNodeType(n, m) == NodeType::Bottom)
            {
                newNodePosition = m_gridNodesCache(n, m) * a + (m_gridNodesCache(n - 1, m) + m_gridNodesCache(n + 1, m) + m_gridNodesCache(n, m + 1)) * constants::numeric::oneThird * b;
            }
            if (m_grid.GetNodeType(n, m) == NodeType::Up)
            {
                newNodePosition = m_gridNodesCache(n, m) * a + (m_gridNodesCache(n - 1, m) + m_gridNodesCache(n + 1, m) + m_gridNodesCache(n, m - 1)) * constants::numeric::oneThird * b;
            }
            if (m_grid.GetNodeType(n, m) == NodeType::Right)
            {
                newNodePosition = m_gridNodesCache(n, m) * a + (m_gridNodesCache(n, m - 1) + m_gridNodesCache(n, m + 1) + m_gridNodesCache(n - 1, m)) * constants::numeric::oneThird * b;
            }
            if (m_grid.GetNodeType(n, m) == NodeType::Left)
            {
                newNodePosition = m_gridNodesCache(n, m) * a + (m_gridNodesCache(n, m - 1) + m_gridNodesCache(n, m + 1) + m_gridNodesCache(n + 1, m)) * constants::numeric::oneThird * b;
            }

            ProjectPointOnClosestGridBoundary(newNodePosition, n, m);
        }
    }
}

void CurvilinearGridSmoothing::ProjectPointOnClosestGridBoundary(Point const& point, UInt n, UInt m)
{
    // Project the new position on the original boundary segment
    Point previousNode;
    Point nextNode;
    if (m_grid.GetNodeType(n, m) == NodeType::Bottom || m_grid.GetNodeType(n, m) == NodeType::Up)
    {
        previousNode = m_gridNodesCache(n - 1, m);
        nextNode = m_gridNodesCache(n + 1, m);
    }
    if (m_grid.GetNodeType(n, m) == NodeType::Right || m_grid.GetNodeType(n, m) == NodeType::Left)
    {
        previousNode = m_gridNodesCache(n, m - 1);
        nextNode = m_gridNodesCache(n, m + 1);
    }

    const auto [firstProjectedPoint, firstRatio, firstProjectedPointOnSegment] = OrthogonalProjectionOnSegment(m_gridNodesCache(n, m), previousNode, point);
    const auto [secondProjectedPoint, secondRatio, secondProjectedPointOnSegment] = OrthogonalProjectionOnSegment(m_gridNodesCache(n, m), nextNode, point);

    if (firstProjectedPointOnSegment && secondProjectedPointOnSegment && secondRatio > firstRatio)
    {
        m_grid.GetNode(n, m) = secondProjectedPoint;
        return;
    }
    if (firstProjectedPointOnSegment && secondProjectedPointOnSegment && secondRatio <= firstRatio)
    {
        m_grid.GetNode(n, m) = firstProjectedPoint;
        return;
    }
    if (firstProjectedPointOnSegment)
    {
        m_grid.GetNode(n, m) = firstProjectedPoint;
        return;
    }
    if (secondProjectedPointOnSegment)
    {
        m_grid.GetNode(n, m) = secondProjectedPoint;
        return;
    }

    m_grid.GetNode(n, m) = (firstProjectedPoint + secondProjectedPoint) * 0.5;
}
