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

#include <MeshKernel/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGridSmoothing.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Operations.hpp>

using meshkernel::CurvilinearGrid;
using meshkernel::CurvilinearGridSmoothing;

CurvilinearGridSmoothing::CurvilinearGridSmoothing(std::shared_ptr<CurvilinearGrid> grid, size_t smoothingIterations) : CurvilinearGridAlgorithm(grid), m_smoothingIterations(smoothingIterations)

{
    // Allocate cache for storing grid nodes values
    m_gridNodesCache.resize(m_grid->m_gridNodes.size(), std::vector<Point>(m_grid->m_gridNodes[0].size()));
    // Compute the grid node types
    m_grid->ComputeGridNodeTypes();
}

std::shared_ptr<CurvilinearGrid> CurvilinearGridSmoothing::Compute()
{
    // Perform smoothing iterations
    for (auto smoothingIterations = 0; smoothingIterations < m_smoothingIterations; ++smoothingIterations)
    {
        Solve();
    }
    return m_grid;
}

std::shared_ptr<CurvilinearGrid> CurvilinearGridSmoothing::Compute(Point const& lowerLeftCornerSmoothingArea, Point const& upperRightCornerSmootingArea)
{
    if (m_lines.empty())
    {
        throw std::invalid_argument("CurvilinearGridSmoothing::Compute No line set for directional refinement.");
    }

    auto const leftPointIndices = m_grid->GetNodeIndices(lowerLeftCornerSmoothingArea);
    auto const rightPointIndices = m_grid->GetNodeIndices(upperRightCornerSmootingArea);

    // Points are coinciding, this no smoothing zone
    if (m_lines[0].m_gridLineType == CurvilinearGrid::GridLine::GridLineType::MGridLine && leftPointIndices.m_n == rightPointIndices.m_n ||
        m_lines[0].m_gridLineType == CurvilinearGrid::GridLine::GridLineType::NGridLine && leftPointIndices.m_m == rightPointIndices.m_m)
    {
        throw std::invalid_argument("CurvilinearGridSmoothing::Compute The points defining the smoothing area have the same direction of the smoothing line.");
    }

    // Compute the smoothing area
    if (m_lines[0].m_gridLineType == CurvilinearGrid::GridLine::GridLineType::MGridLine)
    {
        m_lowerLeft = {m_lines[0].m_startCoordinate, std::min(leftPointIndices.m_n, rightPointIndices.m_n)};
        m_upperRight = {m_lines[0].m_endCoordinate, std::max(leftPointIndices.m_n, rightPointIndices.m_n)};
    }
    else
    {
        m_lowerLeft = {std::min(leftPointIndices.m_m, rightPointIndices.m_m), m_lines[0].m_startCoordinate};
        m_upperRight = {std::max(leftPointIndices.m_m, rightPointIndices.m_m), m_lines[0].m_endCoordinate};
    }

    // compute the box of the smoothing zone, used to determine the smoothing factors
    auto const [lowerLeft, upperRight] = m_grid->ComputeBlockFromCornerPoints(leftPointIndices, rightPointIndices);

    // Perform smoothing iterations
    for (auto smoothingIterations = 0; smoothingIterations < m_smoothingIterations; ++smoothingIterations)
    {
        Solve(lowerLeft, upperRight);
    }

    return m_grid;
}

void CurvilinearGridSmoothing::Solve(CurvilinearGrid::NodeIndices const& lowerLeftCornerRegion,
                                     CurvilinearGrid::NodeIndices const& upperRightCornerSmoothingRegion)
{

    // assign current nodal values to the m_gridNodesCache
    for (auto m = 0; m < m_grid->m_gridNodes.size(); ++m)
    {
        for (auto n = 0; n < m_grid->m_gridNodes[0].size(); ++n)
        {
            m_gridNodesCache[m][n] = m_grid->m_gridNodes[m][n];
        }
    }

    auto isInvalidValidNode = [this](auto const& m, auto const& n) {
        if (m_lines[0].m_gridLineType == CurvilinearGrid::GridLine::GridLineType::MGridLine)
        {
            return m_grid->m_gridNodesMask[m][n] != CurvilinearGrid::NodeType::InternalValid &&
                   m_grid->m_gridNodesMask[m][n] != CurvilinearGrid::NodeType::Bottom &&
                   m_grid->m_gridNodesMask[m][n] != CurvilinearGrid::NodeType::Up;
        }

        return m_grid->m_gridNodesMask[m][n] != CurvilinearGrid::NodeType::InternalValid &&
               m_grid->m_gridNodesMask[m][n] != CurvilinearGrid::NodeType::Left &&
               m_grid->m_gridNodesMask[m][n] != CurvilinearGrid::NodeType::Right;
    };

    // Apply smoothing
    const double smoothingFactor = 0.5;
    for (auto m = m_lowerLeft.m_m; m <= m_upperRight.m_m; ++m)
    {
        for (auto n = m_lowerLeft.m_n; n <= m_upperRight.m_n; ++n)
        {
            // Apply line smoothing only in internal nodes
            if (isInvalidValidNode(m, n))
            {
                continue;
            }

            // Calculate influence radius
            Point firstDelta;
            Point secondDelta;
            if (m_lines[0].m_gridLineType == CurvilinearGrid::GridLine::GridLineType::MGridLine)
            {
                firstDelta = m_gridNodesCache[m][n] - m_gridNodesCache[m - 1][n];
                secondDelta = m_gridNodesCache[m][n] - m_gridNodesCache[m + 1][n];
            }
            else
            {
                firstDelta = m_gridNodesCache[m][n] - m_gridNodesCache[m][n - 1];
                secondDelta = m_gridNodesCache[m][n] - m_gridNodesCache[m][n + 1];
            }

            const auto firstLengthSquared = firstDelta.x * firstDelta.x + firstDelta.y * firstDelta.y;
            const auto secondLengthSquared = secondDelta.x * secondDelta.x + secondDelta.y * secondDelta.y;
            const auto maxlength = std::max(firstLengthSquared, secondLengthSquared);
            const auto characteristicLength = std::abs(secondLengthSquared - firstLengthSquared) * 0.5;
            const auto [mSmoothing, nSmoothing, mixedSmoothing] = CurvilinearGrid::ComputeDirectionalSmoothingFactors({m, n}, m_lines[0].m_startNode, lowerLeftCornerRegion, upperRightCornerSmoothingRegion);

            if (m_lines[0].m_gridLineType == CurvilinearGrid::GridLine::GridLineType::MGridLine)
            {
                // smooth along vertical
                const auto a = maxlength < 1e-8 ? 0.5 : nSmoothing * smoothingFactor * characteristicLength / maxlength;
                const auto maxDelta = firstLengthSquared > secondLengthSquared ? m_gridNodesCache[m - 1][n] - m_grid->m_gridNodes[m][n] : m_gridNodesCache[m + 1][n] - m_grid->m_gridNodes[m][n];
                m_grid->m_gridNodes[m][n] = m_gridNodesCache[m][n] + maxDelta * a;
            }
            else
            {
                // smooth along horizontal
                const auto a = maxlength < 1e-8 ? 0.5 : mSmoothing * smoothingFactor * characteristicLength / maxlength;
                const auto maxDelta = firstLengthSquared > secondLengthSquared ? m_gridNodesCache[m][n - 1] - m_grid->m_gridNodes[m][n] : m_gridNodesCache[m][n + 1] - m_grid->m_gridNodes[m][n];
                m_grid->m_gridNodes[m][n] = m_gridNodesCache[m][n] + maxDelta * a;
            }
        }
    }
}

void CurvilinearGridSmoothing::Solve()
{
    double const a = 0.5;
    double const b = 1.0 - a;

    // assign current nodal values to the m_gridNodesCache
    for (auto m = 0; m < m_grid->m_gridNodes.size(); ++m)
    {
        for (auto n = 0; n < m_grid->m_gridNodes[0].size(); ++n)
        {
            m_gridNodesCache[m][n] = m_grid->m_gridNodes[m][n];
        }
    }

    // Apply smoothing
    for (auto m = m_lowerLeft.m_m; m <= m_upperRight.m_m; ++m)
    {
        for (auto n = m_lowerLeft.m_n; n <= m_upperRight.m_n; ++n)
        {

            // It is invalid or a corner point, skip smoothing
            if (m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeType::Invalid ||
                m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeType::BottomLeft ||
                m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeType::UpperLeft ||
                m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeType::BottomRight ||
                m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeType::UpperRight)
            {
                continue;
            }

            // Compute new position based on a smoothing operator
            if (m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeType::InternalValid)
            {
                m_grid->m_gridNodes[m][n] = m_gridNodesCache[m][n] * a + (m_gridNodesCache[m - 1][n] + m_gridNodesCache[m + 1][n]) * 0.25 * b +
                                            (m_gridNodesCache[m][n - 1] + m_gridNodesCache[m][n + 1]) * 0.25 * b;
                continue;
            }

            // For the point on the boundaries first computed the new position
            Point newNodePosition;
            if (m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeType::Bottom)
            {
                newNodePosition = m_gridNodesCache[m][n] * a + (m_gridNodesCache[m - 1][n] + m_gridNodesCache[m + 1][n] + m_gridNodesCache[m][n + 1]) * oneThird * b;
            }
            if (m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeType::Up)
            {
                newNodePosition = m_gridNodesCache[m][n] * a + (m_gridNodesCache[m - 1][n] + m_gridNodesCache[m + 1][n] + m_gridNodesCache[m][n - 1]) * oneThird * b;
            }
            if (m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeType::Right)
            {
                newNodePosition = m_gridNodesCache[m][n] * a + (m_gridNodesCache[m][n - 1] + m_gridNodesCache[m][n + 1] + m_gridNodesCache[m - 1][n]) * oneThird * b;
            }
            if (m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeType::Left)
            {
                newNodePosition = m_gridNodesCache[m][n] * a + (m_gridNodesCache[m][n - 1] + m_gridNodesCache[m][n + 1] + m_gridNodesCache[m + 1][n]) * oneThird * b;
            }

            ProjectPointOnClosestGridBoundary(newNodePosition, m, n);
        }
    }
}

void CurvilinearGridSmoothing::ProjectPointOnClosestGridBoundary(Point const& point, size_t m, size_t n)
{
    // Project the new position on the original boundary segment
    Point previousNode;
    Point nextNode;
    if (m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeType::Bottom || m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeType::Up)
    {
        previousNode = m_gridNodesCache[m - 1][n];
        nextNode = m_gridNodesCache[m + 1][n];
    }
    if (m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeType::Right || m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeType::Left)
    {
        previousNode = m_gridNodesCache[m][n - 1];
        nextNode = m_gridNodesCache[m][n + 1];
    }

    const auto [firstProjectedPoint, firstRatio, firstProjectedPointOnSegment] = OrthogonalProjectionOnSegment(m_gridNodesCache[m][n], previousNode, point);
    const auto [secondProjectedPoint, secondRatio, secondProjectedPointOnSegment] = OrthogonalProjectionOnSegment(m_gridNodesCache[m][n], nextNode, point);

    if (firstProjectedPointOnSegment && secondProjectedPointOnSegment && secondRatio > firstRatio)
    {
        m_grid->m_gridNodes[m][n] = secondProjectedPoint;
        return;
    }
    if (firstProjectedPointOnSegment && secondProjectedPointOnSegment && secondRatio <= firstRatio)
    {
        m_grid->m_gridNodes[m][n] = firstProjectedPoint;
        return;
    }
    if (firstProjectedPointOnSegment)
    {
        m_grid->m_gridNodes[m][n] = firstProjectedPoint;
        return;
    }
    if (secondProjectedPointOnSegment)
    {
        m_grid->m_gridNodes[m][n] = secondProjectedPoint;
        return;
    }

    m_grid->m_gridNodes[m][n] = (firstProjectedPoint + secondProjectedPoint) * 0.5;
}
