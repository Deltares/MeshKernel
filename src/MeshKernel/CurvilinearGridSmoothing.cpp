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

meshkernel::CurvilinearGridSmoothing::CurvilinearGridSmoothing(std::shared_ptr<CurvilinearGrid> grid,
                                                               size_t smoothingIterations,
                                                               const Point& firstCornerPoint,
                                                               const Point& secondCornerPoint)
    : m_grid(grid),
      m_smoothingIterations(smoothingIterations),
      m_firstCornerPoint(firstCornerPoint),
      m_secondCornerPoint(secondCornerPoint)

{
}

void meshkernel::CurvilinearGridSmoothing::Compute()
{
    // Get the m and n indices from the point coordinates
    auto [mFirstNode, nFirstNode] = m_grid->GetNodeIndices(m_firstCornerPoint);
    auto [mSecondNode, nSecondNode] = m_grid->GetNodeIndices(m_secondCornerPoint);

    // Coinciding corner nodes, no valid area, nothing to do
    if (mFirstNode == mSecondNode && nFirstNode == nSecondNode)
    {
        return;
    }

    // Compute orthogonalization bounding box
    m_minM = std::min(mFirstNode, mSecondNode);
    m_minN = std::min(nFirstNode, nSecondNode);
    m_maxM = std::max(mFirstNode, mSecondNode);
    m_maxN = std::max(nFirstNode, nSecondNode);

    // Compute the grid node types
    m_grid->ComputeGridNodeTypes();

    // Allocate cache for storing grid nodes values
    m_gridNodesCache.resize(m_grid->m_gridNodes.size(), std::vector<Point>(m_grid->m_gridNodes[0].size()));

    // Compute the matrix coefficients
    for (auto smoothingIterations = 0; smoothingIterations < m_smoothingIterations; ++smoothingIterations)
    {
        Solve();
    }
}

void meshkernel::CurvilinearGridSmoothing::Solve()
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
    for (auto m = m_minM; m <= m_maxM; ++m)
    {
        for (auto n = m_minN; n <= m_maxN; ++n)
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
            Point p;
            if (m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeType::Bottom)
            {
                p = m_gridNodesCache[m][n] * a + (m_gridNodesCache[m - 1][n] + m_gridNodesCache[m + 1][n] + m_gridNodesCache[m][n + 1]) * oneThird * b;
            }
            if (m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeType::Up)
            {
                p = m_gridNodesCache[m][n] * a + (m_gridNodesCache[m - 1][n] + m_gridNodesCache[m + 1][n] + m_gridNodesCache[m][n - 1]) * oneThird * b;
            }
            if (m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeType::Right)
            {
                p = m_gridNodesCache[m][n] * a + (m_gridNodesCache[m][n - 1] + m_gridNodesCache[m][n + 1] + m_gridNodesCache[m - 1][n]) * oneThird * b;
            }
            if (m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeType::Left)
            {
                p = m_gridNodesCache[m][n] * a + (m_gridNodesCache[m][n - 1] + m_gridNodesCache[m][n + 1] + m_gridNodesCache[m + 1][n]) * oneThird * b;
            }

            ProjectPointOnClosestGridBoundary(p, m, n);
        }
    }
}

void meshkernel::CurvilinearGridSmoothing::ProjectPointOnClosestGridBoundary(Point const& point, size_t m, size_t n)
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
