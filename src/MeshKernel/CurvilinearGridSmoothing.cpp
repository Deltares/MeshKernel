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
                                                               size_t smoothingIterations)
    : m_grid(grid),
      m_smoothingIterations(smoothingIterations)

{
}

void meshkernel::CurvilinearGridSmoothing::SetBlock(Point const& firstCornerPoint, Point const& secondCornerPoint)
{
    // Get the m and n indices from the point coordinates
    auto const firstNode = m_grid->GetNodeIndices(firstCornerPoint);
    auto const secondNode = m_grid->GetNodeIndices(secondCornerPoint);

    // Coinciding corner nodes, no valid area, nothing to do
    if (firstNode == secondNode)
    {
        return;
    }

    // Compute orthogonalization bounding box
    auto const [upperLeft, lowerRight] = m_grid->ComputeBoundingBoxCornerPoints(firstNode, secondNode);
    m_minM = upperLeft.m;
    m_maxM = lowerRight.m;
    m_minN = lowerRight.n;
    m_maxN = upperLeft.n;
}

void meshkernel::CurvilinearGridSmoothing::Compute()
{
    // Compute the grid node types
    m_grid->ComputeGridNodeTypes();

    // Allocate cache for storing grid nodes values
    m_gridNodesCache.resize(m_grid->m_gridNodes.size(), std::vector<Point>(m_grid->m_gridNodes[0].size()));

    // Perform smoothing iterations
    for (auto smoothingIterations = 0; smoothingIterations < m_smoothingIterations; ++smoothingIterations)
    {
        Solve();
    }
}

void meshkernel::CurvilinearGridSmoothing::ComputeLineSmooth(Point const& firstLinePoint,
                                                             Point const& secondLinePoint,
                                                             Point const& leftPointInfluenceZone,
                                                             Point const& rightPointInfluenceZone)
{
    // Get the m and n indices from the point coordinates
    auto const firstLinePointIndices = m_grid->GetNodeIndices(firstLinePoint);
    auto const secondLinePointIndices = m_grid->GetNodeIndices(secondLinePoint);
    auto const leftPointIndices = m_grid->GetNodeIndices(leftPointInfluenceZone);
    auto const rightPointIndices = m_grid->GetNodeIndices(rightPointInfluenceZone);
    // Points are coinciding, this no influence area
    if (leftPointIndices == rightPointIndices)
    {
        return;
    }

    // Determine if an m-grid line is used
    bool const isMGridLine = firstLinePointIndices.n == secondLinePointIndices.n;

    auto const [upperLeft, lowerRight] = m_grid->ComputeBoundingBoxCornerPoints(leftPointIndices, rightPointIndices);
    const double smoothingFactor = 0.2;

    // Perform smoothing iterations
    for (auto smoothingIterations = 0; smoothingIterations < m_smoothingIterations; ++smoothingIterations)
    {

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
                // Apply line smoothing only in internal nodes
                if (m_grid->m_gridNodesMask[m][n] != CurvilinearGrid::NodeType::InternalValid)
                {
                    continue;
                }

                // Calculate influence radius
                Point firstDelta;
                Point secondDelta;
                if (isMGridLine)
                {
                    firstDelta = m_grid->m_gridNodes[m][n] - m_grid->m_gridNodes[m - 1][n];
                    secondDelta = m_grid->m_gridNodes[m][n] - m_grid->m_gridNodes[m + 1][n];
                }
                else
                {
                    firstDelta = m_grid->m_gridNodes[m][n] - m_grid->m_gridNodes[m][n - 1];
                    secondDelta = m_grid->m_gridNodes[m][n] - m_grid->m_gridNodes[m][n + 1];
                }

                const auto firstLengthSquared = firstDelta.x * firstDelta.x + firstDelta.y * firstDelta.y;
                const auto secondLengthSquared = secondDelta.x * secondDelta.x + secondDelta.y * secondDelta.y;
                const auto maxlength = std::max(firstLengthSquared, secondLengthSquared);
                const auto carateristicLength = std::abs(secondLengthSquared - firstLengthSquared) * 0.5;
                const auto [mSmoothing, nSmoothing, mixedSmoothing] = ComputeSmoothingFactors({m, n}, firstLinePointIndices, secondLinePointIndices, upperLeft, lowerRight);
                const auto a = maxlength < 1e-8 ? 0.5 : mSmoothing * smoothingFactor * carateristicLength / maxlength;

                // smooth
                if (isMGridLine)
                {
                    m_grid->m_gridNodes[m][n] = m_gridNodesCache[m][n] + (m_gridNodesCache[m + 1][n] - m_gridNodesCache[m][n]) * a;
                }
                else
                {
                    m_grid->m_gridNodes[m][n] = m_gridNodesCache[m][n] + (m_gridNodesCache[m][n + 1] - m_gridNodesCache[m][n]) * a;
                }
            }
        }
    }
}

std::tuple<double, double, double> meshkernel::CurvilinearGridSmoothing::ComputeSmoothingFactors(CurvilinearGrid::NodeIndices const& gridpoint,
                                                                                                 const CurvilinearGrid::NodeIndices& firstLinePointIndices,
                                                                                                 const CurvilinearGrid::NodeIndices& secondLinePointIndices,
                                                                                                 const CurvilinearGrid::NodeIndices& leftPointInfluenceZone,
                                                                                                 const CurvilinearGrid::NodeIndices& rightPointInfluenceZone) const
{

    // horizontal smoothing factor
    const auto horizontalDelta = gridpoint.m > firstLinePointIndices.m ? gridpoint.m - firstLinePointIndices.m : firstLinePointIndices.m - gridpoint.m;
    const auto maxHorizontalDelta = gridpoint.m > firstLinePointIndices.m ? rightPointInfluenceZone.m - firstLinePointIndices.m : firstLinePointIndices.m - leftPointInfluenceZone.m;
    const auto horizontalSmoothingFactor = (1.0 + std::cos(M_PI * static_cast<double>(horizontalDelta) / static_cast<double>(maxHorizontalDelta))) * 0.5;

    // vertical smoothing factor
    const auto verticalDelta = gridpoint.n > firstLinePointIndices.n ? gridpoint.n - firstLinePointIndices.n : firstLinePointIndices.n - gridpoint.n;
    const auto maxVerticalDelta = gridpoint.n > firstLinePointIndices.n ? rightPointInfluenceZone.n - firstLinePointIndices.n : firstLinePointIndices.n - leftPointInfluenceZone.n;
    const auto verticalSmoothingFactor = (1.0 + std::cos(M_PI * static_cast<double>(verticalDelta) / static_cast<double>(maxVerticalDelta))) * 0.5;

    // mixed smoothing factor
    const auto mixedSmoothingFactor = std::sqrt(verticalSmoothingFactor * horizontalSmoothingFactor);

    return {horizontalSmoothingFactor, verticalSmoothingFactor, mixedSmoothingFactor};
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
