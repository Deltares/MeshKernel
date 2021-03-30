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
#include <iostream>

meshkernel::CurvilinearGridSmoothing::CurvilinearGridSmoothing(std::shared_ptr<CurvilinearGrid> grid,
                                                               size_t smoothingIterations)
    : m_grid(grid),
      m_smoothingIterations(smoothingIterations)

{
    // Allocate cache for storing grid nodes values
    m_gridNodesCache.resize(m_grid->m_gridNodes.size(), std::vector<Point>(m_grid->m_gridNodes[0].size()));
    // Compute the grid node types
    m_grid->ComputeGridNodeTypes();
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
    auto const [lowerLeft, upperRight] = m_grid->ComputeBoundingBoxCornerPoints(firstNode, secondNode);
    m_minM = lowerLeft.m;
    m_maxM = upperRight.m;
    m_minN = lowerLeft.n;
    m_maxN = upperRight.n;
}

void meshkernel::CurvilinearGridSmoothing::Compute()
{
    // Perform smoothing iterations
    for (auto smoothingIterations = 0; smoothingIterations < m_smoothingIterations; ++smoothingIterations)
    {
        Solve();
    }
}

void meshkernel::CurvilinearGridSmoothing::ComputedDirectionalSmooth(Point const& firstLinePoint,
                                                                     Point const& secondLinePoint,
                                                                     Point const& lowerLeftCornerSmoothingArea,
                                                                     Point const& upperRightCornerSmootingArea)
{
    // Get the m and n indices from the point coordinates
    auto const firstLinePointIndices = m_grid->GetNodeIndices(firstLinePoint);
    auto const secondLinePointIndices = m_grid->GetNodeIndices(secondLinePoint);
    auto const leftPointIndices = m_grid->GetNodeIndices(lowerLeftCornerSmoothingArea);
    auto const rightPointIndices = m_grid->GetNodeIndices(upperRightCornerSmootingArea);

    // Determine if smoothing should occur in the m direction
    bool const isSmoothingAlongM = firstLinePointIndices.n == secondLinePointIndices.n;

    // Points are coinciding, this no smoothing zone
    if (leftPointIndices == rightPointIndices)
    {
        throw std::invalid_argument("CurvilinearGridSmoothing::ComputedDirectionalSmooth The points defining the smoothing area coincides.");
    }
    if (isSmoothingAlongM && leftPointIndices.n == rightPointIndices.n || !isSmoothingAlongM && leftPointIndices.m == rightPointIndices.m)
    {
        throw std::invalid_argument("CurvilinearGridSmoothing::ComputedDirectionalSmooth The points defining the smoothing area have the same direction of the smoothing line.");
    }

    // Compute the smoothing area
    if (isSmoothingAlongM)
    {
        m_minM = std::min(firstLinePointIndices.m, secondLinePointIndices.m);
        m_maxM = std::max(firstLinePointIndices.m, secondLinePointIndices.m);
        m_minN = std::min(leftPointIndices.n, rightPointIndices.n);
        m_maxN = std::max(leftPointIndices.n, rightPointIndices.n);
    }
    else
    {
        m_minM = std::min(leftPointIndices.m, rightPointIndices.m);
        m_maxM = std::max(leftPointIndices.m, rightPointIndices.m);
        m_minN = std::min(firstLinePointIndices.n, secondLinePointIndices.n);
        m_maxN = std::max(firstLinePointIndices.n, secondLinePointIndices.n);
    }

    // compute the box of the smoothing zone, used to determine the smoothing factors
    auto const [lowerLeft, upperRight] = m_grid->ComputeBoundingBoxCornerPoints(leftPointIndices, rightPointIndices);

    // Perform smoothing iterations
    for (auto smoothingIterations = 0; smoothingIterations < m_smoothingIterations; ++smoothingIterations)
    {
        SolveDirectionalSmooth(isSmoothingAlongM, firstLinePointIndices, lowerLeft, upperRight);
    }
}

void meshkernel::CurvilinearGridSmoothing::SolveDirectionalSmooth(bool isSmoothingAlongM,
                                                                  CurvilinearGrid::NodeIndices const& pointOnLineIndices,
                                                                  CurvilinearGrid::NodeIndices const& lowerLeftCornerRegion,
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

    auto isInvalidValidNode = [&isSmoothingAlongM, this](auto const& m, auto const& n) {
        if (isSmoothingAlongM)
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
    for (auto m = m_minM; m <= m_maxM; ++m)
    {
        for (auto n = m_minN; n <= m_maxN; ++n)
        {
            // Apply line smoothing only in internal nodes
            if (isInvalidValidNode(m, n))
            {
                continue;
            }

            // Calculate influence radius
            Point firstDelta;
            Point secondDelta;
            if (isSmoothingAlongM)
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
            const auto [mSmoothing, nSmoothing, mixedSmoothing] = ComputeDirectionalSmoothingFactors({m, n}, pointOnLineIndices, lowerLeftCornerRegion, upperRightCornerSmoothingRegion);

            if (isSmoothingAlongM)
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

std::tuple<double, double, double> meshkernel::CurvilinearGridSmoothing::ComputeDirectionalSmoothingFactors(CurvilinearGrid::NodeIndices const& gridpoint,
                                                                                                            const CurvilinearGrid::NodeIndices& pointOnSmoothingLineIndices,
                                                                                                            const CurvilinearGrid::NodeIndices& lowerLeftIndices,
                                                                                                            const CurvilinearGrid::NodeIndices& upperRightIndices) const
{

    // horizontal smoothing factor
    const auto horizontalDelta = gridpoint.m > pointOnSmoothingLineIndices.m ? gridpoint.m - pointOnSmoothingLineIndices.m : pointOnSmoothingLineIndices.m - gridpoint.m;
    const auto maxHorizontalDelta = gridpoint.m > pointOnSmoothingLineIndices.m ? upperRightIndices.m - pointOnSmoothingLineIndices.m : pointOnSmoothingLineIndices.m - lowerLeftIndices.m;
    const auto horizontalSmoothingFactor = (1.0 + std::cos(M_PI * static_cast<double>(horizontalDelta) / static_cast<double>(maxHorizontalDelta))) * 0.5;

    // vertical smoothing factor
    const auto verticalDelta = gridpoint.n > pointOnSmoothingLineIndices.n ? gridpoint.n - pointOnSmoothingLineIndices.n : pointOnSmoothingLineIndices.n - gridpoint.n;
    const auto maxVerticalDelta = gridpoint.n > pointOnSmoothingLineIndices.n ? upperRightIndices.n - pointOnSmoothingLineIndices.n : pointOnSmoothingLineIndices.n - lowerLeftIndices.n;
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
