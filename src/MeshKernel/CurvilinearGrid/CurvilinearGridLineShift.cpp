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
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLineShift.hpp>
#include <MeshKernel/Entities.hpp>

using meshkernel::CurvilinearGrid;
using meshkernel::CurvilinearGridLineShift;
using meshkernel::Point;

CurvilinearGridLineShift::CurvilinearGridLineShift(std::shared_ptr<CurvilinearGrid> grid) : CurvilinearGridAlgorithm(grid)

{
    // store a deep copy of the grid for computing the displacements
    m_originalGrid = m_grid.CloneCurvilinearGrid();
}

CurvilinearGrid CurvilinearGridLineShift::Compute()
{
    if (m_lines.empty())
    {
        throw std::invalid_argument("CurvilinearGridLineShift::Compute No candidate line to shift has been selected");
    }

    /// The first delta
    auto const previousNodeIndex = m_lines[0].m_startNode;
    auto previousDelta = m_grid.m_gridNodes[previousNodeIndex.m_m][previousNodeIndex.m_n] -
                         m_originalGrid.m_gridNodes[previousNodeIndex.m_m][previousNodeIndex.m_n];

    const double eps = 1e-5;
    auto previousCoordinate = m_lines[0].m_startCoordinate;
    for (auto i = 1; i <= m_lines[0].m_endCoordinate; ++i)
    {
        auto const currentNodeIndex = m_lines[0].GetNodeIndexFromCoordinate(i);

        auto const currentDelta = m_grid.m_gridNodes[currentNodeIndex.m_m][currentNodeIndex.m_n] -
                                  m_originalGrid.m_gridNodes[currentNodeIndex.m_m][currentNodeIndex.m_n];

        if (std::abs(currentDelta.x) < eps && std::abs(currentDelta.y) < eps && i != m_lines[0].m_endCoordinate)
        {
            continue;
        }

        /// On the original algorithm currentDelta is distributed on the nodes above the current i,
        /// except for the last node m_endCoordinate, where currentDelta is distributed on the entire grid line
        const auto currentLastCoordinate = i == m_lines[0].m_endCoordinate ? i : i - 1;
        for (auto j = previousCoordinate; j <= currentLastCoordinate; ++j)
        {

            auto const nodeIndex = m_lines[0].GetNodeIndexFromCoordinate(j);

            auto const firstFactor = static_cast<double>(j - previousCoordinate) / static_cast<double>(i - previousCoordinate);
            auto const secondFactor = 1.0 - firstFactor;

            // Now distribute the shifting
            m_grid.m_gridNodes[nodeIndex.m_m][nodeIndex.m_n] = m_originalGrid.m_gridNodes[nodeIndex.m_m][nodeIndex.m_n] +
                                                               previousDelta * secondFactor + currentDelta * firstFactor;
            // Field transformation on the influence area
            TransformGrid(nodeIndex);
        }
        previousCoordinate = i;
        previousDelta = currentDelta;
    }

    return m_grid;
}

void CurvilinearGridLineShift::TransformGrid(CurvilinearGrid::NodeIndices const& node)
{
    auto delta = m_grid.m_gridNodes[node.m_m][node.m_n] - m_originalGrid.m_gridNodes[node.m_m][node.m_n];
    delta = m_originalGrid.TransformDisplacement(delta, node, true);

    auto const start = m_lines[0].m_gridLineType == CurvilinearGrid::GridLineDirection::MDirection ? m_lowerLeft.m_n : m_lowerLeft.m_m;
    auto const end = m_lines[0].m_gridLineType == CurvilinearGrid::GridLineDirection::MDirection ? m_upperRight.m_n : m_upperRight.m_m;

    for (auto i = start; i <= end; ++i)
    {
        CurvilinearGrid::NodeIndices currentNode{m_lines[0].m_gridLineType == CurvilinearGrid::GridLineDirection::MDirection ? node.m_m : i,
                                                 m_lines[0].m_gridLineType == CurvilinearGrid::GridLineDirection::MDirection ? i : node.m_n};

        if (!m_originalGrid.m_gridNodes[currentNode.m_m][currentNode.m_n].IsValid())
        {
            continue;
        }
        const auto [mSmoothing, nSmoothing, mixedSmoothing] = CurvilinearGrid::ComputeDirectionalSmoothingFactors(currentNode, m_lines[0].m_startNode, m_lowerLeft, m_upperRight);
        Point currentDelta{0.0, 0.0};
        if (m_lines[0].m_gridLineType == CurvilinearGrid::GridLineDirection::MDirection)
        {
            currentDelta = delta * nSmoothing;
        }

        if (m_lines[0].m_gridLineType == CurvilinearGrid::GridLineDirection::NDirection)
        {
            currentDelta = delta * mSmoothing;
        }

        currentDelta = m_originalGrid.TransformDisplacement(currentDelta, currentNode, false);
        m_grid.m_gridNodes[currentNode.m_m][currentNode.m_n] = m_originalGrid.m_gridNodes[currentNode.m_m][currentNode.m_n] + currentDelta;
    }
}

void CurvilinearGridLineShift::MoveNode(Point const& fromPoint, Point const& toPoint)
{
    if (m_lines.empty())
    {
        throw std::invalid_argument("CurvilinearGridLineShift::MoveNode No candidate line to shift has been selected");
    }

    if (!m_lowerLeft.IsValid() && !m_upperRight.IsValid())
    {
        throw std::invalid_argument("CurvilinearGridLineShift::MoveNode No block for smoothing the line shift has been selected");
    }

    // Get the index in the grid of the line to be shifted
    auto const nodeIndex = m_grid.GetNodeIndices(fromPoint);

    //Check the nodes are on the line to shift
    if (!m_lines[0].IsNodeOnLine(nodeIndex))
    {
        throw std::invalid_argument("CurvilinearGridLineShift::MoveNode The selected node does not belong to the line to be shifted");
    }

    m_grid.m_gridNodes[nodeIndex.m_m][nodeIndex.m_n] = toPoint;
}