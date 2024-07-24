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
#include <MeshKernel/UndoActions/CompoundUndoAction.hpp>

using meshkernel::CurvilinearGrid;
using meshkernel::CurvilinearGridLineShift;
using meshkernel::CurvilinearGridNodeIndices;
using meshkernel::Point;

CurvilinearGridLineShift::CurvilinearGridLineShift(CurvilinearGrid& grid) : CurvilinearGridAlgorithm(grid),
                                                                            m_originalGrid(grid)

{
}

meshkernel::UndoActionPtr CurvilinearGridLineShift::Compute()
{
    if (m_lines.empty())
    {
        throw std::invalid_argument("CurvilinearGridLineShift::Compute No candidate line to shift has been selected");
    }

    auto const startN = m_lines[0].IsNGridLine() ? m_lines[0].m_startCoordinate : m_lowerLeft.m_n;
    auto const endN = m_lines[0].IsNGridLine() ? m_lines[0].m_endCoordinate : m_upperRight.m_n;
    auto const startM = m_lines[0].IsMGridLine() ? m_lines[0].m_startCoordinate : m_lowerLeft.m_m;
    auto const endM = m_lines[0].IsMGridLine() ? m_lines[0].m_endCoordinate : m_upperRight.m_m;

    std::unique_ptr<CurvilinearGridBlockUndoAction> undoAction = CurvilinearGridBlockUndoAction::Create(m_grid, {startN, startM}, {endN + 1, endM + 1});

    /// The first delta
    auto const previousNodeIndex = m_lines[0].m_startNode;
    auto previousDelta = m_grid.GetNode(previousNodeIndex.m_n, previousNodeIndex.m_m) -
                         m_originalGrid.GetNode(previousNodeIndex.m_n, previousNodeIndex.m_m);

    const double eps = 1e-5;
    auto previousCoordinate = m_lines[0].m_startCoordinate;
    for (UInt i = 1; i <= m_lines[0].m_endCoordinate; ++i)
    {
        auto const currentNodeIndex = m_lines[0].GetNodeIndexFromCoordinate(i);

        auto const currentDelta = m_grid.GetNode(currentNodeIndex.m_n, currentNodeIndex.m_m) -
                                  m_originalGrid.GetNode(currentNodeIndex.m_n, currentNodeIndex.m_m);

        if (std::abs(currentDelta.x) < eps && std::abs(currentDelta.y) < eps && i != m_lines[0].m_endCoordinate)
        {
            continue;
        }

        /// On the original algorithm currentDelta is distributed on the nodes above the current i,
        /// except for the last node m_endCoordinate, where currentDelta is distributed on the entire grid line
        const auto currentLastCoordinate = i == m_lines[0].m_endCoordinate ? i : i - 1;
        for (UInt j = previousCoordinate; j <= currentLastCoordinate; ++j)
        {

            auto const nodeIndex = m_lines[0].GetNodeIndexFromCoordinate(j);

            auto const firstFactor = static_cast<double>(j - previousCoordinate) / static_cast<double>(i - previousCoordinate);
            auto const secondFactor = 1.0 - firstFactor;

            // Now distribute the shifting
            const auto newPosition = m_originalGrid.GetNode(nodeIndex.m_n, nodeIndex.m_m) +
                                     previousDelta * secondFactor + currentDelta * firstFactor;
            m_grid.GetNode(nodeIndex.m_n, nodeIndex.m_m) = newPosition;
            // Field transformation on the influence area
            TransformGrid(nodeIndex);
        }
        previousCoordinate = i;
        previousDelta = currentDelta;
    }

    return undoAction;
}

meshkernel::UndoActionPtr CurvilinearGridLineShift::Compute(const Point& fromPoint, const Point& toPoint)
{
    std::unique_ptr<CompoundUndoAction> compoundAction = CompoundUndoAction::Create();

    compoundAction->Add(MoveNode(fromPoint, toPoint));
    compoundAction->Add(Compute());

    return compoundAction;
}

void CurvilinearGridLineShift::TransformGrid(CurvilinearGridNodeIndices const& node)
{
    auto delta = m_grid.GetNode(node.m_n, node.m_m) - m_originalGrid.GetNode(node.m_n, node.m_m);
    delta = m_originalGrid.TransformDisplacement(delta, node, true);

    auto const start = m_lines[0].IsMGridLine() ? m_lowerLeft.m_n : m_lowerLeft.m_m;
    auto const end = m_lines[0].IsMGridLine() ? m_upperRight.m_n : m_upperRight.m_m;

    for (auto i = start; i <= end; ++i)
    {
        CurvilinearGridNodeIndices currentNode{
            m_lines[0].IsMGridLine() ? i : node.m_n,
            m_lines[0].IsMGridLine() ? node.m_m : i,
        };

        if (!m_originalGrid.GetNode(currentNode.m_n, currentNode.m_m).IsValid())
        {
            continue;
        }
        const auto [mSmoothing, nSmoothing, mixedSmoothing] = CurvilinearGrid::ComputeDirectionalSmoothingFactors(currentNode, m_lines[0].m_startNode, m_lowerLeft, m_upperRight);
        Point currentDelta{0.0, 0.0};
        if (m_lines[0].IsMGridLine())
        {
            currentDelta = delta * nSmoothing;
        }

        if (m_lines[0].IsNGridLine())
        {
            currentDelta = delta * mSmoothing;
        }

        currentDelta = m_originalGrid.TransformDisplacement(currentDelta, currentNode, false);
        m_grid.GetNode(currentNode.m_n, currentNode.m_m) = m_originalGrid.GetNode(currentNode.m_n, currentNode.m_m) + currentDelta;
    }
}

meshkernel::UndoActionPtr CurvilinearGridLineShift::MoveNode(Point const& fromPoint, Point const& toPoint)
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

    auto const nodePosition = m_grid.FindLocationIndex(fromPoint, Location::Nodes);
    auto const nodeIndex = m_grid.GetCurvilinearGridNodeIndices(nodePosition);

    // Check the nodes are on the line to shift
    if (!m_lines[0].IsNodeOnLine(nodeIndex))
    {
        throw std::invalid_argument("CurvilinearGridLineShift::MoveNode The selected node does not belong to the line to be shifted");
    }

    return m_grid.MoveNode(fromPoint, toPoint);
}
