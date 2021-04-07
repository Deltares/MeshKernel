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
#include <MeshKernel/CurvilinearGridLineShift.hpp>
#include <MeshKernel/Entities.hpp>

using meshkernel::CurvilinearGrid;
using meshkernel::CurvilinearGridLineShift;

CurvilinearGridLineShift::CurvilinearGridLineShift(std::shared_ptr<CurvilinearGrid> grid) : CurvilinearGridAlgorithm(grid)

{
    // Compute the grid node types
    m_grid->ComputeGridNodeTypes();
    m_gridCache = std::make_shared<CurvilinearGrid>(m_grid->m_gridNodes, m_grid->m_projection);
}

std::shared_ptr<CurvilinearGrid> CurvilinearGridLineShift::Compute()
{
    if (m_lines.empty())
    {
        throw std::invalid_argument("CurvilinearGridLineShift::Compute No candidate line to shift has been selected");
    }

    auto firstLineCoordinate = m_lines[0].m_startCoordinate;
    auto firstNodeOnTheLine = m_lines[0].m_startNode;
    auto firstDelta = m_gridCache->m_gridNodes[firstNodeOnTheLine.m_m][firstNodeOnTheLine.m_n] -
                      m_grid->m_gridNodes[firstNodeOnTheLine.m_m][firstNodeOnTheLine.m_n];

    const double eps = 1e-5;
    for (auto i = 1; i <= m_lines[0].m_endCoordinate; ++i)
    {
        auto const externalNodeIndex = m_lines[0].GetNodeindexFromCoordinate(i);

        auto secondDelta = m_grid->m_gridNodes[externalNodeIndex.m_m][externalNodeIndex.m_n] -
                           m_gridCache->m_gridNodes[externalNodeIndex.m_m][externalNodeIndex.m_n];

        if (std::abs(secondDelta.x) > eps || std::abs(secondDelta.y) > eps || i == m_lines[0].m_endCoordinate)
        {

            for (auto j = firstLineCoordinate; j <= i; ++j)
            {

                auto const nodeIndex = m_lines[0].GetNodeindexFromCoordinate(j);

                double const firstFactor = static_cast<double>(j - firstLineCoordinate) / static_cast<double>(i - firstLineCoordinate);
                double const secondFactor = 1.0 - firstFactor;

                // now distribute the shifting
                m_grid->m_gridNodes[nodeIndex.m_m][nodeIndex.m_n] = m_grid->m_gridNodes[nodeIndex.m_m][nodeIndex.m_n] +
                                                                    firstDelta * secondFactor + secondDelta * firstFactor;
                // field transformation in the influence area
            }
            firstLineCoordinate = i;
            firstDelta = secondDelta;
        }
    }

    return nullptr;
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
    auto const nodeIndex = m_gridCache->GetNodeIndices(fromPoint);

    //Check the nodes are on the line to shift

    if (!m_lines[0].IsNodeOnLine(nodeIndex))
    {
        throw std::invalid_argument("CurvilinearGridLineShift::MoveNode The selected node does not belong to the line to be shifted");
    }

    m_gridCache->m_gridNodes[nodeIndex.m_m][nodeIndex.m_n] = toPoint;
}