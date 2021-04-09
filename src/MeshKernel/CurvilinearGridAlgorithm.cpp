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
#include <MeshKernel/CurvilinearGridAlgorithm.hpp>
#include <MeshKernel/CurvilinearGridLine.hpp>
#include <MeshKernel/Exceptions.hpp>

using meshkernel::CurvilinearGrid;
using meshkernel::CurvilinearGridAlgorithm;
using meshkernel::CurvilinearGridLine;

CurvilinearGridAlgorithm::CurvilinearGridAlgorithm(std::shared_ptr<CurvilinearGrid> grid)
    : m_grid(grid)

{
}

void CurvilinearGridAlgorithm::SetBlock(Point const& firstCornerPoint, Point const& secondCornerPoint)
{
    // Get the m_m and m_n indices from the point coordinates
    auto const [lowerLeft, upperRight] = m_grid->ComputeBlockFromCornerPoints(firstCornerPoint, secondCornerPoint);

    // Coinciding corner nodes, no valid area, nothing to do
    if (lowerLeft == upperRight)
    {
        throw std::invalid_argument("CurvilinearGridSmoothing::SetBlock coinciding corner nodes, no valid area to smooth");
    }

    m_lowerLeft = lowerLeft;
    m_upperRight = upperRight;
}

void CurvilinearGridAlgorithm::SetLine(Point const& firstPoint, Point const& secondPoint)
{
    // The selected nodes must be on the vertical or horizontal line
    const auto [newLineLowerLeft, newLineUpperRight] = m_grid->ComputeBlockFromCornerPoints(firstPoint, secondPoint);

    // Coinciding nodes, no valid line, nothing to do
    if (newLineLowerLeft == newLineUpperRight)
    {
        throw std::invalid_argument("CurvilinearGridAlgorithm::SetLine the start and the end points of the line are coinciding");
    }

    // The points of the frozen line, must be on the same grid-line
    if (!newLineLowerLeft.IsOnTheSameGridLine(newLineUpperRight))
    {
        throw std::invalid_argument("CurvilinearGridAlgorithm::SetLine the current line is intersection another line");
    }

    CurvilinearGridLine const newGridline{newLineLowerLeft, newLineUpperRight};

    // Frozen lines cannot cross existing frozen lines
    for (auto const& frozenLine : m_lines)
    {
        for (auto i = frozenLine.m_startCoordinate; i <= frozenLine.m_endCoordinate; ++i)
        {
            for (auto j = newGridline.m_startCoordinate; j <= newGridline.m_endCoordinate; ++j)
            {
                if (j == frozenLine.m_constantCoordinate && i == newGridline.m_constantCoordinate)
                {
                    throw AlgorithmError("CurvilinearGridOrthogonalization::SetFrozenLine the new line to freeze is crossing an existing line");
                }
            }
        }
    }

    // Now a frozen line can be stored
    m_lines.emplace_back(newGridline);
}
