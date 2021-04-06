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
#include <MeshKernel/Exceptions.hpp>

using meshkernel::CurvilinearGrid;
using meshkernel::CurvilinearGridAlgorithm;

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

    // The start-end coordinates of the new line, along m_m or m_n direction
    auto const isNewlineAnMLine = newLineUpperRight.m_m == newLineLowerLeft.m_m;
    auto const startNewLine = isNewlineAnMLine ? newLineLowerLeft.m_n : newLineLowerLeft.m_m;
    auto const endNewLine = isNewlineAnMLine ? newLineUpperRight.m_n : newLineUpperRight.m_m;
    auto const constantCoordinateLine = isNewlineAnMLine ? newLineLowerLeft.m_m : newLineLowerLeft.m_n;

    // Frozen lines cannot cross existing frozen lines
    for (auto const& frozenLine : m_lines)
    {
        auto const& [lowerLeft, upperRight] = frozenLine;

        auto const deltaMCurrentLine = upperRight.m_m - lowerLeft.m_m;
        auto const startCurrentLine = deltaMCurrentLine == 0 ? lowerLeft.m_n : lowerLeft.m_m;
        auto const endCurrentLine = deltaMCurrentLine == 0 ? upperRight.m_n : upperRight.m_m;
        auto const constantCoordinateCurrentLine = deltaMCurrentLine == 0 ? lowerLeft.m_m : lowerLeft.m_n;
        for (auto i = startCurrentLine; i <= endCurrentLine; ++i)
        {
            for (auto j = startNewLine; j <= endNewLine; ++j)
            {
                if (j == constantCoordinateCurrentLine && i == constantCoordinateLine)
                {
                    throw AlgorithmError("CurvilinearGridOrthogonalization::SetFrozenLine the new line to freeze is crossing an existing line");
                }
            }
        }
    }

    // Now a frozen line can be stored
    m_lines.emplace_back(newLineLowerLeft, newLineUpperRight);
}
