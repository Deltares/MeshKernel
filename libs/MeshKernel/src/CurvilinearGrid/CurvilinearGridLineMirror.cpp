//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2021.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation version 3.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without seteven the implied warranty of
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

#include <utility>
#include <iomanip>

#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridLine.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridLineMirror.hpp"

using meshkernel::CurvilinearGrid;
using meshkernel::CurvilinearGridLineMirror;
using meshkernel::Point;

CurvilinearGridLineMirror::CurvilinearGridLineMirror(CurvilinearGrid& grid,
                                                     double mirroringFactor,
                                                     int numRowsToMirror) : CurvilinearGridAlgorithm(grid),
                                                                            m_mirroringFactor(mirroringFactor),
                                                                            m_numRowsToMirror(numRowsToMirror)

{
    if (m_mirroringFactor <= 0)
    {
        throw std::invalid_argument("CurvilinearGridLineMirror::CurvilinearGridLineMirror mirroring factor cannot be less or equal to zero");
    }
    if (m_numRowsToMirror <= 0)
    {
        throw std::invalid_argument("CurvilinearGridLineMirror::CurvilinearGridLineMirror the number of rows to mirror cannot be less or equal to zero");
    }
}

meshkernel::UndoActionPtr CurvilinearGridLineMirror::Compute()
{
    if (m_lines.empty())
    {
        throw std::invalid_argument("CurvilinearGridLineMirror::Compute No candidate line to shift has been selected");
    }
    if (!m_grid.IsValid())
    {
        throw std::invalid_argument("CurvilinearGridLineMirror:: Invalid curvilinear grid");
    }

    std::cout << "CurvilinearGridLineMirror::Compute "<< std::endl;

    CurvilinearGridNodeIndices startOffset = m_grid.StartOffset();
    CurvilinearGridNodeIndices endOffset   = m_grid.EndOffset();

    const auto startNode = m_lines[0].m_startNode;
    const auto endNode = m_lines[0].m_endNode;
    m_grid.ComputeGridNodeTypes();
    auto [numAddedLines, addLinesUndoAction] = m_grid.AddGridLinesAtBoundary(startNode, endNode, m_numRowsToMirror);

    CurvilinearGridNodeIndices startOffset2 = m_grid.StartOffset();
    CurvilinearGridNodeIndices endOffset2   = m_grid.EndOffset();

    std::cout << "start offset: " << startOffset.m_n << ", " << startOffset.m_m << " -- " << startOffset2.m_n << ", " << startOffset2.m_m << std::endl;
    std::cout << "end offset  : " << endOffset.m_n << ", " << endOffset.m_m << " -- " << endOffset2.m_n << ", " << endOffset2.m_m << std::endl;

    auto undoAction = CompoundUndoAction::Create();
    undoAction->Add(std::move(addLinesUndoAction));

    double const a = 1.0 + m_mirroringFactor;
    double const b = -m_mirroringFactor;

    CurvilinearGridNodeIndices lowerLeft;
    CurvilinearGridNodeIndices upperRight;

    if (startOffset.m_n > startOffset2.m_n) {
    }


    using enum CurvilinearGrid::BoundaryGridLineType;

    auto const gridLineType = m_grid.GetBoundaryGridLineType(startNode, endNode);
    switch (gridLineType)
    {
    case Bottom:
        lowerLeft =  {m_lines[0].m_constantCoordinate + numAddedLines - m_numRowsToMirror, m_lines[0].m_startCoordinate};
        upperRight = {m_lines[0].m_constantCoordinate + numAddedLines,                     m_lines[0].m_endCoordinate + 1};
        break;
    case Top:
        lowerLeft =  {m_lines[0].m_constantCoordinate + 1,                     m_lines[0].m_startCoordinate};
        upperRight = {m_lines[0].m_constantCoordinate + m_numRowsToMirror + 1, m_lines[0].m_endCoordinate + 1};
        break;
    case Right:
        lowerLeft  = {m_lines[0].m_startCoordinate,   m_lines[0].m_constantCoordinate + 1};
        upperRight = {m_lines[0].m_endCoordinate + 1, m_lines[0].m_constantCoordinate + m_numRowsToMirror + 1};
        break;
    case Left:
        lowerLeft  = {m_lines[0].m_startCoordinate,   m_lines[0].m_constantCoordinate + numAddedLines - m_numRowsToMirror};
        upperRight = {m_lines[0].m_endCoordinate + 1, m_lines[0].m_constantCoordinate + numAddedLines};
        break;
    default:
        throw std::invalid_argument("CurvilinearGridLineMirror:: Invalid grid line type");
    }

    std::cout << " range " << m_lines[0].m_constantCoordinate << "  " << numAddedLines << "  " << m_numRowsToMirror
              << " {" << lowerLeft.m_n << ", " << lowerLeft.m_m << "} {"
              << upperRight.m_n << ", " << upperRight.m_m << "}" << std::endl;

    undoAction->Add(CurvilinearGridBlockUndoAction::Create(m_grid, lowerLeft, upperRight));
    UInt boundaryCoordinate;
    switch (gridLineType)
    {
    case Bottom:
        std::cout << "Bottom"<< std::endl;
        boundaryCoordinate = m_lines[0].m_constantCoordinate + numAddedLines;
        for (int r = 0; r < m_numRowsToMirror; ++r)
        {
            for (auto i = m_lines[0].m_startCoordinate; i <= m_lines[0].m_endCoordinate; ++i)
            {
                m_grid.GetNode(boundaryCoordinate - 1, i) = m_grid.GetNode(boundaryCoordinate, i) * a +
                                                            m_grid.GetNode(boundaryCoordinate + 1, i) * b;
            std::cout << "{" << std::setw (5) << boundaryCoordinate - 1 << ", " << std::setw (5) << i << "} ";
            }

            std::cout << std::endl;
            boundaryCoordinate--;
        }
        break;
    case Top:
        std::cout << "Top"<< std::endl;
        boundaryCoordinate = m_lines[0].m_constantCoordinate;
        for (int r = 0; r < m_numRowsToMirror; ++r)
        {
            for (auto i = m_lines[0].m_startCoordinate; i <= m_lines[0].m_endCoordinate; ++i)
            {
                m_grid.GetNode(boundaryCoordinate + 1, i) = m_grid.GetNode(boundaryCoordinate, i) * a +
                                                            m_grid.GetNode(boundaryCoordinate - 1, i) * b;
            std::cout << "{" << std::setw (5) << boundaryCoordinate + 1 << ", " << std::setw (5) << i << "} ";
            }

            std::cout << std::endl;
            boundaryCoordinate++;
        }
        break;

    case Right:
        std::cout << "Right"<< std::endl;
        boundaryCoordinate = m_lines[0].m_constantCoordinate;
        for (int r = 0; r < m_numRowsToMirror; ++r)
        {
            for (auto i = m_lines[0].m_startCoordinate; i <= m_lines[0].m_endCoordinate; ++i)
            {
                m_grid.GetNode(i, boundaryCoordinate + 1) = m_grid.GetNode(i, boundaryCoordinate) * a +
                                                            m_grid.GetNode(i, boundaryCoordinate - 1) * b;

                std::cout << "{" << std::setw (5) << i << ", " << std::setw (5) << boundaryCoordinate + 1 << "} ";

            }
            std::cout << std::endl;
            boundaryCoordinate++;
        }
        break;

    case Left:
        std::cout << "Left"<< std::endl;
        boundaryCoordinate = m_lines[0].m_constantCoordinate + numAddedLines;
        for (int r = 0; r < m_numRowsToMirror; ++r)
        {
            for (auto i = m_lines[0].m_startCoordinate; i <= m_lines[0].m_endCoordinate; ++i)
            {
                m_grid.GetNode(i, boundaryCoordinate - 1) = m_grid.GetNode(i, boundaryCoordinate) * a +
                                                            m_grid.GetNode(i, boundaryCoordinate + 1) * b;
                std::cout << "{" << std::setw (5) << i << ", " << std::setw (5) << boundaryCoordinate - 1 << "} ";
            }
            std::cout << std::endl;
            boundaryCoordinate--;
        }
        break;
    default:
        throw std::invalid_argument("CurvilinearGridLineMirror::Compute Invalid grid line type");
    }

    return undoAction;
}
