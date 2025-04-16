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

#include <utility>

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

    auto startNode = m_lines[0].m_startNode;
    auto endNode = m_lines[0].m_endNode;
    m_grid.ComputeGridNodeTypes();
    auto [numAddedLines, gridLineType, addLinesUndoAction] = m_grid.AddGridLinesAtBoundary(startNode, endNode, m_numRowsToMirror);

    auto undoAction = CompoundUndoAction::Create();
    undoAction->Add(std::move(addLinesUndoAction));

    double const a = 1.0 + m_mirroringFactor;
    double const b = -m_mirroringFactor;

    CurvilinearGridNodeIndices lowerLeft;
    CurvilinearGridNodeIndices upperRight;

    using enum CurvilinearGrid::BoundaryGridLineType;

    switch (gridLineType)
    {
    case Bottom:
        lowerLeft = {m_lines[0].m_constantCoordinate + numAddedLines - m_numRowsToMirror, m_lines[0].m_startCoordinate};
        upperRight = {m_lines[0].m_constantCoordinate + numAddedLines, m_lines[0].m_endCoordinate + 1};
        break;
    case Top:
        lowerLeft = {m_lines[0].m_constantCoordinate + 1, m_lines[0].m_startCoordinate};
        upperRight = {m_lines[0].m_constantCoordinate + m_numRowsToMirror + 1, m_lines[0].m_endCoordinate + 1};
        break;
    case Right:
        lowerLeft = {m_lines[0].m_startCoordinate, m_lines[0].m_constantCoordinate + 1};
        upperRight = {m_lines[0].m_endCoordinate + 1, m_lines[0].m_constantCoordinate + m_numRowsToMirror + 1};
        break;
    case Left:
        lowerLeft = {m_lines[0].m_startCoordinate, m_lines[0].m_constantCoordinate + numAddedLines - m_numRowsToMirror};
        upperRight = {m_lines[0].m_endCoordinate + 1, m_lines[0].m_constantCoordinate + numAddedLines};
        break;
    default:
        throw std::invalid_argument("CurvilinearGridLineMirror:: Invalid grid line type");
    }

    undoAction->Add(CurvilinearGridBlockUndoAction::Create(m_grid, lowerLeft, upperRight));

    UInt boundaryCoordinate;

    switch (gridLineType)
    {
    case Bottom:
        boundaryCoordinate = m_lines[0].m_constantCoordinate + numAddedLines;
        for (int r = 0; r < m_numRowsToMirror; ++r)
        {
            for (auto i = m_lines[0].m_startCoordinate; i <= m_lines[0].m_endCoordinate; ++i)
            {
                m_grid.GetNode(boundaryCoordinate - 1, i) = m_grid.GetNode(boundaryCoordinate, i) * a +
                                                            m_grid.GetNode(boundaryCoordinate + 1, i) * b;
            }

            boundaryCoordinate--;
        }
        break;

    case Top:
        boundaryCoordinate = m_lines[0].m_constantCoordinate;
        for (int r = 0; r < m_numRowsToMirror; ++r)
        {
            for (auto i = m_lines[0].m_startCoordinate; i <= m_lines[0].m_endCoordinate; ++i)
            {
                m_grid.GetNode(boundaryCoordinate + 1, i) = m_grid.GetNode(boundaryCoordinate, i) * a +
                                                            m_grid.GetNode(boundaryCoordinate - 1, i) * b;
            }

            boundaryCoordinate++;
        }
        break;

    case Right:
        boundaryCoordinate = m_lines[0].m_constantCoordinate;
        for (int r = 0; r < m_numRowsToMirror; ++r)
        {
            for (auto i = m_lines[0].m_startCoordinate; i <= m_lines[0].m_endCoordinate; ++i)
            {
                m_grid.GetNode(i, boundaryCoordinate + 1) = m_grid.GetNode(i, boundaryCoordinate) * a +
                                                            m_grid.GetNode(i, boundaryCoordinate - 1) * b;
            }

            boundaryCoordinate++;
        }
        break;

    case Left:
        boundaryCoordinate = m_lines[0].m_constantCoordinate + numAddedLines;
        for (int r = 0; r < m_numRowsToMirror; ++r)
        {
            for (auto i = m_lines[0].m_startCoordinate; i <= m_lines[0].m_endCoordinate; ++i)
            {
                m_grid.GetNode(i, boundaryCoordinate - 1) = m_grid.GetNode(i, boundaryCoordinate) * a +
                                                            m_grid.GetNode(i, boundaryCoordinate + 1) * b;
            }

            boundaryCoordinate--;
        }
        break;

    default:
        throw std::invalid_argument("CurvilinearGridLineMirror::Compute Invalid grid line type");
    }

    return undoAction;
}
