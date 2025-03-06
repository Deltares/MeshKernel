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

CurvilinearGridLineMirror::CurvilinearGridLineMirror(CurvilinearGrid& grid, double mirroringFactor) : CurvilinearGridAlgorithm(grid), m_mirroringFactor(mirroringFactor)

{
    if (m_mirroringFactor <= 0)
    {
        throw std::invalid_argument("CurvilinearGridLineMirror::CurvilinearGridLineMirror mirroring factor cannot be less or equal to zero");
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

    std::unique_ptr<CompoundUndoAction> undoAction = CompoundUndoAction::Create();

    const auto startNode = m_lines[0].m_startNode;
    const auto endNode = m_lines[0].m_endNode;

    m_grid.ComputeGridNodeTypes();
    auto const gridLineType = m_grid.GetBoundaryGridLineType(startNode, endNode);

    auto [addedLine, addGridLineAction] = m_grid.AddGridLineAtBoundary(startNode, endNode);
    m_grid.ComputeGridNodeTypes();
    undoAction->Add(std::move(addGridLineAction));

    double const a = 1.0 + m_mirroringFactor;
    double const b = -m_mirroringFactor;

    CurvilinearGridNodeIndices lowerLeft;
    CurvilinearGridNodeIndices upperRight;

    using enum CurvilinearGrid::BoundaryGridLineType;

    switch (gridLineType)
    {
    case Bottom:
        lowerLeft = {0, m_lines[0].m_startCoordinate};
        upperRight = {1, m_lines[0].m_endCoordinate + 1};
        break;
    case Top:
        lowerLeft = {m_grid.NumN() - 1, m_lines[0].m_startCoordinate};
        upperRight = {m_grid.NumN() - 0, m_lines[0].m_endCoordinate + 1};
        break;
    case Right:
        lowerLeft = {m_lines[0].m_startCoordinate, m_grid.NumM() - 1};
        upperRight = {m_lines[0].m_endCoordinate + 1, m_grid.NumM() - 0};
        break;
    case Left:
        lowerLeft = {0, m_lines[0].m_startCoordinate};
        upperRight = {m_lines[0].m_endCoordinate + 1, m_lines[0].m_startCoordinate + 1};
        break;
    default:
        throw std::invalid_argument("CurvilinearGridLineMirror:: Invalid grid line type");
    }

    undoAction->Add(CurvilinearGridBlockUndoAction::Create(m_grid, lowerLeft, upperRight));
    const auto boundaryCoordinate = m_lines[0].m_constantCoordinate == 0 ? 1 : m_lines[0].m_constantCoordinate;

    switch (gridLineType)
    {
    case Bottom:
        for (auto i = m_lines[0].m_startCoordinate; i <= m_lines[0].m_endCoordinate; ++i)
        {
            m_grid.GetNode(boundaryCoordinate - 1, i) = m_grid.GetNode(boundaryCoordinate, i) * a + m_grid.GetNode(boundaryCoordinate + 1, i) * b;
        }
        break;

    case Top:
        for (auto i = m_lines[0].m_startCoordinate; i <= m_lines[0].m_endCoordinate; ++i)
        {
            m_grid.GetNode(boundaryCoordinate + 1, i) = m_grid.GetNode(boundaryCoordinate, i) * a + m_grid.GetNode(boundaryCoordinate - 1, i) * b;
        }
        break;

    case Right:
        for (auto i = m_lines[0].m_startCoordinate; i <= m_lines[0].m_endCoordinate; ++i)
        {
            m_grid.GetNode(i, boundaryCoordinate + 1) = m_grid.GetNode(i, boundaryCoordinate) * a + m_grid.GetNode(i, boundaryCoordinate - 1) * b;
        }
        break;

    case Left:
        for (auto i = m_lines[0].m_startCoordinate; i <= m_lines[0].m_endCoordinate; ++i)
        {
            m_grid.GetNode(i, boundaryCoordinate - 1) = m_grid.GetNode(i, boundaryCoordinate) * a + m_grid.GetNode(i, boundaryCoordinate + 1) * b;
        }
        break;
    default:
        throw std::invalid_argument("CurvilinearGridLineMirror::Compute Invalid grid line type");
    }

    return undoAction;
}
