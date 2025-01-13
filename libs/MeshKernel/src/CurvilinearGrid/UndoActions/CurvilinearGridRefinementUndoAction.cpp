//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2025.
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

#include "MeshKernel/CurvilinearGrid/UndoActions/CurvilinearGridRefinementUndoAction.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"

std::unique_ptr<meshkernel::CurvilinearGridRefinementUndoAction> meshkernel::CurvilinearGridRefinementUndoAction::Create(CurvilinearGrid& grid)
{
    return std::make_unique<CurvilinearGridRefinementUndoAction>(grid);
}

meshkernel::CurvilinearGridRefinementUndoAction::CurvilinearGridRefinementUndoAction(CurvilinearGrid& grid)
    : BaseMeshUndoAction<CurvilinearGridRefinementUndoAction, CurvilinearGrid>(grid), m_nodes(grid.GetNodes()), m_startOffset(grid.StartOffset()), m_endOffset(grid.EndOffset()) {}

void meshkernel::CurvilinearGridRefinementUndoAction::Swap(lin_alg::Matrix<Point>& nodes, CurvilinearGridNodeIndices& startOffset, CurvilinearGridNodeIndices& endOffset)
{
    std::swap(m_startOffset, startOffset);
    std::swap(m_endOffset, endOffset);
    m_nodes.swap(nodes);
}

std::uint64_t meshkernel::CurvilinearGridRefinementUndoAction::MemorySize() const
{
    std::uint64_t result = 0;
    result += sizeof(*this);
    result += static_cast<std::uint64_t>(m_nodes.rows() * m_nodes.cols() * sizeof(Point));
    return result;
}
