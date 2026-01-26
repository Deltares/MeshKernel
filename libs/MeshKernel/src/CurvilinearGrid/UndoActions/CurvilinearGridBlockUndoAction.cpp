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

#include "MeshKernel/CurvilinearGrid/UndoActions/CurvilinearGridBlockUndoAction.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"

std::unique_ptr<meshkernel::CurvilinearGridBlockUndoAction> meshkernel::CurvilinearGridBlockUndoAction::Create(CurvilinearGrid& grid,
                                                                                                               const CurvilinearGridNodeIndices& startOffset,
                                                                                                               const CurvilinearGridNodeIndices& endOffset)
{
    struct DerivedCurvilinearGridBlockUndoAction : public CurvilinearGridBlockUndoAction
    {
        DerivedCurvilinearGridBlockUndoAction(CurvilinearGrid& grid,
                                              const CurvilinearGridNodeIndices& startOffset,
                                              const CurvilinearGridNodeIndices& endOffset) : CurvilinearGridBlockUndoAction(grid, startOffset, endOffset) {}
    };

    UInt endMOffset = endOffset.m_m + (grid.NumM() == endOffset.m_m ? 0 : 1);
    UInt endNOffset = endOffset.m_n + (grid.NumN() == endOffset.m_n ? 0 : 1);

    CurvilinearGridNodeIndices endOffset2(endNOffset, endMOffset);

    return std::make_unique<DerivedCurvilinearGridBlockUndoAction>(grid, startOffset, endOffset2);
}

std::unique_ptr<meshkernel::CurvilinearGridBlockUndoAction> meshkernel::CurvilinearGridBlockUndoAction::Create(CurvilinearGrid& grid)
{
    return Create(grid, CurvilinearGridNodeIndices(0, 0), CurvilinearGridNodeIndices(grid.NumN(), grid.NumM()));
}

meshkernel::CurvilinearGridBlockUndoAction::CurvilinearGridBlockUndoAction(CurvilinearGrid& grid,
                                                                           const CurvilinearGridNodeIndices& startOffset,
                                                                           const CurvilinearGridNodeIndices& endOffset)
    : BaseMeshUndoAction<CurvilinearGridBlockUndoAction, CurvilinearGrid>(grid), m_block(startOffset, endOffset)
{
    m_block.CopyFrom(grid);
}

void meshkernel::CurvilinearGridBlockUndoAction::Swap(CurvilinearGrid& grid)
{
    m_block.Swap(grid);
}

std::uint64_t meshkernel::CurvilinearGridBlockUndoAction::MemorySize() const
{
    return static_cast<std::uint64_t>(sizeof(*this)) + m_block.MemorySize();
}
