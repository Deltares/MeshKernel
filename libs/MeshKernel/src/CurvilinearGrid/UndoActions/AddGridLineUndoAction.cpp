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

#include "MeshKernel/CurvilinearGrid/UndoActions/AddGridLineUndoAction.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"

std::unique_ptr<meshkernel::AddGridLineUndoAction> meshkernel::AddGridLineUndoAction::Create(CurvilinearGrid& grid,
                                                                                             const CurvilinearGridNodeIndices& startOffset,
                                                                                             const CurvilinearGridNodeIndices& endOffset)
{
    return std::make_unique<AddGridLineUndoAction>(grid, startOffset, endOffset);
}

meshkernel::AddGridLineUndoAction::AddGridLineUndoAction(CurvilinearGrid& grid,
                                                         const CurvilinearGridNodeIndices& startOffset,
                                                         const CurvilinearGridNodeIndices& endOffset)
    : BaseMeshUndoAction<AddGridLineUndoAction, CurvilinearGrid>(grid),
      m_startOffset(startOffset),
      m_endOffset(endOffset) {}
