//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2024.
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

#include <memory>
#include <utility>
#include <vector>

#include "MeshKernel/CurvilinearGrid/CurvilinearGridBlock.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridNodeIndices.hpp"
#include "MeshKernel/UndoActions/BaseMeshUndoAction.hpp"

namespace meshkernel
{
    // Forward declaration of the curvilinear mesh
    class CurvilinearGrid;

    class CurvilinearGridBlockUndo : public BaseMeshUndoAction<CurvilinearGridBlockUndo, CurvilinearGrid>
    {
    public:
        static std::unique_ptr<CurvilinearGridBlockUndo> Create(CurvilinearGrid& grid,
                                                                const CurvilinearGridNodeIndices& startOffset,
                                                                const CurvilinearGridNodeIndices& endOffset);

        // Not all is needed
        // which direction, which side and the offset of 1 is implied
        CurvilinearGridBlockUndo(CurvilinearGrid& grid,
                                 const CurvilinearGridNodeIndices& startOffset,
                                 const CurvilinearGridNodeIndices& endOffset);

        void Swap(CurvilinearGrid& grid);

        /// @brief Print the add grid line action to the stream
        void Print(std::ostream& out = std::cout) const override;

    private:
        CurvilinearGridBlock m_block;
    };

} // namespace meshkernel
