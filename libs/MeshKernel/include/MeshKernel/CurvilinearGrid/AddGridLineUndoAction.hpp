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

#include "MeshKernel/CurvilinearGrid/CurvilinearGridNodeIndices.hpp"
#include "MeshKernel/UndoActions/BaseMeshUndoAction.hpp"

namespace meshkernel
{
    /// @brief Forward declaration of the curvilinear mesh
    class CurvilinearGrid;

    /// @brief Undo action for adding grid line to boundary of curvilinear mesh
    class AddGridLineUndoAction : public BaseMeshUndoAction<AddGridLineUndoAction, CurvilinearGrid>
    {
    public:
        /// @brief Return unique pointer to newly created AddGridLineUndoAction object
        static std::unique_ptr<AddGridLineUndoAction> Create(CurvilinearGrid& grid,
                                                             const CurvilinearGridNodeIndices& startOffset,
                                                             const CurvilinearGridNodeIndices& endOffset);

        /// @brief Constructor
        AddGridLineUndoAction(CurvilinearGrid& grid,
                              const CurvilinearGridNodeIndices& startOffset,
                              const CurvilinearGridNodeIndices& endOffset);

        /// @brief Get the start offset
        const CurvilinearGridNodeIndices& StartOffset() const;

        /// @brief Get the end offset
        const CurvilinearGridNodeIndices& EndOffset() const;

        /// @brief Print the add grid line action to the stream
        void Print(std::ostream& out = std::cout) const override;

    private:
        /// @brief Start offset, to be added when restoring action and subtracted when undoing action.
        CurvilinearGridNodeIndices m_startOffset;

        /// @brief End offset, to be added when restoring action and subtracted when undoing action.
        CurvilinearGridNodeIndices m_endOffset;
    };

} // namespace meshkernel

inline const meshkernel::CurvilinearGridNodeIndices& meshkernel::AddGridLineUndoAction::StartOffset() const
{
    return m_startOffset;
}

inline const meshkernel::CurvilinearGridNodeIndices& meshkernel::AddGridLineUndoAction::EndOffset() const
{
    return m_endOffset;
}
