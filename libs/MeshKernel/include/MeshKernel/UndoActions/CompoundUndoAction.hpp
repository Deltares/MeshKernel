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

#include "MeshKernel/UndoActions/UndoAction.hpp"

namespace meshkernel
{
    /// @brief A composite of multiple undo actions
    class CompoundUndoAction : public UndoAction
    {
    public:
        // TODO how to make std::vector<std::unique_ptr<const UndoAction>>
        //
        /// @brief Iterator over composite undo actions.
        using const_iterator = std::vector<std::unique_ptr<UndoAction>>::const_iterator;

        /// @brief Allocate a CompoundUndoAction and return a unique_ptr to the newly create object.
        static std::unique_ptr<CompoundUndoAction> Create();

        /// @brief Add an undo action to the compound action
        void Add(UndoActionPtr&& action);

        /// @brief Iterator to start of composite undo actions
        const_iterator begin() const;

        /// @brief Iterator to one past end of composite undo actions
        const_iterator end() const;

        /// \brief Compute the approximate amount of memory being used, in bytes.
        std::uint64_t MemorySize() const override;

        /// @brief Print the compound undo action to the stream
        void Print(std::ostream& out = std::cout) const override;

    private:
        /// @brief Commit all undo actions.
        void DoCommit() override;

        /// @brief Restore all undo actions, in reverse order.
        void DoRestore() override;

        /// @brief A sequence of all the undo actions
        std::vector<UndoActionPtr> m_undoActions;
    };

} // namespace meshkernel

inline meshkernel::CompoundUndoAction::const_iterator meshkernel::CompoundUndoAction::begin() const
{
    return m_undoActions.begin();
}

inline meshkernel::CompoundUndoAction::const_iterator meshkernel::CompoundUndoAction::end() const
{
    return m_undoActions.end();
}
