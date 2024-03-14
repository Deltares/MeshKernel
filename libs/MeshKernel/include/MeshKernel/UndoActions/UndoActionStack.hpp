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

#include <cstdint>
#include <list>
#include <utility>
#include <vector>

#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/UndoActions/UndoAction.hpp"

namespace meshkernel
{

    /// @brief A stack of UndoActions.
    ///
    /// Actions can be undo and subsequently re-done.
    class UndoActionStack
    {
    public:
        /// @brief Maximum number of undo action items
        static const UInt MaxUndoSize;

        // Add info about the action, both short and long form
        // Perhaps Short for for menu items, long form for tooltips?
        // long form for any exceptions?

        /// @brief Add an UndoAction.
        ///
        /// All added undo-actions must be in the committed state, if not then a ConstraintError
        /// will be raised.
        /// No null undo-actions will be added to the stack.
        /// All restored items will be removed, since after adding a new undo-action they are no
        /// longer restore-able.
        void Add(UndoActionPtr&& transaction);

        /// @brief Undo the action at the top of the committed stack
        ///
        /// The undo-action will be moved to the restored stack in case it should be re-done.
        /// \returns true if an undo-action was performed, false otherwise
        bool Undo();

        // Another name
        /// @brief Redo the action at the top of the restored stack.
        ///
        /// The undo-action will be moved to the committed stack in case it needs to be undone.
        /// \returns true if an redo-action was performed, false otherwise
        bool Commit();

        /// @brief Clear all undo actions.
        void Clear();

        /// @brief Get the number of undo action items
        ///
        /// The total includes both the number of committed and restored actions.
        UInt Size() const;

        /// \brief Compute the approximate amount of memory being used, in bytes, for all undo actions.
        std::uint64_t MemorySize() const;

    private:
        /// @brief Stack of committed undo actions
        std::list<UndoActionPtr> m_committed;

        /// @brief Stack of restored undo actions
        std::list<UndoActionPtr> m_restored;
    };

} // namespace meshkernel

inline meshkernel::UInt meshkernel::UndoActionStack::Size() const
{
    return static_cast<UInt>(m_committed.size() + m_restored.size());
}
