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

#include <utility>
#include <vector>

#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/UndoAction.hpp"

namespace meshkernel
{

    /// @brief A stack of UndoActions.
    ///
    /// Actions can be undo and subsequently re-done.
    class UndoActionStack
    {
    public:
        // When adding new transactions, could check the size of the committed list and remove transactions more than some number ago
        // e.g. keep the undo list no longer than 10

        // Add info about the action, both short and long form
        // Perhaps Short for for menu items, long form for tooltips?
        // long form for any exceptions?

        // Do we need a clear function? If there are operations on a mesh that
        // are not undo-able then the all undo actions should be removed.

        /// @brief Constructor
        UndoActionStack();

        /// @brief Add an UndoAction.
        ///
        /// All added undo-actions must be in the committed state, if not then a ConstraintError
        /// will be raised.
        /// No null undo-actions will be added to the stack.
        /// All actions that have be restored will be deleted.
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

    private:
        /// @brief The initial reserved size of the committed and restored undo-action arrays
        static const UInt DefaultReserveSize = 10;

        /// @brief Stack of committed undo actions
        std::vector<UndoActionPtr> m_committed;

        /// @brief Stack of restored undo actions
        std::vector<UndoActionPtr> m_restored;
    };

} // namespace meshkernel
