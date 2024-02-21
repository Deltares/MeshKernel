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

namespace meshkernel
{

    // All undo actions must be created in the committed state
    // actions should not be copyable/moveable
    // change of state
    // Each action should be applied/restored as a single atomic action.
    class UndoAction
    {
    public:
        enum ActionState
        {
            Committed, ///< The action has been applied.
            Restored   ///< The action has been undone, the state immediately before the action occurred has been restored.
        };

        /// @brief Default constructor
        UndoAction() = default;

        /// @brief Delete the copy constructor
        UndoAction(const UndoAction& copy) = delete;

        /// @brief Delete the move constructor
        UndoAction(UndoAction&& copy) = delete;

        /// @brief Destructor
        virtual ~UndoAction() = default;

        /// @brief Delete the copy assignment operator
        UndoAction& operator=(const UndoAction& copy) = delete;

        /// @brief Delete the move assignment operator
        UndoAction& operator=(UndoAction&& copy) = delete;

        // TODO are the names commit and restore ok?
        // e.g. rollforward, rollback

        /// @brief Apply the changes and change the state to ActionState::Committed.
        void Commit();

        /// @brief Undo the changes and change the state to ActionState::Restored.
        void Restore();

        /// @brief Get the current state of the undo action.
        ActionState State() const;

    private:
        /// @brief Operation to apply the changes required by the UndoAction
        virtual void DoCommit() = 0;

        /// @brief Operation to restore the changes made by the UndoAction.
        ///
        /// The state of the object on which the action is to be restored should
        /// be set to ....
        virtual void DoRestore() = 0;

        /// @brief The current state of the action, all actions are constructed having already been applied
        ActionState m_state = Committed;
    };

    using UndoActionPtr = std::unique_ptr<UndoAction>;

} // namespace meshkernel

inline meshkernel::UndoAction::ActionState meshkernel::UndoAction::State() const
{
    return m_state;
}
