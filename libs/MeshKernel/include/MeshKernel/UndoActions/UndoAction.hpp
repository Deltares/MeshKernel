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
#include <iostream>
#include <memory>
#include <string>

namespace meshkernel
{

    /// @brief A single unit of work performed, usually, on a mesh that can be undone and re-applied.
    ///
    /// Each undo-action and, possible, subsequent restore-action should be atomic, that is: an indivisble
    /// and irreducible series of performed operations, most likely, on a mesh.
    /// All undo-action's must be created in the committed state.
    class UndoAction
    {
    public:
        /// @brief The current state of an UndoAction
        enum class State
        {
            Committed, ///< The action has been applied.
            Restored   ///< The action has been undone, the state immediately before the action occurred has been restored.
        };

        /// @brief Return the string representation of the State enum value
        static const std::string& to_string(const State state);

        /// @brief Default constructor
        UndoAction() = default;

        /// @brief Delete the copy constructor
        UndoAction(const UndoAction& copy) = delete;

        /// @brief Destructor
        virtual ~UndoAction() = default;

        /// @brief Delete the copy assignment operator
        UndoAction& operator=(const UndoAction& copy) = delete;

        /// @brief Apply the changes and change the state to State::Committed.
        void Commit();

        /// @brief Undo the changes and change the state to State::Restored.
        void Restore();

        /// @brief Get the current state of the undo action.
        State GetState() const;

        /// \brief Compute the approximate amount of memory being used, in bytes.
        virtual std::uint64_t MemorySize() const;

        /// @brief Print the undo action to the stream
        virtual void Print(std::ostream& out = std::cout) const;

    private:
        /// @brief Operation to apply the changes required by the UndoAction
        virtual void DoCommit() = 0;

        /// @brief Operation to restore the changes made by the UndoAction.
        ///
        /// The state of the object on which the action is to be restored should
        /// be set to ....
        virtual void DoRestore() = 0;

        /// @brief The current state of the action, all actions are constructed having already been applied
        State m_state = State::Committed;
    };

    /// @brief Typedef of the pointer to an UndoAction.
    using UndoActionPtr = std::unique_ptr<UndoAction>;

} // namespace meshkernel
