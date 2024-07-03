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
        /// @brief Allows for a simplified insertion of an undo action in std::tie
        class StoreExpression
        {
        public:
            /// @brief Constructor
            StoreExpression(CompoundUndoAction& action) : m_undoAction(action) {}

            /// @brief Insert undo action into compound undo action sequence
            void operator=(UndoActionPtr&& action)
            {
                m_undoAction.Add(std::move(action));
            }

        private:
            /// @brief Reference to the compoind undo action object.
            CompoundUndoAction& m_undoAction;
        };

        /// @brief Iterator over composite undo actions.
        using const_iterator = std::vector<std::unique_ptr<UndoAction>>::const_iterator;

        /// @brief Allocate a CompoundUndoAction and return a unique_ptr to the newly create object.
        static std::unique_ptr<CompoundUndoAction> Create();

        /// @brief Constructor
        CompoundUndoAction() : m_storeExpression(*this) {}

        /// @brief Add an undo action to the compound action
        void Add(UndoActionPtr&& action);

        /// @brief Allows for insertion of an undo action into the compund undo action.
        StoreExpression& Insert();

        /// @brief Iterator to start of composite undo actions
        const_iterator begin() const;

        /// @brief Iterator to one past end of composite undo actions
        const_iterator end() const;

        /// \brief Compute the approximate amount of memory being used, in bytes.
        std::uint64_t MemorySize() const override;

    private:
        /// @brief Commit all undo actions.
        void DoCommit() override;

        /// @brief Restore all undo actions, in reverse order.
        void DoRestore() override;

        /// @brief A sequence of all the undo actions
        std::vector<UndoActionPtr> m_undoActions;

        /// @brief A store expresion, enables simple insertion of undo actions in functions returning tuples.
        StoreExpression m_storeExpression;
    };

} // namespace meshkernel

inline meshkernel::CompoundUndoAction::StoreExpression& meshkernel::CompoundUndoAction::Insert()
{
    return m_storeExpression;
}
