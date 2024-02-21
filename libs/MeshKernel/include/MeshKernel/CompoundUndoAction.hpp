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

#include "MeshKernel/UndoAction.hpp"

namespace meshkernel
{
    /// @brief A composite of multiple undo actions
    class CompoundUndoAction : public UndoAction
    {
    public:
        /// @brief Allocate a CompoundUndoAction and return a unique_ptr to the newly create object.
        static std::unique_ptr<CompoundUndoAction> Create();

        /// @brief Add an undo action to the compound action
        void Add(UndoActionPtr&& action);

    private:
        /// @brief Commit all undo actions.
        void DoCommit();

        /// @brief Restore all undo actions, in reverse order.
        void DoRestore();

        /// @brief A sequence of all the undo actions
        std::vector<UndoActionPtr> m_undoActions;
    };

} // namespace meshkernel
