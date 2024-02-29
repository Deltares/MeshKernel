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

    class UndoActionStack
    {
    public:
        // When adding a new transaction, should probably delete all restored transactions
        // When adding new transactions, could check the size of the committed list and remove transactions more than some number ago
        // e.g. keep the undo list no longer than 10

        // Add info about the action, both short and long form
        // Perhaps Short for for menu items, long form for tooltips?
        // long form for any exceptions?

        // Do we need a clear function? If there are operations on a mesh that
        // are not undo-able then the all undo actions should be removed.

        /// @brief Constructor
        UndoActionStack();

        // When adding a new transaction, should probably delete all restored transactions
        // actions should be moved so that there can only be a single reference to the same transaction
        void Add(UndoActionPtr&& transaction);

        /// @brief Undo the action at the top of the committed stack
        bool Undo();

        /// @brief Redo the action at the top of the restored stack
        // Another name
        bool Commit();

    private:
        /// @brief The initial reserved size of the committed and restored undo-action arrays
        static const UInt DefaultReserveSize = 10;

        std::vector<UndoActionPtr> m_committed;
        std::vector<UndoActionPtr> m_restored;
    };

} // namespace meshkernel
