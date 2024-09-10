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

#include "MeshKernel/UndoActions/UndoAction.hpp"

namespace meshkernel
{

    /// @brief A very simple undo action that makes no changes
    ///
    /// This is required because null undo action pointers will not be added to the undo stack.
    class NoActionUndo : public UndoAction
    {
    public:
        /// @brief Allocate a NoActionUndo and return a unique_ptr to the newly created object.
        static std::unique_ptr<NoActionUndo> Create();

        /// @brief Default constructor
        NoActionUndo() = default;

    private:
        /// @brief Perform the undo action, does not change anything
        void DoCommit() override;

        /// @brief Perform the redo action, does not change anything
        void DoRestore() override;
    };

} // namespace meshkernel
