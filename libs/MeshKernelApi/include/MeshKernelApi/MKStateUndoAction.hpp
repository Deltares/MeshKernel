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

#include "MeshKernel/UndoActions/UndoAction.hpp"

#include "MeshKernelApi/State.hpp"

namespace meshkernelapi
{

    /// @brief Undo action for MeshKernelState.
    ///
    /// \note This is for creational api functions, such as mesh generation, mesh set, mesh conversion and state deallocate.
    /// The undo action manages a reference to the active MeshKernelState and keeps a copy of the pointers.
    /// When the pointer is updated in the mk api, undoing is achieved by swapping the pointers.
    class MKStateUndoAction : public meshkernel::UndoAction
    {
    public:
        /// @brief Allocate a CompoundUndoAction and return a unique_ptr to the newly create object.
        static std::unique_ptr<MKStateUndoAction> Create(MeshKernelState& mkState);

        /// @brief Constructor.
        ///
        /// Keeps a reference to the active MeshKernelState and copies its pointers
        explicit MKStateUndoAction(MeshKernelState& mkState);

    private:
        /// @brief Swap the pointer in the active state (referred to by m_mkStateReference) with the copies of the pointers.
        void SwapContents();

        /// @brief Commit undo action.
        void DoCommit() override;

        /// @brief Restore undo action
        void DoRestore() override;

        /// \brief Keeps a copy of the MeshKernelState, so that the pointers can be swapped with the active MeshKernelState on undo and redo.
        MeshKernelState m_mkState;

        /// \brief Reference to the active MeshKernelState.
        MeshKernelState& m_mkStateReference;
    };

} // namespace meshkernelapi
