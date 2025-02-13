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
#include "MeshKernelApi/State.hpp"

namespace meshkernelapi
{

    /// @brief Undo action for deleting a frozen lines
    class CurvilinearFrozenLinesDeleteUndoAction : public meshkernel::UndoAction
    {
    public:
        /// @brief Constructor.
        ///
        /// Keeps a reference to the active MeshKernelState, the last used frozenLineId and frozen line points
        explicit CurvilinearFrozenLinesDeleteUndoAction(MeshKernelState& mkState,
                                                        int frozenLineId,
                                                        const std::pair<meshkernel::Point, meshkernel::Point>& frozenLinePoints) : m_mkStateReference(mkState),
                                                                                                                                   m_frozenLinesCounter(frozenLineId),
                                                                                                                                   m_frozenLinePoints(frozenLinePoints)
        {
        }

    private:
        /// @brief Restoring delete
        void DoCommit() override
        {
            if (!m_mkStateReference.m_frozenLines.contains(m_frozenLinesCounter))
            {
                throw meshkernel::MeshKernelError("Frozen line counter in meshkernel state should exist when committing a deletion of a frozen line");
            }

            m_mkStateReference.m_frozenLines.erase(m_frozenLinesCounter);
        }

        /// @brief Reverting delete
        void DoRestore() override
        {
            if (m_mkStateReference.m_frozenLines.contains(m_frozenLinesCounter))
            {
                throw meshkernel::MeshKernelError("Frozen line counter in meshkernel state should not exist when restoring a deletion frozen line");
            }
            m_mkStateReference.m_frozenLines[m_frozenLinesCounter] = m_frozenLinePoints;
        }

        /// \brief Reference to the active MeshKernelState.
        MeshKernelState& m_mkStateReference;

        /// \brief A copy of the frozen line counter
        int m_frozenLinesCounter;

        /// \brief A copy of the frozen line points
        std::pair<meshkernel::Point, meshkernel::Point> m_frozenLinePoints;
    };

} // namespace meshkernelapi
