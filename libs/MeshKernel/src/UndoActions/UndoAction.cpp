//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2025.
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

#include "MeshKernel/UndoActions/UndoAction.hpp"

const std::string& meshkernel::UndoAction::to_string(const State state)
{
    static std::string committedStr = "Committed";
    static std::string restoredStr = "Restored";
    static std::string unknownStr = "UNKNOWN";

    switch (state)
    {
    case State::Committed:
        return committedStr;
    case State::Restored:
        return restoredStr;
    default:
        return unknownStr;
    }
}

meshkernel::UndoAction::State meshkernel::UndoAction::GetState() const
{
    return m_state;
}

void meshkernel::UndoAction::Commit()
{
    if (m_state == State::Restored)
    {
        DoCommit();
        m_state = State::Committed;
    }
}

void meshkernel::UndoAction::Restore()
{
    if (m_state == State::Committed)
    {
        DoRestore();
        m_state = State::Restored;
    }
}

std::uint64_t meshkernel::UndoAction::MemorySize() const
{
    return sizeof(*this);
}
