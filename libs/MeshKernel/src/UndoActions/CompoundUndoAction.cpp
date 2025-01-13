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

#include "MeshKernel/UndoActions/CompoundUndoAction.hpp"
#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Exceptions.hpp"

#include <iomanip>
#include <ranges>

std::unique_ptr<meshkernel::CompoundUndoAction> meshkernel::CompoundUndoAction::Create()
{
    return std::make_unique<CompoundUndoAction>();
}

void meshkernel::CompoundUndoAction::Add(UndoActionPtr&& action)
{
    if (action != nullptr)
    {
        if (action->GetState() == UndoAction::State::Restored)
        {
            throw ConstraintError("Cannot add an action in the {} state.", UndoAction::to_string(action->GetState()));
        }

        m_undoActions.emplace_back(std::move(action));
    }
}

void meshkernel::CompoundUndoAction::DoCommit()
{
    for (const UndoActionPtr& action : m_undoActions)
    {
        action->Commit();
    }
}

void meshkernel::CompoundUndoAction::DoRestore()
{
    for (const UndoActionPtr& action : m_undoActions | std::views::reverse)
    {
        action->Restore();
    }
}

meshkernel::CompoundUndoAction::const_iterator meshkernel::CompoundUndoAction::begin() const
{
    return m_undoActions.begin();
}

meshkernel::CompoundUndoAction::const_iterator meshkernel::CompoundUndoAction::end() const
{
    return m_undoActions.end();
}

std::uint64_t meshkernel::CompoundUndoAction::MemorySize() const
{
    std::uint64_t size = sizeof(*this) + m_undoActions.capacity() * sizeof(UndoActionPtr);

    for (const UndoActionPtr& action : m_undoActions)
    {
        size += action->MemorySize();
    }

    return size;
}
