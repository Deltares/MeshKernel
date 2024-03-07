#include "MeshKernel/CompoundUndoAction.hpp"
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
        if (action->State() == UndoAction::Restored)
        {
            throw ConstraintError("Cannot add an action in the {} state.", UndoAction::to_string(action->State()));
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

std::uint64_t meshkernel::CompoundUndoAction::MemorySize() const
{
    std::uint64_t size = sizeof(*this) + m_undoActions.size() * sizeof(UndoActionPtr);

    for (const UndoActionPtr& action : m_undoActions)
    {
        size += action->MemorySize();
    }

    return size;
}

void meshkernel::CompoundUndoAction::Print(std::ostream& out) const
{
    out << "CompoundUndoAction: " << m_undoActions.size() << std::endl;

    UInt count = 0;

    for (const UndoActionPtr& action : m_undoActions)
    {
        out << "action: " << std::setw(4) << count << ": ";
        action->Print(out);
        ++count;
    }
}
