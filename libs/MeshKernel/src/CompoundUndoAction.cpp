#include "MeshKernel/CompoundUndoAction.hpp"

#include <ranges>

std::unique_ptr<meshkernel::CompoundUndoAction> meshkernel::CompoundUndoAction::Create ()
{
    return std::make_unique<CompoundUndoAction>();
}

void meshkernel::CompoundUndoAction::Add(UndoActionPtr&& action)
{
    m_undoActions.emplace_back(std::move(action));
}

void meshkernel::CompoundUndoAction::DoCommit()
{
    for (UndoActionPtr& action : m_undoActions)
    {
        action->Commit ();
    }
}

void meshkernel::CompoundUndoAction::DoRestore()
{
    for (UndoActionPtr& action : m_undoActions | std::views::reverse)
    {
        action->Restore ();
    }
}
