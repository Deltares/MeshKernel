#include "MeshKernel/UndoActionStack.hpp"
#include "MeshKernel/Exceptions.hpp"

#include <algorithm>

const meshkernel::UInt meshkernel::UndoActionStack::MaxUndoSize = 10;

void meshkernel::UndoActionStack::Add(UndoActionPtr&& action)
{
    if (action != nullptr)
    {
        if (action->State() != UndoAction::Committed)
        {
            throw ConstraintError("Cannot add an action in the {} state.", UndoAction::to_string(action->State()));
        }

        m_committed.emplace_back(std::move(action));
        // Clear the restored actions. Adding a new undo action means that they cannot be re-done
        m_restored.clear();

        if (m_committed.size() > MaxUndoSize)
        {
            // If the number of undo-actions is greater than the maximum, then remove the first item in the list.
            m_committed.pop_front();
        }
    }
    else
    {
        // Logging message
    }
}

bool meshkernel::UndoActionStack::Undo()
{
    bool didUndo = false;

    if (m_committed.size() > 0)
    {
        // Perform undo operation
        m_committed.back()->Restore();
        // Now move to restored stack
        m_restored.emplace_back(std::move(m_committed.back()));
        m_committed.pop_back();
        didUndo = true;
    }

    return didUndo;
}

bool meshkernel::UndoActionStack::Commit()
{
    bool didCommit = false;

    if (m_restored.size() > 0)
    {
        // Perform commit (redo) operation
        m_restored.back()->Commit();
        // Now move to committed stack
        m_committed.emplace_back(std::move(m_restored.back()));
        m_restored.pop_back();
        didCommit = true;
    }

    return didCommit;
}

void meshkernel::UndoActionStack::Clear()
{
    m_committed.clear();
    m_restored.clear();
}

std::uint64_t meshkernel::UndoActionStack::MemorySize() const
{
    std::uint64_t committedSize = std::accumulate(m_committed.begin(), m_committed.end(), 0, [](const std::uint64_t partialSum, const std::unique_ptr<UndoAction>& action)
                                                  { return partialSum + action->MemorySize(); });
    std::uint64_t restoredSize = std::accumulate(m_restored.begin(), m_restored.end(), 0, [](const std::uint64_t partialSum, const std::unique_ptr<UndoAction>& action)
                                                 { return partialSum + action->MemorySize(); });

    return sizeof(*this) + committedSize + restoredSize;
}
