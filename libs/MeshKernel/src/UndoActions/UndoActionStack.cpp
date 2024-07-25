#include "MeshKernel/UndoActions/UndoActionStack.hpp"
#include "MeshKernel/Exceptions.hpp"

#include <algorithm>
#include <utility>

const meshkernel::UInt meshkernel::UndoActionStack::MaxUndoSize = 50;

void meshkernel::UndoActionStack::Add(UndoActionPtr&& action, const int actionId)
{
    if (action != nullptr)
    {
        if (action->GetState() != UndoAction::State::Committed)
        {
            throw ConstraintError("Cannot add an action in the {} state.", UndoAction::to_string(action->GetState()));
        }

        m_committed.emplace_back(std::move(action), actionId);
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

meshkernel::UInt meshkernel::UndoActionStack::Size() const
{
    return static_cast<UInt>(m_committed.size() + m_restored.size());
}

meshkernel::UInt meshkernel::UndoActionStack::CommittedSize() const
{
    return static_cast<UInt>(m_committed.size());
}

meshkernel::UInt meshkernel::UndoActionStack::RestoredSize() const
{
    return static_cast<UInt>(m_restored.size());
}

bool meshkernel::UndoActionStack::Undo()
{
    bool didUndo = false;

    if (!m_committed.empty())
    {
        // Perform undo operation
        m_committed.back().m_undoAction->Restore();
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

    if (!m_restored.empty())
    {
        // Perform commit (redo) operation
        m_restored.back().m_undoAction->Commit();
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

meshkernel::UInt meshkernel::UndoActionStack::Remove(const int actionId)
{
    UInt removedCount = 0;

    auto hasMatchingActionId = [actionId](const UndoActionForMesh& action)
    { return action.m_actionId == actionId; };

    removedCount = static_cast<UInt>(m_committed.remove_if(hasMatchingActionId));
    removedCount += static_cast<UInt>(m_restored.remove_if(hasMatchingActionId));

    return removedCount;
}

std::uint64_t meshkernel::UndoActionStack::MemorySize() const
{
    const std::uint64_t committedSize = std::accumulate(m_committed.begin(), m_committed.end(), static_cast<std::uint64_t>(0u),
                                                        [](const std::uint64_t partialSum, const UndoActionForMesh& action)
                                                        { return partialSum + action.m_undoAction->MemorySize(); });
    const std::uint64_t restoredSize = std::accumulate(m_restored.begin(), m_restored.end(), static_cast<std::uint64_t>(0u),
                                                       [](const std::uint64_t partialSum, const UndoActionForMesh& action)
                                                       { return partialSum + action.m_undoAction->MemorySize(); });

    return sizeof(*this) + committedSize + restoredSize;
}
